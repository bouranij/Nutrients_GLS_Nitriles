#Script for analysis of 16S data published in Nutrients by Bouranis et al. (2021)
#This script was written by John Bouranis

#Coding Tools
library(tidyverse)
library(magrittr)
#Micorbiome Tools
library(phyloseq)
library(vegan)
library(dunn.test)
#library(DESeq2)
#Visualization Tools
library(DT)
library(cowplot)
library(plotly)
library(RColorBrewer)

#Set up color palette:
paltrt <- RColorBrewer::brewer.pal(name = 'Dark2', n = 6)
palsmp <- viridis::viridis(n = 10)
palgrp <- c('#A6611A', '#018571')

#Load in the data
asvtab <- readRDS('./data/asv_tab.RDS')
taxatab <- readRDS('./data/tax_tab.RDS')
metadata <- read.csv('./data/microbiome_metadata.csv')
metabinterest <- read_csv('./data/metab_interest.csv')
rawdata_atwell <- read_csv('./data/atwell_nitriles_rerun.csv')
volumes <- read_csv('./data/atwell_volumes.csv') 

metadata %<>% mutate(treatment = ifelse(treatment == 'no_veg', 'NC', treatment))

#Factorize and set the levels on the metadata
metadata$treatment %<>% factor(levels = c('fecal_stock', 'NC', 'broc', 'brus', 'combo', 'control_digest'))
metadata$fecal_sample %<>% factor(levels = c('T5631','T5632','T6260','T6291','T4669','T1995','T5627','T5717','T5854','T6382')) 

#Set up colors for each treatment
names(paltrt) <- levels(metadata$treatment)
colTrt <- scale_colour_manual(name = "treatment", values = paltrt)

names(palsmp) <- levels(metadata$fecal_sample)
colSmp <- scale_colour_manual(name = "fecal_sample", values = palsmp)

names(palgrp) <- levels(metadata$group)
colGrp <- scale_color_manual(name = "group", values = palgrp)


#Construct the phyloseq object
rownames(metadata) <- metadata$sample
rownames(asvtab) <-metadata$sample
ps_raw <- phyloseq(otu_table(asvtab, taxa_are_rows = FALSE),
               sample_data(metadata),
               tax_table(taxatab))

#Give arbitrary names to the taxa as opposed to keeping as just DNA-sequences which identify them
taxa_names(ps_raw) <- paste0("ASV", seq(ntaxa(ps_raw)))

#Fill in missing genus names:
renames <- rownames(tax_table(ps_raw)[is.na(tax_table(ps_raw)[, 'Genus'])])
taxdf <- tax_table(ps_raw)[renames,]
renamed_genus <- unname(sapply(taxa_names(taxdf), function(x) paste0('f_', taxdf[x, 'Family'], '_', x)))
tax_table(ps_raw)[renames, 'Genus'] <- renamed_genus

#Remove the control digests, these are not relevant to our analysis
ps_raw <- ps_raw %>% subset_samples(treatment != 'control_digest')
#Agglomerate to the genus level
ps_genera <- ps_raw %>% tax_glom(taxrank = "Genus")
#Remove taxa not seen more than 3 times in at least 20% of the samples
ps_counts <- ps_genera %>% filter_taxa(function(x) sum(x > 3) > (0.2*length(x)), TRUE)
#Convert from counts to relative abundance
ps_relab <- ps_counts %>% transform_sample_counts(function(x) x / sum(x) )
#Filter out low abundance (>1e-5) taxa
ps <- ps_relab %>% filter_taxa(function(x) mean(x) > 1e-5, TRUE)

#Number of ASVs of each
summary <- data.frame(
  list('ps_raw' = dim(tax_table(ps_raw))[1],
  'ps_genera' = dim(tax_table(ps_genera))[1],
  'ps_counts' = dim(tax_table(ps_counts))[1],
  'ps_relab' = dim(tax_table(ps_relab))[1],
  'ps' = dim(tax_table(ps))[1]
))
rownames(summary) <- 'Number of ASVs'

knitr::kable(summary, caption = 'ASVs across data preprocessing')

#Extract relevant sub-data:
ps_stock <- ps %>% subset_samples(treatment == 'fecal_stock')
ps_broc <- ps %>% subset_samples(treatment == 'broc')
ps_bnv <- ps %>% subset_samples(treatment %in% c('NC', 'broc'))
ps_fbnv <- ps %>% subset_samples(treatment %in% c('fecal_stock', 'NC', 'broc'))


datatable(data.frame(sample_data(ps_fbnv)),
           filter = 'top',
           extensions = 'Buttons',
           options = list(pageLength = 12,
                          dom = 'Bfrtip',
                          scrollX = TRUE))      

#Fecal-stock composition:
stocks <- plot_bar(ps_stock, fill = 'Family', x = 'fecal_sample') + 
  guides(col = guide_legend(ncol = 2)) +
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5)) + 
  xlab('\nStock') + 
  ylab('Relative Abundance') +
  ggtitle('Relative Abundance of Bacterial Taxa in Fecal Stocks')

ggplotly(stocks, width = 900, height = 600)


#Alpha diversity metrics:
adiv_stock <- ps_raw %>% 
  subset_samples(treatment == 'fecal_stock') %>%
  estimate_richness(measures = c('Observed','Shannon','Simpson')) %>%
  cbind(., data.frame(sample_data(ps_raw %>% 
                         subset_samples(treatment == 'fecal_stock')))%>%
                         dplyr::select(fecal_sample)) %>%
  remove_rownames() %>%
  column_to_rownames('fecal_sample')

#apply(adiv_stock , 2, mean)
adiv_stock

#Plot of alpha-diversity metrics:

gstock <- adiv_stock %>%
  rownames_to_column('fecal_sample') %>%
  pivot_longer(where(is.numeric), names_to = 'measure')

stk <- ggplot(gstock, aes(x = fecal_sample, y = value, color = fecal_sample)) +
  geom_point() +
  facet_wrap(vars(measure), scales = 'free', ncol = 3) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90), axis.title.x = element_blank()) 
stk

#Significance testing of alpha-diversity metrics:

st <- adiv_stock %>%
  rownames_to_column('fecal_sample') %>%
  left_join(.,data.frame(sample_data(ps_stock))) %>%
  select(fecal_sample, Observed, Shannon, Simpson, birth_year, biological_sex) %>%
  mutate(by_bin = ifelse(birth_year <= 1975, 1, 
                         ifelse(birth_year <= 1985, 2, 
                                ifelse(birth_year <= 1995, 3, 4)))) %>%
  pivot_longer(cols = c('Observed', 'Shannon', 'Simpson'), names_to = 'measure') %>%
  group_by(measure) %>%
  nest() %>%
  mutate(kw_b = map(data, function(x) kruskal.test(value ~ as.factor(by_bin), data = x) %>% broom::tidy() %>% pull('p.value'))) %>%
  mutate(kw_s = map(data, function(x) kruskal.test(value ~ biological_sex, data = x) %>% broom::tidy() %>% pull(2))) %>%
  unnest(cols = c('kw_b', 'kw_s'))

stc <- st %>%
  select(-data) %>%
  rename('Biological Sex' = kw_s, 'Age' = kw_b, 'Measure' = measure)
stc

#Fecal stock beta-diversity:

ordBCstock <- ordinate(ps_stock, method = 'PCoA', distance = 'bray')
sbd <- plot_ordination(ps_stock, ordBCstock, color = 'fecal_sample') +
  geom_point(aes(text = paste0('Fecal Sample: ', fecal_sample, '\n',
                               'Biological Sex: ', biological_sex, '\n',
                               'Year of Bith: ', birth_year)),
             size = 3) +
  theme_cowplot() +
  colSmp +
  ggtitle('Beta Diversity of Fecal Stocks')

sbd

#PERMANOVA by Age:

ordBCstock <- ordinate(ps_stock, method = 'PCoA', distance = 'bray')
sbd <- plot_ordination(ps_stock, ordBCstock, justDF=T) 
distmat <- phyloseq::distance(ps_stock, method = 'bray')
#Extract the metadata
mdata <- data.frame(sample_data(ps_stock)) %>% dplyr::select(fecal_sample, treatment, biological_sex, birth_year)  %>%
  mutate(by_bin = ifelse(birth_year <= 1975, 1, 
                         ifelse(birth_year <= 1985, 2, 
                                ifelse(birth_year <= 1995, 3, 4)))) %>%
  modify_at('by_bin', as.factor)

#Set the seed for reproducibility
set.seed(120)
#Run betadisper to verify the distrubtion of our groups is equal, an underlying assumption of PERMANOVA
bdisp <- betadisper(distmat, mdata$biological_sex)
#Evaluate using a permutation test
permutest(bdisp) 
adonis2(distmat~biological_sex, data = mdata)

#PERMANOVA by Birth Year:

bdisp <- betadisper(distmat, mdata$by_bin)
#Evaluate using a permutation test
permutest(bdisp) 
adonis2(distmat~by_bin, data = mdata)

#Composition of stocks after incubation:

plot_bar(ps_fbnv, fill = 'Family', x = 'treatment') + 
  facet_wrap(~fecal_sample, ncol = 5) + 
  theme_cowplot() + theme(axis.text.x = element_text(angle = 45, hjust = 0.5, size = 12), 
                          axis.text.y = element_text(size = 12),
                          axis.title.x = element_blank(),
                          axis.title.y = element_blank()) +
  ggtitle('Relative Abundance of Bacterial Taxa Following Incubation')


#Dunn Test to compare alpha diversity between pre- and post-incubation treatments
adiv_dunn <- ps_raw %>% 
  subset_samples(treatment %in% c('fecal_stock', 'broc', 'NC')) %>%
  estimate_richness(measures = c('Observed','Shannon','Simpson')) %>%
  cbind(., data.frame(sample_data(ps_raw %>% 
                                    subset_samples(treatment %in% c('fecal_stock', 'broc', 'NC'))))) %>%
  pivot_longer(cols = c('Observed', 'Shannon', 'Simpson'), names_to = 'measure') %>%
  dplyr::select(measure, value, treatment, fecal_sample) %>%
  group_by(measure) %>%
  nest() 
walk2(.x = adiv_dunn$data, .y = list('Observed', 'Shannon', 'Simpson'), function(x, y) {
  print(y)
  dunn.test(x$value, g = x$treatment, kw = T, method = 'bh')})


# Beta-diversity following incubation:

ordBCall <- ordinate(ps_fbnv, method = 'PCoA', distance = 'bray')
allbd <- plot_ordination(ps_fbnv, ordBCall, color = 'treatment') +
  geom_point(aes(text = paste0('Fecal Sample: ', fecal_sample, '\n',
                               'Treatment: ', treatment)),
             size = 3) +
  theme_cowplot() +
  colTrt +
  ggtitle('Beta Diversity of All Samples')
allbd


#Plotting total within-cluster sum of squares:
ugh <- ordinate(ps_fbnv, method = 'PCoA', distance = 'bray') %>%
  plot_ordination(ps_fbnv, ., justDF = TRUE) 


map_df(2:10, function(x) kmeans(ugh[,1:2], x, nstart = 25) %>% extract('tot.withinss')) %>%
  cbind(2:10, .) %>%
    plot(type = 'b', xlab = 'Number of Clusters K', ylab = 'Total Within-Clusters Sum of Squares')


#Plotting with our confirmed clusters from k-means
ughc <- kmeans(ugh[,1:2], 3, nstart = 25) %>% 
  extract('cluster') %>%
  cbind(ugh, .) %>%
  modify_at('cluster', as.factor)


c <- ggplot(ughc, aes(x = Axis.1, y = Axis.2)) + 
  stat_ellipse(aes(group = cluster, fill = cluster), geom = 'polygon', color = 'black', alpha = 0.37, show.legend = FALSE) +
  #geom_point(aes(fill = treatment), color = 'black', pch = 21, size = 4, show.legend = TRUE) +
  xlab('Axis 1 [50.1%]') +
  ylab('Axis 2 [26.5%]') +
  theme_cowplot() + 
  theme(legend.position = 'none') +
  scale_x_continuous(limits = c(-0.8, 0.7)) +
  scale_y_continuous(limits = c(-0.8, 0.5)) +
  scale_fill_manual(values = c('#A50104', '#7DAF9C','#FCBA04')) +
  annotate('text', label = c('Fecal Stock Cluster', 'Mixed Cluster 1', 'Mixed Cluster 2'), 
           x = c(-0.5, 0.45, 0.1), y = c(0.3, 0.35, -0.8),
           size = 5)
  
trt <- RColorBrewer::brewer.pal(n = 3, name = 'Accent')

d <- ggplot(ughc, aes(x = Axis.1, y = Axis.2)) + 
  #stat_ellipse(aes(group = cluster, fill = cluster), geom = 'polygon', color = 'black', alpha = 0.37, show.legend = FALSE) +
  geom_point(aes(fill = treatment), color = 'black', pch = 21, size = 4, show.legend = TRUE) +
  xlab('Axis 1 [50.1%]') +
  ylab('Axis 2 [26.5%]') +
  theme_cowplot() + 
  scale_x_continuous(limits = c(-0.8, 0.7)) +
  scale_y_continuous(limits = c(-0.8, 0.5)) +
  scale_fill_manual(values = c(trt[2], trt[3], trt[1]), 
                    labels = c('Fecal Stock', 'NC', 'Broc'),
                    name = 'Treatment') +
  theme(legend.position = c(0.02,0.15))


aligned_plots <- align_plots(d, c, align = 'hv', axis = 'tbl')

ap <- ggdraw(aligned_plots[[2]]) + draw_plot(aligned_plots[[1]])

ap


#Removing the pre-incubation samples and plotting again:
s1 <- ordinate(ps_bnv, method = 'PCoA', distance = 'bray') %>%
  plot_ordination(ps_bnv, ., justDF = TRUE) 

map_df(1:9, function(x) kmeans(s1[,1:2], x, nstart = 25) %>% extract('tot.withinss')) %>%
  cbind(1:9, .) %>%
    plot(type = 'b', xlab = 'Number of Clusters K', ylab = 'Total Within-Clusters Sum of Squares')

s1c <- kmeans(s1[,1:2], 2, nstart = 25) %>% 
  extract('cluster') %>%
  cbind(s1, .) 
s1c$cluster %<>% as.factor()


sc <- ggplot(s1c, aes(x = Axis.1, y = Axis.2)) + 
  stat_ellipse(aes(group = cluster, fill = cluster), geom = 'polygon', color = 'black', alpha = 0.37, show.legend = FALSE) +
  #geom_point(aes(fill = treatment), color = 'black', pch = 21, size = 4, show.legend = TRUE) +
  xlab('Axis 1 [69.2%]') +
  ylab('Axis 2 [22.3%]') +
  theme_cowplot() + 
  theme(legend.position = 'none') +
  scale_x_continuous(limits = c(-0.5, 1)) +
  scale_y_continuous(limits = c(-1, 0.8)) +
  scale_fill_manual(values = rev(c('#A50104', '#7DAF9C'))) +
  annotate('text', label = c('Mixed Cluster 1', 'Mixed Cluster 2'), 
           x = c(-0.3, 0.3), y = c(0.25, -0.6),
           size = 5)

sd <- ggplot(s1c, aes(x = Axis.1, y = Axis.2)) + 
  #stat_ellipse(aes(group = cluster, fill = cluster), geom = 'polygon', color = 'black', alpha = 0.37, show.legend = FALSE) +
  geom_point(aes(fill = treatment), color = 'black', pch = 21, size = 4, show.legend = TRUE) +
  xlab('Axis 1 [69.2%]') +
  ylab('Axis 2 [22.3%]') +
  theme_cowplot() + 
  scale_x_continuous(limits = c(-0.5, 1)) +
  scale_y_continuous(limits = c(-1, 0.8)) +
  scale_fill_manual(values = c(trt[3], trt[1]), 
                    labels = c('NC', 'Broc'),
                    name = 'Treatment') +
  theme(legend.position = c(0.02,0.1))
sap <- align_plots(sd, sc, align = 'hv', axis = 'tbl')

ggdraw(sap[[2]]) + draw_plot(sap[[1]])


#PERMANOVA by treatment post-incubation

#First calculate the Bray-Curtis distances of the samples
distmat <- phyloseq::distance(ps_bnv, method = 'bray')
#Extract the metadata
mdata <- data.frame(sample_data(ps_bnv)) %>% dplyr::select(fecal_sample, treatment)
#Set the seed for reproducibility
set.seed(120)
#Run betadisper to verify the distrubtion of our groups is equal, an underlying assumption of PERMANOVA
bdisp <- betadisper(distmat, mdata$treatment)
#Evaluate using a permutation test
permutest(bdisp)
adonis2(distmat~treatment, data = mdata)
plot(bdisp, label = F, col = paltrt[2:5]) %$%
legend('bottomright', legend = c('NC', 'broc'), col = paltrt[2:3], pch = 19) %$%
points(bdisp$centroids, pch = 20, cex = 2, col = paltrt[2:3])

#PERMANOVA by fecal donor post-incubation
#First calculate the Bray-Curtis distances of the samples
distmat <- phyloseq::distance(ps_bnv, method = 'bray')
#Extract the metadata
mdata <- data.frame(sample_data(ps_bnv)) %>% dplyr::select(metabolomics_pos_sample, fecal_sample, treatment)
#Set the seed for reproducibility
set.seed(120)
#Run betadisper to verify the distrubtion of our groups is equal, an underlying assumption of PERMANOVA
bdisp <- betadisper(distmat, mdata$fecal_sample) #Significant but due PERMANOVA is robust to this when there is a balanced design (like what we have)

#Evaluate using a permutation test
permutest(bdisp)
adonis2(distmat~fecal_sample, data = mdata)

plot(bdisp, label = F, col = palsmp) %$%
legend('bottomright', legend = names(palsmp), col = palsmp, pch = 19) %$%
points(bdisp$centroids, pch = 20, cex = 2, col = palsmp)

#Relative abundance of NITs in each cluster:

mdataveg <- data.frame(sample_data(ps_broc)) %>% dplyr::select(metabolomics_pos_sample, fecal_sample, treatment)
metabinterest <- read_csv('./data/metab_interest.csv') %>%
  right_join(., mdataveg, by = c('sample' = 'metabolomics_pos_sample'))

nitdata <- metabinterest %>% #Quantified IBN, IBNNIT, and SFNNIT from multiquant
  pivot_longer(where(is.numeric), names_to = 'metabolite', values_to = 'intensity') %>%
  left_join(., data.frame(sample_data(ps_broc)) %>% select(fecal_sample, group)) %>%
  filter(treatment == 'broc')

kt <- nitdata %>%
  group_by(metabolite) %>%
  nest() %>%
  mutate('p-value' = map_dbl(data, function(x) kruskal.test(intensity ~ group, data = x) %>% broom::tidy() %>% pull(p.value))) %>%
  dplyr::select(-data) %>%
  rename('Metabolite' = metabolite)
  
knitr::kable(kt, caption = 'Differences in Means NIT between Clusters')


SEM <- function(x) sd(x)/sqrt(length(x))

nit2 <- nitdata %>%
  group_by(metabolite, group) %>%
  summarise(mean = mean(intensity), se = SEM(intensity)) %>%
  modify_at('metabolite', as.factor) %>%
  modify_at('metabolite', fct_relevel, c('SFNNIT', 'IBNNIT'))

sfnnit <- nit2 %>%
  filter(metabolite == 'SFNNIT') %>%
  modify_at('group', function(x) as.factor(ifelse(x == 'A', 1, 2)))

ibnnit <- nit2 %>%
  filter(metabolite == 'IBNNIT') %>%
  modify_at('group', function(x) as.factor(ifelse(x == 'A', 1, 2)))

s <- ggplot(sfnnit, aes(x = group, y = mean, group = group)) +
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se), width = 0.5) +
  geom_col(fill = 'darkgrey', color = 'black') +
  ylab('SFN-NIT (Relative Intensity)') +
  xlab('Cluster') +
  theme_cowplot() +
  scale_fill_manual(labels = c('Cluster 1', 'Cluster 2'),
                    name = 'Cluster') +
  scale_y_continuous(label = scales::label_comma(),
                     expand = c(0,0),
                     limits = c(0,172000)) +
  annotate('text', label = c('*'), 
           x = 2, y = 35000,
           size = 10) 

i <- ggplot(ibnnit, aes(x = group, y = mean, group = group)) +
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se), width = 0.5) +
  geom_col(fill = 'darkgrey', color = 'black') +
  ylab('IBN-NIT (Relative Intensity)') +
  xlab('Cluster') +
  theme_cowplot() +
  scale_fill_manual(labels = c('Cluster 1', 'Cluster 2'),
                    name = 'Cluster') +
  scale_y_continuous(label = scales::label_comma(),
                     expand = c(0,0),
                     limits = c(0,13000)) +
  annotate('text', label = c('*'), 
           x = 2, y = 6300,
           size = 10) 

plot_grid(s, i, labels = c('A','B'), axis = 'b')


#PCoA with NIT abundance overlaid 

ordBCveg <- ordinate(ps_broc, method = 'PCoA', distance = 'bray') %>%
  plot_ordination(ps_broc, ., justDF = TRUE) %>%
  left_join(., metabinterest, by = c('metabolomics_pos_sample' = 'sample')) %>%
  select(-ends_with('.y'))


vegbd <- ggplot(ordBCveg, aes(x = Axis.1, y = Axis.2, size = SFNNIT, color = IBNNIT)) +
  geom_point(aes(text = paste0('Fecal Sample: ', fecal_sample.x, '\n',
                               'Treatment: ', treatment.x, '\n',
                               'IBNNIT: ', IBNNIT, '\n',
                               'SFNNIT: ', SFNNIT))) +
  scale_color_viridis_c() +
  theme_cowplot() +
  ggtitle('Beta Diversity of Digests Following Incubation')

vegbd

#PERMANOVA on NIT abundance:

set.seed(101)
mdataveg <- data.frame(sample_data(ps_broc)) %>% dplyr::select(metabolomics_pos_sample, fecal_sample, treatment)
distmatveg <- phyloseq::distance(ps_broc, method = 'bray')
metabinterest <- read_csv('./data/metab_interest.csv') %>%
  right_join(., mdataveg, by = c('sample' = 'metabolomics_pos_sample'))

suppressWarnings(    #SuppressWarnings because broom::tidy() does not play well PERMANOVA outputs
perma <- metabinterest %>%
  pivot_longer(where(is.numeric), names_to = 'Metabolite', values_to = 'intensity') %>%
  dplyr::group_by(Metabolite) %>%
  nest() %>%
  mutate(permanova = map(data, function(x){ 
    adonis2(distmatveg~intensity, data = x) %>% 
      broom::tidy() %>%
      magrittr::extract('p.value') %>%
      drop_na() 
    })) %>%
  unnest(cols = 'permanova') %>% 
  select(-data)
)

#Negative Binomial Model:

library(lme4)
library(phia)
set.seed(123)

psglm <- ps_counts %>%
  subset_samples(treatment %in% c('broc', 'NC')) %>%
  rarefy_even_depth() %>%
  psmelt()

glmnb <- psglm %>%
  group_by(OTU) %>%
  nest() %>%
  mutate_at( # change the data in 'lmer' column
    "data",
    purrr::map, # allows iteration over data in 'lmer' column, returns list
    function(x) {
      lmer <- try(glmer.nb(Abundance ~ treatment*group + (1 |fecal_sample), data = x), silent = TRUE) # do the lmer test
      c(lmer, x) # return results of lmer, as well as data for (r)anova later
    } 
  )

#Remove the failed runs which contain error messages
glmnbclean <- glmnb[sapply(glmnb$data, function(x) class(x[[1]]) !=  'character'),]

#Function to calculate the log2FoldChange of each microbial species:
#Helper function
l2fcmin <- function(x, counts, factor, class1, class2){
  x[[counts]] <- log2(x[[counts]]+sqrt(nzmin(x[[counts]])*0.001))
  c1 <- x[[counts]][sapply(x[[factor]], function(x) x == class1)]
  c2 <- x[[counts]][sapply(x[[factor]], function(x) x == class2)]
  lf2 <- mean(c1[sapply(c1,is.finite)]) - mean(c2[sapply(c2, is.finite)])
  return(lf2)
}

nzmin <- function(x){
  min(x[x > 0])
}


lfc <- psglm %>%
  group_by(OTU, treatment) %>%
  nest() %>%
  mutate(L2FC = purrr::map(data, function(x) l2fcmin(x, counts = 'Abundance', factor = 'group', class1 = 'A', class2 = 'B'))) %>%
  dplyr::select(-data) %>%
  unnest(L2FC) %>%
  ungroup()

glmres <- glmnbclean %>%
  group_by(OTU) %>%
  group_modify(~{
    testInteractions(.x$data[[1]][[1]], fixed = 'treatment', across = 'group') %>% broom::tidy()
  }) %>%
  ungroup() %>%
  mutate(padj = p.adjust(p.value, method = 'BH')) %>%
  left_join(lfc, by = c('OTU' = 'OTU', 'term' = 'treatment'))

glmsig <- glmres %>%
  filter(padj <= 0.05) %>%
  left_join(data.frame(tax_table(ps_counts)) %>% rownames_to_column('OTU')) 

glmsig

#Spearman's Correlation Analysis:

psfam <- ps_broc %>%
  tax_glom(taxrank = 'Family') %>%
  psmelt()

cordata <- psfam %>%
  dplyr::rename('ASV' = 'OTU') %>%
  dplyr::select(ASV, Abundance, metabolomics_pos_sample, group, fecal_sample, Family) %>%
  dplyr::rename('sample' = 'metabolomics_pos_sample') %>%
  left_join(metabinterest) %>%
  pivot_longer(cols = c('IBNNIT', 'SFNNIT'), names_to = 'metabolite', values_to = 'intensity') %>%
  group_by(Family, metabolite) %>%
  nest() %>%
  mutate(rho = map(data, function(x){ 
    stats::cor(x$intensity, x$Abundance, method = 'spearman')
    })) %>%
  unnest(cols = 'rho') %>%
  ungroup() %>%
  dplyr::select(-data) %>%
  drop_na()

hmap <- ggplot(cordata, aes(x = Family, y = metabolite, fill = rho)) +
  geom_tile(aes(text = paste0('Metabolite: ', metabolite, '\n',
                                'Family: ', Family, '\n',
                                'Rho:',  round(rho,3)))) +
  coord_fixed() +
  scale_fill_gradientn(colors = c('blue', 'white', 'red')) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  ggtitle('Spearman Correlations Between Families and Metabolites for Broccoli')


#Human NIT data:

ad <- rawdata_atwell%>%
  dplyr::select(-phase, -sample_site) %>%
  group_by(subject, time) %>%
  mutate(SFN_tot = SFN+SFNNIT+SFN_CYS+SFN_NAC+SFN_CG+SFN_GSH) %>%
  modify_at(c('time','subject'), as.factor)

v <- volumes %>%
  #Convert mL to L
  mutate(across(c(2:7), ~.x/1000)) %>%
  pivot_longer(cols = -subject, names_to = 'time', values_to = 'volume') %>%
  modify_at('time', as.integer) %>%
  modify_at(c('time','subject'), as.factor) %>%
  #Add µM to volume data
  left_join(ad) %>%
  group_by(subject, time) %>%
  #Multiple each µM by volume of urine to get µmol
  mutate(across(starts_with('SFN'), ~.x*volume)) 

bleh <- v %>%
  dplyr::select(subject, time, SFNNIT, SFN_NAC) %>%
  pivot_longer(where(is.numeric), names_to = 'metabolite', values_to = 'umol')

plot <- ggplot(bleh, aes(x = time, y = umol, color = subject, group = subject)) +
  geom_bar(aes(group = time),stat = 'summary', fun = mean, fill = 'grey', color = 'black', alpha = 0.5) +
  geom_point(size = 3, aes(text = paste0('Metabolite: ', metabolite, '\n',
                                         'µmol: ', umol, '\n',
                                         'Subject: ', subject, '\n',
                                         'Time: ', time, 'h'))) +
  geom_path(size = 1.5) + 
  theme_minimal_hgrid() +
  ylab('SFN-NAC (µmol)') +
  xlab('Time (Hours)') + 
  scale_y_continuous(limits = c(0, 65), expand = c(0.00,0.00), breaks = seq(0, 70, by = 5)) +
  scale_x_discrete(expand = c(0.005,0.005)) +
  scale_color_manual(name = 'Subject ID', 
                     values = c('#F8766D', '#D89000', '#A3A500', '#39B600', '#00BF7D', '#00BFC4', '#00B0F6', '#9590FF', '#E76BF3', '#FF62BC')) +
  facet_wrap(~metabolite)

plot





