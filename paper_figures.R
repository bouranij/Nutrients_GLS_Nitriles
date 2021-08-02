#Poop Soup Nitriles Figures:
library(tidyverse)
library(magrittr)
library(phyloseq)
library(scales)
library(cowplot)
setwd('~/Documents/Projects/poopsoup/')

#Data Load In
asvtab <- readRDS('./data/microbiome/asv_tab.RDS')
taxatab <- readRDS('./results/data/microbiome/tax_tab.RDS')
metadata <- read.csv('./results/data/microbiome/microbiome_metadata.csv')
metabinterest <- read_csv('./results/metab_interest.csv')
rawdata_atwell <- read_csv('./Atwell_data/atwell_nitriles_rerun.csv')
volumes <- read_csv('./Atwell_data/atwell_volumes.csv') 
#Factorize and set the levels on the metadata
metadata$treatment %<>% factor(levels = c('fecal_stock', 'no_veg', 'broc', 'brus', 'combo', 'control_digest'))
metadata$fecal_sample %<>% factor(levels = c('T5631','T5632','T6260','T6291','T4669','T1995','T5627','T5717','T5854','T6382')) 

#Construct the PhyloSeq Object
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

ps_stock <- ps %>% subset_samples(treatment == 'fecal_stock')
ps_broc <- ps %>% subset_samples(treatment == 'broc')
ps_bnv <- ps %>% subset_samples(treatment %in% c('no_veg', 'broc'))
ps_fbnv <- ps %>% subset_samples(treatment %in% c('fecal_stock', 'no_veg', 'broc'))

#Figure 1
fbnv <- ps_fbnv %>%
  tax_glom('Family') %>%
  psmelt() %>%
  modify_at('treatment', as.character) %>%
  mutate(treatment = gsub('_','', x = treatment)) %>%
  mutate(samptreat = paste0(treatment, '_', fecal_sample)) %>%
  arrange(Family, fecal_sample)
  
goodlevels <- c(c(fbnv$samptreat %>% unique() %>% grep('stock',., value = T)), 
                c(fbnv$samptreat %>% unique() %>% grep('no',., value = T)),
                c(fbnv$samptreat %>% unique() %>% grep('broc',., value = T)))
fbnv$samptreat %<>% factor(., levels = goodlevels)
fss <- sapply(goodlevels, function(x) str_split(x, '_')[[1]][2])


h1 <- ggplot(fbnv,aes(y = Abundance, x = samptreat, fill = Family, group = treatment)) + 
  geom_bar(stat = 'identity', color = 'black') + 
  scale_x_discrete(labels = fss) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, size = 25, vjust = 0.5),
        axis.text.y = element_text(size = 25, vjust = 0.8),
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 28),
        axis.title.x = element_text(size = 30, vjust = 0), 
        axis.title.y = element_text(size = 30 , vjust = 1)) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_percent()) +
  scale_fill_manual(values = colorRamps::primary.colors()) +
  labs(x = 'Human Stool/Fecal Subject ID',
       y = 'Microbial Abundance') 
leg <- plot_grid(get_legend(h1))

h2 <- ggplot(fbnv) + 
  geom_bar(aes(x = samptreat, y = 1, fill = treatment), stat = 'identity', width = 1) +
  theme_void() +
  theme(legend.position = 'none') + 
  annotate('text', x = 5.5, y = 10, label = 'Fecal Stock', size = 12) +
  annotate('text', x = 15.5, y = 10, label = 'Negative Control', size = 12) +
  annotate('text', x = 25.5, y = 10, label = 'Broccoli', size = 12) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 3, name = 'Accent'))

h1 <- h1 + theme(legend.position = 'none')
plot <- plot_grid(h2, h1, align = "v", ncol = 1, axis = "tb", rel_heights = c(0.7, 15))
plot_grid(plot,leg, nrow = 1, rel_widths = c(10, 1.5))

#Figure 2: save as 2100w x 2320h

ugh <- ordinate(ps_fbnv, method = 'PCoA', distance = 'bray') %>%
  plot_ordination(ps_fbnv, ., justDF = TRUE) 

map_df(1:20, function(x) kmeans(ugh[,1:2], x, nstart = 25) %>% extract('tot.withinss')) %>%
  cbind(1:20, .) %>%
    plot(type = 'b', xlab = 'Number of Clusters K', ylab = 'Total Within-Clusters Sum of Squares')

ughc <- kmeans(ugh[,1:2], 3, nstart = 25) %>% 
  extract('cluster') %>%
  cbind(ugh, .) 
ughc$cluster %<>% as.factor()

trt <- RColorBrewer::brewer.pal(n = 3, name = 'Accent')

c <- ggplot(ughc, aes(x = Axis.1, y = Axis.2)) + 
  stat_ellipse(aes(group = cluster, fill = cluster), geom = 'polygon', color = 'black', alpha = 0.37, show.legend = FALSE, size = 1.5) +
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
           size = 15) +
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 35))

d <- ggplot(ughc, aes(x = Axis.1, y = Axis.2)) + 
  #stat_ellipse(aes(group = cluster, fill = cluster), geom = 'polygon', color = 'black', alpha = 0.37, show.legend = FALSE) +
  geom_point(aes(fill = treatment), color = 'black', pch = 21, size = 15, show.legend = TRUE) +
  xlab('Axis 1 [50.1%]') +
  ylab('Axis 2 [26.5%]') +
  theme_cowplot() + 
  scale_x_continuous(limits = c(-0.8, 0.7)) +
  scale_y_continuous(limits = c(-0.8, 0.5)) +
  scale_fill_manual(values = c(trt[2], trt[3], trt[1]), 
                    labels = c('Fecal Stock', 'NC', 'Broc'),
                    name = 'Treatment') +
  theme(legend.position = c(0.02,0.1),
        legend.text = element_text(size = 28),
        legend.title = element_text(size = 30),
        axis.text = element_text(size = 30),
        axis.title = element_text(size = 35))
aligned_plots <- align_plots(d, c, align = 'hv', axis = 'tbl')

ggdraw(aligned_plots[[2]]) + draw_plot(aligned_plots[[1]])

#Figure 3: Save as 1500w x 1981h
mdataveg <- data.frame(sample_data(ps_broc)) %>% dplyr::select(metabolomics_pos_sample, fecal_sample, treatment)
metabinterest <- read_csv('./results/metab_interest.csv') %>%
  right_join(., mdataveg, by = c('sample' = 'metabolomics_pos_sample'))

nitdata <- metabinterest %>% #Quantified IBN, IBNNIT, and SFNNIT from multiquant
  pivot_longer(where(is.numeric), names_to = 'metabolite', values_to = 'intensity') %>%
  left_join(., data.frame(sample_data(ps_broc)) %>% select(fecal_sample, group)) %>%
  filter(treatment == 'broc')

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
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se), width = 0.5, size = 1.5) +
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
           size = 20) +
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 35))

i <- ggplot(ibnnit, aes(x = group, y = mean, group = group)) +
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se), width = 0.5, size = 1.5) +
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
           size = 20) +
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 35))

plot_grid(s, i, labels = c('A','B'), axis = 'b', label_size = 30)

#Figure 4: Save as 3000w x 1285h
ad <- rawdata_atwell%>%
  dplyr::select(-phase, -sample_site) %>%
  group_by(subject, time) %>%
  mutate(SFN_tot = SFN+SFNNIT+SFN_CYS+SFN_NAC+SFN_CG+SFN_GSH) %>%
  modify_at(c('time','subject'), as.factor)

#Convert µM to µmol
#Load in urine volumes
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

#Plot the results
pnit <- ggplot(v, aes(x = time, y = SFNNIT, color = subject, group = subject)) +
  #Calculate the means as solid bars
  geom_bar(aes(group = time),stat = 'summary', fun = mean, fill = 'grey', color = 'black', alpha = 0.5) +
  #Plot each individual as a point and path
  geom_point(size = 9) +
  geom_path(size = 4.5) + 
  theme_minimal_hgrid() +
  ylab('SFNNIT (µmol)') +
  xlab('Time (Hours)')  +
  scale_y_continuous(limits = c(0, 65), expand = c(0.0,0.00), breaks = seq(0, 70, by = 5)) +
  scale_x_discrete(expand = c(0.005,0.005)) +
  scale_color_manual(name = 'Subject ID', 
                     values = c('#F8766D', '#D89000', '#A3A500', '#39B600', '#00BF7D', '#00BFC4', '#00B0F6', '#9590FF', '#E76BF3', '#FF62BC')) +
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 35),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 25),
        plot.margin = margin(2,2,2,2, 'cm'))

pnac <- ggplot(v, aes(x = time, y = SFN_NAC, color = subject, group = subject)) +
  geom_bar(aes(group = time),stat = 'summary', fun = mean, fill = 'grey', color = 'black', alpha = 0.5) +
  geom_point(size = 9) +
  geom_path(size = 4.5) + 
  theme_minimal_hgrid() +
  ylab('SFN-NAC (µmol)') +
  xlab('Time (Hours)') + 
  theme(legend.position = 'none') +
  scale_y_continuous(limits = c(0, 65), expand = c(0.00,0.00), breaks = seq(0, 70, by = 5)) +
  scale_x_discrete(expand = c(0.005,0.005)) +
  scale_color_manual(name = 'Subject ID', 
                     values = c('#F8766D', '#D89000', '#A3A500', '#39B600', '#00BF7D', '#00BFC4', '#00B0F6', '#9590FF', '#E76BF3', '#FF62BC'))+
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 35),
        plot.margin = margin(2,2,2,2, 'cm'))
#Put them together into one plot
plot_grid(pnac, pnit, labels = c('A','B'), nrow = 1, label_size = 30)

#Supplemental Figure 1
s1 <- ordinate(ps_bnv, method = 'PCoA', distance = 'bray') %>%
  plot_ordination(ps_bnv, ., justDF = TRUE) 

map_df(1:9, function(x) kmeans(s1[,1:2], x, nstart = 25) %>% extract('tot.withinss')) %>%
  cbind(1:9, .) %>%
    plot(type = 'b', xlab = 'Number of Clusters K', ylab = 'Total Within-Clusters Sum of Squares')

s1c <- kmeans(s1[,1:2], 2, nstart = 25) %>% 
  extract('cluster') %>%
  cbind(s1, .) 
s1c$cluster %<>% as.factor()

ggplot(s1c, aes(x = Axis.1, y = Axis.2, color = treatment)) + 
  geom_point() +
  stat_ellipse(aes(group = cluster, fill = cluster), geom = 'polygon', color = 'black', alpha = 0.37, show.legend = FALSE) 

ordinate(ps_bnv, method = 'PCoA', distance = 'bray') %>%
  plot_ordination(ps_bnv, ., justDF = F) 

sc <- ggplot(s1c, aes(x = Axis.1, y = Axis.2)) + 
  stat_ellipse(aes(group = cluster, fill = cluster), geom = 'polygon', color = 'black', alpha = 0.37, show.legend = FALSE, size = 1.5) +
  #geom_point(aes(fill = treatment), color = 'black', pch = 21, size = 4, show.legend = TRUE) +
  xlab('Axis 1 [69.2%]') +
  ylab('Axis 2 [22.3%]') +
  theme_cowplot() + 
  theme(legend.position = 'none') +
  scale_x_continuous(limits = c(-0.5, 1)) +
  scale_y_continuous(limits = c(-1, 0.8)) +
  scale_fill_manual(values = rev(c('#A50104', '#7DAF9C'))) +
  annotate('text', label = c('Mixed Cluster 1', 'Mixed Cluster 2'), 
           x = c(-0.3, 0.3), y = c(0.15, -0.6),
           size = 15) +
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 35))

sd <- ggplot(s1c, aes(x = Axis.1, y = Axis.2)) + 
  #stat_ellipse(aes(group = cluster, fill = cluster), geom = 'polygon', color = 'black', alpha = 0.37, show.legend = FALSE) +
  geom_point(aes(fill = treatment), color = 'black', pch = 21, size = 15, show.legend = TRUE) +
  xlab('Axis 1 [69.2%]') +
  ylab('Axis 2 [22.3%]') +
  theme_cowplot() + 
  scale_x_continuous(limits = c(-0.5, 1)) +
  scale_y_continuous(limits = c(-1, 0.8)) +
  scale_fill_manual(values = c(trt[3], trt[1]), 
                    labels = c('NC', 'Broc'),
                    name = 'Treatment') +
  theme(legend.position = c(0.02,0.1),
        legend.text = element_text(size = 28),
        legend.title = element_text(size = 30),
        axis.text = element_text(size = 30),
        axis.title = element_text(size = 35))
sap <- align_plots(sd, sc, align = 'hv', axis = 'tbl')

ggdraw(sap[[2]]) + draw_plot(sap[[1]])

