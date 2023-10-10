library(tidyverse)

#### load data ####
load('data.RData') 

#### Detectability ####
detectability <- npxData_simulatedClinicalData %>% 
  group_by(Assay, Panel, OlinkID, UniProt) %>% 
  summarise(percDetect = 100*sum(NPX>LOD)/n(), .groups = 'drop')

#plot
detectability %>% 
  ggplot(aes(x = percDetect)) + 
  geom_histogram() +
  facet_wrap(~Panel, nrow = 2, ncol = 4) +
  xlab('Percent samples with NPX > LOD') +
  ylab('Number of proteins') +
  theme_bw() +
  theme(strip.text = element_text(size = 15), 
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 13))

#### NPX ranges ####
npxRanges <- 
  npxData_simulatedClinicalData %>%
  group_by(Panel, OlinkID, UniProt, Assay) %>%
  summarise(median = median(NPX, na.rm = T), 
            perc10 = quantile(NPX, na.rm = T, probs = seq(0,1,.1))[2],
            perc90 = quantile(NPX, na.rm = T, probs = seq(0,1,.1))[10]) %>% 
  ungroup() %>% 
  mutate(range = perc90 - perc10) %>% 
  arrange(desc(range))

#plot
npxRanges %>%
  ggplot(aes(x = range)) + 
  geom_histogram() +
  facet_wrap(~Panel, nrow = 2, ncol = 4) +
  xlab('NPX range') +
  ylab('Number of proteins') +
  theme_bw() +
  theme(strip.text = element_text(size = 15), 
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 13)) +
  scale_x_continuous(labels = function(x){signif(x, digits = 2)}, limits = c(0,10.3)) 


#### NPX levels versus baseline characteristics ####
baselineTraits_lm <- npxData_simulatedClinicalData %>% 
  mutate(menopause_status = case_when(menopause_status == 1 ~ '.pre',
                                      menopause_status == 2 ~ '.peri',
                                      menopause_status == 3 ~ '.post') %>% factor(., levels = c('.pre', '.peri', '.post'))) %>% 
  group_by(OlinkID, Assay, UniProt) %>% 
  group_modify(.f = function(x, y){broom::tidy(lm(NPX ~ age + bmi + menopause_status + birth_times + HRT + alcohol_gram_week + smoking, 
                                                  data = x))}) %>% 
  ungroup() %>% 
  filter(term != '(Intercept)') %>% 
  mutate(p.value_adjusted = p.adjust(p.value, method = 'fdr'))

#volcano plot
protsToLabel <- 
  baselineTraits_lm %>%
  filter(p.value_adjusted < .05) %>% 
  group_by(term) %>% 
  slice_min(order_by = p.value_adjusted, n = 10) %>% 
  ungroup() %>% 
  rename(label = Assay) %>% 
  select(OlinkID, term, label)

baselineTraits_lm %>% 
  left_join(protsToLabel, by = c('OlinkID', 'term')) %>% 
  mutate(trait = case_when(term == 'alcohol_gram_week' ~ 'alcohol (gram/week)',
                           term == 'HRT.current' ~ 'hormone replacement therapy',
                           term == 'menopause_status.peri' ~ 'menopause (pre vs peri)',
                           term == 'menopause_status.post' ~ 'menopause (pre vs post)',
                           term == 'smoking.current' ~ 'smoking',
                           term == 'birth_times' ~ 'birth times',
                           TRUE ~ term)) %>% 
  ggplot(aes(x = estimate, y = -log10(p.value), label = label)) +
  geom_point() +
  xlab('Effect Size') +
  ylab('-log10(p-value)') +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 16)) +
  facet_wrap(~trait, scales = 'free', nrow = 2, ncol = 4) +
  ggrepel::geom_label_repel()
