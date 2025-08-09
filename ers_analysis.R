# created by Lingwei Ouyang for Neurohackademy project on August 8th, 2025
# analyze RSA results
getwd()
library(stringr)
library(dplyr)
library(ggplot2)
install.packages("viridis")
library(viridis)
library(lme4)
library(emmeans)
library(rstatix)
library(tidyr)
setwd('/Users/lingwei/Desktop/neurohackademy/')
list.files()

all_rsa <- read.csv('all_subjects_similarity_2.csv', stringsAsFactors = F)
head(all_rsa)
View(all_rsa)

length(unique(all_rsa$X73KID))

colnames(all_rsa)
# plot rsa, for retrieval 1 and retrieval 2
# for retrieval 1, separate by representation type
# for retrieval 2, separate by representation type
# separate by roi
# average over X73KID


#------data examination------
all_rsa
all_rsa_long <- all_rsa %>%
    pivot_longer(names_to = 'rsa_time', values_to = 'rsa', cols = c('enc_ret1', 'enc_ret2', 'ret1_ret2')) %>%
    pivot_longer(names_to = 'button_time', values_to = 'response', cols = c('enc_BUTTON', 'ret1_BUTTON', 'ret2_BUTTON')) %>%
    pivot_longer(names_to = 'rep_time', values_to = 'rep_type', cols = c('enc_presentation_type','ret1_presentation_type', 'ret2_presentation_type'))
all_rsa_long %>%
    group_by(SUBJECT) %>%
    summarise(n = length(unique(X73KID)))

#------RSA------
# averaged across stimuli
rsa_cols <- colnames(all_rsa)[10:12]
all_rsa %>%
    filter(ret2_presentation == 'Long-term' & ret1_presentation == 'Long-term') %>%
    group_by(SUBJECT, ROI) %>%
    summarise(n = n()) %>%
    filter(n < 300)
unique(all_rsa$ROI)
rois <- c('ERC', '35', '36', 'CA1', 'CA3', 'DG', 'PHC')

all_rsa <- all_rsa %>%
    mutate(ret1_presentation = ifelse(ret1_presentation_type == 'easy', 'Short-term', 'Long-term'),
           ret2_presentation = ifelse(ret2_presentation_type == 'easy', 'Short-term', 'Long-term'))

all_rsa$ret1_presentation <- factor(all_rsa$ret1_presentation, levels = c('Short-term', 'Long-term'))
rsa_sub <- all_rsa %>%
    filter(ROI %in% rois) %>%
    group_by(SUBJECT, ROI) %>%
    summarise_at(.funs = mean, .vars = rsa_cols)
View(rsa_sub)


#------plot------
dodge_width_param = 0.7
# plot first retrieval for short-term
rsa_sub <- all_rsa %>%
    filter(ROI %in% rois) %>%
    group_by(SUBJECT, ROI, ret1_presentation) %>%
    summarise_at(.funs = mean, .vars = 'enc_ret1')

range(rsa_sub$enc_ret1)

ggplot(rsa_sub %>% filter(ret1_presentation == 'Long-term'), aes(x = as.factor(ROI), y = enc_ret1, color = as.factor(ROI))) +
   # facet_wrap(~ret1_presentation) +
    geom_point(position = position_jitterdodge(dodge.width = dodge_width_param, jitter.width = 0.1), size = 5) + 
    geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 1) + 
 #   stat_summary(fun = mean, geom = "bar", linewidth = 0.6, position = position_dodge(dodge_width_param), width = 0.5) +
  #  stat_summary(geom = "pointrange", position = position_dodge(dodge_width_param), color = 'black') + 
  #  scale_fill_viridis() + 
    scale_color_brewer(palette = 'Paired') + 
 #   geom_errorbar(aes(ymax = emmean + SE, ymin = emmean - SE), position = position_dodge(dodge_width_param), width = 0.3) +
    ylab('RSA') + ggtitle(paste0('Long-term')) +
    xlab('ROI') +
    coord_cartesian(ylim = c(-0.01,0.13)) +
  #  scale_y_continuous(expand = c(0.01,0)) +
   theme_1 +
   facet_theme +
   theme(legend.position = 'none')
ggsave("long_term_by_roi.png", width = 12, height = 12, dpi = 300)


# plot second retrieval 
rsa_sub <- all_rsa %>%
    filter(ret1_presentation == 'Long-term' & ret2_presentation == 'Long-term') %>%
    filter(ROI %in% rois) %>%
    group_by(SUBJECT, ROI) %>%
    summarise_at(.funs = mean, .vars = 'enc_ret2')

ggplot(rsa_sub, aes(x = as.factor(ROI), y = enc_ret2, color = as.factor(ROI))) +
    geom_point(position = position_jitterdodge(dodge.width = dodge_width_param, jitter.width = 0.1), size = 5) + 
    geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 1) + 
    #   stat_summary(fun = mean, geom = "bar", linewidth = 0.6, position = position_dodge(dodge_width_param), width = 0.5) +
    #  stat_summary(geom = "pointrange", position = position_dodge(dodge_width_param), color = 'black') + 
    #  scale_fill_viridis() + 
    scale_color_brewer(palette = 'Paired') + 
    #   geom_errorbar(aes(ymax = emmean + SE, ymin = emmean - SE), position = position_dodge(dodge_width_param), width = 0.3) +
    ylab('RSA') + ggtitle(paste0('Remote')) +
    xlab('ROI') +
    #  coord_cartesian(ylim = c(-0.2,0.5)) +
    coord_cartesian(ylim = c(-0.01,0.13)) +
    theme_1 +
    facet_theme +
    theme(legend.position = 'none')
ggsave("remote_by_roi.png", width = 12, height = 12, dpi = 300)

rsa_sub %>%
    group_by(SUBJECT, ROI) %>%
    summarise(n = n())

#------logistic regression------
# response ~ RSA * time + (1/subject)
colnames(all_rsa)
all_rsa <- all_rsa %>%
    mutate(first_old = ifelse(ret1_BUTTON == 1, 0, ifelse(ret1_BUTTON == 2, 1, NaN)),
           second_old = ifelse(ret2_BUTTON == 1, 0, ifelse(ret2_BUTTON == 2, 1, NaN)))

# for short-term memory, P = 0.197, b = 0.49
fit.1.phc <- glmer(data = all_rsa %>% filter(ROI == 'PHC' & ret1_presentation == 'Short-term'), formula = first_old ~ enc_ret1 + (1|SUBJECT) + (1|X73KID), family = 'binomial')
check_convergence(fit.1.phc)
summary(fit.1.phc)

# for long-term memory, P = 0.15, b = 0.31
fit.2.phc <- glmer(data = all_rsa %>% filter(ROI == 'PHC' & ret1_presentation == 'Long-term'), formula = first_old ~ enc_ret1 + (1|SUBJECT) + (1|X73KID), family = 'binomial')
check_convergence(fit.2.phc)
summary(fit.2.phc)

# for remote memory, P = 0.30, b = 0.28
fit.3.phc <- glmer(data = all_rsa %>% filter(ROI == 'PHC' & ret1_presentation == 'Long-term' & ret2_presentation == 'Long-term'), 
                   formula = second_old ~ enc_ret2 + (1|SUBJECT) + (1|X73KID), family = 'binomial')
check_convergence(fit.3.phc)
summary(fit.3.phc)

# make plots, just for PHC and long term
all_rsa <- all_rsa %>%
    mutate(recog_first = ifelse(first_old == 1, 'Remember', ifelse(complete.cases(first_old), 'Forget', NaN)),
           recog_second = ifelse(second_old == 1, 'Remember', ifelse(complete.cases(second_old), 'Forget', NaN)))

# short-term
rsa_by_hit <- all_rsa %>%
    filter(ROI == 'PHC') %>%
    filter(ret1_presentation == 'Long-term') %>%
    group_by(SUBJECT, recog_first) %>%
    summarise(mean_first = mean(enc_ret1, na.rm = T))

rsa_by_hit$recog_first <- factor(rsa_by_hit$recog_first, levels = c('Remember', 'Forget'))
ggplot(rsa_by_hit %>% filter(complete.cases(recog_first)), aes(x = recog_first, y = mean_first, color = recog_first)) +
    geom_point(size = 5, position = position_jitter(0.07)) + 
    scale_color_manual(values=c("#CB181D",  "#084594")) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 1) + 
    #   geom_errorbar(aes(ymax = emmean + SE, ymin = emmean - SE), position = position_dodge(dodge_width_param), width = 0.3) +
    xlab('Response') +
    ylab('RSA') + ggtitle('Long-term') +
    coord_cartesian(ylim = c(-0.05,0.12)) +
    #   scale_y_continuous(expand = c(0.01,0)) +
    theme_1 +
    facet_theme +
    theme(legend.position = 'none')
ggsave("Remote_logistic.png", width = 12, height = 12, dpi = 300)




# long-term
rsa_by_hit <- all_rsa %>%
    filter(ROI == 'PHC') %>%
    filter(ret1_presentation == 'Long-term' & ret2_presentation == 'Long-term') %>%
    group_by(SUBJECT, recog_second) %>%
    summarise(mean_first = mean(enc_ret2, na.rm = T))

rsa_by_hit$recog_second <- factor(rsa_by_hit$recog_second, levels = c('Remember', 'Forget'))
ggplot(rsa_by_hit %>% filter(complete.cases(recog_second)), aes(x = recog_second, y = mean_first, color = recog_second)) +
    geom_point(size = 5, position = position_jitter(0.07)) + 
    scale_color_manual(values=c("#CB181D",  "#084594")) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 1) + 
    #   geom_errorbar(aes(ymax = emmean + SE, ymin = emmean - SE), position = position_dodge(dodge_width_param), width = 0.3) +
    xlab('Response') +
    ylab('RSA') + ggtitle('Remote') +
    coord_cartesian(ylim = c(-0.05,0.12)) +
 #   scale_y_continuous(expand = c(0.01,0)) +
    theme_1 +
    facet_theme +
    theme(legend.position = 'none')
ggsave("Remote_logistic.png", width = 12, height = 12, dpi = 300)
