## qPCR resluts
########## Library########## 
library(phyloseq)
library(tidyverse)
library(tidyr)
library(data.table)
require(ggpubr)
require(rstatix)
library(wesanderson)
library(cowplot)
library(ggsci)
library(rtk)
library(factoextra)

ana.red.comp <- read.table(file = "output/ana.red.metadecon.tsv", sep = "\t", header = T, row.names = 1)

vag.red.comp <- read.table(file = "output/vag.red.metadecon.tsv", sep = "\t", header = T, row.names = 1)

qpcr <- read.table(file = "../Figures/march.fig/input/qpcr_final.csv", header = T, sep =",", row.names = 1)

# Treatment == 1, Control == 0
# complinance = baseline/non vs. good 2/3
# 
qpcr$treatment[qpcr$treatment == "1"] <- "treatment"
qpcr$treatment[qpcr$treatment == "0"] <- "control"

qpcr$compliance[qpcr$compliance == "Na"] <- NA

qpcr$timepoint[qpcr$timepoint == "1"] <- "baseline"
qpcr$timepoint[qpcr$timepoint == "2"] <- "timepoint 2"
qpcr$timepoint[qpcr$timepoint == "3"] <- "timepoint 3"
qpcr <- qpcr %>% 
  filter(!compliance == "NA")

## Rank values but keep NAs percentile rank function in R
prank<-function(x){
  r<-rank(x)/sum(!is.na(x))
  r[is.na(x)]<-NA
  r
}
## normal rank without percentile
prank<-function(x) ifelse(is.na(x),NA,rank(x))
for (i in 6:11) {
  qpcr[, i] <- prank(qpcr[, i])
} 

# qpcr_vag <- qpcr[qpcr$swabsite == "Vag_Swabs", ]
# write.table(qpcr_vag, file = "../Figures/march.fig/helper/perecent.cal.vag.only.pcr.csv", sep = ",", quote = F)
# qpcr_ana <- qpcr[qpcr$swabsite == "Ana_Swabs", ]
# write.table(qpcr_ana, file = "../Figures/march.fig/helper/perecent.cal.rec.only.pcr.csv", sep = ",", quote = F)
# write.table(qpcr_ana, file = "../Figures/march.fig/helper/perecent.cal.rec.only.pcr.csv", sep = ",", quote = F)

qpcr_vag <- read.table(file = "../Figures/march.fig/helper/perecent.cal.vag.only.pcr.csv", sep = ",", header =T)
qpcr_rec <- read.table(file = "../Figures/march.fig/helper/perecent.cal.rec.only.pcr.csv", sep = ",", header =T)



### tabel 
qpcr%>% 
  filter(!swabsite == "Ana_Swabs") %>% 
  filter(timepoint == "baseline") %>% 
  filter(compliance == "control") %>% 
  dplyr::group_by(timepoint) -> baseline.vag

qpcr%>% 
  filter(!swabsite == "Ana_Swabs") %>% 
  filter(timepoint == "timepoint 2") %>% 
  filter(compliance == "control") %>% 
  dplyr::group_by(timepoint) -> timep2.vag

qpcr%>% 
  filter(!swabsite == "Ana_Swabs") %>% 
  filter(timepoint == "timepoint 3") %>% 
  filter(compliance == "control") %>% 
  dplyr::group_by(timepoint) -> timep3.vag

qpcr%>% 
  filter(!swabsite == "Ana_Swabs") %>% 
  filter(timepoint == "baseline") %>% 
  filter(treatment == "treatment") %>% 
  filter(compliance == "control") %>% 
  dplyr::group_by(timepoint) -> baseline.vag

qpcr%>% 
  filter(!swabsite == "Ana_Swabs") %>% 
  filter(timepoint == "timepoint 2") %>% 
  filter(treatment == "treatment") %>% 
  filter(compliance == "good") %>%
  dplyr::group_by(timepoint) -> timep2.vag

qpcr%>% 
  filter(!swabsite == "Ana_Swabs") %>% 
  filter(timepoint == "timepoint 3") %>% 
  filter(treatment == "treatment") %>% 
  filter(compliance == "good") %>%
  dplyr::group_by(timepoint) -> timep3.vag


library(rstatix)
library(ggpubr)


##Color palette for collections ##


pal.treat <- c("#803300ff", "#822bc1ff")


########## loop vaginal control ########## 
all.red.na <- qpcr
all.red.na$ID <- as.factor(all.red.na$ID)


all.red.na <- all.red.na[, -c(1,3)]
all.red.na %>% 
  reshape2::melt(id.vars=c("ID", "timepoint", "compliance", "swabsite")) -> all.red.na

all.red.na %>% 
  filter(variable == "total_bacterialload") %>% 
  ggplot(aes(x = timepoint, y = value, color = swabsite, fill = swabsite)) +
  geom_boxplot(color = "black") +
  # geom_point(aes(fill = timepoint, group=ID),size=3,shape=21, position = position_dodge(0.2))+
  facet_wrap(~compliance, ncol = 2) +
  # geom_line(aes(group=ID)) +
  theme_bw() +
  ylab("rank (copynumber)")+
  scale_fill_manual(values = alpha(pal.treat, 0.8)) +
  scale_colour_manual(values = pal.treat) +
  theme (axis.text.x = element_text(angle = 90, 
                                    vjust = 0.5, 
                                    hjust = 1,size = 16), 
         axis.title.x = element_blank(),
         axis.text.y = element_text(size = 16),
         axis.title.y = element_text ( size = 16),
         legend.position = "none") -> qpcr.total



# stats
qpcr %>% 
  filter(compliance == "good") %>% 
  rstatix::wilcox_test(total_bacterialload~ swabsite)%>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance()
# control/good - Vag
# no between total_bacterialload

# control/good - Anal
# no between total_bacterialload

# anal - rectal 
# 1.26e-23


ggsave(qpcr.total, file = "../Figures/march.fig/figure/qpcr.total.svg", width = 5, height = 4, device = "svg")




########## loop vaginal good ########## 
all.red.na <- qpcr[qpcr$swabsite == "Vag_Swabs", ]

# Dear Theda (CC) - could you test mixed effect models for the qPCR time series alone, for these species alone, for non-baseline time points?
# I would propose: m = ranked_copy_number ~ (1 | individual) + ranked_copy_number_baseline + time_point + compliance

all.red.na <- all.red.na[, -c(1:3)]
all.red.na.basline <- all.red.na %>% 
  filter(timepoint == "baseline")

colnames(all.red.na.basline) <- paste0("baseline_",colnames(all.red.na.basline))
all.red.na.basline[all.red.na.basline$baseline_ID]<- "ID"
colnames(all.red.na.basline)[9] <- "ID"

all.red.na.basline %>% 
reshape2::melt(id.vars=c("ID", "baseline_timepoint", "baseline_compliance")) -> all.red.na.basline
all.red.na %>% 
  reshape2::melt(id.vars=c("ID", "timepoint", "compliance")) -> all.red.na



all.red.na <- all.red.na %>% 
  filter(!variable == "total_bacterialload")
all.red.na.basline <- all.red.na.basline %>% 
  filter(!variable == "baseline_total_bacterialload")

write.table(all.red.na, file = "../Figures/march.fig/helper/17423.helper.lmer.qpcr.sofia.csv", sep = ",", quote = F)
write.table(all.red.na.basline, file = "../Figures/march.fig/helper/17423.helper.lmer.qpcr.sofia.baseline.csv", sep = ",", quote = F)

new <- read.table(file = "../Figures/march.fig/helper/17423.helper.lmer.qpcr.sofia.csv", sep = ",", header = T, row.names = 1)
new <- new %>% 
  filter(!timepoint == "baseline")


library(lmtest)
library(lme4)
new$timepoint[new$timepoint == "timepoint 2"] <- 0
new$timepoint[new$timepoint == "timepoint 3"] <- 1
new$timepoint <- as.numeric(new$timepoint)

new$compliance[new$compliance == "control"] <- 0
new$compliance[new$compliance == "good"] <- 1
new$compliance <- as.numeric(new$compliance)
new <- new %>% 
  filter(variable == "L.Iners")

llr.test <- function(x){
  l1 <- lmer(value ~ (1 | ID)+ value_baseline + timepoint + compliance, REML = F, data = x)
  l2 <- lmer(value ~ (1 | ID)+ value_baseline + timepoint, REML = F, data = x)
  return(lrtest(l1,l2)$'Pr(>Chisq)'[[2]])}

llr.test(new) #  is significant, it would indicate higher impact in the probiotic group.
# L.crispatus  0.9084932
# L.Gasseri  0.1210994
# L.Jensenii 0.3047795
# L.Rhamnosus no enough observations
# L.Iners 0.5315707

##### do the same for the anal ########
all.red.na <- qpcr[qpcr$swabsite == "Ana_Swabs", ]

# Dear Theda (CC) - could you test mixed effect models for the qPCR time series alone, for these species alone, for non-baseline time points?
# I would propose: m = ranked_copy_number ~ (1 | individual) + ranked_copy_number_baseline + time_point + compliance

all.red.na <- all.red.na[, -c(1:3)]
all.red.na.basline <- all.red.na %>% 
  filter(timepoint == "baseline")

colnames(all.red.na.basline) <- paste0("baseline_",colnames(all.red.na.basline))
colnames(all.red.na.basline)[9] <- "ID"

all.red.na.basline %>% 
  reshape2::melt(id.vars=c("ID", "baseline_timepoint", "baseline_compliance")) -> all.red.na.basline
all.red.na %>% 
  reshape2::melt(id.vars=c("ID", "timepoint", "compliance")) -> all.red.na



all.red.na <- all.red.na %>% 
  filter(!variable == "total_bacterialload")
all.red.na.basline <- all.red.na.basline %>% 
  filter(!variable == "baseline_total_bacterialload")

write.table(all.red.na, file = "../Figures/march.fig/helper/17423.anal.helper.lmer.qpcr.sofia.csv", sep = ",", quote = F)
write.table(all.red.na.basline, file = "../Figures/march.fig/helper/17423.anal.helper.lmer.qpcr.sofia.baseline.csv", sep = ",", quote = F)

new <- read.table(file = "../Figures/march.fig/helper/17423.anal.helper.lmer.qpcr.sofia.csv", sep = ",", header = T, row.names = 1)
new <- new %>% 
  filter(!timepoint == "baseline")


library(lmtest)
library(lme4)
new$timepoint[new$timepoint == "timepoint 2"] <- 0
new$timepoint[new$timepoint == "timepoint 3"] <- 1
new$timepoint <- as.numeric(new$timepoint)

new$compliance[new$compliance == "control"] <- 0
new$compliance[new$compliance == "good"] <- 1
new$compliance <- as.numeric(new$compliance)
new <- new %>% 
  filter(variable == "L.Iners")

llr.test <- function(x){
  l1 <- lmer(value ~ (1 | ID)+ value_baseline + timepoint + compliance, REML = F, data = x)
  l2 <- lmer(value ~ (1 | ID)+ value_baseline + timepoint, REML = F, data = x)
  return(lrtest(l1,l2)$'Pr(>Chisq)'[[2]])}

llr.test(new) #  is significant, it would indicate higher impact in the probiotic group.
# L.crispatus   0.4269775
# L.Gasseri  0.3841625
# L.Jensenii 0.3637829
# L.Rhamnosus no enough observations 6.342111e-52
# L.Iners  0.96843






all.red.na <- all.red.na %>% 
  filter(compliance == "good")

all.red.na <- all.red.na[, -c(1:4)]
all.red.na %>% 
  reshape2::melt(id.vars=c("ID", "timepoint")) -> all.red.na
all.red.na <- all.red.na %>% 
filter(!variable == "total_bacterialload")


pal.treat <- c("#822bc1ff","#c094dfff", "#ec92aeff")
all.red.na %>% 
  ggplot(aes(x = timepoint, y = value, color = timepoint, fill = timepoint)) +
  geom_boxplot(aes(fill = timepoint, alpha = 0.6),color = "black") +
  geom_point(aes(fill = timepoint, group=ID),size=3,shape=21, position = position_dodge(0.2), color = "black")+
  facet_wrap(~variable, ncol = 3) +
  geom_line(aes(group=ID)) +
  theme_bw() +
  ylab("rank (copynumber)")+
  scale_fill_manual(values = pal.treat) +
  scale_colour_manual(values = pal.treat) +
  theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size = 16), axis.title.x = element_blank(),
         axis.text.y = element_text(size = 16),
         axis.title.y = element_text(size = 16),
         legend.position = "none")  -> qpcr.vag.treatment



qpcr.vagina <- cowplot::plot_grid(qpcr.vag.control, qpcr.vag.treatment, align = "hv")


ggsave(qpcr.vagina, file = "../Figures/march.fig/figure/qpcr.vagina.svg", device = "svg", width = 14, height = 9)

all.red.na <- qpcr[qpcr$swabsite == "Vag_Swabs", ]
all.red.na <- all.red.na %>% 
  filter(compliance == "control")

all.red.na%>% 
  rstatix::wilcox_test(L.Rhamnosus~ timepoint)%>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance()
# good/control
# 16S total no/no
# L.Crispatus no/no
# L.Gasseri no/no
# L.Jensenii no/no
# L.Iners no/no
# L.Rhamnosus no/no





####### loop anal control ########## 
all.red.na <- qpcr[qpcr$swabsite == "Ana_Swabs", ]
all.red.na <- all.red.na %>% 
  filter(compliance == "control")

all.red.na <- all.red.na[, -c(1:4)]
all.red.na %>% 
  reshape2::melt(id.vars=c("ID", "timepoint")) -> all.red.na
all.red.na <- all.red.na %>% 
  filter(!variable == "total_bacterialload")


pal.treat <- c("#803300ff","#c97a20ff", "#e4bc8dff")
all.red.na %>% 
  ggplot(aes(x = timepoint, y = value, color = timepoint, fill = timepoint)) +
  geom_boxplot(aes(fill = timepoint, alpha = 0.6),color = "black") +
  geom_point(aes(fill = timepoint, group=ID),size=3,shape=21, position = position_dodge(0.2), color = "black")+
  facet_wrap(~variable, ncol = 3) +
  geom_line(aes(group=ID)) +
  theme_bw() +
  ylab("rank (copynumber)")+
  scale_fill_manual(values = pal.treat) +
  scale_colour_manual(values = pal.treat) +
  theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size = 16), axis.title.x = element_blank(),
         axis.text.y = element_text(size = 16),
         axis.title.y = element_text(size = 16),
         legend.position = "none")  -> qpcr.vag.control

qpcr.ana <- cowplot::plot_grid(qpcr.vag.control, qpcr.vag.treatment, align = "hv")


ggsave(qpcr.ana, file = "../Figures/march.fig/figure/qpcr.ana.svg", device = "svg", width = 14, height = 9)

all.red.na <- qpcr[qpcr$swabsite == "Ana_Swabs", ]
all.red.na <- all.red.na %>% 
  filter(compliance == "good")

all.red.na%>% 
  rstatix::wilcox_test(L.Rhamnosus~ timepoint)%>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance()
# control
# 16S total no
# L.Crispatus no/no
# L.Gasseri no/no
# L.Jensenii no/no
# L.Iners no/no
# L.Rhamnosus no/no

all.red.c%>% 
  rstatix::wilcox_test(total_bacterialload~ timepoint)%>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance()%>%
  rstatix::add_xy_position(x = "timepoint") -> stats_test

# control   tota… Ana_S… Vag_S…    37    33      1221 2e-20 2e-20 ****            9.82e10

qpcr %>% 
  filter(compliance == "control") %>%  
  ggplot(aes(x =  swabsite, #descending order
             y = log(total_bacterialload), fill = timepoint)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  labs(caption = get_pwc_label(stats_test))+
  # geom_point(aes(y = diversity_shannon, color = origin), 
  # size = 2, alpha = 1) +
  labs(y = "log(total 16S copy)", x = NULL) + #\n adds a new line which creates some space between the axis and axis title
  guides(color = FALSE) + #remove legends
  scale_x_discrete(labels = c("anal", "vaginal")) +
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal) +
  theme_classic()  +
  ggtitle("control")+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=14)) -> bact.contr
  # stat_pvalue_manual(stats_test, label = "{p.adj} {p.adj.signif}")-> bact.contr

grobs <- ggplotGrob(bact.contr)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]


qpcr%>% 
  filter(treatment == "treatment") %>% 
  dplyr::group_by(treatment)%>%
  rstatix::wilcox_test(total_bacterialload ~ swabsite)%>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance()%>%
  rstatix::add_xy_position(x = "timepoint") -> stats_test
# treatment total_bac… Ana_S… Vag_S…    57    57      3249 3.46e-20 3.46e-20 ****  

qpcr %>% 
  filter(compliance == "good") %>%  
  ggplot(aes(x =  swabsite, #descending order
             y = log(total_bacterialload), fill = timepoint)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  labs(caption = get_pwc_label(stats_test))+
  # geom_point(aes(y = diversity_shannon, color = origin), 
  # size = 2, alpha = 1) +
  labs(y = "log(total 16S copy)", x = NULL) + #\n adds a new line which creates some space between the axis and axis title
  guides(color = FALSE) + #remove legends
  scale_x_discrete(labels = c("anal", "vaginal")) +
  scale_fill_manual(values = pal.good) +
  scale_colour_manual(values = pal.good) +
  theme_classic()  +
  ggtitle("treatment")+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=14)) -> bact.treat

bact <- cowplot::plot_grid(bact.contr, bact.treat, rel_widths = c(1, 1))
bact <- cowplot::plot_grid(bact, legend, rel_widths = c(1, .2))


ggsave(bact, width = 8, height = 4, dpi = 300, device = "svg", filename = "../Figures/2023.figure/2322.bacload.compliance. vag.ana.svg")




####### vagina ###########

###### L.Crispatus ########
qpcr_vag %>% 
  filter(!compliance == "NA") %>% 
  # filter(treatment == "control") %>% 
  dplyr::group_by(compliance)%>%
  rstatix::wilcox_test(L.Crispatus ~ timepoint)%>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance()%>%
  rstatix::add_xy_position(x = "timepoint") -> stats_test

# ns both control and treatment over time/ same for compliance 


qpcr %>%
  filter(swabsite == "Vag_Swabs") %>%
  filter(!compliance == "NA") %>% 
  filter(compliance == "control") %>%  
  ggplot(aes(x =  timepoint, #descending order
             y = log(L.Crispatus), fill = timepoint)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  labs(caption = get_pwc_label(stats_test))+
  # geom_point(aes(y = diversity_shannon, color = origin), 
  # size = 2, alpha = 1) +
  labs(y = "log(L.Crispatus copy)", x = NULL) + #\n adds a new line which creates some space between the axis and axis title
  guides(color = FALSE) + #remove legends
  scale_fill_manual(values = pal.treat) +
  scale_colour_manual(values = pal.treat) +
  theme_classic()  +
  ggtitle("control")+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=14)) -> bact.contr
# stat_pvalue_manual(stats_test, label = "{p.adj} {p.adj.signif}")-> bact.contr

grobs <- ggplotGrob(bact.contr)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

qpcr %>% 
  filter(swabsite == "Vag_Swabs") %>%
  filter(!compliance == "NA") %>% 
  filter(compliance == "good") %>%  
  ggplot(aes(x =  timepoint, #descending order
             y = log(L.Crispatus), fill = timepoint)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  labs(caption = get_pwc_label(stats_test))+
  # geom_point(aes(y = diversity_shannon, color = origin), 
  # size = 2, alpha = 1) +
  labs(y = "log(L.Crispatus copy)", x = NULL) + #\n adds a new line which creates some space between the axis and axis title
  guides(color = FALSE) + #remove legends
  scale_fill_manual(values = pal.vag) +
  scale_colour_manual(values = pal.vag) +
  theme_classic()  +
  ggtitle("treatment")+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=14)) -> bact.treat


###### L.Gasseri ########

qpcr_vag %>% 
  filter(!compliance == "NA") %>% 
  # filter(treatment == "control") %>% 
  dplyr::group_by(timepoint)%>%
  rstatix::wilcox_test(L.Gasseri ~ compliance)%>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance()%>%
  rstatix::add_xy_position(x = "timepoint") -> stats_test

# ns both control and treatment over time/ same for compliance 


qpcr %>%
  filter(swabsite == "Vag_Swabs") %>%
  filter(!compliance == "NA") %>% 
  filter(compliance == "control") %>%  
  ggplot(aes(x =  timepoint, #descending order
             y = log(L.Gasseri), fill = timepoint)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  labs(caption = get_pwc_label(stats_test))+
  # geom_point(aes(y = diversity_shannon, color = origin), 
  # size = 2, alpha = 1) +
  labs(y = "log(L.Crispatus copy)", x = NULL) + #\n adds a new line which creates some space between the axis and axis title
  guides(color = FALSE) + #remove legends
  scale_fill_manual(values = pal.vag) +
  scale_colour_manual(values = pal.vag) +
  theme_classic()  +
  ggtitle("control")+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=14)) -> grass.contr
# stat_pvalue_manual(stats_test, label = "{p.adj} {p.adj.signif}")-> bact.contr

grobs <- ggplotGrob(grass.contr)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

qpcr %>% 
  filter(swabsite == "Vag_Swabs") %>%
  filter(!compliance == "NA") %>% 
  filter(compliance == "good") -> helper
  
helper %>% 
  ggplot(aes(x =  timepoint, #descending order
             y = log(L.Gasseri), fill = timepoint)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  labs(caption = get_pwc_label(stats_test))+
  # geom_point(aes(y = diversity_shannon, color = origin), 
  # size = 2, alpha = 1) +
  labs(y = "log(L.Crispatus copy)", x = NULL) + #\n adds a new line which creates some space between the axis and axis title
  guides(color = FALSE) + #remove legends
  scale_fill_manual(values = pal.vag) +
  scale_colour_manual(values = pal.vag) +
  theme_classic()  +
  ggtitle("treatment")+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=14)) -> grassi.treat



###### L.Jensenii ########

qpcr_vag %>% 
  filter(!compliance == "NA") %>% 
  # filter(treatment == "control") %>% 
  dplyr::group_by(timepoint)%>%
  rstatix::wilcox_test(L.Jensenii ~ compliance)%>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance()%>%
  rstatix::add_xy_position(x = "timepoint") -> stats_test

# ns both control and treatment over time/ same for compliance 


qpcr %>%
  filter(swabsite == "Vag_Swabs") %>%
  filter(!compliance == "NA") %>% 
  filter(compliance == "control") %>%  
  ggplot(aes(x =  timepoint, #descending order
             y = log(L.Jensenii), fill = timepoint)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  labs(caption = get_pwc_label(stats_test))+
  # geom_point(aes(y = diversity_shannon, color = origin), 
  # size = 2, alpha = 1) +
  labs(y = "log(L.Crispatus copy)", x = NULL) + #\n adds a new line which creates some space between the axis and axis title
  guides(color = FALSE) + #remove legends
  scale_fill_manual(values = pal.vag) +
  scale_colour_manual(values = pal.vag) +
  theme_classic()  +
  ggtitle("control")+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=14)) -> jensi.contr
# stat_pvalue_manual(stats_test, label = "{p.adj} {p.adj.signif}")-> bact.contr


qpcr %>% 
  filter(swabsite == "Vag_Swabs") %>%
  filter(!compliance == "NA") %>% 
  filter(compliance == "good") -> helper

helper %>% 
  ggplot(aes(x =  timepoint, #descending order
             y = log(L.Jensenii), fill = timepoint)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  labs(caption = get_pwc_label(stats_test))+
  # geom_point(aes(y = diversity_shannon, color = origin), 
  # size = 2, alpha = 1) +
  labs(y = "log(L.Crispatus copy)", x = NULL) + #\n adds a new line which creates some space between the axis and axis title
  guides(color = FALSE) + #remove legends
  scale_fill_manual(values = pal.vag) +
  scale_colour_manual(values = pal.vag) +
  theme_classic()  +
  ggtitle("treatment")+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=14)) -> jensi.treat


###### L.Jensenii ########

qpcr %>%
  filter(swabsite == "Vag_Swabs") %>%
  filter(!compliance == "NA") %>%  
  rstatix::wilcox_test(L.Rhamnosus ~ timepoint)%>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance()%>%
  rstatix::add_xy_position(x = "timepoint") 


# ns both control and treatment over time/ same for compliance 


qpcr %>%
  filter(swabsite == "Vag_Swabs") %>%
  filter(!compliance == "NA") %>% 
  filter(compliance == "control") %>%  
  ggplot(aes(x =  timepoint, #descending order
             y = log(L.Rhamnosus), fill = timepoint)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  labs(caption = get_pwc_label(stats_test))+
  # geom_point(aes(y = diversity_shannon, color = origin), 
  # size = 2, alpha = 1) +
  labs(y = "log(L.Crispatus copy)", x = NULL) + #\n adds a new line which creates some space between the axis and axis title
  guides(color = FALSE) + #remove legends
  scale_fill_manual(values = pal.vag) +
  scale_colour_manual(values = pal.vag) +
  theme_classic()  +
  ggtitle("control")+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=14)) -> rham.contr
# stat_pvalue_manual(stats_test, label = "{p.adj} {p.adj.signif}")-> bact.contr



helper %>% 
  ggplot(aes(x =  timepoint, #descending order
             y = log(L.Rhamnosus), fill = timepoint)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  labs(caption = get_pwc_label(stats_test))+
  # geom_point(aes(y = diversity_shannon, color = origin), 
  # size = 2, alpha = 1) +
  labs(y = "log(L.Crispatus copy)", x = NULL) + #\n adds a new line which creates some space between the axis and axis title
  guides(color = FALSE) + #remove legends
  scale_fill_manual(values = pal.vag) +
  scale_colour_manual(values = pal.vag) +
  theme_classic()  +
  ggtitle("treatment")+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=14)) -> rham.treat



######## L.Iners ########

qpcr %>%
  filter(swabsite == "Vag_Swabs") %>%
  filter(!compliance == "NA") %>%  
  rstatix::wilcox_test(L.Iners ~ treatment)%>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance()%>%
  rstatix::add_xy_position(x = "timepoint") 

# ns both control and treatment over time/ same for compliance 


qpcr %>%
  filter(swabsite == "Vag_Swabs") %>%
  filter(!compliance == "NA") %>% 
  filter(compliance == "control") %>%  
  ggplot(aes(x =  timepoint, #descending order
             y = log(L.Iners), fill = timepoint)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  labs(caption = get_pwc_label(stats_test))+
  # geom_point(aes(y = diversity_shannon, color = origin), 
  # size = 2, alpha = 1) +
  labs(y = "log(L.Crispatus copy)", x = NULL) + #\n adds a new line which creates some space between the axis and axis title
  guides(color = FALSE) + #remove legends
  scale_fill_manual(values = pal.vag) +
  scale_colour_manual(values = pal.vag) +
  theme_classic()  +
  ggtitle("control")+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=14)) -> iner.contr
# stat_pvalue_manual(stats_test, label = "{p.adj} {p.adj.signif}")-> bact.contr



helper %>% 
  ggplot(aes(x =  timepoint, #descending order
             y = log(L.Iners), fill = timepoint)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  labs(caption = get_pwc_label(stats_test))+
  # geom_point(aes(y = diversity_shannon, color = origin), 
  # size = 2, alpha = 1) +
  labs(y = "log(L.Crispatus copy)", x = NULL) + #\n adds a new line which creates some space between the axis and axis title
  guides(color = FALSE) + #remove legends
  scale_fill_manual(values = pal.vag) +
  scale_colour_manual(values = pal.vag) +
  theme_classic()  +
  ggtitle("treatment")+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=14)) -> iner.treat







bact <- cowplot::plot_grid(bact.contr, bact.treat, rel_widths = c(1, 1))
bact <- cowplot::plot_grid(bact, legend, rel_widths = c(1, .2))


ggsave(bact, width = 8, height = 4, dpi = 300, device = "svg", filename = "../Figures/2023.figure/2322.bacload.vag.ana.svg")





######## anal ########
### tabel 
qpcr%>% 
  filter(swabsite == "Ana_Swabs") %>% 
  filter(timepoint == "baseline") %>% 
  filter(treatment == "control") %>% 
  dplyr::group_by(timepoint) -> baseline.ana

qpcr%>% 
  filter(swabsite == "Ana_Swabs") %>% 
  filter(timepoint == "timepoint 2") %>% 
  filter(treatment == "control") %>% 
  dplyr::group_by(timepoint) -> timep2.ana

qpcr%>% 
  filter(swabsite == "Ana_Swabs") %>% 
  filter(timepoint == "timepoint 3") %>% 
  filter(treatment == "control") %>% 
  dplyr::group_by(timepoint) -> timep3.ana

qpcr%>% 
  filter(swabsite == "Ana_Swabs") %>% 
  filter(timepoint == "baseline") %>% 
  filter(treatment == "treatment") %>% 
  filter(compliance == "control") %>% 
  dplyr::group_by(timepoint) -> baseline.ana

qpcr%>% 
  filter(swabsite == "Ana_Swabs") %>% 
  filter(treatment == "treatment") %>% 
  filter(timepoint == "timepoint 2") %>% 
  filter(!compliance == "NA") %>% 
  dplyr::group_by(timepoint) -> timep2.ana

qpcr%>% 
  filter(swabsite == "Ana_Swabs") %>% 
  filter(treatment == "treatment") %>% 
  filter(timepoint == "timepoint 3") %>% 
  filter(!compliance == "NA") %>% 
  dplyr::group_by(timepoint) -> timep3.ana


################### combination with clustertypes vagina 

qpcr <- read.table(file = "../Figures/march.fig/helper/qpcr.ranked.helper.csv", sep = ",", header = T, row.names = 1)
# write.table(otu.gen, file = "../Figures/march.fig/helper/otu.gen.need.id.tsv", sep = "\t", quote = F)
qpcr$ID <- row.names(qpcr)
qpcr$Row.names <- NULL
qpcr$Sampl_ID <- NULL


library(DirichletMultinomial)
load(file = "output/vag.dmn_list.tsv")
Dirichlet_multinomial_1 = mixture(dmn_list[[1]], assign = TRUE)
Dirichlet_multinomial_2 = mixture(dmn_list[[2]], assign = TRUE)
Dirichlet_multinomial_3 = mixture(dmn_list[[3]], assign = TRUE)
Dirichlet_multinomial_4 = mixture(dmn_list[[4]], assign = TRUE)
Dirichlet_multinomial_5 = mixture(dmn_list[[5]], assign = TRUE)
Dirichlet_multinomial_6 = mixture(dmn_list[[6]], assign = TRUE)

Dirichlet_multinomial_all = data.frame(cbind(Dirichlet_multinomial_1,
                                             Dirichlet_multinomial_2,Dirichlet_multinomial_3,
                                             Dirichlet_multinomial_4,Dirichlet_multinomial_5,
                                             Dirichlet_multinomial_6))
colnames(Dirichlet_multinomial_all) = c("DMM_k=1","DMM_k=2","DMM_k=3",
                                        "DMM_k=4","DMM_k=5","DMM_k=6")
cluster_result <- as.data.frame(Dirichlet_multinomial_3)
cluster_result$ID <- row.names(cluster_result)
# qpcr_red <- read.table(file = "qpcr.helper.tsv", sep ="\t", header = T, row.names = 1)



########### vagina #########
qpcr_vag <- qpcr %>% 
  filter(!swabsite == "Ana_Swabs")

comp <- left_join(qpcr_vag, cluster_result, by = "ID")
comp <- comp %>% 
  filter(!compliance == "NA")


### sum 

# MEDIAN
comp_1 <-comp[comp$Dirichlet_multinomial_3 == 1, ]

comp_2 <-comp[comp$Dirichlet_multinomial_3 == 2, ]

comp_3 <-comp[comp$Dirichlet_multinomial_3 == 3, ]

# merging
cluster_merged <-
  rbind(comp_1,
        comp_2,
        comp_3)

write.table(cluster_merged, file = "../Figures/march.fig/helper/cluster.merge.vagina.tsv", sep = "\t", quote = F)




## Plot the Abundance 
cluster_rel <- cluster_merged
cluster_rel <- as.data.frame(t(cluster_rel))
cluster_rel <- apply(cluster_merged, 1, function(x) x/sum(x)) # rows 1 als prozente

## use cluster_imp 
clr <- cluster_rel
clr_import <- merge(clr, cluster_imp, by = 0)
row.names(clr_import) <- clr_import$Row.names
clr_import$Row.names <- NULL






######## treatment #############
comp_treat <- comp %>% 
  filter(!treatment == "control")

comp_treat <- comp_treat[, c("timepoint", "L.Crispatus", "L.Gasseri",
                             "L.Jensenii", "L.Rhamnosus", "L.Iners",
                             "Dirichlet_multinomial_3", "Smplnm") ]



comp_treat <- comp_treat %>% 
  melt(id.vars=c("Smplnm", "Dirichlet_multinomial_3", "timepoint"))
comp_treat$value <- as.numeric(comp_treat$value)

comp_treat %>% 
  filter(timepoint == "baseline") -> vb
vb <- vb %>% 
  filter(!Dirichlet_multinomial_3 == "NA")




helper_rel%>%
  filter(time == "baseline") %>% 
  #dplyr::group_by(Line)%>%
  ggplot(aes(y = rank(value),x= as.factor(variable), fill =as.factor(CST )))+
  geom_boxplot(
    aes(fill = as.factor(CST)), outlier.shape = NA
  ) +
  scale_fill_manual(values=col_vector_4)+
  scale_color_manual(values = col_vector)+
  theme_bw()+
  theme (axis.title.x = element_blank (), 
         axis.text.x = element_text (size = 13), 
         axis.text.y = element_text (size = 13), 
         axis.title.y = element_text (size = 13),
         legend.position = "none")+
  ggtitle("baseline")+
  ylab("ranked taxa count") -> v.1

grobs <- ggplotGrob(v.1)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

comp_treat %>% 
  filter(timepoint == "timepoint 2") -> vb
vb <- vb %>% 
  filter(!Dirichlet_multinomial_3 == "NA")


vb%>%
  #dplyr::group_by(Line)%>%
  ggplot(aes(y = rank(value),x= as.factor(variable), fill =as.factor(Dirichlet_multinomial_3 )))+
  geom_boxplot(
    aes(fill = as.factor(Dirichlet_multinomial_3)), outlier.shape = NA
  ) +
  scale_fill_manual(values=col_vector_4)+
  scale_color_manual(values = col_vector)+
  theme_bw()+
  theme (axis.title.x = element_blank (), 
         axis.text.x = element_text (size = 13), 
         axis.text.y = element_text (size = 13), 
         axis.title.y = element_text (size = 13),
         legend.position = "none")+
  ggtitle("timepoint 2")+
  ylab("ranked taxa count") -> v.2

comp_treat %>% 
  filter(timepoint == "timepoint 3") -> vb
vb <- vb %>% 
  filter(!Dirichlet_multinomial_3 == "NA")

###### Color palette for collections #########
library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


n <- 4
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector_4 = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))




vb%>%
  #dplyr::group_by(Line)%>%
  ggplot(aes(y = rank(value),x= as.factor(variable), fill =as.factor(Dirichlet_multinomial_3 )))+
  geom_boxplot(
    aes(fill = as.factor(Dirichlet_multinomial_3)), outlier.shape = NA
  ) +
  scale_fill_manual(values=col_vector_4)+
  scale_color_manual(values = col_vector)+
  theme_bw()+
  theme (axis.title.x = element_blank (), 
         axis.text.x = element_text (size = 13), 
         axis.text.y = element_text (size = 13), 
         axis.title.y = element_text (size = 13),
         legend.position = "none")+
  ggtitle("timepoint 3")+
  ylab("ranked taxa count") -> v.3

title_gg <- ggplot() + 
  labs(title = "Vagina Swabs", subtitle = "treatment with good compliance")

cowplot::plot_grid(v.1,v.2,v.3, align = "hv") ->vvv
bact <- cowplot::plot_grid(vvv, legend, rel_widths = c(1, .8))
cowplot::plot_grid(title_gg, bact, nrow = 2, rel_heights=c(0.1, 1))->vvv
ggsave(vvv, file = "../Figures/2023.figure/qpcr.vagina.types.svg", device = "svg", width = 22, height = 10)


######## control #############
comp_con <- comp %>% 
  filter(!treatment == "treatment")

comp_con <- comp_con[, c("timepoint", "L.Crispatus", "L.Gasseri",
                             "L.Jensenii", "L.Rhamnosus", "L.Iners",
                             "Dirichlet_multinomial_3", "Smplnm") ]



comp_con <- comp_con %>% 
  melt(id.vars=c("Smplnm", "Dirichlet_multinomial_3", "timepoint"))
comp_con$value <- as.numeric(comp_con$value)

comp_con %>% 
  filter(timepoint == "baseline") -> vb
vb <- vb %>% 
  filter(!Dirichlet_multinomial_3 == "NA")


vb%>%
  #dplyr::group_by(Line)%>%
  ggplot(aes(y = rank(value),x= as.factor(variable), fill =as.factor(Dirichlet_multinomial_3 )))+
  geom_boxplot(
    aes(fill = as.factor(Dirichlet_multinomial_3)), outlier.shape = NA
  ) +
  scale_fill_manual(values=col_vector_4)+
  scale_color_manual(values = col_vector)+
  theme_bw()+
  theme (axis.title.x = element_blank (), 
         axis.text.x = element_text (size = 13), 
         axis.text.y = element_text (size = 13), 
         axis.title.y = element_text (size = 13),
         legend.position = "none")+
  ggtitle("baseline")+
  ylab("ranked taxa count") -> v.1

grobs <- ggplotGrob(v.1)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

comp_con %>% 
  filter(timepoint == "timepoint 2") -> vb
vb <- vb %>% 
  filter(!Dirichlet_multinomial_3 == "NA")


vb%>%
  #dplyr::group_by(Line)%>%
  ggplot(aes(y = rank(value),x= as.factor(variable), fill =as.factor(Dirichlet_multinomial_3 )))+
  geom_boxplot(
    aes(fill = as.factor(Dirichlet_multinomial_3)), outlier.shape = NA
  ) +
  scale_fill_manual(values=col_vector_4)+
  scale_color_manual(values = col_vector)+
  theme_bw()+
  theme (axis.title.x = element_blank (), 
         axis.text.x = element_text (size = 13), 
         axis.text.y = element_text (size = 13), 
         axis.title.y = element_text (size = 13),
         legend.position = "none")+
  ggtitle("timepoint 2")+
  ylab("ranked taxa count") -> v.2

comp_con %>% 
  filter(timepoint == "timepoint 3") -> vb
vb <- vb %>% 
  filter(!Dirichlet_multinomial_3 == "NA")


vb%>%
  #dplyr::group_by(Line)%>%
  ggplot(aes(y = rank(value),x= as.factor(variable), fill =as.factor(Dirichlet_multinomial_3 )))+
  geom_boxplot(
    aes(fill = as.factor(Dirichlet_multinomial_3)), outlier.shape = NA
  ) +
  scale_fill_manual(values=col_vector_4)+
  scale_color_manual(values = col_vector)+
  theme_bw()+
  theme (axis.title.x = element_blank (), 
         axis.text.x = element_text (size = 13), 
         axis.text.y = element_text (size = 13), 
         axis.title.y = element_text (size = 13),
         legend.position = "none")+
  ggtitle("timepoint 3")+
  ylab("ranked taxa count") -> v.3

title_gg <- ggplot() + 
  labs(title = "Vagina Swabs", subtitle = "control samples only")

cowplot::plot_grid(v.1,v.2,v.3, align = "hv") ->vvv
bact <- cowplot::plot_grid(vvv, legend, rel_widths = c(1, .8))
cowplot::plot_grid(title_gg, bact, nrow = 2, rel_heights=c(0.1, 1))->vvv
ggsave(vvv, file = "../Figures/2023.figure/qpcr.vagina.types.control.svg", device = "svg", width = 22, height = 10)



library(DirichletMultinomial)
load(file = "output/ana.dmn_list.tsv")
Dirichlet_multinomial_1 = mixture(dmn_list[[1]], assign = TRUE)
Dirichlet_multinomial_2 = mixture(dmn_list[[2]], assign = TRUE)
Dirichlet_multinomial_3 = mixture(dmn_list[[3]], assign = TRUE)
Dirichlet_multinomial_4 = mixture(dmn_list[[4]], assign = TRUE)
Dirichlet_multinomial_5 = mixture(dmn_list[[5]], assign = TRUE)
Dirichlet_multinomial_6 = mixture(dmn_list[[6]], assign = TRUE)

Dirichlet_multinomial_all = data.frame(cbind(Dirichlet_multinomial_1,
                                             Dirichlet_multinomial_2,Dirichlet_multinomial_3,
                                             Dirichlet_multinomial_4,Dirichlet_multinomial_5,
                                             Dirichlet_multinomial_6))
colnames(Dirichlet_multinomial_all) = c("DMM_k=1","DMM_k=2","DMM_k=3",
                                        "DMM_k=4","DMM_k=5","DMM_k=6")
cluster_result <- as.data.frame(Dirichlet_multinomial_2)
cluster_result$ID <- row.names(cluster_result)
qpcr_red <- read.table(file = "qpcr.helper.tsv", sep ="\t", header = T, row.names = 1)

qpcr <- read.table(file = "../Figures/march.fig/helper/qpcr.ranked.helper.csv", sep = ",", header = T, row.names = 1)
# write.table(otu.gen, file = "../Figures/march.fig/helper/otu.gen.need.id.tsv", sep = "\t", quote = F)
qpcr$ID <- row.names(qpcr)
qpcr$Row.names <- NULL
qpcr$Sampl_ID <- NULL






########### anal #########
qpcr_ana <- qpcr %>% 
  filter(swabsite == "Ana_Swabs")

comp <- left_join(qpcr_ana, cluster_result, by = "ID")
comp <- comp %>% 
  filter(!compliance == "NA")


write.table(comp, file = "../Figures/march.fig/helper/cluster.merge.ana.tsv", sep = "\t", quote = F)





########### anal ###########
qpcr_red <- qpcr_red %>% 
  filter(!swabsite == "Vag_Swabs")

comp <- left_join(qpcr_red, cluster_result, by = "ID")
comp <- comp %>% 
  filter(!compliance == "NA")

######## treatment #############
comp_treat <- comp %>% 
  filter(!treatment == "control")

comp_treat <- comp_treat[, c("timepoint", "L.Crispatus", "L.Gasseri",
                             "L.Jensenii", "L.Rhamnosus", "L.Iners",
                             "Dirichlet_multinomial_4", "Smplnm") ]



comp_treat <- comp_treat %>% 
  melt(id.vars=c("Smplnm", "Dirichlet_multinomial_4", "timepoint"))
comp_treat$value <- as.numeric(comp_treat$value)

comp_treat %>% 
  filter(timepoint == "baseline") -> vb
vb <- vb %>% 
  filter(!Dirichlet_multinomial_4 == "NA")


vb%>%
  #dplyr::group_by(Line)%>%
  ggplot(aes(y = rank(value),x= as.factor(variable), fill =as.factor(Dirichlet_multinomial_4 )))+
  geom_boxplot(
    aes(fill = as.factor(Dirichlet_multinomial_4)), outlier.shape = NA
  ) +
  scale_fill_manual(values=col_vector_4)+
  scale_color_manual(values = col_vector)+
  theme_bw()+
  theme (axis.title.x = element_blank (), 
         axis.text.x = element_text (size = 13), 
         axis.text.y = element_text (size = 13), 
         axis.title.y = element_text (size = 13),
         legend.position = "none")+
  ggtitle("baseline")+
  ylab("ranked taxa count") -> v.1

grobs <- ggplotGrob(v.1)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

comp_treat %>% 
  filter(timepoint == "timepoint 2") -> vb
vb <- vb %>% 
  filter(!Dirichlet_multinomial_4 == "NA")


vb%>%
  #dplyr::group_by(Line)%>%
  ggplot(aes(y = rank(value),x= as.factor(variable), fill =as.factor(Dirichlet_multinomial_4 )))+
  geom_boxplot(
    aes(fill = as.factor(Dirichlet_multinomial_4)), outlier.shape = NA
  ) +
  scale_fill_manual(values=col_vector_4)+
  scale_color_manual(values = col_vector)+
  theme_bw()+
  theme (axis.title.x = element_blank (), 
         axis.text.x = element_text (size = 13), 
         axis.text.y = element_text (size = 13), 
         axis.title.y = element_text (size = 13),
         legend.position = "none")+
  ggtitle("timepoint 2")+
  ylab("ranked taxa count") -> v.2

comp_treat %>% 
  filter(timepoint == "timepoint 3") -> vb
vb <- vb %>% 
  filter(!Dirichlet_multinomial_4 == "NA")

###### Color palette for collections #########
library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


n <- 4
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector_4 = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))




vb%>%
  #dplyr::group_by(Line)%>%
  ggplot(aes(y = rank(value),x= as.factor(variable), fill =as.factor(Dirichlet_multinomial_4 )))+
  geom_boxplot(
    aes(fill = as.factor(Dirichlet_multinomial_4)), outlier.shape = NA
  ) +
  scale_fill_manual(values=col_vector_4)+
  scale_color_manual(values = col_vector)+
  theme_bw()+
  theme (axis.title.x = element_blank (), 
         axis.text.x = element_text (size = 13), 
         axis.text.y = element_text (size = 13), 
         axis.title.y = element_text (size = 13),
         legend.position = "none")+
  ggtitle("timepoint 3")+
  ylab("ranked taxa count") -> v.3

title_gg <- ggplot() + 
  labs(title = "Anal Swabs", subtitle = "treatment with good compliance")

cowplot::plot_grid(v.1,v.2,v.3, align = "hv") ->vvv
bact <- cowplot::plot_grid(vvv, legend, rel_widths = c(1, .8))
cowplot::plot_grid(title_gg, bact, nrow = 2, rel_heights=c(0.1, 1))->vvv
ggsave(vvv, file = "../Figures/2023.figure/qpcr.anal.types.svg", device = "svg", width = 22, height = 10)


######## control #############
comp_con <- comp %>% 
  filter(!treatment == "treatment")

comp_con <- comp_con[, c("timepoint", "L.Crispatus", "L.Gasseri",
                         "L.Jensenii", "L.Rhamnosus", "L.Iners",
                         "Dirichlet_multinomial_4", "Smplnm") ]



comp_con <- comp_con %>% 
  melt(id.vars=c("Smplnm", "Dirichlet_multinomial_4", "timepoint"))
comp_con$value <- as.numeric(comp_con$value)

comp_con %>% 
  filter(timepoint == "baseline") -> vb
vb <- vb %>% 
  filter(!Dirichlet_multinomial_4 == "NA")


vb%>%
  #dplyr::group_by(Line)%>%
  ggplot(aes(y = rank(value),x= as.factor(variable), fill =as.factor(Dirichlet_multinomial_4 )))+
  geom_boxplot(
    aes(fill = as.factor(Dirichlet_multinomial_4)), outlier.shape = NA
  ) +
  scale_fill_manual(values=col_vector_4)+
  scale_color_manual(values = col_vector)+
  theme_bw()+
  theme (axis.title.x = element_blank (), 
         axis.text.x = element_text (size = 13), 
         axis.text.y = element_text (size = 13), 
         axis.title.y = element_text (size = 13),
         legend.position = "none")+
  ggtitle("baseline")+
  ylab("ranked taxa count") -> v.1

grobs <- ggplotGrob(v.1)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

comp_con %>% 
  filter(timepoint == "timepoint 2") -> vb
vb <- vb %>% 
  filter(!Dirichlet_multinomial_4 == "NA")


vb%>%
  #dplyr::group_by(Line)%>%
  ggplot(aes(y = rank(value),x= as.factor(variable), fill =as.factor(Dirichlet_multinomial_4 )))+
  geom_boxplot(
    aes(fill = as.factor(Dirichlet_multinomial_4)), outlier.shape = NA
  ) +
  scale_fill_manual(values=col_vector_4)+
  scale_color_manual(values = col_vector)+
  theme_bw()+
  theme (axis.title.x = element_blank (), 
         axis.text.x = element_text (size = 13), 
         axis.text.y = element_text (size = 13), 
         axis.title.y = element_text (size = 13),
         legend.position = "none")+
  ggtitle("timepoint 2")+
  ylab("ranked taxa count") -> v.2

comp_con %>% 
  filter(timepoint == "timepoint 3") -> vb
vb <- vb %>% 
  filter(!Dirichlet_multinomial_4 == "NA")


vb%>%
  #dplyr::group_by(Line)%>%
  ggplot(aes(y = rank(value),x= as.factor(variable), fill =as.factor(Dirichlet_multinomial_4 )))+
  geom_boxplot(
    aes(fill = as.factor(Dirichlet_multinomial_4)), outlier.shape = NA
  ) +
  scale_fill_manual(values=col_vector_4)+
  scale_color_manual(values = col_vector)+
  theme_bw()+
  theme (axis.title.x = element_blank (), 
         axis.text.x = element_text (size = 13), 
         axis.text.y = element_text (size = 13), 
         axis.title.y = element_text (size = 13),
         legend.position = "none")+
  ggtitle("timepoint 3")+
  ylab("ranked taxa count") -> v.3

title_gg <- ggplot() + 
  labs(title = "Anal Swabs", subtitle = "control samples only")

cowplot::plot_grid(v.1,v.2,v.3, align = "hv") ->vvv
bact <- cowplot::plot_grid(vvv, legend, rel_widths = c(1, .8))
cowplot::plot_grid(title_gg, bact, nrow = 2, rel_heights=c(0.1, 1))->vvv
ggsave(vvv, file = "../Figures/2023.figure/qpcr.anal.types.control.svg", device = "svg", width = 22, height = 10)

# - It looks like essentially CST 1 is replaced by CST 3 as time passes - is that the case? Could one make a "cross-over" matrix heatmap (cells show fraction/percentage having CST X at time point A, CST Y at time point B so that diagonal is no change/persistence)? Is this a general shift towards CST 3 or just that CST 1 becomes CST 3 (meaning, CST 2 typically stays that way)?

# What of a heatmap where X axis is CST, Y axis is qPCR species, and cell values are ranked taxa counts as per below? That might show clearer if there is a skew between CSTs and which species score high in qPCR? Do you see how I mean?

##### heatmap
qpcr_red <- read.table(file = "qpcr.helper.tsv", sep ="\t", header = T, row.names = 1)

########### vagina #########
qpcr_red <- qpcr_red %>% 
  filter(!swabsite == "Ana_Swabs")

comp <- left_join(qpcr_red, cluster_result, by = "ID")
comp <- comp %>% 
  filter(!compliance == "NA")

######## treatment #############
comp_treat <- comp %>% 
  filter(!treatment == "control")

comp_treat <- comp_treat[, c("timepoint", "L.Crispatus", "L.Gasseri",
                             "L.Jensenii", "L.Rhamnosus", "L.Iners",
                             "Dirichlet_multinomial_3", "Smplnm") ]



comp_treat <- comp_treat %>% 
  melt(id.vars=c("Smplnm", "Dirichlet_multinomial_3", "timepoint"))
comp_treat$value <- as.numeric(comp_treat$value)

comp_treat <- comp_treat %>% 
  filter(!Dirichlet_multinomial_3 == "NA")

comp_treat <- comp_treat %>% 
  filter(!value == "NA")

comp_treat$value <- rank(comp_treat$value)

write.table(comp_treat, file = "../helperbaseline.vaginalheatb.tsv", sep = "\t", quote = F)


vaginalheatmap <- read.table(file = "../summed_values_control_vagina.csv", sep = ",", header = T)



###### relatice ####
helper_rel <- read.table(file = "../Figures/march.fig/helper/vagina_sums_relatives.csv", sep =",", header = T)



  helper_rel%>%
  ggplot (aes (x =as.factor(CST) , fill = round(value, digits = 2), y = variable), colour = "black") +
  geom_tile(aes(fill = value), colour = "black") +
  facet_wrap(~time, ncol = 3) +
  scale_fill_gradient(high = "#822bc1ff", low =  "#cea5eb", na.value = "white") +
  geom_text (aes (label = round(value, digits = 2)), size = 6)+
  labs(x = "CST", y = "", fill = "relative counts")+
  scale_x_discrete(labels = c("CST 1","CST 2", "CST 3"))+
  theme_classic()+
  theme (axis.title.x = element_blank (), 
         text = element_text(size = 15),
         axis.text.x = element_text (size = 15, angle = 90), 
         axis.text.y = element_text (size = 15), 
         axis.title.y = element_text (size = 15),
         legend.position = "none")-> vag
  

helper_rel <- read.table(file = "../Figures/march.fig/helper/ana_sum_relatives.csv", sep =",", header = T)
  
helper_rel %>%
    ggplot (aes (x =as.factor(CST) , fill = round(value, digits = 2), y = variable), colour = "black") +
    geom_tile(aes(fill = value), colour = "black") +
    facet_wrap(~time, ncol = 3) +
    scale_fill_gradient(high =  "#803300ff", low =  "#ffa366", na.value = "white") +
    geom_text (aes (label = round(value, digits = 2)), size = 6)+
    labs(x = "CST", y = "", fill = "relative counts")+
    scale_x_discrete(labels = c("community type 1","community type 2"))+
    theme_classic()+
    theme (axis.title.x = element_blank (), 
           text = element_text(size = 15),
           axis.text.x = element_text (size = 15, angle = 90), 
           axis.text.y = element_text (size = 15), 
           axis.title.y = element_text (size = 15),
           legend.position = "none")-> anal
  
cowplot::plot_grid(vag, anal, nrow = 1) ->totalheat
ggsave(totalheat, file = "../Figures/march.fig/figure/total_rel_heatmap_qpcr.svg",device = "svg", height = 8, width = 12)  


comp_treat <- comp %>% 
  filter(!treatment == "treatment")

comp_treat <- comp_treat[, c("timepoint", "L.Crispatus", "L.Gasseri",
                             "L.Jensenii", "L.Rhamnosus", "L.Iners",
                             "Dirichlet_multinomial_3", "Smplnm") ]



comp_treat <- comp_treat %>% 
  melt(id.vars=c("Smplnm", "Dirichlet_multinomial_3", "timepoint"))
comp_treat$value <- as.numeric(comp_treat$value)

comp_treat <- comp_treat %>% 
  filter(!Dirichlet_multinomial_3 == "NA")

comp_treat <- comp_treat %>% 
  filter(!value == "NA")

comp_treat$value <- rank(comp_treat$value)

write.table(comp_treat, file = "../helpertreat.vaginalheatb.tsv", sep = "\t", quote = F)


vaginalheatmap <- read.table(file = "../summed_values_treatment_vagina.csv", sep = ",", header = T)
vaginalheatmap <- vaginalheatmap %>% 
  filter(!Dirichlet_multinomial_3 == "NA")

vaginalheatmap %>% 
  ggplot (aes (x =as.factor(Dirichlet_multinomial_3) , y = variable, fill = summed_values), colour = "black") +
  geom_tile(aes(fill = summed_values), colour = "black") +
  facet_wrap(~timepoint, ncol = 3) +
  scale_fill_gradient(low = "#80CDC1", high = "#F1B6DA", na.value = "white") +
  geom_text (aes (label = summed_values), size = 6)+
  labs(x = "community type", y = "", fill = "ranked bacterial\ncopynumbers")+
  scale_x_discrete(labels = c("community type 1","community type 2", "community type 3"))+
  theme_classic()+
  theme (axis.title.x = element_blank (), 
         text = element_text(size = 15),
         axis.text.x = element_text (size = 15, angle = 90), 
         axis.text.y = element_text (size = 15), 
         axis.title.y = element_text (size = 15),
         legend.position = "none")+
  ggtitle("vaginal swabs treatment group") -> treatment.vag


vagina <- cowplot::plot_grid(contr.vag, treatment.vag, align = "hv")
ggsave(vagina, file = "../Figures/2023.figure/heatmap.qpcr.vaginatype.svg", width = 14, height = 5, device = "svg")



########### Anal #########
qpcr_red <- qpcr_red %>% 
  filter(swabsite == "Ana_Swabs")

comp <- left_join(qpcr_red, cluster_result, by = "ID")
comp <- comp %>% 
  filter(!compliance == "NA")

######## treatment #############
comp_treat <- comp %>% 
  filter(treatment == "control")

comp_treat <- comp_treat[, c("timepoint", "L.Crispatus", "L.Gasseri",
                             "L.Jensenii", "L.Rhamnosus", "L.Iners",
                             "Dirichlet_multinomial_4", "Smplnm") ]



comp_treat <- comp_treat %>% 
  melt(id.vars=c("Smplnm", "Dirichlet_multinomial_4", "timepoint"))
comp_treat$value <- as.numeric(comp_treat$value)

comp_treat <- comp_treat %>% 
  filter(!Dirichlet_multinomial_4 == "NA")

comp_treat <- comp_treat %>% 
  filter(!value == "NA")

comp_treat$value <- rank(comp_treat$value)

write.table(comp_treat, file = "../control.anal.helper.tsv", sep = "\t", quote = F)


analheatmap <- read.table(file = "../summed_control.anal.csv", sep = ",", header = T)

analheatmap %>% 
  ggplot (aes (x =as.factor(Dirichlet_multinomial_4) , y = variable, fill = summed_values), colour = "black") +
  geom_tile(aes(fill = summed_values), colour = "black") +
  facet_wrap(~timepoint, ncol = 3) +
  scale_fill_gradient(low = "#80CDC1", high = "#F1B6DA", na.value = "white") +
  geom_text (aes (label = summed_values))+
  labs(x = "community type", y = "", fill = "ranked bacterial\ncopynumbers")+
  theme_classic()+
  theme (axis.title.x = element_blank (), 
         axis.text.x = element_text (size = 13),
         text = element_text(size = 15),
         axis.text.y = element_text (size = 13), 
         axis.title.y = element_text (size = 13),
         legend.position = "none")+
  ggtitle("rectal swabs control group") -> contr.ana

comp_treat <- comp %>% 
  filter(treatment == "treatment")

comp_treat <- comp_treat[, c("timepoint", "L.Crispatus", "L.Gasseri",
                             "L.Jensenii", "L.Rhamnosus", "L.Iners",
                             "Dirichlet_multinomial_4", "Smplnm") ]



comp_treat <- comp_treat %>% 
  melt(id.vars=c("Smplnm", "Dirichlet_multinomial_4", "timepoint"))
comp_treat$value <- as.numeric(comp_treat$value)

comp_treat <- comp_treat %>% 
  filter(!Dirichlet_multinomial_4 == "NA")

comp_treat <- comp_treat %>% 
  filter(!value == "NA")

comp_treat$value <- rank(comp_treat$value)

write.table(comp_treat, file = "../helpertreat.anal.tsv", sep = "\t", quote = F)


analheatmap <- read.table(file = "../summed_values_treatment_anal.csv", sep = ",", header = T)
analheatmap <- analheatmap %>% 
  filter(!Dirichlet_multinomial_4 == "NA")

analheatmap %>% 
  ggplot (aes (x =as.factor(Dirichlet_multinomial_4) , y = variable, fill = summed_values), colour = "black") +
  geom_tile(aes(fill = summed_values), colour = "black") +
  facet_wrap(~timepoint, ncol = 3) +
  scale_fill_gradient(low = "#80CDC1", high = "#F1B6DA", na.value = "white") +
  geom_text (aes (label = summed_values))+
  labs(x = "community type", y = "", fill = "ranked bacterial\ncopynumbers")+
  theme_classic()+
  theme (axis.title.x = element_blank (), 
         text = element_text(size = 15),
         axis.text.x = element_text (size = 15), 
         axis.text.y = element_text (size = 15), 
         axis.title.y = element_text (size = 15),
         legend.position = "none")+
  ggtitle("rectal swabs treatment group") -> treatment.ana


anal <- cowplot::plot_grid(contr.ana, treatment.ana, align = "hv")
ggsave(anal, file = "../Figures/2023.figure/heatmap.qpcr.analtype.svg", width = 18, height = 4, device = "svg")



########## privious preterm birth ##########
########## check for both groups had high L. Crisp or high L. Iners GA<259
qpcr_red <- read.table(file = "qpcr.helper.tsv", sep ="\t", header = T, row.names = 1)

qpcr_ga <- left_join(qpcr_red, meta_Time, by = "Sampl_ID")
qpcr_ga <- qpcr_ga %>% 
  filter(Swabsite == "Vag_Swabs")

qpcr_ga$GA_at_birth_days[qpcr_ga$GA_at_birth_days < 259] <- 1
qpcr_ga$GA_at_birth_days[qpcr_ga$GA_at_birth_days > 259] <- 0
qpcr_ga$GA_at_birth_days <- as.factor(qpcr_ga$GA_at_birth_days )
###### L.Crispatus ########
qpcr_ga %>% 
  filter(GA_at_birth_days == "1") %>% 
  # filter(treatment == "control") %>% 
  # dplyr::group_by(timepoint)%>%
  rstatix::wilcox_test(L.Iners ~ timepoint)%>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance()%>%
  rstatix::add_xy_position(x = "timepoint") -> stats_test

qpcr_ga %>% 
  filter(!GA_at_birth_days == "NA") %>% 
  # filter(treatment == "control") %>% 
  dplyr::group_by(timepoint)%>%
  rstatix::wilcox_test(L.Iners ~ GA_at_birth_days)%>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance()%>%
  rstatix::add_xy_position(x = "GA_at_birth_days") -> stats_test

# L.chrispatus/L.Gasseri/L.Jensenii/Iners and Rh. ns both between 0/1 at the different timepoint and betwenn the different timepoint when 1
qpcr_ga %>%
  filter(GA_at_birth_days == "1") -> ga.pcr.1
row.names(ga.pcr.1)<- ga.pcr.1$Sampl_ID
ga.pcr.1 <- ga.pcr.1[, c(1,6:12)]
gen.red <- otu.gen %>%
  select((
    c(
      "OTU_26_Prevotella.9",
      "OTU_18_Prevotella.6",
      "OTU_38_Prevotella",
      "OTU_198_Prevotella.7",
      "OTU_8_Gardnerella",
      "OTU_45_Atopobium",
      "OTU_117_Sneathia",
      "OTU_392_Ureaplasma",
      "OTU_785_Mycoplasma",
      "OTU_1_Lactobacillus"
    )
  ))


ga.pcr.1 %>% 
  reshape2::melt(id.vars=c("Smplnm", "timepoint")) -> ga.pcr.1

ga.pcr.1 %>% 
  ggplot(aes(x = timepoint, y = value, color = timepoint, fill = timepoint)) +
  geom_boxplot(color = "black") +
  geom_point(aes(fill = timepoint, group=Smplnm),size=3,shape=21, position = position_dodge(0.2))+
  facet_wrap(~variable, ncol = 3) +
  geom_line(aes(group=Smplnm)) +
  theme_classic() +
  ylab("log (copynumber)")+
  scale_fill_manual(values = alpha(pal.vag, 0.8)) +
  scale_colour_manual(values = pal.vag) +
  theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size = 15), axis.title.x = element_blank(),
         axis.text.y = element_text(size = 15),
         axis.title.y = element_text ( size = 15),
         legend.position = "none")  -> qpcr.vag.ga.1



#### qpcr #####
ga_otu.gen <- merge(meta_Time, gen.red, by = 0)
ga_otu.gen <- ga_otu.gen %>% 
  filter(Swabsite == "Vag_Swabs")

ga_otu.gen$GA_at_birth_days[ga_otu.gen$GA_at_birth_days < 259] <- 1
ga_otu.gen$GA_at_birth_days[ga_otu.gen$GA_at_birth_days > 259] <- 0
ga_otu.gen$GA_at_birth_days <- as.factor(ga_otu.gen$GA_at_birth_days )

ga_otu.gen %>% 
  filter(!GA_at_birth_days == "NA") %>% 
  # filter(treatment == "control") %>% 
  dplyr::group_by(meta_Time)%>%
  rstatix::wilcox_test(OTU_1_Lactobacillus ~ GA_at_birth_days)%>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance()%>%
  rstatix::add_xy_position(x = "GA_at_birth_days") -> stats_test
# OTU_26_Prevotella.9 not working/OTU_18_Prevotella.6/OTU_38_Prevotella

ga_otu.gen %>% 
  filter(GA_at_birth_days == "1") %>% 
  # filter(treatment == "control") %>% 
  # dplyr::group_by(timepoint)%>%
  rstatix::wilcox_test(OTU_8_Gardnerella ~ meta_Time)%>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance()%>%
  rstatix::add_xy_position(x = "meta_Time") -> stats_test




#### qpcr #####
ga.pcr.1 <- ga.pcr.1[, c(1,6:12)]

ga_otu.gen %>% 
  reshape2::melt(id.vars=c("Smplnm", "timepoint")) -> ga_otu.gen


ga_otu.gen %>% 
  ggplot(aes(x = timepoint, y = value, color = timepoint, fill = timepoint)) +
  geom_boxplot(color = "black") +
  geom_point(aes(fill = timepoint, group=Smplnm),size=3,shape=21, position = position_dodge(0.2))+
  facet_wrap(~variable, ncol = 3) +
  geom_line(aes(group=Smplnm)) +
  theme_classic() +
  ylab("read counts")+
  scale_fill_manual(values = alpha(pal.vag, 0.8)) +
  scale_colour_manual(values = pal.vag) +
  theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size = 15), axis.title.x = element_blank(),
         axis.text.y = element_text(size = 15),
         axis.title.y = element_text ( size = 15),
         legend.position = "none")  -> otu.vag.ga.1




bact <- cowplot::plot_grid(otu.vag.ga.1, qpcr.vag.ga.1, align = "hv")


ggsave(bact, width = 10, height = 8, dpi = 300, device = "svg", filename = "../Figures/2023.figure/ga.1.request.leonoar.svg")


qpcr <- read.table(file = "../Figures/march.fig/input/qpcr_final.csv", header = T, sep =",", row.names = 1)

# Treatment == 1, Control == 0
# complinance = baseline/non vs. good 2/3
# 
qpcr$treatment[qpcr$treatment == "1"] <- "treatment"
qpcr$treatment[qpcr$treatment == "0"] <- "control"

qpcr$compliance[qpcr$compliance == "Na"] <- NA

qpcr$timepoint[qpcr$timepoint == "1"] <- "baseline"
qpcr$timepoint[qpcr$timepoint == "2"] <- "timepoint 2"
qpcr$timepoint[qpcr$timepoint == "3"] <- "timepoint 3"
qpcr <- qpcr %>% 
  filter(!compliance == "NA")

## Rank values but keep NAs percentile rank function in R
prank<-function(x){
  r<-rank(x)/sum(!is.na(x))
  r[is.na(x)]<-NA
  r
}
## normal rank without percentile
prank<-function(x) ifelse(is.na(x),NA,rank(x))
for (i in 6:11) {
  qpcr[, i] <- prank(qpcr[, i])
} 




qpcr.presence <- qpcr[, c(7:11)]
x <- qpcr[, c(1:5,12)]
qpcr.presence[qpcr.presence > 1] <- 1
qpcr.presence[is.na(qpcr.presence)] <- 0

qpcr.presence <- merge(qpcr.presence, x, by = 0)
qpcr.presence$Row.names <- NULL
qpcr.presence$Samplename <- NULL
qpcr.presence$treatment <- NULL


qpcr.presence %>% 
  reshape2::melt(id.vars=c("ID", "timepoint","swabsite", "compliance", "treatment")) -> qpcr.presence

library(dplyr)

# Filter the data to include only the control group
control_data <- qpcr.presence %>% filter(compliance == "control")
good_data <- qpcr.presence %>% filter(compliance == "good")
# Group the data by taxa and summarize the total bacteria count
taxa_summary <- control_data %>% 
  group_by(variable, timepoint, swabsite) %>% 
  summarize(total_bacteria = sum(value))

taxa_summary.g <- good_data %>% 
  group_by(variable, timepoint, swabsite) %>% 
  summarize(total_bacteria = sum(value))


# If the test is for presence/absence, then for each species which is in the probiotic, you can check for each sample if it is there in the qPCR or not. You can then tally this up as a contingency table: found/not found on one axis, compliant/not compliant on another axis. If one does this separately for time point 2 and time point 3, there is no pseudoreplication, and then this can be tested just with a Chi-squared test. Does this achieve significance? 
#   
#   It could also be tested with a logistic regression mixed effects model:
#   
#   m1 = presence_qpcr ~ compliance * time_point + (1 | individual) + time_point + compliance, error model for logistic regression
# m0 = presence_qpcr ~ (1 | individual) + time_point, error model for logistic regression
# 
# and compare those models.
# 
# Ultimately, unless you can show significance statistically, you can only report trends which will depend on what angle you are looking from.


qpcr.presence <- qpcr.presence %>% 
  filter(!variable == "total_bacterialload")

qpcr.presence.vag <- qpcr.presence %>% 
  filter(!swabsite == "Ana_Swabs")


control_data <- qpcr.presence.vag %>% filter(compliance == "control")
good_data <- qpcr.presence %>% filter(compliance == "good")
# Group the data by taxa and summarize the total bacteria count
taxa_summary <- qpcr.presence.vag %>% 
  group_by(variable, timepoint, compliance) %>% 
  summarize(total_bacteria = sum(value))


ggplot(taxa_summary, aes(x = compliance, y = variable, fill = total_bacteria), color = "black") +
  geom_tile(color = "black") +
  facet_wrap(~ timepoint, ncol = 1) +
  scale_fill_gradient(high = "#822bc1ff", low =  "#cea5eb", na.value = "white")+
  geom_text (aes (label = round(total_bacteria, digits = 2)), size = 6)+
  theme_classic()+
  theme(text = element_text(size = 15),
        axis.text.x = element_text (size = 15), 
        axis.text.y = element_text (size = 15, face = "italic"), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) -> pres.vag



qpcr.presence.23 <- qpcr.presence.vag %>% 
  filter(!timepoint == "baseline")

qpcr.presence.b2 <- qpcr.presence.vag %>% 
  filter(!timepoint == "timepoint 3")

qpcr.presence.b3 <- qpcr.presence.vag %>% 
  filter(!timepoint == "timepoint 2")

library(lmtest)
library(lme4)
qpcr.presence.b3$timepoint[qpcr.presence.b3$timepoint == "timepoint 3"] <- 1
qpcr.presence.b3$timepoint[qpcr.presence.b3$timepoint == "baseline"] <- 0
qpcr.presence.b3$timepoint <- as.numeric(qpcr.presence.b3$timepoint)

qpcr.presence.b3$compliance[qpcr.presence.b3$compliance == "control"] <- 0
qpcr.presence.b3$compliance[qpcr.presence.b3$compliance == "good"] <- 1
qpcr.presence.b3$compliance <- as.numeric(qpcr.presence.b3$compliance)
qpcr.presence.b3 <- qpcr.presence.b3 %>% 
  filter(variable == "L.Iners")

logistic_regression <- function(data) {
  model1 <- glmer(value ~ timepoint + compliance + (1 | ID), family = binomial, data = data)
  model2 <- glmer(value ~ timepoint + (1 | ID), family = binomial, data = data)
  return(lrtest(model1, model2))$'Pr(>Chisq)'[[1]]
}

logistic_regression(qpcr.presence.23) 
# L.iners   0.8659
# L.Crispatus 0.9417
# L.Gasseri    0.939
# L.Jensenii 0.008144 **
# L.Rhamnosus 0.3466

logistic_regression(qpcr.presence.b2) 
# L.iners   0.8659
# L.Crispatus 0.9155
# L.Gasseri     0.9934
# L.Jensenii  0.09977
# L.Rhamnosus 0.6418

logistic_regression(qpcr.presence.b3) 
# L.iners   0.8893
# L.Crispatus 0.9951
# L.Gasseri     0.9786
# L.Jensenii  0.1501
# L.Rhamnosus 0.000109 ***




qpcr.presence.rec <- qpcr.presence %>% 
  filter(!swabsite == "Vag_Swabs")

taxa_summary <- qpcr.presence.rec %>% 
  group_by(variable, timepoint, compliance) %>% 
  summarize(total_bacteria = sum(value))

ggplot(taxa_summary, aes(x = compliance, y = variable, fill = total_bacteria), color = "black") +
  geom_tile(color = "black") +
  facet_wrap(~ timepoint, ncol = 1) +
  scale_fill_gradient(high =  "#803300ff", low =  "#ffa366", na.value = "white") +
  geom_text (aes (label = round(total_bacteria, digits = 2)), size = 6)+
  theme_classic()+
  theme(text = element_text(size = 15),
        axis.text.x = element_text (size = 15), 
        axis.text.y = element_text (size = 15, face = "italic"), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) -> pres.rec

cowplot::plot_grid(pres.vag, pres.rec) -> pres.abs
ggsave(pres.abs, file = "../Figures/march.fig/figure/pres.abs.qpcr.svg", device = "svg", width = 8, height = 8)
qpcr.presence.23 <- qpcr.presence.rec %>% 
  filter(!timepoint == "baseline")

qpcr.presence.b2 <- qpcr.presence.rec %>% 
  filter(!timepoint == "timepoint 3")

qpcr.presence.b3 <- qpcr.presence.rec %>% 
  filter(!timepoint == "timepoint 2")

library(lmtest)
library(lme4)
qpcr.presence.b3$timepoint[qpcr.presence.b3$timepoint == "timepoint 3"] <- 1
qpcr.presence.b3$timepoint[qpcr.presence.b3$timepoint == "baseline"] <- 0
qpcr.presence.b3$timepoint <- as.numeric(qpcr.presence.b3$timepoint)

qpcr.presence.b3$compliance[qpcr.presence.b3$compliance == "control"] <- 0
qpcr.presence.b3$compliance[qpcr.presence.b3$compliance == "good"] <- 1
qpcr.presence.b3$compliance <- as.numeric(qpcr.presence.b3$compliance)
qpcr.presence.b3 <- qpcr.presence.b3 %>% 
  filter(variable == "L.Rhamnosus")

logistic_regression <- function(data) {
  model1 <- glmer(value ~ timepoint + compliance + (1 | ID), family = binomial, data = data)
  model2 <- glmer(value ~ timepoint + (1 | ID), family = binomial, data = data)
  return(lrtest(model1, model2))$'Pr(>Chisq)'[[1]]
}

logistic_regression(qpcr.presence.23) 
# L.iners    0.8314
# L.Crispatus 0.9706
# L.Gasseri    0.07867
# L.Jensenii 0.005644 **
# L.Rhamnosus 0.02277 *

logistic_regression(qpcr.presence.b2) 
# L.iners   0.8482
# L.Crispatus .08195
# L.Gasseri     2.279e-05 ***
# L.Jensenii  0.247
# L.Rhamnosus 0.01472 *

logistic_regression(qpcr.presence.b3) 
# L.iners     0.7368
# L.Crispatus 0.1196
# L.Gasseri     8.754e-06 ***
# L.Jensenii  0.2654
# L.Rhamnosus 0.000109 ***