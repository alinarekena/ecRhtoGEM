# This script takes the enzyme usage data from ec-models and plots this in various ways
#install.packages("tidyverse") # Install tidyverse if required
library(tidyverse)             # for write_delim
#install.packages("magrittr")
library(dplyr)                 # for mutate_if
#install.packages("dplyr")
#install.packages("readr")
#library(readr)
#library(tidyr)
#library(ggplot2)
# Adjust to the correct directory 
setwd("~/Documents/GitHub/ecRhtoGEM/results/modelSimulation")  # directory path for Mac computers

# Load usage information and remove proteins with always zero usage
capUse <- read.delim('enzymeUsages.txt')
capUse <- capUse[,1:5]
capSum <- rowSums(capUse[,4:5])
capUse <- capUse[!capSum==0,]
capUse[,4:5] <- capUse[,4:5]*100

# Plot based on GO term annotation as obtained from Uniprot (first rearrange data)
GO <- read.delim('../../data/selectedAnnotation.txt', stringsAsFactors = F)

# Only keep data from relevant GO terms
capUse <- capUse[capUse$protID %in% GO$Entry,]
idx <- match(capUse$protID, GO$Entry)
capUse$GOterm <- GO$system[idx]
colnames(capUse) <- gsub('capUse_','',colnames(capUse))

capUse <- capUse %>% mutate_if(is.numeric, round, digits = 3)
write_delim(capUse,'../../results/modelSimulation/capUsage.txt',delim = '\t')

capUse <- gather(capUse, 'Condition', 'Usage', 4:5)
capUse$GOterm <- factor(capUse$GOterm, levels=c('glycolysis','TCA cycle','ETC','PPP','Ribosome'))
capUse$Condition <- factor(capUse$Condition, levels=c('Xexp','XNlim'))

plot1<-capUse[capUse$GOterm %in% c('glycolysis','TCA cycle','ETC','PPP','Ribosome'),]

#plot1_annot <- plot1 %>%
#  count(GOterm) %>%
#  mutate(
#    label = paste0("n = ", n)
#  )

ggplot(plot1, aes(x = Condition, y = Usage, color=GOterm)) +
  geom_boxplot(lwd = 0.35) +
#  geom_text(aes(Condition, y = Inf, label = label), data = plot1_annot, hjust = 0.05, vjust = +1) +
  scale_color_manual(values=c('#CBBBA0','#1D1D1B','#1D71B8','#878787','#FF0000')) +
  facet_grid(. ~ GOterm) +
  labs(x = '', y = 'Capacity usage (%)') + 
  theme_classic() +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5), text = element_text(size=7), 
        line = element_line(size=0.15), strip.background = element_blank(),
        axis.line = element_line(size=0.15), legend.position='none')
ggsave("selectedGOtermUsage.pdf", width=10, height=4.5, units='cm')

#plot2<-capUse[capUse$GOterm %in% c('PP shunt','THF cycle'),]
#ggplot(plot2, aes(x = Condition, y = Usage, color=GOterm)) +
#  geom_boxplot(lwd = 0.35) +
#  scale_color_manual(values=c('#CBBBA0','#1D1D1B','#1D71B8','#878787')) +
#  facet_grid(. ~ GOterm) +
#  labs(x = '', y = 'Capacity usage (%)') + 
#  theme_classic() +
#  theme(axis.text.x=element_text(angle = 90, vjust = 0.5), text = element_text(size=7), 
#        line = element_line(size=0.15), strip.background = element_blank(),
#        axis.line = element_line(size=0.15), legend.position='none')
#ggsave("supplementGOtermUsage.pdf", width=6, height=4.5, units='cm')

