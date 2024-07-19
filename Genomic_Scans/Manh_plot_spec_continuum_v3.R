rm(list = ls())
library(qqman)
library(tidyverse)
library(beepr)
library(dplyr)
library(VennDiagram)
library(ggplot2)
library(ggpubr)
library(gridExtra)

setwd("/mnt/DATA/B_JUN_reseq/B_PopIndexes/H_PLOTS_2.0/A_data")

# Save the entire workspace
save.image(file = "DEJUYEJU_spec_cont_v3.RData")
## Load the saved workspace
load("DEJUYEJU_spec_cont_v3.RData")

#pdf(file="/mnt/DATA/B_JUN_reseq/B_PopIndexes/A_Plots/Man_plot_CANDOR.pdf")
#jpeg(file="CANDOR.jpeg")
#png(file="CANDOR.png",    width=600, height=350)

#Carga datos y limpia NAs. El input es es output de Pixy
fst_sp<-read_tsv("sp_DE_YE_50kb_fst.txt")
fst_lin<-read_tsv("lin_50kb_fst.txt")
fst_pop<-read_tsv("POP_50kb_fst.txt")
dxy_sp<-read_tsv("sp_DE_YE_50kb_dxy.txt")
dxy_lin<-read_tsv("lin_50kb_dxy.txt")
dxy_pop<-read_tsv("POP_50kb_dxy.txt")
pi_sp <- read_tsv("sp_DE_YE_50kb_pi.txt")
pi_lin <- read_tsv("lin_50kb_pi.txt")
pi_pop <- read_tsv("POP_50kb_pi.txt")
xpehh_SP<-read_tsv("xpehh_DEJUYEJU.txt")
rmap<-read_table("/mnt/DATA/B_JUN_reseq/N_genetic_map/C_final_LDhat/A_HYE/Recombination_map_rho_HYE.txt")

fst_sp<-na.omit(fst_sp)
fst_lin<-na.omit(fst_lin)
fst_pop<-na.omit(fst_pop)
rmap<-na.omit(rmap)
xpehh<-na.omit(xpehh_SP)

#Utils:
win<- (dxy_sp$window_pos_2[1]-dxy_sp$window_pos_1[1])/2
ChrOrder<-c("ScoVZU6_963__HRSCAF___984","ScoVZU6_1238__HRSCAF___1266", "ScoVZU6_2668__HRSCAF___2714", "ScoVZU6_1491__HRSCAF___1522", "ScoVZU6_2655__HRSCAF___2701",  "ScoVZU6_1225__HRSCAF___1252","ScoVZU6_3544__HRSCAF___3601", "ScoVZU6_982__HRSCAF___1003", "ScoVZU6_3251__HRSCAF___3306", "ScoVZU6_2192__HRSCAF___2234", "ScoVZU6_201__HRSCAF___207", "ScoVZU6_4447__HRSCAF___4518", "ScoVZU6_3582__HRSCAF___3639", "ScoVZU6_4448__HRSCAF___4521", "ScoVZU6_2057__HRSCAF___2097", "ScoVZU6_3803__HRSCAF___3864", "ScoVZU6_4455__HRSCAF___4539", "ScoVZU6_1929__HRSCAF___1968", "ScoVZU6_582__HRSCAF___593", "ScoVZU6_2973__HRSCAF___3026", "ScoVZU6_1043__HRSCAF___1065", "ScoVZU6_3047__HRSCAF___3100", "ScoVZU6_1118__HRSCAF___1143", "ScoVZU6_300__HRSCAF___308", "ScoVZU6_957__HRSCAF___978", "ScoVZU6_1653__HRSCAF___1687", "ScoVZU6_3291__HRSCAF___3346", "ScoVZU6_3522__HRSCAF___3579", "ScoVZU6_4457__HRSCAF___4541", "ScoVZU6_4456__HRSCAF___4540", "ScU1KuS_2__HRSCAF___21")
ChrNum<- c(seq(1:31))
ChrName<-c("1","1A","2","3","4","4A","5","6", "7","8","9","10","11","12","13","14","15","17","18","19","20","21","22","23","24","25","26","27","28","29","Z")
ChrTable<-tibble(ChrOrder,ChrNum,ChrName)

#Ampl?a datos y pone todo lo necesario para plotear tanto Fst como Dxy
#### Fst and Dxy SP
FstDxy_sp <- dxy_sp %>%
  left_join(fst_sp) %>% #Junta Fst y Dxy en un solo df 
  select(-count_diffs,-count_comparisons,-count_missing,-no_sites,-no_snps) %>% #Quita cols que no se usan
  merge(ChrTable, by.x = "chromosome", by.y="ChrOrder", all=TRUE) %>% #Junta df de útiles
  arrange(factor(chromosome,levels = ChrOrder), window_pos_1, pop1, pop2) %>% #Ordena
  mutate(Comps=paste(pop1,pop2,sep="_")) %>% #Columna para agrupar por comparaciones
  mutate(WinNum= seq(1, length(chromosome))) %>% #Mete col de id de ventanas
  mutate(ChrPos= window_pos_2 - win) %>% # Para colocar el nombre del crm en el gr?fico
  group_by(chromosome,Comps) %>% # Hace que lo siguiente ocurra por cromosoma y no todo junto
  mutate(FstZscore= scale(avg_wc_fst)) %>% # Halla zscore de Fst
  mutate(DxyZscore= scale(avg_dxy)) %>% #De Dxy
  mutate(FstPvalue2sided= pnorm(-abs(FstZscore))) %>%  # Halla p-valores de Fst
  mutate(FstPcorrected= p.adjust(FstPvalue2sided, "BH")) %>% #Halla p-valor corregido de Fst
  mutate(DxyPvalue2sided= pnorm(-abs(DxyZscore))) %>% #Igual para Dxy
  mutate(DxyPcorrected= p.adjust(DxyPvalue2sided, "BH"))

#### Fst and Dxy LIN
FstDxy_lin <- dxy_lin %>%
  left_join(fst_lin) %>% #Junta Fst y Dxy en un solo df 
  select(-count_diffs,-count_comparisons,-count_missing,-no_sites,-no_snps) %>% #Quita cols que no se usan
  merge(ChrTable, by.x = "chromosome", by.y="ChrOrder", all=TRUE) %>% #Junta df de útiles
  arrange(factor(chromosome,levels = ChrOrder), window_pos_1, pop1, pop2) %>% #Ordena
  mutate(Comps=paste(pop1,pop2,sep="_")) %>% #Columna para agrupar por comparaciones
  mutate(WinNum= seq(1, length(chromosome))) %>% #Mete col de id de ventanas
  mutate(ChrPos= window_pos_2 - win) %>% # Para colocar el nombre del crm en el gr?fico
  group_by(chromosome,Comps) %>% # Hace que lo siguiente ocurra por cromosoma y no todo junto
  mutate(FstZscore= scale(avg_wc_fst)) %>% # Halla zscore de Fst
  mutate(DxyZscore= scale(avg_dxy)) %>% #De Dxy
  mutate(FstPvalue2sided= pnorm(-abs(FstZscore))) %>%  # Halla p-valores de Fst
  mutate(FstPcorrected= p.adjust(FstPvalue2sided, "BH")) %>% #Halla p-valor corregido de Fst
  mutate(DxyPvalue2sided= pnorm(-abs(DxyZscore))) %>% #Igual para Dxy
  mutate(DxyPcorrected= p.adjust(DxyPvalue2sided, "BH"))

#### Fst and Dxy POP
FstDxy_pop <- dxy_pop %>%
  left_join(fst_pop) %>% #Junta Fst y Dxy en un solo df 
  select(-count_diffs,-count_comparisons,-count_missing,-no_sites,-no_snps) %>% #Quita cols que no se usan
  merge(ChrTable, by.x = "chromosome", by.y="ChrOrder", all=TRUE) %>% #Junta df de útiles
  arrange(factor(chromosome,levels = ChrOrder), window_pos_1, pop1, pop2) %>% #Ordena
  mutate(Comps=paste(pop1,pop2,sep="_")) %>% #Columna para agrupar por comparaciones
  mutate(WinNum= seq(1, length(chromosome))) %>% #Mete col de id de ventanas
  mutate(ChrPos= window_pos_2 - win) %>% # Para colocar el nombre del crm en el gr?fico
  group_by(chromosome,Comps) %>% # Hace que lo siguiente ocurra por cromosoma y no todo junto
  mutate(FstZscore= scale(avg_wc_fst)) %>% # Halla zscore de Fst
  mutate(DxyZscore= scale(avg_dxy)) %>% #De Dxy
  mutate(FstPvalue2sided= pnorm(-abs(FstZscore))) %>%  # Halla p-valores de Fst
  mutate(FstPcorrected= p.adjust(FstPvalue2sided, "BH")) %>% #Halla p-valor corregido de Fst
  mutate(DxyPvalue2sided= pnorm(-abs(DxyZscore))) %>% #Igual para Dxy
  mutate(DxyPcorrected= p.adjust(DxyPvalue2sided, "BH"))

rmapDF <- rmap %>%
  merge(ChrTable, by.x = "CHR", by.y="ChrOrder", all=TRUE) %>%
  arrange(factor(CHR,levels = ChrOrder), POS) %>%
  mutate(SNPNum= seq(1, length(CHR))) %>%
  mutate(ChrPos= POS) %>%
  group_by(CHR)

### XPEHH
xpehh_sp_DF <- xpehh_SP %>%
  rename(chromosome = CHR.name) %>%  # Rename CHR.name to chromosome
  merge(ChrTable, by.x = "chromosome", by.y = "ChrOrder", all = TRUE) %>%
  arrange(factor(chromosome, levels = ChrOrder), POSITION, POP) %>%
  mutate(SNPNum = seq(1, n())) %>%  # Using n() instead of length(CHR.name)
  mutate(ChrPos = POSITION) %>%
  group_by(chromosome)

### Voy a crear un dataset de rmap que solo tenga las posiciones del dataset maf 0.05
filt_rmapDF <- semi_join(rmapDF, xpehh_sp_DF, by = c("ChrName", "ChrPos"))
# Remove the variables
rm(rmapDF, rmap)

combined_data <- bind_rows(FstDxy_sp, FstDxy_lin, FstDxy_pop)

POP1_sp <- "DEJU"
POP2_sp <- "YEJU"
POP1_lin <- "CANsp"
POP2_lin <- "HYEsp"
POP3_lin <- "OREsp"
POP1_pop <- "DOE"
POP2_pop <- "DOW"
POP3_pop <- "THU"
POP4_pop <- "TOW"
POP5_pop <- "CAR"
POP6_pop <- "HYE"
POP7_pop <- "PAL"
POP8_pop <- "AIK"
POP9_pop <- "MON"
POP10_pop <- "ORE"
POP11_pop <- "PIN"
POP12_pop <- "CAN"
POP13_pop <- "MEA"
POP14_pop <- "PHA"

MIX1_sp <- paste(POP1_sp,POP2_sp,sep="")
MIX1_lin <- paste(POP1_lin,POP2_lin,sep="")
MIX2_lin <- paste(POP1_lin,POP3_lin,sep="")
MIX3_lin <- paste(POP2_lin,POP3_lin,sep="")

MIX1_pop <- paste(POP1_pop,POP2_pop,sep="") # DOE DOW
MIX2_pop <- paste(POP3_pop,POP4_pop,sep="") # THU TOW
MIX3_pop <- paste(POP5_pop,POP6_pop,sep="") # CAR HYE
MIX4_pop <- paste(POP2_pop,POP7_pop,sep="") # DOW PAL
MIX5_pop <- paste(POP12_pop,POP2_pop,sep="") # CAN DOW
MIX6_pop <- paste(POP8_pop,POP6_pop,sep="") # AIK HYE
MIX7_pop <- paste(POP9_pop,POP11_pop,sep="") # MON PIN
MIX8_pop <- paste(POP9_pop,POP10_pop,sep="") # MON ORE
MIX9_pop <- paste(POP12_pop,POP13_pop,sep="" ) # CAN MEA
MIX10_pop <- paste(POP8_pop,POP13_pop,sep="" ) # MEA AIK
MIX11_pop <- paste(POP13_pop,POP9_pop,sep="" ) # MEA MON
MIX12_pop <- paste(POP9_pop,POP3_pop, sep="") #MON THU
MIX13_pop <- paste(POP7_pop,POP14_pop, sep="") # PAL PHA

CHRS <- 1
CHRE <- 31

NAME_SAVE<-paste0("/mnt/DATA/B_JUN_reseq/B_PopIndexes/H_PLOTS_2.0/C_figures/", "Rrate_fst.tiff", sep="")
tiff(NAME_SAVE, units="in", width=16, height=12, res=600)    #CHANGE THE POPS PLZ
tiff(NAME_SAVE, units="in", width=16, height=2, res=600)    #CHANGE THE POPS PLZ
dev.off()
#### OUTLIERS SPECIATION CONTINUUM ----

# outliers DEJUYEJU
outliersFst_sp_DEJUYEJU<-subset(combined_data,pop1==POP1_sp&pop2==POP2_sp&ChrNum>=CHRS&ChrNum<=CHRE&
                              FstPcorrected<0.05) # YEJU DEJU
outliersFst_sp_MIX1 <- outliersFst_sp_DEJUYEJU

# Outliers MIX1
outliersFst_lin_MIX1 <- combined_data %>%
  filter(pop1 == POP1_lin & pop2 == POP2_lin & ChrNum >= CHRS & ChrNum <= CHRE & FstPcorrected < 0.05)
# Unique outliers MIX1
unique_outliersFst_lin_MIX1 <- outliersFst_lin_MIX1 %>%
  anti_join(outliersFst_sp_DEJUYEJU, by = c("window_pos_1", "ChrName"))
# Only shared outliers MIX1 - DEJUYEJU
shared_outliersFst_sp_lin_MIX1 <- outliersFst_lin_MIX1 %>%
  anti_join(unique_outliersFst_lin_MIX1, by = c("window_pos_1", "ChrName"))

# Outliers MIX2
outliersFst_lin_MIX2 <- combined_data %>%
  filter(pop1 == POP1_lin & pop2 == POP3_lin & ChrNum >= CHRS & ChrNum <= CHRE & FstPcorrected < 0.05)
# Unique outliers CANsOREsp
unique_outliersFst_lin_MIX2 <- outliersFst_lin_MIX2 %>%
  anti_join(outliersFst_sp_DEJUYEJU, by = c("window_pos_1", "ChrName"))
# Only shared outliers CANsOREsp - DEJUYEJU
shared_outliersFst_sp_lin_MIX2 <- outliersFst_lin_MIX2 %>%
  anti_join(unique_outliersFst_lin_MIX2, by = c("window_pos_1", "ChrName"))

# Outliers MIX3
outliersFst_lin_MIX3 <- combined_data %>%
  filter(pop1 == POP2_lin & pop2 == POP3_lin & ChrNum >= CHRS & ChrNum <= CHRE & FstPcorrected < 0.05)
# Unique outliers MIX3
unique_outliersFst_lin_MIX3 <- outliersFst_lin_MIX3 %>%
  anti_join(outliersFst_sp_DEJUYEJU, by = c("window_pos_1", "ChrName"))
# Only shared outliers MIX3 - DEJUYEJU
shared_outliersFst_sp_lin_MIX3 <- outliersFst_lin_MIX3 %>%
  anti_join(unique_outliersFst_lin_MIX3, by = c("window_pos_1", "ChrName"))

# Outliers DOEDOW
outliersFst_pop_MIX1 <- combined_data %>%
  filter(pop1 == POP1_pop & pop2 == POP2_pop & ChrNum >= CHRS & ChrNum <= CHRE & FstPcorrected < 0.05)
# Unique outliers DOEDOW
unique_outliersFst_pop_MIX1 <- outliersFst_pop_MIX1 %>%
  anti_join(outliersFst_sp_MIX1, by = c("window_pos_1", "ChrName"))
# Only shared outliers DOEDOW - DEJUYEJU
shared_outliersFst_sp_pop_MIX1 <- outliersFst_pop_MIX1 %>%
  anti_join(unique_outliersFst_pop_MIX1, by = c("window_pos_1", "ChrName"))

# Outliers THUTOW
outliersFst_pop_MIX2 <- combined_data %>%
  filter(pop1 == POP3_pop & pop2 == POP4_pop & ChrNum >= CHRS & ChrNum <= CHRE & FstPcorrected < 0.05)
# Unique outliers THUTOW
unique_outliersFst_pop_MIX2 <- outliersFst_pop_MIX2 %>%
  anti_join(outliersFst_sp_MIX1, by = c("window_pos_1", "ChrName"))
# Only shared outliers THUTOW - DEJUYEJU
shared_outliersFst_sp_pop_MIX2 <- outliersFst_pop_MIX2 %>%
  anti_join(unique_outliersFst_pop_MIX2, by = c("window_pos_1", "ChrName"))
nrow(shared_outliersFst_sp_pop_MIX2)

# Outliers HYECAR
outliersFst_pop_MIX3 <- combined_data %>%
  filter(pop1 == POP5_pop & pop2 == POP6_pop & ChrNum >= CHRS & ChrNum <= CHRE & FstPcorrected < 0.05)
# Unique outliers HYECAR
unique_outliersFst_pop_MIX3 <- outliersFst_pop_MIX3 %>%
  anti_join(outliersFst_sp_MIX1, by = c("window_pos_1", "ChrName"))
# Only shared outliers HYECAR - DEJUYEJU
shared_outliersFst_sp_pop_MIX3 <- outliersFst_pop_MIX3 %>%
  anti_join(unique_outliersFst_pop_MIX3, by = c("window_pos_1", "ChrName"))

# Outliers DOWPAL
outliersFst_pop_MIX4 <- combined_data %>%
  filter(pop1 == POP2_pop & pop2 == POP7_pop & ChrNum >= CHRS & ChrNum <= CHRE & FstPcorrected < 0.05)
# Unique outliers DOWPAL
unique_outliersFst_pop_MIX4 <- outliersFst_pop_MIX4 %>%
  anti_join(outliersFst_sp_MIX1, by = c("window_pos_1", "ChrName"))
# Only shared outliers DOWPAL - DEJUYEJU
shared_outliersFst_sp_pop_MIX4 <- outliersFst_pop_MIX4 %>%
  anti_join(unique_outliersFst_pop_MIX4, by = c("window_pos_1", "ChrName"))

# Outliers CANDOW
outliersFst_pop_MIX5 <- combined_data %>%
  filter(pop1 == POP12_pop & pop2 == POP2_pop & ChrNum >= CHRS & ChrNum <= CHRE & FstPcorrected < 0.05)
# Unique outliers CANDOW
unique_outliersFst_pop_MIX5 <- outliersFst_pop_MIX5 %>%
  anti_join(outliersFst_sp_MIX1, by = c("window_pos_1", "ChrName"))
# Only shared outliers CANDOW - DEJUYEJU
shared_outliersFst_sp_pop_MIX5 <- outliersFst_pop_MIX5 %>%
  anti_join(unique_outliersFst_pop_MIX5, by = c("window_pos_1", "ChrName"))

# Outliers AIKHYE
outliersFst_pop_MIX6 <- combined_data %>%
  filter(pop1 == POP8_pop & pop2 == POP6_pop & ChrNum >= CHRS & ChrNum <= CHRE & FstPcorrected < 0.05)
# Unique outliers CANDOW
unique_outliersFst_pop_MIX6 <- outliersFst_pop_MIX6 %>%
  anti_join(outliersFst_sp_MIX1, by = c("window_pos_1", "ChrName"))
# Only shared outliers CANDOW - DEJUYEJU
shared_outliersFst_sp_pop_MIX6 <- outliersFst_pop_MIX6 %>%
  anti_join(unique_outliersFst_pop_MIX6, by = c("window_pos_1", "ChrName"))

# Outliers MON PIN 
outliersFst_pop_MIX7 <- combined_data %>%
  filter(pop1 == POP9_pop & pop2 == POP11_pop & ChrNum >= CHRS & ChrNum <= CHRE & FstPcorrected < 0.05)
# Unique outliers CANDOW
unique_outliersFst_pop_MIX7 <- outliersFst_pop_MIX7 %>%
  anti_join(outliersFst_sp_MIX1, by = c("window_pos_1", "ChrName"))
# Only shared outliers CANDOW - DEJUYEJU
shared_outliersFst_sp_pop_MIX7 <- outliersFst_pop_MIX7 %>%
  anti_join(unique_outliersFst_pop_MIX7, by = c("window_pos_1", "ChrName"))

# Outliers  ORE MON
outliersFst_pop_MIX8 <- combined_data %>%
  filter(pop1 == POP9_pop & pop2 == POP10_pop & ChrNum >= CHRS & ChrNum <= CHRE & FstPcorrected < 0.05)
# Unique outliers CANDOW
unique_outliersFst_pop_MIX8 <- outliersFst_pop_MIX8 %>%
  anti_join(outliersFst_sp_MIX1, by = c("window_pos_1", "ChrName"))
# Only shared outliers CANDOW - DEJUYEJU
shared_outliersFst_sp_pop_MIX8 <- outliersFst_pop_MIX8 %>%
  anti_join(unique_outliersFst_pop_MIX8, by = c("window_pos_1", "ChrName"))

# Outliers  CAN MEA
outliersFst_pop_MIX9 <- combined_data %>%
  filter(pop1 == POP12_pop & pop2 == POP13_pop & ChrNum >= CHRS & ChrNum <= CHRE & FstPcorrected < 0.05)
# Unique outliers CANDOW
unique_outliersFst_pop_MIX9 <- outliersFst_pop_MIX9 %>%
  anti_join(outliersFst_sp_MIX1, by = c("window_pos_1", "ChrName"))
# Only shared outliers CANDOW - DEJUYEJU
shared_outliersFst_sp_pop_MIX9 <- outliersFst_pop_MIX9 %>%
  anti_join(unique_outliersFst_pop_MIX9, by = c("window_pos_1", "ChrName"))

# Outliers  AIK MEA
outliersFst_pop_MIX10 <- combined_data %>%
  filter(pop1 == POP8_pop & pop2 == POP13_pop & ChrNum >= CHRS & ChrNum <= CHRE & FstPcorrected < 0.05)
# Unique outliers CANDOW
unique_outliersFst_pop_MIX10 <- outliersFst_pop_MIX10 %>%
  anti_join(outliersFst_sp_MIX1, by = c("window_pos_1", "ChrName"))
# Only shared outliers CANDOW - DEJUYEJU
shared_outliersFst_sp_pop_MIX10 <- outliersFst_pop_MIX10 %>%
  anti_join(unique_outliersFst_pop_MIX10, by = c("window_pos_1", "ChrName"))

# Outliers  MEAMON
outliersFst_pop_MIX11 <- combined_data %>%
  filter(pop1 == POP13_pop & pop2 == POP9_pop & ChrNum >= CHRS & ChrNum <= CHRE & FstPcorrected < 0.05)
# Unique outliers CANDOW
unique_outliersFst_pop_MIX11 <- outliersFst_pop_MIX11 %>%
  anti_join(outliersFst_sp_MIX1, by = c("window_pos_1", "ChrName"))
# Only shared outliers CANDOW - DEJUYEJU
shared_outliersFst_sp_pop_MIX11 <- outliersFst_pop_MIX11 %>%
  anti_join(unique_outliersFst_pop_MIX11, by = c("window_pos_1", "ChrName"))

############################################# Overlapping outliers for the venn diagram - LINEAGES

outliersFst_lin_MIX1_2 <- outliersFst_lin_MIX1 %>%
  semi_join(outliersFst_lin_MIX2, by = c("window_pos_1", "ChrName"))
outliersFst_lin_MIX2_1 <- outliersFst_lin_MIX2 %>%
  semi_join(outliersFst_lin_MIX1, by = c("window_pos_1", "ChrName"))
outliersFst_lin_MIX1_3 <- outliersFst_lin_MIX1 %>%
  semi_join(outliersFst_lin_MIX3, by = c("window_pos_1", "ChrName"))
outliersFst_lin_MIX3_1 <- outliersFst_lin_MIX3 %>%
  semi_join(outliersFst_lin_MIX1, by = c("window_pos_1", "ChrName"))
outliersFst_lin_MIX2_3 <- outliersFst_lin_MIX2 %>%
  semi_join(outliersFst_lin_MIX3, by = c("window_pos_1", "ChrName"))
outliersFst_lin_MIX3_2 <- outliersFst_lin_MIX3 %>%
  semi_join(outliersFst_lin_MIX2, by = c("window_pos_1", "ChrName"))
outliersFst_sp_lin1_2 <- outliersFst_lin_MIX1_2 %>%
  semi_join(outliersFst_sp_MIX1, by = c("window_pos_1", "ChrName"))
outliersFst_sp_lin1_3 <- outliersFst_lin_MIX1_3 %>%
  semi_join(outliersFst_sp_MIX1, by = c("window_pos_1", "ChrName"))
outliersFst_sp_lin2_3 <- outliersFst_lin_MIX2_3 %>%
  semi_join(outliersFst_sp_MIX1, by = c("window_pos_1", "ChrName"))
outliersFst_lin1_2_3 <- outliersFst_lin_MIX1_2 %>%
  semi_join(outliersFst_lin_MIX3, by = c("window_pos_1", "ChrName"))
outliersFst_lin3_1_2 <- outliersFst_lin_MIX2_3 %>%
  semi_join(outliersFst_lin_MIX1, by = c("window_pos_1", "ChrName"))
outliersFst_lin2_3_1 <- outliersFst_lin_MIX3_1 %>%
  semi_join(outliersFst_lin_MIX1, by = c("window_pos_1", "ChrName"))
outliersFst_sp_lin1_2_3 <- outliersFst_lin1_2_3 %>%
  semi_join(outliersFst_sp_MIX1, by = c("window_pos_1", "ChrName"))

############################################# Overlapping outliers for the venn diagram - POPULATIONS

outliersFst_pop_MIX1_2 <- outliersFst_pop_MIX1 %>%
  semi_join(outliersFst_pop_MIX2, by = c("window_pos_1", "ChrName"))
outliersFst_pop_MIX1_3 <- outliersFst_pop_MIX1 %>%
  semi_join(outliersFst_pop_MIX3, by = c("window_pos_1", "ChrName"))
outliersFst_pop_MIX2_3 <- outliersFst_pop_MIX2 %>%
  semi_join(outliersFst_pop_MIX3, by = c("window_pos_1", "ChrName"))
outliersFst_sp_pop1_2 <- outliersFst_pop_MIX1_2 %>%
  semi_join(outliersFst_sp_MIX1, by = c("window_pos_1", "ChrName"))
outliersFst_sp_pop1_3 <- outliersFst_pop_MIX1_3 %>%
  semi_join(outliersFst_sp_MIX1, by = c("window_pos_1", "ChrName"))
outliersFst_sp_pop2_3 <- outliersFst_pop_MIX2_3 %>%
  semi_join(outliersFst_sp_MIX1, by = c("window_pos_1", "ChrName"))
outliersFst_pop1_2_3 <- outliersFst_pop_MIX1_2 %>%
  semi_join(outliersFst_pop_MIX3, by = c("window_pos_1", "ChrName"))
outliersFst_sp_pop1_2_3 <- outliersFst_pop1_2_3 %>%
  semi_join(outliersFst_sp_MIX1, by = c("window_pos_1", "ChrName"))


##### MANHATTAN PLOT ----
manhattanJAVI <- function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("gray10", 
                                                                                    "gray60"), chrlabs = NULL, suggestiveline = -log10(1e-05), 
                           genomewideline = -log10(5e-08), highlight1 = NULL, highlight2 = NULL, 
                           highlight3 = NULL, highlight4 = NULL, highlight5 = NULL, highlight6 = NULL, highlight7 = NULL, logp = TRUE, 
                           annotatePval = NULL, annotateTop = TRUE, space = 5e6, ...) 
{
  CHR = BP = P = index = NULL
  if (!(chr %in% names(x))) 
    stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x))) 
    stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x))) 
    stop(paste("Column", p, "not found!"))
  if (!(snp %in% names(x))) 
    warning(paste("No SNP column found. OK unless you're trying to highlight."))
  if (!is.numeric(x[[chr]])) 
    stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]])) 
    stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]])) 
    stop(paste(p, "column should be numeric."))
  if (!is.null(x[[snp]])) 
    d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]], 
                   pos = NA, index = NA, SNP = x[[snp]], stringsAsFactors = FALSE)
  else d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]], 
                      pos = NA, index = NA)
  d <- d[order(d$CHR, d$BP), ]
  if (logp) {
    d$logp <- -log10(d$P)
  }
  else {
    d$logp <- d$P
  }
  d$index = rep.int(seq_along(unique(d$CHR)), times = tapply(d$SNP, 
                                                             d$CHR, length))
  nchr = length(unique(d$CHR))
  if (nchr == 1) {
    d$pos = d$BP
    xlabel = paste("Chromosome", unique(d$CHR), "position")
  }
  else {
    lastbase = 0
    ticks = NULL
    for (i in unique(d$index)) {
      if (i == 1) {
        d[d$index == i, ]$pos = d[d$index == i, ]$BP
      }
      else {
        lastbase = lastbase + max(d[d$index == (i - 1), "BP"]) + space
        d[d$index == i, "BP"] = d[d$index == i, "BP"] - 
          min(d[d$index == i, "BP"]) + 1
        d[d$index == i, "pos"] = d[d$index == i, "BP"] + 
          lastbase
      }
    }
    ticks <- tapply(d$pos, d$index, quantile, probs = 0.5)
    xlabel = "Chromosome"
    labs <- unique(d$CHR)
  }
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
                   las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0, 
                                                                     ceiling(max(d$logp))), xlab = xlabel, ylab = expression(-log[10](italic(p))))
  dotargs <- list(...)
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% 
                                            names(dotargs)]))
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs) == length(labs)) {
        labs <- chrlabs
      }
      else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    }
    else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }
  if (nchr == 1) {
    axis(1, ...)
  }
  else {
    axis(1, at = ticks, labels = labs, ...)
  }
  col = rep_len(col, max(d$index))
  if (nchr == 1) {
    with(d, points(pos, logp, pch = 20, col = col[1], ...))
  }
  else {
    icol = 1
    for (i in unique(d$index)) {
      points(d[d$index == i, "pos"], d[d$index == i, "logp"], 
             col = col[icol], pch = 20, ...)
      icol = icol + 1
    }
  }
  if (suggestiveline) 
    abline(h = suggestiveline, col = "blue")
  if (genomewideline) 
    abline(h = genomewideline, col = "red")
  if (!is.null(highlight1)) {
    if (any(!(highlight1 %in% d$SNP))) 
      warning("You're trying to highlight1 SNPs that don't exist in your results.")
    d.highlight1 = d[which(d$SNP %in% highlight1), ]
    with(d.highlight1, points(pos, logp, col = "green",  # #F9C825
                              pch = 20, ...))
    
  }
  if (!is.null(highlight2)) {
    if (any(!(highlight2 %in% d$SNP))) 
      warning("You're trying to highlight2 SNPs that don't exist in your results.")
    d.highlight2 = d[which(d$SNP %in% highlight2), ]
    with(d.highlight2, points(pos, logp, col = "purple", #"#009E73"
                              pch = 20, ...))
    
  }
  if (!is.null(highlight3)) {
    if (any(!(highlight3 %in% d$SNP))) 
      warning("You're trying to highlight3 SNPs that don't exist in your results.")
    d.highlight3 = d[which(d$SNP %in% highlight3), ]
    with(d.highlight3, points(pos, logp, col = "red",  #BC79A7
                              pch = 20, ...))
    
  }
  if (!is.null(highlight4)) {
    if (any(!(highlight4 %in% d$SNP))) 
      warning("You're trying to highlight4 SNPs that don't exist in your results.")
    d.highlight4 = d[which(d$SNP %in% highlight4), ]
    with(d.highlight4, points(pos, logp, col = "blue", # #56B4E9
                              pch = 20, ...))
  }
  if (!is.null(highlight5)) {
    if (any(!(highlight5 %in% d$SNP))) 
      warning("You're trying to highlight5 SNPs that don't exist in your results.")
    d.highlight5 = d[which(d$SNP %in% highlight5), ]
    with(d.highlight5, points(pos, logp, col = "red", # EF005E
                              pch = 20, ...))
    
  }
  if (!is.null(annotatePval)) {
    if (logp) {
      topHits = subset(d, P <= annotatePval)
    }
    else topHits = subset(d, P >= annotatePval)
    par(xpd = TRUE)
    if (annotateTop == FALSE) {
      if (logp) {
        with(subset(d, P <= annotatePval), textxy(pos, 
                                                  -log10(P), offset = 0.625, labs = topHits$SNP, 
                                                  cex = 0.45), ...)
      }
      else with(subset(d, P >= annotatePval), textxy(pos, 
                                                     P, offset = 0.625, labs = topHits$SNP, cex = 0.45), 
                ...)
    }
    else {
      topHits <- topHits[order(topHits$P), ]
      topSNPs <- NULL
      for (i in unique(topHits$CHR)) {
        chrSNPs <- topHits[topHits$CHR == i, ]
        topSNPs <- rbind(topSNPs, chrSNPs[1, ])
      }
      if (logp) {
        textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.625, 
               labs = topSNPs$SNP, cex = 0.5, ...)
      }
      else textxy(topSNPs$pos, topSNPs$P, offset = 0.625, 
                  labs = topSNPs$SNP, cex = 0.5, ...)
    }
  }
  par(xpd = FALSE)
}


# Reset to default
par()

n=3

# Set graphical parameters
par(mfrow=c(n,1), oma=c(2,0,0,0), mar=c(1,4,1,1))

manhattanJAVI(subset(FstDxy_sp,pop1==POP1_sp & pop2==POP2_sp & ChrNum>=CHRS & ChrNum<=CHRE),#Fst SP
              chr="ChrNum", bp="ChrPos", p="FstZscore",logp="FALSE", snp="WinNum",
              suggestiveline=FALSE,genomewideline=FALSE,
              highlight1 = outliersFst_sp_MIX1$WinNum,
              ylim = c(-2,20),col = c("grey40", "grey"), cex=1,
              ylab ="Fst Z-score\nDEJU YEJU", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI(subset(FstDxy_lin,pop1==POP1_lin & pop2==POP2_lin & ChrNum>=CHRS & ChrNum<=CHRE),#Fst LIN1-2 (MIX1)
              chr="ChrNum", bp="ChrPos", p="FstZscore",logp="FALSE", snp="WinNum",
              suggestiveline=FALSE,genomewideline=FALSE,
              highlight4 = c(outliersFst_lin1_2_3$WinNum, outliersFst_lin_MIX1_2$WinNum, outliersFst_lin_MIX1_3$WinNum),
              highlight3 = outliersFst_lin_MIX1$WinNum,
              ylim = c(-2,20),col = c("grey40", "grey"), cex=1,
              ylab ="Fst Z-score\nCANsp HYEsp", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI(subset(FstDxy_lin,pop1==POP1_lin & pop2==POP3_lin & ChrNum>=CHRS & ChrNum<=CHRE),#Fst LIN 1-3 (MIX 2)
              chr="ChrNum", bp="ChrPos", p="FstZscore",logp="FALSE", snp="WinNum",
              suggestiveline=FALSE,genomewideline=FALSE,
              highlight4 = c(outliersFst_lin3_1_2$WinNum, outliersFst_lin_MIX2_1$WinNum, outliersFst_lin_MIX2_3$WinNum),
              highlight3 = outliersFst_lin_MIX2$WinNum,
              ylim = c(-2,20),col = c("grey40", "grey"), cex=1,
              ylab ="Fst Z-score\nCANsp OREsp", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI(subset(FstDxy_lin,pop1==POP2_lin & pop2==POP3_lin & ChrNum>=CHRS & ChrNum<=CHRE),#Fst LIN 2-3 (MIX 3)
              chr="ChrNum", bp="ChrPos", p="FstZscore",logp="FALSE", snp="WinNum",
              suggestiveline=FALSE,genomewideline=FALSE,
              highlight4 = c(outliersFst_lin2_3_1$WinNum, outliersFst_lin_MIX3_1$WinNum, outliersFst_lin_MIX3_2$WinNum),
              highlight3 = outliersFst_lin_MIX3$WinNum,
              ylim = c(-2,20),col = c("grey40", "grey"), cex=1,
              ylab ="Fst Z-score\nHYEsp OREsp", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI(subset(FstDxy_pop,pop1==POP2_pop & pop2==POP7_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst DOWPAL
              chr="ChrNum", bp="ChrPos", p="FstZscore",logp="FALSE", snp="WinNum",
              suggestiveline=FALSE,genomewideline=FALSE,
              highlight1 = shared_outliersFst_sp_pop_MIX4$WinNum,
              highlight4 = unique_outliersFst_pop_MIX4$WinNum,
              ylim = c(-2,20),col = c("grey40", "grey"), cex=1,
              ylab ="Fst Z-score\nDOW - PAL", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI(subset(FstDxy_pop,pop1==POP1_pop & pop2==POP2_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst DOEDOW
              chr="ChrNum", bp="ChrPos", p="FstZscore",logp="FALSE", snp="WinNum",
              suggestiveline=FALSE,genomewideline=FALSE,
              highlight1 = shared_outliersFst_sp_pop_MIX1$WinNum,
              highlight4 = unique_outliersFst_pop_MIX1$WinNum,
              ylim = c(-2,20),col = c("grey40", "grey"), cex=1,
              ylab ="Fst Z-score\nDOE - DOW", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI(subset(FstDxy_pop,pop1==POP12_pop & pop2==POP2_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst CAN DOW
              chr="ChrNum", bp="ChrPos", p="FstZscore",logp="FALSE", snp="WinNum",
              suggestiveline=FALSE,genomewideline=FALSE,
              highlight1 = shared_outliersFst_sp_pop_MIX5$WinNum,
              highlight4 = unique_outliersFst_pop_MIX5$WinNum,
              ylim = c(-2,20),col = c("grey40", "grey"), cex=1,
              ylab ="Fst Z-score\nDOW - CAN", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI(subset(FstDxy_pop,pop1==POP12_pop & pop2==POP13_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst CANMEA
              chr="ChrNum", bp="ChrPos", p="FstZscore",logp="FALSE", snp="WinNum",
              suggestiveline=FALSE,genomewideline=FALSE,
              highlight1 = shared_outliersFst_sp_pop_MIX9$WinNum,
              highlight4 = unique_outliersFst_pop_MIX9$WinNum,
              ylim = c(-2,20),col = c("grey40", "grey"), cex=1,
              ylab ="Fst Z-score\nCAN - MEA", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI(subset(FstDxy_pop,pop1==POP13_pop & pop2==POP9_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst MEA MON
              chr="ChrNum", bp="ChrPos", p="FstZscore",logp="FALSE", snp="WinNum",
              suggestiveline=FALSE,genomewideline=FALSE,
              highlight1 = shared_outliersFst_sp_pop_MIX11$WinNum,
              highlight4 = unique_outliersFst_pop_MIX11$WinNum,
              ylim = c(-2,20),col = c("grey40", "grey"), cex=1,
              ylab ="Fst Z-score\nMEA - MON", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI(subset(FstDxy_pop,pop1==POP9_pop & pop2==POP11_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst MONPIN
              chr="ChrNum", bp="ChrPos", p="FstZscore",logp="FALSE", snp="WinNum",
              suggestiveline=FALSE,genomewideline=FALSE,
              highlight1 = shared_outliersFst_sp_pop_MIX7$WinNum,
              highlight4 = unique_outliersFst_pop_MIX7$WinNum,
              ylim = c(-2,20),col = c("grey40", "grey"), cex=1,
              ylab ="Fst Z-score\nCAN - DOW", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI(subset(FstDxy_pop,pop1==POP9_pop & pop2==POP10_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst MONORE
              chr="ChrNum", bp="ChrPos", p="FstZscore",logp="FALSE", snp="WinNum",
              suggestiveline=FALSE,genomewideline=FALSE,
              highlight1 = shared_outliersFst_sp_pop_MIX8$WinNum,
              highlight4 = unique_outliersFst_pop_MIX8$WinNum,
              ylim = c(-2,20),col = c("grey40", "grey"), cex=1,
              ylab ="Fst Z-score\nMON - ORE", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI(subset(FstDxy_pop,pop1==POP3_pop & pop2==POP4_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst THU TOW
              chr="ChrNum", bp="ChrPos", p="FstZscore",logp="FALSE", snp="WinNum",
              suggestiveline=FALSE,genomewideline=FALSE,
              highlight1 = shared_outliersFst_sp_pop_MIX2$WinNum,
              highlight4 = unique_outliersFst_pop_MIX2$WinNum,
              ylim = c(-2,20),col = c("grey40", "grey"), cex=1,
              ylab ="Fst Z-score\nTHU - TOW", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI(subset(FstDxy_pop,pop1==POP8_pop & pop2==POP13_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst AIK MEA
              chr="ChrNum", bp="ChrPos", p="FstZscore",logp="FALSE", snp="WinNum",
              suggestiveline=FALSE,genomewideline=FALSE,
              highlight1 = shared_outliersFst_sp_pop_MIX10$WinNum,
              highlight4 = unique_outliersFst_pop_MIX10$WinNum,
              ylim = c(-2,20),col = c("grey40", "grey"), cex=1,
              ylab ="Fst Z-score\nAIK - MEA", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI(subset(FstDxy_pop,pop1==POP8_pop & pop2==POP6_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst AIKHYE
              chr="ChrNum", bp="ChrPos", p="FstZscore",logp="FALSE", snp="WinNum",
              suggestiveline=FALSE,genomewideline=FALSE,
              highlight1 = shared_outliersFst_sp_pop_MIX6$WinNum,
              highlight4 = unique_outliersFst_pop_MIX6$WinNum,
              ylim = c(-2,20),col = c("grey40", "grey"), cex=1,
              ylab ="Fst Z-score\nAIK - HYE", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI(subset(FstDxy_pop,pop1==POP5_pop & pop2==POP6_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst CAR HYE
              chr="ChrNum", bp="ChrPos", p="FstZscore",logp="FALSE", snp="WinNum",
              suggestiveline=FALSE,genomewideline=FALSE,
              highlight1 = shared_outliersFst_sp_pop_MIX3$WinNum,
              highlight4 = unique_outliersFst_pop_MIX3$WinNum,
              ylim = c(-2,20),col = c("grey40", "grey"), cex=1,
              ylab ="Fst Z-score\nHYE - CAR", chrlabs = ChrName[CHRS:CHRE])


manhattanJAVI(subset(filt_rmapDF, ChrNum>=CHRS & ChrNum<=CHRE),#Rmap
              chr="ChrNum", bp="ChrPos", p="Mean_rho",logp="FALSE", snp="SNPNum",
              suggestiveline=FALSE,genomewideline=FALSE,
              #highlight1 = shared_outliersFst_sp_pop_MIX3$WinNum,
              #highlight3 = unique_outliersFst_pop_MIX3$WinNum,
              ylim = c(-5,300),col = c("grey40", "grey"), cex=0.3,
              ylab ="Rho", chrlabs = ChrName[CHRS:CHRE])

dev.off()
beep()


##### Reference four-set diagram ----

venn.plot_LIN <- draw.quad.venn(
  area1 = nrow(outliersFst_sp_MIX1),
  area2 = nrow(outliersFst_lin_MIX1),
  area3 = nrow(outliersFst_lin_MIX2),
  area4 = nrow(outliersFst_lin_MIX3),
  n12 = nrow(shared_outliersFst_sp_lin_MIX1),
  n13 = nrow(shared_outliersFst_sp_lin_MIX2),
  n14 = nrow(shared_outliersFst_sp_lin_MIX3),
  n23 = nrow(outliersFst_lin_MIX1_2),
  n24 = nrow(outliersFst_lin_MIX1_3),
  n34 = nrow(outliersFst_lin_MIX2_3),
  n123 = nrow(outliersFst_sp_lin1_2),
  n124 = nrow(outliersFst_sp_lin1_3),
  n134 = nrow(outliersFst_sp_lin2_3),
  n234 = nrow(outliersFst_lin1_2_3),
  n1234 = nrow(outliersFst_sp_lin1_2_3),
  category = c("DEJU - YEJU", "CANsp - HYEsp", "CANsp - OREsp", "HYEsp - OREsp"),
  fill = c("gold", "#CC79A7", "#009E73", "#56B4E9"),
  lty = 1,
  lwd = 0.7, 
  cex = 2,
  alpha = 0.5,
  cat.cex = 1,
  cat.col = "black"
)

# Reference four-set diagram
venn.plot_LIN <- draw.quad.venn(
  area1 = nrow(outliersFst_sp_MIX1),
  area2 = nrow(outliersFst_pop_MIX1),
  area3 = nrow(outliersFst_pop_MIX2),
  area4 = nrow(outliersFst_pop_MIX3),
  n12 = nrow(shared_outliersFst_sp_pop_MIX1),
  n13 = nrow(shared_outliersFst_sp_pop_MIX2),
  n14 = nrow(shared_outliersFst_sp_pop_MIX3),
  n23 = nrow(outliersFst_pop_MIX1_2),
  n24 = nrow(outliersFst_pop_MIX1_3),
  n34 = nrow(outliersFst_pop_MIX2_3),
  n123 = nrow(outliersFst_sp_pop1_2),
  n124 = nrow(outliersFst_sp_pop1_3),
  n134 = nrow(outliersFst_sp_pop2_3),
  n234 = nrow(outliersFst_pop1_2_3),
  n1234 = nrow(outliersFst_sp_pop1_2_3),
  category = c("DEJU - YEJU", "DOW - DOE", "THU - TOW", "CAR - HYE"),
  fill = c("gold", "#CC79A7", "#009E73", "#56B4E9"),
  lty = 1,
  lwd = 0.7, 
  cex = 2,
  alpha = 0.5,
  cat.cex = 1,
  cat.col = "black"
)

##### PLOT recomb rate vs Fst species ----
# Define the function for rounding down to the nearest multiple
round_down_to_nearest <- function(x, base) {
  return(base * floor(x / base) + 1)
}

### AUTOSOMES

# Initialize an empty data frame to store the results
rmapDF_fst <- data.frame()
# Loop over each ChrNum value from 1 to 31
for (i in 1:30) {
  # Subset data for the current ChrNum value
  A <- subset(FstDxy_sp, ChrNum == i)
  B <- subset(filt_rmapDF, ChrNum == i)
  
  # Apply rounding to the nearest multiple of 50000 for window_pos_1
  B$window_pos_1 <- round_down_to_nearest(B$POS, 50000)
  
  # Perform the left join
  merged_AB <- dplyr::left_join(B, A[, c("window_pos_1", "avg_wc_fst")], by = "window_pos_1")
  
  # Concatenate the result into the final data frame
  rmapDF_fst <- rbind(rmapDF_fst, merged_AB)
}

avg_mean_rho_df <- data.frame()
avg_mean_rho_df <- rmapDF_fst %>%
  group_by(ChrNum, ChrName, window_pos_1) %>%
  summarise(avg_mean_rho = mean(Mean_rho, na.rm = TRUE)) %>%
  ungroup()

# Perform the left join to include avg_mean_rho in FstDxy_sp
FstDxy_sp_rrate <- dplyr::left_join(FstDxy_sp, avg_mean_rho_df, by = c("ChrNum", "ChrName", "window_pos_1"))
# Remove rows with NA in avg_mean_rho
FstDxy_sp_rrate <- FstDxy_sp_rrate %>% filter(!is.na(avg_mean_rho))

# Calculate the threshold for the top 0.5% FST
fst_threshold <- quantile(FstDxy_sp_rrate$avg_wc_fst, 0.995, na.rm = TRUE)

# Subset the data for regions with elevated differentiation and moderate to high recombination rate
high_diff_regions_AUTO <- subset(FstDxy_sp_rrate, avg_wc_fst >= fst_threshold & avg_mean_rho > 25)

# Subset the data for regions with elevated differentiation and low recombination rate
fake_high_diff_regions_AUTO <- subset(FstDxy_sp_rrate, avg_wc_fst >= fst_threshold & avg_mean_rho <= 25)

# 3rd step: scatterplot

# Create the scatter plot
scatter_plot <- ggplot(FstDxy_sp_rrate, aes(x = avg_mean_rho, y = avg_wc_fst)) +
  geom_point(alpha = 0.5, size = 2, stroke = 0, color = "gray60") +
  geom_point(data = fake_high_diff_regions_AUTO, aes(x = avg_mean_rho, y = avg_wc_fst), color = "red", alpha = 1, size = 3, stroke = 0) +
  geom_point(data = high_diff_regions_AUTO, aes(x = avg_mean_rho, y = avg_wc_fst), color = "green4", alpha = 1, size = 4, stroke = 0) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  labs(title = "DEJU vs YEJU Autosomes", x = NULL, y = NULL) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),         # No major gridlines
    panel.grid.minor = element_blank(),         # No minor gridlines
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white")) +
  stat_cor(method = "pearson", 
           aes(label = after_stat(r)), 
           label.x = max(FstDxy_sp_rrate$avg_mean_rho, na.rm = TRUE) * 0.6, 
           label.y = max(FstDxy_sp_rrate$avg_wc_fst, na.rm = TRUE) * 0.9, 
           size = 9, 
           label.sep = "")

# Save the plot as a TIFF file
ggsave("/mnt/DATA/B_JUN_reseq/B_PopIndexes/H_PLOTS_2.0/F_rrate_vs_fst/scatter_plot_DEJUYEJU_auto.tiff",  plot = scatter_plot, width = 3, height = 4, dpi = 280, compression = "lzw", bg ="white")

### SEX Z CHROMOSOME

# Initialize an empty data frame to store the results
rmapDF_fst <- data.frame()
# Loop over each ChrNum value from 1 to 31
for (i in 31:31) {
  # Subset data for the current ChrNum value
  A <- subset(FstDxy_sp, ChrNum == i)
  B <- subset(filt_rmapDF, ChrNum == i)
  
  # Apply rounding to the nearest multiple of 50000 for window_pos_1
  B$window_pos_1 <- round_down_to_nearest(B$POS, 50000)
  
  # Perform the left join
  merged_AB <- dplyr::left_join(B, A[, c("window_pos_1", "avg_wc_fst")], by = "window_pos_1")
  
  # Concatenate the result into the final data frame
  rmapDF_fst <- rbind(rmapDF_fst, merged_AB)
}

avg_mean_rho_df <- data.frame()
avg_mean_rho_df <- rmapDF_fst %>%
  group_by(ChrNum, ChrName, window_pos_1) %>%
  summarise(avg_mean_rho = mean(Mean_rho, na.rm = TRUE)) %>%
  ungroup()

# Perform the left join to include avg_mean_rho in FstDxy_sp
FstDxy_sp_rrate <- dplyr::left_join(FstDxy_sp, avg_mean_rho_df, by = c("ChrNum", "ChrName", "window_pos_1"))
# Remove rows with NA in avg_mean_rho
FstDxy_sp_rrate <- FstDxy_sp_rrate %>% filter(!is.na(avg_mean_rho))

# Calculate the threshold for the top 0.5% FST
fst_threshold <- quantile(FstDxy_sp_rrate$avg_wc_fst, 0.995, na.rm = TRUE)

# Subset the data for regions with elevated differentiation and moderate to high recombination rate &  elevated differentiation and low recombination rate
high_diff_regions_Zchr <- subset(FstDxy_sp_rrate, avg_wc_fst >= fst_threshold & avg_mean_rho > 25)
fake_high_diff_regions_Zchr <- subset(FstDxy_sp_rrate, avg_wc_fst >= fst_threshold & avg_mean_rho <= 25)

# 3rd step: scatterplot

# Create the scatter plot
scatter_plot <- ggplot(FstDxy_sp_rrate, aes(x = avg_mean_rho, y = avg_wc_fst)) +
  geom_point(alpha = 0.5, size = 2, stroke = 0, color = "gray60") +
  geom_point(data = fake_high_diff_regions_Zchr, aes(x = avg_mean_rho, y = avg_wc_fst), color = "red", alpha = 1, size = 3, stroke = 0) +
  geom_point(data = high_diff_regions_Zchr, aes(x = avg_mean_rho, y = avg_wc_fst), color = "green4", alpha = 1, size = 4, stroke = 0) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  labs(title = "DEJU vs YEJU Z-Chr", x = NULL, y = NULL) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),         # No major gridlines
    panel.grid.minor = element_blank(),         # No minor gridlines
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white")) +
  stat_cor(method = "pearson", 
           aes(label = after_stat(r)), 
           label.x = max(FstDxy_sp_rrate$avg_mean_rho, na.rm = TRUE) * 0.6, 
           label.y = max(FstDxy_sp_rrate$avg_wc_fst, na.rm = TRUE) * 0.9, 
           size = 9, 
           label.sep = "")

# Save the plot as a TIFF file
ggsave("/mnt/DATA/B_JUN_reseq/B_PopIndexes/H_PLOTS_2.0/F_rrate_vs_fst/scatter_plot_DEJUYEJU_Zchr.tiff",  plot = scatter_plot, width = 3, height = 4, dpi = 280, compression = "lzw", bg ="white")

high_diff_regions_SP <- dplyr::bind_rows(high_diff_regions_Zchr, high_diff_regions_AUTO)
fake_high_diff_regions_SP <- dplyr::bind_rows(fake_high_diff_regions_Zchr, fake_high_diff_regions_AUTO)


##### PLOT recomb rate vs Fst lineages ----

# Initialize an empty data frame to store the results
FstDxy_lin <- FstDxy_lin %>%
  unite(combined_pop, pop1, pop2, sep = "_", remove = FALSE)

lineages <- c("CANsp_HYEsp", "CANsp_OREsp", "HYEsp_OREsp" )

### AUTOSOMES
# Loop over each ChrNum value from 1 to 30
for (j in lineages) {
  assign(paste0("rmapDF_fst_", j), data.frame())
  LIN <- subset(FstDxy_lin,combined_pop == j )
  for (i in 1:30) {
    # Subset data for the current ChrNum value
    A <- subset(LIN, ChrNum == i)
    B <- subset(filt_rmapDF, ChrNum == i)
    
    # Apply rounding to the nearest multiple of 50000 for window_pos_1
    B$window_pos_1 <- round_down_to_nearest(B$POS, 50000)
    
    # Perform the left join
    merged_AB <- dplyr::left_join(B, A[, c("window_pos_1", "avg_wc_fst")], by = "window_pos_1")
    
    # Concatenate the result into the final data frame
    assign(paste0("rmapDF_fst_", j), rbind(get(paste0("rmapDF_fst_", j)), merged_AB))
  }
}
beepr::beep()

for (x in lineages) {
  # Calculate avg_mean_rho_df
  avg_mean_rho_df <- get(paste0("rmapDF_fst_", x)) %>%
    group_by(ChrNum, ChrName, window_pos_1) %>%
    summarise(avg_mean_rho = mean(Mean_rho, na.rm = TRUE), .groups = "drop")  # Ensure to drop grouping
  
  # Perform left join
  FstDxy_lin_rrate <- dplyr::left_join(FstDxy_lin, avg_mean_rho_df, by = c("ChrNum", "ChrName", "window_pos_1"))
  
  # Filter FstDxy_lin_rrate by combined_pop == x and remove NA values in avg_mean_rho
  assign(paste0("FstDxy_lin_rrate_", x),
         FstDxy_lin_rrate %>%
           filter(!is.na(avg_mean_rho) & combined_pop == x))
  
  # Calculate the threshold for the top 0.5% FST
  fst_threshold <- quantile(get(paste0("FstDxy_lin_rrate_", x))$avg_wc_fst, 0.995, na.rm = TRUE)
  
  # Subset the data for regions with elevated differentiation and moderate to high recombination rate & elevated differentiation and low rrate
  assign(paste0("high_diff_regions_AUTO_", x), subset(get(paste0("FstDxy_lin_rrate_", x)), avg_wc_fst >= fst_threshold & avg_mean_rho > 25))
  assign(paste0("fake_high_diff_regions_AUTO_", x), subset(get(paste0("FstDxy_lin_rrate_", x)), avg_wc_fst >= fst_threshold & avg_mean_rho <= 25))
  
  # Create the scatter plot
  scatter_plot <- ggplot(get(paste0("FstDxy_lin_rrate_", x)), aes(x = avg_mean_rho, y = avg_wc_fst)) +
    geom_point(alpha = 0.5, size = 2, stroke = 0,  color = "gray60") +
    geom_point(data = get(paste0("fake_high_diff_regions_AUTO_", x)), aes(x = avg_mean_rho, y = avg_wc_fst), color = "red", alpha = 1, size = 3, stroke = 0) +
    geom_point(data = get(paste0("high_diff_regions_AUTO_", x)), aes(x = avg_mean_rho, y = avg_wc_fst), color = "green4", alpha = 1, size = 4, stroke = 0) +
    geom_smooth(method = "lm", color = "black", se = FALSE) +
    labs(title = paste0(gsub("_", " vs ", x), " Autosomes"), x = NULL, y = NULL) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),         # No major gridlines
      panel.grid.minor = element_blank(),         # No minor gridlines
      panel.background = element_rect(fill = "white", color = "white"),
      plot.background = element_rect(fill = "white", color = "white")) +
    stat_cor(method = "pearson", 
             aes(label = after_stat(r)), 
             label.x = max(get(paste0("FstDxy_lin_rrate_", x))$avg_mean_rho, na.rm = TRUE) * 0.6, 
             label.y = max(get(paste0("FstDxy_lin_rrate_", x))$avg_wc_fst, na.rm = TRUE) * 0.9, 
             size = 9, 
             label.sep = "")
  #ggsave(paste0("/mnt/DATA/B_JUN_reseq/B_PopIndexes/H_PLOTS_2.0/F_rrate_vs_fst/scatter_plot_",x, "_AUTO.tiff"),  plot = scatter_plot, width = 3, height = 4, dpi = 280, compression = "lzw", bg ="white")
  beepr::beep()
}

##### Z Chr 

# Loop over each ChrNum value from 1 to 31
for (j in lineages) {
  assign(paste0("rmapDF_fst_", j), data.frame())
  LIN <- subset(FstDxy_lin,combined_pop == j )
  for (i in 31:31) {
    # Subset data for the current ChrNum value
    A <- subset(LIN, ChrNum == i)
    B <- subset(filt_rmapDF, ChrNum == i)
    
    # Apply rounding to the nearest multiple of 50000 for window_pos_1
    B$window_pos_1 <- round_down_to_nearest(B$POS, 50000)
    
    # Perform the left join
    merged_AB <- dplyr::left_join(B, A[, c("window_pos_1", "avg_wc_fst")], by = "window_pos_1")
    
    # Concatenate the result into the final data frame
    assign(paste0("rmapDF_fst_", j), rbind(get(paste0("rmapDF_fst_", j)), merged_AB))
  }
}
beepr::beep()

for (x in lineages) {
  # Calculate avg_mean_rho_df
  avg_mean_rho_df <- get(paste0("rmapDF_fst_", x)) %>%
    group_by(ChrNum, ChrName, window_pos_1) %>%
    summarise(avg_mean_rho = mean(Mean_rho, na.rm = TRUE), .groups = "drop")  # Ensure to drop grouping
  
  # Perform left join
  FstDxy_lin_rrate <- dplyr::left_join(FstDxy_lin, avg_mean_rho_df, by = c("ChrNum", "ChrName", "window_pos_1"))
  
  # Filter FstDxy_lin_rrate by combined_pop == x and remove NA values in avg_mean_rho
  assign(paste0("FstDxy_lin_rrate_", x),
         FstDxy_lin_rrate %>%
           filter(!is.na(avg_mean_rho) & combined_pop == x))
  
  # Calculate the threshold for the top 0.5% FST
  fst_threshold <- quantile(get(paste0("FstDxy_lin_rrate_", x))$avg_wc_fst, 0.995, na.rm = TRUE)
  
  # Subset the data for regions with elevated differentiation and moderate to high recombination rate & elevated differentiation and low rrate
  assign(paste0("high_diff_regions_Zchr_", x), subset(get(paste0("FstDxy_lin_rrate_", x)), avg_wc_fst >= fst_threshold & avg_mean_rho > 25))
  assign(paste0("fake_high_diff_regions_Zchr_", x), subset(get(paste0("FstDxy_lin_rrate_", x)), avg_wc_fst >= fst_threshold & avg_mean_rho <= 25))
  
  # Create the scatter plot
  scatter_plot <- ggplot(get(paste0("FstDxy_lin_rrate_", x)), aes(x = avg_mean_rho, y = avg_wc_fst)) +
    geom_point(alpha = 0.5, size = 2, stroke = 0,  color = "gray60") +
    geom_point(data = get(paste0("fake_high_diff_regions_Zchr_", x)), aes(x = avg_mean_rho, y = avg_wc_fst), color = "red", alpha = 1, size = 3, stroke = 0) +
    geom_point(data = get(paste0("high_diff_regions_Zchr_", x)), aes(x = avg_mean_rho, y = avg_wc_fst), color = "green4", alpha = 1, size = 4, stroke = 0) +
    geom_smooth(method = "lm", color = "black", se = FALSE) +
    labs(title = paste0(gsub("_", " vs ", x), " Z Chr"), x = NULL, y = NULL) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),         # No major gridlines
      panel.grid.minor = element_blank(),         # No minor gridlines
      panel.background = element_rect(fill = "white", color = "white"),
      plot.background = element_rect(fill = "white", color = "white")) +
    stat_cor(method = "pearson", 
             aes(label = after_stat(r)), 
             label.x = max(get(paste0("FstDxy_lin_rrate_", x))$avg_mean_rho, na.rm = TRUE) * 0.6, 
             label.y = max(get(paste0("FstDxy_lin_rrate_", x))$avg_wc_fst, na.rm = TRUE) * 0.9, 
             size = 9, 
             label.sep = "")
  #ggsave(paste0("/mnt/DATA/B_JUN_reseq/B_PopIndexes/H_PLOTS_2.0/F_rrate_vs_fst/scatter_plot_",x, "_Zchr.tiff"),  plot = scatter_plot, width = 3, height = 4, dpi = 280, compression = "lzw", bg ="white")
  beepr::beep()
}

for (x in lineages) {
  assign(paste0("high_diff_regions_", x), dplyr::bind_rows(get(paste0("high_diff_regions_Zchr_", x)),  get(paste0("high_diff_regions_AUTO_", x))))
  assign(paste0("fake_high_diff_regions_", x), dplyr::bind_rows(get(paste0("fake_high_diff_regions_Zchr_", x)),  get(paste0("fake_high_diff_regions_AUTO_", x))))
}

##### PLOT recomb rate vs Fst population ----

# Initialize an empty data frame to store the results
FstDxy_pop <- FstDxy_pop %>%
  unite(combined_pop, pop1, pop2, sep = "_", remove = FALSE)

population <- c("DOE_DOW", "THU_TOW", "CAR_HYE","DOW_PAL" ,"CAN_DOW", "AIK_HYE", "MON_PIN", "MON_ORE", "CAN_MEA", "AIK_MEA", "MEA_MON", "MON_THU", "PAL_PHA" )

MIX1_pop <- paste(POP1_pop,POP2_pop,sep="") # DOE DOW
MIX2_pop <- paste(POP3_pop,POP4_pop,sep="") # THU TOW
MIX3_pop <- paste(POP5_pop,POP6_pop,sep="") # CAR HYE
MIX4_pop <- paste(POP2_pop,POP7_pop,sep="") # DOW PAL
MIX5_pop <- paste(POP12_pop,POP2_pop,sep="") # CAN DOW
MIX6_pop <- paste(POP8_pop,POP6_pop,sep="") # AIK HYE
MIX7_pop <- paste(POP9_pop,POP11_pop,sep="") # MON PIN
MIX8_pop <- paste(POP9_pop,POP10_pop,sep="") # MON ORE
MIX9_pop <- paste(POP12_pop,POP13_pop,sep="" ) # CAN MEA
MIX10_pop <- paste(POP8_pop,POP13_pop,sep="" ) # AIK MEA
MIX11_pop <- paste(POP13_pop,POP9_pop,sep="" ) # MEA MON
MIX12_pop <- paste(POP9_pop,POP3_pop, sep="") #MON THU
MIX13_pop <- paste(POP7_pop,POP14_pop, sep="") # PAL PHA

### AUTOSOMES
# Loop over each ChrNum value from 1 to 30
for (j in population) {
  assign(paste0("rmapDF_fst_", j), data.frame())
  POPUL <- subset(FstDxy_pop,combined_pop == j )
  for (i in 1:30) {
    # Subset data for the current ChrNum value
    A <- subset(POPUL, ChrNum == i)
    B <- subset(filt_rmapDF, ChrNum == i)
    
    # Apply rounding to the nearest multiple of 50000 for window_pos_1
    B$window_pos_1 <- round_down_to_nearest(B$POS, 50000)
    
    # Perform the left join
    merged_AB <- dplyr::left_join(B, A[, c("window_pos_1", "avg_wc_fst")], by = "window_pos_1")
    
    # Concatenate the result into the final data frame
    assign(paste0("rmapDF_fst_", j), rbind(get(paste0("rmapDF_fst_", j)), merged_AB))
  }
  beepr::beep()
}

for (x in population) {
  # Calculate avg_mean_rho_df
  avg_mean_rho_df <- get(paste0("rmapDF_fst_", x)) %>%
    group_by(ChrNum, ChrName, window_pos_1) %>%
    summarise(avg_mean_rho = mean(Mean_rho, na.rm = TRUE), .groups = "drop")  # Ensure to drop grouping
  
  # Perform left join
  FstDxy_pop_rrate <- dplyr::left_join(FstDxy_pop, avg_mean_rho_df, by = c("ChrNum", "ChrName", "window_pos_1"))
  
  # Filter FstDxy_pop_rrate by combined_pop == x and remove NA values in avg_mean_rho
  assign(paste0("FstDxy_pop_rrate_", x),
         FstDxy_pop_rrate %>%
           filter(!is.na(avg_mean_rho) & combined_pop == x))
  
  # Calculate the threshold for the top 0.5% FST
  fst_threshold <- quantile(get(paste0("FstDxy_pop_rrate_", x))$avg_wc_fst, 0.995, na.rm = TRUE)
  
  # Subset the data for regions with elevated differentiation and moderate to high recombination rate & elevated differentiation and low rrate
  assign(paste0("high_diff_regions_AUTO_", x), subset(get(paste0("FstDxy_pop_rrate_", x)), avg_wc_fst >= fst_threshold & avg_mean_rho > 25))
  assign(paste0("fake_high_diff_regions_AUTO_", x), subset(get(paste0("FstDxy_pop_rrate_", x)), avg_wc_fst >= fst_threshold & avg_mean_rho <= 25))
  
  # Create the scatter plot
  scatter_plot <- ggplot(get(paste0("FstDxy_pop_rrate_", x)), aes(x = avg_mean_rho, y = avg_wc_fst)) +
    geom_point(alpha = 0.5, size = 2, stroke = 0,  color = "gray60") +
    geom_point(data = get(paste0("fake_high_diff_regions_AUTO_", x)), aes(x = avg_mean_rho, y = avg_wc_fst), color = "red", alpha = 1, size = 3, stroke = 0) +
    geom_point(data = get(paste0("high_diff_regions_AUTO_", x)), aes(x = avg_mean_rho, y = avg_wc_fst), color = "green4", alpha = 1, size = 4, stroke = 0) +
    geom_smooth(method = "lm", color = "black", se = FALSE) +
    labs(title = paste0(gsub("_", " vs ", x), " Autosomes"), x = NULL, y = NULL) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),         # No major gridlines
      panel.grid.minor = element_blank(),         # No minor gridlines
      panel.background = element_rect(fill = "white", color = "white"),
      plot.background = element_rect(fill = "white", color = "white")) +
    stat_cor(method = "pearson", 
             aes(label = after_stat(r)), 
             label.x = max(get(paste0("FstDxy_pop_rrate_", x))$avg_mean_rho, na.rm = TRUE) * 0.6, 
             label.y = max(get(paste0("FstDxy_pop_rrate_", x))$avg_wc_fst, na.rm = TRUE) * 0.9, 
             size = 9, 
             label.sep = "")
  #ggsave(paste0("/mnt/DATA/B_JUN_reseq/B_PopIndexes/H_PLOTS_2.0/F_rrate_vs_fst/scatter_plot_",x, "_AUTO.tiff"),  plot = scatter_plot, width = 3, height = 4, dpi = 280, compression = "lzw", bg ="white")
  beepr::beep()
}

### SEX Z chr
# Loop over each ChrNum value from 31 to 31
for (j in population) {
  assign(paste0("rmapDF_fst_", j), data.frame())
  POPUL <- subset(FstDxy_pop,combined_pop == j )
  for (i in 31:31) {
    # Subset data for the current ChrNum value
    A <- subset(POPUL, ChrNum == i)
    B <- subset(filt_rmapDF, ChrNum == i)
    
    # Apply rounding to the nearest multiple of 50000 for window_pos_1
    B$window_pos_1 <- round_down_to_nearest(B$POS, 50000)
    
    # Perform the left join
    merged_AB <- dplyr::left_join(B, A[, c("window_pos_1", "avg_wc_fst")], by = "window_pos_1")
    
    # Concatenate the result into the final data frame
    assign(paste0("rmapDF_fst_", j), rbind(get(paste0("rmapDF_fst_", j)), merged_AB))
  }
}
beepr::beep()

for (x in population) {
  # Calculate avg_mean_rho_df
  avg_mean_rho_df <- get(paste0("rmapDF_fst_", x)) %>%
    group_by(ChrNum, ChrName, window_pos_1) %>%
    summarise(avg_mean_rho = mean(Mean_rho, na.rm = TRUE), .groups = "drop")  # Ensure to drop grouping
  
  # Perform left join
  FstDxy_pop_rrate <- dplyr::left_join(FstDxy_pop, avg_mean_rho_df, by = c("ChrNum", "ChrName", "window_pos_1"))
  
  # Filter FstDxy_pop_rrate by combined_pop == x and remove NA values in avg_mean_rho
  assign(paste0("FstDxy_pop_rrate_", x),
         FstDxy_pop_rrate %>%
           filter(!is.na(avg_mean_rho) & combined_pop == x))
  
  # Calculate the threshold for the top 0.5% FST
  fst_threshold <- quantile(get(paste0("FstDxy_pop_rrate_", x))$avg_wc_fst, 0.995, na.rm = TRUE)
  
  # Subset the data for regions with elevated differentiation and moderate to high recombination rate & elevated differentiation and low rrate
  assign(paste0("high_diff_regions_Zchr_", x), subset(get(paste0("FstDxy_pop_rrate_", x)), avg_wc_fst >= fst_threshold & avg_mean_rho > 25))
  assign(paste0("fake_high_diff_regions_Zchr_", x), subset(get(paste0("FstDxy_pop_rrate_", x)), avg_wc_fst >= fst_threshold & avg_mean_rho <= 25))
  
  # Create the scatter plot
  scatter_plot <- ggplot(get(paste0("FstDxy_pop_rrate_", x)), aes(x = avg_mean_rho, y = avg_wc_fst)) +
    geom_point(alpha = 0.5, size = 2, stroke = 0,  color = "gray60") +
    geom_point(data = get(paste0("fake_high_diff_regions_Zchr_", x)), aes(x = avg_mean_rho, y = avg_wc_fst), color = "red", alpha = 1, size = 3, stroke = 0) +
    geom_point(data = get(paste0("high_diff_regions_Zchr_", x)), aes(x = avg_mean_rho, y = avg_wc_fst), color = "green4", alpha = 1, size = 4, stroke = 0) +
    geom_smooth(method = "lm", color = "black", se = FALSE) +
    labs(title = paste0(gsub("_", " vs ", x), " Z Chr"), x = NULL, y = NULL) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),         # No major gridlines
      panel.grid.minor = element_blank(),         # No minor gridlines
      panel.background = element_rect(fill = "white", color = "white"),
      plot.background = element_rect(fill = "white", color = "white")) +
    stat_cor(method = "pearson",
             aes(label = after_stat(r)), 
             label.x = max(get(paste0("FstDxy_pop_rrate_", x))$avg_mean_rho, na.rm = TRUE) * 0.6, 
             label.y = max(get(paste0("FstDxy_pop_rrate_", x))$avg_wc_fst, na.rm = TRUE) * 0.9, 
             size = 9, 
             label.sep = "")
  #ggsave(paste0("/mnt/DATA/B_JUN_reseq/B_PopIndexes/H_PLOTS_2.0/F_rrate_vs_fst/scatter_plot_",x, "_Zchr.tiff"),  plot = scatter_plot, width = 3, height = 4, dpi = 280, compression = "lzw", bg ="white")
  beepr::beep()
}

for (x in population) {
  assign(paste0("high_diff_regions_", x), dplyr::bind_rows(get(paste0("high_diff_regions_Zchr_", x)),  get(paste0("high_diff_regions_AUTO_", x))))
  assign(paste0("fake_high_diff_regions_", x), dplyr::bind_rows(get(paste0("fake_high_diff_regions_Zchr_", x)),  get(paste0("fake_high_diff_regions_AUTO_", x))))
}

##### Density plot Populations fst pi dxy ----

# Define the populations you want to filter (7 comparisons)
pop1_filter <- c("CAN", "AIK", "MON", "MON", "DOE", "THU", "CAR", "INS") 
pop2_filter <- c("DOW", "HYE", "PIN", "ORE", "DOW", "TOW", "HYE", "THU") 

# Create a data frame of the population pairs
population_pairs <- data.frame(pop1 = pop1_filter, pop2 = pop2_filter)

# Filter the data for specific population pairs
filtered_data_auto <- FstDxy_pop %>%
  filter(ChrNum != 31) %>%
  inner_join(population_pairs, by = c("pop1", "pop2")) %>%
  mutate(population_pair = paste(pop1, pop2, sep = "_")) %>%
  mutate(population_pair = factor(population_pair, levels = paste(pop1_filter, pop2_filter, sep = "_")))

filtered_data_Z <- FstDxy_pop %>%
  filter(ChrNum == 31) %>%
  inner_join(population_pairs, by = c("pop1", "pop2")) %>%
  mutate(population_pair = paste(pop1, pop2, sep = "_")) %>%
  mutate(population_pair = factor(population_pair, levels = paste(pop1_filter, pop2_filter, sep = "_")))

# Define a manual color palette
color_palette <- c("CAN_DOW" = "purple", "AIK_HYE" = "cyan", "MON_PIN" = "forestgreen", 
                   "MON_ORE" = "palegreen1", "DOE_DOW" = "purple4", "THU_TOW" = "olivedrab3", 
                   "CAR_HYE" = "skyblue2", "INS_THU" = "orange")

# Define a manual linetype palette
linetype_palette <- c("CAN_DOW" = "solid", "AIK_HYE" = "solid", "MON_PIN" = "solid", 
                      "MON_ORE" = "solid", "DOE_DOW" = "dashed", "THU_TOW" = "dashed", 
                      "CAR_HYE" = "solid", "INS_THU" = "dashed")

# Define the x-axis limits for Fst plot
x_axis_limits_fst <- c(-0.05, 0.3) 

# Plot the density plot for avg_wc_fst values per filtered population pair
fst_plot_auto <- ggplot(filtered_data_auto, aes(x = avg_wc_fst, fill = population_pair, color = population_pair, linetype = population_pair)) +
  geom_density(alpha = 0.05, linewidth = 1) +
  labs(x = "Fst Values", y = "Density") +
  theme_minimal() +
  xlim(x_axis_limits_fst) +
  scale_fill_manual(values = color_palette) +
  scale_color_manual(values = color_palette) +
  scale_linetype_manual(values = linetype_palette)+
  theme(legend.position = "none")  # Remove legend for now

fst_plot_Z <- ggplot(filtered_data_Z, aes(x = avg_wc_fst, fill = population_pair, color = population_pair, linetype = population_pair)) +
  geom_density(alpha = 0.05, linewidth = 1) +
  labs(x = "Fst Values", y = "Density") +
  theme_minimal() +
  xlim(x_axis_limits_fst) +
  scale_fill_manual(values = color_palette) +
  scale_color_manual(values = color_palette) +
  scale_linetype_manual(values = linetype_palette)  +
  theme(legend.position = "none")  # Remove legend for now

fst_plot_auto
fst_plot_Z

# Define the x-axis limits for Dxy plot
x_axis_limits_dxy <- c(0, 0.0105)

# Plot the density plot for avg_dxy values per filtered population pair
dxy_plot_auto <- ggplot(filtered_data_auto, aes(x = avg_dxy, fill = population_pair, color = population_pair, linetype = population_pair)) +
  geom_density(alpha = 0.05, size = 1) +
  labs(x = "Dxy Values", y = "Density") +
  theme_minimal() +
  xlim(x_axis_limits_dxy) +
  scale_fill_manual(values = color_palette) +
  scale_color_manual(values = color_palette) +
  scale_linetype_manual(values = linetype_palette)+
  theme(legend.position = "none")  # Remove legend for now

dxy_plot_Z <- ggplot(filtered_data_Z, aes(x = avg_dxy, fill = population_pair, color = population_pair, linetype = population_pair)) +
  geom_density(alpha = 0.05, size = 1) +
  labs(x = "Dxy Values", y = "Density") +
  theme_minimal() +
  xlim(x_axis_limits_dxy) +
  scale_fill_manual(values = color_palette) +
  scale_color_manual(values = color_palette) +
  scale_linetype_manual(values = linetype_palette)+
  theme(legend.position = "none")  # Remove legend for now

dxy_plot_auto
dxy_plot_Z

# Define the populations you want to filter
pop_filter <- c("AIK", "CAN", "CAR", "DOE", "DOW", "MON", "ORE", "HYE", "PIN", "THU", "TOW", "INS")

# Filter the data for specific populations and create the version without ChrNum == 31
filtered_pi_auto <- pi_pop %>%
  filter(pop %in% pop_filter & chromosome != "ScU1KuS_2__HRSCAF___21") %>%
  mutate(pop = factor(pop, levels = pop_filter))

# Filter the data for specific populations and create the version with ChrNum == 31
filtered_pi_Z <- pi_pop %>%
  filter(pop %in% pop_filter & chromosome == "ScU1KuS_2__HRSCAF___21") %>%
  mutate(pop = factor(pop, levels = pop_filter))

# Define the x-axis limits
x_axis_limits <- c(-0.0001, 0.01) # Replace with your desired x-axis limits

# Define a manual color palette
color_palette_pi <- c("AIK" = "cyan", "CAN" = "orchid1", "CAR" = "skyblue2", "DOE" = "purple4", 
                      "DOW" = "purple", "MON" = "forestgreen", "ORE" = "palegreen1", "HYE" = "#F781BF", 
                      "PIN" = "springgreen2", "THU" = "darkseagreen", "TOW" = "olivedrab3", "INS" = "orange")

# Define a manual linetype palette with denser dashing for DOE and TOW
linetype_palette_pi <- c("AIK" = "solid", "CAN" = "solid", "CAR" = "solid", "DOE" = "dashed", 
                         "DOW" = "solid", "MON" = "solid", "ORE" = "solid", "HYE" = "solid", 
                         "PIN" = "solid", "THU" = "solid", "TOW" = "dashed", "INS" = "dashed")

# Plot the density plot for avg_pi values per filtered population
pi_plot_auto <- ggplot(filtered_pi_auto, aes(x = avg_pi, fill = pop, color = pop, linetype = pop)) +
  geom_density(alpha = 0.05, size = 1) +
  labs(x = "Pi Values", y = "Density") +
  theme_minimal() +
  xlim(x_axis_limits) +
  scale_fill_manual(values = color_palette_pi) + # Use manual color palette for fill
  scale_color_manual(values = color_palette_pi) + # Use manual color palette for lines
  scale_linetype_manual(values = linetype_palette_pi)+
  theme(legend.position = "none")  # Remove legend for now

# Plot the density plot for avg_pi values per filtered population
pi_plot_Z <- ggplot(filtered_pi_Z, aes(x = avg_pi, fill = pop, color = pop, linetype = pop)) +
  geom_density(alpha = 0.05, size = 1) +
  labs(x = "Pi Values", y = "Density") +
  theme_minimal() +
  xlim(x_axis_limits) +
  scale_fill_manual(values = color_palette_pi) + # Use manual color palette for fill
  scale_color_manual(values = color_palette_pi) + # Use manual color palette for lines
  scale_linetype_manual(values = linetype_palette_pi)+
  theme(legend.position = "none")  # Remove legend for now

# Arrange the three plots into a grid
grid.arrange(fst_plot_auto,fst_plot_Z, dxy_plot_auto, dxy_plot_Z, pi_plot_auto, pi_plot_Z, ncol = 2)

##### Incorporate xpehh outliers into fst plots SP ----

# Define the functions for rounding down and up to the nearest multiple
round_down_to_nearest <- function(x, base) {
  return(base * floor(x / base) + 1)
}

###  SP

xpehh_SP<-read_tsv("xpehh_DEJUYEJU.txt") # Already read

### XPEHH

xpehh_sp_DF <- xpehh_SP %>%
  rename(chromosome = CHR.name) %>%  # Rename CHR.name to chromosome
  merge(ChrTable, by.x = "chromosome", by.y = "ChrOrder", all = TRUE) %>%
  arrange(factor(chromosome, levels = ChrOrder), POSITION, POP) %>%
  mutate(SNPNum = seq(1, n())) %>%  # Using n() instead of length(CHR.name)
  mutate(ChrPos = POSITION) %>%
  group_by(chromosome)

# Combine all the steps into one pipeline
outliersXpehh_sp_MIX1 <- xpehh_sp_DF %>%
  filter(ChrNum >= CHRS & ChrNum <= CHRE & LOGPVALUE > 8) %>% # Apply the filter condition
  mutate(window_pos_1 = round_down_to_nearest(POSITION, 50000)) %>% # Create new column with rounded values
  distinct(ChrNum, window_pos_1, .keep_all = TRUE) # Remove duplicates based on ChrNum and window_pos_1

rm(xpehh_SP, xpehh_sp_DF)

FstDxy_XPEHH_MIX1_sp <- combined_data %>%
  inner_join(outliersXpehh_sp_MIX1 %>% select(chromosome,window_pos_1, ChrNum, ChrName),
             by = c("chromosome","window_pos_1", "ChrNum", "ChrName" )) %>%
  filter(pop1 == POP1_sp & pop2 == POP2_sp) %>%
  distinct() %>%
  subset(FstPcorrected<0.05)

# Only outliers high fst & xpehh
unique_outliersFstDxy_XPEHH_MIX1_sp <- FstDxy_XPEHH_MIX1_sp %>%
  anti_join(high_diff_regions_SP, by = c("window_pos_1", "ChrName"))

# Only outliers high rrate+fst & xpehh
shared_outliersFstDxy_XPEHH_MIX1_sp <- FstDxy_XPEHH_MIX1_sp %>%
  anti_join(unique_outliersFstDxy_XPEHH_MIX1_sp, by = c("window_pos_1", "ChrName"))

# Create new DF with outliers Fst + rrate that erases coincidenced with fst + rrate + xpehh for plotting purposes
real_high_diff_regions_SP <- high_diff_regions_SP %>%
  anti_join(shared_outliersFstDxy_XPEHH_MIX1_sp, by = c("window_pos_1", "ChrName"))

#### LINEAGES

xpehh_LIN <- read_tsv("xpehh_concat_LIN.txt")

### XPEHH

xpehh_lin_DF <- xpehh_LIN %>%
  rename(chromosome = CHR.name) %>%  # Rename CHR.name to chromosome
  merge(ChrTable, by.x = "chromosome", by.y = "ChrOrder", all = TRUE) %>%
  arrange(factor(chromosome, levels = ChrOrder), POSITION, POP) %>%
  mutate(SNPNum = seq(1, n())) %>%  # Using n() instead of length(CHR.name)
  mutate(ChrPos = POSITION) %>%
  group_by(chromosome)

lineages2 <- c("CANspHYEsp", "CANspOREsp", "HYEspOREsp" )

# Combine all the steps into one pipeline for each lineage
for (x in lineages2) {
  outliersXpehh <- xpehh_lin_DF %>%
    filter(ChrNum >= CHRS & ChrNum <= CHRE & LOGPVALUE > 8 & POP == x) %>% # Apply the filter condition
    mutate(window_pos_1 = round_down_to_nearest(POSITION, 50000)) %>% # Create new column with rounded values
    distinct(ChrNum, window_pos_1, .keep_all = TRUE) # Remove duplicates based on ChrNum and window_pos_1
  assign(paste0("outliersXpehh_", x), outliersXpehh)
}

rm(xpehh_lin_DF, xpehh_LIN)

for (x in lineages2) {
  # Split x into POP1 and POP2 based on the pattern
  lineage_parts <- unlist(strsplit(x, "(?<=[a-z])(?=[A-Z])", perl=TRUE))
  POP1 <- lineage_parts[1]
  POP2 <- lineage_parts[2]
  assign(paste0("FstDxy_XPEHH_", x), 
         combined_data %>%
           inner_join(get(paste0("outliersXpehh_", x)), 
                      by = c("chromosome", "window_pos_1", "ChrNum", "ChrName")) %>%
           filter(pop1 == POP1 & pop2 == POP2) %>%
           distinct() %>%
           subset(FstPcorrected < 0.05))
  # Only outliers high fst & xpehh
  outliers <- anti_join(get(paste0("FstDxy_XPEHH_", x)), get(paste0("high_diff_regions_", POP1,"_" ,POP2)), by = c("window_pos_1", "ChrNum"))
  assign(paste0("unique_outliersFstDxy_XPEHH_", x), outliers)
  # Only outliers high rrate+fst & xpehh
  assign(paste0("shared_outliersFstDxy_XPEHH_", x),
         get(paste0("FstDxy_XPEHH_", x)) %>%
           anti_join(get(paste0("unique_outliersFstDxy_XPEHH_", POP1, POP2)), by = c("window_pos_1", "ChrNum")))
  # Create new DF with outliers Fst + rrate that erases coincidenced with fst + rrate + xpehh for plotting purposes
  assign(paste0("real_high_diff_regions_", x),
         get(paste0("high_diff_regions_", POP1, "_", POP2)) %>%
           anti_join(get(paste0("shared_outliersFstDxy_XPEHH_", POP1, POP2)), by = c("window_pos_1", "ChrNum")))
}

#### POP
xpehh_POP <- read_tsv("xpehh_concat_POP.txt")

### XPEHH
xpehh_pop_DF <- xpehh_POP %>%
  rename(chromosome = CHR.name) %>%  # Rename CHR.name to chromosome
  merge(ChrTable, by.x = "chromosome", by.y = "ChrOrder", all = TRUE) %>%
  arrange(factor(chromosome, levels = ChrOrder), POSITION, POP) %>%
  mutate(SNPNum = seq(1, n())) %>%  # Using n() instead of length(CHR.name)
  mutate(ChrPos = POSITION) %>%
  group_by(chromosome)

population2 <- c("DOEDOW", "THUTOW", "CARHYE","DOWPAL" ,"CANDOW", "AIKHYE", "MONPIN", "MONORE", "CANMEA", "AIKMEA", "MEAMON", "MONTHU", "PALPHA" )
# Combine all the steps into one pipeline for each lineage
for (x in population2) {
  outliersXpehh <- xpehh_pop_DF %>%
    filter(ChrNum >= CHRS & ChrNum <= CHRE & LOGPVALUE > 8 & POP == x) %>% # Apply the filter condition
    mutate(window_pos_1 = round_down_to_nearest(POSITION, 50000)) %>% # Create new column with rounded values
    distinct(ChrNum, window_pos_1, .keep_all = TRUE) # Remove duplicates based on ChrNum and window_pos_1
  assign(paste0("outliersXpehh_", x), outliersXpehh)
}


rm(xpehh_pop_DF, xpehh_POP)

for (x in population2) {
  # Split the element into two parts at the specific position (after the third character)
  POP1 <- substr(x, 1, 3)
  POP2 <- substr(x, 4, nchar(x))
  assign(paste0("FstDxy_XPEHH_", x), 
         combined_data %>%
           inner_join(get(paste0("outliersXpehh_", x)), 
                      by = c("chromosome", "window_pos_1", "ChrNum", "ChrName")) %>%
           filter(pop1 == POP1 & pop2 == POP2) %>%
           distinct() %>%
           subset(FstPcorrected < 0.05))
  # Only outliers high fst & xpehh
  outliers <- anti_join(get(paste0("FstDxy_XPEHH_", x)), get(paste0("high_diff_regions_", POP1,"_" ,POP2)), by = c("window_pos_1", "ChrNum"))
  assign(paste0("unique_outliersFstDxy_XPEHH_", x), outliers)
  # Only outliers high rrate+fst & xpehh
  assign(paste0("shared_outliersFstDxy_XPEHH_", x),
         get(paste0("FstDxy_XPEHH_", x)) %>%
           anti_join(get(paste0("unique_outliersFstDxy_XPEHH_", POP1, POP2)), by = c("window_pos_1", "ChrNum")))
  # Create new DF with outliers Fst + rrate that erases coincidenced with fst + rrate + xpehh for plotting purposes
  assign(paste0("real_high_diff_regions_", x),
         get(paste0("high_diff_regions_", POP1, "_", POP2)) %>%
           anti_join(get(paste0("shared_outliersFstDxy_XPEHH_", POP1, POP2)), by = c("window_pos_1", "ChrNum")))
}


##### Plotting fst continuum with rrate and xpehh consideration ----

manhattanJAVI_2 <- function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", 
                             col = c("gray10", "gray60"), chrlabs = NULL, 
                             suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08), 
                             highlight1 = NULL, highlight2 = NULL, highlight3 = NULL, 
                             highlight4 = NULL, highlight5 = NULL, highlight6 = NULL, 
                             highlight7 = NULL, logp = TRUE, annotatePval = NULL, 
                             annotateTop = TRUE, cex.highlight = 1, ...) 
{
  CHR = BP = P = index = NULL
  if (!(chr %in% names(x))) 
    stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x))) 
    stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x))) 
    stop(paste("Column", p, "not found!"))
  if (!(snp %in% names(x))) 
    warning(paste("No SNP column found. OK unless you're trying to highlight."))
  if (!is.numeric(x[[chr]])) 
    stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]])) 
    stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]])) 
    stop(paste(p, "column should be numeric."))
  if (!is.null(x[[snp]])) 
    d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]], 
                   pos = NA, index = NA, SNP = x[[snp]], stringsAsFactors = FALSE)
  else d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]], 
                      pos = NA, index = NA)
  d <- d[order(d$CHR, d$BP), ]
  if (logp) {
    d$logp <- -log10(d$P)
  }
  else {
    d$logp <- d$P
  }
  d$index = rep.int(seq_along(unique(d$CHR)), times = tapply(d$SNP, 
                                                             d$CHR, length))
  nchr = length(unique(d$CHR))
  if (nchr == 1) {
    d$pos = d$BP
    xlabel = paste("Chromosome", unique(d$CHR), "position")
  }
  else {
    lastbase = 0
    ticks = NULL
    for (i in unique(d$index)) {
      if (i == 1) {
        d[d$index == i, ]$pos = d[d$index == i, ]$BP
      }
      else {
        lastbase = lastbase + max(d[d$index == (i - 
                                                  1), "BP"])
        d[d$index == i, "BP"] = d[d$index == i, "BP"] - 
          min(d[d$index == i, "BP"]) + 1
        d[d$index == i, "pos"] = d[d$index == i, "BP"] + 
          lastbase
      }
    }
    ticks <- tapply(d$pos, d$index, quantile, probs = 0.5)
    xlabel = "Chromosome"
    labs <- unique(d$CHR)
  }
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
                   las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0, 
                                                                     ceiling(max(d$logp))), xlab = xlabel, ylab = expression(-log[10](italic(p))))
  dotargs <- list(...)
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% 
                                            names(dotargs)]))
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs) == length(labs)) {
        labs <- chrlabs
      }
      else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    }
    else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }
  if (nchr == 1) {
    axis(1, ...)
  }
  else {
    axis(1, at = ticks, labels = labs, ...)
  }
  col = rep_len(col, max(d$index))
  if (nchr == 1) {
    with(d, points(pos, logp, pch = 20, col = col[1], ...))
  }
  else {
    icol = 1
    for (i in unique(d$index)) {
      points(d[d$index == i, "pos"], d[d$index == i, "logp"], 
             col = col[icol], pch = 20, ...)
      icol = icol + 1
    }
  }
  if (suggestiveline) 
    abline(h = suggestiveline, col = "blue")
  if (genomewideline) 
    abline(h = genomewideline, col = "red")
  if (!is.null(highlight1)) {
    if (any(!(highlight1 %in% d$SNP))) 
      warning("You're trying to highlight1 SNPs that don't exist in your results.")
    d.highlight1 = d[which(d$SNP %in% highlight1), ]
    with(d.highlight1, points(pos, logp, col = "blue", pch = 20, cex = cex.highlight, ...))
  }
  if (!is.null(highlight2)) {
    if (any(!(highlight2 %in% d$SNP))) 
      warning("You're trying to highlight2 SNPs that don't exist in your results.")
    d.highlight2 = d[which(d$SNP %in% highlight2), ]
    with(d.highlight2, points(pos, logp, col = "purple", pch = 20, ...))
  }
  if (!is.null(highlight3)) {
    if (any(!(highlight3 %in% d$SNP))) 
      warning("You're trying to highlight3 SNPs that don't exist in your results.")
    d.highlight3 = d[which(d$SNP %in% highlight3), ]
    with(d.highlight3, points(pos, logp, col = "red", pch = 20, cex = cex.highlight, ...))
  }
  if (!is.null(highlight4)) {
    if (any(!(highlight4 %in% d$SNP))) 
      warning("You're trying to highlight4 SNPs that don't exist in your results.")
    d.highlight4 = d[which(d$SNP %in% highlight4), ]
    with(d.highlight4, points(pos, logp, col = "gold", pch = 20, ...))
  }
  if (!is.null(highlight5)) {
    if (any(!(highlight5 %in% d$SNP))) 
      warning("You're trying to highlight5 SNPs that don't exist in your results.")
    d.highlight5 = d[which(d$SNP %in% highlight5), ]
    with(d.highlight5, points(pos, logp, col = "green", pch = 20, cex = cex.highlight, ...))
  }
  if (!is.null(annotatePval)) {
    if (logp) {
      topHits = subset(d, P <= annotatePval)
    }
    else topHits = subset(d, P >= annotatePval)
    par(xpd = TRUE)
    if (annotateTop == FALSE) {
      if (logp) {
        with(subset(d, P <= annotatePval), textxy(pos, 
                                                  -log10(P), offset = 0.625, labs = topHits$SNP, 
                                                  cex = 0.45), ...)
      }
      else with(subset(d, P >= annotatePval), textxy(pos, 
                                                     P, offset = 0.625, labs = topHits$SNP, cex = 0.45), 
                ...)
    }
    else {
      topHits <- topHits[order(topHits$P), ]
      topSNPs <- NULL
      for (i in unique(topHits$CHR)) {
        chrSNPs <- topHits[topHits$CHR == i, ]
        topSNPs <- rbind(topSNPs, chrSNPs[1, ])
      }
      if (logp) {
        textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.625, 
               labs = topSNPs$SNP, cex = 0.5, ...)
      }
      else textxy(topSNPs$pos, topSNPs$P, offset = 0.625, 
                  labs = topSNPs$SNP, cex = 0.5, ...)
    }
  }
  par(xpd = FALSE)
}

n=13

par(mfrow=c(n,1),oma=c(3,0,0,0),mar=c(1,5,1,1))

manhattanJAVI_2(subset(FstDxy_sp,pop1==POP1_sp & pop2==POP2_sp & ChrNum>=CHRS & ChrNum<=CHRE),#Fst SP
              chr="ChrNum", bp="ChrPos", p="avg_wc_fst",logp="FALSE", snp="WinNum",
              suggestiveline=FALSE,genomewideline=FALSE,
              highlight1 = FstDxy_XPEHH_MIX1_sp$WinNum,   # yellow shared outliers xpehh and fst (no rrate)
              highlight3 = shared_outliersFstDxy_XPEHH_MIX1_sp$WinNum,   # 
              highlight5 = real_high_diff_regions_SP$WinNum,
              #highlight4 = fake_high_diff_regions_SP$WinNum,
              ylim = c(0,1),col = c("grey40", "grey"), cex=1,
              ylab ="Fst \nDEJU YEJU", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI_2(subset(FstDxy_lin,pop1==POP1_lin & pop2==POP2_lin & ChrNum>=CHRS & ChrNum<=CHRE),#Fst LIN1-2
              chr="ChrNum", bp="ChrPos", p="avg_wc_fst",logp="FALSE", snp="WinNum",
              suggestiveline=FALSE,genomewideline=FALSE,
              highlight1 = FstDxy_XPEHH_CANspHYEsp$WinNum,   # yellow shared outliers xpehh and fst (no rrate)
              highlight5 = shared_outliersFstDxy_XPEHH_CANspHYEsp$WinNum,   # 
              highlight3 = real_high_diff_regions_CANspHYEsp$WinNum,
              ylim = c(0,0.6),col = c("grey40", "grey"), cex=1,
              ylab ="Fst \nCANsp HYEsp", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI_2(subset(FstDxy_lin,pop1==POP1_lin & pop2==POP3_lin & ChrNum>=CHRS & ChrNum<=CHRE),#Fst LIN 1-3
              chr="ChrNum", bp="ChrPos", p="avg_wc_fst",logp="FALSE", snp="WinNum",
              suggestiveline=FALSE,genomewideline=FALSE,
              highlight1 = FstDxy_XPEHH_CANspOREsp$WinNum,   # yellow shared outliers xpehh and fst (no rrate)
              highlight5 = shared_outliersFstDxy_XPEHH_CANspOREsp$WinNum,   # 
              highlight3 = real_high_diff_regions_CANspOREsp$WinNum,
              ylim = c(0,0.6),col = c("grey40", "grey"), cex=1,
              ylab ="Fst \nCANsp OREsp", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI_2(subset(FstDxy_lin,pop1==POP2_lin & pop2==POP3_lin & ChrNum>=CHRS & ChrNum<=CHRE),#Fst LIN 2-3
              chr="ChrNum", bp="ChrPos", p="avg_wc_fst",logp="FALSE", snp="WinNum",
              suggestiveline=FALSE,genomewideline=FALSE,
              highlight1 = FstDxy_XPEHH_HYEspOREsp$WinNum,   # yellow shared outliers xpehh and fst (no rrate)
              highlight5 = shared_outliersFstDxy_XPEHH_HYEspOREsp$WinNum,   # 
              highlight3 = real_high_diff_regions_HYEspOREsp$WinNum,
              ylim = c(0,0.6),col = c("grey40", "grey"), cex=1,
              ylab ="Fst Z-score\nHYEsp OREsp", chrlabs = ChrName[CHRS:CHRE])




par(mfrow=c(n,1),oma=c(3,0,0,0),mar=c(1,5,1,1))

manhattanJAVI_2(subset(FstDxy_pop,pop1==POP7_pop & pop2==POP14_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst PAL PHA
              chr="ChrNum", bp="ChrPos", p="avg_wc_fst",logp="FALSE", snp="WinNum",
              suggestiveline=FALSE,genomewideline=FALSE,
              highlight1 = FstDxy_XPEHH_PALPHA$WinNum,   # yellow shared outliers xpehh + fst (no rrate)
              highlight3 = shared_outliersFstDxy_XPEHH_PALPHA$WinNum,   # Red high rrate+fst & xpehh
              highlight4 = real_high_diff_regions_PALPHA$WinNum,   # Green ONLY outliers fst + rrate
              ylim = c(0,0.65),col = c("grey40", "grey"), cex=1,
              ylab ="Fst \nPHA - PAL", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI_2(subset(FstDxy_pop,pop1==POP2_pop & pop2==POP7_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst  DOW PAL
              chr="ChrNum", bp="ChrPos", p="avg_wc_fst",logp="FALSE", snp="WinNum",
              suggestiveline=FALSE,genomewideline=FALSE,
              highlight1 = FstDxy_XPEHH_DOWPAL$WinNum,   # Yellow shared outliers xpehh + fst (no rrate)
              highlight3 = shared_outliersFstDxy_XPEHH_DOWPAL$WinNum,   # Red high rrate+fst & xpehh
              highlight4 = real_high_diff_regions_DOWPAL$WinNum,   # Green ONLY outliers fst + rrate (rm xpehh)
              ylim = c(0,0.65),col = c("grey40", "grey"), cex=1,
              ylab ="Fst \nPAL - DOW", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI_2(subset(FstDxy_pop,pop1==POP1_pop & pop2==POP2_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst  DOE DOW
              chr="ChrNum", bp="ChrPos", p="avg_wc_fst",logp="FALSE", snp="WinNum",
              suggestiveline=FALSE,genomewideline=FALSE,
              highlight1 = FstDxy_XPEHH_DOEDOW$WinNum,   # Yellow shared outliers xpehh + fst (no rrate)
              highlight3 = shared_outliersFstDxy_XPEHH_DOEDOW$WinNum,   # Red high rrate+fst & xpehh
              highlight4 = real_high_diff_regions_DOEDOW$WinNum,   # Green ONLY outliers fst + rrate (rm xpehh)
              ylim = c(0,0.3),col = c("grey40", "grey"), cex=1,
              ylab ="Fst \nDOW - DOE", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI_2(subset(FstDxy_pop,pop1==POP12_pop & pop2==POP2_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst  CAN DOW
              chr="ChrNum", bp="ChrPos", p="avg_wc_fst",logp="FALSE", snp="WinNum",
              suggestiveline=FALSE,genomewideline=FALSE,
              highlight1 = FstDxy_XPEHH_CANDOW$WinNum,   # Yellow shared outliers xpehh + fst (no rrate)
              highlight3 = shared_outliersFstDxy_XPEHH_CANDOW$WinNum,   # Red high rrate+fst & xpehh
              highlight5 = real_high_diff_regions_CANDOW$WinNum,   # Green ONLY outliers fst + rrate (rm xpehh)
              ylim = c(0,0.3),col = c("grey40", "grey"), cex=1,
              ylab ="Fst \nDOW - CAN", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI_2(subset(FstDxy_pop,pop1==POP12_pop & pop2==POP13_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst  CAN MEA
                chr="ChrNum", bp="ChrPos", p="avg_wc_fst",logp="FALSE", snp="WinNum",
                suggestiveline=FALSE,genomewideline=FALSE,
                highlight1 = FstDxy_XPEHH_CANMEA$WinNum,   # Yellow shared outliers xpehh + fst (no rrate)
                highlight3 = shared_outliersFstDxy_XPEHH_CANMEA$WinNum,   # Red high rrate+fst & xpehh
                highlight4 = real_high_diff_regions_CANMEA$WinNum,   # Green ONLY outliers fst + rrate (rm xpehh)
                ylim = c(0,0.3),col = c("grey40", "grey"), cex=1,
                ylab ="Fst \nCAN - MEA", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI_2(subset(FstDxy_pop,pop1==POP13_pop & pop2==POP9_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst  MEA MON
                chr="ChrNum", bp="ChrPos", p="avg_wc_fst",logp="FALSE", snp="WinNum",
                suggestiveline=FALSE,genomewideline=FALSE,
                highlight1 = FstDxy_XPEHH_MEAMON$WinNum,   # Yellow shared outliers xpehh + fst (no rrate)
                highlight3 = shared_outliersFstDxy_XPEHH_MEAMON$WinNum,   # Red high rrate+fst & xpehh
                highlight4 = real_high_diff_regions_MEAMON$WinNum,   # Green ONLY outliers fst + rrate (rm xpehh)
                ylim = c(0,0.3),col = c("grey40", "grey"), cex=1,
                ylab ="Fst \nMEA - MON", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI_2(subset(FstDxy_pop,pop1==POP9_pop & pop2==POP10_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst  MON ORE
              chr="ChrNum", bp="ChrPos", p="avg_wc_fst",logp="FALSE", snp="WinNum",
              suggestiveline=FALSE,genomewideline=FALSE,
              highlight1 = FstDxy_XPEHH_MONORE$WinNum,   # Yellow shared outliers xpehh + fst (no rrate)
              highlight3 = shared_outliersFstDxy_XPEHH_MONORE$WinNum,   # Red high rrate+fst & xpehh
              highlight4 = real_high_diff_regions_MONORE$WinNum,   # Green ONLY outliers fst + rrate (rm xpehh)
              ylim = c(0,0.3),col = c("grey40", "grey"), cex=1,
              ylab ="Fst \nMON - ORE", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI_2(subset(FstDxy_pop,pop1==POP9_pop & pop2==POP11_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst MON PIN
              chr="ChrNum", bp="ChrPos", p="avg_wc_fst",logp="FALSE", snp="WinNum",
              suggestiveline=FALSE,genomewideline=FALSE,
              highlight1 = FstDxy_XPEHH_MONPIN$WinNum,   # Yellow shared outliers xpehh + fst (no rrate)
              highlight3 = shared_outliersFstDxy_XPEHH_MONPIN$WinNum,   # Red high rrate+fst & xpehh
              highlight4 = real_high_diff_regions_MONPIN$WinNum,   # Green ONLY outliers fst + rrate (rm xpehh)
              ylim = c(0,0.3),col = c("grey40", "grey"), cex=1,
              ylab ="Fst \nMON - PIN", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI_2(subset(FstDxy_pop,pop1==POP9_pop & pop2==POP3_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst MON THU
                chr="ChrNum", bp="ChrPos", p="avg_wc_fst",logp="FALSE", snp="WinNum",
                suggestiveline=FALSE,genomewideline=FALSE,
                highlight1 = FstDxy_XPEHH_MONTHU$WinNum,   # Yellow shared outliers xpehh + fst (no rrate)
                highlight3 = shared_outliersFstDxy_XPEHH_MONTHU$WinNum,   # Red high rrate+fst & xpehh
                highlight4 = real_high_diff_regions_MONTHU$WinNum,   # Green ONLY outliers fst + rrate (rm xpehh)
                ylim = c(0,0.3),col = c("grey40", "grey"), cex=1,
                ylab ="Fst \nMON - THU", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI_2(subset(FstDxy_pop,pop1==POP3_pop & pop2==POP4_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst THU TOW
              chr="ChrNum", bp="ChrPos", p="avg_wc_fst",logp="FALSE", snp="WinNum",
              suggestiveline=FALSE,genomewideline=FALSE,
              highlight1 = FstDxy_XPEHH_THUTOW$WinNum,   # Yellow shared outliers xpehh + fst (no rrate)
              highlight3 = shared_outliersFstDxy_XPEHH_THUTOW$WinNum,   # Red high rrate+fst & xpehh
              highlight5 = real_high_diff_regions_THUTOW$WinNum,   # Green ONLY outliers fst + rrate (rm xpehh)
              ylim = c(0,0.6),col = c("grey40", "grey"), cex=1,
              ylab ="Fst \nTHU - TOW", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI_2(subset(FstDxy_pop,pop1==POP8_pop & pop2==POP13_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst AIK MEA
                chr="ChrNum", bp="ChrPos", p="avg_wc_fst",logp="FALSE", snp="WinNum",
                suggestiveline=FALSE,genomewideline=FALSE,
                highlight1 = FstDxy_XPEHH_AIKMEA$WinNum,   # Yellow shared outliers xpehh + fst (no rrate)
                highlight3 = shared_outliersFstDxy_XPEHH_AIKMEA$WinNum,   # Red high rrate+fst & xpehh
                highlight4 = real_high_diff_regions_AIKMEA$WinNum,   # Green ONLY outliers fst + rrate (rm xpehh)
                ylim = c(0,0.3),col = c("grey40", "grey"), cex=1,
                ylab ="Fst \nMEA - AIK", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI_2(subset(FstDxy_pop,pop1==POP8_pop & pop2==POP6_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst AIK HYE
              chr="ChrNum", bp="ChrPos", p="avg_wc_fst",logp="FALSE", snp="WinNum",
              suggestiveline=FALSE,genomewideline=FALSE,
              highlight1 = FstDxy_XPEHH_AIKHYE$WinNum,   # Yellow shared outliers xpehh + fst (no rrate)
              highlight3 = shared_outliersFstDxy_XPEHH_AIKHYE$WinNum,   # Red high rrate+fst & xpehh
              highlight4 = real_high_diff_regions_AIKHYE$WinNum,   # Green ONLY outliers fst + rrate (rm xpehh)
              ylim = c(0,0.3),col = c("grey40", "grey"), cex=1,
              ylab ="Fst \nAIK - HYE", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI_2(subset(FstDxy_pop,pop1==POP5_pop & pop2==POP6_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst CAR HYE
              chr="ChrNum", bp="ChrPos", p="avg_wc_fst",logp="FALSE", snp="WinNum",
              suggestiveline=FALSE,genomewideline=FALSE,
              highlight1 = FstDxy_XPEHH_CARHYE$WinNum,   # Yellow shared outliers xpehh + fst (no rrate)
              highlight3 = shared_outliersFstDxy_XPEHH_CARHYE$WinNum,   # Red high rrate+fst & xpehh
              highlight4 = real_high_diff_regions_CARHYE$WinNum,   # Green ONLY outliers fst + rrate (rm xpehh)
              ylim = c(0,0.3),col = c("grey40", "grey"), cex=1,
              ylab ="Fst \nHYE - CAR", chrlabs = ChrName[CHRS:CHRE])


overlap_rrate_mix1_2 <- inner_join(high_diff_regions_CANsp_HYEsp, high_diff_regions_CANsp_OREsp, by = c("window_pos_1", "ChrNum"))
overlap_rrate_mix1_3 <- inner_join(high_diff_regions_CANsp_HYEsp, high_diff_regions_HYEsp_OREsp, by = c("window_pos_1", "ChrNum"))

overlap_rrate_mix2_1 <- inner_join(high_diff_regions_CANsp_OREsp, high_diff_regions_CANsp_HYEsp, by = c("window_pos_1", "ChrNum"))
overlap_rrate_mix2_3 <- inner_join(high_diff_regions_CANsp_OREsp, high_diff_regions_HYEsp_OREsp, by = c("window_pos_1", "ChrNum"))

overlap_rrate_mix3_2 <- inner_join(high_diff_regions_HYEsp_OREsp, high_diff_regions_CANsp_OREsp, by = c("window_pos_1", "ChrNum"))
overlap_rrate_mix3_1 <- inner_join(high_diff_regions_HYEsp_OREsp, high_diff_regions_CANsp_HYEsp, by = c("window_pos_1", "ChrNum"))

# Plot with Zscore
n=4

par(mfrow=c(n,1),oma=c(3,0,0,0),mar=c(1,5,1,1))

manhattanJAVI(subset(combined_data,pop1==POP1_lin & pop2==POP2_lin & ChrNum>=CHRS & ChrNum<=CHRE),#Fst canSP_hyeSP
              chr="ChrNum", bp="ChrPos", p="FstZscore",logp="FALSE", snp="WinNum",
              suggestiveline=FALSE,genomewideline=FALSE, 
              #highlight1 = FstDxy_XPEHH_CANspHYEsp$WinNum,   # yellow shared outliers xpehh + fst (no rrate)
              #highlight3 = shared_outliersFstDxy_XPEHH_CANspHYEsp$WinNum,   # Red high rrate+fst & xpehh
              highlight5 = c(overlap_rrate_mix1_2$WinNum.x, overlap_rrate_mix1_3$WinNum.x),   # Green ONLY outliers fst + rrate
              highlight4 = high_diff_regions_CANsp_HYEsp$WinNum,   # Green ONLY outliers fst + rrate
              ylim = c(-2,26),col = c("grey40", "grey"), cex=1,
              ylab ="Fst \nCANsp - HYEsp", chrlabs = ChrName[CHRS:CHRE])
              
manhattanJAVI(subset(combined_data,pop1==POP1_lin & pop2==POP3_lin & ChrNum>=CHRS & ChrNum<=CHRE),#Fst CANspOREsp
              chr="ChrNum", bp="ChrPos", p="FstZscore",logp="FALSE", snp="WinNum",
              suggestiveline=FALSE,genomewideline=FALSE, 
              #highlight1 = FstDxy_XPEHH_CANspHYEsp$WinNum,   # yellow shared outliers xpehh + fst (no rrate)
              highlight5 = c(overlap_rrate_mix2_3$WinNum.x, overlap_rrate_mix2_1$WinNum.x),   # Red high rrate+fst & xpehh
              highlight4 = high_diff_regions_CANsp_OREsp$WinNum,   # Green ONLY outliers fst + rrate
              ylim = c(-2,26),col = c("grey40", "grey"), cex=1,
              ylab ="Fst \nCANsp - OREsp", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI(subset(combined_data,pop1==POP2_lin & pop2==POP3_lin & ChrNum>=CHRS & ChrNum<=CHRE),#Fst HYEspOREsp
              chr="ChrNum", bp="ChrPos", p="FstZscore",logp="FALSE", snp="WinNum",
              suggestiveline=FALSE,genomewideline=FALSE, 
              #highlight1 = FstDxy_XPEHH_CANspHYEsp$WinNum,   # yellow shared outliers xpehh + fst (no rrate)
              highlight5 = c(overlap_rrate_mix3_1$WinNum.x, overlap_rrate_mix3_2$WinNum.x),   # Red high rrate+fst & xpehh
              highlight4 = high_diff_regions_HYEsp_OREsp$WinNum,   # Green ONLY outliers fst + rrate
              ylim = c(-2,26),col = c("grey40", "grey"), cex=1,
              ylab ="Fst \nHYEsp - OREsp", chrlabs = ChrName[CHRS:CHRE])


manhattanJAVI_2(subset(FstDxy_pop,pop1==POP7_pop & pop2==POP14_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst PAL PHA
                chr="ChrNum", bp="ChrPos", p="FstZscore",logp="FALSE", snp="WinNum",
                suggestiveline=FALSE,genomewideline=FALSE,
                highlight1 = FstDxy_XPEHH_PALPHA$WinNum,   # yellow shared outliers xpehh + fst (no rrate)
                highlight3 = shared_outliersFstDxy_XPEHH_PALPHA$WinNum,   # Red high rrate+fst & xpehh
                highlight4 = real_high_diff_regions_PALPHA$WinNum,   # Green ONLY outliers fst + rrate
                ylim = c(-2,26),col = c("grey40", "grey"), cex=1,
                ylab ="Fst \nPHA - PAL", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI_2(subset(FstDxy_pop,pop1==POP2_pop & pop2==POP7_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst  DOW PAL
                chr="ChrNum", bp="ChrPos", p="FstZscore",logp="FALSE", snp="WinNum",
                suggestiveline=FALSE,genomewideline=FALSE,
                highlight1 = FstDxy_XPEHH_DOWPAL$WinNum,   # Yellow shared outliers xpehh + fst (no rrate)
                highlight3 = shared_outliersFstDxy_XPEHH_DOWPAL$WinNum,   # Red high rrate+fst & xpehh
                highlight4 = real_high_diff_regions_DOWPAL$WinNum,   # Green ONLY outliers fst + rrate (rm xpehh)
                ylim = c(-2,26),col = c("grey40", "grey"), cex=1,
                ylab ="Fst \nPAL - DOW", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI(subset(FstDxy_pop,pop1==POP1_pop & pop2==POP2_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst  DOE DOW
                chr="ChrNum", bp="ChrPos", p="FstZscore",logp="FALSE", snp="WinNum",
                suggestiveline=FALSE,genomewideline=FALSE,
                highlight1 = FstDxy_XPEHH_DOEDOW$WinNum,   # Yellow shared outliers xpehh + fst (no rrate)
                highlight3 = shared_outliersFstDxy_XPEHH_DOEDOW$WinNum,   # Red high rrate+fst & xpehh
                highlight4 = real_high_diff_regions_DOEDOW$WinNum,   # Green ONLY outliers fst + rrate (rm xpehh)
                ylim = c(-2,20),col = c("grey40", "grey"), cex=1,
                ylab ="Fst \nDOW - DOE", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI(subset(FstDxy_pop,pop1==POP12_pop & pop2==POP2_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst  CAN DOW
                chr="ChrNum", bp="ChrPos", p="FstZscore",logp="FALSE", snp="WinNum",
                suggestiveline=FALSE,genomewideline=FALSE,
                highlight1 = FstDxy_XPEHH_CANDOW$WinNum,   # Yellow shared outliers xpehh + fst (no rrate)
                highlight3 = shared_outliersFstDxy_XPEHH_CANDOW$WinNum,   # Red high rrate+fst & xpehh
                highlight4 = real_high_diff_regions_CANDOW$WinNum,   # Green ONLY outliers fst + rrate (rm xpehh)
                ylim = c(-2,26),col = c("grey40", "grey"), cex=1,
                ylab ="Fst \nDOW - CAN", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI_2(subset(FstDxy_pop,pop1==POP12_pop & pop2==POP13_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst  CAN MEA
                chr="ChrNum", bp="ChrPos", p="FstZscore",logp="FALSE", snp="WinNum",
                suggestiveline=FALSE,genomewideline=FALSE,
                highlight1 = FstDxy_XPEHH_CANMEA$WinNum,   # Yellow shared outliers xpehh + fst (no rrate)
                highlight3 = shared_outliersFstDxy_XPEHH_CANMEA$WinNum,   # Red high rrate+fst & xpehh
                highlight4 = real_high_diff_regions_CANMEA$WinNum,   # Green ONLY outliers fst + rrate (rm xpehh)
                ylim = c(-2,26),col = c("grey40", "grey"), cex=1,
                ylab ="Fst \nCAN - MEA", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI_2(subset(FstDxy_pop,pop1==POP13_pop & pop2==POP9_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst  MEA MON
                chr="ChrNum", bp="ChrPos", p="FstZscore",logp="FALSE", snp="WinNum",
                suggestiveline=FALSE,genomewideline=FALSE,
                highlight1 = FstDxy_XPEHH_MEAMON$WinNum,   # Yellow shared outliers xpehh + fst (no rrate)
                highlight3 = shared_outliersFstDxy_XPEHH_MEAMON$WinNum,   # Red high rrate+fst & xpehh
                highlight4 = real_high_diff_regions_MEAMON$WinNum,   # Green ONLY outliers fst + rrate (rm xpehh)
                ylim = c(-2,26),col = c("grey40", "grey"), cex=1,
                ylab ="Fst \nMEA - MON", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI_2(subset(FstDxy_pop,pop1==POP9_pop & pop2==POP10_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst  MON ORE
                chr="ChrNum", bp="ChrPos", p="FstZscore",logp="FALSE", snp="WinNum",
                suggestiveline=FALSE,genomewideline=FALSE,
                highlight1 = FstDxy_XPEHH_MONORE$WinNum,   # Yellow shared outliers xpehh + fst (no rrate)
                highlight3 = shared_outliersFstDxy_XPEHH_MONORE$WinNum,   # Red high rrate+fst & xpehh
                highlight4 = real_high_diff_regions_MONORE$WinNum,   # Green ONLY outliers fst + rrate (rm xpehh)
                ylim = c(-2,26),col = c("grey40", "grey"), cex=1,
                ylab ="Fst \nMON - ORE", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI_2(subset(FstDxy_pop,pop1==POP9_pop & pop2==POP11_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst MON PIN
                chr="ChrNum", bp="ChrPos", p="FstZscore",logp="FALSE", snp="WinNum",
                suggestiveline=FALSE,genomewideline=FALSE,
                highlight1 = FstDxy_XPEHH_MONPIN$WinNum,   # Yellow shared outliers xpehh + fst (no rrate)
                highlight3 = shared_outliersFstDxy_XPEHH_MONPIN$WinNum,   # Red high rrate+fst & xpehh
                highlight4 = real_high_diff_regions_MONPIN$WinNum,   # Green ONLY outliers fst + rrate (rm xpehh)
                ylim = c(-2,20),col = c("grey40", "grey"), cex=1,
                ylab ="Fst \nMON - PIN", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI(subset(FstDxy_pop,pop1==POP9_pop & pop2==POP3_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst MON THU
                chr="ChrNum", bp="ChrPos", p="FstZscore",logp="FALSE", snp="WinNum",
                suggestiveline=FALSE,genomewideline=FALSE,
                highlight1 = FstDxy_XPEHH_MONTHU$WinNum,   # Yellow shared outliers xpehh + fst (no rrate)
                highlight3 = shared_outliersFstDxy_XPEHH_MONTHU$WinNum,   # Red high rrate+fst & xpehh
                highlight4 = real_high_diff_regions_MONTHU$WinNum,   # blue ONLY outliers fst + rrate (rm xpehh)
                ylim = c(-2,20),col = c("grey40", "grey"), cex=1,
                ylab ="Fst \nMON - THU", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI(subset(FstDxy_pop,pop1==POP3_pop & pop2==POP4_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst THU TOW
                chr="ChrNum", bp="ChrPos", p="FstZscore",logp="FALSE", snp="WinNum",
                suggestiveline=FALSE,genomewideline=FALSE,
                highlight1 = FstDxy_XPEHH_THUTOW$WinNum,   # Yellow shared outliers xpehh + fst (no rrate)
                highlight3 = shared_outliersFstDxy_XPEHH_THUTOW$WinNum,   # Red high rrate+fst & xpehh
                highlight4 = real_high_diff_regions_THUTOW$WinNum,   # blue ONLY outliers fst + rrate (rm xpehh)
                ylim = c(-2,20),col = c("grey40", "grey"), cex=1,
                ylab ="Fst \nTHU - TOW", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI_2(subset(FstDxy_pop,pop1==POP8_pop & pop2==POP13_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst AIK MEA
                chr="ChrNum", bp="ChrPos", p="FstZscore",logp="FALSE", snp="WinNum",
                suggestiveline=FALSE,genomewideline=FALSE,
                highlight1 = FstDxy_XPEHH_AIKMEA$WinNum,   # Yellow shared outliers xpehh + fst (no rrate)
                highlight3 = shared_outliersFstDxy_XPEHH_AIKMEA$WinNum,   # Red high rrate+fst & xpehh
                highlight4 = real_high_diff_regions_AIKMEA$WinNum,   # Green ONLY outliers fst + rrate (rm xpehh)
                ylim = c(-2,26),col = c("grey40", "grey"), cex=1,
                ylab ="Fst \nMEA - AIK", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI_2(subset(FstDxy_pop,pop1==POP8_pop & pop2==POP6_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst AIK HYE
                chr="ChrNum", bp="ChrPos", p="FstZscore",logp="FALSE", snp="WinNum",
                suggestiveline=FALSE,genomewideline=FALSE,
                highlight1 = FstDxy_XPEHH_AIKHYE$WinNum,   # Yellow shared outliers xpehh + fst (no rrate)
                highlight3 = shared_outliersFstDxy_XPEHH_AIKHYE$WinNum,   # Red high rrate+fst & xpehh
                highlight4 = real_high_diff_regions_AIKHYE$WinNum,   # Green ONLY outliers fst + rrate (rm xpehh)
                ylim = c(-2,26),col = c("grey40", "grey"), cex=1,
                ylab ="Fst \nAIK - HYE", chrlabs = ChrName[CHRS:CHRE])

manhattanJAVI(subset(FstDxy_pop,pop1==POP5_pop & pop2==POP6_pop & ChrNum>=CHRS & ChrNum<=CHRE),#Fst CAR HYE
                chr="ChrNum", bp="ChrPos", p="FstZscore",logp="FALSE", snp="WinNum",
                suggestiveline=FALSE,genomewideline=FALSE,
                highlight1 = FstDxy_XPEHH_CARHYE$WinNum,  # Green shared outliers xpehh + fst (no rrate)
                highlight3 = shared_outliersFstDxy_XPEHH_CARHYE$WinNum,   # Red high rrate+fst & xpehh
                highlight4 = real_high_diff_regions_CARHYE$WinNum,   # blue ONLY outliers fst + rrate (rm xpehh)
                ylim = c(-2,20),col = c("grey40", "grey"), cex=1,
                ylab ="Fst \nHYE - CAR", chrlabs = ChrName[CHRS:CHRE])
