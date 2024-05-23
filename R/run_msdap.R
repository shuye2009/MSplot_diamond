scriptdir <- dirname(rstudioapi::getSourceEditorContext()$path)

source(file.path(scriptdir,"msdap_lib.R"))

fastas <- c(human = "C:/GREENBLATT/resource/human_uniprot_reference_proteome_9606.fasta/UP000005640_9606_061820.fasta",
            sars2 = "C:/GREENBLATT/resource/sars2_uniprot_reference_proteome_2697049.fasta/REFSEQ_Wuhan_Hu_1_nr.fasta",
            ebv = "C:/GREENBLATT/resource/EB_virus_fasta/uniprot-proteome_UP000153037_strain_B95-8.fasta",
            HA = "C:/data/raw/EDYTA/PROTEIN/240117_RobertBonin_DIA/HA-C1.fasta",
            clostridium = "C:/GREENBLATT/resource/Clostridium_perfringens_fasta/UP000000818_195102_reference.fasta")

### 240410_EM_DIA ####
wd <- "C:/data/raw/EDYTA/PROTEIN/240410_EM_DIA/DIANN"

mydataset <- prepare_msdap(wd)

#comps = "group|C"
comps = "group|group_6"
contrasts <- make_contrasts(wd, comparisons = comps, allpairs = TRUE)

res <- run_msdap(wd, fastas[c(1,5)], mydataset, contrasts=contrasts)
res$dir
head(res$data)

wd <- file.path(wd,res$dir)

plot_volcano_for_msdap(wd, padj_cutoff = 0.05,
                       lfc_cutoffs = c("group_6.vs.group_X" = 1.812),
                       col_factor=NULL, gsea=TRUE)

# for the contrasts, the first of the pair is the background
#

## 240311_AL_DIA ####
wd <- "C:/data/raw/EDYTA/PROTEIN/240311_AL_DIA/DIANN"

mydataset <- prepare_msdap(wd)

#comps = "group|C"
comps = "group|Group1"
contrasts <- make_contrasts(wd, comparisons = comps, allpairs = TRUE)

res <- run_msdap(wd, fastas[c(1)], mydataset, contrasts=contrasts)
res$dir
head(res$data)

wd <- file.path(wd,res$dir)

plot_volcano_for_msdap(wd, padj_cutoff=0.05,
                       lfc_cutoffs = c("Group1.vs.Group2" = 2.142,
                                       "Group1.vs.Group3" = 2.189,
                                       "Group1.vs.Group4" = 2.089,
                                       "Group2.vs.Group3" = 2.06,
                                       "Group2.vs.Group4" = 2.098,
                                       "Group3.vs.Group4" = 1.92),
                       col_factor=NULL, gsea=TRUE)
## 240305_EM_DIA ####
wd <- "C:/data/raw/EDYTA/PROTEIN/240305_EM_DIA/DIANN"

mydataset <- prepare_msdap(wd)

#comps = "group|C"
comps = "group|Group1"
contrasts <- make_contrasts(wd, comparisons = comps, allpairs = TRUE)

res <- run_msdap(wd, fastas[c(1)], mydataset, contrasts=contrasts)
res$dir
head(res$data)

wd <- file.path(wd,res$dir)

plot_volcano_for_msdap(wd, padj_cutoff=0.05,
                       lfc_cutoffs = NULL,
                       col_factor=NULL, gsea=TRUE)
#


## 240209_SGC_AL_DIA ####
wd <- "C:/data/raw/EDYTA/PROTEIN/240209_SGC_AL_DIA/DIANN"

mydataset <- prepare_msdap(wd)

comps = "group|wt"
contrasts <- make_contrasts(wd, comparisons = comps, allpairs = TRUE)

res <- run_msdap(wd, fastas[c(1)], mydataset, contrasts=contrasts)
res$dir
head(res$data)

wd <- file.path(wd,res$dir)

plot_volcano_for_msdap(wd, padj_cutoff=0.05, col_factor=NULL, gsea=FALSE)

## 240117_RobertBonin_DIA ####
wd <- "C:/data/raw/EDYTA/PROTEIN/240117_RobertBonin_DIA/DIANN"
fastas <- c(fastas, file.path("C:/data/raw/EDYTA/PROTEIN/240117_RobertBonin_DIA", "HA-C1.fasta"))

mydataset <- prepare_msdap(wd)
comps = "group|pbs_wt"
contrasts <- make_contrasts(wd, comparisons = comps, allpairs = TRUE)

res <- run_msdap(wd, fastas[c(1,4)], mydataset, contrasts=contrasts)
res$dir
head(res$data)

wd <- file.path(wd,res$dir)
lfc_cutoffs <- c("pep1_c11.vs.pep2_c12" = 1.027,
                 "pbs_wt.vs.pep1_c11" = 0.977,
                 "pbs_wt.vs.pep2_c12" = 1.362)
plot_volcano_for_msdap(wd, padj_cutoff = 0.05, lfc_cutoffs = lfc_cutoffs,
                       col_factor = NULL, gsea = FALSE)

## 231215_EM_EBV_DIA ####
wd <- "C:/data/raw/EDYTA/PROTEIN/231215_EM_EBV_DIA/DIANN"

mydataset <- prepare_msdap(wd)
contrasts <- make_contrasts(wd, comparisons = "Treatment|UI,Time|0,Phenotype|MOCK")

res <- run_msdap(wd, fastas[c(1,3)], mydataset, contrasts=contrasts)
res$dir

wd <- file.path(wd,res$dir)
geneNames <- gene_names_from_fasta(fastas[3])
plot_for_msdap_drugTreatment(wd, geneNames, impute="halfmin", comparisons = "Treatment|UI,Time|0,Phenotype|MOCK", data_name="231215_EM_EBV_DIA")
plot_volcano_for_msdap(wd, padj_cutoff=0.05, col_factor=NULL, gsea=FALSE)

## 231206_EM_EBV #####
wd <- "C:/data/raw/EDYTA/PROTEIN/231206_EM_EBV/combined/txt"

mydataset <- prepare_msdap(wd, searchType = "MAXQUANT")
contrasts <- make_contrasts(wd, comparisons = "Treatment|UI,Time|0,Phenotype|MOCK")

res <- run_msdap(wd, fastas, mydataset, contrasts=contrasts, plot_data = FALSE)
res$dir

wd <- file.path(wd,res$dir)
geneNames <- gene_names_from_fasta(fastas[3])
plot_for_msdap_drugTreatment(wd, geneNames, impute="halfmin", comparisons = "Treatment|UI,Time|0,Phenotype|MOCK", data_name="231206_EM_EBV")
plot_volcano_for_msdap(wd, padj_cutoff=0.05, col_factor=NULL, gsea=FALSE)
#
## 221026_CALU3_SARS2_INF8_vcp_DIA ####
wd <- "C:/data/raw/EDYTA/PROTEIN/221026_CALU3_SARS2_INF8_vcp_DIA/DIANN"


mydataset <- prepare_msdap(wd)
contrasts <- make_contrasts(wd, treatRef="mock", groupRef="nd")

res <- run_msdap(wd, fastas, mydataset, contrasts=contrasts)
res$dir

wd <- file.path(wd,res$dir)
geneNames <- gene_names_from_fasta(fastas[2])
plot_for_msdap_drugTreatment(wd, geneNames, impute="halfmin", treatRef="mock", drugRef="nd")
plot_volcano_for_msdap(wd, padj_cutoff=0.05, col_factor=NULL, gsea=FALSE)

## 221017_calu3_sars2_inf8_vcp_dda #####
wd <- "C:/data/raw/EDYTA/PROTEIN/221017_calu3_sars2_inf8_vcp_dda/combined/txt"

mydataset <- prepare_msdap(wd, searchType = "MAXQUANT")
contrasts <- make_contrasts(wd, treatRef="mock", groupRef="nd")

res <- run_msdap(wd, fastas, mydataset, contrasts=contrasts)
res$dir

wd <- file.path(wd,res$dir)
geneNames <- gene_names_from_fasta(fastas[2])
plot_for_msdap_drugTreatment(wd, geneNames, impute="halfmin", treatRef="mock", drugRef="nd")
plot_volcano_for_msdap(wd, padj_cutoff=0.05, col_factor=NULL, gsea=TRUE)
#
## 220913_calu3_sars2_vcp_DDA #####
wd <- "C:/data/raw/EDYTA/PROTEIN/220913_calu3_sars2_vcp_DDA/combined/txt"

mydataset <- prepare_msdap(wd, searchType = "MAXQUANT")
contrasts <- make_contrasts(wd, treatRef="mock", groupRef="nd")

res <- run_msdap(wd, fastas, mydataset, contrasts=contrasts)
res$dir

wd <- file.path(wd,res$dir)
geneNames <- gene_names_from_fasta(fastas[2])
plot_for_msdap_drugTreatment(wd, geneNames, impute="halfmin", treatRef="mock", drugRef="nd")
plot_volcano_for_msdap(wd, padj_cutoff=0.05, col_factor=NULL, gsea=FALSE)
#
## 220913_calu3_sars2_vcp_DIA #####
wd <- "C:/data/raw/EDYTA/PROTEIN/220913_calu3_sars2_vcp_DIA/DIANN"

mydataset <- prepare_msdap(wd)
contrasts <- make_contrasts(wd, treatRef="mock", groupRef="nd")

res <- run_msdap(wd, fastas, mydataset, contrasts=contrasts)
res$dir

wd <- file.path(wd,res$dir)
geneNames <- gene_names_from_fasta(fastas[2])
plot_for_msdap_drugTreatment(wd, geneNames, impute="halfmin", treatRef="mock", drugRef="nd")
plot_volcano_for_msdap(wd, padj_cutoff=0.05, col_factor=NULL, gsea=FALSE)


## 220902_NASAL_basal_apical_day1&2_180min basal #####
wd <- "C:/data/raw/EDYTA/PROTEIN/220902_NASAL_basal_apical_day1&2_180min/basal_180_D1&D2/DIANN"

mydataset <- prepare_msdap(wd)
contrasts <- make_contrasts(wd, treatRef="mock", groupRef="dmso")

res <- run_msdap(wd, fastas, mydataset, contrasts=contrasts)
res$dir

wd <- file.path(wd,res$dir)
geneNames <- gene_names_from_fasta(fastas[2])
plot_for_msdap_drugTreatment(wd, geneNames, impute="halfmin", treatRef="mock", drugRef="dmso")
plot_volcano_for_msdap(wd, padj_cutoff=0.05, col_factor=NULL, gsea=FALSE)

## 220902_NASAL_basal_apical_day1&2_180min apical #####
wd <- "C:/data/raw/EDYTA/PROTEIN/220902_NASAL_basal_apical_day1&2_180min/apical_180_D1&D2/DIANN"

mydataset <- prepare_msdap(wd)
contrasts <- make_contrasts(wd, treatRef="mock", groupRef="dmso")

res <- run_msdap(wd, fastas, mydataset, contrasts=contrasts)
res$dir

wd <- file.path(wd,res$dir)
geneNames <- gene_names_from_fasta(fastas[2])
plot_for_msdap_drugTreatment(wd, geneNames, impute="halfmin", treatRef="mock", drugRef="dmso")
plot_volcano_for_msdap(wd, padj_cutoff=0.05, col_factor=NULL, gsea=FALSE)

## 220504_nasal_basal_apical_Day1_180min_apical DIANN #####
wd <- "C:/data/raw/EDYTA/PROTEIN/220504_nasal_basal_apical_Day1_180min/apical/DIANN"

mydataset <- prepare_msdap(wd)
contrasts <- make_contrasts(wd, treatRef="mock", groupRef="dmso")

res <- run_msdap(wd, fastas, mydataset, contrasts=contrasts)
res$dir

wd <- file.path(wd,res$dir)
geneNames <- gene_names_from_fasta(fastas[2])
plot_for_msdap_drugTreatment(wd, geneNames, impute="halfmin", treatRef="mock", drugRef="dmso")
plot_volcano_for_msdap(wd, padj_cutoff=0.05, col_factor=NULL, gsea=FALSE)

## 220504_nasal_basal_apical_Day1_180min_base DIANN #####
wd <- "C:/data/raw/EDYTA/PROTEIN/220504_nasal_basal_apical_Day1_180min/basal/DIANN"

mydataset <- prepare_msdap(wd)
contrasts <- make_contrasts(wd, treatRef="mock", groupRef="dmso")

res <- run_msdap(wd, fastas, mydataset, contrasts=contrasts)
res$dir

wd <- file.path(wd,res$dir)
geneNames <- gene_names_from_fasta(fastas[2])
plot_for_msdap_drugTreatment(wd, geneNames, impute="halfmin", treatRef="mock", drugRef="dmso")
plot_volcano_for_msdap(wd, padj_cutoff=0.05, col_factor=NULL, gsea=FALSE)

## 220804_nasal_basal_apical_Day1_90min_basal_Day1 DIANN #####
wd <- "C:/data/raw/EDYTA/PROTEIN/220804_nasal_basal_apical_Day1_90min/basal_Day1/DIANN"

mydataset <- prepare_msdap(wd)
contrasts <- make_contrasts(wd, treatRef="mock", groupRef="dmso")

res <- run_msdap(wd, fastas, mydataset, contrasts=contrasts)
res$dir

wd <- file.path(wd,res$dir)
geneNames <- gene_names_from_fasta(fastas[2])
plot_for_msdap_drugTreatment(wd, geneNames, impute="halfmin", treatRef="mock", drugRef="dmso")
plot_volcano_for_msdap(wd, padj_cutoff=0.05, col_factor=NULL, gsea=FALSE)

## 220804_nasal_basal_apical_Day1_90min_apical_Day1 DIANN #####
wd <- "C:/data/raw/EDYTA/PROTEIN/220804_nasal_basal_apical_Day1_90min/apical_Day1/DIANN"

mydataset <- prepare_msdap(wd)
contrasts <- make_contrasts(wd, treatRef="mock", groupRef="dmso")

res <- run_msdap(wd, fastas, mydataset, contrasts=contrasts)
res$dir

wd <- file.path(wd,res$dir)
geneNames <- gene_names_from_fasta(fastas[2])
plot_for_msdap_drugTreatment(wd, geneNames, impute="halfmin", treatRef="mock", drugRef="dmso")
plot_volcano_for_msdap(wd, padj_cutoff=0.05, col_factor=NULL, gsea=FALSE)

## 220705_nasal_basal_EXP3 DIANN #####
wd <- "C:/data/raw/EDYTA/PROTEIN/220705_nasal_basal_EXP3/apical_exp3/DIANN"

mydataset <- prepare_msdap(wd)
contrasts <- make_contrasts(wd, treatRef="mock", groupRef="dmso")

res <- run_msdap(wd, fastas, mydataset, contrasts=contrasts)
res$dir

wd <- file.path(wd,res$dir)
geneNames <- gene_names_from_fasta(fastas[2])
plot_for_msdap_drugTreatment(wd, geneNames, impute="halfmin", treatRef="mock", drugRef="dmso")

### apical DIANN_D2 ####
wd <- "C:/data/raw/EDYTA/PROTEIN/220705_nasal_basal_EXP3/apical_exp3/DIANN_D2"


mydataset <- prepare_msdap(wd)
contrasts <- make_contrasts(wd, treatRef="mock", groupRef="dmso")

res <- run_msdap(wd, fastas, mydataset, contrasts=contrasts)
res$dir

wd <- file.path(wd,res$dir)
geneNames <- gene_names_from_fasta(fastas[2])
plot_for_msdap_drugTreatment(wd, geneNames, impute="halfmin", treatRef="mock", drugRef="dmso")
### apical DIANN_D3 ####
wd <- "C:/data/raw/EDYTA/PROTEIN/220705_nasal_basal_EXP3/apical_exp3/DIANN_D3"

# the first of the pair is the background

mydataset <- prepare_msdap(wd)
contrasts <- make_contrasts(wd, treatRef="mock", groupRef="dmso")

res <- run_msdap(wd, fastas, mydataset, contrasts=contrasts)
res$dir

wd <- file.path(wd,res$dir)
geneNames <- gene_names_from_fasta(fastas[2])
plot_for_msdap_drugTreatment(wd, geneNames, impute="halfmin", treatRef="mock", drugRef="dmso")

### basal DIANN_D2 ####
wd <- "C:/data/raw/EDYTA/PROTEIN/220705_nasal_basal_EXP3/basal_exp3/DIANN_D2"

mydataset <- prepare_msdap(wd)
contrasts <- make_contrasts(wd, treatRef="mock", groupRef="dmso")

res <- run_msdap(wd, fastas, mydataset, contrasts=contrasts)
res$dir

wd <- file.path(wd,res$dir)
geneNames <- gene_names_from_fasta(fastas[2])
plot_for_msdap_drugTreatment(wd, geneNames, impute="halfmin", treatRef="mock", drugRef="dmso")

### basal DIANN_D3 ####
wd <- "C:/data/raw/EDYTA/PROTEIN/220705_nasal_basal_EXP3/basal_exp3/DIANN_D3"

mydataset <- prepare_msdap(wd)
contrasts <- make_contrasts(wd, treatRef="mock", groupRef="dmso")

res <- run_msdap(wd, fastas, mydataset, contrasts=contrasts)
res$dir

wd <- file.path(wd,res$dir)
geneNames <- gene_names_from_fasta(fastas[2])
plot_for_msdap_drugTreatment(wd, geneNames, impute="halfmin", treatRef="mock", drugRef="dmso")

## 220516_INF2_INF5_DIA DIA-NN #####
wd <- "C:/data/raw/EDYTA/PROTEIN/220516_INF2_INF5_DIA"

preprocess_diann_report(wd)
mydataset <- prepare_msdap(wd)

contrasts <- list(c("sars2_nd", "sars2_b5"), c("mock_b5", "sars2_b5"), c("mock_nd", "sars2_nd"), c("mock_nd", "mock_b5"))

contrasts <- make_contrasts(wd, treatRef="mock", groupRef="nd")
contrasts <- contrasts[]

last_dir <- run_msdap(wd, fastas, mydataset, contrasts=contrasts)
last_dir <- "2022-05-18_10;39;20"

wd <- file.path(wd,last_dir)

geneNames <- gene_names_from_fasta(fastas[2])
plot_for_msdap(wd, geneNames)
col_order <- c("mock_nd.vs.mock_b1", "mock_nd.vs.mock_b5", "mock_nd.vs.mock_b10", "mock_nd.vs.mock_b50",
               "mock_nd.vs.sars2_nd", "mock_b1.vs.sars2_b1", "mock_b5.vs.sars2_b5", "mock_b10.vs.sars2_b10", "mock_b50.vs.sars2_b50",
               "sars2_nd.vs.sars2_b1", "sars2_nd.vs.sars2_b5", "sars2_nd.vs.sars2_b10", "sars2_nd.vs.sars2_b50")
col_factor <- factor(c(rep("Mock",4), rep("MockSars2",5), rep("Sars2", 4)), levels=c("Mock", "MockSars2", "Sars2"))
names(col_factor) <- col_order
gsea_for_msdap(wd, padj_cutoff = 0.2, col_factor=col_factor)

### 220324_basal DIA-NN #####
wd <- "C:/data/raw/EDYTA/PROTEIN/220324_nasal_EXP1/220324_basal/DIANN_05_04_2022"

# the first of the pair is the background

mydataset <- prepare_msdap(wd)
contrasts <- make_contrasts(wd, treatRef="mock", groupRef="dmso")
contrasts <- contrasts[6:7]

last_dir <- run_msdap(wd, fastas, mydataset, contrasts=contrasts)
last_dir

wd <- file.path(wd,last_dir)
geneNames <- gene_names_from_fasta(fastas[2])
plot_for_msdap(wd, geneNames)


### 220324_apical DIA-NN #####
wd <- "C:/data/raw/EDYTA/PROTEIN/220324_nasal_EXP1/220324_apical/DIANN_04_04_2022"

# the first of the pair is the background

mydataset <- prepare_msdap(wd)
contrasts <- make_contrasts(wd, treatRef="mock", groupRef="dmso")

last_dir <- run_msdap(wd, fastas, mydataset, contrasts=contrasts)
last_dir

wd <- file.path(wd,last_dir)
geneNames <- gene_names_from_fasta(fastas[2])
plot_for_msdap(wd, geneNames)

### 220222_INF5_DIA DIA-NN #####
wd <- "C:/data/raw/EDYTA/PROTEIN/220222_INF5_DIA/DIANN_02_03_22"
# the first of the pair is the background
contrasts <- list(c("sars2_nd", "sars2_b5"), c("mock_b5", "sars2_b5"), c("mock_nd", "sars2_nd"), c("mock_nd", "mock_b5"))

mydataset <- prepare_msdap(wd)

run_msdap(wd, fastas, mydataset, contrasts=contrasts)

wd <- "C:/data/raw/EDYTA/PROTEIN/220222_INF5_DIA/DIANN_02_03_22/2022-03-03_16;24;58"
gsea_for_msdap(wd, padj_cutoff = 0.2)
geneNames <- gene_names_from_fasta(fastas[2])
plot_for_msdap(wd, geneNames)

### 220222_INF5_DIA mzDIA-NN #####
wd <- "C:/data/raw/EDYTA/PROTEIN/220222_INF5_DIA/mzDIANN_21_03_2022"
# the first of the pair is the background
contrasts <- list(c("sars2_nd", "sars2_b5"), c("mock_b5", "sars2_b5"), c("mock_nd", "sars2_nd"), c("mock_nd", "mock_b5"))

mydataset <- prepare_msdap(wd)

out_dir <- run_msdap(wd, fastas, mydataset, contrasts=contrasts)

wd <- file.path("C:/data/raw/EDYTA/PROTEIN/220222_INF5_DIA/mzDIANN_21_03_2022", out_dir)

geneNames <- gene_names_from_fasta(fastas[2])
plot_for_msdap(wd, geneNames)
gsea_for_msdap(wd)

### 220222_INF5_DIA FraggerPipe DIA-NN #####
wd <- "C:/data/raw/EDYTA/PROTEIN/220222_INF5_DIA/FragPipe_17_03_2022"
# the first of the pair is the background
contrasts <- list(c("sars2_nd", "sars2_b5"), c("mock_b5", "sars2_b5"), c("mock_nd", "sars2_nd"), c("mock_nd", "mock_b5"))

mydataset <- prepare_msdap(wd, searchType = "msfragger")

run_msdap(wd, fastas, mydataset, contrasts=contrasts)

wd <- "C:/data/raw/EDYTA/PROTEIN/220222_INF5_DIA/FragPipe_17_03_2022/2022-03-18_08;46;07"
gsea_for_msdap(wd)
geneNames <- gene_names_from_fasta(fastas[2])
plot_for_msdap(wd, geneNames)

### 220222_INF5_DIA FraggerPipe Umpire DIA-NN #####
wd <- "C:/data/raw/EDYTA/PROTEIN/220222_INF5_DIA/FragPipe_Umpire_18_03_2022"
# the first of the pair is the background
contrasts <- list(c("sars2_nd", "sars2_b5"), c("mock_b5", "sars2_b5"), c("mock_nd", "sars2_nd"), c("mock_nd", "mock_b5"))

mydataset <- prepare_msdap(wd, searchType = "DIA-NN")

run_msdap(wd, fastas, mydataset, contrasts=contrasts)

wd <- "C:/data/raw/EDYTA/PROTEIN/220222_INF5_DIA/FragPipe_Umpire_18_03_2022/2022-03-18_18;18;30"
gsea_for_msdap(wd)
geneNames <- gene_names_from_fasta(fastas[2])
plot_for_msdap(wd, geneNames)

### 220222_INF5_DIA Maxquant DIA#####
wd <- "C:/data/raw/EDYTA/PROTEIN/220222_INF5_DIA/combined/txt"
# the first of the pair is the background
contrasts <- list(c("sars2_nd", "sars2_b5"), c("mock_b5", "sars2_b5"), c("mock_nd", "sars2_nd"), c("mock_nd", "mock_b5"))

mydataset <- prepare_msdap(wd, searchType = "maxquant")

run_msdap(wd, fastas, mydataset, contrasts=contrasts)

wd <- "C:/data/raw/EDYTA/PROTEIN/220222_INF5_DIA/combined/txt/2022-03-07_09;53;49"
gsea_for_msdap(wd)
plot_for_msdap(wd, geneNames)

### 211221_INF5_2 Maxquant DDA #####
wd <- "C:/data/raw/EDYTA/PROTEIN/211221_INF5_2/combined/txt"
# the first of the pair is the background
contrasts <- list(c("sars2_nd", "sars2_b5"), c("mock_b5", "sars2_b5"), c("mock_nd", "sars2_nd"), c("mock_nd", "mock_b5"))

mydataset <- prepare_msdap(wd, searchType = "maxquant")

run_msdap(wd, fastas, mydataset, contrasts=contrasts)

wd <- "C:/data/raw/EDYTA/PROTEIN/211221_INF5_2/combined/txt/2022-03-07_10;18;19"
gsea_for_msdap(wd)
plot_for_msdap(wd, geneNames)

### 220301_INF2_DIA DIA-NN #####
wd <- "C:/data/raw/EDYTA/PROTEIN/220301_INF2_DIA/DIANN_10_03_2022"
mydataset <- prepare_msdap(wd)
mydataset$acquisition_mode

contrasts <- list(c("inf_nd", "inf_b5"),
                  c("inf_nd", "inf_b10"),
                  c("inf_nd", "inf_b50"),
                  c("mock_nd", "mock_b5"),
                  c("mock_nd", "mock_b10"),
                  c("mock_nd", "mock_b50"),
                  c("mock_nd", "inf_nd"),
                  c("mock_b5", "inf_b5"),
                  c("mock_b10", "inf_b10"),
                  c("mock_b50", "inf_b50"))
run_msdap(wd, fastas, mydataset, contrasts=contrasts)


wd <- "C:/data/raw/EDYTA/PROTEIN/220301_INF2_DIA/DIANN_10_03_2022/2022-03-11_14;06;19"
plot_for_msdap(wd, geneNames)
gsea_for_msdap(wd)

### 210727_Calu3_drugs_INF2_regular #####
wd <- "C:/data/raw/EDYTA/PROTEIN/210727_Calu3_drugs_INF2_regular/combined/txt"
mydataset <- prepare_msdap(wd, searchType = "maxquant")
mydataset$acquisition_mode

contrasts <- list(c("inf_nd", "inf_b5"),
                  c("inf_nd", "inf_b10"),
                  c("inf_nd", "inf_b50"),
                  c("mock_nd", "mock_b5"),
                  c("mock_nd", "mock_b10"),
                  c("mock_nd", "mock_b50"),
                  c("mock_nd", "inf_nd"),
                  c("mock_b5", "inf_b5"),
                  c("mock_b10", "inf_b10"),
                  c("mock_b50", "inf_b50"))
run_msdap(wd, fastas, mydataset, contrasts=contrasts)


wd <- "C:/data/raw/EDYTA/PROTEIN/210727_Calu3_drugs_INF2_regular/combined/txt/2022-03-11_14;57;05"
plot_for_msdap(wd, geneNames)
gsea_for_msdap(wd)

### 210929_SARS2_INF4 Maxquant DDA #####
wd <- "C:/data/raw/EDYTA/PROTEIN/210929_SARS2_INF4/combined/txt"
# the first of the pair is the background
contrasts <- list(c("sars2_nd", "sars2_mtx"),
                  c("sars2_nd", "sars2_ri1"),
                  c("sars2_nd", "sars2_tg"),
                  c("mock_nd", "mock_mtx"),
                  c("mock_nd", "mock_ri1"),
                  c("mock_nd", "mock_tg"),
                  c("mock_nd", "sars2_nd"),
                  c("mock_mtx", "sars2_mtx"),
                  c("mock_ri1", "sars2_ri1"),
                  c("mock_tg", "sars2_tg")
)

mydataset <- prepare_msdap(wd, searchType = "maxquant")

run_msdap(wd, fastas, mydataset, contrasts=contrasts)

wd <- "C:/data/raw/EDYTA/PROTEIN/210929_SARS2_INF4/combined/txt/2022-03-10_14;15;02"
plot_for_msdap(wd, geneNames)

### 220310_INF4_DIA DIA-NN #####
wd <- "C:/data/raw/EDYTA/PROTEIN/220310_INF4_DIA/DIANN_11_03_2022"
# the first of the pair is the background
contrasts <- make_contrasts(wd)
mydataset <- prepare_msdap(wd)

last_dir <- run_msdap(wd, fastas, mydataset, contrasts=contrasts)
last_dir

wd <- "C:/data/raw/EDYTA/PROTEIN/220310_INF4_DIA/DIANN_11_03_2022/last_dir"
plot_for_msdap(wd, geneNames)

### 220307_INF1_DIA DIA-NN #####
wd <- "C:/data/raw/EDYTA/PROTEIN/220307_INF1_DIA/DIANN_10_03_2022"
# the first of the pair is the background
mydataset <- prepare_msdap(wd)
contrasts <- make_contrasts(wd)

run_msdap(wd, fastas, mydataset, contrasts=contrasts)

wd <- "C:/data/raw/EDYTA/PROTEIN/220307_INF1_DIA/DIANN_10_03_2022/2022-03-12_12;31;47"
plot_for_msdap(wd, geneNames)

### 210617_SARS2_DRUGS_EXP1 DDA Maxquant #####
wd <- "C:/data/raw/EDYTA/PROTEIN/210617_SARS2_DRUGS_EXP1/combined/txt"
# the first of the pair is the background
mydataset <- prepare_msdap(wd, searchType = "maxquant")

contrasts <- make_contrasts(wd)

run_msdap(wd, fastas, mydataset, contrasts=contrasts)

wd <- "C:/data/raw/EDYTA/PROTEIN/210617_SARS2_DRUGS_EXP1/combined/txt/2022-03-12_15;14;23"
plot_for_msdap(wd, geneNames)

### 210803_Calu3_drugs_Inf3 DDA Maxquant #####
wd <- "C:/data/raw/EDYTA/PROTEIN/210803_Calu3_drugs_Inf3/combined/txt"
# the first of the pair is the background
mydataset <- prepare_msdap(wd, searchType = "maxquant")

contrasts <- make_contrasts(wd, treatRef="mock", groupRef="aavs12")

run_msdap(wd, fastas, mydataset, contrasts=contrasts)

wd <- "C:/data/raw/EDYTA/PROTEIN/210803_Calu3_drugs_Inf3/combined/txt/2022-03-14_15;08;29"
plot_for_msdap(wd, geneNames, groupRef = "aavs12")

### 201105_service_Meena Maxquant #####
wd <- "C:/data/raw/EDYTA/PROTEIN/201105_service_Meena/combined/txt"
fastas <- c("C:/GREENBLATT/resource/uniprot-proteome_BolivianSquirrelMonkey.fasta/uniprot-proteome_UP000233220.fasta")

mydataset <- prepare_msdap(wd, searchType = "MAXQUANT")
# for the contrasts, the first of the pair is the background
contrasts <- list(c("Control", "THC"), c("Control", "THC_CBD"), c("THC", "THC_CBD"))

results <- run_msdap(wd, fastas, mydataset, contrasts=contrasts, plot_data = T)
last_dir <- results$dir
mydataset <- results$data

lwd <- file.path(wd,last_dir)

geneNames <- gene_names_from_fasta(fastas[1])
geneNames <- geneNames[geneNames %in% c("GFAP", "STMN1", "NRCAM", "FBXO2")]
plot_for_msdap_group(lwd, geneNames, groupRef="Control")
plot_peptide_data_group(lwd, mydataset, geneNames)
plot_volcano_for_msdap(lwd, padj_cutoff = 0.2, col_factor=NULL)
