preprocess_diann_report <- function(wd){
   df <- read.delim(file.path(wd, "report.tsv"))
   head(df)
   folder <- sapply(df$File.Name, function(x){
      unlist(strsplit(x, split="\\", fixed=T))[6]})
   names(folder) <- NULL
   df$Run <- paste(folder, df$Run, sep="_")
   df$File.Name <- paste(df$File.Name, folder, sep="_")
   write.table(df, file.path(wd, "report.tsv"), sep="\t", row.names=F, quote=F)
}

# Produce samples.xlsx file with these columns:
#
# sample_id, shortname, Treatment, Time, Phenotype, Cell, group, exclude
#
# edit this file manually according to your experimental design before proceeding
# with msdap analysis
#
prepare_msdap <- function(wd, overwrite=FALSE, searchType="DIA-NN"){
   library(iq)
   library(msdap)
   library(openxlsx)

   searchType <- toupper(searchType)
   if(searchType == "DIA-NN"){
      mydataset <- import_dataset_diann(file.path(wd, "report.tsv"))
   }else if(searchType == "MSFRAGGER"){
      mydataset <- import_dataset_diann(file.path(wd, "diann-output.tsv"))
   }else if(searchType == "MAXQUANT"){
      mydataset <- import_dataset_maxquant_evidencetxt(wd)
   }


   if(!file.exists(file.path(wd, "samples.xlsx")) | overwrite){
      write_template_for_sample_metadata(mydataset, file.path(wd, "samples.xlsx"), TRUE)
      ## make groups out of shortNames
      xlsxfile <- file.path(wd, "samples.xlsx")
      sheets <- openxlsx::getSheetNames(xlsxfile)
      annot <- lapply(sheets, openxlsx::read.xlsx, xlsxFile=xlsxfile)[[1]]
      instruc <- lapply(sheets, openxlsx::read.xlsx, xlsxFile=xlsxfile)[[2]]

      annot <- annot %>%
         dplyr::mutate(Treatment	= shortname,
                       Time = shortname,
                       Phenotype = shortname,
                       Cell = shortname,
                       group = shortname) %>%
         dplyr::select(sample_id, shortname, Treatment, Time, Phenotype, Cell,
                       group, exclude)


      print(annot)

      dataset_names <- list("samples" = annot, "instructions" = instruc)
      write.xlsx(dataset_names, file = file.path(wd, "samples.xlsx"))

      print(paste("!!! Edit", file.path(wd, "samples.xlsx"), "now"))
   }else{
      print(paste(file.path(wd, "samples.xlsx"), "exists already!"))
      print("If you want to make a new one, set 'overwrite' to TRUE and rerun this function!")
   }


   return(mydataset)
}

# Run functions in msdap R package
run_msdap <- function(wd, fastas, mydataset, contrasts=NULL, plot_data=FALSE){

   mydataset <- import_fasta(mydataset, files=fastas)
   mydataset = remove_proteins_by_name(
      mydataset,
      #remove keratins and IGGs using regular expression against uniprot fasta headers
      #(particularly useful for IP experiments);
      regular_expression = "keratin|GN=(krt|try|igk|igg|igkv|ighv|ighg)",
      remove_irt_peptides = FALSE,
      gene_symbols = c("DMKN", "ALB") ## in maxquant contaminants list
   )
   mydataset <- import_sample_metadata(mydataset, file.path(wd, "samples.xlsx"))

   if(is.null(contrasts)){
      mydataset <- setup_contrasts(mydataset, contrast_list = combn(unique(mydataset$samples$group), 2, simplify=F))
   }else{
      mydataset <- setup_contrasts(mydataset, contrast_list = contrasts)
   }


   names(mydataset)
   print_dataset_summary(mydataset)

   mydataset <- analysis_quickstart(mydataset,
                       filter_min_quant=1, #zero to fully rely on MBR
                       filter_by_contrast=T,
                       dea_log2foldchange_threshold=NA,
                       dea_algorithm = c("deqms"),
                       dea_qvalue_threshold = 0.05,
                       pca_sample_labels="shortname",
                       norm_algorithm = c("vwmb", "modebetween_protein"),
                       output_dir=wd)

   df <- file.info(list.files(wd, full.names = T))
   last_dir <- basename(rownames(df)[which.max(df$mtime)])

   if(plot_data){
      plot_peptide_data(mydataset,
                        select_all_proteins = TRUE,
                        select_diffdetect_candidates = TRUE,
                        select_dea_signif = TRUE,
                        output_dir = file.path(wd, last_dir),
                        show_unused_datapoints = TRUE,
                        norm_algorithm = c("vwmb", "modebetween_protein"))
   }

   print_dataset_summary(mydataset)

   return(list("dir"=last_dir, "data"=mydataset))
}

# Plot enhanced volcano and perform GSEA analysis
plot_volcano_for_msdap <- function(wd, padj_cutoff=0.05, lfc_cutoffs=NULL,
                                   col_factor=NULL, gsea=FALSE){
   #wd <- "C:/data/raw/EDYTA/PROTEIN/211221_INF5_2/combined/txt/2022-03-07_10;18;19"
   library(openxlsx)
   library(EnhancedVolcano)
   if(gsea){
      source("C:/GREENBLATT/Rscripts/RNAseq/Function_analysis_for_RNAseq_lib.R")
      PATHWAY_file <- "C:/GREENBLATT/resource/GSEA_gmt/ReactomePathways.gmt"
   }

   setwd(wd)
   xlsxfile <- file.path(wd, "differential_abundance_analysis.xlsx")

   # getting data from sheets
   sheets <- openxlsx::getSheetNames(xlsxfile)
   data_frame <- lapply(sheets, openxlsx::read.xlsx, xlsxFile=xlsxfile)

   # assigning names to data frame
   names(data_frame) <- sheets

   # printing the data
   # print (data_frame)
   my_df <- data_frame[["statistics"]]

   head(my_df)


   lfc_df <- my_df %>%
      select(c("gene_symbols_or_id", starts_with("foldchange.log2_deqms_contrast"), starts_with("qvalue_deqms_contrast"))) %>%
      filter(!grepl(";", gene_symbols_or_id)) %>% # filter out protein groups with more than two proteins
      filter(!duplicated(gene_symbols_or_id)) # filter out protein groups with duplicated names
   lfc_mat <- lfc_df %>%
      select(starts_with("foldchange.log2_deqms_contrast"))
   rownames(lfc_mat) <- lfc_df$gene_symbols_or_id
   colnames(lfc_mat) <- gsub("foldchange.log2_deqms_contrast:.", "", colnames(lfc_mat))

   ## choose gene with less than 25% missing values in all contrasts
   lfc_mat <- lfc_mat[apply(lfc_mat,1, function(x)sum(is.na(x))<0.25*ncol(lfc_mat)),]

   if(!is.null(col_factor)){
      lfc_mat <- lfc_mat[,names(col_factor)]
      cols <- rep("#c6dcff", length(col_factor))
      names(cols) <- col_factor
      ha <- ComplexHeatmap::HeatmapAnnotation(df = data.frame(treat = col_factor), col=list(treat=cols), which="column", show_legend=FALSE)
   }else{
      ha <- NULL
   }

   pdf("Volcano_plots.pdf", width=10, height=8)

   h <- ComplexHeatmap::Heatmap(lfc_mat,
                                #col=circlize::colorRamp2(c(-3, 0, 3), c("green", "black", "red"), space="RGB"),
                                col = circlize::colorRamp2(quantile(lfc_mat, seq(0, 1, by = 0.25), na.rm=TRUE), viridis::viridis(5)),
                                cluster_columns = F,
                                top_annotation = ha,
                                column_names_side = "bottom",
                                show_row_names = F,
                                column_split = col_factor,
                                name="fold change")
   ComplexHeatmap::draw(h)

   lapply(colnames(lfc_df)[grepl("foldchange.log2_deqms_contrast", colnames(lfc_df))], function(x){
      print(x)
      cn <- gsub("foldchange.log2_deqms_contrast:.", "", x)
      print(cn)

      ## for volcano plot
      y <- paste0("qvalue_deqms_contrast:.", cn)
      if(is.null(lfc_cutoffs)){
         lfc_cutoff <- 1
      }else{
         lfc_cutoff <- lfc_cutoffs[cn]
      }

      p <- EnhancedVolcano(lfc_df,
                           lab = lfc_df[,"gene_symbols_or_id"],
                           x = x,
                           y = y,
                           title = cn,
                           subtitle = paste("N =", nrow(lfc_df)),
                           caption = paste("log2 fold change cutoff = +/-", lfc_cutoff,  " adjusted pvalue cutoff = ", padj_cutoff, sep=""),
                           pCutoff = padj_cutoff,
                           FCcutoff = lfc_cutoff,
                           ylim = c(0, max(-log10(lfc_df[[y]]), na.rm = TRUE) + 1.5),
                           labSize = 4,
                           pointSize = 3.0,
                           #drawConnectors = TRUE,
                           widthConnectors = 0.75)
      print(p)

      ## for GSEA
      if(gsea){
         lfc_list <- lfc_df[, x]
         names(lfc_list) <- lfc_df[,1]
         lfc_list <- lfc_list[!is.na(lfc_list)]

         try(run_gseGO_simpleList(gene_list=lfc_list, dataName=cn, adjp_cutoff=padj_cutoff, ont="REACTOME_Pathway", GO_file=PATHWAY_file))
         try(run_gseGO_simpleList(gene_list=lfc_list, dataName=cn, adjp_cutoff=padj_cutoff))
      }

   })
   dev.off()
}

gene_names_from_fasta <- function(fasta){
   uniprot_fasta <- readLines(fasta)
   headers <- uniprot_fasta[grep(">", uniprot_fasta)]
   geneNames <- sapply(headers, function(header){
      pattern <- "GN=.+?\\s"
      uniprot_name <- unlist(lapply(header, function(x){y <- regmatches(x, regexpr(pattern, x)); if(length(y) == 0) y <- "-"; y}))
      uniprot_name <- unlist(lapply(uniprot_name, function(x){gsub("GN=", "", x)}))
      uniprot_name <- unlist(lapply(uniprot_name, function(x){trimws(x)}))
   })

   geneNames

}

#' @title plot two factors boxplot
#' @description This function extracts protein abundance values from output of 'run_msdap' for selected genes, plot
#' boxplots for treat factor and drug factor separately, with statistical significance annotations. The statistical analysis
#' is one-way anova followed by t-test on log transformed intensity values.
#'
#' @param wd working directory where sample information file 'sample.xlsx', and abundance data
#' 'protein_abundance__input data as-is.tsv' are located
#' @param geneNames a vector of characters denoting the selected gene names whose data will be analyzed
#' @param impute a string in ("halfmin", "wrproteo") denoting the imputation methods for missing values, default 'halfmin'
#' @param comparisons reference group for each factor
#'
#' @return NULL
#' @author Shuye Pu

plot_for_msdap_drugTreatment <- function(wd, geneNames, impute="halfmin", comparisons="Treatment|UI,Time|0,Group|MOCK", data_name="msdap"){

   library(ComplexHeatmap)
   library(openxlsx)
   library(ggpubr)
   library(wrProteo)

   setwd(wd)
   xlsxfile <- file.path(wd, "samples.xlsx")
   sheets <- openxlsx::getSheetNames(xlsxfile)
   annot <- lapply(sheets, openxlsx::read.xlsx, xlsxFile=xlsxfile)[[1]]

   NAME <- colnames(annot)[2]
   EXPERIMENT <- colnames(annot)[3]
   TREATMENT <- colnames(annot)[4]
   DOSE <- colnames(annot)[5]
   PHENOTYPE <- colnames(annot)[6]
   CELL <- colnames(annot)[7]
   GROUP <- colnames(annot)[9]

   categs <- unlist(strsplit(comparisons, split=","))
   categ <- unlist(lapply(categs, function(x) unlist(strsplit(x, split="|", fixed=T))[1]))
   categRef <- unlist(lapply(categs, function(x) unlist(strsplit(x, split="|", fixed=T))[2]))
   names(categRef) <- tolower(categ)


   tsvfile <- file.path(wd, "protein_abundance__input data as-is.tsv")

   df_all <- read.delim2(tsvfile)
   if(length(intersect(df_all$gene_symbols_or_id, toupper(geneNames))) == 0){
      message("No target protein is found in the protein aboundance data!")
      return()
   }
   mat <- NULL

   if(impute=="wrproteo"){
      dfImp <- matrixNAneighbourImpute(type.convert(as.matrix(df_all[, 4:ncol(df_all)]), as.is=T), annot[, CELL])
      df <- dfImp$data
      df <- df[df_all$fasta_headers %in% names(geneNames),]
      rownames(df) <- df_all$gene_symbols_or_id[df_all$fasta_headers %in% names(geneNames)]
      mat <- t(df)
   }else{
      df <- df_all[df_all$fasta_headers %in% names(geneNames),]
      rownames(df) <- df_all$gene_symbols_or_id[df_all$fasta_headers %in% names(geneNames)]
      mat <- t(type.convert(as.matrix(df[, 4:ncol(df)]), as.is=T))
      NAs <- length(mat[is.na(mat)])
      MIN <- min(mat[!is.na(mat)])/2
      SD <- sqrt(sd(mat[!is.na(mat)]))
      mat[is.na(mat)] <- rnorm(n=NAs, mean=MIN, sd=SD)
      print(paste("half of min(nonNA values)", MIN, "square root of sd(nonNA values)", SD))
   }

   rownames(mat) <- annot[, EXPERIMENT]

   p <- ComplexHeatmap::Heatmap(mat, name="log(intensity)", cluster_rows=T, cluster_columns=T,
                           column_names_rot = 45,
                           row_split = 3, column_split=3,
                           right_annotation = c(rowAnnotation(treat=annot[,TREATMENT]),
                                                rowAnnotation(time=as.factor(annot[,DOSE])),
                                                rowAnnotation(phenotype=annot[,PHENOTYPE])),
                           cell_fun = function(j, i, x, y, width, height, fill) {
                                 grid.text(sprintf("%.1f", mat[i, j]), x, y, gp = gpar(fontsize = 8))}
                           )
   pdf("heatmap_of_selected_proteins.pdf", width=10, height=8)
   draw(p)
   dev.off()

   df_wide <- as.data.frame(cbind(annot[, c(EXPERIMENT, TREATMENT, DOSE, PHENOTYPE, CELL)], mat))
   df_wide[,DOSE] <- as.factor(df_wide[,DOSE])
   stats_df <- pivot_longer(df_wide, cols=colnames(mat), names_to = "Protein", values_to = "Log2Intensity") %>%
      mutate(Log2Intensity=as.numeric(Log2Intensity))

   for(cell in unique(annot[, CELL])){
      cell_df <- stats_df[stats_df[,CELL] == cell,]
      if(nrow(cell_df) <= 0) next

      print(unique(cell_df[,CELL]))
      print(cell)
      print(head(cell_df))

      pdf(paste(data_name, cell, "allGene_boxplot.pdf", sep="_"), height=10, width=8)
      #ggexport(p, filename=figname)

      p1a <- ggboxplot(cell_df, x=DOSE, y="Log2Intensity",
                       color = "black", palette = "jco", title=data_name, fill=DOSE,
                       xlab=FALSE, add="jitter", facet.by=c(TREATMENT, PHENOTYPE)) +
         geom_hline(yintercept = mean(cell_df$Log2Intensity), linetype = 2,
                    color="red3") + # Add horizontal line at base mean
         stat_compare_means(label = "p.signif", ref.group=categRef[DOSE],
                            method="t.test")   # Pairwise comparison against dose 0
      print(p1a)


      p1b <- ggboxplot(cell_df, x=TREATMENT, y="Log2Intensity", fill=TREATMENT,
                       color = TREATMENT, palette = "jco", title=data_name,
                       xlab=FALSE, add="jitter", facet.by=c(PHENOTYPE, DOSE)) +
         geom_hline(yintercept = mean(cell_df$Log2Intensity), linetype = 2, color="red3") + # Add horizontal line at base mean
         stat_compare_means(method = "anova", label.y = max(cell_df$Log2Intensity)+1)+        # Add global annova p-value
         stat_compare_means(label = "p.signif", ref.group=categRef[TREATMENT], method="t.test") # Pairwise comparison against NT
      print(p1b)

      p1c <- ggboxplot(cell_df, x=PHENOTYPE, y="Log2Intensity", fill=PHENOTYPE,
                       color = PHENOTYPE, palette = "jco", title=data_name,
                       xlab=FALSE, add="jitter", facet.by=c(TREATMENT, DOSE)) +
         geom_hline(yintercept = mean(cell_df$Log2Intensity), linetype = 2, color="red3") + # Add horizontal line at base mean
         stat_compare_means(method = "anova", label.y = max(cell_df$Log2Intensity)+1)+        # Add global annova p-value
         stat_compare_means(label = "p.signif", ref.group=categRef[PHENOTYPE], method="t.test") # Pairwise comparison against MOCK
      print(p1c)

      dev.off()


      if(length(unique(cell_df$Gene)) < 50){
         pdf(paste(data_name, cell, "perGene_boxplot.pdf", sep="_"), height=10, width=8)

         for(mod in unique(cell_df$Protein)){
            #mod <- "GPR89B|203"
            sub_df <- cell_df[cell_df$Protein == mod,]
            p5a <- ggboxplot(sub_df, x=DOSE, y="Log2Intensity",
                             color = DOSE, palette = "ucscgb", title=data_name, subtitle=mod,
                             xlab=FALSE, add="jitter", facet.by=c(TREATMENT, PHENOTYPE)) +
               stat_compare_means(label = "p.signif", ref.group=categRef[DOSE], method="t.test")
            print(p5a)

            p5b <- ggboxplot(sub_df, x=TREATMENT, y="Log2Intensity",
                             color = TREATMENT, palette = "ucscgb", title=data_name, subtitle=mod,
                             xlab=FALSE, add="jitter", facet.by=c(PHENOTYPE, DOSE)) +
               stat_compare_means(label = "p.signif", ref.group=categRef[TREATMENT], method="t.test")
            print(p5b)
            p5c <- ggboxplot(sub_df, x=PHENOTYPE, y="Log2Intensity",
                             color = PHENOTYPE, palette = "ucscgb", title=data_name, subtitle=mod,
                             xlab=FALSE, add="jitter", facet.by=c(TREATMENT, DOSE)) +
               stat_compare_means(label = "p.signif", ref.group=categRef[PHENOTYPE], method="t.test")
            print(p5c)
         }
         dev.off()
      }
   }
}

#' @title plot one factor boxplot
#' @description This function extracts protein abundance values from output of 'run_msdap' for selected genes, plot
#' boxplots for the group vaiable, with statistical significance annotations. The statistical analysis
#' is one-way anova followed by t-test on log transformed intensity values.
#'
#' @param wd working directory where sample information file 'sample.xlsx', and abundance data
#' 'protein_abundance__input data as-is.tsv' are located
#' @param geneNames a vector of characters denoting the selected gene names whose data will be analyzed
#' @param impute a string in ("halfmin", "wrproteo") denoting the imputation methods for missing values, default 'halfmin'
#' @param groupRef reference group for the group factor
#'
#' @return NULL
#' @author Shuye Pu

plot_for_msdap_group <- function(wd, geneNames, impute="halfmin", groupRef="Control"){
   if(0){
      impute="halfmin"
      treatRef="mock"
      drugRef="nd"
   }

   library(ComplexHeatmap)
   library(openxlsx)
   library(ggpubr)
   library(wrProteo)

   setwd(wd)
   xlsxfile <- file.path(wd, "samples.xlsx")
   sheets <- openxlsx::getSheetNames(xlsxfile)
   annot <- lapply(sheets, openxlsx::read.xlsx, xlsxFile=xlsxfile)[[1]]


   tsvfile <- file.path(wd, "protein_abundance__input data as-is.tsv")

   df_all <- read.delim2(tsvfile) %>%
      filter(!duplicated(gene_symbols_or_id)) # filter out protein groups with duplicated names
   rownames(df_all) <- df_all$gene_symbols_or_id

   if(length(intersect(df_all$gene_symbols_or_id, toupper(geneNames))) == 0){
      message("No target protein is found in the protein aboundance data!")
      return()
   }

   if(impute=="wrproteo"){
      dfImp <- matrixNAneighbourImpute(type.convert(as.matrix(df_all[, 4:ncol(df_all)]), as.is=T), annot$group)
      mat_all <- dfImp$data
      rownames(mat_all) <- df_all$gene_symbols_or_id
      colnames(mat_all) <- annot[, "shortname"]
   }else{
      mat_all <- type.convert(as.matrix(df_all[, 4:ncol(df_all)]), as.is=T)
      rownames(mat_all) <- df_all$gene_symbols_or_id
      colnames(mat_all) <- annot[, "shortname"]

      NAs <- length(mat_all[is.na(mat_all)])
      MIN <- min(mat_all[!is.na(mat_all)])/2
      SD <- sqrt(sd(mat_all[!is.na(mat_all)]))
      mat_all[is.na(mat_all)] <- rnorm(n=NAs, mean=MIN, sd=SD)

      print(paste("half of min(nonNA values", MIN, "square root of sd(nonNA values)", SD))
   }

   mat <- t(mat_all[toupper(unique(geneNames)),])

   pdf("protein_boxplot.pdf", width=10, height=8)
   p <- ComplexHeatmap::Heatmap(mat_all, name="log(intensity)", cluster_rows=T, cluster_columns=T,
                                column_names_rot = 45,
                                top_annotation = c(columnAnnotation(treat=annot[,"group"])),
                                column_names_side = "bottom",
                                show_row_names = F
   )

   draw(p)
   p <- ComplexHeatmap::Heatmap(mat, name="log(intensity)", cluster_rows=T, cluster_columns=T,
                                column_names_rot = 45,
                                right_annotation = c(rowAnnotation(treat=annot[,"group"])),
                                cell_fun = function(j, i, x, y, width, height, fill) {
                                   grid.text(sprintf("%.1f", mat[i, j]), x, y, gp = gpar(fontsize = 8))}
   )

   draw(p)
   #dev.off()

   df_wide <- as.data.frame(cbind(SampleID=annot$shortname, Group=annot$group, mat))
   df_long <- pivot_longer(df_wide, cols=colnames(mat), names_to = "Protein", values_to = "Log2Intensity") %>%
      mutate(Log2Intensity=as.numeric(Log2Intensity))
   comp <- combn(seq_along(unique(annot$group)),2, simplify=F)
   #pdf("Viral_protein_boxplot.pdf", height=8, width=6)
   p1 <- ggboxplot(df_long, x="Group", y="Log2Intensity",
                   color = "black", palette = "jco", title="selected proteins", fill="Group",
                   xlab=FALSE, add="jitter", add.params = list(size=1, color="darkred")) +
      geom_hline(yintercept = mean(df_long$Log2Intensity), linetype = 2, color="black") + # Add horizontal line at base mean
      stat_compare_means(label = "p.signif", ref.group=groupRef, method="t.test")   # Pairwise comparison against dose 0
   #print(p1)

   ps1 <- ggplot(df_long, aes(x=Group, y=Log2Intensity, fill=Group)) + scale_fill_manual(values=rainbow(length(unique(annot$group)))) +
      geom_boxplot(notch=FALSE) +
      geom_jitter() +
      theme_classic() +
      theme(axis.text.x = element_text(face="bold", size=10, color="black")) +
      theme(axis.title.x = element_blank(), axis.title.y = element_text(face="bold", size=10)) +
      theme(legend.position = "none") +
      labs(y=expression(paste(Log[2], "(Intensity)"))) +
      ggtitle("Selected proteins") +
      #stat_summary(fun=mean, geom="point", shape=23, size=4, fill="black") +
      geom_signif(comparisons = comp, test="t.test", map_signif_level=TRUE,  step_increase = 0.1)
   print(ps1)

   proteins <- unique(df_long$Protein)
   print(proteins)
   for(protein in proteins){
      protein_df <- df_long %>%
         filter(Protein == protein)

      p1 <- ggboxplot(protein_df, x="Group", y="Log2Intensity",
                      color = "black", palette = "jco", title=protein, fill="Group",
                      xlab=FALSE, add="jitter", add.params = list(size=1, color="darkred")) +
         geom_hline(yintercept = mean(df_long$Log2Intensity), linetype = 2, color="black") + # Add horizontal line at base mean
         stat_compare_means(label = "p.signif", ref.group=groupRef, method="t.test")   # Pairwise comparison against dose 0
      #print(p1)
      ps1 <- ggplot(protein_df, aes(x=Group, y=Log2Intensity, fill=Group)) + scale_fill_manual(values=rainbow(length(unique(annot$group)))) +
         geom_boxplot(notch=FALSE) +
         geom_jitter() +
         theme_classic() +
         theme(axis.text.x = element_text(face="bold", size=10, color="black")) +
         theme(axis.title.x = element_blank(), axis.title.y = element_text(face="bold", size=10)) +
         theme(legend.position = "none") +
         labs(y=expression(paste(Log[2], "(Intensity)"))) +
         ggtitle(protein) +
         #stat_summary(fun=mean, geom="point", shape=23, size=4, fill="black") +
         geom_signif(comparisons = comp, test="t.test", map_signif_level=TRUE,  step_increase = 0.1)
      print(ps1)
   }

   dev.off()
}

# Specify the comparisons for msdap, the format of the column|reference, e.g.
# comparisons <- "Treatment|UI,Time|0,Phenotype|MOCK" if allpairs = FALSE
# comparisons <- "group|UI+0+MOCK" is allpairs = TRUE
make_contrasts <- function(wd, comparisons, allpairs = FALSE){
   stopifnot(any(!is.null(comparisons), allpairs))
   library(readxl)

   setwd(wd)
   xlfile <- file.path(wd, "samples.xlsx")
   groupDesign <- as.data.frame(readxl::read_excel(xlfile, 1)) %>%
      dplyr::filter(!exclude)


   NAME <- colnames(groupDesign)[1]
   EXPERIMENT <- colnames(groupDesign)[2]
   TREATMENT <- colnames(groupDesign)[3]
   DOSE <- colnames(groupDesign)[4]
   PHENOTYPE <- colnames(groupDesign)[5]
   cells <- unique(groupDesign[,6])

   conditions <- groupDesign %>%
      dplyr::pull(group)
   names(conditions) <- groupDesign[, EXPERIMENT]

   categs <- unlist(strsplit(comparisons, split=","))

   categ <- unlist(lapply(categs, function(x)
      unlist(strsplit(x, split="|", fixed=T))[1]))
   categRef <- unlist(lapply(categs, function(x)
      unlist(strsplit(x, split="|", fixed=T))[2]))
   names(categRef) <- categ

   all_pairs <- combn(unique(conditions), 2, simplify = F)

   if(allpairs) {
      if(is.na(categRef["group"])) stop("Must be like 'comps = group|pbs_wt'")
      for(i in seq_along(all_pairs)){
         apair <- all_pairs[[i]]
         if(apair[2] == categRef["group"])
            all_pairs[[i]] <- rev(apair)
      }
      return(all_pairs)
   }

   contrasts <- lapply(all_pairs, function(factors){
      treatments <- unlist(lapply(factors, function(x)unlist(strsplit(x, split="+",fixed=T))[1]))
      doses <- unlist(lapply(factors, function(x)unlist(strsplit(x, split="+",fixed=T))[2]))
      phenotypes <- unlist(lapply(factors, function(x)unlist(strsplit(x, split="+",fixed=T))[3]))

      contrastN <- NULL
      ## reorder group contrast such that background group is fixed as the first element of the pair
      if(treatments[1] == treatments[2] && doses[1] == doses[2] && categRef[PHENOTYPE] %in% phenotypes){
         contrastN <- ifelse(phenotypes[1] == categRef[PHENOTYPE], list(factors), list(rev(factors)))
      }

      if(phenotypes[1] == phenotypes[2] && doses[1] == doses[2] && categRef[TREATMENT] %in% treatments){
         contrastN <- ifelse(treatments[1] == categRef[TREATMENT], list(factors), list(rev(factors)))
      }

      if(phenotypes[1] == phenotypes[2] && treatments[1] == treatments[2] && categRef[DOSE] %in% doses){
         contrastN <- ifelse(doses[1] == categRef[DOSE], list(factors), list(rev(factors)))
      }
      names(contrastN) <- NULL
      return(unlist(contrastN))
   })
   return(contrasts[!sapply(contrasts, is.null)])

}

make_contrasts_old <- function(wd, treatRef="mock", groupRef="nd"){

      library(openxlsx)

      setwd(wd)
      xlsxfile <- file.path(wd, "samples.xlsx")
      sheets <- openxlsx::getSheetNames(xlsxfile)
      annot <- lapply(sheets, openxlsx::read.xlsx, xlsxFile=xlsxfile)[[1]]

      grouping <- annot$group

      all_pairs <- combn(unique(grouping), 2, simplify = F)
      contrasts <- lapply(all_pairs, function(x){
         g1 <- unlist(strsplit(x[1], split="_"))
         g2 <- unlist(strsplit(x[2], split="_"))

         if(g1[1]==g2[1] && g1[2]!=g2[2]){
            if(g1[2]==groupRef){
               return(x)
            }else if(g2[2]==groupRef){
               return(rev(x))
            }
         }
         if(g1[2]==g2[2] && g1[1]!=g2[1]){
            if(g1[1]==treatRef){
               return(x)
            }else if(g2[1]==treatRef){
               return(rev(x))
            }
         }

      })
      return(contrasts[!sapply(contrasts, is.null)])
}

# plot peptide intensity for each gene in geneNames by group
plot_peptide_data_group <- function(lwd, mydataset, geneNames){

   sample_data <- mydataset$samples
   group_data <- mydataset$groups
   protein_data <- mydataset$proteins %>%
      filter(gene_symbols_or_id %in% geneNames)
   if(nrow(protein_data) == 0){
      message("No target protein is found in the protein aboundance data!")
      return()
   }

   peptide_data <- mydataset$peptides %>%
      filter(isdecoy==FALSE) %>%
      filter(detect==TRUE) %>%
      filter(protein_id %in% protein_data$protein_id) %>%
      merge(group_data) %>%
      merge(protein_data) %>%
      select(peptide_id, protein_id, sample_id, confidence, intensity_all_group, group, gene_symbols_or_id)

   pdf(file.path(lwd, "peptide_boxplot.pdf"), height=6, width=8)
   comp <- combn(group_data$key_group ,2, simplify=F)
   for(protein in unique(geneNames)){
      protein_df <- peptide_data %>%
         filter(gene_symbols_or_id == protein)

      ps1 <- ggplot(protein_df, aes(x=group, y=intensity_all_group, color=group)) + scale_fill_manual(values=rainbow(nrow(group_data))) +
         geom_boxplot(notch=FALSE, outlier.shape = NA, size=2) +
         geom_jitter(size=2) +
         theme_classic() +
         theme(axis.text.x = element_text(face="bold", size=10, color="black")) +
         theme(axis.title.x = element_blank(), axis.title.y = element_text(face="bold", size=10)) +
         theme(legend.position = "none") +
         labs(y=expression(paste(Log[2], "(Peptide intensity)"))) +
         ggtitle(protein) +
         #stat_summary(fun=mean, geom="point", shape=23, size=4, fill="black") +
         geom_signif(comparisons = comp, test="t.test", map_signif_level=TRUE,  step_increase = 0.1)
      print(ps1)
   }

   dev.off()

}

