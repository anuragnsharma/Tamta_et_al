# Load packages ----
library(tidyverse)
library(ggpubr)
library(ggsci)
library(ggthemes)
library(extrafont)
library(extrafontdb)
library(ggtext)
library(edgeR)
library(conflicted)
library(dendextend)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
# Functions involving one or more repeats of borders around variable length of strings ----
F02_get_between_strings <- function(colname,left_border, right_border){
  # This does it only for one pattern match
  str_match(colname,str_c(left_border,"\\s*(.*?)\\s*",right_border)) %>% 
    as_tibble(.name_repair = "unique") %>% 
    dplyr::pull(2)
}

F02_get_between_strings_all <- function(colname,left_border, right_border){
  # This does it only for all pattern match
  str_match_all(colname,str_c(left_border,"\\s*(.*?)\\s*",right_border)) %>% 
    purrr::pluck(1) %>% 
    as_tibble(.name_repair = "unique") %>% 
    dplyr::pull(2)
}

F02_get_including_and_between_strings <- function(colname,left_border, right_border){
  # This does it only for one pattern match
  str_match(colname,str_c(left_border,"\\s*(.*?)\\s*",right_border)) %>% 
    as_tibble(.name_repair = "unique") %>% 
    dplyr::pull(1)
}

F02_get_including_and_between_strings_all <- function(colname,left_border, right_border){
  # This does it only for all pattern match
  str_match_all(colname,str_c(left_border,"\\s*(.*?)\\s*",right_border)) %>% 
    purrr::pluck(1) %>% 
    as_tibble(.name_repair = "unique") %>% 
    dplyr::pull(1)
}

F02_get_left_border_exclude_right_border <- function(colname,left_border, right_border){
  # This does it only for one pattern match
  str_match(colname,str_c(left_border,"\\s*(.*?)\\s*",right_border)) %>% 
    as_tibble(.name_repair = "unique") %>% 
    dplyr::pull(1) %>% 
    str_replace(right_border,"") %>% 
    str_squish()
}

F02_get_left_border_exclude_right_border_all <- function(colname,left_border, right_border){
  # This does it only for many pattern matches in the same string
  str_match_all(colname,str_c(left_border,"\\s*(.*?)\\s*",right_border)) %>% 
    purrr::pluck(1) %>% 
    as_tibble(.name_repair = "unique") %>% 
    dplyr::pull(1) %>% 
    str_replace(right_border,"") %>% 
    str_squish()
}

F02_get_right_border_exclude_left_border <- function(colname,left_border, right_border){
  # This does it only for one pattern match
  str_match(colname,str_c(left_border,"\\s*(.*?)\\s*",right_border)) %>% 
    as_tibble(.name_repair = "unique") %>% 
    dplyr::pull(1) %>% 
    str_replace(left_border,"") %>% 
    str_squish()
}

F02_get_right_border_exclude_left_border_all <- function(colname,left_border, right_border){
  # This does it only for many pattern matches in the same string
  str_match_all(colname,str_c(left_border,"\\s*(.*?)\\s*",right_border)) %>% 
    purrr::pluck(1) %>% 
    as_tibble(.name_repair = "unique") %>% 
    dplyr::pull(1) %>% 
    str_replace(left_border,"") %>% 
    str_squish()
}

# Set working directory ----
setwd("/Ankit/Ankit_RNA_seq_June_2024")
# Save locations of each file ----
Location_of_fastq.gz_files <- getwd()
Location_genome_and_gtf_file <- "/General/ncbi_dataset/Mouse/GRCm39/GCF_000001635.27"
Zero_Location_of_index <- "/General/ncbi_dataset/Mouse/GRCm39/GCF_000001635.27/Index"
Zero_location_of_code <- getwd()

# Ensure GTF file has only exons, and not other features (such as UTR or promoter)
# Filter GTF file to contain only exons and write
# rtracklayer::readGFFAsGRanges(str_c(Location_genome_and_gtf_file,"enter name of .gtf file containing only the exons")) %>%
#   plyranges::filter(type == "exon") %>%
#   rtracklayer::export(.,str_c(Location_genome_and_gtf_file,"enter name of .gtf file containing only the exons"))

# Step 1: fastp adapter trimming ----
Step_01_fastp_adapter_trimming <-
  list.files(Location_of_fastq.gz_files,full.names = T,recursive = T) %>%
  # Pick only the fastq files
  str_subset(".fastq.gz$") %>% 
  str_subset("fastp|sha256sum|html|json|\\.bam|\\.sam",negate = T) %>% 
  # Put into a tibble
  tibble::enframe(name = NULL, value = "Directory") %>% 
  # Create an ID column 
  dplyr::mutate(`Read 1 or Read 2` = if_else(str_detect(Directory, "\\_R1"),"Read 1","Read 2"),
                ID = str_replace(Directory,"_R[12]","")) %>% 
  # Arrange the tibble to prepare the code, pivot the Read column
  tidyr::pivot_wider(id_cols = ID, names_from = `Read 1 or Read 2`, values_from = Directory) %>% 
  dplyr::mutate(
    Code = str_c(
      "fastp",
      "--in1", `Read 1`,
      "--in2", `Read 2`,
      "--out1", str_replace(`Read 1`, "\\.fastq","_fastp_trimmed.fastq"),
      "--out2", str_replace(`Read 2`, "\\.fastq","_fastp_trimmed.fastq"),
      "--json", str_replace(ID, "\\.fastq\\.gz",".json"),
      "--html", str_replace(ID, "\\.fastq\\.gz",".html"),
      "--detect_adapter_for_pe",
      # "--qualified_quality_phred 30", this is what I've used before
      "--qualified_quality_phred 20", # reducing this threshold because there were only 2 DEGs with above condition
      # "--unqualified_percent_limit 50", this is what I've used before
      # "--average_qual 30", this is what I've used before
      "--correction",
      "--thread 10",
      sep = " "
    )
  ) %>% 
  dplyr::mutate(`Order of code` =  1,
                `Order within subset` = 1:10)

# Step 2: star build index ----
# 
# Step_02_star_build_index <- str_c(
#   "star ",
#   "--runThreadN 10",
#   "--runMode genomeGenerate",
#   "--genomeDir",
#   Zero_Location_of_index,
#   "--genomeFastaFiles",
#   str_c(Location_genome_and_gtf_file,"/enter name of fasta/.fna/.fa genome file here"),
#   "--sjdbGTFfeatureExon",
#   str_c(Location_genome_and_gtf_file,"/enter name of .gtf file containing only the exons"),
#   sep = " "
# )

# Step 3: star alignment ----

Step_03_star_alignment <-
  list.files(Location_of_fastq.gz_files,full.names = T,recursive = T) %>% 
  # Pick only the fastq files
  str_subset(".fastq.gz$") %>% 
  str_subset("fastp|sha256sum|html|json|\\.bam|\\.sam",negate = T) %>% 
  # Create input names 
  str_replace("\\.fastq\\.gz","_fastp_trimmed.fastq.gz") %>% 
  # Put into a tibble
  tibble::enframe(name = NULL, value = "Directory") %>% 
  # Create an ID column 
  dplyr::mutate(`Read 1 or Read 2` = if_else(str_detect(Directory, "\\_R1"),"Read 1","Read 2"),
                Input = str_replace(Directory, ".fastq","_fastp_trimmed.fastq"),
                ID = str_replace(Directory,"\\_R[12]","")) %>% 
  # Arrange the tibble to prepare the code, pivot the Read column
  tidyr::pivot_wider(id_cols = ID, names_from = `Read 1 or Read 2`, values_from = Directory) %>% 
  dplyr::mutate(Output = 
                  str_replace(`Read 1`,"_fastp_trimmed.fastq","_fastp_trimmed_and_star_") %>% 
                  str_replace_all("\\_R1|\\.gz","")) %>% 
  dplyr::mutate(
    Code = str_c(
      "star",
      "--runThreadN 10",
      "--genomeDir",
      Zero_Location_of_index,
      "--readFilesIn",
      `Read 1`,
      `Read 2`,
      "--readFilesCommand gunzip -c",
      "--outFileNamePrefix",
      Output,
      sep = " "
    )
  ) %>% 
  dplyr::mutate(`Order of code` =  2,
                `Order within subset` = 1:10)


# Step 4: samtools operations 01 - samtools sam to bam ----

Step_04_samtools_01_convert_sam_to_bam <-
  list.files(Location_of_fastq.gz_files,full.names = T,recursive = T) %>% 
  # Pick only the fastq files
  str_subset(".fastq.gz") %>% 
  str_subset("fastp|sha256sum|html|json|\\.bam|\\.sam",negate = T) %>% 
  # Put into a tibble
  tibble::enframe(name = NULL, value = "Directory") %>% 
  # Create an ID column 
  dplyr::mutate(ID = str_replace(Directory,"\\_R[12]","")) %>% 
  dplyr::distinct(ID, .keep_all = T) %>% 
  dplyr::select(Directory) %>% 
  dplyr::mutate(Input = 
                  str_replace(Directory, "\\_R[12]","") %>% 
                  str_replace("\\.fastq\\.gz","_fastp_trimmed_and_star_Aligned.out.sam"),
                Output = str_replace(Input, "\\.sam",".bam")) %>% 
  # Prepare the code
  dplyr::mutate(
    Code = str_c(
      "samtools view -bS",
      Input,
      ">",
      Output,
      "--threads 10",
      sep = " "
    )
  ) %>% 
  dplyr::mutate(`Order of code` =  3,
                `Order within subset` = 1:10)

# Step 5: remove sam files, because they are taking up too much space ----
Step_05_remove_sam_files <-
  list.files(Location_of_fastq.gz_files,full.names = T,recursive = T) %>% 
  # Pick only the fastq files
  str_subset(".fastq.gz") %>% 
  str_subset("fastp|sha256sum|html|json|\\.bam|\\.sam",negate = T) %>% 
  # Put into a tibble
  tibble::enframe(name = NULL, value = "Directory") %>% 
  # Create an ID column 
  dplyr::mutate(ID = str_replace(Directory,"\\_R[12]","")) %>% 
  dplyr::distinct(ID, .keep_all = T) %>% 
  dplyr::select(Directory) %>% 
  dplyr::mutate(Code = 
                  str_replace(Directory, "\\_R[12]","") %>% 
                  str_replace("\\.fastq\\.gz","_fastp_trimmed_and_star_Aligned.out.sam") %>% 
                  str_c("rm ",.)) %>% 
  dplyr::mutate(`Order of code` =  4,
                `Order within subset` = 1:10)
# Step 6: samtools operations 02 - samtools keep paired end reads only, discard single end reads ----

Step_06_samtools_02_convert_bam_to_paired_bam <-
  list.files(Location_of_fastq.gz_files,full.names = T,recursive = T) %>% 
  # Pick only the fastq files
  str_subset(".fastq.gz") %>% 
  str_subset("fastp|sha256sum|html|json|\\.bam|\\.sam",negate = T) %>% 
  # Put into a tibble
  tibble::enframe(name = NULL, value = "Directory") %>% 
  # Create an ID column 
  dplyr::mutate(ID = str_replace(Directory,"\\_R[12]","")) %>% 
  dplyr::distinct(ID, .keep_all = T) %>% 
  dplyr::select(Directory) %>% 
  dplyr::mutate(Input = 
                  str_replace(Directory, "\\_R[12]","") %>% 
                  str_replace("\\.fastq\\.gz","_fastp_trimmed_and_star_Aligned.out.bam"),
                Output = str_replace(Input, "\\.bam","_paired_end_only.bam")) %>% 
  # Prepare the code
  dplyr::mutate(
    Code = str_c(
      "samtools view -bf 1",
      Input,
      ">",
      Output,
      "--threads 10",
      sep = " "
    )
  ) %>% 
  dplyr::mutate(`Order of code` =  5,
                `Order within subset` = 1:10)

# Step 7: featureCounts ----
Step_07_featureCounts <-
  list.files(Location_of_fastq.gz_files,full.names = T,recursive = T) %>% 
  # Pick only the fastq files
  str_subset("\\.fastq.gz") %>% 
  str_subset("fastp|sha256sum|html|json|\\.bam|\\.sam",negate = T) %>% 
  # Put into a tibble
  tibble::enframe(name = NULL, value = "Directory") %>% 
  # Create an ID column 
  dplyr::mutate(ID = str_replace(Directory,"\\_R[12]","")) %>% 
  dplyr::distinct(ID, .keep_all = T) %>% 
  dplyr::select(Directory) %>% 
  dplyr::mutate(Input = 
                  str_replace(Directory, "\\_R[12]","") %>% 
                  str_replace("\\.fastq\\.gz","_fastp_trimmed_and_star_Aligned.out_paired_end_only.bam"),
                Output = 
                  str_replace(Input, "\\.bam","_counts.txt") %>% 
                  # Save in the previous cd ../ location
                  str_replace("augmet-download\\/","")) %>% 
  # Prepare the code
  dplyr::mutate(
    Code = str_c(
      "featureCounts -p -O",
      "-T 10",
      "-a", str_c(Location_genome_and_gtf_file,"/exons.gtf"),
      "-o", Output,
      Input,
      sep = " "
    )
  ) %>% 
  dplyr::mutate(`Order of code` =  6,
                `Order within subset` = 1:10)

# Step final: Put all the code together ----
list(
  Step_01_fastp_adapter_trimming, # Step 1
  # Step_02_star_build_index, # Step 2
  Step_03_star_alignment, # Step 3
  Step_04_samtools_01_convert_sam_to_bam, # Step 4
  Step_05_remove_sam_files, # Step 6
  Step_06_samtools_02_convert_bam_to_paired_bam, # Step 5
  Step_07_featureCounts # Step 7
) %>% 
  purrr::map(~{
    .x %>% 
      dplyr::select(Code, `Order of code`, `Order within subset`)
  }) %>% 
  purrr::reduce(bind_rows) %>% 
  dplyr::arrange(`Order within subset`, `Order of code`) %>% 
  dplyr::group_by(`Order within subset`) %>% 
  dplyr::group_split() %>% 
  purrr::map(~{
    .x %>% 
      dplyr::pull(Code) %>% 
      str_c(collapse = "\n\n")
  }) %>%
  unlist %>% 
  str_c(collapse = "\n\n#\t\t------------------------------------------\t\t# \n\n") %>% 
  write_lines(str_c(Zero_location_of_code,"/00_Complete_script.sh"))

# Concatenate all code together with a one empty line space

# Step 8: Remove all extra files and keep only the final bam file ----
list.files(getwd(), full.names = T) %>% 
  str_subset("paired_end_only|html|json|R[12]\\.fastq\\.gz|fastp_trimmed\\.fastq\\.gz|\\.sh|PCA|[qQ][cC]|Code", negate = T) %>% 
  str_c("rm ", .) %>% 
  str_c(collapse = "\n\n") %>% 
  write_lines("Remove_extras.sh")
# # Optional step 8 when I want to restart the analysis ----
#   list.files(Location_of_fastq.gz_files, full.names = T) %>%
#   str_subset("R[12]\\.fastq\\.gz|\\.sh|Code", negate = T) %>%
#   str_c("rm ", .) %>% 
#   str_c(collapse = "\n\n") %>% 
#   write_lines("Restart.sh")
# Read count data ----

Sample_names <- 
  list.files(Location_of_fastq.gz_files,
             full.names = F, recursive = T) %>% 
  # Pick only the fastq files
  str_subset(".fastq.gz") %>% 
  str_subset("fastp|sha256sum|html|json|\\.bam|\\.sam",negate = T) %>% 
  # Put into a tibble
  tibble::enframe(name = NULL, value = "Directory") %>% 
  # Create an ID column 
  dplyr::mutate(ID = str_replace(Directory,"\\_R[12]","")) %>% 
  dplyr::distinct(ID, .keep_all = T) %>% 
  dplyr::pull(Directory) %>% 
  str_replace("\\_R[12]","") %>% 
  str_replace("\\.fastq\\.gz","")

Counts <- seqUtils::loadCounts(
  sample_dir = Location_of_fastq.gz_files,
  sample_names = Sample_names,
  counts_suffix = "_fastp_trimmed_and_star_Aligned.out_paired_end_only_counts.txt",
  sub_dir = F
) %>% 
  as_tibble %>% 
  dplyr::rename(`Gene name` = gene_id)

Counts %>% 
  write_csv(str_c(Location_of_fastq.gz_files,"/Readcounts.csv"))

Names <-
  Counts %>%
  colnames %>% 
  tibble::enframe(name = NULL, value = "Colnames") %>% 
  dplyr::mutate(Names = 
                  str_extract(Colnames,"T2[WK][TO]-[0-9]"),
                Names = if_else(is.na(Names),Colnames, Names)) %>% 
  dplyr::pull(Names)

Counts <- Counts %>% 
  setNames(Names) %>% 
  dplyr::select(-length)

# ----
# Set working directory ----
setwd("/Ankit/Ankit_RNA_seq_June_2024")
# Reading in readcount data ----
Counts <- read_csv("Readcounts.csv")
# Other functions ---- 
F01_get_DEGs <- function(Input_data){
  # Combinations of WT and KO selected for this analysis
  # Preparing the groups
  Group <- 
    Input_data %>% 
    dplyr::select(matches("KO|WT")) %>% 
    colnames %>% 
    str_replace("[1234567]","") %>% 
    as_factor()
  
  # Preparing the DGEset object for analysis 
  DGEset <- DGEList(
    counts = as.matrix(Input_data %>% 
                         dplyr::select(matches("KO|WT"))),
    group = Group,
    genes = Input_data %>% select("Gene name")
  )
  colnames(DGEset) <- colnames(Input_data %>% 
                                 dplyr::select(matches("KO|WT")))
  
  # Filtering out low expressed genes
  keep_DGEset <- filterByExpr(DGEset)
  DGEset <- DGEset[keep_DGEset,,
                   keep.lib.sizes = F]
  
  # Calculating library size and normalization factor
  DGEset <- calcNormFactors(DGEset)
  
  # plotMDS(DGEset)
  
  # Obtaining CPM values that is to be used in the analysis later
  # CPM <- cpmByGroup(DGEset, log = F, normalized.lib.sizes = TRUE) %>%
  #   apply(.,1,function(x) (x - mean(x)) / sd(x)) %>%
  #   as_tibble(.name_repair = "universal") %>%
  #   data.table::transpose(.) %>%
  #   dplyr::rename_all(., ~levels(Group)) %>%
  #   bind_cols(DGEset %>% pluck("genes")) %>%
  #   dplyr::select(SystematicName,Col0,J,JKM,JKD)
  
  # Obtaining the differentially expressed genes
  Design <- model.matrix(~0 + Group)
  colnames(Design) <- levels(Group)
  DGEsetDisp <- estimateDisp(DGEset,
                             Design,
                             robust = T)
  Fit <- glmQLFit(DGEsetDisp,
                  Design,
                  robust = T)
  getDEGs <- function(Contrast){
    glmQLFTest(Fit,contrast = Contrast) %>%
      topTags(.,DGEset %>% pluck(3) %>% nrow()) %>%
      pluck(1) %>%
      dplyr::select(-"logCPM",-"F",-"PValue") %>%
      dplyr::mutate(logFC = round(logFC,2)) %>%
      as_tibble
  }
  
  KO_WT <- makeContrasts(KO - WT, levels = Design) %>%
    getDEGs %>%
    dplyr::filter(FDR < 0.05) %>%
    dplyr::left_join(Input_data) %>% 
    dplyr::rename(`Log fold change`  = logFC,
                  `Adjusted p-value` = FDR) %>% 
    dplyr::arrange(`Adjusted p-value`)
}
F01_get_DEGs_unfiltered <- function(Input_data){
  # Combinations of WT and KO selected for this analysis
  # Preparing the groups
  Group <- 
    Input_data %>% 
    dplyr::select(matches("KO|WT")) %>% 
    colnames %>% 
    str_replace("[1234567]","") %>% 
    as_factor()
  
  # Preparing the DGEset object for analysis 
  DGEset <- DGEList(
    counts = as.matrix(Input_data %>% 
                         dplyr::select(matches("KO|WT"))),
    group = Group,
    genes = Input_data %>% select("Gene name")
  )
  colnames(DGEset) <- colnames(Input_data %>% 
                                 dplyr::select(matches("KO|WT")))
  
  # Filtering out low expressed genes
  keep_DGEset <- filterByExpr(DGEset)
  DGEset <- DGEset[keep_DGEset,,
                   keep.lib.sizes = F]
  
  # Calculating library size and normalization factor
  DGEset <- calcNormFactors(DGEset)
  
  # plotMDS(DGEset)
  
  # Obtaining CPM values that is to be used in the analysis later
  # CPM <- cpmByGroup(DGEset, log = F, normalized.lib.sizes = TRUE) %>%
  #   apply(.,1,function(x) (x - mean(x)) / sd(x)) %>%
  #   as_tibble(.name_repair = "universal") %>%
  #   data.table::transpose(.) %>%
  #   dplyr::rename_all(., ~levels(Group)) %>%
  #   bind_cols(DGEset %>% pluck("genes")) %>%
  #   dplyr::select(SystematicName,Col0,J,JKM,JKD)
  
  # Obtaining the differentially expressed genes
  Design <- model.matrix(~0 + Group)
  colnames(Design) <- levels(Group)
  DGEsetDisp <- estimateDisp(DGEset,
                             Design,
                             robust = T)
  Fit <- glmQLFit(DGEsetDisp,
                  Design,
                  robust = T)
  getDEGs <- function(Contrast){
    glmQLFTest(Fit,contrast = Contrast) %>%
      topTags(.,DGEset %>% pluck(3) %>% nrow()) %>%
      pluck(1) %>%
      dplyr::select(-"logCPM",-"F",-"PValue") %>%
      dplyr::mutate(logFC = round(logFC,2)) %>%
      as_tibble
  }
  
  makeContrasts(KO - WT, levels = Design) %>%
    getDEGs() %>%
    # dplyr::filter(FDR < 0.05) %>%
    dplyr::left_join(Input_data) %>% 
    dplyr::rename(`Log fold change`  = logFC,
                  `Adjusted p-value` = FDR) %>% 
    dplyr::arrange(`Adjusted p-value`)
}
F01_plot_MDS <- function(Input_data){
  # Combinations of WT and KO selected for this analysis
  # Preparing the groups
  Group <- 
    Input_data %>% 
    dplyr::select(matches("KO|WT")) %>% 
    colnames %>% 
    str_replace("[1234567]","") %>% 
    as_factor()
  
  # Preparing the DGEset object for analysis 
  DGEset <- DGEList(
    counts = as.matrix(Input_data %>% 
                         dplyr::select(matches("KO|WT"))),
    group = Group,
    genes = Input_data %>% select("Gene name")
  )
  colnames(DGEset) <- colnames(Input_data %>% 
                                 dplyr::select(matches("KO|WT")))
  
  # Filtering out low expressed genes
  keep_DGEset <- filterByExpr(DGEset)
  DGEset <- DGEset[keep_DGEset,,
                   keep.lib.sizes = F]
  
  # Calculating library size and normalization factor
  DGEset <- calcNormFactors(DGEset)
  
  MDS <- plotMDS(DGEset)
  
  X_explained <- (((MDS %>% purrr::pluck("var.explained"))[1])*100) %>% round(2) %>% 
    str_c("Leading logFC dimension 1 (",.,"%)")
  Y_explained <- (((MDS %>% purrr::pluck("var.explained"))[2])*100) %>% round(2) %>% 
    str_c("Leading logFC dimension 2 (",.,"%)")
  
  Data <- tibble(
    X = MDS %>%
      purrr::pluck("x"),
    Y = MDS %>%
      purrr::pluck("y"),
    Text = colnames(DGEset),
    `Text category` = Group
  )
  Data %>%
    ggplot(aes(x = X, y = Y)) +
    geom_point(size = 2) +
    geom_text_repel(aes(label = Text, color = `Text category`), family = "Times New Roman", size = 10) +
    labs(x = X_explained, y = Y_explained) +
    theme_classic() +
    theme(
      panel.background = element_rect(fill = "white",colour = "gray10"),
      axis.title = element_text(family = "Times New Roman",size = 20),
      axis.text = element_text(family = "Times New Roman",size = 20),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    ) +
    scale_color_manual(values = c("brown","black")) +
    ggsave("PCA plot.png", device = "png",
           height = 10,
           width = 10)
}
F01_get_CPM <- function(Input_data){
  # Combinations of WT and KO selected for this analysis
  # Preparing the groups
  Group <- 
    Input_data %>% 
    dplyr::select(matches("KO|WT")) %>% 
    colnames %>% 
    str_replace("[1234567]","") %>% 
    as_factor()
  
  # Preparing the DGEset object for analysis 
  DGEset <- DGEList(
    counts = as.matrix(Input_data %>% 
                         dplyr::select(matches("KO|WT"))),
    group = Group,
    genes = Input_data %>% select("Gene name")
  )
  colnames(DGEset) <- colnames(Input_data %>% 
                                 dplyr::select(matches("KO|WT")))
  
  # Filtering out low expressed genes
  keep_DGEset <- filterByExpr(DGEset)
  DGEset <- DGEset[keep_DGEset,,
                   keep.lib.sizes = F]
  
  # Calculating library size and normalization factor
  DGEset <- calcNormFactors(DGEset)
  
  # plotMDS(DGEset)
  
  # Obtaining CPM values that is to be used in the analysis later
  CPM <- cpmByGroup(DGEset, log = F, normalized.lib.sizes = TRUE) %>% 
    as_tibble() %>% 
    dplyr::rename_all(~levels(Group)) %>%
    bind_cols(DGEset %>% pluck("genes")) 
  # tidyr::pivot_longer(cols = !matches("Gene name"),names_to = "Genotype",values_to = "CPM") %>% 
  # dplyr::group_by(`Gene name`) %>% 
  # dplyr::mutate(CPM = (CPM - mean(CPM)) / sd(CPM)) %>% 
  # tidyr::pivot_wider(id_cols = `Gene name`,names_from = Genotype, values_from = `CPM`)
}
# MDS plot ----
Counts %>% 
  setNames(colnames(.) %>% str_replace("T2KO-","KO") %>% str_replace("T2WT-","WT") %>% str_replace("_L3","")) %>% 
  F01_plot_MDS()

# ----
# Get differentially expressed genes ----
DEGs <-
  Counts %>% 
  dplyr::select(-length) %>% 
  setNames(colnames(.) %>% str_extract("Gene name|KO-\\d|WT-\\d") %>% str_replace("-","")) %>% 
  dplyr::select(matches("G|KO[134]|WT[567]")) %>% 
  F01_get_DEGs()

DEGs %>% 
  dplyr::arrange(`Adjusted p-value`, `Log fold change`) %>% 
  write_csv("KO1-KO3-KO4 vs WT5-WT6-WT7 (3783 DEGs - 1607 down and 2176 up).csv")

DEGs %>% 
  dplyr::arrange(`Gene name`) %>% 
  write_csv("KO1-KO3-KO4 vs WT5-WT6-WT7 (3783 DEGs - 1607 down and 2176 up) - alphabetically ordered.csv")

DEGs %>% 
  dplyr::mutate(Down = if_else(`Log fold change` < 0, "Down","Up")) %>% 
  group_by(Down) %>% 
  dplyr::group_split() %>% 
  purrr::map(~{
    .x %>% 
      dplyr::select(-Down) %>% 
      dplyr::arrange(`Gene name`) 
  }) %>% 
  openxlsx::write.xlsx("KO1-KO3-KO4 vs WT5-WT6-WT7 (3783 DEGs - 1607 down and 2176 up) - split by up or downregulated into two sheets.xlsx")


# ----
# Volcano plot ----
Counts %>% 
  dplyr::select(-length) %>% 
  setNames(colnames(.) %>% str_extract("Gene name|KO-\\d|WT-\\d") %>% str_replace("-","")) %>% 
  dplyr::select(matches("G|KO[134]|WT[567]")) %>% 
  F01_get_DEGs_unfiltered() %>% 
  dplyr::mutate(Less = if_else(`Adjusted p-value` < 0.05, "Yes", "No"),
                `-log10 (Adjusted p-value)` = -log10(`Adjusted p-value`)) %>% 
  ggplot(aes(y = `-log10 (Adjusted p-value)`, x = `Log fold change`)) +
  geom_point(aes(color = Less), size = 0.8) +
  geom_hline(yintercept = -log10(0.05), linewidth = 1.5, color = "grey25", lty = "dotted") +
  scale_color_manual(values = c("grey","red")) +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white"),
    axis.line = element_line(color = "grey"),
    axis.text = element_text(family = "Helvetica"),
    axis.title = element_text(family = "Helvetica"),
    legend.position = "none"
  ) +
  ggsave("Volcano plot.png",
         height = 5,
         width = 5)

# ----
# Generating CPM of DEGs (to be used for heatmaps) ----

CPM <- Counts %>%
  F01_get_CPM() %>%
  setNames(colnames(.) %>% str_replace("TKO","KO") %>% str_replace("TWT","WT") %>% str_replace("_L3","") %>% str_replace("-","")) %>%
  dplyr::select(KO1,KO3,KO4,WT5,WT6,WT7,`Gene name`)

DEGs_CPM <-
  CPM %>%
  dplyr::filter(`Gene name` %in% (DEGs %>% dplyr::pull(`Gene name`)))

# Generating z-score of CPM of DEGs (to be used for heatmaps) ----
DEGs_CPM_z_score <-
  DEGs_CPM %>%
  tidyr::pivot_longer(cols = matches("WT|KO"),names_to = "Genotype",values_to = "Original") %>%
  dplyr::group_by(`Gene name`) %>%
  dplyr::mutate(Normalized = (Original - mean(Original)) / sd(Original)) %>%
  dplyr::select(-Original) %>%
  tidyr:::pivot_wider(names_from = Genotype, values_from = Normalized) %>%
  dplyr::ungroup()

remove(DEGs_CPM)
# ----
# Latest code ----
# Arathi-list - extracting genes enriched in the chosen pathways - data loading ----
All_metascape <-
  list.files("Metascape/",recursive = T,full.names = T) %>% 
  str_subset("xlsx|csv") %>% 
  str_subset("Arathi",negate = T) %>% 
  purrr::map(~{
    if(str_detect(.x,"csv")){
      Data =  read_csv(.x)
    }else{
      Data =  readxl::read_excel(.x)
    }
    Data %>% 
      dplyr::mutate(
        `Up or downregulated GO set` = str_extract(.x,"Upregulated|Downregulated")
      ) %>% 
      dplyr::relocate(`Up or downregulated GO set`)
  })
# Arathi-list - preparing data for plots - heatmap ----
All_metascape_heatmap <-
  All_metascape %>% 
  purrr::reduce(bind_rows) %>% 
  list() %>% 
  purrr::map2(.x = .,
              .y = list("Upregulated","Downregulated"),
              ~{
                .x %>% 
                  dplyr::filter(`Up or downregulated GO set` == .y)
              }) %>% 
  purrr::map2(.x = .,
              .y = list("MAPK6\\/MAPK4 signaling|^apoptotic signaling pathway|^regulation of autophagy|regulation of muscle atrophy|skeletal muscle contraction",
                        "AMPK signaling pathway - Mus musculus \\(house mouse\\)|Negative regulation of MAPK pathway|mTOR signaling pathway - Mus musculus \\(house mouse\\)|regulation of ATP metabolic process"),
              ~{
                List = .x %>% 
                  dplyr::filter(str_detect(Description,.y)) %>% 
                  dplyr::filter(!is.na(Category)) %>% 
                  dplyr::filter(!is.na(`_MEMBER_MyList`)) %>% 
                  dplyr::select(`Up or downregulated GO set`,
                                Description,
                                LogP,
                                Category,
                                Hits) %>% 
                  list()
                names(List) = .x %>% 
                  dplyr::pull(`Up or downregulated GO set`) %>% 
                  unique
                List
              }) %>% 
  purrr::flatten()

All_metascape_heatmap %>% 
  openxlsx::write.xlsx("Metascape/Heatmap related gene list.xlsx")

All_metascape_heatmap <- 
  All_metascape_heatmap %>% 
  purrr::reduce(bind_rows)


# Arathi-list - preparing data for plots - pathways ----
Filter <- list(
  readxl::read_excel("Metascape/Arathi - Enriched GO related genes.xlsx",sheet = "upregulated pathway"),
  readxl::read_excel("Metascape/Arathi - Enriched GO related genes.xlsx", sheet = "downregulated pathway")
) %>% purrr::map(~{.x %>% dplyr::pull(`GO description`) %>% str_c("^",.,"$") %>% str_c(collapse = "|") %>% str_replace_all("\\(","\\\\(") %>% str_replace_all("\\)","\\\\)") %>% str_replace_all("\\,","\\\\,")})

All_metascape_pathways <-
  All_metascape %>% 
  purrr::reduce(bind_rows) %>% 
  list() %>% 
  purrr::map2(.x = .,
              .y = list("Upregulated","Downregulated"),
              ~{
                .x %>% 
                  dplyr::filter(`Up or downregulated GO set` == .y)
              }) %>% 
  purrr:: map2(.x = .,
               .y = Filter,
               ~{
                 List = .x %>% 
                   dplyr::filter(str_detect(Description,.y)) %>%
                   dplyr::filter(!is.na(Category)) %>% 
                   dplyr::filter(!is.na(`_MEMBER_MyList`)) %>% 
                   dplyr::select(`Up or downregulated GO set`,
                                 Description,
                                 LogP,
                                 Hits,
                                 Category) %>% 
                   list()
                 names(List) = .x %>% dplyr::pull(`Up or downregulated GO set`) %>% unique()
                 List
               }) %>% 
  purrr::flatten()

All_metascape_pathways %>% 
  openxlsx::write.xlsx("Metascape/Description plot related gene list.xlsx")

All_metascape_pathways <- 
  All_metascape_pathways %>% 
  purrr::reduce(bind_rows)


# Arathi-list - Plotting description and p-values only ----
All_metascape_pathways %>% 
  dplyr::arrange(desc(`Up or downregulated GO set`),desc(LogP)) %>% 
  dplyr::mutate(Description = as_factor(Description),
                `Up or downregulated GO set` = as_factor(`Up or downregulated GO set`)) %>% 
  ggplot(aes(y = Description, x = LogP)) +
  geom_col(aes(fill = `Up or downregulated GO set`)) +
  facet_wrap(`Up or downregulated GO set`~.,nrow = 2, scales = "free") +
  theme_minimal() +
  scale_fill_manual(values = c("brown","green4")) +
  theme(
    plot.background = element_rect(fill = "white"),
    axis.line = element_line(color = "gray"),
    strip.background = element_rect(fill = "white"),
    axis.text = element_text(family = "Helvetica"),
    legend.title = element_blank(),
    legend.position = "none",
    legend.text = element_text(family = "Helvetica"),
    axis.title.x = element_text(family = "Helvetica"),
    axis.title.y = element_blank()
  ) + 
  ggsave("Metascape/Description plot.png",
         width = 10,
         height = 5)
# Arathi-list - preparing CPM values for heatmap ----
Arathi_DEGs_CPM_z_score <-
  All_metascape_heatmap %>% 
  dplyr::group_by_all() %>% 
  dplyr::group_split() %>% 
  purrr::map(~{
    .x %>% 
      dplyr::pull(Hits) %>% 
      str_split("\\|") %>% 
      unlist %>% 
      tibble::enframe(name = NULL, value = "Gene name") %>% 
      dplyr::bind_cols(
        .x %>% 
          dplyr::select(-Hits)
      ) %>% 
      dplyr::inner_join(DEGs_CPM_z_score)
  })
# Arathi-list - Renaming columns and rearranging order ----
Arathi_DEGs_CPM_z_score <- 
  Arathi_DEGs_CPM_z_score %>% 
  purrr::map(~{
    .x %>% 
      dplyr::select(matches("WT"),matches("KO"),1:ncol(.)) %>% 
      dplyr::rename(WT1 = WT5,
                    WT2 = WT6,
                    WT3 = WT7,
                    
                    KO2 = KO3,
                    KO3 = KO4)
  })
# Arathi-list - Function to plot heatmaps ----
F01_plot_heatmap <- function(Data, Pluck, Width = 5, Height = 15){
  data = Data %>%
    purrr::pluck(Pluck)
  # prepare matrix ----
  mat <- 
    data %>% 
    dplyr::select(matches("KO|WT")) %>%
    as.matrix()
  
  rownames(mat) <- data %>% 
    dplyr::pull(`Gene name`)
  
  # prepare highlight values ----
  # Rows to highlight
  Highlight_rows <- "Sirt2"
  
  # Set stylings for row names and make our selected rows unique
  Highlight_rows_id <- which(rownames(mat) %in% Highlight_rows)
  Fontsizes <- rep(10, nrow(mat))
  Fontsizes[Highlight_rows_id] <- 14
  Fontcolors <- rep('black', nrow(mat))
  Fontcolors[Highlight_rows_id] <- 'blue3'
  Fontfaces <- rep('plain',nrow(mat))
  Fontfaces[Highlight_rows_id] <- 'bold'
  
  # Create text annotation object for displaying row names
  Row_annotation <- rowAnnotation(rows = anno_text(rownames(mat),
                                                   gp = gpar(
                                                     fontsize = Fontsizes,
                                                     fontface = Fontfaces,
                                                     col = Fontcolors)))
  
  # Plot Heatmap ----
  col_fun = colorRamp2(c(-2, 0, 2), c("green","black","red"),transparency = 0)
  name <-
    str_c(
      "Metascape/Heatmap of ",
      (data %>% 
         dplyr::pull(`Up or downregulated GO set`) %>%
         unique()),
      " DEGs belonging to the GO: ",
      (data %>%
         dplyr::pull(Description) %>%
         unique),
      ".png"
    ) %>% 
    str_replace("6\\/M","6 or M")
  png(name,
      width = Width,
      height = Height,
      units = "in",
      res = 1200)
  Heatmap_GR <- ComplexHeatmap::Heatmap(mat,
                                        col = col_fun,
                                        name = "Z-score of CPM",
                                        border = T,
                                        show_row_names = F,
                                        cluster_columns = F,
                                        column_names_rot = 0,
                                        column_names_centered = T,
                                        column_names_gp = gpar(fontsize = 20, fontfamily = "Times"),
                                        row_names_gp = gpar(fontsize = 2, fontfamily = "Times"),
                                        column_gap = unit(0.4,"mm"),
                                        row_gap = unit(3,"mm"),
                                        # row_names_gp = gpar(fontfamily = "Times"),
                                        right_annotation = Row_annotation,
                                        heatmap_legend_param = list(
                                          title_gp = gpar(fontfamily = "Times"),
                                          lables_gp = gpar(fontfamily = "Times"),
                                          legend_direction = "horizontal",
                                          title_position = "leftcenter",
                                          legend_gp = gpar(fontsize = 2, fontfamily = "Times")
                                        ))
  
  draw(Heatmap_GR, heatmap_legend_side = "top")
  dev.off()
}

# Run this only to see how all of them will look at default dimensions
# 1:length(Arathi_DEGs_CPM_z_score) %>% 
#   purrr::map(~{
#     F01_plot_heatmap(Arathi_DEGs_CPM_z_score, .x)
#   })


# Arathi-list - Plotting custom sized heatmaps ----

# Apopotic signalling pathway
F01_plot_heatmap(Arathi_DEGs_CPM_z_score, 6)
# Autophagy
F01_plot_heatmap(Arathi_DEGs_CPM_z_score, 7)
# MAPK6/MAPK4
F01_plot_heatmap(Arathi_DEGs_CPM_z_score, 5, Height = 8)
# Skeletal muscle contraction
F01_plot_heatmap(Arathi_DEGs_CPM_z_score, 9, Height = 5)
# muscle atrophy
F01_plot_heatmap(Arathi_DEGs_CPM_z_score, 8, Height = 2)
# mTOR
F01_plot_heatmap(Arathi_DEGs_CPM_z_score, 3, Height = 5)
# AMPK
F01_plot_heatmap(Arathi_DEGs_CPM_z_score, 1, Height = 5)
# ----
# Arathi-list - heatmap of genes in the chosen plots - downregulated only - Plot heatmap of all DEGs - prepare matrix ----
mat <- Arathi_DEGs_CPM_z_score %>%
  dplyr::filter(`Up or downregulated GO set` == "Downregulated") %>% 
  dplyr::select(-`Gene name`,-`Up or downregulated GO set`) %>%
  as.matrix()

rownames(mat) <- Arathi_DEGs_CPM_z_score %>%
  dplyr::filter(`Up or downregulated GO set` == "Downregulated") %>% 
  dplyr::pull(`Gene name`)

# Arathi-list - heatmap of genes in the chosen plots - downregulated only - Plot heatmap of all DEGs - prepare highlight values 
----
  # Rows to highlight
  Highlight_rows <- "Sirt2"

# Set stylings for row names and make our selected rows unique
Highlight_rows_id <- which(rownames(mat) %in% Highlight_rows)
Fontsizes <- rep(10, nrow(mat))
Fontsizes[Highlight_rows_id] <- 14
Fontcolors <- rep('black', nrow(mat))
Fontcolors[Highlight_rows_id] <- 'blue3'
Fontfaces <- rep('plain',nrow(mat))
Fontfaces[Highlight_rows_id] <- 'bold'

# Create text annotation object for displaying row names
Row_annotation <- rowAnnotation(rows = anno_text(rownames(mat),
                                                 gp = gpar(
                                                   fontsize = Fontsizes,
                                                   fontface = Fontfaces,
                                                   col = Fontcolors)))

# Arathi-list - heatmap of genes in the chosen plots - downregulated only - Plot heatmap of all DEGs - Plot Heatmap ----
col_fun = colorRamp2(c(-2, 0, 2), c("green","black","red"),transparency = 0)
png("Metascape/Heatmap of Downregulated GO enriched DEGs.png",
    width = 5,
    height = 15,
    units = "in",
    res = 1200)
Heatmap_GR <- ComplexHeatmap::Heatmap(mat,
                                      col = col_fun,
                                      name = "Z-score of CPM",
                                      border = T,
                                      show_row_names = F,
                                      # row_split = Readcounts_CPM_GR %>%
                                      #   pull(`GO Term`) %>%
                                      #   as.factor(),
                                      cluster_columns = F,
                                      column_names_rot = 0,
                                      column_names_centered = T,
                                      column_names_gp = gpar(fontsize = 20, fontfamily = "Times"),
                                      row_names_gp = gpar(fontsize = 2, fontfamily = "Times"),
                                      column_gap = unit(0.4,"mm"),
                                      row_gap = unit(3,"mm"),
                                      # row_names_gp = gpar(fontfamily = "Times"),
                                      right_annotation = Row_annotation,
                                      heatmap_legend_param = list(
                                        title_gp = gpar(fontfamily = "Times"),
                                        lables_gp = gpar(fontfamily = "Times"),
                                        legend_direction = "horizontal",
                                        title_position = "leftcenter",
                                        legend_gp = gpar(fontsize = 2, fontfamily = "Times")
                                      ))

draw(Heatmap_GR, heatmap_legend_side = "top")
dev.off()

# Arathi-list - heatmap of genes in the chosen plots - downregulated only - Plot heatmap of all DEGs - cleaning up ----
remove(mat, Heatmap_GR, Fontcolors, Fontfaces, Fontsizes, Highlight_rows, Highlight_rows_id, Row_annotation)
# ----
# Bhoomika-list - data loading ----

# There is some discrepancy? Figure it out 3799 in mine vs 3783 in Bhoomika sheet I sent when?
# There is some discrepancy? 
# For now not using the Bhoomika list

Bhoomika_list <- '/Ankit/Ankit_RNA_seq_June_2024/Metascape/Bhoomika heat maps KO1-KO3-KO4 vs WT5-WT6-WT7 (3783 DEGs - 1607 down and 2176 up) - split by up or downregulated into two sheets.xlsx' %>% 
  purrr::map(~{
    Path = .x
    
    Path %>% 
      readxl::excel_sheets() %>% 
      str_subset("down|up", negate = T) %>%
      purrr::map2(.x = .,
                  .y = c(4,5,4),
                  ~{
                    Path %>% 
                      readxl::read_excel(sheet = .x,skip = .y) %>% 
                      dplyr::mutate(Description = .x)
                  })
  }) %>% 
  purrr::flatten() %>% 
  purrr::map(~{
    .x %>% 
      dplyr::filter(!is.na(`Gene name`))
  })

Bhoomika_list <- 
  Bhoomika_list %>% 
  purrr::map(~{
    .x %>% 
      dplyr::select(Description, `Gene name`) %>% 
      dplyr::left_join(DEGs_CPM_z_score)
  }
  )

Bhoomika_list <- Bhoomika_list %>% 
  purrr::map(~{
    .x %>% 
      dplyr::rename_all(~c("Description", "Gene name",
                           "KO1", "KO2", "KO3", 
                           "WT1", "WT2", "WT3")) %>% 
      dplyr::select(Description, `Gene name`, matches("WT"), matches("KO"))
  })

Bhoomika_list %>% 
  pluck(1) %>% 
  view


# Bhoomika-list - prepare matrix ----
mat <- Bhoomika_list %>%
  purrr::map(~{
    .x %>% 
      dplyr::select(-`Gene name`,-`Description`) %>%
      as.matrix()
  })

mat <- Bhoomika_list %>%
  purrr::map2(.x = .,
              .y = mat,
              ~{
                Matrix = .y
                rownames(Matrix) = .x %>% 
                  dplyr::pull(`Gene name`)
                Matrix
              })

names(mat) <- Bhoomika_list %>% 
  purrr::map(~{
    .x %>% 
      dplyr::pull(Description) %>% 
      unique
  }) %>% 
  unlist

# Bhoomika-list - plot heatmap ----

F99_Bhoomika_heatmaps <- function(Mat, Png){
  mat <- Mat
  
  # Rows to highlight
  Highlight_rows <- "Sirt2"
  
  # Set stylings for row names and make our selected rows unique
  Highlight_rows_id <- which(rownames(mat) %in% Highlight_rows)
  Fontsizes <- rep(10, nrow(mat))
  Fontsizes[Highlight_rows_id] <- 14
  Fontcolors <- rep('black', nrow(mat))
  Fontcolors[Highlight_rows_id] <- 'blue3'
  Fontfaces <- rep('plain',nrow(mat))
  Fontfaces[Highlight_rows_id] <- 'bold'
  
  # Create text annotation object for displaying row names
  Row_annotation <- rowAnnotation(rows = anno_text(rownames(mat),
                                                   gp = gpar(
                                                     fontsize = Fontsizes,
                                                     fontface = Fontfaces,
                                                     col = Fontcolors)))
  
  # Bhoomika-list - Plot Heatmap 
  col_fun = colorRamp2(c(-2, 0, 2), c("green","black","red"),transparency = 0)
  png(str_c("Metascape/Heatmap of ",Png,".png"),
      width = 5,
      height = 15,
      units = "in",
      res = 1200)
  Heatmap_GR <- ComplexHeatmap::Heatmap(mat,
                                        col = col_fun,
                                        name = "Z-score of CPM",
                                        border = T,
                                        show_row_names = F,
                                        # row_split = Readcounts_CPM_GR %>%
                                        #   pull(`GO Term`) %>%
                                        #   as.factor(),
                                        cluster_columns = F,
                                        column_names_rot = 0,
                                        column_names_centered = T,
                                        column_names_gp = gpar(fontsize = 20, fontfamily = "Times"),
                                        row_names_gp = gpar(fontsize = 2, fontfamily = "Times"),
                                        column_gap = unit(0.4,"mm"),
                                        row_gap = unit(3,"mm"),
                                        # row_names_gp = gpar(fontfamily = "Times"),
                                        right_annotation = Row_annotation,
                                        heatmap_legend_param = list(
                                          title_gp = gpar(fontfamily = "Times"),
                                          lables_gp = gpar(fontfamily = "Times"),
                                          legend_direction = "horizontal",
                                          title_position = "leftcenter",
                                          legend_gp = gpar(fontsize = 2, fontfamily = "Times")
                                        ))
  
  draw(Heatmap_GR, heatmap_legend_side = "top")
  dev.off()
}

mat %>% 
  purrr::map2(.x = .,
              .y = names(mat),
              ~{
                F99_Bhoomika_heatmaps(.x, .y)
              })

# ----
# ----