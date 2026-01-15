## simple-ARG-duplication-analysis.R by Rohan Maddamsetti.
## analyse the distribution of antibiotic resistance genes (ARGs)
## on chromosomes versus plasmids in fully-sequenced genomes and plasmids
## in the NCBI RefSeq database (dated March 26 2021).

library(tidyverse)
library(cowplot)
library(data.table)
library(stringi)
library(ggalluvial)

################################################################################
## Regular expressions used in this analysis.

## antibiotic-specific keywords.
chloramphenicol.keywords <- "chloramphenicol|Chloramphenicol"
tetracycline.keywords <- "tetracycline efflux|Tetracycline efflux|TetA|Tet(A)|tetA|tetracycline-inactivating"
MLS.keywords <- "macrolide|lincosamide|streptogramin"
multidrug.keywords <- "Multidrug resistance|multidrug resistance|antibiotic resistance"
beta.lactam.keywords <- "lactamase|LACTAMASE|beta-lactam|oxacillinase|carbenicillinase|betalactam\\S*"
glycopeptide.keywords <- "glycopeptide resistance|VanZ|vancomycin resistance|VanA|VanY|VanX|VanH|streptothricin N-acetyltransferase"
polypeptide.keywords <- "bacitracin|polymyxin B|phosphoethanolamine transferase|phosphoethanolamine--lipid A transferase"
diaminopyrimidine.keywords <- "trimethoprim|dihydrofolate reductase|dihydropteroate synthase"
sulfonamide.keywords <- "sulfonamide|Sul1|sul1|sulphonamide"
quinolone.keywords <- "quinolone|Quinolone|oxacin|qnr|Qnr"
aminoglycoside.keywords <- "Aminoglycoside|aminoglycoside|streptomycin|Streptomycin|kanamycin|Kanamycin|tobramycin|Tobramycin|gentamicin|Gentamicin|neomycin|Neomycin|16S rRNA (guanine(1405)-N(7))-methyltransferase|23S rRNA (adenine(2058)-N(6))-methyltransferase|spectinomycin 9-O-adenylyltransferase|Spectinomycin 9-O-adenylyltransferase|Rmt"
macrolide.keywords <- "macrolide|ketolide|Azithromycin|azithromycin|Clarithromycin|clarithromycin|Erythromycin|erythromycin|Erm|EmtA"
antimicrobial.keywords <- "QacE|Quaternary ammonium|quaternary ammonium|Quarternary ammonium|quartenary ammonium|fosfomycin|ribosomal protection|rifampin ADP-ribosyl|azole resistance|antimicrob\\S*"

antibiotic.keywords <- paste(chloramphenicol.keywords, tetracycline.keywords, MLS.keywords, multidrug.keywords,
    beta.lactam.keywords, glycopeptide.keywords, polypeptide.keywords, diaminopyrimidine.keywords,
    sulfonamide.keywords, quinolone.keywords, aminoglycoside.keywords, macrolide.keywords, antimicrobial.keywords, sep="|")


################################################################################
## Functions

## return the first column for several tables.
## shows the number of isolates in each category.
make.isolate.totals.col <- function(gbk.annotation) {
    isolate.totals <- gbk.annotation %>%
        group_by(Annotation) %>%
        summarize(total_isolates = n()) %>%
        arrange(desc(total_isolates))
    return(isolate.totals)
}


read.LLM.gbk.annotation.csv <- function(gbk.annotation.path, ground.truth.gbk.annotation) {
    
    ## only select Annotation_Accession, Organism, Strain, Genus columns from ground_truth_gbk.annotation.
    relevant.ground.truth <- ground.truth.gbk.annotation %>%
        select(Annotation_Accession, Organism, Strain, Genus)
    
    gbk.annotation.path %>%
        read.csv() %>%
        as_tibble() %>%
        ## filter based on ground truth genomes
        inner_join(relevant.ground.truth) %>%
        ## refer to NA annotations as "Unannotated".
        mutate(Annotation = replace_na(Annotation,"Unannotated")) %>%
        mutate(Annotation = replace(Annotation, Annotation == "NoAnnotation", "Unannotated")) %>%
        ## collapse Annotations into a smaller number of categories as follows:
        ## Marine, Freshwater --> Water
        ## Sediment, Soil, Terrestrial --> Earth
        ## Plants, Agriculture, Animals --> Plants & Animals
        ## Anthropogenic -> Human-impacted
        mutate(Annotation = replace(Annotation, Annotation == "Marine", "Water")) %>%
        mutate(Annotation = replace(Annotation, Annotation == "Freshwater", "Water")) %>%
        mutate(Annotation = replace(Annotation, Annotation == "Sediment", "Earth")) %>%
        mutate(Annotation = replace(Annotation, Annotation == "Soil", "Earth")) %>%
        mutate(Annotation = replace(Annotation, Annotation == "Terrestrial", "Earth")) %>%
        mutate(Annotation = replace(Annotation, Annotation == "Plants", "Plants & Animals")) %>%
        mutate(Annotation = replace(Annotation, Annotation == "Agriculture", "Plants & Animals")) %>%
        mutate(Annotation = replace(Annotation, Annotation == "Animals", "Plants & Animals")) %>%
        mutate(Annotation = replace(Annotation, Annotation == "Anthropogenic", "Human-impacted")) %>%
        ## And now remove all Unannotated genomes, since these are not analyzed
        ## at all in this first paper.
        filter(Annotation != "Unannotated") %>%
        ## and remove any strains (although none should fall in this category)
        ## that were not annotated by annotate-ecological-category.py.
        filter(Annotation != "blank")
}


##########################################################################
## Functions for Figure 2.

## See Wikipedia reference:
## https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval

## Make Z-distributed confidence intervals for the fraction of isolates with
## duplicated ARGs (panel A),
## the fraction of isolates with single-copy ARGs (panel B),
## the fraction of isolates with duplicated genes (panel C).

## Count data for isolates with duplicated ARGs
## goes into Supplementary Table S1.

calc.isolate.confints <- function(df) {
    df %>%
        ## use the normal approximation for binomial proportion conf.ints
        mutate(se = sqrt(p*(1-p)/total_isolates)) %>%
        ## See Wikipedia reference:
        ## https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
        mutate(Left = p - 1.96*se) %>%
        mutate(Right = p + 1.96*se) %>%
        ## truncate confidence limits to interval [0,1].
        rowwise() %>% mutate(Left = max(0, Left)) %>%
        rowwise() %>% mutate(Right = min(1, Right)) %>%
        ## Sort every table by the total number of isolates.
        arrange(desc(total_isolates))
}


make.TableS1 <- function(gbk.annotation, duplicate.ARGs) {

    ## count the number of isolates with duplicated ARGs in each category.
    ARG.category.counts <- duplicate.ARGs %>%
        ## next two lines is to count isolates rather than genes
        select(Annotation_Accession, Annotation) %>%
        distinct() %>%
        count(Annotation, sort = TRUE) %>%
        rename(isolates_with_duplicated_ARGs = n)
    
    ## join columns to make Table S1.
    TableS1 <- make.isolate.totals.col(gbk.annotation) %>%
        left_join(ARG.category.counts) %>%
        mutate(isolates_with_duplicated_ARGs =
                   replace_na(isolates_with_duplicated_ARGs,0)) %>%
        mutate(p = isolates_with_duplicated_ARGs/total_isolates) %>%
        calc.isolate.confints()
    
    return(TableS1)
}


make.confint.figure.panel <- function(Table, order.by.total.isolates, title,
                                      no.category.label = FALSE) {    
    Fig.panel <- Table %>%
        mutate(Annotation = factor(
                   Annotation,
                   levels = rev(order.by.total.isolates))) %>%
        ggplot(aes(y = Annotation, x = p)) +
        geom_point(size=1) +
        ylab("") +
        xlab("Proportion of Isolates") +
        theme_classic() +
        ggtitle(title) +
        ## plot CIs.
        geom_errorbar(aes(xmin=Left,xmax=Right), height=0.2, linewidth=0.2, orientation = "y")
    
    if (no.category.label)
        Fig.panel <- Fig.panel +
            theme(axis.text.y=element_blank())
    
    return(Fig.panel)
}


make.TableS2 <- function(gbk.annotation, singleton.ARGs) {

## count the number of isolates with singleton AR genes in each category.
    ARG.category.counts <- singleton.ARGs %>%
        ## next two lines is to count isolates rather than genes
        select(Annotation_Accession, Annotation) %>%
        distinct() %>%
        group_by(Annotation) %>%
        summarize(isolates_with_singleton_ARGs = n()) %>%
        arrange(desc(isolates_with_singleton_ARGs))
    gc() ## free memory.
    
    ## join columns to make Table S2.
    TableS2 <- make.isolate.totals.col(gbk.annotation) %>%
        left_join(ARG.category.counts) %>%
        mutate(isolates_with_singleton_ARGs =
                   replace_na(isolates_with_singleton_ARGs, 0)) %>%
        mutate(p = isolates_with_singleton_ARGs/total_isolates) %>%
        calc.isolate.confints()
    return(TableS2)
}


make.TableS3 <- function(gbk.annotation, duplicate.proteins) {
    ## count the number of isolates with duplicated genes in each category.
    category.counts <- duplicate.proteins %>%
        ## next two lines is to count isolates rather than genes
        select(Annotation_Accession, Annotation) %>%
        distinct() %>%
        group_by(Annotation) %>%
        summarize(isolates_with_duplicated_genes = n()) %>%
        arrange(desc(isolates_with_duplicated_genes))
    
    ## join columns to make Table S3.
    TableS3 <- make.isolate.totals.col(gbk.annotation) %>%
        left_join(category.counts) %>%
        mutate(isolates_with_duplicated_genes =
                   replace_na(isolates_with_duplicated_genes, 0)) %>%
        mutate(p = isolates_with_duplicated_genes/total_isolates) %>%
        calc.isolate.confints()
    return(TableS3)
}


################################################################################
## Set up the key data structures for the analysis.

ground.truth.gbk.annotation <- read.csv(
    "../data/Maddamsetti2024/FileS3-Complete-Genomes-with-Duplicated-ARG-annotation.csv") %>%
    as_tibble()

## This vector is used for ordering axes in figures and tables.
order.by.total.isolates <- make.isolate.totals.col(ground.truth.gbk.annotation)$Annotation

ground.truth.duplicate.proteins <- data.table::fread("../data/Maddamsetti2024/duplicate-proteins.csv",
                                                     drop="sequence") %>%
    ## now merge with gbk annotation.
    inner_join(ground.truth.gbk.annotation, by="Annotation_Accession")

ground.truth.duplicate.ARGs <- ground.truth.duplicate.proteins %>%
    filter(str_detect(product, antibiotic.keywords))

## Import llama3.2 ecological annotations
llama3.2.gbk.annotation <- read.LLM.gbk.annotation.csv(
    "../results/llama3.2_latest_gbk-annotation-table.csv",
    ground.truth.gbk.annotation)

llama3.2.duplicate.proteins <- data.table::fread("../data/Maddamsetti2024/duplicate-proteins.csv",
                                                 drop="sequence") %>%
    ## now merge with gbk annotation.
    inner_join(llama3.2.gbk.annotation, by="Annotation_Accession")

llama3.2.duplicate.ARGs <- llama3.2.duplicate.proteins %>%
    filter(str_detect(product, antibiotic.keywords))


################################################################################
## For speed, reload files on disk if they exist, otherwise recreate from scratch.

if (file.exists("../results/filtered-singleton-ARGs.csv")) {
    ## then load prefiltered data on disk.
    ground.truth.singleton.ARGs <- data.table::fread("../results/filtered-singleton-ARGs.csv")

} else { ## filter using regexes in R.

## get singleton proteins and filter.
ground.truth.singleton.proteins <- data.table::fread(
                                                   "../data/Maddamsetti2024/all-proteins.csv",
                                                   drop="sequence") %>%
    filter(count == 1) %>%
    inner_join(ground.truth.gbk.annotation, by="Annotation_Accession")
    
    ## THIS LINE IS SLOW. That said, it's as faster than my naive attempt to  filter with grep on the command-line...
    ground.truth.singleton.ARGs <- ground.truth.singleton.proteins[stri_detect_regex(product, antibiotic.keywords)]

    ## save data frames to disk for next time.
    write.csv(ground.truth.singleton.ARGs, "../results/filtered-singleton-ARGs.csv", row.names=FALSE, quote=FALSE)
}


####################################################################

## get LifestyleAnnotation, rename to Annotation for existing code to work nicely.
gbk.reannotation <- read.csv("../results/llama3.2_latest_Complete-Genomes-with-lifestyle-annotation.csv") %>%
    select(-hasDuplicatedARGs, -Annotation) %>%
    rename(Annotation = LifestyleAnnotation) %>%
    ## refer to NA annotations as "Unannotated".
    mutate(Annotation = replace_na(Annotation,"Unannotated")) %>%
    ## And now remove all Unannotated genomes, since these are not analyzed
    ## at all in this first paper.
    filter(Annotation != "Unannotated")

## and merge.
reannotated.singleton.ARGs <- ground.truth.singleton.ARGs %>%
    select(-Annotation) %>%
    left_join(gbk.reannotation)

reannotated.duplicate.ARGs <- ground.truth.duplicate.ARGs %>%
    select(-Annotation) %>%
    left_join(gbk.reannotation)

reannotated.duplicate.proteins <- ground.truth.duplicate.proteins %>%
    select(-Annotation) %>%
    left_join(gbk.reannotation)

##########################################################################
## Data structure for alluvial diagram of reannotations.

manual.annotation.df <- ground.truth.gbk.annotation %>%
    select(Annotation_Accession, Annotation) %>%
    rename(Original = "Annotation")

llama.annotation.df <- llama3.2.gbk.annotation %>%
        select(Annotation_Accession, Annotation) %>%
    rename(LLM = "Annotation")

lifestyle.reannotation.df <- gbk.reannotation %>%
        select(Annotation_Accession, Annotation) %>%
    rename(Lifestyle = "Annotation")

alluvial.plot.df <- manual.annotation.df %>%
    inner_join(llama.annotation.df) %>%
    inner_join(lifestyle.reannotation.df) %>%
    ## aggregate the counts in each category
    group_by(Original, LLM, Lifestyle) %>%
    summarize(Count = n()) %>%
    as.data.frame()

## validate the data.frame
is_alluvia_form(alluvial.plot.df, axes = 1:3, silent = TRUE)


################################################################################
## count the number of matches between the Original and LLM columns.

precision.df <- full_join(manual.annotation.df, llama.annotation.df) %>%
    group_by(Original, LLM) %>%
    summarize(Count = n())

## get the total number of annotated genomes.
## 18938 genomes here.
sum(precision.df$Count)

matched.precision.df <- precision.df %>%
    select(Original, LLM, Count) %>%
    filter(Original == LLM)

## 15062 match.
sum(matched.precision.df$Count)

## 15062/18938 = 0.795.
print(sum(matched.precision.df$Count)/sum(precision.df$Count))

## repeat the calculation, excluding Unannotated Genomes
precision.df2 <- precision.df %>%
    filter(!is.na(LLM))

## 18062 genomes here.
sum(precision.df2$Count)

matched.precision.df2 <- precision.df2 %>%
    select(Original, LLM, Count) %>%
    filter(Original == LLM)

## 15062 match.
sum(matched.precision.df2$Count)

## 15062/18062 = 0.834 precision.
print(sum(matched.precision.df2$Count)/sum(precision.df2$Count))


################################################################################
## Data structure for Figure 1BC:
## normal-approximation confidence intervals for the percentage
## of isolates with duplicated ARGs.
ground.truth.TableS1 <- make.TableS1(ground.truth.gbk.annotation, ground.truth.duplicate.ARGs)

llama3.2.TableS1 <- make.TableS1(llama3.2.gbk.annotation, llama3.2.duplicate.ARGs)

reannotated.TableS1 <- make.TableS1(gbk.reannotation, reannotated.duplicate.ARGs)

################################################################################
## Table S2. Control: does the distribution of ARG singletons
## (i.e. genes that have NOT duplicated) follow the distribution
## of sampled isolates?

## No categories are enriched with ARG singletons,
## as most isolates have a gene that matches an antibiotic keyword.
## Animal-host isolates are depleted (perhaps due to aphid bacteria isolates?)


## This data frame will be used for Figure 2A.
TableS2 <- make.TableS2(ground.truth.gbk.annotation, ground.truth.singleton.ARGs)

## This data frame will be used for Figure 2B.
reannotated.TableS2 <- make.TableS2(gbk.reannotation, reannotated.singleton.ARGs)

#########################################################################
## Table S3. Control: does the number of isolates with duplicate genes
## follow the sampling distribution of isolates?

## Data structure for Figure 2C
TableS3 <- make.TableS3(ground.truth.gbk.annotation, ground.truth.duplicate.proteins)

## Data structure for Figure 2D
reannotated.TableS3 <- make.TableS3(gbk.reannotation, reannotated.duplicate.proteins)


################################################################################
## CRITICAL TODO: DOUBLE-CHECK THAT THESE ARE CORRECT FOR THE PAPER!

## Save alluvial.plot.df and Tables S1, S2, S3 as Source Data.
##write.csv(alluvial.plot.df, "../results/Source-Data/Fig1A-Source-Data.csv", row.names=FALSE, quote=FALSE)
##write.csv(TableS1, "../results/Source-Data/Fig1BCD-Source-Data.csv", row.names=FALSE, quote=FALSE)
##write.csv(reannotated.TableS1, "../results/Source-Data/Fig1BCD-Source-Data.csv", row.names=FALSE, quote=FALSE)
##write.csv(TableS2, "../results/Source-Data/Fig2A-Source-Data.csv", row.names=FALSE, quote=FALSE)
##write.csv(reannotatedTableS2, "../results/Source-Data/Fig2B-Source-Data.csv", row.names=FALSE, quote=FALSE)
##write.csv(TableS3, "../results/Source-Data/Fig2C-Source-Data.csv", row.names=FALSE, quote=FALSE)
##write.csv(reannotated.TableS3, "../results/Source-Data/Fig2D-Source-Data.csv", row.names=FALSE, quote=FALSE)

## make Figures.
## Throughout, add special scales for panels as needed.

## This vector is used for ordering axes in figures and tables in the genomes reannotated by lifestyle
new.order.by.total.isolates <- make.isolate.totals.col(gbk.reannotation)$Annotation

## make Figure 1A.
Fig1A <- ggplot(alluvial.plot.df,
       aes(y = Count, axis1 = Original, axis2 = LLM, axis3 = Lifestyle)) +
    geom_alluvium(aes(fill = Original), width = 1/12) +
    geom_stratum(width = 1/12, fill = "black", color = "grey") +
    geom_label(stat = "stratum", size= 3.75, aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c("Original", "LLM", "Lifestyle"), expand = c(0.09, 0.09)) +
    scale_fill_brewer(type = "qual", palette = "Set1", name="Original annotation") +
    ggtitle("Rapid ecological reannotation of microbial genomes with a large language model") +
    theme_minimal() +
    theme(legend.position="top")

## add a bottom row of panels for D-ARGs.
Fig1B <- make.confint.figure.panel(ground.truth.TableS1, order.by.total.isolates, "D-ARGs, original\nannotation") +
    scale_x_continuous(breaks = c(0, 0.15), limits = c(0,0.16))

Fig1C <- make.confint.figure.panel(llama3.2.TableS1, order.by.total.isolates, "D-ARGs, llama3.2\nannotation") +
    scale_x_continuous(breaks = c(0, 0.15), limits = c(0,0.16))

Fig1D <- make.confint.figure.panel(reannotated.TableS1, new.order.by.total.isolates,
                                   "D-ARGs in genomes\nreannotated by lifestyle")

Fig1BCD <- plot_grid(Fig1B, Fig1C, Fig1D, labels=c("B", "C", "D"), nrow=1)

Fig1 <- plot_grid(Fig1A, Fig1BCD, labels=c("A", ""), nrow=2, rel_heights=c(3,1))
ggsave("../results/Fig1.pdf", Fig1, height=7.5, width=9)


Fig2A <- make.confint.figure.panel(TableS3, order.by.total.isolates,
                                   "All D-genes, original\nannotation")

Fig2B <- make.confint.figure.panel(TableS2, order.by.total.isolates,
                                   "S-ARGs, original\nannotation")

Fig2C <- make.confint.figure.panel(reannotated.TableS3, new.order.by.total.isolates,
                                   "All D-genes in genomes\nreannotated by lifestyle")

Fig2D <- make.confint.figure.panel(reannotated.TableS2, new.order.by.total.isolates,
                                   "S-ARGs in genomes\nreannotated by lifestyle")


Fig2 <- plot_grid(Fig2A, Fig2B, Fig2C, Fig2D,
                  labels=c('A', 'B', 'C', 'D'), nrow=2)    
ggsave("../results/Fig2.pdf", Fig2, height=4, width=7)
