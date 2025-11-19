## simple-ARG-duplication-analysis.R by Rohan Maddamsetti.
## analyse the distribution of antibiotic resistance genes (ARGs)
## on chromosomes versus plasmids in  fully-sequenced genomes and plasmids
## in the NCBI RefSeq database (dated March 26 2021).

library(tidyverse)
library(cowplot)
library(ggrepel)
library(data.table)

################################################################################

fancy_scientific <- function(x) {
    ## function for plotting better axis labels.
    ## see solution here for nice scientific notation on axes.
    ## https://stackoverflow.com/questions/10762287/how-can-i-format-axis-labels-with-exponents-with-ggplot2-and-scales
    ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scales::scientific_format()(x)))))
}

################################################################################
## Regular expressions used in this analysis.

## match MGE genes using the following keywords in the "product" annotation
transposon.keywords <- "IS|transpos\\S*|insertion|Tra[A-Z]|Tra[0-9]|tra[A-Z]|conjugate transposon|Transpos\\S*|Tn[0-9]|tranposase|Tnp|Ins|ins"
plasmid.keywords <- "relax\\S*|conjug\\S*|mob\\S*|plasmid|type IV|chromosome partitioning|chromosome segregation|Mob\\S*|Plasmid|Rep|Conjug\\S*"
phage.keywords <- "capsid|phage|Tail|tail|head|tape measure|antiterminatio|Phage|virus|Baseplate|baseplate|coat|entry exclusion"
other.HGT.keywords <- "Integrase|integrase|excision\\S*|exonuclease|recomb|toxin|restrict\\S*|resolv\\S*|topoisomerase|reverse transcrip|intron|antitoxin|toxin|Toxin|Reverse transcriptase|hok|Hok|competence|addiction"

MGE.keywords <- paste(transposon.keywords, plasmid.keywords, phage.keywords, other.HGT.keywords, sep="|")

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

antibiotic.or.MGE.keywords <- paste(MGE.keywords,antibiotic.keywords,sep="|")


categorize.as.MGE.ARG.or.other <- function(product) {
    if (is.na(product))
        return("Other function")
    else if (str_detect(product, antibiotic.keywords))
        return("ARG")
    else if (str_detect(product, MGE.keywords))
        return("MGE")
    else
        return("Other function")
}


################################################################################
## Set up the key data structures for the analysis:
## gbk.annotation, in particular.

gbk.annotation <- read.csv("../data/Maddamsetti2024/FileS3-Complete-Genomes-with-Duplicated-ARG-annotation.csv") %>%
  as_tibble()

#gbk.annotation <- read.csv(
#    "../results/computationally-annotated-gbk-annotation-table.csv") %>%
#    as_tibble() %>%
    ## filter based on QCed.genomes.
#    filter(Annotation_Accession %in% QCed.genomes$Annotation_Accession) %>%
    ## remove the overlaps with the Hawkey et al. (2022) clinical validation data.
#    filter(!(Annotation_Accession %in% Hawkey.overlaps.in.gbk.annotation.vec)) %>%
    ## refer to NA annotations as "Unannotated".
#    mutate(Annotation = replace_na(Annotation,"Unannotated")) %>%
    ## collapse Annotations into a smaller number of categories as follows:
    ## Marine, Freshwater --> Water
    ## Sediment, Soil, Terrestrial --> Earth
    ## Plants, Agriculture, Animals --> Plants & Animals
#    mutate(Annotation = replace(Annotation, Annotation == "Marine", "Water")) %>%
#    mutate(Annotation = replace(Annotation, Annotation == "Freshwater", "Water")) %>%
#    mutate(Annotation = replace(Annotation, Annotation == "Sediment", "Earth")) %>%
#    mutate(Annotation = replace(Annotation, Annotation == "Soil", "Earth")) %>%
#    mutate(Annotation = replace(Annotation, Annotation == "Terrestrial", "Earth")) %>%
#    mutate(Annotation = replace(Annotation, Annotation == "Plants", "Plants & Animals")) %>%
#    mutate(Annotation = replace(Annotation, Annotation == "Agriculture", "Plants & Animals")) %>%
#    mutate(Annotation = replace(Annotation, Annotation == "Animals", "Plants & Animals")) %>%
    ## get species name annotation from episome.database.
#    left_join(episome.database) %>%
    ## Annotate the genera.
#    mutate(Genus = stringr::word(Organism, 1)) %>%
    ## CRITICAL STEP: remove the NCBI_Nucleotide_Accession and SequenceType columns.
    ## This is absolutely critical, otherwise each row is duplicated for every
    ## chromosome and plasmid, breaking the invariant that each row refers to one sequence,
    ## when we add this annotation to duplicate.proteins and singleton.proteins.
#    select(-NCBI_Nucleotide_Accession, -SequenceType) %>%
    ## and we have to explicitly remove redundant rows now.
#    distinct() %>%
    ## And now remove all Unannotated genomes, since these are not analyzed
    ## at all in this first paper.
#    filter(Annotation != "Unannotated") %>%
    ## and remove any strains (although none should fall in this category)
    ## that were not annotated by annotate-ecological-category.py.
#    filter(Annotation != "blank")

## return the first column for several tables.
## shows the number of isolates in each category.
make.isolate.totals.col <- function(gbk.annotation) {
    isolate.totals <- gbk.annotation %>%
        group_by(Annotation) %>%
        summarize(total_isolates = n()) %>%
        arrange(desc(total_isolates))
    return(isolate.totals)
}


## This vector is used for ordering axes in figures and tables.
order.by.total.isolates <- make.isolate.totals.col(gbk.annotation)$Annotation


## read in duplicate proteins with sequences, using a separate file.
## I want the sequence column for the duplicate genes,
## but not for the singletons, to save memory.
duplicate.proteins <- read.csv("../data/Maddamsetti2024/duplicate-proteins.csv") %>%
    ## now merge with gbk annotation.
    inner_join(gbk.annotation) %>%
    mutate(Category = sapply(product, categorize.as.MGE.ARG.or.other))


## import the 37GB file containing all proteins, including singletons.
## I can save a ton of memory if I don't import the sequence column,
## and by using the data.table package for import.
singleton.proteins <- data.table::fread("../data/Maddamsetti2024/all-proteins.csv",
                                        drop="sequence") %>%
    filter(count == 1) %>%
    inner_join(gbk.annotation) %>%
    mutate(Category = sapply(product, categorize.as.MGE.ARG.or.other))
        
duplicate.ARGs <- duplicate.ARGs.by.keyword
duplicate.MGE.genes <- duplicate.MGE.genes.by.keyword

singleton.ARGs <- singleton.proteins %>%
    filter(str_detect(.$product, antibiotic.keywords))

singleton.MGE.genes <- singleton.proteins %>%
    filter(str_detect(.$product, MGE.keywords))

##########################################################################
## Code and data structures for Figure 4ABC.

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


## generic version of make.TableS1, for examining classes of genes other than
## antibiotic resistance genes.
make.IsolateEnrichmentTable <- function(gbk.annotation, duplicate.proteins, keywords) {
    ## count the number of isolates with duplicated genes of interest in each category.
    category.counts <- duplicate.proteins %>%
        filter(str_detect(.$product, keywords)) %>%
        ## next two lines is to count isolates rather than genes
        select(Annotation_Accession, Annotation) %>%
        distinct() %>%
        count(Annotation, sort = TRUE) %>%
        rename(isolates_with_duplicated_function = n)
    
    ## join columns to make the Table.
    Table <- make.isolate.totals.col(gbk.annotation) %>%
        left_join(category.counts) %>%
        mutate(isolates_with_duplicated_function =
                   replace_na(isolates_with_duplicated_function,0)) %>%
        mutate(p = isolates_with_duplicated_function/total_isolates) %>%
        calc.isolate.confints()
    return(Table)
}


make.IsolateEnrichmentControlTable <- function(gbk.annotation, singleton.proteins, keywords) {
    ## count the number of isolates with singleton genes of interest in each category.
    category.counts <- singleton.proteins %>%
        filter(str_detect(.$product, keywords)) %>%
        ## next two lines is to count isolates rather than genes
        select(Annotation_Accession, Annotation) %>%
        distinct() %>%
        count(Annotation, sort = TRUE) %>%
        rename(isolates_with_singleton_function = n)
    
    ## join columns to make the Table.
    Table <- make.isolate.totals.col(gbk.annotation) %>%
        left_join(category.counts) %>%
        mutate(isolates_with_singleton_function =
                   replace_na(isolates_with_singleton_function,0)) %>%
        mutate(p = isolates_with_singleton_function/total_isolates) %>%
        calc.isolate.confints()
    return(Table)
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
        geom_errorbarh(aes(xmin=Left,xmax=Right), height=0.2, size=0.2)
    
    if (no.category.label)
        Fig.panel <- Fig.panel +
            theme(axis.text.y=element_blank())
    
    return(Fig.panel)
}

## Data structure for Figure 4A:
## normal-approximation confidence intervals for the percentage
## of isolates with duplicated ARGs.
TableS1 <- make.TableS1(gbk.annotation, duplicate.ARGs)


######################
## Table S2. Control: does the distribution of ARG singletons
## (i.e. genes that have NOT duplicated) follow the distribution
## of sampled isolates?

## No categories are enriched with ARG singletons,
## as most isolates have a gene that matches an antibiotic keyword.
## Animal-host isolates are depleted (perhaps due to aphid bacteria isolates?)

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
                   replace_na(isolates_with_singleton_ARGs,0)) %>%
        mutate(p = isolates_with_singleton_ARGs/total_isolates) %>%
        calc.isolate.confints()
    return(TableS2)
}

## This data frame will be used for Figure 4B.
TableS2 <- make.TableS2(gbk.annotation, singleton.ARGs)

#########################################################################
## Table S3. Control: does the number of isolates with duplicate genes
## follow the sampling distribution of isolates?

## Most follow the expected distribution.
## however, isolates from animal-hosts are signficantly depleted
## in duplicate genes: FDR-corrected p = 0.0000314
## while isolates from anthropogenic environments are weakly enriched
## in multi-copy genes: FDR-corrected p = 0.0212.

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

## Data structure for Figure 4C.
TableS3 <- make.TableS3(gbk.annotation, duplicate.proteins)

            
################################################################################
## Save Tables S1, S2, and S3 as Source Data for Fig4ABC.
write.csv(TableS1, "../results/Source-Data/Fig4A-Source-Data.csv", row.names=FALSE, quote=FALSE)
write.csv(TableS2, "../results/Source-Data/Fig4B-Source-Data.csv", row.names=FALSE, quote=FALSE)
write.csv(TableS3, "../results/Source-Data/Fig4C-Source-Data.csv", row.names=FALSE, quote=FALSE)

## Finally -- make Figure 4ABC.
## Throughout, add special scales for Figure 4 but not for Supplementary Figure S14.
Fig4A <- make.confint.figure.panel(TableS1, order.by.total.isolates, "D-ARGs") +
    scale_x_continuous(breaks = c(0, 0.15), limits = c(0,0.16))


Fig4B <- make.confint.figure.panel(TableS2, order.by.total.isolates,
                                   "S-ARGs", no.category.label=TRUE) +
    scale_x_continuous(breaks = c(0.85, 1.0), limits = c(0.85,1))


Fig4C <- make.confint.figure.panel(TableS3, order.by.total.isolates,
                                   "All D-genes", no.category.label=TRUE) +
        scale_x_continuous(breaks = c(0.75, 0.95), limits = c(0.75, 1.0))


Fig4ABC.title <- title_theme <- ggdraw() +
    draw_label("Isolate-level analysis",fontface="bold")

Fig4ABC <- plot_grid(Fig4A, Fig4B, Fig4C, labels=c('A','B','C'),
                     rel_widths = c(1.5,1,1), nrow=1)

Fig4ABC.with.title <- plot_grid(Fig4ABC.title, Fig4ABC, ncol = 1, rel_heights = c(0.1, 1))


## show the plot
Fig4ABC.with.title




