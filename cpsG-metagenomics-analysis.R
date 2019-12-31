## cpsG-metagenomics-analysis.R by Rohan Maddamsetti.
## analysis of molecular evolution following cpsG and manB gene fusion events.

library(tidyverse)
library(DescTools) ## for G-test implementation in the multinomial selection test.

## get the lengths of all genes in REL606.
## This excludes genes in repetitive regions of the genome.
## See Section 4.3.1 "Removing mutations in repetitive regions of the genome"
## in Ben Good's LTEE metagenomics paper for more details.
## This filtering is done in my python script printEcoliIDs.py.
##Do by running:
##python printEcoliIDs.py -i ../data/REL606.7.gbk > ../results/REL606_IDs.csv.
REL606.genes <- read.csv('../results/REL606_IDs.csv',as.is=TRUE) %>%
mutate(gene_length=strtoi(gene_length))

## Order nonmutator pops, then hypermutator pops by converting Population to
## factor type and setting the levels.
nonmutator.pops <- c("Ara-5", "Ara-6", "Ara+1", "Ara+2", "Ara+4", "Ara+5")
hypermutator.pops <- c("Ara-1", "Ara-2", "Ara-3", "Ara-4", "Ara+3", "Ara+6")


## import genes in the table created by:
## conda activate ltee-metagenomix
## cd LTEE-metagenomic-repo
## python rohan-write-csv.py > ../results/LTEE-metagenome-mutations.csv
mutation.data <- read.csv(
    '../results/LTEE-metagenome-mutations.csv',
    header=TRUE,as.is=TRUE) %>%
    mutate(Generation=t0/10000) %>%
    ## This for changing the ordering of populations in plots.
    mutate(Population=factor(Population,levels=c(nonmutator.pops,hypermutator.pops)))


gene.mutation.data <- inner_join(mutation.data,REL606.genes)

## It turns out that some gene names map to multiple genes!!!
duplicate.genes <- gene.mutation.data %>%
    filter(Gene!='intergenic') %>%
    group_by(Gene) %>%
    summarize(checkme=length(unique(gene_length))) %>%
    filter(checkme>1)

## filter those duplicates.
gene.mutation.data <- filter(gene.mutation.data,
                             !(Gene %in% duplicate.genes$Gene))

########################################################################
## TODO: import and analyze cpsG and manB mutations in the 264 LTEE Genomes.
## Some of these mutations appear to be rare variants that were not sampled in
## Ben Good's metagenomics paper.
## These were downloaded from Jeff Barrick's shiny app:
## https://barricklab.org/shiny/LTEE-Ecoli.

cpsG.LTEE.genome.mutations <- read.csv("../results/cpsG-LTEE-Ecoli-query.csv",
                              header=TRUE,
                              as.is=TRUE) %>%
    rename(Gene='cpsG') %>%
    inner_join(REL606.genes)


#################################################################################

## For every gene, calculate my multinomial test for fit to a "neutral" model of
## hits per gene, based on the number of mutations in each population per mutation class.
## combine this in a dataframe.

## Then,
## Calculate the total density of mutations for each gene
## (for different classes of mutations),
## and combine in a dataframe.

## I want to use both sets of these statistics, on a per-gene basis, in order to
## see if any of these features relate to the rates at which various classes of
## mutations occur in different sets of genes.

################
## calculate the probability of a given configuration of parallel mutations
## in the metagenomes across populations.
## This seems to be a neat test for positive selection!
## Seems to give the same answer as results in Tenaillon Nature paper.
## TODO: compare to Supplementary Table S3 in Good et al. Nature paper.
gene.multinom.probability <- function(data, gene, mut.class.vec) {

    gene.data <- filter(data,Gene==gene, Annotation %in% mut.class.vec)
    
    gene.data.counts <- gene.data %>% group_by(Population) %>%
        summarize(muts=n())

    ## if there are no mutations in this gene, return NA.
    if(sum(gene.data.counts$muts) == 0) return(NA)
    
    ## calculate the bin size for each population as the fraction of mutations
    ## found in that population.
    population.probs <- data %>% group_by(Population) %>%
        summarize(total.muts=n()) %>% mutate(prob=total.muts/sum(total.muts))
  
    ## hack to get lengths of vectors matching for the GTest.
    full.pop.column <- data.frame(Population=population.probs$Population,
                                  stringsAsFactors=FALSE)
    stat.df <- full_join(full.pop.column,gene.data.counts,by='Population')
    stat.df[is.na(stat.df)] <- 0

    ## Use a G-test from the DescTools package because
    ## an exact multinomial test is too slow.
    GTest(x=stat.df$muts, p=population.probs$prob)$p.value
}

## draw every gene in the genome. calculate log probability and rank.
genes.vec <- unique(REL606.genes$Gene)
all.muts.pval.vec <- sapply(genes.vec, function(gene) {gene.multinom.probability(gene.mutation.data,gene,c("missense", "synonymous", "sv", "indel", "nonsense"))})
dN.pval.vec <- sapply(genes.vec, function(gene) {gene.multinom.probability(gene.mutation.data,gene,c("missense"))})
dS.pval.vec <- sapply(genes.vec,function(gene) {gene.multinom.probability(gene.mutation.data,gene,c("synonymous"))})
nonsense.sv.indel.pval.vec <- sapply(genes.vec, function(gene) {gene.multinom.probability(gene.mutation.data,gene,c("sv", "indel", "nonsense"))})

## IMPORTANT: KEEP NA VALUES.
multinom.result <- data.frame(Gene=genes.vec,
                              dN.pvalue=dN.pval.vec,
                              dS.pvalue=dS.pval.vec,
                              all.mut.pvalue=all.muts.pval.vec,
                              nonsense.sv.indel.pvalue=nonsense.sv.indel.pval.vec) %>%
    arrange(dN.pvalue) %>%
    mutate(index=1:n()) %>%
    mutate(dS.qvalue=p.adjust(dS.pvalue,'fdr')) %>%
    mutate(all.mut.qvalue=p.adjust(all.mut.pvalue,'fdr')) %>%
    mutate(nonsense.sv.indel.qvalue=p.adjust(nonsense.sv.indel.pvalue,'fdr')) %>%
    mutate(dN.qvalue=p.adjust(dN.pvalue,'fdr')) %>%
    ## keep the annotation in REL606.genes.
    full_join(REL606.genes)

## write multinomial selection test results to file.
write.csv(multinom.result,file='../results/multinomial_test_for_selection.csv')

## take a look at these distributions. The top genes are all under strong
## positive selection in the LTEE, as reported by Tenaillon et al.
multinom.plot <- ggplot(multinom.result,aes(x=index,y=-log10(dN.pvalue),label=Gene,color)) + geom_point() + theme_classic()

multinom.plot2 <- ggplot(arrange(multinom.result,all.mut.pvalue),aes(x=index,y=-log10(all.mut.pvalue),label=Gene,color)) + geom_point() + theme_classic()

## Now, let's examine the results of this test.

## interesting! cpsG is significant in terms of dS!
sig.multinom.dS <- filter(multinom.result,dS.pvalue<0.05)

## most of these genes are already known. Some aren't.
sig.multinom.dN <- filter(multinom.result,dN.pvalue<0.05)

sig.multinom.all.muts <- filter(multinom.result,all.mut.pvalue<0.05)
sig.multinom.nonsense.sv.indels <- filter(multinom.result,nonsense.sv.indel.pvalue<0.05)

sig.5percent.FDR <- filter(multinom.result,all.mut.qvalue<0.05)

## look at genes which are significant based on this test, but
## not based on parallelism alone, from Ben Good Table S3.
## these may be candidates for contingency.
Good.significant.genes <- read.csv('../ltee-metagenomics-paper/nature24287-s5.csv.txt')
contingency.candidates <- filter(sig.5percent.FDR,!(Gene %in% Good.significant.genes$Gene))
## only one new hit here: the significant dS gene cpsG.
## I queried this gene in the Good metagenomic data and in Jeff Barrick's
## LTEE-Ecoli shiny app. Extraordinarily interesting!
## cpsG looks like a case study for
## 1) historical contigency
## 2) epistasis
## 3) beneficial synonymous mutations.

## the last domain of cpsG was deleted by a parallel large deletion.
## also note cpsG mutations in the LTEE genomics data that are not found
## in the metagenomic data.
## These look like rare mutations that were missed in the metagenomics
## because they didn't rise to high enough frequency, that are driven by selection.


## Examine results of multinomial selection/contingency test in further detail.
## TODO: examine mutations in genes beyond cpsG. Anything worthwhile following up?

## NOTE: The multinom analysis here is based ONLY on distribution of dN.

## NOTE: rate of accumulation depends on the size of the gene set,
## when looking at the left and right tails.
## (I can't remember what this comment means-- take a look at the data to figure out.)

## also compare genes in the tails of the multinomial hit distribution.
## left tail is strong selection. right tail fits the null well.

## ODD! way more mutations in top 100 or top 200 genes.
## but way fewer in the top 50! Take a look at 1-50, and 51-100 genes.

multinom.1to50.genes <- multinom.result$Gene[1:50]
multinom.50to100.genes <- multinom.result$Gene[51:100]

multinom.1to100.genes <- multinom.result$Gene[1:100]
multinom.1to200.genes <- multinom.result$Gene[1:200]


multinom.1to50.mutation.data <- filter(mutation.data,
                                  Gene %in% multinom.1to50.genes)

multinom.50to100.mutation.data <- filter(mutation.data,
                                  Gene %in% multinom.50to100.genes)

multinom.1to100.mutation.data <- filter(mutation.data,
                                  Gene %in% multinom.1to100.genes)

multinom.1to200.mutation.data <- filter(mutation.data,
                                  Gene %in% multinom.1to200.genes)

## for the bottom, get the number of genes in the multinom.result
## (for now, including NA values).

num.multinom.genes <- length(multinom.result$Gene)

multinom.bottom50to1.genes <- multinom.result$Gene[num.multinom.genes-50:num.multinom.genes]
multinom.bottom100to50.genes <- multinom.result$Gene[(num.multinom.genes-100):(num.multinom.genes-50)]
multinom.bottom100to1.genes <- multinom.result$Gene[(num.multinom.genes-100):num.multinom.genes]
multinom.bottom200to1.genes <- multinom.result$Gene[num.multinom.genes-200:num.multinom.genes]

multinom.bottom200to1.mutation.data <- filter(mutation.data,
                                              Gene %in% multinom.bottom200to1.genes)

## length of tails of multinomial distibution.
multinom.1to50.length <- sum(filter(REL606.genes, Gene %in% multinom.1to50.genes)$gene_length, na.rm=TRUE)
multinom.50to100.length <- sum(filter(REL606.genes, Gene %in% multinom.50to100.genes)$gene_length, na.rm=TRUE)
multinom.1to100.length <- sum(filter(REL606.genes, Gene %in% multinom.1to100.genes)$gene_length, na.rm=TRUE)
multinom.1to200.length <- sum(filter(REL606.genes, Gene %in% multinom.1to200.genes)$gene_length, na.rm=TRUE)

multinom.bottom50to1.length <- sum(filter(REL606.genes, Gene %in% multinom.bottom50to1.genes)$gene_length, na.rm=TRUE)
multinom.bottom100to50.length <- sum(filter(REL606.genes, Gene %in% multinom.bottom100to50.genes)$gene_length, na.rm=TRUE)
multinom.bottom100to1.length <- sum(filter(REL606.genes, Gene %in% multinom.bottom100to1.genes)$gene_length, na.rm=TRUE)
multinom.bottom200to1.length <- sum(filter(REL606.genes, Gene %in% multinom.bottom200to1.genes)$gene_length, na.rm=TRUE)


#############################################################################
## examining the rates that genes in the left-hand tail and right-hand tail
## of the positive selection distribution get hit by mutations.
## TODO: what is the relationship, if any, between genes with dN (or not)
## and those with nonsense.indels.sv's, (or not), and rates of accumulation?
## Is there a relationship between genes under positive or purifying selection?

### make plots for the top, median, and bottom genes by multinomial test.
## IMPORTANT NOTE: THE UNDERLYING SETS OF GENES NEED TO BE THE SAME.

#### Plots for multinom1to50 genes.
c.multinom1to50.dN <- filter(multinom.1to50.mutation.data, Annotation=='missense') %>%
    calc.cumulative.muts(multinom.1to50.length)
c.multinom1to50.dS <- filter(multinom.1to50.mutation.data, Annotation=='synonymous') %>%
    calc.cumulative.muts(multinom.1to50.length)
c.multinom1to50.nonsense <- filter(multinom.1to50.mutation.data, Annotation=='nonsense') %>%
    calc.cumulative.muts(multinom.1to50.length)
c.multinom1to50.nonsense.indel.sv <- filter(multinom.1to50.mutation.data,
                                           Annotation %in% c("nonsense","indel","sv")) %>%
    calc.cumulative.muts(multinom.1to50.length)

## Now make the plot.
log.multinom1to50.plot <- plot.cumulative.muts(c.multinom1to50.dN, my.color="purple") %>%
    add.cumulative.mut.layer(c.multinom1to50.dS, my.color="green") %>%
    add.cumulative.mut.layer(c.multinom1to50.nonsense, my.color="gray") %>%
    add.cumulative.mut.layer(c.multinom1to50.nonsense.indel.sv, my.color="black")

log.multinom1to50.plot
ggsave(log.multinom1to50.plot,filename="../results/figures/log-multinom-1to50.pdf")

multinom1to50.plot <- plot.cumulative.muts(c.multinom1to50.dN, my.color="purple",logscale=FALSE) %>%
    add.cumulative.mut.layer(c.multinom1to50.dS, my.color="green",logscale=FALSE) %>%
    add.cumulative.mut.layer(c.multinom1to50.nonsense, my.color="gray",logscale=FALSE) %>%
    add.cumulative.mut.layer(c.multinom1to50.nonsense.indel.sv, my.color="black",logscale=FALSE)

log.multinom1to50.plot
ggsave(multinom1to50.plot,filename="../results/figures/multinom-1to50.pdf")




################################################################################
## cpsG is not essential in the KEIO collection.
## there is a hit to growth in MOPS minimal media when knocked out.
## notice how manB is absent from K-12.
KEIO.data <- read.csv("../data/KEIO_Essentiality.csv", header=TRUE,as.is=TRUE) %>%
    select(-JW_id)


################################################################################
## Statistical test for purifying selection.
## Look at accumulation of stars over time.
## in other words, look at the rates at which the mutations occur over time.

## To normalize, we need to supply the number of sites at risk
## (such as sum of gene length)

cumsum.per.pop.helper.func <- function(df) {
    df %>%
        arrange(t0) %>%
        group_by(Population,Generation) %>%
        summarize(count=n()) %>%
        mutate(cs=cumsum(count)) %>%
        ungroup() 
}

calc.cumulative.muts <- function(data, normalization.constant=NA) {

    ## if normalization.constant is not provided, then
    ## calculate based on gene length by default.
    if (is.na(normalization.constant)) {
        my.genes <- data %>% select(Gene,gene_length) %>% distinct()
        normalization.constant <- sum(my.genes$gene_length)
    }
    
    c.dat <- data %>%
        split(.$Population) %>%
        map_dfr(.f=cumsum.per.pop.helper.func) %>%
        mutate(normalized.cs=cs/normalization.constant) %>%
        ## remove any NA values.
        na.omit()
    return(c.dat)
}

###############################################################
## Examine the rate of cumulative accumulation of mutations over time in
## cpsG and manB compared to random genes.


plot.cumulative.muts <- function(mut.data,logscale=TRUE, my.color="black") {
    if (logscale) {
        p <- ggplot(mut.data,aes(x=Generation,y=log10(normalized.cs))) +
            ylim(-7,-2) +
            ylab('log[Cumulative number of mutations (normalized)]')
    } else {
        p <- ggplot(mut.data,aes(x=Generation,y=normalized.cs)) +
            ylim(0,0.003) +
            ylab('Cumulative number of mutations (normalized)')
    }
    p <- p +
        theme_classic() +
        geom_point(size=0.2, color=my.color) +
        geom_step(size=0.2, color=my.color) +
        facet_wrap(.~Population,scales='fixed') +
        xlab('Generations (x 10,000)') +
        xlim(0,6.3)
    return(p)
}

## take a ggplot object output by plot.cumulative.muts, and add an extra layer.
add.cumulative.mut.layer <- function(p, layer.df, my.color, logscale=TRUE) {
    if (logscale) {
        p <- p +
            geom_point(data=layer.df, aes(x=Generation,y=log10(normalized.cs)), color=my.color, size=0.2) +
            geom_step(data=layer.df, aes(x=Generation,y=log10(normalized.cs)), color=my.color, size=0.2)
        } else {
            p <- p +
                geom_point(data=layer.df, aes(x=Generation,y=normalized.cs), color=my.color, size=0.2) +
                geom_step(data=layer.df, aes(x=Generation,y=normalized.cs), color=my.color, size=0.2)
        }
    return(p)
}


c.mutations <- calc.cumulative.muts(gene.mutation.data)

c.mutation.plot <- plot.cumulative.muts(c.mutations, logscale=TRUE)
c.mutation.plot2 <- plot.cumulative.muts(c.mutations, logscale=FALSE)

##########################################################################
## Bootstrap a distribution for the distribution of the accumulation of stars
## over time for random subsets of genes.
## Say, 1000 or 10000 random subsets of genes.
## This idea could be used to bootstrap p-values for statistics, by
## comparing subsets of mutations in genes of interest to this random null model.
## Distributions greater or lower than all 1000 curves is significantly greater or lesser
## at a p = 0.001/2 (I think? check this calculation more rigorously.)

plot.random.subsets <- function(data, subset.size=300, N=1000,logscale=TRUE) {
    
  ## set up an empty plot then add random trajectories, one by one.
    my.plot <- ggplot(data) +
        theme_classic() +
        facet_wrap(.~Population,scales='fixed',nrow=4) +
        xlab('Generations (x 10,000)') +
        xlim(0,6.3)

    if (logscale) {
        my.plot <- my.plot +
            ylim(-7,-2) +
            ylab('log[Cumulative number of mutations (normalized)]')
    } else {
        my.plot <- my.plot +
            ylim(0,0.003) +
            ylab('Cumulative number of mutations (normalized)')
    }
    
    for (i in 1:N) {
        rando.genes <- sample(unique(data$Gene),subset.size)
        mut.subset <- filter(data,Gene %in% rando.genes)
        c.mut.subset <- calc.cumulative.muts(mut.subset)
        
        if (logscale) {
            my.plot <- my.plot +
                geom_point(data=c.mut.subset,aes(x=Generation,y=log10(normalized.cs)), color='gray',size=0.2,alpha = 0.1)                
                } else {
                    my.plot <- my.plot +
                        geom_point(data=c.mut.subset,aes(x=Generation,y=normalized.cs), color='gray',size=0.2,alpha = 0.1)
                }
    }
    return(my.plot)
}

## Base plots here of null distributions: add the data lines on top to compare.
log.all.rando.plot <- plot.random.subsets(gene.mutation.data, subset.size=1, logscale=TRUE)
log.sv.indel.nonsen.rando.plot <- plot.random.subsets(sv.indel.nonsense.gene.mutation.data,
                                                      subset.size=1, logscale=TRUE)

all.rando.plot <- plot.random.subsets(gene.mutation.data, subset.size=1, logscale=FALSE)
sv.indel.nonsen.rando.plot <- plot.random.subsets(sv.indel.nonsense.gene.mutation.data,
                                                  subset.size=1, logscale=FALSE)

## for rapid testing.
log.small.rando.plot <- plot.random.subsets(gene.mutation.data, subset.size=1,
                                            logscale=TRUE,N=10)
log.small.sv.indel.nonsen.rando.plot <- plot.random.subsets(sv.indel.nonsense.gene.mutation.data, subset.size=1, logscale=TRUE,N=10)


## TODO: plot curves for manB and cpsG in order to do hypothesis test for purifying
## selection on these loci.
