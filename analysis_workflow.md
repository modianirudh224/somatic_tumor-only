

# Somatic variant calling
This workflow is for the analysis of tumor-only unmatched Whole Exome Sequencing(WES) samples. 

The general idea is for doing two types of analysis:

1.  **Finding recurrent mutations in a given set of samples.**
> Using GATK tools, bwa, bcftools etc.
2. **Finding the most significantly mutated genes given a background mutation rate.**
> Using MutSigCV
  
  ## Setup GATK and bwa
  GATK 4.1.8.1 and bwa 0.7.17 is used herein example. Download the latest package as present at the time. 
 A tutorial to install all GATK related tools can be found [on the GATK website](https://gatk.broadinstitute.org/hc/en-us/articles/360041320571--How-to-Install-all-software-packages-required-to-follow-the-GATK-Best-Practices)


    wget https://github.com/broadinstitute/gatk/releases/download/4.1.8.1/gatk-4.1.8.1.zip
    unzip gatk-4.1.8.1.zip
   > Although the jars themselves cannot simply be added to your PATH, you can do so with the gatk wrapper script. Please look up instructions depending on the terminal shell you use; in bash the typical syntax is `export PATH="/path/to/gatk-package/:$PATH"` where /path/to/gatk-package/ is the path to the location of the gatk executable. Note that the jars must remain in the same directory as gatk for it to work. Be sure to include the final / in your path.
   
    wget https://excellmedia.dl.sourceforge.net/project/bio-bwa/bwa-0.7.17.tar.bz2
    bzip2 -d bwa-0.7.17.tar.bz2
    tar -xvf bwa-0.7.17.tar
    cd bwa-0.7.17/
    make
    cd ..
    
To add both the tools in path and apply the changes

    echo "export PATH='~/gatk-4.1.8.1:\$PATH' " >> .bash_profile
    echo "export PATH='~/bwa-0.7.17:\$PATH' " >> .bash_profile
    source ~/.bash_profile

## Getting all the initial references
Create a folder named ref to store all the references needed. We will be using gsutils to download resources from gatk cloud buckets. The resources files used are given below and if any other needed can be downloaded from the following :
1.  [Resource files for hg38 reference genome](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?pli=1&prefix=&forceOnObjectsSortingFiltering=false)
    
2.  [Resource files for calling something mutations using hg38 as a reference](https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38;tab=objects?project=broad-dsde-outreach&prefix=&forceOnObjectsSortingFiltering=false)

To download the references we need for now:
   
    gsutil -m cp -r gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict /path/to/ref/
    gsutil -m cp -r gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta* /path/to/ref/
    gsutil -m cp -r gs://genomics-public-data/resources/broad/hg38/v0/wgs_calling_regions.hg38.interval_list /path/to/ref/
    gsutil -m cp -r gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz* /path/to/ref/
    gsutil -m cp -r gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz /path/to/ref/
    gsutil -m cp -r gs://gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz* /path/to/ref/

## Workflow to call somatic variants in tumor only mode
>The rough plan is to map and clean short read sequences as given by GATK best practices [here](https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently). Followed by using a workflow from the GATK Mutect2 tutorial on [terra](https://app.terra.bio/#workspaces/help-gatk/GATKTutorials-Somatic/notebooks/launch/1-somatic-mutect2-tutorial.ipynb).

**1. Change directory to the folder containing the two paired-end WES fastq files R1 and R2. 
2. Considering the reads are named as X_R1.fastq.gz and X_R2.fastq.gz.**

Unzip both the files 

    gunzip X*.gz

Convert the fastq files to a single unaligned bam file. SM option adds the sample name X to the files in downstream 

    gatk FastqToSam -F1 X_R1.fastq -F2 X_R2.fastq -O X_fastqtosam.bam -SM X

Mark Illumina adapters and clip them 

    gatk MarkIlluminaAdapters \
    -I X_fastqtosam.bam \
    -O X_markilluminaadapters.bam \
    -M X_illumina_metrics.txt
 
 Convert the bam file back to a single fastq file. The options are selected as per the best practices mentioned on the GATK website given above. 

    gatk SamToFastq \
    -I X_markilluminaadapters.bam \
    -FASTQ X_samtofastq.fq \
    -CLIPPING_ATTRIBUTE XT \
    -CLIPPING_ACTION 2 \
    -INTERLEAVE true

Align this file to the reference genome (preferably hg38). Bwa supports multi-threading on a single node with 16/32 nodes being optimal. No of threads can be changed by the option -t.

    bwa mem -M -p -t 16 /path/to/ref/Homo_sapiens_assembly38.fasta X_samtofastq.fq > X_bwamem.sam

Merge the unaligned bam with aligned bam file, while clipping adapters and coordinate sorting the files ready to mark duplicates.

    gatk MergeBamAlignment \
    -R /path/to/ref/Homo_sapiens_assembly38.fasta \
    -UNMAPPED_BAM X_fastqtosam.bam \
    -ALIGNED_BAM X_bwamem.sam \
    -O X_mergebam.bam \
    -CREATE_INDEX true \
    -CLIP_ADAPTERS false \
    -MAX_INSERTIONS_OR_DELETIONS -1 \
    -ATTRIBUTES_TO_RETAIN XS

Mark the remove the duplicates arose due to PCR which might create a bias in results or cause a change in confidence levels.

    gatk MarkDuplicates \
    --CREATE_INDEX true \
    --REMOVE_DUPLICATES true \
    -I X_mergebam.bam \
    -O X_marked_duplicates.bam \
    -M X_marked_dup_metrics.txt 

Call mutect2 along with the panel of normals, germline resource and interval list.

    gatk Mutect2 \
    --native-pair-hmm-threads 16 \
    -R /path/to/ref/Homo_sapiens_assembly38.fasta \
    -I X_marked_duplicates.bam \
    -pon /path/to/ref/1000g_pon.hg38.vcf.gz \
    --germline-resource /path/to/ref/af-only-gnomad.hg38.vcf.gz \
    -L /path/to/ref/wgs_calling_regions.hg38.interval_list \
    -O X_somaticm2.vcf.gz \
    -bamout X_ashap_m2.bam

Calculate pileup summaries to summarizes counts of reads that support reference, alternate and other alleles for given sites.

    gatk GetPileupSummaries \
    -I X_marked_duplicates.bam \
    -V /path/to/ref/small_exac_common_3.hg38.vcf.gz \
    -L /path/to/ref/small_exac_common_3.hg38.vcf.gz \
    -O X_getpileupsummaries.table

Using pileup summaries created, calculate the fraction of reads coming from cross-sample calculation.

    gatk CalculateContamination \
    -I X_getpileupsummaries.table \
    -tumor-segmentation segments.table \
    -O X_calculatecontamination.table

Finally apply 20 hard filters to mark unwanted variants and retain somatic variants with high confidence.

    gatk FilterMutectCalls \
    -R /path/to/ref/Homo_sapiens_assembly38.fasta \
    -V X_somaticm2.vcf.gz \
    --contamination-table X_calculatecontamination.table \
    --stats X_somaticm2.vcf.gz.stats \
    --tumor-segmentation segments.table \
    -O X_somatic_oncefiltered.vcf.gz


Remove the filtered variants for downstream analysis and index it.

    gatk SelectVariants \
    --exclude-filtered \
    -R /path/to/ref/Homo_sapiens_assembly38.fasta \
    -V X_somatic_oncefiltered.vcf.gz \
    -O X_somatic_oncefiltered_passed.vcf.gz

    gatk IndexFeatureFile -F X_somatic_oncefiltered_passed.vcf.gz

Liftover to hg19 reference if needed for any tools that require, for example, MutSig in this case.

>Reference file needed for this process are chain files and the reference sequence of the genome that is being liftover to. 
>
>Chain files can be found at [http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/).
>hg19 reference genome at [http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz).

  
  

To download the chain file for liftover for hg38 to hg19 and hg19 reference genome

    cd /path/to/ref/
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

  Proceed with Liftover

    gatk LiftoverVcf \
    -I X_somatic_oncefiltered_passed.vcf.gz \
    -O X_somatic_oncefiltered_passed_hg19.vcf.gz \
    -C /path/to/ref/hg38ToHg19.over.chain 
    --REJECT reject_X.vcf \
    -R /path/to/ref/hg19.fa
    
Annotate variants using funcotator. 
The format can be chosen from VCF or MAF. MAF is needed for downstream analysis in MutSig. Reference version is chosen as hg19 as MutSig supports hg19 version only as of now. This can similarly be easily used with hg38 VCF files.

Download data sources needed to use funcotator. Change directory to ref to download it store in the folder containing all references.

    gatk FuncotatorDataSourceDownloader --somatic --validate-integrity --extract-after-download

  >This creates a folder called funcotator_dataSources.v1.6.20190124s (version at the time of writing). This folder contains multiple data sources. The only required data source is gencode and the rest can be deleted for the improving speed and reducing file size. For example, in our analysis using MutSig downstream we don't need any of those.

    gatk Funcotator \
    --variant X_somatic_oncefiltered_passed.vcf.gz \
    --reference /path/to/ref/hg19.fa \
    --ref-version hg19 \
    --data-sources-path /path/to/ref/funcotator_dataSources.v1.6.20190124s \
    --output X_somatic_oncefiltered_passed_hg19_funcotate.maf \
    --output-file-format MAF \
    --annotation-default Tumor_Sample_Barcode:X

  Simply add all mafs together to form a multi-sample maf. Copy mafs for all samples in a folder and concatenate.

    cat * > final_combined.maf

## Calling out recurrent mutations in multiple samples

We will be using bcftools for finding the intersection between multiple samples.

To setup bcftools:

    git clone git://github.com/samtools/bcftools.git
    cd bcftools
    make
>Add bcftools to the path like GATK and bwa in previous examples

To find the common mutations present in atleast n samples - 

*Copy all **X_somatic_oncefiltered_passed.vcf.gz** files into one folder and change the directory to that folder. Use bcftools to find the intersection.* 
>Given example is for common mutations in atleast 10 samples. 10 can be replaced with any number n for n samples.

    bcftools isec -c all -n 10 -p common_10 -O z *_somatic_oncefiltered_passed.vcf.gz

This will create a folder called `common_10` with one file corresponding to each sample containing only the mutations which were present in atleast 10 files. Change directory to this folder. 

Merge these using bcftools merge.

    bcftools merge -O v -o common_merged_10.vcf *.vcf.gz

**This will produce a VCF file with high confidence somatic variants common in atleast 10 samples among all the samples given.** 

1. Then these files were annotated using VEP web interface and sorted descending according to the CADD scores given. 
2. Manually, variants with mentions in COSMIC and low gnomad allele frequencies (<0.1) were checked for pathogenicity. 
3. SIFT and polyphen scores were also considered to consider their probable pathogenicity.

## Extracting cancer driving genes using MutSigCV
Download latest MutSig from the [website](https://software.broadinstitute.org/cancer/cga/mutsig).

    wget https://software.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutsig/MutSigCV_1.41.zip

1. The installation process is given [here](https://software.broadinstitute.org/cancer/cga/mutsig_download)
>This requires Matlab MCR version 2016a which can be downloaded and installed using using the [following manual](https://nl.mathworks.com/help/compiler/install-the-matlab-runtime.html)
2. The tutorial to run MutSigCV is given [here](https://software.broadinstitute.org/cancer/cga/mutsig_run)

The reference files needed are taken from the tutorial itself.  The following code downloads and unzips the required files.
    
    wget http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutsig/reference_files/chr_files_hg19.zip
    unzip chr_files_hg19.zip
    
    wget http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutsig/reference_files/mutation_type_dictionary_file.txt
    
    wget http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutsig/reference_files/exome_full192.coverage.zip
    unzip exome_full192.coverage.zip
    
    wget http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutsig/reference_files/gene.covariates.txt

After the setup is done. We have the given four reference files, a multi-sample maf file, MutSigCV installed and Matlab MCR installed. Given all these are are at their right places. Run the following :
 
    run_MutSigCV.sh <path_to_MCR> final_combined.maf exome_full192.coverage.txt gene.covariates.txt results mutation_type_dictionary_file.txt chr_files_hg19

The most important file generated will be the one called **results.sig_genes.txt**. 
>This file contains a list of most significantly mutated genes in descending order as a table, with one of the columns being the p-value. The lower the p-value, the higher the significance.

## Plotting the results using Maftools
 
 Maftools is a powerful R package used to summarise and plot results from multi-sample MAF files.
>Use the latest R version possible for the running maftools, otherwise some functions might not work.

>The example below uses a sample maf called *SQCC.maf*

Start the R program by typing R and install maftools

    BiocManager::install("maftools")

Read the file

    library(maftools)
    SQCC.maf = file.path('SQCC_combined.maf') 
    SQCC = read.maf(maf = SQCC.maf)

Make a set of genes recognised as top 20 in the file **results.sig_genes.txt created using MutSigCV**. The genes given are an example. The top 20 most significantly mutated genes can be viewed along with the header by `head -21 results.sig_genes.txt`,

    SQCC_genes = c('F11' ,'TP53' ,'CXCL9' ,'RBM17' ,'IMMP1L' ,'TSPAN2' ,'TRPC1' ,'BMPR1B' ,'CLUAP1' ,'C10orf88' ,'FGF14' ,'UBA2' ,'COX18' ,'GNG11' ,'RASL12' ,'HMGCLL1' ,'DARS' ,'YES1' ,'ST6GAL1' ,'GPR137C')

Create a subset of the combined maf with entries only from these genes.

    SQCC_top = subsetMaf(maf = SQCC, genes = c('F11' ,'TP53' ,'CXCL9' ,'RBM17' ,'IMMP1L' ,'TSPAN2' ,'TRPC1' ,'BMPR1B' ,'CLUAP1' ,'C10orf88' ,'FGF14' ,'UBA2' ,'COX18' ,'GNG11' ,'RASL12' ,'HMGCLL1' ,'DARS' ,'YES1' ,'ST6GAL1' ,'GPR137C'))

Read the table in results.sig_genes.txt

    SQCC.mutsig = file.path('results.sig_genes.txt')
    SQCC.mutsig = data.table::fread(input = SQCC.mutsig)[,.(gene, p)]
    SQCC.mutsig[,p := -log10(p)]

Create a summary plot of all genes

    pdf('SQCC_summary_complete.pdf')
    plotmafSummary(maf = SQCC, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
    dev.off()

Create a summary plot for top 20 mutated genes

    pdf('SQCC_summary_top20.pdf')
    plotmafSummary(maf = SQCC_top, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
    dev.off()


Create a complete oncoplot for the top 20 most significantly mutated genes

    pdf('SQCC_oncoplot.pdf')
    oncoplot(maf = SQCC, genes = SQCC_genes, rightBarData = SQCC.mutsig, rightBarLims = c(0,4), keepGeneOrder = TRUE, draw_titv = TRUE, keepGeneOrder, cohortSize = 9)
    dev.off()

Create a lollipop plot for all genes in the gene list

    for (gn in SQCC_genes) {
    pdf(paste(gn, ".pdf" , sep = ""))
    lollipopPlot(maf = SQCC, gene = gn, AACol= 'Protein_Change', showMutationRate = TRUE)
    dev.off()
    }

Createa a plot with mutation burden compared against datasets from TCGA

    pdf('SQCC_mutload.pdf')
    SQCC.mutload = tcgaCompare(maf = SQCC, cohortName = 'SQCC')
    dev.off()

Create a graph showing somatic interactions between genes being mutated

    pdf('SQCC_somaticInteractions.pdf')
    somaticInteractions(maf = SQCC, genes = SQCC_genes, pvalue = c(0.05, 0.1))
    dev.off()

Get drug interactions associated with the mutations

    pdf('SQCC_druginteract.pdf')
    dgi = drugInteractions(maf = SQCC_top, fontSize = 0.75)
    dev.off()

Get pathways affected by the mutated genes

    pdf('SQCC_pathway.pdf')
    OncogenicPathways(maf = SQCC)
    dev.off()

Get pathways affected by  the top 20 mutated genes

    pdf('SQCC_top_pathway.pdf')
    OncogenicPathways(maf = SQCC_top)
    dev.off()

>For all the other information about maftools : 
>1. A detailed guide can be found [here](https://www.bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/maftools.html)
>2. A detailed manual with all the tools and its usage options can be found [here](https://rdrr.io/bioc/maftools/man/)






