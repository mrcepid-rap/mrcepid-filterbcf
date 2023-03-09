# FilterBCF (DNAnexus Platform App)

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://documentation.dnanexus.com/.

### Table of Contents

- [Introduction](#introduction)
  * [Background](#background)
  * [Dependencies](#dependencies)
    + [Docker](#docker)
    + [Resource Files](#resource-files)
- [Methodology](#methodology)
  * [Potential Caveats to this Approach](#potential-caveats-to-this-approach)
  * [1. Split multiallelic variants and normalise all variants](#1-split-multiallelic-variants-and-normalise-all-variants)
  * [2. Perform variant filtering](#2-perform-variant-filtering)
  * [3. VEP, gnomAD, and CADD annotation](#3-vep--gnomad--and-cadd-annotation)
  * [4. Parsing VEP consequences](#4-parsing-vep-consequences)
- [Running on DNANexus](#running-on-dnanexus)
  * [Inputs](#inputs)
  * [Outputs](#outputs)
  * [Command line example](#command-line-example)
    
## Introduction

This applet filters and annotates a target BCF/VCF according to parameters set by MRC Epidemiology.
This applet makes heavy use of bcftools. Please see the bcftools [manual](https://samtools.github.io/bcftools/bcftools.html)
for more information on how individual commands/tools used in this applet work.

This README makes use of DNANexus file and project naming conventions. Where applicable, an object available on the DNANexus 
platform has a hash ID like:

* file – `file-1234567890ABCDEFGHIJKLMN`
* project – `project-1234567890ABCDEFGHIJKLMN`

Information about files and projects can be queried using the `dx describe` tool native to the DNANexus SDK:

```shell
dx describe file-1234567890ABCDEFGHIJKLMN
```

**Note:** This README pertains to data included as part of the DNANexus project "MRC - Variant Filtering" (project-G2XK5zjJXk83yZ598Z7BpGPk)

### Background

The current proposal for variant quality control is to use a missingness-based approach for variant-level filters.
We use this approach as UK Biobank does not provide more fine-tuned parameters that we can use for variant-level
filtering. In brief, UK Biobank used the “OQFE” [calling approach](https://www.nature.com/articles/s41588-021-00885-0) for
the 200k exomes, which involves alignment with [BWA mem](http://bio-bwa.sourceforge.net/bwa.shtml) and variant calling
with [DeepVariant](https://github.com/google/deepvariant). Variants were restricted to ±100bps from exome capture regions
and then filtered using the following parameters:

1. Hardy-Weinberg Equil. p.value < 1x10<sup>-15</sup>
2. Minimum read coverage depth > 7 for SNVs and > 10 for InDels
3. One sample per variant passed allele balance > 0.15 and > 0.20 for InDels

This applet takes this data as input and performs additional filtering as outlined in the [Methology](#methodology) section.

### Dependencies

#### Docker

This applet uses [Docker](https://www.docker.com/) to supply dependencies to the underlying AWS instance
launched by DNANexus. The Dockerfile used to build dependencies is available as part of the MRCEpid organisation at:

https://github.com/mrcepid-rap/dockerimages/blob/main/filterbcf.Dockerfile

This Docker image is built off of the primary 20.04 Ubuntu distribution available via [dockerhub](https://hub.docker.com/layers/ubuntu/library/ubuntu/20.04/images/sha256-644e9b64bee38964c4d39b8f9f241b894c00d71a932b5a20e1e8ee8e06ca0fbd?context=explore).
This image is very light-weight and only provides basic OS installation. Other basic software (e.g. wget, make, and gcc) need
to be installed manually. For more details on how to build a Docker image for use on the UKBiobank RAP, please see:

https://github.com/mrcepid-rap#docker-images

In brief, the primary **bioinformatics software** dependencies required by this Applet (and provided in the associated Docker image)
are:

* [htslib and samtools](http://www.htslib.org/)
* [bcftools](https://samtools.github.io/bcftools/bcftools.html)
* [Variant Effect Predictor](https://www.ensembl.org/info/docs/tools/vep/index.html), which requires the plugins:
    * [LoFTEE](https://github.com/konradjk/loftee)
    * [REVEL](https://github.com/Ensembl/VEP_plugins/blob/release/104/REVEL.pm)
* [plink1.9](https://www.cog-genomics.org/plink2)
* [plink2](https://www.cog-genomics.org/plink/2.0/)
* [cadd](https://cadd.gs.washington.edu/)
    * For more details on how to install CADD, please see the [CADD github repo](https://github.com/kircherlab/CADD-scripts/)

This list is not exhaustive and does not include dependencies of dependencies and software needed
to acquire other resources (e.g. wget). See the referenced Dockerfile for more information.

#### Resource Files

This app also makes use of several resource files that provide various annotations for individual variants. To download
these files, we recommend using the [URL Fetcher App](https://ukbiobank.dnanexus.com/app/url_fetcher) on the DNANexus platform (DNANexus login required). Specific
instructions for obtaining these files is listed below each bullet point:

1. The hg38 human reference genome.

  These files are available for free on the DNANexus platform:

    .fa file: `file-Fx2x270Jx0j17zkb3kbBf6q2`
    .fa.fai file: `file-Fx2x21QJ06f47gV73kZPjkQQ`

2. VEP hg38 cache for VEP108. **NOTE:** This applet is designed to work with VEP108, errors may arise if the correct version is not used!

       http://ftp.ensembl.org/pub/release-108/variation/indexed_vep_cache/homo_sapiens_vep_108_GRCh38.tar.gz

3. Files required for the loftee vep plugin. These files are downloaded from a variety of sources. See [the LOFTEE documentation](https://github.com/konradjk/loftee/tree/grch38) for more details. The LOFTEE files require a very specific file-structure to be used by this applet. Please follow the directions below carefully:

    ```{commandline}
    mkdir loftee_hg38/
    cd loftee_hg38/
    wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw
    wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz
    wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.fai
    wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.gzi
    wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/loftee.sql.gz
    gunzip loftee.sql.gz
    cd ../
    tar -czf loftee_hg38.tar.gz loftee_hg38/
    dx upload loftee_hg38.tar.gz
    ```

4. gnomAD***v3*** Non-finish European MAF files for **Hg38** generated from [gnomAD](https://gnomad.broadinstitute.org/downloads). This is a simple tabix indexed file for use with `bcftools annotate`. For instructions on generating this file, see below:

    ```{commandline}
    # You will need to run the following for all chromosomes available, but an example for chromosome 1 is shown below.
    bcftools query -f '%CHROM\t%POS\t$REF\t%ALT\t%FILTER\t%AF_nfe\n' https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr1.vcf.bgz > chr1.gnomad.tsv
    
    # Once all 23 chromosomes have completed, run the following:
    cat chr*.gnomad.tsv | sort -k 1,1 -k 2,2n > gnomad.tsv
    bgzip gnomad.tsv
    tabix -s 1 -b 2 -e 2 gnomad.tsv.gz
    dx upload gnomad.tsv.gz
    ```

5. REVEL tabix indexed files from [REVEL](https://sites.google.com/site/revelgenomics/downloads) for use with the REVEL vep plugin:

    ```{commandline}
    # Download the original REVEL file:
    wget https://zenodo.org/record/7072866/files/revel-v1.3_all_chromosomes.zip
   
    # Process according to instructions here: https://github.com/Ensembl/VEP_plugins/blob/release/109/REVEL.pm
    unzip revel-v1.3_all_chromosomes.zip
    cat revel_with_transcript_ids | tr "," "\t" > tabbed_revel.tsv
    sed '1s/.*/#&/' tabbed_revel.tsv > new_tabbed_revel.tsv
    bgzip new_tabbed_revel.tsv
    zcat new_tabbed_revel.tsv.gz | head -n1 > h
    zgrep -h -v ^#chr new_tabbed_revel.tsv.gz | awk '$3 != "." ' | sort -k1,1 -k3,3n - | cat h - | bgzip -c > new_tabbed_revel_grch38.tsv.gz
    tabix -f -s 1 -b 3 -e 3 new_tabbed_revel_grch38.tsv.gz
    dx upload new_tabbed_revel_grch38.tsv.gz
    ```

6. CADD resource files downloaded from the [CADD Downloads website](https://cadd.gs.washington.edu/download). **NOTE**: These files are VERY large:

       CADD Annotations: https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/annotationsGRCh38_v1.6.tar.gz
       Precomputed SNVs: https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz
       Precomputed SNVs index: https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz.tbi
       Precomputed InDels: https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.indel.tsv.gz
       Precomputed InDels index: https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.indel.tsv.gz.tbi

## Methodology

This applet is step 2 (mrc-filterbcf) of the rare variant testing pipeline developed by Eugene Gardner for the UKBiobank RAP at the MRC
Epidemiology Unit:

![](https://github.com/mrcepid-rap/.github/blob/main/images/RAPPipeline.v3.png)

As transparency regarding how variant filtering is performed is crucial, I have documented each step that this applet performs below.
Code in this section is meant for example purposes only. For more details, please see the commented source code available at
`src/mrcepid-filterbcf.py` of this repository.

### Potential Caveats to this Approach

The filtering approach used has several caveats:

1. Potential issues with multi-allelic variants – splitting multi-allelics is not fullproof. In particular, sites that are multiallelic
   within one individual will not be considered when performing rare variant association testing.
2. X-chromosome calls in males are not adjusted for male hemizygosity and will typically be ‘1/1’ genotypes (excluding the PAR).
3. Unsure of how well this QC works on the Y-chromosome.
4. The above quality control does not include any sample-level variant QC. This is assuming that UKB (a.k.a. Regeneron)
   does reasonable sample-level QC and excludes any samples that fail various filters.

### 1. Split multiallelic variants and normalise all variants

I use `bcftools norm` to first split all multi-allelic variants (e.g. any variant with a REF/ALT VCF field like A T;C into
two separate fields):

```shell
bcftools norm --threads 4 -Oz -o variants.norm.vcf.gz -m - -f reference.fasta variants.vcf.gz
```

Here, the variants.vcf.gz is the raw input VCF. See [inputs](#inputs) for what this means.

The purpose of the fasta file is to make sure bcftools checks any resulting left-correction against the human reference 
to ensure there are no annotation errors.

### 2. Perform variant filtering

Next, we use bcftools to perform genotype and variant-level filtering. This is performed using `bcftools filter` in a 
two-step process: 

1. Check per-sample genotypes:

```shell
bcftools filter -Oz -o /test/variants.norm.filtered.vcf.gz -S . \
          -i '(TYPE="snp" & FMT/DP >= 7 & ((FMT/GT="RR" & FMT/GQ >= 20) |  \
          (FMT/GT="RA" & FMT/GQ >= 20 & binom(FMT/AD) > 0.001) |  \
          (FMT/GT="AA"))) | \
          (TYPE="indel" & FMT/DP >= 10 & FMT/GQ >= 20)' \
          /test/variants.norm.vcf.gz
```

This is a complex command-line and filtering is done differently based on genotype (e.g. het/hom) and variant class (e.g. 
SNP/InDel). Thus, I am going to break it down here. In brief, the `-i` parameter combined with `-S .` tells bcftools to 
**i**nclude per **S**ample genotypes that pass the described filters. Otherwise, set to missing (./.):

* SNP genotypes are retained if:
  * Depth ≥ 7
  * A homozyogous ref genotype (0/0) with genotype quality ≥ 20
  * A heterozygous alt genotype (0/1) with genotype quality ≥ 20 **AND** [binomial test](https://en.wikipedia.org/wiki/Binomial_test)
  p. value > 1x10<sup>-3</sup>
    * If we did this in R, binomial p. is calculated as if we run the function `binom.test(x, y)`, where x is the number
    of reads supporting the alternate allele (AD[1]), and y is the read depth (DP or AD[0] + AD[1])
  * **ALL** homozygous alt genotypes (1/1)
    * This is due to a [known issue](https://gist.github.com/23andme-jaredo/dd356e97dc55ced4cee1050c892915b8) with 
      homozygous ALT genotypes in the 200k data. I have also checked the new 450k data and this still seems to be a problem.
* InDel genotypes are retained if:
  * Depth is ≥ 10 **AND** genotype quality ≥ 20
  
2. Check the proportion of missing genotypes for each variant:

First we calculate per-allele missingness using the bcftools plugin `bcftools +fill-tags`:

```shell
bcftools +fill-tags variants.norm.filtered.vcf.gz -Oz -o \
          variants.norm.filtered.tagged.vcf.gz -- -t F_MISSING,AC,AF,AN
```

The above syntax is different than standard due to the use of a plugin, please see the [bcftools](https://samtools.github.io/bcftools/bcftools.html) 
documentation on running plugins for more details. `-t` here is simply telling bcftools to calculate the INFO tags listed
for each allele. Now we perform the filtering. Here you will notice that we are using lower-case `-s` with the flag `FAIL`.
This just tells bcftools to set the vcf FILTER tag to `FAIL` if the allele does not:

* Having missingness ≤ 50%
* Have an allele count > 0

```shell
# First calculate missingness
bcftools filter -i 'F_MISSING<=0.50 & AC!=0' -s 'FAIL' -Oz -o \
          variants.norm.filtered.tagged.missingness_filtered.vcf.gz \
          variants.norm.filtered.tagged.vcf.gz
```

### 3. VEP, gnomAD, and CADD annotation

This step is relatively straightforward from a "running VEP" standpoint. Please see VEP documentation on [input options](https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html#basic)
for more information on what each flag here means. VEP annotation runs in three separate steps:

1. Conversion of the filtered VCF into a sites vcf to speedup VEP annotation
2. Running of VEP
3. Running of CADD
4. Using `bcftools annotate` to add gnomAD MAF and CADD score information for all alleles

```shell
# Generate a sites file from our filtered VCF
bcftools view -G -Oz -o variants.norm.filtered.tagged.missingness_filtered.sites.vcf.gz \
          variants.norm.filtered.tagged.missingness_filtered.vcf.gz

# Run VEP
perl -Iensembl-vep/cache/Plugins/loftee/ -Iensembl-vep/cache/Plugins/loftee/maxEntScan/ \
          ensembl-vep/vep --offline --cache --assembly GRCh38 --dir_cache vep_caches/ --everything --allele_num \
          -i variants.norm.filtered.tagged.missingness_filtered.sites.vcf.gz --format vcf --fasta reference.fasta \
          -o variants.norm.filtered.tagged.missingness_filtered.sites.vep.vcf.gz --compress_output bgzip --vcf \
          --dir_plugins ensembl-vep/cache/Plugins/ \
          --plugin LoF,loftee_path:ensembl-vep/cache/Plugins/loftee,human_ancestor_fa:loftee_files/loftee_hg38/human_ancestor.fa.gz,conservation_file:loftee_files/loftee_hg38/loftee.sql,gerp_bigwig:loftee_files/loftee_hg38/gerp_conservation_scores.homo_sapiens.GRCh38.bw \
          --plugin REVEL,revel_files/new_tabbed_revel_grch38.tsv.gz

# Run CADD
CADD-scripts/CADD.sh -g GRCh38 -c 2 -o variants.cadd.tsv.gz variants.vcf
          
# gnomAD MAF information per-allele is added:
bcftools annotate -a gnomad_files/gnomad.tsv.gz -c CHROM,POS,REF,ALT,-,gnomAD_MAF \
          -h gnomad_files/gnomad.header.txt -Oz -o variants.norm.filtered.tagged.missingness_filtered.sites.vep.gnomad.vcf.gz \
          variants.norm.filtered.tagged.missingness_filtered.sites.vep.vcf.gz
```

VEP annotations are then parsed into a tab-delimited format using the bcftools plugin `bcftools +split-vep` to generate an easy-to-parse
file for later. Please see the [documentation on split-vep](https://samtools.github.io/bcftools/howtos/plugin.split-vep.html) for more detailed
information on how this tool works. In brief `-d` causes each **d**uplicate VEP record for a single allele (i.e. those split 
by `,` and typically representing transcripts) to be split into a separate row in the tab-delimited file with the
information requested with `-f`.

```shell
bcftools +split-vep -df '%CHROM\t%POS\t%REF\t%ALT\t%ID\t%FILTER\t%INFO/AF\t%F_MISSING\t%AN\t%AC\t%MANE_SELECT\t%Feature\t%Gene\t%BIOTYPE\t%CANONICAL\t%SYMBOL\t%Consequence\t%gnomAD_MAF\t%REVEL\t%SIFT\t%PolyPhen\t%LoF\n' \
          -o variants.vep_table.tsv variants.norm.filtered.tagged.missingness_filtered.sites.vep.gnomad.vcf.gz
```

### 4. Parsing VEP consequences

I have written custom python code to then prioritize **ONE** annotation per variant. This code iterates through the `variants.vep_table.tsv`
generated at the end of the prior step and compares all consequence annotations based on the following, ordered, criteria.
Where the two annotations being compared both fulfill a given criteria, the code moves to the next tier. 

1. Is the consequence for a protein-coding *gene*?
2. Is the consequence for the [MANE transcript](https://www.ncbi.nlm.nih.gov/refseq/MANE/)
3. Is the consequence for the VEP canonical transcript?
4. Which consequence is more severe? For this, I have generated a priority score that largely matches that [used by VEP](https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences).
   This also gives an annotation 'type' that is used later in the variant testing workflow to group variants with very 
   similar consequences together. Variants with a score closer to 0 are considered more damaging:

| VEP Consequence (CSQ)             | score | type       |
|-----------------------------------|-------|------------|
| stop_gained                       | 1     | PTV        |
| frameshift_variant                | 2     | PTV        |
| splice_acceptor_variant           | 3     | PTV        |
| splice_donor_variant              | 4     | PTV        |
| stop_lost                         | 5     | STOP_LOST  |
| start_lost                        | 6     | START_LOST |
| inframe_insertion                 | 7     | INFRAME    |
| inframe_deletion                  | 8     | INFRAME    |
| missense_variant                  | 9     | MISSENSE   |
| protein_altering_variant          | 10    | INFRAME    |
| splice_region_variant             | 11    | NONCODING  |
| incomplete_terminal_codon_variant | 12    | INFRAME    |
| start_retained_variant            | 13    | SYN        |
| stop_retained_variant             | 14    | SYN        |
| synonymous_variant                | 15    | SYN        |
| 5_prime_UTR_variant               | 16    | UTR        |
| 3_prime_UTR_variant               | 17    | UTR        |
| intron_variant                    | 18    | INTRONIC   |
| intergenic_variant                | 19    | INTERGENIC |
| upstream_gene_variant             | 20    | INTERGENIC |
| downstream_gene_variant           | 21    | INTERGENIC |
| no_score                          | 22    | ERROR      |

The end-result of this iterative process is each variant in the VCF file gets a collection of INFO fields that contain annotations.
For information on these annotations, please see the README for the next applet in this pipeline [mrcepid-annotatecadd](https://github.com/mrcepid-rap/mrcepid-annotatecadd).

5. Which consequence appears first in the file?

Filtered VCFs are then annotated with information provided by VEP using `bcftools annotate`:

```shell
bcftools annotate -a variants.vep_table.annote.tsv.gz \ 
          -c "CHROM,POS,REF,ALT,MANE,ENST,ENSG,BIOTYPE,SYMBOL,CSQ,gnomAD_AF,REVEL,SIFT,POLYPHEN,LOFTEE,PARSED_CSQ,MULTI,INDEL,MINOR,MAJOR,MAF,MAC " \
          -h variants.header.txt -Oz -o variants.norm.filtered.tagged.missingness_filtered.annotated.vcf.gz variants.norm.filtered.tagged.missingness_filtered.vcf.gz
```

My hope is the fields that are included in the final annotation are named in a straightforward way. Please open an issue 
on this repository if you feel otherwise.

This final, annotated vcf represents the output of this applet. See below in [outputs](#outputs) for more details.

## Running on DNANexus

### Inputs

| input      | description                                                                                                     |
|------------|-----------------------------------------------------------------------------------------------------------------|
| input_vcfs | List of files from [mrcepid-bcfsplitter](https://github.com/mrcepid-rap/mrcepid-bcfsplitter) to annotate filter |

`input_vcfs` is a file list that **MUST** contain DNANexus file hash keys (e.g. like file-1234567890ABCDEFGHIJ). A simple
way to generate such a list is with the following bash/perl one-liner:

```shell
dx ls -l filtered_vcfs/ukb23148_c7_b*_v1_chunk*.bcf | perl -ane 'chomp $_; if ($F[6] =~ /^\((\S+)\)$/) {print "$1\n";}' > bcf_list.txt
```

This command will:

1. Find all filtered vcf files on chromosome 7 and print in dna nexus "long" format which includes a column for file hash (column 7)
2. Extract the file hash using a perl one-liner and print one file hash per line

The final input file will look something like:

```text
file-1234567890ABCDEFGHIJ
file-2345678901ABCDEFGHIJ
file-3456789012ABCDEFGHIJ
file-4567890123ABCDEFGHIJ
```

This file then needs to be uploaded to the DNANexus platform, so it can be provided as input:

```shell
dx upload bcf_list.txt
```

There are also several command-line inputs that should not need to be changed if running from within application 9905. These
mostly have to do with the underlying inputs to models that are generated by other tools in this pipeline. We have set
sensible defaults for these files and only change them if running from a different set of filtered data.

| input                   | description                                                                                                                                                                                                                                                      | 
|-------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| human_reference         | GRCh38 human reference genome. This file is available as part of standard DNANexus resources in project `project-Fx2x0fQJ06KfqV7Y3fFZq1jp`.                                                                                                                      |
| human_reference_index   | GRCh38 human reference genome .fai index. This file is available as part of standard DNANexus resources in project `project-Fx2x0fQJ06KfqV7Y3fFZq1jp`.                                                                                                           |
| vep_cache               | Pointer to the VEP cache stored on DNA nexus. Current default is v108                                                                                                                                                                                            |
| loftee_libraries        | Pointer to the tar.gz file of loftee databases on DNA nexus. See the [loftee github page](https://github.com/konradjk/loftee/tree/grch38) for more information.                                                                                                  |
| gnomad_maf_db           | Pointer to a precompiled .tsv format file of gnomAD AFs. Calculated from gnomAD v2 vcf files for Non-Finish European individuals. A corresponding `.tbi` index MUST be located at the same location as this file!                                                |
| revel_db                | Precompiled REVEL score database. See the [REVEL](https://sites.google.com/site/revelgenomics/) website for more information. A corresponding `.tbi` index MUST be located at the same location as this file!                                                    |
| cadd_annotations        | Pointer to annotations used by CADD to compute variant deleteriousness scores. See the [CADD downloads website](https://cadd.bihealth.org/download) for more information.                                                                                        |
| precomputed_cadd_snvs   | `tsv.gz` file of all possible SNVs in the human genome with pre-computed CADD scores. See the [CADD downloads website](https://cadd.bihealth.org/download) for more information. A corresponding `.tbi` index MUST be located at the same location as this file! |
| precomputed_cadd_indels | `tsv.gz` file of all gnomAD3.0 InDels with pre-computed CADD scores. See the [CADD downloads website](https://cadd.bihealth.org/download) for more information. A corresponding `.tbi` index MUST be located at the same location as this file!                  |

### Outputs

| output             | description                                                                      |
|--------------------|----------------------------------------------------------------------------------|
| output_bcfs        | Output VCF(s) with filtered genotypes and sites, annotated with CADD             |
| output_bcf_idxs    | .csi index files for `output_vcfs`                                               |
| output_veps        | tabix indexed TSV file with all variants from `outputvcfs` with annotations      |
| output_vep_idxs    | .tbi index files for `output_veps`                                               |
| output_per_samples | .tsv file of variants per-individual to collate number of variants per person    |
| coordinates_file   | tab-delimited chromosome/start/stop and DNANexus file IDs for each bcf processed |

`output_vcfs` is a standard vcf.gz format file with INFO fields derived from VEP annotations. These fields are identical 
to those in the paired .tsv.gz file provided with the `output_veps` output. These are the following (note that `varID` is
added during the following `makebgen` step):

| INFO Field  | row number in .tsv.gz | Field dtype (pandas) | Possible Levels (If Factor)                                                                                       | Short Description                                                                                                                                                                                           |
|-------------|-----------------------|----------------------|-------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| CHROM       | 1                     | object               | chr1-chr22,chrX, chrY                                                                                             | chromosome of variant                                                                                                                                                                                       |
| POS         | 2                     | int64                | NA                                                                                                                | position of variant                                                                                                                                                                                         |
| REF         | 3                     | object               | NA                                                                                                                | reference allele (Note that this is not major/minor allele!)                                                                                                                                                |
| ALT         | 4                     | object               | NA                                                                                                                | alternate allele (Note that this is not major/minor allele!)                                                                                                                                                |
| varID       | 5                     | object               | NA                                                                                                                | variant ID. Simple a concatenation like CHROM:POS:REF:ALT. This ID should be used to match on per-marker association test results.                                                                          |
| ogVarID     | 6                     | object               | NA                                                                                                                | original variant ID in UKBiobank VCF files. **Note:** This ID SHOULD NOT be used to match on association test results. It is included for posterity and to allow matching on original UK Biobank VCF files. |
| FILTER      | 7                     | object               | PASS, FAIL                                                                                                        | Did this variant pass QC?                                                                                                                                                                                   |
| AF          | 8                     | float64              | NA                                                                                                                | Allele Frequency of the allele listed in `ALT`                                                                                                                                                              |
| F_MISSING   | 9                     | float64              | NA                                                                                                                | Proportion of missing (./.) genotypes                                                                                                                                                                       |
| AN          | 10                    | int64                | NA                                                                                                                | Number of possible alleles ([sample size - n.missing] * 2)                                                                                                                                                  |
| AC          | 11                    | int64                | NA                                                                                                                | Number of non-reference alleles of the allele listed in `ALT`                                                                                                                                               |
| MANE        | 12                    | object               | NA                                                                                                                | The [MANE](https://www.ncbi.nlm.nih.gov/refseq/MANE/) transcript for this variant                                                                                                                           |
| ENST        | 13                    | object               | NA                                                                                                                | ENSEMBL transcript (e.g. ENST) for this variant                                                                                                                                                             |
| ENSG        | 14                    | object               | NA                                                                                                                | ENSEMBL gene (e.g. ENSG) for this variant                                                                                                                                                                   | 
| BIOTYPE     | 15                    | object               | protein_coding                                                                                                    | VEP [biotype](https://m.ensembl.org/info/genome/genebuild/biotypes.html) of the ENST.                                                                                                                       |
| SYMBOL      | 16                    | object               | NA                                                                                                                | [HGNC Gene](https://www.genenames.org/) symbol                                                                                                                                                              | 
| CSQ         | 17                    | object               | All possible values in [this table](https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html).  | VEP annotated [CSQ](https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html)                                                                                                             |
| gnomAD_AF   | 18                    | float64              | NA                                                                                                                | Non-Finnish European allele frequency of this variant in gnomAD exomes. 0 if not in gnomAD                                                                                                                  |
| CADD        | 19                    | float64              | NA                                                                                                                | [CADDv1.6](https://cadd.gs.washington.edu/) phred score                                                                                                                                                     |
| REVEL       | 20                    | float64              | NA                                                                                                                | [REVEL](https://sites.google.com/site/revelgenomics/) score for this variant if CSQ is "missense", else nan                                                                                                 |
| SIFT        | 21                    | object               | NA                                                                                                                | [SIFT](https://sift.bii.a-star.edu.sg/) score for this variant if CSQ is "missense", else NA                                                                                                                |
| POLYPHEN    | 22                    | object               | NA                                                                                                                | [POLYPHEN](http://genetics.bwh.harvard.edu/pph2/) score for this variant if CSQ is "missense", else NA                                                                                                      |
| LOFTEE      | 23                    | object               | HC,LC                                                                                                             | [LOFTEE score](https://github.com/konradjk/loftee) for this variant if PARSED_CSQ is "PTV", else NA. Variants with a "HC" value are high-confidence and should be retained for testing                      |
| AA          | 24                    | object               | NA                                                                                                                | The amino acid change if a missense, InDel, or PTV variant                                                                                                                                                  |
| AApos       | 25                    | object               | NA                                                                                                                | The position of this variant in the translated protein of this ENST (based on associated ENSP ID)                                                                                                           |
| PARSED_CSQ  | 26                    | object               | All possible values in [this table](https://github.com/mrcepid-rap/mrcepid-filterbcf#4-parsing-vep-consequences). | Eugene-determined consequence. See the README for [mrcepid-filterbcf](https://github.com/mrcepid-rap/mrcepid-filterbcf) for more information                                                                |
| MULTI       | 27                    | bool                 | True, False                                                                                                       | Was this variant original multiallelic? Determined based on presence of at least one ";" in the ID field                                                                                                    |
| INDEL       | 28                    | bool                 | True, False                                                                                                       | Is this variant an InDel? True if len(REF) != len(ALT)                                                                                                                                                      | 
| MINOR       | 29                    | object               | NA                                                                                                                | The minor allele for this variant. Will be the same as ALT if AF < 0.5                                                                                                                                      |
| MAJOR       | 30                    | object               | NA                                                                                                                | The minor allele for this variant. Will be the same as REF if AF < 0.5                                                                                                                                      |
| MAF         | 31                    | float64              | NA                                                                                                                | The minor allele frequency for this variant. Will be the same as AF if AF < 0.5                                                                                                                             |
| MAC         | 32                    | int64                | NA                                                                                                                | The minor allele count for this variant. Will be the same as AC if AF < 0.5                                                                                                                                 |

### Command line example

If this is your first time running this applet within a project other than "MRC - Variant Filtering", please see our 
organisational documentation on how to download and build this app on the DNANexus Research Access Platform:

https://github.com/mrcepid-rap

Running this command is fairly straightforward using the DNANexus SDK toolkit. For the input vcf (provided with the flag 
`-iinputvcf`) one can use either the file hash OR the full path:

```shell
# Using file hash
dx run mrcepid-filterbcf --priority low --destination filtered_vcfs/ -iinputvcfs=file-Fz7JXxjJYy0zPf5VFJGGgzBP

# Using full path
dx run mrcepid-filterbcf --priority low --destination filtered_vcfs/ -iinputvcfs=bcf_list.txt
```

Brief I/O information can also be retrieved on the command line:

```shell
dx run mrcepid-filterbcf --help
```

Some notes here regarding execution:
1. For ease of execution, I prefer using the file hash. This is mostly because DNANexus has put lots of spaces in their 
   filepaths, AND it is easier to programmatically access many files at once using hashes.

2. Outputs are automatically named based on the prefix of the input vcf full path (this is regardless of if you use hash or full path). So 
   the primary VCF output for the above command-line will be `ukb23156_c1_b0_v1.filtered_annotated.bcf`. All outputs 
   will be named using a similar convention.

3. I have set a sensible (and tested) default for compute resources on DNANexus that is baked into the json used for building the app (at `dxapp.json`) 
   so setting an instance type is unnecessary. This current default is for a mem2_ssd1_v2_x32 instance (32 CPUs, 128 Gb RAM, 1200Gb storage). This
   instance prioritises more RAM over other types of instances, which is required for the [normalisation step](#1-split-multiallelic-variants-and-normalise-all-variants)
   outlined above. **Please note** that this applet is set up for the parallelisation of many files. To run one file, one needs much less 
   memory. If necessary to adjust compute resources, one can provide a flag like `--instance-type mem3_ssd1_v2_x8` to 
   `dx run`.