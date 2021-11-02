# FilterBCF (DNAnexus Platform App)

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://documentation.dnanexus.com/.

Table of Contents:

- [FilterBCF (DNAnexus Platform App)](#filterbcf--dnanexus-platform-app-)
  * [Introduction](#introduction)
    + [Dependencies](#dependencies)
      - [Docker](#docker)
        * [Building The Docker Image](#building-the-docker-image)
      - [Resource Files](#resource-files)
  * [Methodology](#methodology)
    + [1. Split multiallelic variants and normalise all variants](#1-split-multiallelic-variants-and-normalise-all-variants)
    + [2. Perform variant filtering](#2-perform-variant-filtering)
    + [3.  VEP Annotation](#3--vep-annotation)
    + [4. Parsing VEP Consequences](#4-parsing-vep-consequences)
  * [Running on DNANexus](#running-on-dnanexus)
    + [Inputs](#inputs)
    + [Outputs](#outputs)
    + [Command Line Example](#command-line-example)
        * [Batch Running](#batch-running)

## Introduction

This applet filters and annotates a target BCF/VCF according to parameters set by MRC Epidemiology.
This applet makes heavy use of bcftools. Please see the bcftools [manual](https://samtools.github.io/bcftools/bcftools.html)
for more information on how individual commands/tools used in this applet work.

This README makes use of DNANexus file and project naming conventions. Where applicable, an object available on the DNANexus 
platform has a hash ID like:

* file – `file-1234567890ABCDEFGHIJKLMN`
* project – `project-1234567890ABCDEFGHIJKLMN`

Information about files and projects can be queried using the `dx describe` tool native to the DNANexus SDK:

```commandline
dx describe file-1234567890ABCDEFGHIJKLMN
```

**Note:** This README pertains to data included as part of the DNANexus project "MRC - Variant Filtering" (project-G2XK5zjJXk83yZ598Z7BpGPk)
          
### Dependencies

#### Docker

This applet uses [Docker](https://www.docker.com/) to supply dependencies to the underlying AWS instance
launched by DNANexus. All Docker images produced for this and other MRC projects are available as part of Eugene Gardner's 
Dockerhub profile: **[egardner413](https://hub.docker.com/u/egardner413)**. See [below](#building-the-docker-image) for
specifics regarding the image produced for this specific project. The Dockerfile used to build dependencies is available 
as part of this repository at:

`resources/Dockerfile`

This Docker image is built off of the primary 20.04 Ubuntu distribution available via [dockerhub](https://hub.docker.com/layers/ubuntu/library/ubuntu/20.04/images/sha256-644e9b64bee38964c4d39b8f9f241b894c00d71a932b5a20e1e8ee8e06ca0fbd?context=explore).
This image is very light-weight and only provides basic OS installation. Other basic software (e.g. wget, make, and gcc) need
to be installed manually.

In brief, the primary **bioinformatics software** dependencies required by this Dockerfile are:

* [htslib and samtools](http://www.htslib.org/)
* [bcftools](https://samtools.github.io/bcftools/bcftools.html)
* [Variant Effect Predictor](https://www.ensembl.org/info/docs/tools/vep/index.html), which requires the plugins:
    * [LoFTEE](https://github.com/konradjk/loftee)
    * [REVEL](https://github.com/Ensembl/VEP_plugins/blob/release/104/REVEL.pm)
* [plink1.9](https://www.cog-genomics.org/plink2)
* [plink2](https://www.cog-genomics.org/plink/2.0/)

This list is not exhaustive and does not include dependencies of dependencies and software needed
to acquire other resources (e.g. wget). See the referenced Dockerfile for more information. 

**Note:** For most of the above have chosen **not** to hardcode programme version numbers. This means that Docker will install the latest 
version of any software provided as part of the Dockerfile! The exception to this is plink/plink2 as the developers do not provide static links to
the latest version. If running this Dockerfile you will need to change  

This docker image is stored for public use on Dockerhub:

https://hub.docker.com/r/egardner413/mrcepid-filtering

For how this Docker image is used to perform specific tasks as part of this applet, please see below.

##### Building The Docker Image

The docker image was built on an AWS Instance launched via DNANexus. I include here a brief example workflow to enable reproducibility.
Please note that uploading to dockerhub using Eugene Gardner's dockerhub account *will not work* as described in step 5 below.
If you want to store this docker image for yourself, you will need to change the commands accordingly.
All commands are run either locally using the [dx-toolkit](https://documentation.dnanexus.com/downloads) or on a DNANexus AWS instance. 

1. Launch a cloud-workstation:

```commandline
dx run app-cloud_workstation --ssh --instance-type mem1_ssd1_v2_x4 -imax_session_length=2h
```

This command will generate a series of prompts and eventually lead to a query for your DNANexus ssh password. This will then
log you into an AWS instance.

**Note:** One can adjust the amount of time the instance will exist for by changed the `max_session_length`
parameter. 2 hours should be enough to build most Docker instances.

2. Set your permissions to `root` on the instance:

```commandline
dx-su-contrib
```

3. Make and enter a fresh directory for building a Docker image and create a Dockerfile:

```commandline
mkdir dockerbuild
cd dockerbuild
```

After entering this directory, you can either copy-and-paste the text from the provided Dockerfile (using something like 
`vi`), or upload the file using a combination of `dx upload`/`dx download`. Regardless of method, the following commands
require that this file is named `Dockerfile`.

4. Build the Docker image:

```commandline
docker build -t egardner413/mrcepid-filtering:latest .
```

5. Log in to dockerhub and upload the docker image:

```commandline
docker login
# Enter prompted credentials
docker push egardner413/mrcepid-filtering:latest
```

After this last command, the Docker image will be available to pull with the address `egardner413/mrcepid-filtering:latest`.

#### Resource Files

This app also makes use of several resource files that provide various annotations for individual variants. These are 
provided as part of the DNANexus project "MRC - Variant Filtering" (project-G2XK5zjJXk83yZ598Z7BpGPk) in the folder `project_resources/`
with specific directories listed here:

* The hg38 human reference genome. These files are available for free on the DNANexus platform.
    * .fa file: `file-Fx2x270Jx0j17zkb3kbBf6q2`
    * .fa.fai file: `file-Fx2x21QJ06f47gV73kZPjkQQ`
* VEP hg38 cache. Available [here](http://ftp.ensembl.org/pub/release-104/variation/indexed_vep_cache/) – `project_resources/vep_caches`
* Files required for the loftee vep plugin. These files are downloaded from a variety of sources. See [here](https://github.com/konradjk/loftee/tree/grch38) for more details – `project_resources/loftee_files/`
* gnomADv3 MAF files available from [gnomAD](https://gnomad.broadinstitute.org/downloads). This is a simple tabix indexed file for use with `bcftools annotate`  – `project_resources/gnomad_files/`
* REVEL tabix indexed files from [REVEL](https://sites.google.com/site/revelgenomics/downloads) for use with the REVEL vep plugin – `project_resources/revel_files/`

## Methodology

This applet is step 1 (mrc-filterbcf) of the rare variant testing pipeline developed by Eugene Gardner for the UKBiobank RAP at the MRC
Epidemiology Unit:

![](https://github.com/MRCEpid-DNANexus/mrcepid-filterbcf/blob/main/resources/RAPPipeline.png)

As transparency regarding how variant filtering is performed is crucial, I have documented each step that this applet performs below.
Code in this section is meant for example purposes only. For more details, please see the commented source code available at
`src/mrcepid-filterbcf.py` of this repository.

### 1. Split multiallelic variants and normalise all variants

I use `bcftools norm` to first split all multi-allelic variants (e.g. any variant with a REF/ALT VCF field like A T;C into
two separate fields):

```commandline
bcftools norm --threads 4 -Oz -o variants.norm.vcf.gz -m - -f reference.fasta variants.vcf.gz
```

Here, the variants.vcf.gz is the raw input VCF. See [inputs](#inputs) for what this means.

The purpose of the fasta file is to make sure bcftools checks any resulting left-correction against the human reference 
to ensure there are no annotation errors.

### 2. Perform variant filtering

Next, we use bcftools to perform genotype and variant-level filtering. This is performed using `bcftools filter` in a 
two-step process: 

1. Check per-sample genotypes:

```commandline
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
      homozygous ALT genotypes in the 200k data. This will be revisited in the 450k data.
* InDel genotypes are retained if:
  * Depth is ≥ 10 **AND** genotype quality ≥ 20
  
2. Check the proportion of missing genotypes for each variant:

First we calculate per-allele missingness using the bcftools plugin `bcftools +fill-tags`:

```commandline
"bcftools +fill-tags variants.norm.filtered.vcf.gz -Oz -o \
          variants.norm.filtered.tagged.vcf.gz -- -t F_MISSING,AC,AF,AN
```

The above syntax is different than standard due to the use of a plugin, please see the [bcftools](https://samtools.github.io/bcftools/bcftools.html) 
documentation on running plugins for more details. `-t` here is simply telling bcftools to calculate the INFO tags listed
for each allele. Now we perform the filtering. Here you will notice that we are using lower-case `-s` with the flag `FAIL`.
This just tells bcftools to set the vcf FILTER tag to `FAIL` if the allele does not:

* Having missingness ≤ 50%
* Have an allele count > 0

```commandline
# First calculate missingness
bcftools filter -i 'F_MISSING<=0.50 & AC!=0' -s 'FAIL' -Oz -o \
          variants.norm.filtered.tagged.missingness_filtered.vcf.gz \
          variants.norm.filtered.tagged.vcf.gz
```

### 3.  VEP Annotation

This step is relatively straightforward from a "running VEP" standpoint. Please see VEP documentation on [input options](https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html#basic)
for more information on what each flag here means. VEP annotation runs in three separate steps:

1. Conversion of the filtered VCF into a sites vcf to speedup VEP annotation
2. Running of VEP
3. Using `bcftools annotate` to add gnomAD MAF information for all alleles

```commandline
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

```commandline
bcftools +split-vep -df '%CHROM\t%POS\t%REF\t%ALT\t%ID\t%FILTER\t%INFO/AF\t%F_MISSING\t%AN\t%AC\t%MANE_SELECT\t%Feature\t%Gene\t%BIOTYPE\t%CANONICAL\t%SYMBOL\t%Consequence\t%gnomAD_MAF\t%REVEL\t%SIFT\t%PolyPhen\t%LoF\n' \
          -o variants.vep_table.tsv variants.norm.filtered.tagged.missingness_filtered.sites.vep.gnomad.vcf.gz
```

### 4. Parsing VEP Consequences

I have written custom python code to then prioritize **ONE** annotation per variant. This code iterates through the `variants.vep_table.tsv`
generated at the end of the prior step and compares all consequence annotations based on the following, ordered, criteria.
Where the two annotations being compared both fulfill a given criteria, the code moves to the next tier. 

1. Is the consequence for a protein-coding gene?
2. Is the consequence for the [MANE transcript](https://www.ncbi.nlm.nih.gov/refseq/MANE/)
3. Is the consequence for the VEP canonical transcript?
4. Which consequence is more severe? For this, I have generated a priority score that largely matches that [used by VEP](https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences).
   This also gives an annotation 'type' that is used later in the variant testing workflow to group variants with very 
   similar consequences together. Variants with a score closer to 0 are considered more damaging:

| VEP Consequence (CSQ)             | score | type       |
| --------------------------------- | ----- | ---------- |
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

5. Which consequence appears first in the file?

Filtered VCFs are then annotated with information provided by VEP using `bcftools annotate`:

```commandline
bcftools annotate -a variants.vep_table.annote.tsv.gz \ 
          -c "CHROM,POS,REF,ALT,MANE,ENST,ENSG,BIOTYPE,SYMBOL,CSQ,gnomAD_AF,REVEL,SIFT,POLYPHEN,LOFTEE,PARSED_CSQ,MULTI,INDEL,MINOR,MAJOR,MAF,MAC " \
          -h variants.header.txt -Oz -o variants.norm.filtered.tagged.missingness_filtered.annotated.vcf.gz variants.norm.filtered.tagged.missingness_filtered.vcf.gz
```

My hope is the fields that are included in the final annotation are named in a straightforward way. Please open an issue 
on this repository if you feel otherwise.

This final, annotated vcf represents the output of this applet. See below in [outputs](#outputs) for more details.

## Running on DNANexus

### Inputs

|input|description             |
|---- |------------------------|
|vcf  |Input vcf file to filter|

### Outputs

|output                 |description       |
|-----------------------|------------------|
|output_vcf             |  Output VCF with filtered genotypes and sites, annotated with VEP |
|output_vcf_idx         |  .csi index file for `output_vcf`                                 |
|site_stats_pre_filter  |  site_stats.tsv file from pre-filtered VCF with the columns: chrom, pos, filter, %missing, allele count, allele number |
|site_stats_post_filter |  site_stats.tsv file from filtered VCF with the columns: chrom, pos, filter, %missing, allele count, allele number     |
|indv_stats_pre_filter  |  indv_stats.tsv file from pre-filtered VCF generated by `bcftools stats` |
|indv_stats_post_filter |  indv_stats.tsv file from pre-filtered VCF generated by `bcftools stats` |

### Command Line Example

Running this command is fairly straightforward using the DNANexus SDK toolkit. If you are using a project OTHER THAN 
"MRC - Variant Filtering" (project-G2XK5zjJXk83yZ598Z7BpGPk), you will need to compile this applet for your respective DNANexus project:

1. Clone this github repo:

```commandline
git clone https://github.com/MRCEpid-DNANexus/mrcepid-filterbcf.git
```

This will create a folder named mrcepid-filterbcf, you can then:

2. Compile the source code:

```commandline
dx build -f mrcepid-filterbcf
```

the `-f` flag just tells DNANexus to overwrite older versions of the applet within the same project.

**BIG NOTE:** This has not been tested in any project other than "MRC - Variant Filtering". It MAY BREAK due to hardcoded
paths to [Resource Files](#resource-files).

You can then run the following to run this applet. For the input vcf (provided with the flag `-ivcf`) one can use either the file hash OR the full path:

```commandline
# Using file hash
dx run mrcepid-filterbcf --priority low --destination filtered_vcfs/ -ivcf=file-Fz7JXxjJYy0zPf5VFJGGgzBP

# Using full path
dx run mrcepid-filterbcf --priority low --destination filtered_vcfs/ -ivcf="Bulk/Exome sequences/Population level exome OQFE variants, pVCF format/ukb23156_c1_b0_v1.vcf.gz"
```

Some notes here regarding execution:
1. For ease of execution, I prefer using the file hash. This is mostly because DNANexus has put lots of spaces in their filepaths 
   AND it is easier to programmatically access many files at once using hashes as [described below](#batch-running).

2. Outputs are automatically named based on the prefix of the input vcf full path (this is regardless of if you use hash or full path). So 
   the primary VCF output for the above command-line will be `ukb23156_c1_b0_v1.norm.filtered.tagged.missingness_filtered.annotated.vcf.gz`. 
   All outputs will be named using a similar convention.   

3. I have set a sensible (and tested) default for compute resources on DNANexus that is baked into the json used for building the app (at `dxapp.json`) 
   so setting an instance type is unnecessary. This current default is for a mem1_ssd1_v2_x4 instance (4 CPUs, 8 Gb RAM, 100Gb storage). 
   If necessary to adjust compute resources, one can provide a flag like `--instance-type mem1_ssd1_v2_x8`.
   
#### Batch Running

DNANexus have generated a tool to create batch input for running multiple files through this pipeline:

```commandline
dx generate_batch_inputs --path "/Bulk/Exome sequences/Population level exome OQFE variants, pVCF format/" -ivcf='ukb23156_(.*)_v1.vcf.gz$'
```

This will generate files with 500 lines each that have has IDs for all VCF files. One can then modify the above command-line
to run all VCFs in one of these files like so:

```commandline
dx run mrcepid-filterbcf --priority low --destination filtered_vcfs/ --batch-tsv dx_batch.0000.tsv
```