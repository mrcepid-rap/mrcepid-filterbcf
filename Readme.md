# FilterBCF (DNAnexus Platform App)

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://documentation.dnanexus.com/.

### Table of Contents

- [Introduction](#introduction)
  * [Changelog](#changelog)
  * [Background](#background)
  * [Dependencies](#dependencies)
    + [Docker](#docker)
    + [Resource Files](#resource-files)
- [Methodology](#methodology)
  * [Potential Caveats to this Approach](#potential-caveats-to-this-approach)
  * [1. Perform Variant Filtering](#1-perform-variant-filtering)
  * [2. VEP and Additional Annotations](#2-vep-and-additional-annotations)
  * [3. Parsing VEP Consequences](#3-parsing-vep-consequences)
- [Running on DNANexus](#running-on-dnanexus)
  * [Inputs](#inputs)
    + [input_vcfs](#input_vcfs)
    + [additional_annotations](#additional_annotations)
  * [Outputs](#outputs)
    + [output_bcfs](#output_bcfs)
    + [output_veps](#output_veps)
    + [coordinates_file](#coordinates_file)
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

### Changelog

* v2.0.0
  * Applet has gone through a major refactor to support WGS analysis. Do not expect backwards compatability with previous versions.
  * Splitting of multiallelics has been refactored to the [mrcepid-bcfsplitter](https://github.com/mrcepid-rap/mrcepid-bcfsplitter).
  * Removed the per-sample output file as this was computationally expensive and not used in downstream analysis.
  * Added a full testing suite to ensure that the applet is working as expected.
    * Please see `Readme.developer.md` for more information on how to run these tests.
  * The tool no longer specifies default annotations (e.g., REVEL / gnomAD). It now uses a more flexible approach to allow for
    the addition of any annotation that can be added to a VCF file. Please see this README for information on preparing these files.
    * LOFTEE is excluded from this list, as it is run via VEP and requires a specific file structure.
  * It is now possible to set custom values when filtering on GQ (`-igq`).
  * WES and WGS are now supported, please use the `-iwes true` flag to specify that the input is WGS data. Please see the README on why this is required.
  * Variant IDs have been modified to be in `CHR_POS_REF_ALT` format. This is to avoid issues with `:` in the variant ID.
  * The code has been refactored to be more modular and easier to read. Files / methods may be in different locations than before.
  * The default instance type has been changed to an mem1_ssd1_v2_x72. Please update your workflows and cost expectations accordingly.

* v1.0.1
  * Added support for .bcf or .vcf.gz input to this applet

* v1.0.0
  * Initial numbered release. See git changes for history

### Background

The current proposal for variant quality control is to use a missingness-based approach for variant-level filters.
We use this approach as UK Biobank does not provide more fine-tuned parameters that we can use for variant-level
filtering. This applet takes this data as input and performs additional filtering as outlined in the 
[Methology](#methodology) section. WES and WGS are processed slightly differently due to differences in the data; however, 
the general methodolgy is the same.

#### Whole Exome Data

UK Biobank used the “OQFE” [calling approach](https://www.nature.com/articles/s41588-021-00885-0) for the 500k exomes, 
which involves alignment with [BWA mem](http://bio-bwa.sourceforge.net/bwa.shtml) and variant calling with [DeepVariant](https://github.com/google/deepvariant). Variants were restricted to 
±100bps from exome capture regions and then filtered using the following parameters:

1. Hardy-Weinberg Equil. p.value < 1x10<sup>-15</sup>
2. Minimum read coverage depth > 7 for SNVs and > 10 for InDels
3. One sample per variant passed allele balance > 0.15 and > 0.20 for InDels

#### Whole Genome Data

UK Biobank used both a DRAGEN-based and GraphTyper-based variant calling pipeline for the full 500k WGS release. This is 
documented in more detail [here](https://www.medrxiv.org/content/10.1101/2023.12.06.23299426v1). We have chosen the DRAGEN
pipeline for our analysis as it is the preferred datatype, with GraphTyper being deprecated in the future.

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
* [Variant Effect Predictor](https://www.ensembl.org/info/docs/tools/vep/index.html), which requires the plugin:
    * [LoFTEE](https://github.com/konradjk/loftee)

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

## Methodology

This applet is step 2 (mrcepid-filterbcf) of the rare variant testing pipeline:

![](https://github.com/mrcepid-rap/.github/blob/main/images/RAPPipeline.v3.png)

As transparency regarding how variant filtering is performed is crucial, I have documented each step that this applet performs below.
Code in this section is meant for example purposes only. For more details, please see the commented source code available at
`src/mrcepid-filterbcf.py` of this repository.

### Potential Caveats to this Approach

The filtering approach used has several caveats:

1. Potential issues with multiallelic variants – splitting multiallelics is not fullproof. In particular, sites that are multiallelic
   within one individual will not be considered when performing rare variant association testing.
2. X-chromosome calls in males are not adjusted for male hemizygosity and will be ‘1/1’ genotypes (excluding the PAR).
3. Unsure of how well this QC works on the Y-chromosome.
4. The above quality control does not include any sample-level QC. This is assuming that files provided contain reasonable 
   sample-level QC.
5. We prioritise the most severe consequence per variant. This ensures a 1:1 relationship between variant and consequence and thus each variant
6. appears once in the resulting consequences table; however, this could result in rare cases where a more severe 
   consequence is mapped to one of two overlapping genes, and a less severe consequence is ignored. This is a known limitation of this approach.

### 1. Perform Variant Filtering

We 1st use bcftools to perform genotype and variant-level filtering. This is performed using `bcftools filter` in a 
two-step process: 

1. Check per-sample genotypes:

```shell
# gq is provided on the command line at start; by default it is 20
bcftools filter -Ob -o /test/<PREFIX>.filtered.bcf --threads 4 -S . \
          -i '(TYPE="snp" & sSUM(FMT/LAD) >= 7 & (
          (FMT/GT="RR" & FMT/GQ >= gq) |  \
          (FMT/GT="RA" & FMT/GQ >= gq & binom(FMT/LAD) > 0.001) |  \
          (FMT/GT="AA"))) | \
          (TYPE="indel" & sSUM(FMT/DP) >= 10 & FMT/GQ >= gq)' \
          /test/<PREFIX>.bcf
```

This is a complex command-line and filtering is done differently based on genotype (e.g. het/hom) and variant class (e.g. 
SNP/InDel). In brief, the `-i` parameter combined with `-S .` tells bcftools to **i**nclude per **S**ample genotypes 
that pass the described filters. Otherwise, set to missing (./.):

* SNP genotypes are retained if:
  * Depth ≥ 7
  * A homozyogous ref genotype (0/0) with genotype quality ≥ `gq`
  * A heterozygous alt genotype (0/1) with genotype quality ≥ `gq` **AND** [binomial test](https://en.wikipedia.org/wiki/Binomial_test)
  p. value > 1x10<sup>-3</sup>
    * If we did this in R, binomial p. is calculated as if we run the function `binom.test(x, y)`, where x is the number
    of reads supporting the alternate allele (AD[1]), and y is the read depth (DP or AD[0] + AD[1])
  * **ALL** homozygous alt genotypes (1/1)
    * This is due to a [known issue](https://gist.github.com/23andme-jaredo/dd356e97dc55ced4cee1050c892915b8) with homozygous ALT genotypes in the 200k data. I have also checked the new 
      450k data and this still seems to be a problem.
* InDel genotypes are retained if:
  * Depth is ≥ 10 **AND** genotype quality ≥ `gq`
  
**Note**: The process is slightly modified if using WES data (using the `-iwes true` flag at runtime); due to different
processing pipelines, FILTER flags in the VCFs are different. The LAD flag becomes the AD flag in WES data, and the DP flag 
is used rather than summing the LAD flag at each genotype.

2. Check the proportion of missing genotypes and totals for each genotype for each variant

Perform calculations using the bcftools plugin `bcftools +fill-tags`:

```shell
bcftools +fill-tags <PREFIX>.filtered.bcf -Ob \
          -o <PREFIX>.tagged.bcf -- 
          -t F_MISSING,AC,AF,AN,GTM=count(FORMAT/GT == "./."),GT0=count(FORMAT/GT == "0/0"),GT1=count(FORMAT/GT == "0/1"),GT2=count(FORMAT/GT == "1/1")
```

The above syntax is different than standard due to the use of a plugin, please see the [bcftools](https://samtools.github.io/bcftools/bcftools.html) 
documentation on running plugins for more details. `-t` here is telling bcftools to calculate the INFO tags listed
for each allele. The `GTM` tag is the total number of missing genotypes for a given variant. The `GT0`, `GT1`, and `GT2`
tags are the total number of homozygous ref, heterozygous, and homozygous alt genotypes, respectively.

3. Set IDs according to a standard format.

Ensure standard formatting for variant IDs using `bcftools annotate`:

```shell
bcftools annotate -Ob -I '%CHROM\_%POS\_%REF\_%ALT' \
         -Ob --threads 4
         -o <PREFIX>.id_fixed.bcf
         <PREFIX>.tagged.bcf
```

4. Filter variants based on missingness

Perform filtering. We use the flag `-s` with the flag `FAIL`. This tells bcftools to set the vcf FILTER tag to `FAIL` 
if the allele does not:

* Having missingness ≤ 50%
* Have an allele count > 0

```shell
# First calculate missingness
bcftools filter -i 'F_MISSING<=0.50 & AC!=0' -s 'FAIL' -Ob --threads 4 
         -o <PREFIX>.missingness_filtered.bcf \
          <PREFIX>.id_fixed.bcf
```

5. Create a .csi index for the filtered BCF

### 2. VEP and Additional Annotations

This step is relatively straightforward from a "running VEP" standpoint. Please see VEP documentation on [input options](https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html#basic)
for more information on what each flag here means. VEP annotation runs in four separate steps:

1. Conversion of the filtered VCF into a sites vcf to speed up VEP annotation
2. Running of VEP
3. Parsing of VEP consequences
4. Using `annot-tsv` to add additional annotations

1. Convert VEP into a sites file:

```shell
# Generate a sites file from our filtered VCF
bcftools view -G -Oz -o <PREFIX>.missingness_filtered.sites.vcf.gz \
          <PREFIX>.missingness_filtered.bcf
```

2. Run VEP

```shell
# Run VEP
perl -Iensembl-vep/cache/Plugins/loftee/ -Iensembl-vep/cache/Plugins/loftee/maxEntScan/ \
          ensembl-vep/vep --offline --cache --assembly GRCh38 --dir_cache vep_caches/ --everything --allele_num --fork 4 \
          -i <PREFIX>.missingness_filtered.sites.vcf.gz --format vcf --fasta reference.fasta \
          -o <PREFIX>.missingness_filtered.sites.vep.vcf.gz --compress_output bgzip --vcf \
          --dir_plugins ensembl-vep/cache/Plugins/ \
          --plugin LoF,loftee_path:ensembl-vep/cache/Plugins/loftee,human_ancestor_fa:loftee_files/loftee_hg38/human_ancestor.fa.gz,conservation_file:loftee_files/loftee_hg38/loftee.sql,gerp_bigwig:loftee_files/loftee_hg38/gerp_conservation_scores.homo_sapiens.GRCh38.bw
```

3. Converting VEP outputs into TSV 

VEP annotations are then converted into a tab-delimited format using the bcftools plugin `bcftools +split-vep` to generate a 
parseable file for later. Please see the [documentation on split-vep](https://samtools.github.io/bcftools/howtos/plugin.split-vep.html) for more detailed information on how this 
tool works. In brief `-d` causes each **d**uplicate VEP record for a single allele (i.e. those split by `,` and
typically representing transcripts) to be split into a separate row in the tab-delimited file with the information requested with `-f`.

```shell
bcftools +split-vep -df '%CHROM\t%POS\t%REF\t%ALT\t%ID\t%FILTER\t%INFO/AF\t%F_MISSING\t%AN\t%AC\t%MANE_SELECT\t%Feature\t%Gene\t%BIOTYPE\t%CANONICAL\t%SYMBOL\t%Consequence\t%SIFT\t%PolyPhen\t%LoF\tAmino_acids\tProtein_position\tGTM\tGT0\tGT1\tGT2\n' \
          -o variants.vep_table.tsv variants.norm.filtered.tagged.missingness_filtered.sites.vep.gnomad.vcf.gz
```

4. Adding any additional annotations to the VEP tsv file.

The [annot-tsv tool provided with htslib](https://www.htslib.org/doc/annot-tsv.html) is used to add additional annotations:

```shell
annot-tsv 
  -c CHROM,POS,POS
  -f {annotation["annotation_name"]} ' \
  -m <MATCH_STRING> 
  -t /test/<PREFIX>.vep_table.tsv 
  -s /test/<ANNOTATION_NAME>.tsv.gz 
  -o /test/<PREFIX>.<ANNOTATION_NAME>.vep_table.tsv
```

<MATCH_STRING> is a formatted string that describes how to match columns between the two files. For matching purely on coordinates, 
this string is `REF,ALT:REF,ALT` (CHROM, POS must always be matched on with `-c`). For transcript- or symbol-aware
matching this string would be `REF,ALT,SYMBOL:REF,ALT,SYMBOL`. For more information, please see the annot-tsv documentation linked above.

For the format of `<ANNOTATION_NAME>.tsv`, please see the [Additional Annotations Format](#additional-annotations-format) section below.

### 3. Parsing VEP Consequences

I have written custom python code to then prioritize **ONE** annotation per variant. This code iterates through the `<PREFIX>.vep_table.tsv`
generated at the end of the prior step and compares all consequence annotations based on the following, ordered, criteria.
Where the two annotations being compared both fulfill a given criteria, the code moves to the next tier to prioritize the annotation.

1. Is the consequence protein-coding?
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

5. Which consequence appears first in the file?

The final .tsv & filtered .bcf represents the output of this applet. See below in [outputs](#outputs) for more details.

## Running on DNANexus

### Inputs

| input                  | default                                                          | description                                                                                                                                                     |
|------------------------|------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------|
| input_vcfs             | None                                                             | List of files from [mrcepid-bcfsplitter](https://github.com/mrcepid-rap/mrcepid-bcfsplitter) to annotate filter                                                 |
| vep_cache              | None                                                             | Pointer to a VEP cache stored on DNA nexus.                                                                                                                     |
| loftee_libraries       | None                                                             | Pointer to the tar.gz file of loftee databases on DNA nexus. See the [loftee github page](https://github.com/konradjk/loftee/tree/grch38) for more information. |
| coordinates_name       | `coordinates.tsv`                                                | Name of the output storing coordinate information for all output files.                                                                                         |
| human_reference        | `project-Fx2x0fQJ06KfqV7Y3fFZq1jp:file-Fx2x270Jx0j17zkb3kbBf6q2` | GRCh38 human reference genome.                                                                                                                                  |
| human_reference_index  | `project-Fx2x0fQJ06KfqV7Y3fFZq1jp:file-Fx2x21QJ06f47gV73kZPjkQQ` | GRCh38 human reference genome .fai index.                                                                                                                       |
| gq                     | 20                                                               | What genotype quality threshold to use when filtering genotypes                                                                                                 |
| wes                    | false                                                            | Is the data derived from Whole Exome Sequencing?                                                                                                                |
| additional_annotations | None                                                             | A list of additional annotations to add to the vep.tsv.gz output. See below for formatting.                                                                     |

#### input_vcfs

`input_vcfs` is a file list that **MUST** contain DNANexus file hash keys (e.g. like file-1234567890ABCDEFGHIJ). A simple
way to generate a list is with the following bash/perl one-liner; this command will find all filtered vcf files on 
chromosome 7 and print in dna nexus "brief" format which includes a column for file hash (column 7):

```shell
dx ls --brief filtered_vcfs/ukb23148_c7_b*_v1_chunk*.bcf > bcf_list.txt
```

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

#### additional_annotations

Additional annotations are .tsv files that contain the following columns:

```text
#CHROM  POS  REF  ALT  ANNOTATION_NAME
chr1    100  A    T    0.1
chr1    200  C    G    0.2
chr1    300  G    A    0.3
```

An optional column 6 can be added to allow matching on gene. This is useful for annotations that are gene-specific. The
format of this column is:

```text
#CHROM  POS  REF  ALT  ANNOTATION_NAME  SYMBOL
chr1    100  A    T    0.1              GENE1
chr1    200  C    G    0.2              GENE1
chr1    300  G    A    0.3              GENE2
```

This column **must** be named either ENST or SYMBOL, which will match on the Feature and SYMBOL columns output by VEP, respectively.

### Outputs

| output             | description                                                                      |
|--------------------|----------------------------------------------------------------------------------|
| output_bcfs        | Output VCF(s) with filtered genotypes and sites                                  |
| output_bcf_idxs    | .csi index files for `output_vcfs`                                               |
| output_veps        | tabix indexed TSV file with all variants from `outputvcfs` with annotations      |
| output_vep_idxs    | .tbi index files for `output_veps`                                               |
| coordinates_file   | tab-delimited chromosome/start/stop and DNANexus file IDs for each bcf processed |

#### output_bcfs

`output_bcfs` are standard bcf format file. It *does not* contain any annotations (e.g., INFO fields) added by this tool. 
Annotations are stored in tab-delimited output as part of the `output_veps` output.

#### output_veps

Annotations for each variant stored in `output_bcfs`. The format of this table is the following:

| INFO Field               | column  number in .tsv.gz | Field dtype (pandas) | Possible Levels (If Factor)                                                                                       | Short Description                                                                                                                                                                                                              |
|--------------------------|---------------------------|----------------------|-------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| CHROM                    | 1                         | object               | chr1-chr22,chrX, chrY                                                                                             | chromosome of variant                                                                                                                                                                                                          |
| POS                      | 2                         | int64                | NA                                                                                                                | position of variant                                                                                                                                                                                                            |
| REF                      | 3                         | object               | NA                                                                                                                | reference allele (Note that this is not major/minor allele!)                                                                                                                                                                   |
| ALT                      | 4                         | object               | NA                                                                                                                | alternate allele (Note that this is not major/minor allele!)                                                                                                                                                                   |
| ID                       | 5                         | object               | NA                                                                                                                | variant ID. Simple a concatenation like CHROM:POS:REF:ALT. This ID should be used to match on per-marker association test results                                                                                              |
| FILTER                   | 6                         | object               | PASS, FAIL                                                                                                        | Did this variant pass QC?                                                                                                                                                                                                      |
| AF                       | 7                         | float64              | NA                                                                                                                | Allele Frequency of the allele listed in `ALT`                                                                                                                                                                                 |
| F_MISSING                | 8                         | float64              | NA                                                                                                                | Proportion of missing (./.) genotypes                                                                                                                                                                                          |
| AN                       | 9                         | int64                | NA                                                                                                                | Number of possible alleles ([sample size - n.missing] * 2)                                                                                                                                                                     |
| AC                       | 10                        | int64                | NA                                                                                                                | Number of non-reference alleles of the allele listed in `ALT`                                                                                                                                                                  |
| MA                       | 11                        | object               | NA                                                                                                                | Will list the original variant ID for this variant if it was originally part of a multiallelic in the format CHR\|POS\|REF\|ALT_alleles\|NUM where ALT_alleles are comma separated and NUM represents the order of the alleles |
| MANE                     | 12                        | object               | NA                                                                                                                | The [MANE](https://www.ncbi.nlm.nih.gov/refseq/MANE/) transcript for this variant                                                                                                                                              |
| ENST                     | 13                        | object               | NA                                                                                                                | ENSEMBL transcript (e.g. ENST) for this variant                                                                                                                                                                                |
| ENSG                     | 14                        | object               | NA                                                                                                                | ENSEMBL gene (e.g. ENSG) for this variant                                                                                                                                                                                      | 
| BIOTYPE                  | 15                        | object               | All possible values in [this table](https://useast.ensembl.org/info/genome/genebuild/biotypes.html).              | VEP [biotype](https://m.ensembl.org/info/genome/genebuild/biotypes.html) of the ENST.                                                                                                                                          |
| CANONICAL                | 16                        | object               | YES, NO                                                                                                           | Is this transcript the [VEP Canonical](https://useast.ensembl.org/info/genome/genebuild/canonical.html) transcript for this gene?                                                                                              | 
| SYMBOL                   | 17                        | object               | All possible values in [this table](https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html).  | [HGNC Gene](https://www.genenames.org/) symbol                                                                                                                                                                                 |
| CSQ                      | 18                        | object               | NA                                                                                                                | VEP annotated [CSQ](https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html)                                                                                                                                |
| SIFT                     | 19                        | object               | NA                                                                                                                | [SIFT](https://sift.bii.a-star.edu.sg/) score for this variant if CSQ is "missense", else NA                                                                                                                                   |
| POLYPHEN                 | 20                        | object               | HC,LC                                                                                                             | [POLYPHEN](http://genetics.bwh.harvard.edu/pph2/) score for this variant if CSQ is "missense", else NA                                                                                                                         |
| LOFTEE                   | 21                        | object               | NA                                                                                                                | [LOFTEE score](https://github.com/konradjk/loftee) for this variant if PARSED_CSQ is "PTV", else NA. Variants with a "HC" value are high-confidence and should be retained for testing                                         |
| AA                       | 22                        | object               | NA                                                                                                                | The amino acid change if a missense, InDel, or PTV variant                                                                                                                                                                     |
| AApos                    | 23                        | object               | NA                                                                                                                | The position of this variant in the translated protein of this ENST (based on associated ENSP ID)                                                                                                                              |
| GTM                      | 24                        | int64                | NA                                                                                                                | Number of individuals with a missing genotype                                                                                                                                                                                  |
| GT0                      | 25                        | int64                | NA                                                                                                                | Number of individuals with a homozygous REF genotype                                                                                                                                                                           |
| GT1                      | 26                        | int64                | NA                                                                                                                | Number of individuals with a heterozyous genotype                                                                                                                                                                              |
| GT2                      | 27                        | int64                | NA                                                                                                                | Number of individuals with a homozygous ALT genotype                                                                                                                                                                           |
| <Additional Annotations> | 28 to (27 + num. add.)    | varies               | varies                                                                                                            | Type and description varies depending on inputs, see below.                                                                                                                                                                    |
| PARSED_CSQ               | 28 + num. add.            | object               | All possible values in [this table](https://github.com/mrcepid-rap/mrcepid-filterbcf#4-parsing-vep-consequences). | See the section on [parsing consequences](#3-parsing-vep-consequences) for more information                                                                                                                                    |
| IS_MULTIALLELIC          | 29 + num. add.            | bool                 | True, False                                                                                                       | Was this variant original multiallelic? Determined based on presence of at least one ";" in the ID field                                                                                                                       |
| IS_INDEL                 | 30 + num. add.            | bool                 | True, False                                                                                                       | Is this variant an InDel? True if len(REF) != len(ALT)                                                                                                                                                                         | 
| MINOR_ALLELE             | 31 + num. add.            | object               | NA                                                                                                                | The minor allele for this variant. Will be the same as ALT if AF < 0.5                                                                                                                                                         |
| MAJOR_ALLELE             | 32 + num. add.            | object               | NA                                                                                                                | The minor allele for this variant. Will be the same as REF if AF < 0.5                                                                                                                                                         |
| MAF                      | 33 + num. add.            | float64              | NA                                                                                                                | The minor allele frequency for this variant. Will be the same as AF if AF < 0.5                                                                                                                                                |
| MAC                      | 34 + num. add.            | int64                | NA                                                                                                                | The minor allele count for this variant. Will be the same as AC if AF < 0.5                                                                                                                                                    |

#### coordinates_file

`coordinates_file` is a tab-delimited file that contains the chromosome, start, and stop positions of all bcfs processed
via the current applet. The columns are as follows:

| column           | description                           |
|------------------|---------------------------------------|
| chrom            | chromosome                            |
| start            | start position                        |
| end              | end position                          |
| vcf_prefix       | chunk prefix                          |
| output_bcf       | DNANexus ID for the output bcf        |
| output_bcf_idx   | DNANexus ID for the output bcf.csi    |
| output_vep       | DNANexus ID for the output vep.tsv    |
| output_vep_idx   | DNANexus ID for the output vep.tsv.gz |

All coordinates are inclusive; i.e., they are 1-based and include the position listed.

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
   so setting an instance type is unnecessary. This current default is for a mem1_ssd1_v2_x72 instance (72 CPUs, 144 Gb RAM, 1.8Tb storage). 
   **Please note** that this applet is set up for the parallelisation of many files. To run one file, one needs much less 
   memory. If necessary to adjust compute resources, one can provide a flag like `--instance-type mem3_ssd1_v2_x8` to 
   `dx run`.