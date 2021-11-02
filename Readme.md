<!-- dx-header -->
# FilterBCF (DNAnexus Platform App)

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://documentation.dnanexus.com/.

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
Dockerhub profile: **[egardner413](https://hub.docker.com/u/egardner413)**. See below for specifics regarding the image 
produced for this specific project. The Dockerfile used to build dependencies is available as part of this repository at:

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

This applet is STEP 1 of the rare variant testing pipleine developed by Eugene Gardner for the UKBiobank RAP at the MRC
Epidemiology Unit:

![](https://github.com/MRCEpid-DNANexus/mrcepid-filterbcf/tree/main/resources/RAPPipeline.png)



## Running on DNANexus

### Inputs

|input|description|
|---|---|
|vcf|Input vcf file to filter|

### Outputs

|output|description|
|------|-----------|
|output_vcf|           |
|output_vcf_idx|           |


### Command Line Examples
