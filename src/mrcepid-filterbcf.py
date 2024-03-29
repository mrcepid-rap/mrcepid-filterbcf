#!/usr/bin/env python
# mrcepid-filterbcf 1.0.0
# Generated by dx-app-wizard.
#
# Author: Eugene Gardner (eugene.gardner at mrc.epid.cam.ac.uk)
#
# DNAnexus Python Bindings (dxpy) documentation:
#   http://autodoc.dnanexus.com/bindings/python/current/
import csv
import sys
import dxpy

from pathlib import Path
from time import sleep
from os.path import exists
from typing import TypedDict, Tuple

from general_utilities.association_resources import download_dxfile_by_name
from general_utilities.job_management.thread_utility import ThreadUtility

# We have to do this to get modules to run properly on DNANexus while still enabling easy editing in PyCharm
sys.path.append('/')
sys.path.append('/filterbcf/')

from filterbcf.ingest_data import IngestData
from filterbcf.vcf_filter.vcf_filter import VCFFilter
from filterbcf.vcf_annotate.vcf_annotate import VCFAnnotate


class ProcessedReturn(TypedDict):
    chrom: str
    start: int
    end: int
    vcf_prefix: str
    output_bcf: dxpy.DXFile
    output_bcf_idx: dxpy.DXFile
    output_vep: dxpy.DXFile
    output_vep_idx: dxpy.DXFile
    output_per_sample: dxpy.DXFile


def identify_vcf_format(vcf_name: str) -> Tuple[str, str]:
    """Decide which variant data type we are using

    Identifies .bcf or .vcf.gz binary vcf formats and reports the prefix and suffix. Also throws an error if the suffix
    is not correct. Thus

    test.vcf.gz or test.bcf

    becomes

    test, .vcf.gz or test, .bcf

    :param vcf_name: A vcf file in str format
    :return: A tuple of vcf_prefix, vcf_suffix
    """

    if vcf_name.endswith('.vcf.gz'):
        vcf_prefix = vcf_name[:-7]
        vcf_suffix = '.vcf.gz'
    elif vcf_name.endswith('.bcf'):
        vcf_prefix = vcf_name[:-4]
        vcf_suffix = '.bcf'
    else:
        raise ValueError('VCF does not end with .vcf.gz or .bcf and may indicate formatting issues!')

    return vcf_prefix, vcf_suffix


# This is a method that will execute all steps necessary to process one VCF file
# It is the primary unit that is executed by individual threads from the 'main()' method
def process_vcf(vcf: str) -> ProcessedReturn:

    # Create a DXFile instance of the given file:
    vcf = dxpy.DXFile(vcf)

    # Get the suffix on the vcf file
    vcfname = vcf.describe()['name']
    vcfprefix, vcfsuffix = identify_vcf_format(vcfname)

    # Download the VCF file chunk to the instance
    download_dxfile_by_name(vcf, project_id=dxpy.PROJECT_CONTEXT_ID)

    # 1. Do normalisation and filtering
    print(f'Filtering bcf: {vcfname}')
    VCFFilter(vcfprefix, vcfsuffix)

    # We need to pause here in each thread to make sure that CADD and VEP files have downloaded in separate threads...
    # We know that when the original .tar.gz files are gone, that it is safe as deleting those files is the final step
    # of the download process.
    vep_tar = 'vep_caches/vep_cache.tar.gz'
    cadd_tar = 'cadd_files/annotationsGRCh38_v1.6.tar.gz'
    precomputed_index = 'cadd_precomputed/gnomad.genomes.r3.0.indel.tsv.gz.tbi'
    while (exists(vep_tar) or exists(cadd_tar)) and exists(precomputed_index) is False:
        sleep(5)

    # 2. Do annotation
    print(f'Annotating bcf: {vcfname}')
    vcf_annotater = VCFAnnotate(vcfprefix)

    print(f'Finished bcf: {vcfname}')

    return {'chrom': vcf_annotater.chunk_chrom,
            'start': vcf_annotater.chunk_start,
            'end': vcf_annotater.chunk_end,
            'vcf_prefix': vcfprefix,
            'output_bcf': vcf_annotater.finalbcf,
            'output_bcf_idx': vcf_annotater.finalbcf_index,
            'output_vep': vcf_annotater.finalvep,
            'output_vep_idx': vcf_annotater.finalvep_index,
            'output_per_sample': vcf_annotater.output_per_sample}


@dxpy.entry_point('main')
def main(input_vcfs, coordinates_name, human_reference, human_reference_index, vep_cache, loftee_libraries,
         gnomad_maf_db, revel_db,
         cadd_annotations, precomputed_cadd_snvs, precomputed_cadd_indels):

    # Separate function to acquire necessary resource files
    ingested_data = IngestData(input_vcfs, human_reference, human_reference_index, vep_cache, loftee_libraries,
                               gnomad_maf_db, revel_db,
                               cadd_annotations, precomputed_cadd_snvs, precomputed_cadd_indels)

    # Now build a thread worker that contains as many threads, divided by 2 that have been requested since each bcftools
    # 1 thread for monitoring threads
    # 2 threads for downloading (1 each for CADD and VEP)
    # 2 threads for each BCF
    thread_utility = ThreadUtility(thread_factor=2, error_message='A bcffiltering thread failed', incrementor=5)

    # And launch the requested threads
    for input_vcf in ingested_data.input_vcfs:
        thread_utility.launch_job(process_vcf,
                                  vcf=input_vcf)
    print("All threads submitted...")

    # And add the resulting futures to relevant output arrays / file
    output_bcfs = []
    output_bcf_idxs = []
    output_veps = []
    output_vep_idxs = []
    output_per_samples = []

    # This file contains information about the 'chunks' that have been processed. It DOES NOT have a header to
    # make it easier to concatenate coordinate files from multiple runs. The columns are as follows:
    # [0] #chrom
    # [1] start
    # [2] end
    # [3] chunk_prefix
    # [4] bcf_dxpy
    # [5] vep_dxpy

    with Path(coordinates_name).open('w') as coordinate_writer:

        coordinate_csv = csv.DictWriter(coordinate_writer,
                                        fieldnames=['chrom', 'start', 'end', 'vcf_prefix', 'output_bcf', 'output_vep'],
                                        delimiter="\t")
        # And gather the resulting futures
        for result in thread_utility:
            output_bcfs.append(result['output_bcf'])
            output_bcf_idxs.append(result['output_bcf_idx'])
            output_veps.append(result['output_vep'])
            output_vep_idxs.append(result['output_vep_idx'])
            output_per_samples.append(result['output_per_sample'])
            writer_dict = {
                'chrom': result['chrom'],
                'start': result['start'],
                'end': result['end'],
                'vcf_prefix': result['vcf_prefix'],
                'output_bcf': result['output_bcf'].describe()['id'],
                'output_vep': result['output_vep'].describe()['id']}
            coordinate_csv.writerow(writer_dict)

    print("All threads completed...")

    # Getting files back into your project directory on DNAnexus is a two-step process:
    # 1. uploading the local file to the DNA nexus platform to assign it a file-ID (looks like file-ABCDEFGHIJKLMN1234567890)
    # 2. linking this file ID to your project and placing it within your project's directory structure
    # (the subdirectory can be controlled on the command-line by adding a flag to `dx run` like: --destination test/)
    output = {"output_bcfs": [dxpy.dxlink(item) for item in output_bcfs],
              "output_bcf_idxs": [dxpy.dxlink(item) for item in output_bcf_idxs],
              "output_veps": [dxpy.dxlink(item) for item in output_veps],
              "output_vep_idxs": [dxpy.dxlink(item) for item in output_vep_idxs],
              "output_per_samples": [dxpy.dxlink(item) for item in output_per_samples],
              "coordinates_file": dxpy.dxlink(dxpy.upload_local_file(coordinates_name))}

    # This returns all the information about your exit files to the work managing your job via DNANexus:
    return output


dxpy.run()
