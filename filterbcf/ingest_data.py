import gzip
import os
from typing import List, TypedDict

import dxpy

from pathlib import Path

from general_utilities.association_resources import find_index, download_dxfile_by_name
from general_utilities.job_management.command_executor import build_default_command_executor
from general_utilities.job_management.thread_utility import ThreadUtility


class AdditionalAnnotation(TypedDict):
    file: Path
    index: Path
    header_file: Path
    annotation_name: str


class IngestData:

    def __init__(self, input_vcfs: dict, human_reference: dict, human_reference_index: dict, vep_cache: dict,
                 loftee_libraries: dict,
                 cadd_annotations: dict, precomputed_cadd_snvs: dict, precomputed_cadd_indels: dict,
                 additional_annotations: List[dict]):

        # Get a Docker image with executables
        self.cmd_executor = build_default_command_executor()

        # Set variables we want to store:
        self._set_vcf_list(input_vcfs)

        # Here we are downloading & unpacking resource files that are required for the annotation engine, they are:
        self._ingest_human_reference(human_reference, human_reference_index)
        self._ingest_loftee_files(loftee_libraries)

        # And here we are downloading and processing additional resource files (if required):
        self.annotations: List[AdditionalAnnotation] = []
        for annotation_file in additional_annotations:
            self.annotations.append(self._ingest_annotation(annotation_file))

        # These two downloads are submitted to a ThreadPoolExecutor so that they download in the background while
        # we perform initial filtering. The monitoring of these threads is handled by the parent class
        # (unsure if good programming practice or not...?)
        thread_utility = ThreadUtility(thread_factor=2, error_message='A data ingest thread failed...', incrementor=1)
        thread_utility.launch_job(self._ingest_vep_cache, vep_cache=vep_cache)
        thread_utility.launch_job(self._ingest_cadd_files, cadd_annotations=cadd_annotations)
        thread_utility.launch_job(self._ingest_precomputed_cadd_files,
                                  precomputed_cadd_snvs=precomputed_cadd_snvs,
                                  precomputed_cadd_indels=precomputed_cadd_indels)
        thread_utility.collect_futures()

    def _set_vcf_list(self, input_vcfs):

        self.input_vcfs = []
        input_vcfs = dxpy.DXFile(input_vcfs)
        dxpy.download_dxfile(input_vcfs.get_id(), "vcf_list.txt")  # Actually download the file
        input_vcf_reader = open("vcf_list.txt", 'r')
        for input_vcf in input_vcf_reader:
            self.input_vcfs.append(input_vcf.rstrip())
        input_vcf_reader.close()

    # Human reference files – default dxIDs are the location of the GRCh38 reference file on AWS London
    def _ingest_human_reference(self, human_reference: dict, human_reference_index: dict) -> None:
        dxpy.download_dxfile(dxpy.DXFile(human_reference).get_id(), "reference.fasta.gz")
        dxpy.download_dxfile(dxpy.DXFile(human_reference_index).get_id(), "reference.fasta.fai")
        cmd = "gunzip reference.fasta.gz"  # Better to unzip the reference for most commands for some reason...
        self.cmd_executor.run_cmd(cmd)

    # loftee reference files:
    def _ingest_loftee_files(self, loftee_libraries) -> None:

        os.mkdir("loftee_files/")
        dxpy.download_dxfile(dxpy.DXFile(loftee_libraries), 'loftee_files/loftee_hg38.tar.gz')
        cmd = "tar -zxf loftee_files/loftee_hg38.tar.gz -C loftee_files/"
        self.cmd_executor.run_cmd(cmd)
        Path('loftee_files/loftee_hg38.tar.gz').unlink()

    # VEP cache file
    def _ingest_vep_cache(self, vep_cache: dict) -> None:
        os.mkdir("vep_caches/")  # This is for legacy reasons to make sure all tests work...

        dxpy.download_dxfile(dxpy.DXFile(vep_cache).get_id(), 'vep_caches/vep_cache.tar.gz')

        cmd = "tar -zxf vep_caches/vep_cache.tar.gz -C vep_caches/"
        self.cmd_executor.run_cmd(cmd)

        # Write the header for use with bcftools annotate
        self._write_vep_header()

        # And purge the large tarball
        Path('vep_caches/vep_cache.tar.gz').unlink()
        print('VEP resources finished downloading and unpacking...')

    # CADD reference files – These are the resource files so InDel CADD scores can be calculated from scratch
    # CADD known reference files - pre-computed sites files so we don't have to recompute CADD for SNVs/known InDels
    def _ingest_cadd_files(self, cadd_annotations: dict) -> None:

        # First the annotations
        os.mkdir("cadd_files/")
        dxpy.download_dxfile(dxpy.DXFile(cadd_annotations).get_id(), 'cadd_files/annotationsGRCh38_v1.6.tar.gz')
        cmd = "tar -zxf cadd_files/annotationsGRCh38_v1.6.tar.gz -C cadd_files/"
        self.cmd_executor.run_cmd(cmd)
        # And finally remove the large annotations tar ball
        Path('cadd_files/annotationsGRCh38_v1.6.tar.gz').unlink()
        print('CADD resources finished downloading and unpacking...')

    def _ingest_precomputed_cadd_files(self, precomputed_cadd_snvs: dict, precomputed_cadd_indels: dict):

        # Now precomputed files...
        os.mkdir("cadd_precomputed/")

        # SNVs...
        cadd_snvs_dx_file = dxpy.DXFile(precomputed_cadd_snvs)
        cadd_snvs_dx_index = find_index(cadd_snvs_dx_file, 'tbi')
        dxpy.download_dxfile(cadd_snvs_dx_file.get_id(), 'cadd_precomputed/whole_genome_SNVs.tsv.gz')
        dxpy.download_dxfile(cadd_snvs_dx_index.get_id(), 'cadd_precomputed/whole_genome_SNVs.tsv.gz.tbi')
        # InDels...
        cadd_indels_dx_file = dxpy.DXFile(precomputed_cadd_indels)
        cadd_index_dx_index = find_index(cadd_indels_dx_file, 'tbi')
        dxpy.download_dxfile(cadd_indels_dx_file.get_id(), 'cadd_precomputed/gnomad.genomes.r3.0.indel.tsv.gz')
        dxpy.download_dxfile(cadd_index_dx_index.get_id(), 'cadd_precomputed/gnomad.genomes.r3.0.indel.tsv.gz.tbi')

        # Write a header:
        self._write_cadd_header()
        print('Precomputed CADD resources finished downloading and unpacking...')

    def _ingest_annotation(self, annotation_dxfile: dict) -> AdditionalAnnotation:

        annotation_path = download_dxfile_by_name(annotation_dxfile)
        index_dxfile = find_index(dxpy.DXFile(annotation_dxfile), 'tbi')
        annotation_name = ''
        index_path = download_dxfile_by_name(index_dxfile)

        with gzip.open(annotation_path, 'rt') as annotation_reader:
            file_header = annotation_reader.readlines(1)[0].rstrip()

        # File header should only have 5 items:
        # 0 = CHROM
        # 1 = POS
        # 2 = REF
        # 3 = ALT
        # 4 = Annotation itself
        if len(file_header) == 5 and file_header[0] == 'CHROM' and file_header[1] == 'POS' and file_header[2] == 'REF' and file_header[3] == 'ALT':
            annotation_name = file_header[4]
            header_file = Path(f'{annotation_name}.header.txt')
            with header_file.open('w') as header_writer:
                # Because I am never using the VCF to do any filtering on these fields, we use the most unbounded
                # 'Number' and 'Type' fields so essentially anything can be entered by the user
                header_writer.write(f'##INFO=<ID={annotation_name},Number=.,Type=String,Description="{annotation_name} annotation.">\n')

        else:
            raise dxpy.AppError(f'File format of annotations file {annotation_path} is incorrect – '
                                f'header = {file_header}')

        return {'file': annotation_path, 'index': index_path, 'annotation_name': annotation_name,
                'header_file': header_file}


    # Writes a VCF style header that is compatible with bcftools annotate for adding VEP info back into our filtered VCF
    @staticmethod
    def _write_vep_header() -> None:
        header_writer = open('vep_vcf.header.txt', 'w')
        header_writer.writelines('##INFO=<ID=MANE,Number=1,Type=String,Description="Canonical MANE Transcript">' + "\n")
        header_writer.writelines('##INFO=<ID=ENST,Number=1,Type=String,Description="Canonical Ensembl '
                                 'Transcript">' + "\n")
        header_writer.writelines('##INFO=<ID=ENSG,Number=1,Type=String,Description="Canonical Ensembl Gene">' + "\n")
        header_writer.writelines('##INFO=<ID=BIOTYPE,Number=1,Type=String,Description="Biotype of ENSG as in '
                                 'VEP">' + "\n")
        header_writer.writelines('##INFO=<ID=SYMBOL,Number=1,Type=String,Description="HGNC Gene ID">' + "\n")
        header_writer.writelines('##INFO=<ID=CSQ,Number=1,Type=String,Description="Most severe VEP CSQ for this '
                                 'variant">' + "\n")
        header_writer.writelines('##INFO=<ID=gnomAD_AF,Number=1,Type=Float,Description="gnomAD v3.0 Exomes AF. If 0,'
                                 ' variant does not exist in gnomAD">' + "\n")
        header_writer.writelines('##INFO=<ID=CADD,Number=1,Type=Float,Description="CADD Phred Score">' + "\n")
        header_writer.writelines('##INFO=<ID=REVEL,Number=1,Type=Float,Description="REVEL Score. '
                                 'NaN if not a missense variant.">' + "\n")
        header_writer.writelines('##INFO=<ID=SIFT,Number=1,Type=String,Description="SIFT Score. '
                                 'NA if not a missense variant.">' + "\n")
        header_writer.writelines('##INFO=<ID=POLYPHEN,Number=1,Type=String,Description="PolyPhen Score. '
                                 'NA if not a missense variant.">' + "\n")
        header_writer.writelines('##INFO=<ID=LOFTEE,Number=1,Type=String,Description="LOFTEE annotation if LoF CSQ. '
                                 'NA if not a PTV.">' + "\n")
        header_writer.writelines('##INFO=<ID=AA,Number=1,Type=String,Description="Amino acid change for this '
                                 'variant.">' + "\n")
        header_writer.writelines('##INFO=<ID=AApos,Number=1,Type=String,Description="Amino acid location in target '
                                 'protein for change indicated by AA">' + "\n")
        header_writer.writelines('##INFO=<ID=PARSED_CSQ,Number=1,Type=String,Description="Parsed simplified '
                                 'CSQ">' + "\n")
        header_writer.writelines('##INFO=<ID=MULTI,Number=1,Type=String,Description="Is variant multiallelic?">' + "\n")
        header_writer.writelines('##INFO=<ID=INDEL,Number=1,Type=String,Description="Is variant an InDel?">' + "\n")
        header_writer.writelines('##INFO=<ID=MINOR,Number=1,Type=String,Description="Minor Allele">' + "\n")
        header_writer.writelines('##INFO=<ID=MAJOR,Number=1,Type=String,Description="Major Allele">' + "\n")
        header_writer.writelines('##INFO=<ID=MAF,Number=1,Type=Float,Description="Minor Allele Frequency">' + "\n")
        header_writer.writelines('##INFO=<ID=MAC,Number=1,Type=Float,Description="Minor Allele Count">' + "\n")
        header_writer.close()

    # And write a brief gnomad VCF header file so that bcftools knows how to process this data:
    @staticmethod
    def _write_gnomad_header() -> None:
        gnomad_header_writer = open('gnomad_files/gnomad.header.txt', 'w')
        gnomad_header_writer.writelines(
            '##INFO=<ID=gnomAD_MAF,Number=1,Type=Float,Description="gnomAD Exomes AF">' + "\n")
        gnomad_header_writer.close()

    # Write a header for cadd annotation with bcftools annotate
    @staticmethod
    def _write_cadd_header() -> None:
        with Path('cadd.header.txt').open('w') as header_writer:
            header_writer.write(f'##INFO=<ID=CADD,Number=1,Type=Float,Description="CADD Phred Score">\n"')
