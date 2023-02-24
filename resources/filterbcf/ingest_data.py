import os
import dxpy

from pathlib import Path

from general_utilities.association_resources import run_cmd
from general_utilities.job_management.thread_utility import ThreadUtility


class IngestData:

    def __init__(self, input_vcfs: dict, human_reference: dict, human_reference_index: dict, vep_cache: dict,
                 loftee_libraries: dict, gnomad_maf_db: dict, revel_db: dict,
                 cadd_annotations: dict, precomputed_cadd_snvs: dict, precomputed_cadd_indels: dict):

        # Set variables we want to store:
        self._set_vcf_list(input_vcfs)

        # Here we are downloading & unpacking resource files that are required for the annotation engine, they are:
        self._ingest_docker_file()
        self._ingest_human_reference(human_reference, human_reference_index)
        self._ingest_loftee_files(loftee_libraries)
        self._ingest_gnomad_files(gnomad_maf_db)
        self._ingest_revel_files(revel_db)

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

    # This function will locate an associated tbi/csi index:
    @staticmethod
    def find_index(parent_file: dxpy.DXFile, index_suffix: str) -> dxpy.DXFile:

        # Describe the file to get attributes:
        file_description = parent_file.describe(fields={'folder': True, 'name': True, 'project': True})

        # First set the likely details of the corresponding index:
        project_id = file_description['project']
        index_folder = file_description['folder']
        index_name = file_description['name'] + '.' + index_suffix

        # Run a dxpy query.
        # This will fail if no or MULTIPLE indices are found
        index_object = dxpy.find_one_data_object(more_ok=False, classname='file', project=project_id,
                                                 folder=index_folder,
                                                 name=index_name, name_mode='exact')

        # Set a dxfile of the index itself:
        found_index = dxpy.DXFile(dxid=index_object['id'], project=index_object['project'])

        return found_index

    # Bring a prepared docker image into our environment so that we can run commands we need:
    # The Dockerfile to build this image is located at resources/Dockerfile
    @staticmethod
    def _ingest_docker_file() -> None:
        cmd = "docker pull egardner413/mrcepid-burdentesting:latest"
        run_cmd(cmd, is_docker=False)

    # Human reference files – default dxIDs are the location of the GRCh38 reference file on AWS London
    @staticmethod
    def _ingest_human_reference(human_reference: dict, human_reference_index: dict) -> None:
        dxpy.download_dxfile(dxpy.DXFile(human_reference).get_id(), "reference.fasta.gz")
        dxpy.download_dxfile(dxpy.DXFile(human_reference_index).get_id(), "reference.fasta.fai")
        cmd = "gunzip reference.fasta.gz"  # Better to unzip the reference for most commands for some reason...
        run_cmd(cmd, is_docker=False)

    # loftee reference files:
    @staticmethod
    def _ingest_loftee_files(loftee_libraries) -> None:

        os.mkdir("loftee_files/")
        dxpy.download_dxfile(dxpy.DXFile(loftee_libraries), 'loftee_files/loftee_hg38.tar.gz')
        cmd = "tar -zxf loftee_files/loftee_hg38.tar.gz -C loftee_files/"
        run_cmd(cmd, is_docker=False)
        Path('loftee_files/loftee_hg38.tar.gz').unlink()

    # gnomAD MAF files:
    def _ingest_gnomad_files(self, gnomad_maf_db: dict) -> None:

        os.mkdir("gnomad_files/")
        gnomad_dx_file = dxpy.DXFile(gnomad_maf_db)
        dxpy.download_dxfile(gnomad_dx_file.get_id(), 'gnomad_files/gnomad.tsv.gz')
        dxpy.download_dxfile(self.find_index(gnomad_dx_file, 'tbi'), 'gnomad_files/gnomad.tsv.gz.tbi')

        self._write_gnomad_header()

    # REVEL files
    def _ingest_revel_files(self, revel_db):

        os.mkdir("revel_files/")
        revel_dx_file = dxpy.DXFile(revel_db)
        dxpy.download_dxfile(revel_dx_file.get_id(), 'revel_files/new_tabbed_revel_grch38.tsv.gz')
        dxpy.download_dxfile(self.find_index(revel_dx_file, 'tbi').get_id(),
                             'revel_files/new_tabbed_revel_grch38.tsv.gz.tbi')

    # VEP cache file
    def _ingest_vep_cache(self, vep_cache: dict) -> None:
        os.mkdir("vep_caches/")  # This is for legacy reasons to make sure all tests work...

        dxpy.download_dxfile(dxpy.DXFile(vep_cache).get_id(), 'vep_caches/vep_cache.tar.gz')

        cmd = "tar -zxf vep_caches/vep_cache.tar.gz -C vep_caches/"
        run_cmd(cmd, is_docker=False)

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
        run_cmd(cmd, is_docker=False)
        # And finally remove the large annotations tar ball
        Path('cadd_files/annotationsGRCh38_v1.6.tar.gz').unlink()
        print('CADD resources finished downloading and unpacking...')

    def _ingest_precomputed_cadd_files(self, precomputed_cadd_snvs: dict, precomputed_cadd_indels: dict):

        # Now precomputed files...
        os.mkdir("cadd_precomputed/")

        # SNVs...
        cadd_snvs_dx_file = dxpy.DXFile(precomputed_cadd_snvs)
        dxpy.download_dxfile(cadd_snvs_dx_file.get_id(), 'cadd_precomputed/whole_genome_SNVs.tsv.gz')
        dxpy.download_dxfile(self.find_index(cadd_snvs_dx_file, 'tbi').get_id(),
                             'cadd_precomputed/whole_genome_SNVs.tsv.gz.tbi')
        # InDels...
        cadd_indels_dx_file = dxpy.DXFile(precomputed_cadd_indels)
        dxpy.download_dxfile(cadd_indels_dx_file.get_id(), 'cadd_precomputed/gnomad.genomes.r3.0.indel.tsv.gz')
        dxpy.download_dxfile(self.find_index(cadd_indels_dx_file, 'tbi').get_id(),
                             'cadd_precomputed/gnomad.genomes.r3.0.indel.tsv.gz.tbi')

        # Write a header:
        self._write_cadd_header()
        print('Precomputed CADD resources finished downloading and unpacking...')

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
        header_writer = open('cadd.header.txt', 'w')
        header_writer.writelines('##INFO=<ID=CADD,Number=1,Type=Float,Description="CADD Phred Score">' + "\n")
        header_writer.close()
