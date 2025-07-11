import csv
import gzip
from pathlib import Path
from typing import List, TypedDict, Tuple, Optional

import dxpy
from general_utilities.association_resources import find_index
from general_utilities.import_utils.file_handlers.input_file_handler import InputFileHandler
from general_utilities.job_management.command_executor import build_default_command_executor
from general_utilities.job_management.thread_utility import ThreadUtility
from general_utilities.mrc_logger import MRCLogger


class AdditionalAnnotation(TypedDict):
    """A :func:`TypedDict` class to store a tabix-indexed annotation file alongside information about it. See the method
    :func:`IngestData._ingest_annotation()` for more information.
    """
    file: Path
    index: Path
    header_file: Path
    annotation_name: str
    symbol_mode: Optional[str]


class IngestData:

    def __init__(self, input_vcfs: dict, human_reference: dict, human_reference_index: dict, vep_cache: dict,
                 loftee_libraries: dict, additional_annotations: List[dict], thread_utility: ThreadUtility):
        """A class to download resources required for annotation of UKBiobank BCF files. This class is simple and in
        most cases only downloads inputs. It is designed to be used in a multi-threaded environment so that the
        download of very large resource files (e.g. VEP cache) can be done in tandem.

        :param input_vcfs: A txt file containing input VCF/BCFs to be annotated. Individual variant are NOT downloaded
            here, but are instead downloaded in parallel by the calling class.
        :param human_reference: A dxlink to the human reference file in .fasta.gz format
        :param human_reference_index: A dxlink to the human reference index file in .fasta.fai format
        :param vep_cache: A dxlink to a vep cache downloaded from ENSEMBLE
        :param loftee_libraries: A dxlink to a tarball containing LoFTEE libraries
        :param additional_annotations: A list of dxlinks to additional annotation files that are added to the final
            VCF file. These most follow a strict format (see :func:`_ingest_annotation()` for more information)
        """

        self._logger = MRCLogger(__name__).get_logger()

        # Get a Docker image with executables
        self.cmd_executor = build_default_command_executor()

        # Set variables we want to store:
        self.input_vcfs = self._set_vcf_list(input_vcfs)

        # And here we are downloading and processing additional resource files (if required):
        self.annotations: List[AdditionalAnnotation] = []
        for annotation_file in additional_annotations:
            self.annotations.append(self._ingest_annotation(annotation_file))

        # Here we are downloading & unpacking resource files that are required for the annotation engine:
        # These downloads are submitted to a ThreadPoolExecutor so that they download in the background while
        # we perform initial filtering. The monitoring of these threads is handled by the parent class
        thread_utility.launch_job(self._ingest_human_reference, human_reference=human_reference,
                                  human_reference_index=human_reference_index)
        thread_utility.launch_job(self._ingest_loftee_files, loftee_libraries=loftee_libraries)
        thread_utility.launch_job(self._ingest_vep_cache, vep_cache=vep_cache)

    def _set_vcf_list(self, input_vcfs: dict) -> List[str]:
        """Download the input VCF list.

        :param input_vcfs: A txt file containing input VCF/BCFs to be annotated. Individual variant are NOT downloaded
            here, but are instead downloaded in parallel by the calling class.
        :return: A List containing string representations of VCF/BCF files to download
        """

        input_vcf_list = []
        input_vcf_path = InputFileHandler(input_vcfs).get_file_handle()
        with input_vcf_path.open('r') as input_vcf_reader:
            for input_vcf in input_vcf_reader:
                input_vcf_list.append(input_vcf.rstrip())

        self._logger.info(f'Input VCF list contains {len(input_vcf_list)} files...')
        return input_vcf_list

    # Human reference files –
    def _ingest_human_reference(self, human_reference: dict, human_reference_index: dict) -> None:
        """Download the human reference and index files.

        Default dxIDs in dxapp.json are the location of the GRCh38 reference file on AWS London.

        :param human_reference: A dxlink to the human reference file in .fasta.gz format
        :param human_reference_index: A dxlink to the human reference index file in .fasta.fai format
        :return: None
        """
        human_reference_path = InputFileHandler(human_reference, download_now=True).get_file_handle()
        human_reference_index_path = InputFileHandler(human_reference_index, download_now=True).get_file_handle()
        human_reference_path.rename('reference.fasta.gz')
        human_reference_index_path.rename('reference.fasta.fai')

        # Better to unzip the reference for most tools for some reason...
        self.cmd_executor.run_cmd("gunzip reference.fasta.gz")
        self._logger.info('Human reference files finished downloading and unpacking...')

    def _ingest_loftee_files(self, loftee_libraries: dict) -> None:
        """Download the LoFTEE libraries required for the LoFTEE VEP plugin.

        :param loftee_libraries: A dxlink to a tarball containing LoFTEE libraries.
        :return: None
        """

        loftee_dir = Path('loftee_files/')
        loftee_dir.mkdir(exist_ok=True)
        loftee_path = InputFileHandler(loftee_libraries, download_now=True).get_file_handle()

        cmd = f'tar -zxf {loftee_path} -C {loftee_dir}/'
        self.cmd_executor.run_cmd(cmd)
        loftee_path.unlink()
        self._logger.info('LoFTEE resources finished downloading and unpacking...')

    def _ingest_vep_cache(self, vep_cache: dict) -> None:
        """Download the VEP cache from ENSEMBLE and unpack it. This function also writes a VEP BCF header for use with
        BCFtools annotate.

        :param vep_cache: A dxlink to a vep cache downloaded from ENSEMBLE
        :return: None
        """
        vep_dir = Path('vep_caches/')
        vep_dir.mkdir(exist_ok=True)  # This is for legacy reasons to make sure all tests work...
        vep_path = InputFileHandler(vep_cache, download_now=True).get_file_handle()

        cmd = f'tar -zxf {vep_path} -C {vep_dir}/'
        self.cmd_executor.run_cmd(cmd)

        # And purge the large tarball
        vep_path.unlink()
        self._logger.info('VEP resources finished downloading and unpacking...')

    @staticmethod
    def _determine_annotation_type(values: List[str]) -> Tuple[str, str]:
        """Take a list of values extracted from an annotation .tsv and determines the corresponding VCF format.

        As specified by the VCF format, Format can only be one of Integer, Float, Flag, or String. The number of values
        (if delimited by ',') is also determined. If the number of values is not consistent, the number is set to '.'.

        :param values: A list of values from an annotation file.
        :return: A tuple containing the VCF format and the number of values.
        """

        # First check if the value is actually a list of values
        lengths = set()
        new_values = []
        found_comma = False
        for value in values:
            if ',' in value:
                found_comma = True
                split_values = value.split(',')
                new_values.extend(split_values)
                lengths.add(len(split_values))

        # Determine if the number of values is consistent, if no comma was identified, then length MUST be 1
        if found_comma and len(lengths) == 1:
            annotation_number = str(lengths.pop())
            values = new_values
        elif found_comma and len(lengths) > 1:
            annotation_number = '.'
            values = new_values
        else:
            annotation_number = '1'

        # Now determine the format of the annotation
        types = set()
        data_types = [(int, 'Integer'), (float, 'Float'), (bool, 'Flag'), (str, 'String')]
        for value in values:
            # Missing values do not determine annotation type
            if value.lower() in ['na', 'nan', '.']:
                continue

            for dt, label in data_types:
                try:
                    dt(value)
                    types.add(label)
                    break
                except ValueError:
                    continue

        # Determine if annotation types are consistent, if not (except for int + float combo), then we default to string
        if len(types) == 1:
            annotation_type = types.pop()
        elif len(types) == 2 and 'Integer' in types and 'Float' in types:
            annotation_type = 'Float'
        else:
            annotation_type = 'String'

        return annotation_type, annotation_number

    def _ingest_annotation(self, annotation_dxfile: dict) -> AdditionalAnnotation:
        """Download and process an additional annotation file.

        This function also writes a VCF header specific to this annotation for use with bcftools annotate by taking
        the first 1000 records and inferring the data type of the annotation using :func:`_determine_annotation_type()`.

        :param annotation_dxfile: A dxlink to an additional annotation file
        :return: An AdditionalAnnotation object
        """

        index_dxfile = find_index(annotation_dxfile, 'tbi')
        annotation_path = InputFileHandler(annotation_dxfile, download_now=True).get_file_handle()
        index_path = InputFileHandler(index_dxfile, download_now=True).get_file_handle()
        annotation_name = ''

        with gzip.open(annotation_path, 'rt') as annotation_reader:

            annotation_csv = csv.DictReader(annotation_reader, delimiter='\t')
            file_header = annotation_csv.fieldnames

            # File header should only have 5 items:
            # 0 = CHROM
            # 1 = POS
            # 2 = REF
            # 3 = ALT
            # 4 = Annotation itself
            # 5 = Optional SYMBOL tag for gene-specific annotation
            if (file_header[:4] == ['CHROM', 'POS', 'REF', 'ALT'] or file_header[:4] == ['#CHROM', 'POS', 'REF', 'ALT']) and len(file_header) in range(5, 7):

                # Set the annotation name from the column name
                annotation_name = file_header[4]

                # Check if we have an SYMBOL value for each score
                symbol_mode = None
                if len(file_header) == 6:
                    if file_header[5] == 'SYMBOL':
                        symbol_mode = 'SYMBOL'
                    elif file_header[5] == 'ENST':
                        symbol_mode = 'ENST'
                    else:
                        raise dxpy.AppError(f'Sixth column included in annotations file – {annotation_path} – but it '
                                            f'is NOT a SYMBOL / ENST column. Please check the file.')

                # Attempt to infer the datatype of the annotation in the annotation file:
                values = []
                for i, record in enumerate(annotation_csv):
                    if i > 1000:
                        break
                    else:
                        values.append(record[annotation_name])

                # Get information about the information in the annotation file
                annotation_type, annotation_number = self._determine_annotation_type(values)

                header_file = Path(f'{annotation_name}.header.txt')
                with header_file.open('w') as header_writer:
                    # Because I am never using the VCF to do any filtering on these fields, we use the most unbounded
                    # 'Number' and 'Type' fields so essentially anything can be entered by the user
                    header_writer.write(f'##INFO=<ID={annotation_name},'
                                        f'Number={annotation_number},'
                                        f'Type={annotation_type},'
                                        f'Description="{annotation_name} annotation.">\n')
            else:
                raise dxpy.AppError(f'Number / name of columns in annotations file {annotation_path} is incorrect – '
                                    f'header = {file_header}')

        return {'file': annotation_path, 'index': index_path, 'annotation_name': annotation_name,
                'header_file': header_file, 'symbol_mode': symbol_mode}
