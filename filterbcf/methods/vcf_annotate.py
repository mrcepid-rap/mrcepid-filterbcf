import csv
import gzip
import itertools

from pathlib import Path
from typing import TypedDict, Tuple, List
from pathlib import Path

from general_utilities.association_resources import generate_linked_dx_file, bgzip_and_tabix, replace_multi_suffix, \
    check_gzipped
from filterbcf.methods.ingest_data import AdditionalAnnotation
from general_utilities.association_resources import generate_linked_dx_file, bgzip_and_tabix, replace_multi_suffix, \
    check_gzipped
from general_utilities.job_management.command_executor import CommandExecutor
from general_utilities.mrc_logger import MRCLogger

from filterbcf.methods.ingest_data import AdditionalAnnotation


class VCFAnnotate:

    def __init__(self, vcf_path: Path, filtered_vcf: Path, additional_annotations: List[AdditionalAnnotation],
                 cmd_executor: CommandExecutor):

        self._logger = MRCLogger(__name__).get_logger()
        self._cmd_executor = cmd_executor
        self._cadd_executor = cadd_executor

        # Set the prefix for the VCF file. This allows us to just use the prefix for all file names throughout this
        # class
        self.vcfprefix = vcf_path.stem

        # 1. Generate sites files that can be run through the various annotation parts
        sites_vcf = self._generate_sites_file(filtered_vcf)

        # 1a. Get information about this VCF from the VEP sites file since it is much quicker than parsing the full VCF
        # with bcftools:
        (self.chunk_chrom, self.chunk_start, self.chunk_end) = self._get_bcf_information(sites_vcf)

        # 2. This does the initial VEP run
        annotated_vcf = self._run_vep(sites_vcf)

        # 3. Generate an annotations TSV for the VCF
        vep_tsv = self._generate_annotations_tsv(annotated_vcf)

        # 4. Add annotations that need to be done separately from VEP (e.g., gnomAD)
        # Note that 'vep_tsv' changes each time we add an annotation, but the final file should be the same. This was
        # done for simplicity and does not represent a 'vep_tsv' per se.
        annotation_names = []
        for annotation in additional_annotations:
            vep_tsv, annotation_name = self._add_additional_annotation(vep_tsv, annotation)
            annotation_names.append(annotation_name)

        # 5. Generating merged/final files:
        # This function parses the information from the raw VEP run (via run_vep()) and adds it to our filtered vcf
        vep_gz, vep_gz_tbi = self._parse_vep(vep_tsv, annotation_names)

        # 6. Set final file outputs for this entire process:
        self.finalbcf = generate_linked_dx_file(filtered_vcf)
        self.finalbcf_index = generate_linked_dx_file(f'{filtered_vcf}.csi')
        self.finalvep = generate_linked_dx_file(vep_gz)
        self.finalvep_index = generate_linked_dx_file(vep_gz_tbi)

    # Just a class to ensure typing hint compatibility when defining consequence scores
    class ConsequenceSeverity(TypedDict):
        score: int
        type: str

    # Decide how "severe" a CSQ is for a given annotation record
    @staticmethod
    def _define_score(csqs: str) -> ConsequenceSeverity:

        csqs_split = csqs.split('&')

        # This is a Eugene-decided level of severity partially decided based on VEP severity score.
        # This should contain all possible CSQ annotations for a variant other than those reserved
        # for large structural variants
        vep_consequences = {'stop_gained': {'score': 1, 'type': 'PTV'},
                            'frameshift_variant': {'score': 2, 'type': 'PTV'},
                            'splice_acceptor_variant': {'score': 3, 'type': 'PTV'},
                            'splice_donor_variant': {'score': 4, 'type': 'PTV'},
                            'stop_lost': {'score': 5, 'type': 'STOP_LOST'},
                            'start_lost': {'score': 6, 'type': 'START_LOST'},
                            'inframe_insertion': {'score': 7, 'type': 'INFRAME'},
                            'inframe_deletion': {'score': 8, 'type': 'INFRAME'},
                            'missense_variant': {'score': 9, 'type': 'MISSENSE'},
                            'protein_altering_variant': {'score': 10, 'type': 'INFRAME'},
                            'splice_region_variant': {'score': 11, 'type': 'NONCODING'},
                            'incomplete_terminal_codon_variant': {'score': 12, 'type': 'INFRAME'},
                            'start_retained_variant': {'score': 13, 'type': 'SYN'},
                            'stop_retained_variant': {'score': 14, 'type': 'SYN'},
                            'synonymous_variant': {'score': 15, 'type': 'SYN'},
                            '5_prime_UTR_variant': {'score': 16, 'type': 'UTR'},
                            '3_prime_UTR_variant': {'score': 17, 'type': 'UTR'},
                            'intron_variant': {'score': 18, 'type': 'INTRONIC'},
                            'intergenic_variant': {'score': 19, 'type': 'INTERGENIC'},
                            'upstream_gene_variant': {'score': 20, 'type': 'INTERGENIC'},
                            'downstream_gene_variant': {'score': 21, 'type': 'INTERGENIC'},
                            'no_score': {'score': 22, 'type': 'ERROR'}}

        ret_csq = vep_consequences['no_score']
        for c in csqs_split:
            if c in vep_consequences:
                if vep_consequences[c]['score'] < ret_csq['score']:
                    ret_csq = vep_consequences[c]

        return ret_csq

    # Prepares a record for final printing by adding some additional pieces of information and
    # finalising the names of some columns for easy printing.
    @staticmethod
    def _final_process_record(rec: dict, severity: ConsequenceSeverity, annotation_names: List[str]) -> dict:

        # Rename some columns for printing purposes
        rec['PARSED_CSQ'] = severity['type']

        # I use the MA flag set in 'splitvcf' to determine if a variant is multiallelic or not
        if rec['MA'] == '.':
            rec['IS_MULTIALLELIC'] = False
        else:
            rec['IS_MULTIALLELIC'] = True

        # Records who do not have equal REF/ALT length are assigned as InDels
        if len(rec['REF']) != len(rec['ALT']):
            rec['IS_INDEL'] = True
        else:
            rec['IS_INDEL'] = False

        # This corrects an issue with sites with 100% missingness that BCFtools doesn't handle correctly
        if rec['AF'] == '.':
            rec['AF'] = 0
            rec['FILTER'] = 'FAIL'
        else:
            rec['AF'] = float(rec['AF'])

        # Set sensible defaults for a variety of fields:
        # "gnomad_maf", "REVEL" "SIFT", "PolyPhen", "LoF"
        # rec['gnomad_maf'] = rec['gnomad_maf'] if rec['gnomad_maf'] != '.' else '0' # Sites w/o gnomAD don't exist in gnomAD so a MAF of 0 seems appropriate
        # rec['REVEL'] = rec['REVEL'] if rec['REVEL'] != '.' else 'NaN' # NaN is default VCF spec for missing floats
        rec['SIFT'] = rec['SIFT'] if rec['SIFT'] != '.' else 'NA' # NA is default VCF spec for missing strings
        rec['POLYPHEN'] = rec['POLYPHEN'] if rec['POLYPHEN'] != '.' else 'NA' # NA is default VCF spec for missing strings
        rec['LOFTEE'] = rec['LOFTEE'] if rec['LOFTEE'] != '.' else 'NA' # NA is default VCF spec for missing strings
        rec['AA'] = rec['AA'] if rec['AA'] != '.' else 'NA'  # NA is default VCF spec for missing strings
        rec['AApos'] = rec['AApos'] if rec['AApos'] != '.' else 'NA'  # NA is default VCF spec for missing strings

        # Setting additional tags requested by GWAS-y people
        if float(rec['AF']) < 0.5:
            rec['MINOR_ALLELE'] = rec['ALT']
            rec['MAJOR_ALLELE'] = rec['REF']
            rec['MAF'] = '%s' % rec['AF']
            rec['MAC'] = '%s' % rec['AC']
        else:
            rec['MINOR_ALLELE'] = rec['REF']
            rec['MAJOR_ALLELE'] = rec['ALT']
            rec['MAF'] = '%s' % (1 - float(rec['AF']))
            rec['MAC'] = '%s' % (int(rec['AN']) - int(rec['AC']))

        # Just make sure that additional annotations are set to NA rather than '.' for downstream compatability:
        for name in annotation_names:
            rec[name] = rec[name] if rec[name] != '.' else 'NA'

        return rec

    def _generate_sites_file(self, input_vcf: Path) -> Path:
        """This function generates sites files for VEP using bcftools view -G.

        :param input_vcf: The final filtered VCF from VCFFilter
        :return: Path to the sites file
        """

        # Generate a file without genotypes for VEP
        # This should be the ONLY file that is .vcf.gz format as it is required by VEP
        # -G : strips genotypes
        output_sites = Path(f'{self.vcfprefix}.sites.vcf.gz')
        cmd = f'bcftools view --threads 4 -G -Oz ' \
              f'-o /test/{output_sites} ' \
              f'/test/{input_vcf}'
        self._cmd_executor.run_cmd_on_docker(cmd)

        return output_sites

    @staticmethod
    def _get_bcf_information(input_vcf: Path) -> Tuple[str, int, int]:
        """Retrieve the chromosome and start and end coordinates of a BCF file.

        :param input_vcf: Input VCF sites file. Here we use the VEP vcf sites file so we don't have to parse GTs
        :return: A Tuple containing the chromosome, start, and end coordinates
        """

        with check_gzipped(input_vcf) as sites_reader:
            line_num = 1
            for line in sites_reader:
                if line.startswith('#'):
                    continue
                else:
                    site_data = line.rstrip().split('\t')
                    if line_num == 1:
                        chunk_chrom = site_data[0]
                        chunk_start = int(site_data[1])
                    chunk_end = int(site_data[1])
                    line_num += 1

        return chunk_chrom, chunk_start, chunk_end

    def _run_vep(self, input_vcf: Path) -> Path:
        """Run VEP on the sites file generated by _generate_sites_files()

        I am not going to document each individual thing for VEP, but all are available on the VEP webiste

        :param input_vcf: The sites file generated by _generate_sites_files()
        :return: The VEP-annotated VCF file
        """

        output_vcf = Path(f'{self.vcfprefix}.sites.vep.vcf.gz')
        cmd = f'perl -Iensembl-vep/cache/Plugins/loftee/ -Iensembl-vep/cache/Plugins/loftee/maxEntScan/ ' \
              f'ensembl-vep/vep --offline --cache --assembly GRCh38 --dir_cache /test/vep_caches/ ' \
              f'--everything --allele_num --fork 4 ' \
              f'-i /test/{input_vcf} --format vcf --fasta /test/reference.fasta ' \
              f'-o /test/{output_vcf} --compress_output bgzip --vcf ' \
              f'--dir_plugins ensembl-vep/cache/Plugins/ ' \
              f'--plugin LoF,loftee_path:ensembl-vep/cache/Plugins/loftee,human_ancestor_fa:/test/loftee_files/loftee_hg38/human_ancestor.fa.gz,conservation_file:/test/loftee_files/loftee_hg38/loftee.sql,gerp_bigwig:/test/loftee_files/loftee_hg38/gerp_conservation_scores.homo_sapiens.GRCh38.bw'
        self._cmd_executor.run_cmd_on_docker(cmd)
        input_vcf.unlink()

        return output_vcf

    def _run_cadd(self, cadd_vcf: Path) -> AdditionalAnnotation:
        """This method runs CADD on provided variants. It DOES NOT add the annotations to the VCF.

        After this method is called, the AdditionalAnnotation is then provided to :func:`_add_additional_annotation` to
        add scores to the VCF with bcftools annotate.

        :param cadd_vcf: The VCF file in a CADD-friendly format
        :return: An AdditionalAnnotation object containing the CADD annotation information
        """

        # Run CADD on the sites file from generate_sites_files():
        cadd_tsv = Path(f'{self.vcfprefix}.cadd.tsv')
        cmd = f'CADD-scripts/CADD.sh -c 4 -g GRCh38 ' \
              f'-o /test/{cadd_tsv} ' \
              f'/test/{cadd_vcf}'

        self._cadd_executor.run_cmd_on_docker(cmd)

        # Add chr back so BCFtools can understand for reannotation, simplify the output, and then bgzip and tabix
        # index
        cadd_chr = Path(f'{self.vcfprefix}.cadd.chr.tsv')
        with gzip.open(cadd_tsv, 'rt') as cadd_reader, \
                cadd_chr.open('w') as cadd_out:

            cadd_csv = csv.DictWriter(cadd_out, delimiter='\t', fieldnames=['#CHROM', 'POS', 'REF', 'ALT', 'CADD'])
            cadd_csv.writeheader()

            for line in cadd_reader:
                data = line.rstrip().split('\t')
                if data[0].startswith('#'):
                    continue
                cadd_csv.writerow({'#CHROM': f'chr{data[0]}',
                                   'POS': data[1],
                                   'REF': data[2],
                                   'ALT': data[3],
                                   'CADD': data[5]})

        cadd_gz, cadd_idx = bgzip_and_tabix(cadd_chr, end_row=2, comment_char='"#"')

        cadd_annotation: AdditionalAnnotation = {'annotation_name': 'CADD',
                                                 'file': cadd_gz,
                                                 'index': cadd_idx,
                                                 'header_file': Path('cadd.header.txt'),
                                                 'symbol_mode': None}

        # Remove CADD intermediate files to save space
        cadd_vcf.unlink()
        cadd_tsv.unlink()

        return cadd_annotation

    def _add_additional_annotation(self, input_tsv, annotation: AdditionalAnnotation) -> Tuple[Path, str]:
        """This method will take a tab-delimited, tabix-indexed .tsv file with a single annotation, and add it to the
        end of the provided annotations tsv while matching on gene SYMBOL.

        Documenting 'things I know' about how this works for future me:

        * bcftools annotate will not work here in it's current form. This is because it does not allow for matching on
        additional fields. This causes issues when the annotation I care about is specific to a transcript / protein
        and is NOT just a global annotation for that given variant. Examples included AlphaMissense & EVE

        * annot-tsv is a tool found within htslib. Importantly, this tool DOES NOT CARE about the wierd column
        header generated by +split-vep, where it adds the column number to the header. The following method
        (:func:`_parse_vep`) ALSO does not care since I manually set the header.

        * I have also tested if annot-tsv is fast with the entire annotation file, and NOT subset to the chromosome /
        gene region. It is fast!

        * This method is hard-coded to use gene SYMBOL as the matching field. This is because gene SYMBOL is the most
        flexible, even if 'transcript' is more accurate.

        :param input_tsv: The input TSV to be annotated with the additional annotation.
        :param annotation: A AdditionalAnnotation 'object' that contains information about a given supplemental
            annotation.
        :return: A Tuple containing the output TSV with additional annotation and the name of the annotation
        """

        # Memory requires for annot-tsv appear to be very high for large files. I am slicing the annotation file down
        # to the same region as the vcf to see if I can reduce memory requirements.
        sliced_tsv = replace_multi_suffix(annotation["file"], f'.{self.vcfprefix}.sliced.tsv')
        cmd = f'tabix --threads 4 -h /test/{annotation["file"]} ' \
              f'"{self.chunk_chrom}:{self.chunk_start - 1}-{self.chunk_end}" > {sliced_tsv}' # -1 due to 0-based idx
        self._cmd_executor.run_cmd_on_docker(cmd)

        sliced_bgzip, _ = bgzip_and_tabix(sliced_tsv, comment_char='"#"', end_row=2)

        if annotation['symbol_mode']:
            # VEP uses 'Feature' to denote the gene transcript, but we need to match on ENST
            if annotation['symbol_mode'] == 'ENST':
                target = 'Feature'
            else:
                target = 'SYMBOL'
            match_string = f'REF,ALT,{annotation["symbol_mode"]}:REF,ALT,{target}'
        else:
            match_string = 'REF,ALT:REF,ALT'

        output_tsv = Path(f'{self.vcfprefix}.{annotation["annotation_name"]}.vep_table.tsv')
        cmd = f'annot-tsv ' \
              f'-c CHROM,POS,POS ' \
              f'-f {annotation["annotation_name"]} ' \
              f'-m {match_string} ' \
              f'-t /test/{input_tsv} ' \
              f'-s /test/{sliced_bgzip} ' \
              f'-o /test/{output_tsv}'
        self._cmd_executor.run_cmd_on_docker(cmd)
        input_tsv.unlink()

        return output_tsv, annotation['annotation_name']

    def _generate_annotations_tsv(self, input_vcf: Path) -> Path:
        """Generate an output TSV after doing all annotations that we can parse for most deleterious consequence.
        
        The purpose of this is to generate a TSV file of annotations that we care about to parse later :func:`_parse_vep`
        +split-vep is a bcftools plugin that iterates through VEP fields provided in a VCF via the CSQ INFO field.
        -d :  outputs duplicate transcripts on separate lines. In other words, a gene may have multiple transcripts,
        and put each transcript and CSQ on a separate line in the tsv file
        -f : CSQ fields that we want to output into the tsv file
        
        :param input_vcf: VCF file that contains all annotations to add
        :return: Pathlike object of the annotation tsv
        """

        # Set field names that we need to extract from the annotations file
        to_extract = ['CHROM', 'POS', 'REF', 'ALT', 'ID', 'FILTER', 'INFO/AF', 'F_MISSING', 'AN', 'AC', 'MA',
                      'MANE_SELECT', 'Feature', 'Gene', 'BIOTYPE', 'CANONICAL', 'SYMBOL', 'Consequence', 'SIFT',
                      'PolyPhen', 'LoF', 'Amino_acids', 'Protein_position', 'GTM', 'GT0', 'GT1', 'GT2']
        # Have to do this b/c '%' is req'd and join won't append to the front of the first element
        to_extract_string = '\\t'.join([f'%{x}' for x in to_extract])

        vep_table = Path(f'{self.vcfprefix}.vep_table.tsv')
        cmd = f'bcftools +split-vep -df \'{to_extract_string}\\n\' -H ' \
              f'-o /test/{vep_table} ' + \
              f'/test/{input_vcf}'
        self._cmd_executor.run_cmd_on_docker(cmd)
        input_vcf.unlink()

        return vep_table

    def _parse_vep(self, raw_vep: Path, annotation_names: List[str]) -> Tuple[Path, Path]:
        """This function parses VEP output for most severe CSQ for each variant. See individual comments in this code to
         understand how that is done.

        :param raw_vep: The raw VEP output file directly from bcftools split-vep
        :param annotation_names: A list of additional annotation names to append to the end of our required annotations
        :return: A Tuple containing the path to the gzipped and tabixed VEP output file and index file
        """

        # First we need to set the I/O headers for this process. This is for two reasons:
        # 1. BCFTools can only print headers based on field names, and we want to change these to fieldnames that are
        #    more human-readable.
        # 2. We need to add additional fields documenting the processing done in this method.

        # Reader:
        reader_header = ['CHROM', 'POS', 'REF', 'ALT', 'ID', 'FILTER', 'AF', 'F_MISSING', 'AN', 'AC', 'MA',
                         'MANE', 'ENST', 'ENSG', 'BIOTYPE', 'CANONICAL', 'SYMBOL', 'CSQ', 'SIFT',
                         'POLYPHEN', 'LOFTEE', 'AA', 'AApos', 'GTM', 'GT0', 'GT1', 'GT2']
        reader_header.extend(annotation_names)

        # Writer:
        writer_header = reader_header.copy()
        writer_header.extend(['PARSED_CSQ', 'IS_MULTIALLELIC', 'IS_INDEL', 'MINOR_ALLELE', 'MAJOR_ALLELE',
                              'MAF', 'MAC'])

        # Define files:
        annote_file = Path(f'{self.vcfprefix}.vep.tsv')

        with raw_vep.open('r') as annote_reader,\
                annote_file.open('w') as annote_writer:

            reader_csv = csv.DictReader(annote_reader, delimiter='\t', fieldnames=reader_header, quoting=csv.QUOTE_NONE)
            writer_csv = csv.DictWriter(annote_writer, delimiter="\t", fieldnames=writer_header, extrasaction='ignore')
            writer_csv.writeheader()

            # Now we need to iterate through the .tsv file that we made
            # A single variant can be spread across multiple rows as it can be annotated for multiple transcripts in the same gene
            # The annotations ARE always one after another (all entries for one variant are sequential), so we don't have to
            # worry about the order of the file.
            # But we do have to collect multiple entries for one variant, and then decide which one is the most "important". So we:
            # 1. Iterate through some records until we find a record that is NOT the same ref/alt
            # 2. Decide if the severity of the current record is "worse" than the currently held_rec
            # 3. Write the record (function: final_process_record())
            # 4. Repeat steps 1 - 3 for the next record until we reach the end of the file
            held_rec_name = None
            held_rec = None
            held_severity_score = None

            # Iterate through records (step 1)
            for rec in reader_csv:

                # Skip the original header
                if rec['CHROM'].startswith('#'):
                    continue

                # Set a unique record ID
                current_rec_name = f'{rec["CHROM"]}_{rec["POS"]}_{rec["REF"]}_{rec["ALT"]}'

                if current_rec_name != held_rec_name: # If ID is not the same, write currently held record and reset (steps 3 - 4)
                    if held_rec_name != None: # Obviously, don't print if going through the first rec since there is no stored INFO yet
                        # Write the record with the most severe consequence (step 3)
                        held_rec = self._final_process_record(held_rec, held_severity_score, annotation_names)
                        writer_csv.writerow(held_rec)

                    # Reset to a new record (step 4)
                    held_rec_name = current_rec_name
                    held_rec = rec  # Make sure at least one record makes it through
                    # This function decides how "severe" a given CSQ annotation for a record is. See the function for
                    # more details
                    held_severity_score = self._define_score(held_rec['CSQ'])
                else:
                    # Calculate severity of this record
                    current_severity_score = self._define_score(rec['CSQ'])

                    # check to see if we should prioritise the new record based on the following ordered criteria (step 2):
                    # Below are named in DECREASING selection importance
                    # 1. protein_coding transcript
                    # 2. MANE Transcript
                    # 3. VEP Canonical Transcript
                    # 4. CSQ Severity
                    # All "else" statements are when two records have identical annotations for the given category above
                    if rec['BIOTYPE'] == 'protein_coding' and held_rec['BIOTYPE'] != 'protein_coding':
                        held_rec = rec
                        held_severity_score = self._define_score(held_rec['CSQ'])
                    elif rec['BIOTYPE'] != 'protein_coding' and held_rec['BIOTYPE'] == 'protein_coding':
                        held_rec = held_rec
                    else:
                        if (rec['MANE'] != '.' and current_severity_score['score'] <= 18) and held_rec['MANE'] == '.':
                            held_rec = rec
                            held_severity_score = self._define_score(held_rec['CSQ'])
                        elif rec['MANE'] == '.' and (held_rec['MANE'] != '.' and held_severity_score['score'] <= 18):
                            # This doesn't actually do anything, just to keep things obvious / Python happy
                            held_rec = held_rec
                        else:
                            if rec['CANONICAL'] == 'YES' and held_rec['CANONICAL'] == '.':
                                held_rec = rec
                                held_severity_score = self._define_score(held_rec['CSQ'])
                            elif rec['CANONICAL'] == '.' and held_rec['CANONICAL'] == 'YES':
                                # This doesn't actually do anything, just to keep things obvious / Python happy
                                held_rec = held_rec
                            else:
                                if current_severity_score['score'] < held_severity_score['score']:
                                    held_rec = rec
                                    held_severity_score = self._define_score(held_rec['CSQ'])

            # And print the last record since it cannot be compared to the next record:
            held_rec = self._final_process_record(held_rec, held_severity_score, annotation_names)
            writer_csv.writerow(held_rec)

        vep_gz, vep_gz_idx = bgzip_and_tabix(annote_file, comment_char='C', end_row=2)

        return vep_gz, vep_gz_idx
