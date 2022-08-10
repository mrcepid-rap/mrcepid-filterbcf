import csv

from typing import TypedDict, Tuple
from ..filterbcf_resources import *


class VCFAnnotate:

    def __init__(self, vcfprefix: str):

        self.vcfprefix = vcfprefix

        # 1. Generate sites files that can be run through the various annotation parts
        self._generate_sites_files()

        # 1a. Get information about this VCF from the sites file since it is much quicker:
        # chromosome, start, end
        (self.chunk_chrom, self.chunk_start, self.chunk_end) = self._get_bcf_information()

        # 2. This does the initial VEP run
        self._run_vep()

        # 3. Add annotations that need to be done separately from VEP (gnomAD & CADD)
        self._add_gnomad_annotation()
        self._add_cadd_annotation()

        # 4. Generate an annotations TSV for each VCF
        self._generate_annotations_tsv()

        # 5. Generating merged/final files:
        # This function parses the information from the raw VEP run (via run_vep()) and adds it to our filtered vcf
        self._parse_vep()

        # 6. Annotate the final, filtered VCF with our new VEP fields to make it easy to go back and generate files
        # for association testing
        self._annotate_vcf_with_vep()

        # 7. Generate a file of states for QC purposes to make sure this didn't go south:
        self._generate_qc_file()

        # 8. Set final file outputs for this entire process:
        self.finalbcf = generate_linked_dx_file(self.vcfprefix + ".filtered_annotated.bcf")
        self.finalbcf_index = generate_linked_dx_file(self.vcfprefix + ".filtered_annotated.bcf.csi")
        self.finalvep = generate_linked_dx_file(self.vcfprefix + ".vep.tsv.gz")
        self.finalvep_index = generate_linked_dx_file(self.vcfprefix + ".vep.tsv.gz.tbi")
        self.output_per_sample = generate_linked_dx_file(self.vcfprefix + '.per_indv.tsv')

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
    def _final_process_record(rec: dict, severity: ConsequenceSeverity) -> dict:

        # Has to have a "#" to be compatible with VCF I/O
        rec['#CHROM'] = rec['CHROM']

        # Rename some columns for printing purposes
        rec['parsed_csq'] = severity['type']

        # Records who do not have equal REF/ALT length are assigned as InDels
        if len(rec['REF']) != len(rec['ALT']):
            rec['is_indel'] = True
        else:
            rec['is_indel'] = False

        # By UKBB convention, only variants that are multiallelic have a ";" in the name
        if ';' in rec['ID']:
            rec['is_multiallelic'] = True
        else:
            rec['is_multiallelic'] = False

        # This corrects an issue with sites with 100% missingness that BCFtools doesn't handle correctly
        if rec['AF'] == '.':
            rec['AF'] = 0
            rec['FILTER'] = 'FAIL'
        else:
            rec['AF'] = float(rec['AF'])

        # Set sensible defaults for a variety of fields:
        # "gnomad_maf", "REVEL" "SIFT", "PolyPhen", "LoF"
        rec['gnomad_maf'] = rec['gnomad_maf'] if rec['gnomad_maf'] != '.' else '0' # Sites w/o gnomAD don't exist in gnomAD so a MAF of 0 seems appropriate
        rec['REVEL'] = rec['REVEL'] if rec['REVEL'] != '.' else 'NaN' # NaN is default VCF spec for missing floats
        rec['SIFT'] = rec['SIFT'] if rec['SIFT'] != '.' else 'NA' # NA is default VCF spec for missing strings
        rec['PolyPhen'] = rec['PolyPhen'] if rec['PolyPhen'] != '.' else 'NA' # NA is default VCF spec for missing strings
        rec['LoF'] = rec['LoF'] if rec['LoF'] != '.' else 'NA' # NA is default VCF spec for missing strings
        rec['AA'] = rec['AA'] if rec['AA'] != '.' else 'NA'  # NA is default VCF spec for missing strings
        rec['AApos'] = rec['AApos'] if rec['AApos'] != '.' else 'NA'  # NA is default VCF spec for missing strings

        # Setting additional tags requested by GWAS-y people
        if float(rec['AF']) < 0.5:
            rec['minor_allele'] = rec['ALT']
            rec['major_allele'] = rec['REF']
            rec['MAF'] = '%s' % rec['AF']
            rec['MAC'] = '%s' % rec['AC']
        else:
            rec['minor_allele'] = rec['REF']
            rec['major_allele'] = rec['ALT']
            rec['MAF'] = '%s' % (1 - float(rec['AF']))
            rec['MAC'] = '%s' % (int(rec['AN']) - int(rec['AC']))

        return rec

    # This function generates sites files for both VEP and CADD, which require slightly different formats
    def _generate_sites_files(self):

        # Generate a file without genotypes for VEP
        # This should be the ONLY file that is .vcf.gz format as it is required by VEP
        # -G : strips genotypes
        cmd = "bcftools view --threads 2 -G -Oz " \
              "-o /test/" + self.vcfprefix + ".norm.filtered.tagged.missingness_filtered.sites.vcf.gz /test/" + \
              self.vcfprefix + ".norm.filtered.tagged.missingness_filtered.bcf"
        run_cmd(cmd, True)

        # Generate a sites file that is in the correct format for CADD from the above
        cmd = "bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\n' " \
              "-o /test/" + self.vcfprefix + ".norm.filtered.tagged.missingness_filtered.sites.cadd.vcf /test/" + \
              self.vcfprefix + ".norm.filtered.tagged.missingness_filtered.sites.vcf.gz"
        run_cmd(cmd, True)

        # CADD doesn't like the 'chr' prefix..., so remove it!
        cmd = "sed -i \'s_chr__\' " + self.vcfprefix + ".norm.filtered.tagged.missingness_filtered.sites.cadd.vcf"
        run_cmd(cmd)

    def _get_bcf_information(self) -> Tuple[str, int, int]:

        sites_reader = open(self.vcfprefix + '.norm.filtered.tagged.missingness_filtered.sites.cadd.vcf', 'r')
        linenum = 1
        for site in sites_reader:
            site = site.rstrip()
            site_data = site.split('\t')
            if linenum == 1:
                chunk_chrom = 'chr' + site_data[0]  # We are using the stripped CADD file - have to add chr back in
                chunk_start = int(site_data[1])
            chunk_end = int(site_data[1])
            linenum += 1

        return chunk_chrom, chunk_start, chunk_end

    # This function runs VEP
    def _run_vep(self) -> None:

        # Run VEP on the sites file generated by _generate_sites_files():
        # Not going to document each individual thing for VEP, but all are available on the VEP webiste
        cmd = "perl -Iensembl-vep/cache/Plugins/loftee/ -Iensembl-vep/cache/Plugins/loftee/maxEntScan/ " \
              "ensembl-vep/vep --offline --cache --assembly GRCh38 --dir_cache /test/vep_caches/ --everything --allele_num " \
              "-i /test/" + self.vcfprefix + ".norm.filtered.tagged.missingness_filtered.sites.vcf.gz --format vcf --fasta /test/reference.fasta " \
              "-o /test/" + self.vcfprefix + ".norm.filtered.tagged.missingness_filtered.sites.vep.vcf.gz --compress_output bgzip --vcf " \
              "--dir_plugins ensembl-vep/cache/Plugins/ " \
              "--plugin LoF,loftee_path:ensembl-vep/cache/Plugins/loftee,human_ancestor_fa:/test/loftee_files/loftee_hg38/human_ancestor.fa.gz,conservation_file:/test/loftee_files/loftee_hg38/loftee.sql,gerp_bigwig:/test/loftee_files/loftee_hg38/gerp_conservation_scores.homo_sapiens.GRCh38.bw " \
              "--plugin REVEL,/test/revel_files/new_tabbed_revel_grch38.tsv.gz"
        run_cmd(cmd, True)
        purge_file(self.vcfprefix + ".norm.filtered.tagged.missingness_filtered.sites.vcf.gz")

    # Do gnomAD annotation
    def _add_gnomad_annotation(self):

        # Add better gnomAD MAF information
        # bcftools just takes a tabix-formated tsv file (in this case gnomad.tsv.gz) and adds non-coordinate fields as INFO fields
        # In this case, the gnomad.tsv.gz file has 6 columns, and I ignore column 5
        # Thus, an INFO field named "gnomAD_MAF" is added to the VEP-annotated VCF
        # There is also a file (gnomad.header.txt) generated in ingest_data() that contains the INFO field annotation.
        cmd = "bcftools annotate --threads 2 -a /test/gnomad_files/gnomad.tsv.gz -c CHROM,POS,REF,ALT,-,gnomAD_MAF -Ob " \
              "-h /test/gnomad_files/gnomad.header.txt " \
              "-o /test/" + self.vcfprefix + ".norm.filtered.tagged.missingness_filtered.sites.vep.gnomad.bcf /test/" + \
              self.vcfprefix + ".norm.filtered.tagged.missingness_filtered.sites.vep.vcf.gz"
        run_cmd(cmd, True)
        purge_file(self.vcfprefix + ".norm.filtered.tagged.missingness_filtered.sites.vep.vcf.gz")

    # Do CADD annotation
    def _add_cadd_annotation(self):

        # Run CADD on the sites file from generate_sites_files():
        cmd = "CADD-scripts/CADD.sh -c 2 -g GRCh38 " \
              "-o /test/" + self.vcfprefix + ".cadd.tsv.gz /test/" + \
              self.vcfprefix + ".norm.filtered.tagged.missingness_filtered.sites.cadd.vcf"
        run_cmd(cmd, True)

        # Add chr back so BCFtools can understand for reannotation and then bgzip and tabix index
        cmd = "zcat " + self.vcfprefix + ".cadd.tsv.gz | tail -n+3 | sed \'s_^_chr_\' > " + self.vcfprefix + ".cadd.chr.tsv"
        run_cmd(cmd)
        cmd = "bgzip /test/" + self.vcfprefix + ".cadd.chr.tsv"
        run_cmd(cmd, True)
        cmd = "tabix -p vcf /test/" + self.vcfprefix + ".cadd.chr.tsv.gz"
        run_cmd(cmd, True)

        # Now annotate the gnomAD VCF with CADD scores:
        cmd = "bcftools annotate --threads 2 -c CHROM,POS,REF,ALT,-,CADD -h /test/cadd.header.txt -Ob " \
              "-a /test/" + self.vcfprefix + ".cadd.chr.tsv.gz " + \
              "-o /test/" + self.vcfprefix + ".norm.filtered.tagged.missingness_filtered.sites.vep.gnomad.cadd.bcf " + \
              "/test/" + self.vcfprefix + ".norm.filtered.tagged.missingness_filtered.sites.vep.gnomad.bcf"
        run_cmd(cmd, True)

        # Remove CADD annotation files to save space
        purge_file(self.vcfprefix + ".cadd.tsv.gz")
        purge_file(self.vcfprefix + ".cadd.chr.tsv.gz")
        purge_file(self.vcfprefix + ".cadd.chr.tsv.gz.tbi")

    # This function extracts individual annotations using the +split-vep tool in bcftools
    def _generate_annotations_tsv(self):

        # Then generate an output TSV that we can parse:
        # The purpose of this is to generate a TSV file of annotations that we care about to parse later (function parse_vep())
        # +split-vep is a bcftools plugin that iterates through VEP fields provided in a VCF via the CSQ INFO field.
        # -d :  outputs duplicate transcripts on separate lines. In other words, a gene may have multiple transcripts,
        # and put each transcript and CSQ on a separate line in the tsv file
        # -f : CSQ fields that we want to output into the tsv file
        cmd = "bcftools +split-vep -df '%CHROM\\t%POS\\t%REF\\t%ALT\\t%ID\\t%FILTER\\t%INFO/AF\\t%F_MISSING\\t%AN\\t%AC\\t%MANE_SELECT\\t%Feature\\t%Gene\\t%BIOTYPE\\t%CANONICAL\\t%SYMBOL\\t%Consequence\\t%gnomAD_MAF\\t%CADD\\t%REVEL\\t%SIFT\\t%PolyPhen\\t%LoF\\t%Amino_acids\\t%Protein_position\\n' " \
              "-o /test/" + self.vcfprefix + ".vep_table.tsv " + \
              "/test/" + self.vcfprefix + ".norm.filtered.tagged.missingness_filtered.sites.vep.gnomad.cadd.bcf"
        run_cmd(cmd, True)
        purge_file(self.vcfprefix + ".norm.filtered.tagged.missingness_filtered.sites.vep.gnomad.bcf")

    # This function parses VEP output for most severe CSQ for each variant.
    # See individual comments in this code to understand how that is done.
    def _parse_vep(self) -> None:

        # These are all possible fields from the vep table that we generated in run_vep()
        # And then read it in as a csv.DictReader()
        csv_reader_header = ("CHROM", "POS", "REF", "ALT", "ID", "FILTER", "AF", "prop_missing", "AN", "AC",
                             "mane_transcript", "ENST_ID", "ENSG_ID", "biotype", "is_canonical", "symbol", "csq",
                             "gnomad_maf", "CADD", "REVEL", "SIFT", "PolyPhen", "LoF", "AA", "AApos")
        vep_reader = csv.DictReader(open(self.vcfprefix + '.vep_table.tsv', 'r', newline='\n'), delimiter="\t",
                                    fieldnames=csv_reader_header, quoting=csv.QUOTE_NONE)

        # Next, open a file that will contain (in tabix format tsv) the info we want to add back to the vcf
        # And these are all possible output fields that we want
        annote_writer_header = ("#CHROM", "POS", "REF", "ALT", "mane_transcript", "ENST_ID", "ENSG_ID", "biotype", # OG Fields
                                "symbol", "csq", "gnomad_maf", "CADD", "REVEL", "SIFT", "PolyPhen", "LoF", "AA", "AApos", # OG Fields
                                "parsed_csq", "is_multiallelic", "is_indel", "minor_allele", "major_allele", "MAF", "MAC") # New Fields

        # And open the csv DictWriter to put annotations into
        annote_file = open(self.vcfprefix + '.vep_table.annote.tsv', 'w')
        annote_writer = csv.DictWriter(annote_file, delimiter="\t", fieldnames=annote_writer_header, extrasaction='ignore')

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
        for rec in vep_reader:

            # Set a unique record ID
            current_rec_name = '%s_%s_%s_%s' % (rec['CHROM'], rec['POS'], rec['REF'], rec['ALT'])

            if current_rec_name != held_rec_name: # If ID is not the same, write currently held record and reset (steps 3 - 4)
                if (held_rec_name != None): # Obviously, don't print if going through the first rec since there is no stored INFO yet
                    # Write the record with the most severe consequence (step 3)
                    held_rec = self._final_process_record(held_rec, held_severity_score)
                    annote_writer.writerow(held_rec)

                # Reset to a new record (step 4)
                held_rec_name = current_rec_name
                held_rec = rec # Make sure at least one record makes it through
                # This function decides how "severe" a given CSQ annotation for a record is. See the function for more details
                held_severity_score = self._define_score(held_rec['csq'])
            else:
                # Calculate severity of this record
                current_severity_score = self._define_score(rec['csq'])

                # check to see if we should prioritise the new record based on the following ordered criteria (step 2):
                # Below are named in DECREASING selection importance
                # 1. protein_coding transcript
                # 2. MANE Transcript
                # 3. VEP Canonical Transcript
                # 4. CSQ Severity
                # All "else" statements are when two records have identical annotations for the given category above
                if rec['biotype'] == 'protein_coding' and held_rec['biotype'] != 'protein_coding':
                    held_rec = rec
                    held_severity_score = self._define_score(held_rec['csq'])
                elif rec['biotype'] != 'protein_coding' and held_rec['biotype'] == 'protein_coding':
                    held_rec = held_rec
                else:
                    if (rec['mane_transcript'] != '.' and current_severity_score['score'] <= 18) and held_rec['mane_transcript'] == '.':
                        held_rec = rec
                        held_severity_score = self._define_score(held_rec['csq'])
                    elif rec['mane_transcript'] == '.' and (held_rec['mane_transcript'] != '.' and held_severity_score['score'] <= 18):
                        # This doesn't actually do anything, just to keep things obvious / Python happy
                        held_rec = held_rec
                    else:
                        if rec['is_canonical'] == 'YES' and held_rec['is_canonical'] == '.':
                            held_rec = rec
                            held_severity_score = self._define_score(held_rec['csq'])
                        elif rec['is_canonical'] == '.' and held_rec['is_canonical'] == 'YES':
                            # This doesn't actually do anything, just to keep things obvious / Python happy
                            held_rec = held_rec
                        else:
                            if current_severity_score['score'] < held_severity_score['score']:
                                held_rec = rec
                                held_severity_score = self._define_score(held_rec['csq'])

        # And print the last record since it cannot be compared to an old record above:
        held_rec = self._final_process_record(held_rec, held_severity_score)
        annote_writer.writerow(held_rec)

        # This flushes all data & closes the output since this is not done by default in csv.DictWriter
        annote_file.close()

        # bgzip/tabix the output(s) to save space on DNAnexus / allow postprocessing
        cmd = "bgzip /test/" + self.vcfprefix + ".vep_table.annote.tsv"
        run_cmd(cmd, True)
        cmd = "tabix -p vcf /test/" + self.vcfprefix + ".vep_table.annote.tsv.gz"
        run_cmd(cmd, True)

    # This function simply runs bcftools annotate to add VEP information back to our QCd VCF
    def _annotate_vcf_with_vep(self) -> None:

        # This is similar to how BCFtools annotate was run above to add gnomAD MAF but for A LOT more fields that we got via VEP
        cmd = "bcftools annotate --threads 2 -a /test/" + self.vcfprefix + ".vep_table.annote.tsv.gz -c " \
              "CHROM,POS,REF,ALT,MANE,ENST,ENSG,BIOTYPE,SYMBOL,CSQ,gnomAD_AF,CADD,REVEL,SIFT,POLYPHEN,LOFTEE,AA,AApos,PARSED_CSQ,MULTI,INDEL,MINOR,MAJOR,MAF,MAC " \
              "-h /test/vep_vcf.header.txt -Ob -o /test/" + self.vcfprefix + ".filtered_annotated.bcf /test/" + self.vcfprefix + ".norm.filtered.tagged.missingness_filtered.bcf"
        run_cmd(cmd, True)
        purge_file(self.vcfprefix + ".norm.filtered.tagged.missingness_filtered.bcf")

        # Index this file
        cmd = "bcftools index /test/" + self.vcfprefix + ".filtered_annotated.bcf"
        run_cmd(cmd, True)

        # Generate a final .tsv of annotations:
        cmd = 'bcftools query -f ' \
              '"%CHROM\\t%POS\\t%REF\\t%ALT\\t%ID\\t%FILTER\\t%INFO/AF\\t%F_MISSING\\t%AN\\t%AC\\t%MANE\\t%ENST\\t%ENSG\\t%BIOTYPE\\t' \
              '%SYMBOL\\t%CSQ\\t%gnomAD_AF\\t%CADD\\t%REVEL\\t%SIFT\\t%POLYPHEN\\t%LOFTEE\\t%AA\\t%AApos\\t%PARSED_CSQ\\t%MULTI\\t%INDEL\\t%MINOR\\t' \
              '%MAJOR\\t%MAF\\t%MAC\\n" -o /test/' + self.vcfprefix + '.vep.tsv ' + \
              '/test/' + self.vcfprefix + '.filtered_annotated.bcf'
        run_cmd(cmd, True)

        # And bgzip and tabix this file
        cmd = "bgzip /test/" + self.vcfprefix + ".vep.tsv"
        run_cmd(cmd, True)
        cmd = "tabix -p vcf /test/" + self.vcfprefix + ".vep.tsv.gz"
        run_cmd(cmd, True)

    # Just generates per-sample QC information
    def _generate_qc_file(self):

        cmd = 'bcftools query -i \'FILTER == "PASS" && INFO/AF < 0.001 && ((INFO/PARSED_CSQ == "PTV" && INFO/LOFTEE == "HC") || INFO/PARSED_CSQ == "SYN" || (INFO/PARSED_CSQ == "MISSENSE" && INFO/CADD > 25)) && GT="alt"\' ' \
              '-f "[%SAMPLE\\t%PARSED_CSQ\\t%AF\\t%GT\\n]" ' \
              '-o /test/' + self.vcfprefix + '.per_indv.tsv ' + \
              '/test/' + self.vcfprefix + '.filtered_annotated.bcf'
        run_cmd(cmd, True)

