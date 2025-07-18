from pathlib import Path

from general_utilities.job_management.command_executor import CommandExecutor


class VCFFilter:
    """A class that functions as a workflow to filter a vcf based on missingness.

    This class has three steps:

    1. Genotype level filtering using bcftools filter
    2. Adding missingness values and allele counts using bcftools +fill-tags
    3. Setting pass/fail filters using bcftools filter

    :param vcf_path: The original vcf file to be filtered, that has ALREADY been normalised.
    :param cmd_executor: An instance of the CommandExecutor class to run commands on the docker container.
    :param gq: Genotype quality filter to use.
    :param ad_binom: Binomial test p-value for filtering heterozygous genotypes.
    :param snp_depth: Depth filter for snps.
    :param indel_depth: Depth filter for indels.
    :param missingness: Missingness filter to use.
    :param wes: A boolean flag to indicate if the vcf is from whole exome sequencing data.
    :param testing: A boolean flag to indicate if the class is being run as part of a test.
    """

    def __init__(self, vcf_path: Path, cmd_executor: CommandExecutor, gq: int, ad_binom: float, snp_depth: int, indel_depth: int,
                 missingness: float, wes: bool, testing: bool = False):

        self._cmd_executor = cmd_executor
        self._vcf_prefix = vcf_path.stem

        self.testing = testing
        self._files_to_close = []

        filter_out = self._genotype_filter(vcf_path, gq, ad_binom, snp_depth, indel_depth, wes)
        flag_out = self._set_missingness_values(filter_out)
        id_out = self._set_id(flag_out)
        self.filtered_vcf = self._set_filter_flags(id_out, missingness)
        self.filtered_idx = self._write_index(self.filtered_vcf)
        # Final file should have a name like:
        # self._vcf_prefix.missingness_filtered.bcf

        self.close()

    def close(self) -> None:
        """
        If running the applet, close the class and delete any temporary files that were created
        If running a test, do nothing (the temporary directory will be deleted once the tests are complete) - see
        pytest documentation for more information on this.
        """
        if self.testing is False:
            for file in self._files_to_close:
                file.unlink()

    def _genotype_filter(self, input_vcf: Path, gq: int, ad_binom: float, snp_depth: int, indel_depth: int, wes: bool) -> Path:
        """Do genotype level filtering using bcftools filter.

        BCFTools flags used in this function:

        -S : sets genotypes which fail -i to missing (./.)
        -i : filtering expression to decide whether to set to missing or not.
            Note this is to INCLUDE pass genotypes, not EXCLUDE fail genotypes

        Filtering individual genotypes is split up based on variant type (snp / indel) and genotype due to issues in
        variant calling in the UK Biobank WES VCFs.

        For snps:
        - all genotypes are filtered on depth (DP) ≥ 7
        - homozygous REF ("RR") are only filtered on genotype quality (GQ) ≥ 20
        - heterozygous ("RA") are filtered on GQ ≥ 20 and a ref/alt balance binomial test [binom()] p.value of 1e-3
        - homozygous ALT ("AA") are only filtered on DP due to issues with GQ for hom alt alleles being improperly assigned

        For InDels:
        - All variants regardless of genotype are filtered on DP ≥ 10 and GQ ≥ 20

        :param input_vcf: Path to the input vcf file. This file will be deleted on completion of this method.
        :param gq: Genotype quality filter to use.
        :param ad_binom: Binomial test p-value for filtering heterozygous genotypes.
        :param snp_depth: Depth filter for snps.
        :param indel_depth: Depth filter for indels.
        :param wes: A boolean flag to indicate if the vcf is from whole exome sequencing data.
        :return: Path to the filtered bcf file
        """

        output_vcf = Path(f'{self._vcf_prefix}.filtered.bcf')
        if wes:
            filtering_string = f'"(TYPE=\'snp\' & FMT/DP >= {snp_depth} & (' \
                               f'(FMT/GT=\'RR\' & FMT/GQ >= {gq}) | ' \
                               f'(FMT/GT=\'RA\' & FMT/GQ >= {gq} & binom(FMT/AD) > {ad_binom}) | ' \
                               f'(FMT/GT=\'AA\'))) | ' \
                               f'(TYPE=\'indel\' & FMT/DP >= {indel_depth} & FMT/GQ >= {gq})"'
        else:
            filtering_string = f'"(TYPE=\'snp\' & sSUM(FMT/LAD) >= {snp_depth} & (' \
                               f'(FMT/GT=\'RR\' & FMT/GQ >= {gq}) | ' \
                               f'(FMT/GT=\'RA\' & FMT/GQ >= {gq} & binom(FMT/LAD) > {ad_binom}) | ' \
                               f'(FMT/GT=\'AA\' & FMT/GQ >= {gq}))) | ' \
                               f'(TYPE=\'indel\' & sSUM(FMT/LAD) >= {indel_depth} & FMT/GQ >= {gq})"'

        cmd = f'bcftools filter -Ob -o /test/{output_vcf.name} --threads 4 -S . -i {filtering_string} /test/{input_vcf.name}'

        self._cmd_executor.run_cmd_on_docker(cmd)

        # append a list of files to close (to save space)
        self._files_to_close.append(input_vcf)

        return output_vcf

    def _set_missingness_values(self, input_vcf: Path) -> Path:
        """Add values for missingness and AC/AF

        +fill-tags is a bcftools nonstandard plugin which calculates INFO fields from available genotypes. Note
        the non-standard '--' at then end which provides commands to the plugin rather than bcftools itself
        INFO fields added are:

        F_MISSING: fraction of missing genotypes
        AC       : allele count of ALT allele
        AF       : allele frequency of ALT allele
        AN       : number of possible alleles (accounting for missingness)

        :param input_vcf: Path to the input vcf file. This file will be deleted on completion of this method.
        :return: Path to the missingness flagged bcf file
        """

        # Do not change the command line order here. This plugin has a very specific command-line.
        output_vcf = Path(f'{self._vcf_prefix}.tagged.bcf')
        cmd = f'bcftools +fill-tags /test/{input_vcf.name} -Ob ' \
              f'-o /test/{output_vcf.name} -- -t \'F_MISSING,AC,AF,AN,GTM=count(FORMAT/GT == "./."),GT0=count(FORMAT/GT == "0/0"),GT1=count(FORMAT/GT == "0/1"),GT2=count(FORMAT/GT == "1/1")\''
        self._cmd_executor.run_cmd_on_docker(cmd)

        # append a list of files to close (to save space)
        self._files_to_close.append(input_vcf)

        return output_vcf

    def _set_id(self, input_vcf: Path) -> Path:
        """Set ID of all variants to a predictable format for downstream use.

        Method is a wrapper for bcftools annotate -I. Sets all IDs to 'CHROM_POS_REF_ALT'

        :param input_vcf: Path to the input vcf file. This file will be deleted on completion of this method.
        :return: Path to a VCF with set ID
        """
        output_vcf = Path(f'{self._vcf_prefix}.id_fixed.bcf')
        cmd = f'bcftools annotate -I "%CHROM\_%POS\_%REF\_%ALT" -Ob --threads 4 ' \
              f'-o /test/{output_vcf.name} ' \
              f'/test/{input_vcf.name}'
        self._cmd_executor.run_cmd_on_docker(cmd)

        # append a list of files to close (to save space)
        self._files_to_close.append(input_vcf)

        return output_vcf

    def _set_filter_flags(self, input_vcf: Path, missingness: float) -> Path:
        """Set pass/fail filters within the filtered VCF

        BCFTools flags used in this function:

        -s : sets SITES that fail the filtering expression from -i are set to FAIL
        -i : only include sites as PASS if they meet these requirements
            F_MISSING : Only include sits with less than 50% missing genotypes
            AC        : Only include biallelic sites (i.e. no monomorphic)

        :param input_vcf: Path to the input vcf file. This file will be deleted on completion of this method.
        :param missingness: Missingness filter to use.
        :return: Path to the missingness filtered bcf file
        """
        output_vcf = Path(f'{self._vcf_prefix}.missingness_filtered.bcf')
        cmd = f'bcftools filter -i \'F_MISSING<={missingness} & AC!=0\' -s \'FAIL\' -Ob --threads 4 ' \
              f'-o /test/{output_vcf.name} ' \
              f'/test/{input_vcf.name}'
        self._cmd_executor.run_cmd_on_docker(cmd)

        # append a list of files to close (to save space)
        self._files_to_close.append(input_vcf)

        return output_vcf

    def _write_index(self, input_vcf: Path) -> Path:
        """Generate a .csi index file for the final bcf file.

        :param input_vcf: The final, filtered bcf file to be indexed
        :return: The path to the indexed bcf file (in .csi format)
        """

        output_index = Path(f'{input_vcf}.csi')
        cmd = f'bcftools index --threads 4 /test/{input_vcf.name}'
        self._cmd_executor.run_cmd_on_docker(cmd)

        return output_index
