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
    """

    def __init__(self, vcf_path: Path, cmd_executor: CommandExecutor):

        self._cmd_executor = cmd_executor
        self._vcf_prefix = vcf_path.stem

        filter_out = self._genotype_filter(vcf_path)
        flag_out = self._set_missingness_values(filter_out)
        self.filtered_vcf = self._set_filter_flags(flag_out)
        # Final file should have a name like:
        # self._vcf_prefix.missingness_filtered.bcf

    def _genotype_filter(self, input_vcf: Path) -> Path:
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
        :return: Path to the filtered bcf file
        """

        output_vcf = Path(f'{self._vcf_prefix}.filtered.bcf')
        cmd = f'bcftools filter -Ob -o /test/{output_vcf} --threads 2 -S . ' \
              f'-i "(TYPE=\'snp\' & sSUM(FMT/LAD) >= 7 & (' \
              f'(FMT/GT=\'RR\' & FMT/GQ >= 20) | ' \
              f'(FMT/GT=\'RA\' & FMT/GQ >= 20 & binom(FMT/LAD) > 0.001) | ' \
              f'(FMT/GT=\'AA\'))) | ' \
              f'(TYPE=\'indel\' & sSUM(FMT/LAD) >= 10 & FMT/GQ >= 20)" ' \
              f'/test/{input_vcf}'

        # cmd = f'bcftools filter -Ob -o /test/{output_vcf} --threads 2 -S . ' \
        #       f'-i "(TYPE=\'snp\' & FMT/DP >= 7 & (' \
        #       f'(FMT/GT=\'RR\' & FMT/GQ >= 20) | ' \
        #       f'(FMT/GT=\'RA\' & FMT/GQ >= 20 & binom(FMT/AD) > 0.001) | ' \
        #       f'(FMT/GT=\'AA\'))) | ' \
        #       f'(TYPE=\'indel\' & FMT/DP >= 10 & FMT/GQ >= 20)" ' \
        #       f'/test/{input_vcf}'
        self._cmd_executor.run_cmd_on_docker(cmd)

        input_vcf.unlink()  # Purge the original file to save HDD space
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
        cmd = f'bcftools +fill-tags /test/{input_vcf} -Ob ' \
              f'-o /test/{output_vcf} -- -t F_MISSING,AC,AF,AN'
        self._cmd_executor.run_cmd_on_docker(cmd)
        input_vcf.unlink()

        return output_vcf

    def _set_filter_flags(self, input_vcf: Path) -> Path:
        """Set pass/fail filters within the filtered VCF

        BCFTools flags used in this function:
        
        -s : sets SITES that fail the filtering expression from -i are set to FAIL
        -i : only include sites as PASS if they meet these requirements
            F_MISSING : Only include sits with less than 50% missing genotypes
            AC        : Only include biallelic sites (i.e. no monomorphic)

        :param input_vcf: Path to the input vcf file. This file will be deleted on completion of this method.
        :return: Path to the missingness filtered bcf file
        """
        output_vcf = Path(f'{self._vcf_prefix}.missingness_filtered.bcf')
        cmd = f'bcftools filter -i \'F_MISSING<=0.50 & AC!=0\' -s \'FAIL\' -Ob --threads 2 ' \
              f'-o /test/{output_vcf} ' \
              f'/test/{input_vcf}'
        self._cmd_executor.run_cmd_on_docker(cmd)
        input_vcf.unlink()

        return output_vcf

    def _write_index(self, input_vcf: Path) -> None:

        cmd = f'bcftools index /test/{input_vcf}'
        self._cmd_executor.run_cmd_on_docker(cmd)
