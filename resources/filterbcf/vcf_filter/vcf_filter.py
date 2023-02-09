from pathlib import Path

from general_utilities.association_resources import run_cmd


class VCFFilter:

    def __init__(self, vcfprefix: str):

        self.vcfprefix = vcfprefix

        self._normalise_and_left_correct()
        self._genotype_filter()
        self._set_missingness_values()
        self._set_filter_flags()

    # Just a wrapper for BCFtools norm for each bcf file
    # Generate a normalised bcf file for all downstream processing:
    # -m : splits all multiallelics into separate records
    # -f : provides a reference file so bcftools can left-normalise and check records against the reference genome
    def _normalise_and_left_correct(self) -> None:
        cmd = "bcftools norm --threads 2 -Ob -m - -f /test/reference.fasta " \
              "-o /test/" + self.vcfprefix + ".norm.bcf /test/" + self.vcfprefix + ".bcf"
        run_cmd(cmd, is_docker=True, docker_image='egardner413/mrcepid-burdentesting')
        Path(self.vcfprefix + ".bcf").unlink()

    # Do genotype level filtering
    # -S : sets genotypes which fail -i to missing (./.)
    # -i : filtering expression to decide whether to set to missing or not.
    # Note this is to INCLUDE pass genotypes, not EXCLUDE fail genotypes
    # Filtering is split up based on variant type (snp / indel) and genotype due to issues in variant
    # calling in the UK Biobank WES VCFs
    # For snps:
    # - all genotypes are filtered on depth (DP) ≥ 7
    # - homozygous REF ("RR") are only filtered on genotype quality (GQ) ≥ 20
    # - heterozygous ("RA") are filtered on GQ ≥ 20 and a ref/alt balance binomial test [binom()] p.value of 1e-3
    # - homozygous ALT ("AA") are only filtered on DP due to issues with GQ for hom alt alleles being improperly assigned
    # For InDels:
    # - All variants regardless of genotype are filtered on DP ≥ 10 and GQ ≥ 20
    def _genotype_filter(self) -> None:
        cmd = "bcftools filter -Ob -o /test/" + self.vcfprefix + ".norm.filtered.bcf --threads 2 -S . " \
                "-i '(TYPE=\"snp\" & FMT/DP >= 7 & ((FMT/GT=\"RR\" & FMT/GQ >= 20) | " \
                "(FMT/GT=\"RA\" & FMT/GQ >= 20 & binom(FMT/AD) > 0.001) | " \
                "(FMT/GT=\"AA\"))) | " \
                "(TYPE=\"indel\" & FMT/DP >= 10 & FMT/GQ >= 20)' " \
                "/test/" + self.vcfprefix + ".norm.bcf"
        run_cmd(cmd, is_docker=True, docker_image='egardner413/mrcepid-burdentesting')
        Path(self.vcfprefix + ".norm.bcf").unlink()  # Purge the original file to save memory

    # Add values for missingness and AC/AF
    # +fill-tags is a bcftools nonstandard plugin which calculates INFO fields from available genotypes. Note the non-standard '--' at then end which provides
    # commands to the plugin rather than bcftools itself
    # INFO fields added are:
    # F_MISSING: fraction of missing genotypes
    # AC       : allele count of ALT allele
    # AF       : allele frequency of ALT allele
    # AN       : number of possible alleles (accounting for missingness)
    def _set_missingness_values(self):

        # Do not change the command line order here. This plugin has a very specific command-line.
        cmd = "bcftools +fill-tags /test/" + self.vcfprefix + ".norm.filtered.bcf -Ob -o /test/" + self.vcfprefix + ".norm.filtered.tagged.bcf -- -t F_MISSING,AC,AF,AN"
        run_cmd(cmd, is_docker=True, docker_image='egardner413/mrcepid-burdentesting')
        Path(self.vcfprefix + ".norm.filtered.bcf").unlink()

    # Set pass/fail filters within the filtered VCF
    # -s : sets SITES that fail the filtering expression from -i are set to FAIL
    # -i : only include sites as PASS if they meet these requirements
    # F_MISSING : Only include sits with less than 50% missing genotypes
    # AC        : Only include biallelic sites (i.e. no monomorphic)
    def _set_filter_flags(self):
        cmd = "bcftools filter -i \'F_MISSING<=0.50 & AC!=0\' -s \'FAIL\' -Ob --threads 2 " \
              "-o /test/" + self.vcfprefix + ".norm.filtered.tagged.missingness_filtered.bcf " \
              "/test/" + self.vcfprefix + ".norm.filtered.tagged.bcf"
        run_cmd(cmd, is_docker=True, docker_image='egardner413/mrcepid-burdentesting')
        Path(self.vcfprefix + ".norm.filtered.tagged.bcf").unlink()
