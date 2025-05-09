# Developer Readme for mrcepid-filterbcf 

## Testing

Tests can be run using data in the `test_data` directory. Note that if you run tests with `KEEP_TEMP=False` then any data generated
will be kept in a temporary directory by pytest (pytest keeps data from the past 3 tests runs, and deletes data that is older). If
you run the tests with `KEEP_TEMP=True` then the newly generated data will appear in your project. This would be useful if you need to
use the outputs in any downstream analyses. 

If you wish to re-generate the test data, please use the following commands:

### Test data

```
bcftools view -Oz https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr7.recalibrated_variants.vcf.gz "chr7:100679507-100694250" > test_input1.vcf.gz
bcftools view -Oz https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr13.recalibrated_variants.vcf.gz "chr13:36432495-36442870" > test_input2.vcf.gz
bcftools index -t test_input1.vcf.gz
bcftools index -t test_input2.vcf.gz
```

### External files

You also need to follow these commands to download the `reference.fasta` file:

```
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai
```

Rename the files to `reference.fasta` and `reference.fasta.fai` and make sure they are stored in the `/test/` directory.

Finally, the `test_input1.vcf.gz` and `test_input2.vcf.gz` files were altered slightly to also capture the `LAD` tag in the `FORMAT` field, for whole genome 
sequencing data. 

```
bcftools view test_input1.vcf.gz | \
    sed 's/##FORMAT=<ID=AD,/##FORMAT=<ID=LAD,/' | \
    awk -F'\t' 'BEGIN{OFS="\t"} {if ($0 ~ /^#/) print; else {gsub(":AD", ":LAD", $9); gsub(":AD", ":LAD", $10); print}}' | \
    bgzip > test_input1_wgs.vcf.gz
tabix -p vcf test_input1_wgs.vcf.gz

bcftools view test_input2.vcf.gz | \
    sed 's/##FORMAT=<ID=AD,/##FORMAT=<ID=LAD,/' | \
    awk -F'\t' 'BEGIN{OFS="\t"} {if ($0 ~ /^#/) print; else {gsub(":AD", ":LAD", $9); gsub(":AD", ":LAD", $10); print}}' | \
    bgzip > test_input2_wgs.vcf.gz
tabix -p vcf test_input2_wgs.vcf.gz
```

The files may also be missing the MA tag from the header and INFO field. If so, run the following command:

```
# Normalize with bcftools and overwrite input
for file in test_input1.vcf.gz test_input2.vcf.gz; do
    tmp_out="tmp_${file}"
    bcftools norm --threads 8 -w 100 -Oz -m - -f reference.fasta \
        --old-rec-tag MA \
        -o "$tmp_out" "$file"
    mv "$tmp_out" "$file"
done
```

Note that the dependency files (LOFTEE files, VEP files & reference files) should be downloaded
separately and arranged in the way that is outlined below, to ensure successful test runs.

One important thing to note is the folder structure. We are using a local image of Docker that expects
certain dependency files to be in a certain order. The correct folder structure should look like this:

```
- test/
  - test_data/
    - expected_outputs/
      - *all expected outputs*
    - test_input1.vcf.gz
    - test_input2.vcf.gz
  - reference.fasta
  - reference.fasta.fai
  - loftee_files/
    - loftee_hg38/
  - vep_caches/
    - homo_sapiens/
```

### Generating Test Values

We have generated a script (at `test/scripts/generate_test_values.py`) that will generate expected values following 
filtering for the test VCFs in `test/test_data`. This is mostly for sanity checking purposes, but provides a simple
way to generate filtered genotypes & missingness values for use in the `test_genotype_filter` and `test_set_id` functions
of `test/test_vcf_filter.py`. An example command line is:

```shell
python3 test/scripts/generate_test_values.py test/test_data/test_input2.vcf.gz 0.001 7 10 0.5
```

Which should give an output like:

```text
Number of ./. genotypes pre-filtering: 16050
Number of ./. genotypes post-filtering: 289806
Number of missing sites with missingness <= 0.5 and AC != 0: 67
```