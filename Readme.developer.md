# mrcepid-filterbcf Developer Readme

## Testing

Tests can be run using data in the `test_data` directory. Note that if you run tests with `KEEP_TEMP=False` then any data generated
will be kept in a temporary directory by pytest (pytest keeps data from the past 3 tests runs, and deletes data that is older). If
you run the tests with `KEEP_TEMP=True` then the newly generated data will appear in your project. This would be useful if you need to
use the outputs in any downstream analyses. 

If you wish to re-generate the test data, please use the following commands:

#### Test data

```
bcftools view -Oz https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr7.recalibrated_variants.vcf.gz "chr7:100679507-100694250" > test_input1.vcf.gz
bcftools view -Oz https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr13.recalibrated_variants.vcf.gz "chr13:36432495-36442870" > test_input2.vcf.gz
bcftools index -t test_input1.vcf.gz
bcftools index -t test_input2.vcf.gz
```

#### External files

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

## Running this app with additional computational resources

This app has the following entry points:

* main

When running this app, you can override the instance type to be used by
providing the ``systemRequirements`` field to ```/applet-XXXX/run``` or
```/app-XXXX/run```, as follows:

    {
      systemRequirements: {
        "main": {"instanceType": "mem2_hdd2_x2"}
      },
      [...]
    }

See <a
href="https://documentation.dnanexus.com/developer/api/running-analyses/io-and-run-specifications#run-specification">Run
Specification</a> in the API documentation for more information about the
available instance types.

