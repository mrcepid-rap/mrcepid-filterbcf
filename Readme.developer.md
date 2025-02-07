# Developer Readme for mrcepid-filterbcf 

## Testing

Note for development: this applet has now been set to work with local unit tests, to ensure
robust CI/CD.

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
    - loftee_hg38.tar.gz
  - vep_caches/
    - homo_sapiens/
    - homo_sapiens_vep_108_GRCh38.tar.gz
```