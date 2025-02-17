"""
NOTE: running these tests requires quite a bit of setup, namely downloading the correct files
and ensuring correct folder structure. Please refer to the developers README for a folder schema that will work with
these tests.

In brief, you will need to download the following files:
* LOFTEE hg38 files
* VEP cache files (homo_sapiens_vep_108_GRCh38)
* 1000G reference fasta

These should live in the following directories, respectively:
* test/loftee_files/loftee_hg38
* test/vep_caches/homo_sapiens
* test/reference.fasta
* test/reference.fasta.fai

Additional note: these tests take a little while (~1 minute) as we are dealing with very large files.
If you want to keep the output of these tests, change the flag KEEP_TEMP to True.

The recommendation is to run tests for the whole file, to ensure the flow of data.
"""
import os
import shutil
from pathlib import Path

import pandas as pd
import pysam
import pytest
from general_utilities.job_management.command_executor import DockerMount, CommandExecutor

from filterbcf.methods.vcf_annotate import VCFAnnotate
from filterbcf.methods.vcf_filter import VCFFilter

test_data_dir = Path(__file__).parent

# Set this flag to True if you want to keep (copy) the temporary output files
KEEP_TEMP = False

# ensure the necessary test files exist
assert Path("test/loftee_files/loftee_hg38/").exists
assert Path("test/vep_caches/homo_sapiens/").exists
assert Path("test/reference.fasta").exists
assert Path("test/reference.fasta.fai").exists

@pytest.fixture
def temporary_path(tmp_path, monkeypatch):
    """
    Prepare a temporary working directory by copying everything from the project root
    (i.e. the test directory where this file lives) into the temporary directory.

    This way, every file and folder (including reference.fasta, vep_caches, loftee_files,
    test_data, etc.) is available in the temporary directory. The working directory is
    then changed to the temporary directory.

    If KEEP_TEMP is True, after the test the entire temporary directory is copied to a
    folder 'temp_test_outputs' in the project root.
    """
    # Define the project root (where this file lives)
    project_root = Path(__file__).parent

    # Iterate over every item in the project root and copy it to the temporary directory.
    for item in project_root.iterdir():
        dest = tmp_path / item.name
        if item.is_dir():
            shutil.copytree(item, dest)
        else:
            shutil.copy2(item, dest)

    # Change the current working directory to the temporary directory.
    monkeypatch.chdir(tmp_path)

    # Yield the temporary directory to the test.
    yield tmp_path

    # After the test, if KEEP_TEMP is True, copy the temporary directory to a persistent location.
    if KEEP_TEMP:
        persistent_dir = project_root / "temp_test_outputs" / tmp_path.name
        persistent_dir.parent.mkdir(exist_ok=True)
        shutil.copytree(tmp_path, persistent_dir, dirs_exist_ok=True)
        print(f"Temporary output files have been copied to: {persistent_dir}")


@pytest.mark.parametrize(
    "vcf_filename, expected_chrom, expected_start, expected_end, expected_vep,"
    "final_vep_file",
    [
        (
                Path('/test_data/test_input1.vcf.gz'), 'chr7', 100679512, 100694238,
                Path('test_input1.vcf.sites.vcf.gz'),
                Path('test_input1.vcf.vep_table.tsv'),
        ),
        (
                Path('/test_data/test_input2.vcf.gz'), 'chr13', 36432507, 36442739,
                Path('test_input2.vcf.sites.vcf.gz'),
                Path('test_input2.vcf.vep_table.tsv'),
        ),
    ]
)
def test_vcf_annotator(temporary_path: Path, vcf_filename: Path, expected_chrom: str, expected_start: int, expected_end: int,
                       expected_vep: Path, final_vep_file: Path):
    """
    Test the VCFAnnotate class and its methods.

    Note: I am testing all the functions within the test_vcf_annotator function. This
    is because calling the VCFAnnotate class takes time, so this saves having to
    re-do it again and again. Calling it as a fixture also seems to have this
    repetitiveness.

    This test function performs the following:
    1. Loads the VCFFilter and VCFAnnotate classes.
    2. Tests the _define_score method of VCFAnnotate.
    3. Tests the _final_process_record method of VCFAnnotate.
    4. Tests the _generate_sites_file method of VCFAnnotate.
    5. Tests the _get_bcf_information method of VCFAnnotate.
    6. Tests the _run_vep method of VCFAnnotate.
    7. Tests the _parse_vep method of VCFAnnotate.
    8. Deletes test files after testing (or comment this out to keep the output).

    :param vcf_filename: The filename of the VCF file to be tested.
    :type vcf_filename: str
    """

    # load the class for testing
    test_mount = DockerMount(Path(os.getcwd()), Path('/test/'))
    cmd_exec = CommandExecutor(docker_image='egardner413/mrcepid-burdentesting', docker_mounts=[test_mount])

    class_loaded = VCFFilter(vcf_filename, cmd_exec, gq=20, wes=True, testing=True)

    vcf_annotator = VCFAnnotate(vcf_filename, class_loaded.filtered_vcf, additional_annotations=[],
                                cmd_executor=cmd_exec, dna_nexus_run=False)
    print("VCF Annotator class loaded successfully")

    # test for defining the VEP scores
    csq_expected_pairs = [
        ("stop_gained", {'score': 1, 'type': 'PTV'}),
        ("frameshift_variant", {'score': 2, 'type': 'PTV'}),
        ("splice_acceptor_variant", {'score': 3, 'type': 'PTV'}),
        ("splice_donor_variant", {'score': 4, 'type': 'PTV'}),
        ("stop_lost", {'score': 5, 'type': 'STOP_LOST'}),
        ("start_lost", {'score': 6, 'type': 'START_LOST'}),
        ("inframe_insertion", {'score': 7, 'type': 'INFRAME'}),
        ("inframe_deletion", {'score': 8, 'type': 'INFRAME'}),
        ("missense_variant", {'score': 9, 'type': 'MISSENSE'}),
        ("protein_altering_variant", {'score': 10, 'type': 'INFRAME'}),
        ("splice_region_variant", {'score': 11, 'type': 'NONCODING'}),
        ("incomplete_terminal_codon_variant", {'score': 12, 'type': 'INFRAME'}),
        ("start_retained_variant", {'score': 13, 'type': 'SYN'}),
        ("stop_retained_variant", {'score': 14, 'type': 'SYN'}),
        ("synonymous_variant", {'score': 15, 'type': 'SYN'}),
        ("5_prime_UTR_variant", {'score': 16, 'type': 'UTR'}),
        ("3_prime_UTR_variant", {'score': 17, 'type': 'UTR'}),
        ("intron_variant", {'score': 18, 'type': 'INTRONIC'}),
        ("intergenic_variant", {'score': 19, 'type': 'INTERGENIC'}),
        ("upstream_gene_variant", {'score': 20, 'type': 'INTERGENIC'}),
        ("downstream_gene_variant", {'score': 21, 'type': 'INTERGENIC'}),
        ("no_score", {'score': 22, 'type': 'ERROR'})
    ]

    for csq, expected in csq_expected_pairs:
        assert vcf_annotator._define_score(csq) == expected
    print("VEP scores defined correctly")

    # Test final_process_record
    rec = {
        'CHROM': 'chr13',
        'POS': '36432507',
        'REF': 'C',
        'ALT': 'T',
        'ID': 'chr13_36432507_C_T',
        'FILTER': 'PASS',
        'AF': '0.000626763',
        'F_MISSING': '0.00343535',
        'AN': '6382',
        'AC': '4',
        'MA': '.',
        'MANE': 'NM_003914.4',
        'ENST': 'ENST00000255465',
        'ENSG': 'ENSG00000133101',
        'BIOTYPE': 'protein_coding',
        'CANONICAL': 'YES',
        'SYMBOL': 'CCNA1',
        'CSQ': 'stop_gained',
        'SIFT': '.',
        'POLYPHEN': '.',
        'LOFTEE': '.',
        'AA': '.',
        'AApos': '.',
        'GTM': '11',
        'GT0': '3187',
        'GT1': '4',
        'GT2': '0'
    }

    expected = {
        'CHROM': 'chr13',
        'POS': '36432507',
        'REF': 'C',
        'ALT': 'T',
        'ID': 'chr13_36432507_C_T',
        'FILTER': 'PASS',
        'AF': 0.000626763,
        'F_MISSING': '0.00343535',
        'AN': '6382',
        'AC': '4',
        'MA': '.',
        'MANE': 'NM_003914.4',
        'ENST': 'ENST00000255465',
        'ENSG': 'ENSG00000133101',
        'BIOTYPE': 'protein_coding',
        'CANONICAL': 'YES',
        'SYMBOL': 'CCNA1',
        'CSQ': 'stop_gained',
        'SIFT': 'NA',
        'POLYPHEN': 'NA',
        'LOFTEE': 'NA',
        'AA': 'NA',
        'AApos': 'NA',
        'GTM': '11',
        'GT0': '3187',
        'GT1': '4',
        'GT2': '0',
        'PARSED_CSQ': 'PTV',
        'IS_MULTIALLELIC': False,
        'IS_INDEL': False,
        'MINOR_ALLELE': 'T',
        'MAJOR_ALLELE': 'C',
        'MAF': '0.000626763',
        'MAC': '4'
    }

    severity = {'score': 1, 'type': 'PTV'}
    annotation_names = ['SIFT', 'POLYPHEN', 'LOFTEE', 'AA', 'AApos']

    result = vcf_annotator._final_process_record(rec, severity, annotation_names)
    assert result == expected
    print("Record processing working as expected")

    # test for generating the sites file
    vcf_annotator._generate_sites_file(vcf_filename)
    assert Path(f'{vcf_annotator.vcfprefix}.sites.vcf.gz').exists()
    print("Sites file exists")

    # test for getting bcf information
    # Call the method
    chrom, start, end = vcf_annotator._get_bcf_information(Path(f'{vcf_annotator.vcfprefix}.sites.vcf.gz'))

    # Assertions
    assert chrom == expected_chrom
    assert start == expected_start
    assert end == expected_end
    print("BCF information parsed correctly")

    # test running VEP site file generation
    # make sure the file exists
    output_vep = Path(f'{vcf_annotator.vcfprefix}.sites.vcf.gz')
    assert output_vep.exists()
    # read in the expected file
    expected_vep = test_data_dir / 'test_data/expected_outputs' / expected_vep
    vcf1 = pysam.VariantFile(expected_vep)
    # read in the file that we just generated
    vcf2 = pysam.VariantFile(output_vep)
    # make sure the file records are the same
    for rec1, rec2 in zip(vcf1.fetch(), vcf2.fetch()):
        if rec1 != rec2:
            print(f"Record 1: {rec1}")
            print(f"Record 2: {rec2}")
            assert rec1 == rec2, f"Records do not match: {rec1} != {rec2}"
    assert sum(1 for _ in vcf1.fetch()) == sum(1 for _ in vcf2.fetch()), "Number of records do not match"
    print("VEP summary generated")

    # test parsing the final VEP file
    # make sure it exists
    assert Path(f'{vcf_annotator.vcfprefix}.vep_table.tsv').exists()
    # make sure the final VEP files are being properly created
    # read in the expected file
    final_vep_file = test_data_dir / 'test_data/expected_outputs' / final_vep_file
    df1 = pd.read_csv(final_vep_file, sep='\t')
    # read in the file we just created
    df2 = pd.read_csv(Path(f'{vcf_annotator.vcfprefix}.vep_table.tsv'), sep='\t')
    # make sure they are the same
    pd.testing.assert_frame_equal(df1, df2, check_like=True)
    print("VEP TSV files have been created successfully")
