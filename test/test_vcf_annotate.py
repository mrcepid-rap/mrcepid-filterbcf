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
import subprocess
from pathlib import Path

import pandas as pd
import pysam
import pytest
from general_utilities.job_management.command_executor import DockerMount, CommandExecutor

from filterbcf.methods.vcf_annotate import VCFAnnotate
from filterbcf.methods.vcf_filter import VCFFilter

test_data_dir = Path(__file__).parent

# ensure the necessary test files exist
assert Path("test/loftee_files/loftee_hg38/").exists
assert Path("test/vep_caches/homo_sapiens/").exists
assert Path("test/reference.fasta").exists
assert Path("test/reference.fasta.fai").exists

# Set this flag to True if you want to keep (copy) the temporary output files
KEEP_TEMP = False  # or True if needed


@pytest.fixture
def temporary_path(tmp_path, monkeypatch):
    project_root = Path(__file__).parent

    # Skip copying `temp_test_outputs` directory to avoid recursion
    for item in project_root.iterdir():
        if item.name == "temp_test_outputs":
            continue  # skip this
        dest = tmp_path / item.name
        if item.is_dir():
            shutil.copytree(item, dest)
        else:
            shutil.copy2(item, dest)

    monkeypatch.chdir(tmp_path)

    yield tmp_path

    if KEEP_TEMP:
        persistent_dir = project_root / "temp_test_outputs" / tmp_path.name
        persistent_dir.parent.mkdir(parents=True, exist_ok=True)
        # Only copy selected outputs if needed
        try:
            shutil.copytree(tmp_path, persistent_dir, dirs_exist_ok=True)
        except OSError as e:
            print(f"Warning: could not copy temp files: {e}")


@pytest.fixture(scope="session", autouse=True)
def cleanup_temp_outputs():
    """
    Delete temp_test_outputs folder at the end of the test session (only if KEEP_TEMP is False).
    """
    yield  # Run tests first

    # given how large some of the temporary files are, we will delete them after testing
    subprocess.run('./clean_pytest_temp.sh', shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    print("Temporary files deleted successfully")


@pytest.mark.parametrize(
    "vcf_filename, expected_chrom, expected_start, expected_end, expected_vep,"
    "final_vep_file, number_of_variants, vep_table_length, vep_unique_variants_length,"
    "vep_unique_pass_variants, annotation_dict, final_vep_annotation_path,"
    "annotation_name, expected_number_of_annotations",
    [
        (
                Path('/test_data/test_input1.vcf.gz'), 'chr7', 100679512, 100694238,
                Path('test_input1.vcf.sites.vcf.gz'),
                Path('test_input1.vcf.vep_table.tsv'),
                907,
                9597,
                875,
                777,
                {
                    'file': Path('test_data/anno1.tsv.gz'),
                    'index': Path('test_data/anno1.tsv.gz.tbi'),
                    'annotation_name': 'some_annotation',
                    'header_file': "##INFO=<ID=SANNO,Number=1,Type=String,Description='SANNO annotation.'>",
                    'symbol_mode': ''
                },
                Path('test_data/expected_outputs/test_input1.vcf.vep_table.tsv'),
                'some_annotation',
                9136,
        ),
        (
                Path('/test_data/test_input2.vcf.gz'), 'chr13', 36432507, 36442739,
                Path('test_input2.vcf.sites.vcf.gz'),
                Path('test_input2.vcf.vep_table.tsv'),
                558,
                2179,
                504,
                461,
                {
                    'file': Path('test_data/anno2.tsv.gz'),
                    'index': Path('test_data/anno2.tsv.gz.tbi'),
                    'annotation_name': 'some_other_annotation',
                    'header_file': "##INFO=<ID=SMANNO,Number=1,Type=String,Description='SMANNO annotation.'>",
                    'symbol_mode': 'ENST'
                },
                Path('test_data/expected_outputs/test_input2.vcf.vep_table.tsv'),
                'some_other_annotation',
                1089,
        ),
    ]
)
def test_vcf_annotator(temporary_path: Path, vcf_filename: Path, expected_chrom: str, expected_start: int,
                       expected_end: int, expected_vep: Path, final_vep_file: Path, number_of_variants: int,
                       vep_table_length: int, vep_unique_variants_length: int, vep_unique_pass_variants: int,
                       annotation_dict: dict, final_vep_annotation_path: Path, annotation_name: str,
                       expected_number_of_annotations: int):
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
                                cmd_executor=cmd_exec)
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

    variants_df = read_vcf_as_dataframe(Path(f'{vcf_annotator.vcfprefix}.sites.vcf.gz'))
    assert len(variants_df) == number_of_variants

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
    print("VEP summary generated")

    # test parsing the final VEP file
    # make sure it exists
    assert Path(f'{vcf_annotator.vcfprefix}.vep_table.tsv').exists()
    vep_table = pd.read_csv(Path(f'{vcf_annotator.vcfprefix}.vep_table.tsv'), sep='\t')
    # assert table length
    assert len(vep_table) == vep_table_length
    # assert number of unique variants
    assert len(vep_table[['#[1]CHROM', '[2]POS', '[3]REF', '[4]ALT']].drop_duplicates()) == vep_unique_variants_length
    # assert number of unique variants that have passed the quality filter
    assert len(vep_table[vep_table['[6]FILTER'] == 'PASS'][
                   ['#[1]CHROM', '[2]POS', '[3]REF', '[4]ALT']].drop_duplicates()) == vep_unique_pass_variants

    # make sure the final VEP files are being properly created
    # read in the expected file
    final_vep_file = test_data_dir / 'test_data/expected_outputs' / final_vep_file
    df1 = pd.read_csv(final_vep_file, sep='\t')
    # read in the file we just created
    df2 = pd.read_csv(Path(f'{vcf_annotator.vcfprefix}.vep_table.tsv'), sep='\t')
    # make sure they are the same
    # Normalize Consequence column: sort &-separated terms
    df1['[18]Consequence'] = df1['[18]Consequence'].apply(
        lambda x: '&'.join(sorted(str(x).split('&'))) if pd.notna(x) else x
    )
    df2['[18]Consequence'] = df2['[18]Consequence'].apply(
        lambda x: '&'.join(sorted(str(x).split('&'))) if pd.notna(x) else x
    )
    # Compare DataFrames (column order agnostic)
    pd.testing.assert_frame_equal(df1, df2, check_like=True)
    print("VEP TSV files have been created successfully and match expected output.")

    # check for additional annotations
    vep_tsv, annotation_name = vcf_annotator._add_additional_annotation(final_vep_annotation_path,
                                                                        annotation_dict)

    # make sure the file exists
    assert Path(vep_tsv).exists()

    # make sure the new column is in the file
    df = pd.read_csv(vep_tsv, sep='\t')
    assert any(annotation_name in col for col in df.columns)

    # pull out the column indices containing the annotation_name
    annotation_indices = [i for i, col in enumerate(df.columns) if annotation_name in col]
    # make sure number of annotations is correct
    # Count number of cells in those columns that have value equal to annotation_name
    count = sum(
        (df.iloc[:, idx] == annotation_name).sum()
        for idx in annotation_indices
    )
    # Assert that this count matches your expected number
    assert count == expected_number_of_annotations, f"Expected {expected_number_of_annotations}, but got {count}"

    # make sure the annotations match
    anno = pd.read_csv(annotation_dict['file'], sep="\t")
    if annotation_dict['symbol_mode'] == 'ENST':
        # Create a set of keys from the annotation file
        anno_keys = set(zip(anno["#CHROM"], anno["POS"], anno["REF"], anno["ALT"], anno["ENST"]))
        # Create a list of keys from your main df for matching
        df_keys = list(zip(df["#[1]CHROM"], df["[2]POS"], df["[3]REF"], df["[4]ALT"], df["[13]Feature"]))
    else:
        # Create a set of keys from the annotation file
        anno_keys = set(zip(anno["#CHROM"], anno["POS"], anno["REF"], anno["ALT"]))
        # Create a list of keys from your main df for matching
        df_keys = list(zip(df["#[1]CHROM"], df["[2]POS"], df["[3]REF"], df["[4]ALT"]))
    # Generate a boolean mask indicating matching rows
    matching_rows_mask = pd.Series(df_keys).isin(anno_keys)
    # Identify the annotation columns in your df
    annotation_indices = [i for i, col in enumerate(df.columns) if annotation_name in col]
    # Extract just the annotation columns from df
    annotation_df = df.iloc[:, annotation_indices]
    # Assert: All values in matching rows must be 'some_other_annotation'
    assert (annotation_df[matching_rows_mask] == annotation_name).all().all(), \
        "Mismatch: Not all matching rows contain 'some_other_annotation' in annotation columns."
    # Assert: All values in non-matching rows must NOT be 'some_other_annotation'
    assert (annotation_df[~matching_rows_mask] != annotation_name).all().all(), \
        "Mismatch: Some non-matching rows contain 'some_other_annotation' in annotation columns."




def read_vcf_as_dataframe(vcf_file_path: Path):
    """
    A way to read VCF files as pandas dataframes

    :param vcf_file_path: filepath to a VCF file
    :return: a pandas dataframe
    """
    # Open the VCF file using pysam
    vcf_in = pysam.VariantFile(vcf_file_path)

    # Initialize a list to store the records
    records = []

    # Iterate through each record in the VCF file
    for record in vcf_in.fetch():
        # Extract the necessary information from each record
        rec = {
            'CHROM': record.chrom,
            'POS': record.pos,
            'ID': record.id,
            'REF': record.ref,
            'ALT': ','.join(str(alt) for alt in record.alts),
            'QUAL': record.qual,
            'FILTER': ';'.join(record.filter.keys()),
            'INFO': ';'.join(f'{key}={value}' for key, value in record.info.items())
        }
        records.append(rec)

    # Convert the list of records to a pandas DataFrame
    df = pd.DataFrame(records)

    return df
