from pathlib import Path
from typing import Tuple

import pytest
from general_utilities.job_management.command_executor import DockerMount, CommandExecutor
from pysam import VariantFile

from filterbcf.methods.vcf_filter import VCFFilter

test_data_dir = Path(__file__).parent / 'test_data'

EXPECTED_VCF_VALUES = [{'final': 845, 'missing': 6, 'original': 835,
                        'vcf': test_data_dir / 'test_input1.vcf.gz',
                        'index': test_data_dir / 'test_input1.vcf.gz.tbi',
                        'fail': 1},
                       {'final': 417, 'missing': 32, 'original': 410,
                        'vcf': test_data_dir / 'test_input2.vcf.gz',
                        'index': test_data_dir / 'test_input2.vcf.gz.tbi',
                        'fail': 0},
                       ]

for vcf_data in EXPECTED_VCF_VALUES:
    assert vcf_data['vcf'].exists()
    assert vcf_data['index'].exists()


@pytest.fixture(scope='function')
def tmp_data_dir(tmp_path_factory) -> Path:
    """Is a fixture to create a tmp_dir so tests do not clutter the test_data directory.

    This fixture just symlinks the fasta files in test_data/ to a tmp_dir.

    :param tmp_path_factory: A tmp_path_factory object from pytest.
    :return: A Pathlike to the temporary directory.
    """

    tmp_data = tmp_path_factory.mktemp('data')

    return tmp_data


def make_vcf_link(tmp_dir, vcf, idx) -> Tuple[Path, Path]:
    """Helper function to get the VCF file and index being tested into a tmp directory.

    :param tmp_dir: The tmp_dir provided by tmp_path_factory in pytest
    :param vcf: A Path to the VCF file being tested (should be in the test_data directory)
    :param idx: A Path to the index file for the VCF file being tested (should be in the test_data directory)
    :return: A tuple of the VCF and index files in the tmp directory.
    """

    tmp_vcf = tmp_dir / 'test_input.vcf.gz'
    tmp_idx = tmp_dir / 'test_input.vcf.gz.tbi'
    vcf.link_to(tmp_vcf)
    idx.link_to(tmp_idx)

    return tmp_vcf, tmp_idx


@pytest.mark.parametrize(argnames=['sites_suffix', 'vcf_info'],
                         argvalues=zip(['.sites.tsv', '.sites.tsv'], EXPECTED_VCF_VALUES)
                         )
def test_genotype_filter(tmp_data_dir, sites_suffix, vcf_info):
    """
    Test for the genotype filter of the VCF files

    :param tmp_data_dir: temporary data directory made by the tmp_data_dir() function
    :param vcf_info: vcf_file attributes (file path, expected count etc.)
    :return: output
    """
    test_mount = DockerMount(tmp_data_dir, Path('/test/'))
    cmd_exec = CommandExecutor(docker_image='egardner413/mrcepid-burdentesting', docker_mounts=[test_mount])
    tmp_vcf, tmp_idx = make_vcf_link(tmp_data_dir, vcf_info["vcf"], vcf_info["index"])

    assert tmp_vcf.exists()

    class_loaded = VCFFilter(Path(tmp_vcf.name), cmd_exec, gq=20, wes=True, testing=True)

    outfile = class_loaded._genotype_filter(Path(tmp_vcf.name), gq=20, wes=True)
    outfile_path = tmp_vcf.parent / outfile

    assert outfile_path.exists()

    vcf_in = VariantFile(outfile_path, "r")
    num_variants = sum(1 for _ in vcf_in)
    print(num_variants)
    assert num_variants == vcf_info['original']


@pytest.mark.parametrize(argnames=['sites_suffix', 'vcf_info'],
                         argvalues=zip(['.sites.tsv', '.sites.tsv'], EXPECTED_VCF_VALUES)
                         )
def test_set_missingness_values(tmp_data_dir, sites_suffix, vcf_info):
    """
    Create and set missingness values in the VCF file

    :param tmp_data_dir: temporary data directory made by the tmp_data_dir() function
    :param vcf_info: vcf_file attributes (file path, expected count etc.)
    :return: output
    """
    test_mount = DockerMount(tmp_data_dir, Path('/test/'))
    cmd_exec = CommandExecutor(docker_image='egardner413/mrcepid-burdentesting', docker_mounts=[test_mount])
    tmp_vcf, tmp_idx = make_vcf_link(tmp_data_dir, vcf_info["vcf"], vcf_info["index"])

    assert tmp_vcf.exists()

    class_loaded = VCFFilter(Path(tmp_vcf.name), cmd_exec, gq=20, wes=True, testing=True)

    outfile = class_loaded._set_missingness_values(Path(tmp_vcf.name))

    outfile_path = tmp_vcf.parent / outfile
    assert outfile_path.exists()

    for record in VariantFile(outfile_path, "r"):
        assert 'F_MISSING' in record.info
        assert 'GTM' in record.info
        assert 'GT0' in record.info
        assert 'GT1' in record.info
        assert 'GT2' in record.info


@pytest.mark.parametrize(argnames=['sites_suffix', 'vcf_info'],
                         argvalues=zip(['.sites.tsv', '.sites.tsv'], EXPECTED_VCF_VALUES)
                         )
def test_set_id(tmp_data_dir, sites_suffix, vcf_info):
    """
    Ensure variant IDs are properly formatted

    :param tmp_data_dir: temporary data directory made by the tmp_data_dir() function
    :param vcf_info: vcf_file attributes (file path, expected count etc.)
    :return: output
    """
    test_mount = DockerMount(tmp_data_dir, Path('/test/'))
    cmd_exec = CommandExecutor(docker_image='egardner413/mrcepid-burdentesting', docker_mounts=[test_mount])
    tmp_vcf, tmp_idx = make_vcf_link(tmp_data_dir, vcf_info["vcf"], vcf_info["index"])

    assert tmp_vcf.exists()

    class_loaded = VCFFilter(Path(tmp_vcf.name), cmd_exec, gq=20, wes=True, testing=True)

    outfile = class_loaded._set_id(Path(tmp_vcf.name))

    outfile_path = tmp_vcf.parent / outfile
    assert outfile_path.exists()

    with VariantFile(outfile_path) as vcf_in:
        # Print the variant records
        for record in vcf_in:
            assert 'chr' in record.id


@pytest.mark.parametrize(argnames=['sites_suffix', 'vcf_info'],
                         argvalues=zip(['.sites.tsv', '.sites.tsv'], EXPECTED_VCF_VALUES)
                         )
def test_set_filter_flags(tmp_data_dir, sites_suffix, vcf_info):
    """
    Ensure the filter flags are working correctly

    :param tmp_data_dir: temporary data directory made by the tmp_data_dir() function
    :param vcf_info: vcf_file attributes (file path, expected count etc.)
    :return: output
    """
    test_mount = DockerMount(tmp_data_dir, Path('/test/'))
    cmd_exec = CommandExecutor(docker_image='egardner413/mrcepid-burdentesting', docker_mounts=[test_mount])
    tmp_vcf, tmp_idx = make_vcf_link(tmp_data_dir, vcf_info["vcf"], vcf_info["index"])

    assert tmp_vcf.exists()

    class_loaded = VCFFilter(Path(tmp_vcf.name), cmd_exec, gq=20, wes=True, testing=True)

    outfile = class_loaded._set_filter_flags(Path(tmp_vcf.name))

    outfile_path = tmp_vcf.parent / outfile
    assert outfile_path.exists()

    count_pass = 0

    with VariantFile(outfile_path) as vcf_in:
        # Print the variant records
        # Each variant record can be converted to a string, which resembles the original VCF line.
        for record in vcf_in:
            # If 'PASS' is in the filter set, increment the counter
            if "FAIL" in record.filter.keys():
                count_pass += 1
        print(f"Number of FAIL records: {count_pass}")
        assert count_pass == vcf_info['fail']


@pytest.mark.parametrize(argnames=['sites_suffix', 'vcf_info'],
                         argvalues=zip(['.sites.tsv', '.sites.tsv'], EXPECTED_VCF_VALUES)
                         )
def test_write_index(tmp_data_dir, sites_suffix, vcf_info):
    """
    Ensure the index file gets created

    :param tmp_data_dir: temporary data directory made by the tmp_data_dir() function
    :param vcf_info: vcf_file attributes (file path, expected count etc.)
    :return: output
    """
    test_mount = DockerMount(tmp_data_dir, Path('/test/'))
    cmd_exec = CommandExecutor(docker_image='egardner413/mrcepid-burdentesting', docker_mounts=[test_mount])
    tmp_vcf, tmp_idx = make_vcf_link(tmp_data_dir, vcf_info["vcf"], vcf_info["index"])

    assert tmp_vcf.exists()

    class_loaded = VCFFilter(Path(tmp_vcf.name), cmd_exec, gq=20, wes=True, testing=True)

    outfile = class_loaded._write_index(Path(tmp_vcf.name))

    outfile_path = tmp_vcf.parent / outfile
    assert outfile_path.exists()
    assert outfile_path.suffix == '.csi'
