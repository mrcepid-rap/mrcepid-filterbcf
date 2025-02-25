import shutil
from collections import Counter
from pathlib import Path
from typing import Tuple

import pytest
from general_utilities.job_management.command_executor import DockerMount, CommandExecutor
from pysam import VariantFile

from filterbcf.methods.vcf_filter import VCFFilter

test_data_dir = Path(__file__).parent / 'test_data'

# Set this flag to True if you want to keep (copy) the temporary output files.
KEEP_TEMP = False

EXPECTED_VCF_VALUES = [{'final': 845, 'original': 835,
                        'vcf': test_data_dir / 'test_input1.vcf.gz',
                        'index': test_data_dir / 'test_input1.vcf.gz.tbi',
                        'fail': 1},
                       {'final': 417, 'original': 410,
                        'vcf': test_data_dir / 'test_input2.vcf.gz',
                        'index': test_data_dir / 'test_input2.vcf.gz.tbi',
                        'fail': 0},
                       ]

for vcf_data in EXPECTED_VCF_VALUES:
    assert vcf_data['vcf'].exists()
    assert vcf_data['index'].exists()


@pytest.fixture
def temporary_path(tmp_path, monkeypatch):
    """
    Prepare a temporary working directory that contains a copy of the test_data
    directory, then change the working directory to it.

    If KEEP_TEMP is True, after the test the entire temporary directory will be copied
    to a folder 'temp_test_outputs' in the project root.
    """
    # Determine where the original test_data directory is located.
    # (Assumes it is at <project_root>/test_data)
    project_root = Path(__file__).parent
    test_data_source = project_root / "test_data"

    # Create the destination folder inside the tmp_path.
    destination = tmp_path / "test_data"
    destination.parent.mkdir(parents=True, exist_ok=True)

    # Copy the entire test_data directory into the temporary directory.
    shutil.copytree(test_data_source, destination)

    # Change the current working directory to the temporary directory.
    monkeypatch.chdir(tmp_path)

    # Yield the temporary directory to the test.
    yield tmp_path

    # After the test, if KEEP_TEMP is True, copy the temporary directory to a persistent location.
    if KEEP_TEMP:
        persistent_dir = project_root / "temp_test_outputs" / tmp_path.name
        persistent_dir.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(tmp_path, persistent_dir, dirs_exist_ok=True)
        print(f"Temporary output files have been copied to: {persistent_dir}")


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


@pytest.mark.parametrize(argnames=['vcf_info', 'wes', 'gq', 'gt_none'],
                         argvalues=[
                             ({'final': 845, 'original': 835,
                               'vcf': test_data_dir / 'test_input1.vcf.gz',
                               'index': test_data_dir / 'test_input1.vcf.gz.tbi',
                               'fail': 1}, True, 20, 212132),
                             ({'final': 417, 'original': 410,
                               'vcf': test_data_dir / 'test_input2_wgs.vcf.gz',
                               'index': test_data_dir / 'test_input2_wgs.vcf.gz.tbi',
                               'fail': 0}, False, 20, 213456),
                             ({'final': 845, 'original': 835,
                               'vcf': test_data_dir / 'test_input1_wgs.vcf.gz',
                               'index': test_data_dir / 'test_input1_wgs.vcf.gz.tbi',
                               'fail': 1}, False, 20, 212163),
                             ({'final': 417, 'original': 410,
                               'vcf': test_data_dir / 'test_input2.vcf.gz',
                               'index': test_data_dir / 'test_input2.vcf.gz.tbi',
                               'fail': 0}, True, 90, 762698),
                             ({'final': 845, 'original': 835,
                               'vcf': test_data_dir / 'test_input1.vcf.gz',
                               'index': test_data_dir / 'test_input1.vcf.gz.tbi',
                               'fail': 1}, True, 10, 199157),
                             ({'final': 417, 'original': 410,
                               'vcf': test_data_dir / 'test_input2.vcf.gz',
                               'index': test_data_dir / 'test_input2.vcf.gz.tbi',
                               'fail': 0}, True, 20, 213337),
                         ])
def test_genotype_filter(temporary_path: Path, vcf_info, wes: bool, gq: int, gt_none: int) -> None:
    """
    Test for the genotype filter of the VCF files. We are testing to ensure that the genotype filter is working.

    :param temporary_path: temporary data directory made by the temporary_path() function
    :param vcf_info: vcf_file attributes (file path, expected count etc.)
    :param wes: boolean flag for WES
    :param gq: genotype quality threshold
    :return: output
    """
    test_mount = DockerMount(temporary_path, Path('/test/'))
    cmd_exec = CommandExecutor(docker_image='egardner413/mrcepid-burdentesting', docker_mounts=[test_mount])
    tmp_vcf, tmp_idx = make_vcf_link(temporary_path, vcf_info["vcf"], vcf_info["index"])

    assert tmp_vcf.exists()

    class_loaded = VCFFilter(Path(tmp_vcf.name), cmd_exec, gq, wes, testing=True)

    outfile = class_loaded._genotype_filter(Path(tmp_vcf.name), gq, wes)
    outfile_path = tmp_vcf.parent / outfile

    assert outfile_path.exists()

    # Count the number of variants
    vcf_in = VariantFile(outfile_path, "r")

    num_variants = sum(1 for _ in vcf_in)
    assert num_variants == vcf_info['original']

    # Count the number of missing genotypes
    vcf_in = VariantFile(outfile_path, "r")
    genotype_counts = Counter()

    for site in vcf_in:
        for sample in site.samples:
            gt = site.samples[sample]['GT']

            if gt is None:
                gt_str = "./."
            elif None in gt:
                gt_str = "/".join(map(lambda x: "." if x is None else str(x), gt))
            else:
                gt_str = "/".join(map(str, gt))

            genotype_counts[gt_str] += 1

    # Debug assertion failures
    assert genotype_counts.get('./.', 0) == gt_none, f"Expected 212132, but got {genotype_counts.get('./.', 0)}"


@pytest.mark.parametrize(
    argnames='vcf_info, expected_fmissing, expected_gt0, expected_gt1, expected_gt2, expected_ac, expected_af, expected_an',
    argvalues=[
        (EXPECTED_VCF_VALUES[0],
         2.553091820795089,
         2575246,
         31329,
         52639,
         147273,
         23.41906427865615,
         5330990),
        (EXPECTED_VCF_VALUES[1],
         1.8032479549874552,
         1236021,
         32924,
         16160,
         92933,
         14.743444535706658,
         2614092)
    ])
def test_set_missingness_values(temporary_path, vcf_info, expected_fmissing, expected_gt0, expected_gt1, expected_gt2,
                                expected_ac, expected_af, expected_an):
    """
    Create and set missingness values in the VCF file

    :param tmp_data_dir: temporary data directory made by the tmp_data_dir() function
    :param vcf_info: vcf_file attributes (file path, expected count etc.)
    :return: output
    """
    test_mount = DockerMount(temporary_path, Path('/test/'))
    cmd_exec = CommandExecutor(docker_image='egardner413/mrcepid-burdentesting', docker_mounts=[test_mount])
    tmp_vcf, tmp_idx = make_vcf_link(temporary_path, vcf_info["vcf"], vcf_info["index"])

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

    vcf_in = VariantFile(outfile_path, "r")
    # Initialize counters
    genotype_counts = Counter()
    f_missing = 0
    ac = 0
    af = 0
    an = 0

    # Iterate through VCF records
    for site in vcf_in:
        # Extract values safely, summing them if they're lists/tuples
        f_missing += (
            sum(site.info.get('F_MISSING', (0,)))
            if isinstance(site.info.get('F_MISSING', (0,)), (tuple, list))
            else site.info.get('F_MISSING', 0)
        )

        ac += (
            sum(site.info.get('AC', (0,)))
            if isinstance(site.info.get('AC', (0,)), (tuple, list))
            else site.info.get('AC', 0)
        )

        af += (
            sum(site.info.get('AF', (0.0,)))
            if isinstance(site.info.get('AF', (0.0,)), (tuple, list))
            else site.info.get('AF', 0.0)
        )

        an += (
            sum(site.info.get('AN', (0,)))
            if isinstance(site.info.get('AN', (0,)), (tuple, list))
            else site.info.get('AN', 0)
        )

        # Count genotypes
        for sample in site.samples:
            gt = site.samples[sample]['GT']
            gt_str = (
                "./."
                if gt is None or (None in gt)
                else "/".join(map(str, gt))
            )
            genotype_counts[gt_str] += 1

    # Run assertions (example values, replace with your expected values)
    assert f_missing == expected_fmissing, f"Expected {expected_fmissing}, got {f_missing}"
    assert ac == expected_ac, f"Expected {expected_ac}, got {ac}"
    assert af == expected_af, f"Expected {expected_af}, got {af}"
    assert an == expected_an, f"Expected {expected_an}, got {an}"
    assert genotype_counts.get("0/0", 0) == expected_gt0, \
        f"Expected {expected_gt0}, got {genotype_counts.get('0/0', 0)}"
    assert genotype_counts.get("0/1", 0) == expected_gt1, \
        f"Expected {expected_gt1}, got {genotype_counts.get('0/1', 0)}"
    assert genotype_counts.get("1/1", 0) == expected_gt2, \
        f"Expected {expected_gt2}, got {genotype_counts.get('1/1', 0)}"


@pytest.mark.parametrize(argnames='vcf_info',
                         argvalues=EXPECTED_VCF_VALUES)
def test_set_id(temporary_path, vcf_info):
    """
    Ensure variant IDs are properly formatted

    :param tmp_data_dir: temporary data directory made by the tmp_data_dir() function
    :param vcf_info: vcf_file attributes (file path, expected count etc.)
    :return: output
    """
    test_mount = DockerMount(temporary_path, Path('/test/'))
    cmd_exec = CommandExecutor(docker_image='egardner413/mrcepid-burdentesting', docker_mounts=[test_mount])
    tmp_vcf, tmp_idx = make_vcf_link(temporary_path, vcf_info["vcf"], vcf_info["index"])

    assert tmp_vcf.exists()

    class_loaded = VCFFilter(Path(tmp_vcf.name), cmd_exec, gq=20, wes=True, testing=True)

    outfile = class_loaded._set_id(Path(tmp_vcf.name))

    outfile_path = tmp_vcf.parent / outfile
    assert outfile_path.exists()

    with VariantFile(outfile_path) as vcf_in:
        # Print the variant records
        for record in vcf_in:
            assert 'chr' in record.id
            # ensure the record.id is formatted as chr_pos_ref_alt
            assert '_' in record.id
            assert f"{record.chrom}_{record.pos}_{record.ref}_{record.alts[0]}" in record.id


@pytest.mark.parametrize(argnames='vcf_info',
                         argvalues=EXPECTED_VCF_VALUES)
def test_set_filter_flags(temporary_path, vcf_info):
    """
    Ensure the filter flags are working correctly

    :param tmp_data_dir: temporary data directory made by the tmp_data_dir() function
    :param vcf_info: vcf_file attributes (file path, expected count etc.)
    :return: output
    """
    test_mount = DockerMount(temporary_path, Path('/test/'))
    cmd_exec = CommandExecutor(docker_image='egardner413/mrcepid-burdentesting', docker_mounts=[test_mount])
    tmp_vcf, tmp_idx = make_vcf_link(temporary_path, vcf_info["vcf"], vcf_info["index"])

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


@pytest.mark.parametrize(argnames='vcf_info',
                         argvalues=EXPECTED_VCF_VALUES)
def test_write_index(temporary_path, vcf_info):
    """
    Ensure the index file gets created

    :param tmp_data_dir: temporary data directory made by the tmp_data_dir() function
    :param vcf_info: vcf_file attributes (file path, expected count etc.)
    :return: output
    """
    test_mount = DockerMount(temporary_path, Path('/test/'))
    cmd_exec = CommandExecutor(docker_image='egardner413/mrcepid-burdentesting', docker_mounts=[test_mount])
    tmp_vcf, tmp_idx = make_vcf_link(temporary_path, vcf_info["vcf"], vcf_info["index"])

    assert tmp_vcf.exists()

    class_loaded = VCFFilter(Path(tmp_vcf.name), cmd_exec, gq=20, wes=True, testing=True)

    outfile = class_loaded._write_index(Path(tmp_vcf.name))

    outfile_path = tmp_vcf.parent / outfile
    assert outfile_path.exists()
    assert outfile_path.suffix == '.csi'
    # assert that the index is not empty
    assert outfile_path.stat().st_size > 0
