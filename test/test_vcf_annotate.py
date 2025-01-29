# import os
#
# import pytest
#
# from pathlib import Path
# from typing import Tuple
# from pysam import VariantFile
#
# from general_utilities.job_management.command_executor import DockerMount, CommandExecutor
# from filterbcf.methods.vcf_annotate import VCFAnnotate
#
# test_data_dir = Path(__file__).parent / 'test_data'
#
# def test_vcf_annotate():
#
#     test_mount = DockerMount(Path(os.getcwd()), Path('/test/'))
#     cmd_exec = CommandExecutor(docker_image='egardner413/mrcepid-burdentesting', docker_mounts=[test_mount])
#
#     annotate_class = VCFAnnotate(
#         vcf_path=Path("test_data/test_input1.vcf.gz"),
#         filtered_vcf=Path("test_data/test_input1.vcf.filtered.bcf"),
#         additional_annotations=(),
#         cmd_executor=cmd_exec,
#         cadd_executor=cmd_exec
#     )
