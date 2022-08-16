import os
import dxpy
import traceback
import subprocess


# Get threads available to this instance
def get_thread_number() -> int:
    threads = os.cpu_count()
    print('Number of threads available: %i' % threads)
    return threads


# This function runs a command on an instance, either with or without calling the docker instance we downloaded
# By default, commands are not run via Docker, but can be changed by setting is_docker = True
def run_cmd(cmd: str, is_docker: bool = False) -> None:

    if is_docker:
        # -v here mounts a local directory on an instance (in this case the home dir) to a directory internal to the
        # Docker instance named /test/. This allows us to run commands on files stored on the AWS instance within Docker
        cmd = "docker run " \
              "-v /home/dnanexus:/test " \
              "-v /home/dnanexus/cadd_files/:/CADD-scripts/data/annotations/ " \
              "-v /home/dnanexus/cadd_precomputed/:/CADD-scripts/data/prescored/GRCh38_v1.6/incl_anno/ " \
              "egardner413/mrcepid-burdentesting " + cmd

    # Standard python calling external commands protocol
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()

    # If the command doesn't work, print the error stream and close the AWS instance out with 'dxpy.AppError'
    if proc.returncode != 0:
        print("The following cmd failed:")
        print(cmd)
        print("STDERROR follows\n")
        traceback.format_exc()
        print(stderr.decode('utf-8'))
        raise dxpy.AppError("Failed to run properly...")


# Utility function to delete files no longer needed from the AWS instance to save space
def purge_file(file: str) -> None:

    cmd = "rm " + file
    run_cmd(cmd)


# This is a helper function to upload a local file and then remove it from the instance.
# This is different than other applets I have written since CADD takes up so much space.
# I don't want to have to use a massive instance costing lots of £s!
def generate_linked_dx_file(file: str) -> dxpy.DXFile:

    linked_file = dxpy.upload_local_file(file)
    purge_file(file)
    return linked_file

