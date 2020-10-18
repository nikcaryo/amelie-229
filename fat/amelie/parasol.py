import tempfile

import subprocess

import os


def submit_jobs(conda_path, command_to_run, args,
                max_jobs=None, num_cpu=1, num_mem=6,
                do_not_submit=False):
    tempdir = tempfile.mkdtemp(dir='/cluster/tmp')
    current_directory = os.getcwd()
    print(tempdir)
    output_files = [tempfile.mkstemp(dir=tempdir, prefix="out_") for _ in args]
    cpu_mem_str = '{use cpu ' + str(num_cpu) + '} {use ram ' + str(num_mem) + 'g}'

    with tempfile.NamedTemporaryFile(prefix='joblist_', mode='w', dir=tempdir, delete=False) as job_file:
        for arg, output_file in zip(args, output_files):
            job_file.write(conda_path + '/' + command_to_run + ' ' + output_file[1] +
                           ' ' + arg + ' ' + cpu_mem_str + '  \n')

        job_file.flush()
        job_file_path = job_file.name

    print(job_file_path)

    if do_not_submit:
        return {'output': 'dummy', 'run_directory': tempdir}

    if max_jobs is None:
        para_command = "bash --login -ic 'cd {} && para -batch={} make {}'".format(
            current_directory,
            tempdir,
            job_file_path,
        )
    else:
        para_command = "bash --login -ic 'cd {} && para -batch={} make {} maxJob={}'".format(
            current_directory, 
            tempdir,
            job_file_path,
            max_jobs,
        )

    print(para_command)
    try:
        subprocess.run(['ssh', 'hoxa', para_command])
    except:
        print("Caught exception")
        para_command = "bash --login -ic 'cd {} && para -batch={} stop'".format(
            current_directory, 
            tempdir,
        )
        subprocess.run(['ssh', 'hoxa', para_command])
        raise

    result = []

    for output_file in output_files:
        with open(output_file[1]) as file:
            result.append(file.read())

    return {'output': result, 'run_directory': tempdir}
