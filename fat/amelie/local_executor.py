
import tempfile
import subprocess


def submit_jobs(conda_path, command_to_run, args,
                do_not_submit=False):
    tempdir = tempfile.mkdtemp(dir='/tmp')
    print(tempdir)
    output_files = [tempfile.mkstemp(dir=tempdir, prefix="out_") for _ in args]

    with tempfile.NamedTemporaryFile(prefix='joblist_', mode='w', dir=tempdir, delete=False) as job_file:
        for arg, output_file in zip(args, output_files):
            job_file.write(conda_path + '/' + command_to_run + ' ' + output_file[1] + ' ' + arg + '  \n')

        job_file.flush()
        job_file_path = job_file.name

    print(job_file_path)

    if do_not_submit:
        return {'output': 'dummy', 'run_directory': tempdir}

    with open(job_file_path) as f:
        for line in f:
            command = line.strip().split()
            print("Executing %s" % str(command))
            completed = subprocess.run(command)
            print(completed)

    result = []

    for output_file in output_files:
        with open(output_file[1]) as file:
            result.append(file.read())

    return {'output': result, 'run_directory': tempdir}
