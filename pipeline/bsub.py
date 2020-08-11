#! /home/byrne/.venv/bin/python3
'''
bsub.py

script source: https://github.com/slowkow/snakefiles/blob/master/bsub.py

This script checks a Snakemake job's properties (threads, resources) and chooses
an appropriate LSF queue that meets the requirements. It also automatically
chooses the queue that is least busy unless you already specified a queue.
Usage
-----
Add 'threads' and 'resources' to your resource-intensive rules:
    rule my_rule:
        input: ...
        output ...
        threads: 4
        resources:
            mem = 8000                    # megabytes
            runtime = 35                  # minutes
            queue = 'my_favorite_queue'   # queue name
Invoke snakemake with the path to bsub.py:
    snakemake --jobs 999 --cluster "path/to/bsub.py -o bsub.stdout"
Consider adding bsub.py to a folder in your $PATH, so you can do:
    snakemake --jobs 999 --cluster "bsub.py -o bsub.stdout"
Note
----
For your cluster at your institution, you'll have to modify this script.
'''

import os
import sys
import json
import argparse

from subprocess import check_output

from snakemake.utils import read_job_properties

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("jobscript")
    args = parser.parse_args()

    job_properties = read_job_properties(args.jobscript)

    # By default, we use 1 thread.
    threads = job_properties.get('threads', 1)

    # Set Juno defaults mem (GB) and runtime (minutes)
    mem = int(job_properties['resources'].get('mem', '2'))
    runtime = int(job_properties['resources'].get('runtime', '360'))
    out = job_properties['cluster'].get('stdout', 'stdout.txt').replace("=", "_")
    err = job_properties['cluster'].get('stderr', 'stdout.txt').replace("=", "_")
    name = job_properties['cluster'].get('name', 'melanoma').replace("=", "_")

    # Let the user specify the queue.
    queue = job_properties['resources'].get('queue', None)

    # Otherwise, set to default: general 
    if not queue:
        queue = 'general'    

    # Submit the job to the queue.
    run_bsub(name, queue, threads, mem, runtime, out, err, args.jobscript)

def run_bsub(name, queue, threads, mem, runtime, stdout, stderr, script):
    cmd = 'bsub -J {nm} -n {t}'.format(t=threads, nm = name)
    if mem:
        cmd += ' -R "rusage[mem={}]"'.format(mem)
    if runtime:
        cmd += ' -W {}'.format(runtime)
    if stdout:
        cmd += ' -o {}'.format(stdout)
    if stderr:
        cmd += ' -e {}'.format(stderr)
    cmd += ' {s}'.format(s=script)
    return os.system(cmd)

def get_queue(threads, mem, runtime):
    # All the ERISone queues.
    queues = ['vshort', 'short', 'medium',
              'normal', 'long', 'vlong',
              'big', 'big-multi']
    # Find valid queues for this job's requirements.
    retval = []
    # Only consider 'vshort' if we specify a nonzero runtime.
    if threads == 1 and mem <= 1000 and 0 < runtime <= 15:
        retval.append('vshort')
    # The other queues are all ok if we leave runtime=0.
    if threads == 1 and mem <= 4000 and runtime <= 60:
        retval.append('short')
    if threads <= 4 and mem <= 8000 and runtime <= 60 * 24:
        retval.append('medium')
    if threads <= 6 and mem <= 8000 and runtime <= 60 * 24 * 3:
        retval.append('normal')
    if threads <= 4 and mem <= 8000 and runtime <= 60 * 24 * 7:
        retval.append('long')
    if threads <= 4 and mem <= 4000 and runtime <= 60 * 24 * 7 * 4:
        retval.append('vlong')
    if threads <= 6 and mem > 8000:
        retval.append('big')
    if 8 <= threads <= 16 and mem > 8000:
        retval.append('big-multi')
    # Make sure we have at least one valid queue.
    if not len(retval):
        return None
    # Get the number of currently running jobs on each queue.
    lines = check_output('bqueues').split(b'\n')[1:-1]
    lines = [line.decode('utf-8').split() for line in lines]
    njobs = {x[0]: int(x[7]) for x in lines}
    # Among valid queues, choose the one with fewest running jobs.
    return min(retval, key=lambda j: njobs[j])

if __name__ == '__main__':
    main()
