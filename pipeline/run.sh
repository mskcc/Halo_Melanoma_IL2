export PYTHONPATH=~/.venv/lib64/python3.7/site-packages

## activate virtual environment
~/.venv/bin/python3 ~/.venv/bin/activate_this.py 

bsub -J 'melanoma_pipeline' -e 'melanoma_pipeline.log' -o 'melanoma_pipeline.log' \
  snakemake --jobs 999 --cluster '/home/byrne/halo/dev/melanoma/pipeline/bsub.py' \
    --configfile input/config/study_config.yaml \
    --cluster-config /home/byrne/halo/dev/melanoma/pipeline/lsf.yaml \
    --snakefile /home/byrne/halo/dev/melanoma/pipeline/Snakefile \
    -p -r "$@" \
    all

