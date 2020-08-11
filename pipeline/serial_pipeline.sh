export PYTHONPATH=~/.venv/lib64/python3.7/site-packages

## activate virtual environment
~/.venv/bin/python3 ~/.venv/bin/activate_this.py 

nohup snakemake --configfile input/config/study_config.yaml \
  --snakefile /home/byrne/halo/dev/melanoma/pipeline/Snakefile \
  -p \
  all > pipeline.log 2>&1 & 

