export PYTHONPATH=~/.venv/lib64/python3.7/site-packages

## activate virtual environment
~/.venv/bin/python3 ~/.venv/bin/activate_this.py 

snakemake --snakefile /home/byrne/halo/dev/melanoma/pipeline/Snakefile -n --forceall --rulegraph all | dot -Tsvg > melanoma_rules.svg
snakemake --snakefile /home/byrne/halo/dev/melanoma/pipeline/Snakefile -n --forceall --dag all | dot -Tsvg > melanoma_dag.svg

