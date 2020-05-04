# star-bulk-tpn

assess detection of transposable elements from bulk rnaseq

# getting started

```
snakemake --configfile config.yaml --use-conda --use-singularity \
	  --cores 8 -d ~/scratch/200504 --keep-remote \ 
	  --config billing_project=YOUR_GCP_PROJ \ 
	  -u ~/configs/star-bulk-tpn.json \
	  --cluster "sbatch --export=ALL --partition {cluster.partition} --nodes {cluster.n} --time {cluster.time} --ntasks {cluster.tasks} --mem {cluster.mem} \
	   --cpus-per-task={cluster.cpus}" -kp
```
