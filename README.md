# Scripts for palaeoproteomic analysis of dental enamel proteins from ancient Iberian equids

### General
- Scripts are organized in directories named according to their intended use.
- Scripts are numbered in order of recommended usage starting from 00.
- Unnumbered scripts are supporting scripts that are called by other scripts.


### Usage
- Nearly all scripts ending in "_SLURM.sh" are designed to be submitted to the 
[Slurm workload manager](https://slurm.schedmd.com/sbatch.html) 
on a high-performance computational cluster.
- The scripts were run as arrays (several small jobs, usually one for each individual or protein). This helped avoid the use of loops.
- The arrays were started via parameter files that provide all paths and specify the number of arrays. An example template of a parameter file can be found in `diverse/example_parameter_file.sh`.
- Dependencies are usually activated via `module load NAME_OF_SOFTWARE`.
- File paths given in any script are just placeholders. All file paths must be provided by the user.

### Transparency
If code was created with the aid of artificial intelligence, this is stated at the top of the script as "written with support of large language models". 