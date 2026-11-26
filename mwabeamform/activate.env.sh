# code for activating the environment

module use /software/projects/pawsey0272/setonix/2025.08/modules/zen3/gcc/14.2.0
module load singularity/4.1.0-psr

export psrsrc="singularity exec /software/projects/pawsey0272/zwan4817/container/psr-search.sif"
export psraly="singularity exec /software/projects/pawsey0272/zwan4817/container/psr-analysis.sif" 