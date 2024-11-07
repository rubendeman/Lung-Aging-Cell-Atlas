#####################################################################################################################################################
## argument 1: path to parent directory for the project
## eg: Rscript 10X_basecount.R /gpfs/gibbs/pi/kaminski/public/Backup/Ruben/BaylorFASTQoutput
#####################################################################################################################################################

args <- commandArgs(trailingOnly=TRUE)
# args[1] is input parent directory

sample.parent.dir <- file.path(args[1], "sample_out")
scripts.parent.dir <- file.path(args[1], "scripts")
  basecount.scripts.dir <- file.path(scripts.parent.dir, "basecount")
    if(file.exists(basecount.scripts.dir)==F){
                dir.create(basecount.scripts.dir)
    }

soup.sample.names <- list.files(file.path(sample.parent.dir))
soup.sample.dir.paths <- file.path(sample.parent.dir, soup.sample.names)
jobsub.filepath <- file.path(basecount.scripts.dir,"basecount.jobsub.bat")

py.path <- "/home/rd796/project/SComatic-main/scripts/BaseCellCounter/BaseCellCounter.py"
ref.path <- "/home/rd796/project/GRCh38.primary_assembly.genome.fa"

## Wipe the old jobsub.bat file
cat("",sep="",file=jobsub.filepath, append=FALSE)
for(i in 1:length(soup.sample.names)){
        script.filepath <- file.path(basecount.scripts.dir, paste(soup.sample.names[i], "basecount.sh", sep = ""))
##### print the script ##  -N ",num.reads.per.cell,";
        cmd.out <- NULL
        cmd.out <- paste("#!/bin/bash\n")
        cmd.out <- paste(cmd.out,"#SBATCH --time=15:00:00 -p bigmem --ntasks=1 --cpus-per-task=8",
          " --mem=150G --job-name=", paste0(soup.sample.names[i], "basecount"),
          " -o ",script.filepath,".o%J -e ",script.filepath,".e%J\n",sep="")
        cmd.out <- paste(cmd.out,
        "mkdir ", soup.sample.dir.paths[i], "/basecount\n",
"module load miniconda\n",
"conda activate SComatic\n",sep="")
        cmd.out <- paste(cmd.out, "for bam in $(ls -d ",soup.sample.dir.paths[i],"/split_bam2/*bam);do\n", #NOTE THE 2
        "cell_type=$(basename $bam | awk -F'.' '{print $(NF-1)}')\n",
        "temp=",soup.sample.dir.paths[i],"/temp_${cell_type}\n",
        "mkdir -p $temp\n",
        "python ", py.path, 
                    " --bam $bam",
                    " --ref ", ref.path,
                    " --chrom all",
                    " --out_folder ", soup.sample.dir.paths[i], "/basecount",
                    " --min_bq 30 --tmp_dir $temp --nprocs 1\n",
                    "rm -rf $temp\n","done",sep="")
cat(cmd.out,file=script.filepath,append=F)
cat("sbatch ",script.filepath,"\n",sep="",file=jobsub.filepath, append=TRUE)
}
system(paste("chmod 700 ",jobsub.filepath,sep=""))
