#####################################################################################################################################################
## argument 1: path to parent directory for the project
## eg: Rscript 10X_split.R /gpfs/gibbs/pi/kaminski/public/Backup/Ruben/BaylorFASTQoutput
#####################################################################################################################################################

args <- commandArgs(trailingOnly=TRUE)
# args[1] is input parent directory

sample.parent.dir <- file.path(args[1], "sample_out")
scripts.parent.dir <- file.path(args[1], "scripts")
  split.scripts.dir <- file.path(scripts.parent.dir, "split")
    if(file.exists(split.scripts.dir)==F){
                dir.create(split.scripts.dir)
    }

soup.sample.names <- list.files(file.path(sample.parent.dir))
soup.sample.dir.paths <- file.path(sample.parent.dir, soup.sample.names)
jobsub.filepath <- file.path(split.scripts.dir,"split.jobsub.bat")

py.path <- "/home/rd796/project/SComatic-main/scripts/SplitBam/SplitBamCellTypes.py"
meta.path <- "/home/rd796/project/ageproj/meta_scomatic.tsv"

## Wipe the old jobsub.bat file
cat("",sep="",file=jobsub.filepath, append=FALSE)
for(i in 1:length(soup.sample.names)){
        script.filepath <- file.path(split.scripts.dir, paste(soup.sample.names[i], "_split.sh", sep = ""))
##### print the script ##  -N ",num.reads.per.cell,";
        cmd.out <- NULL
        cmd.out <- paste("#!/bin/bash\n")
        cmd.out <- paste(cmd.out,"#SBATCH --time=8:00:00 -p bigmem --ntasks=1 --cpus-per-task=8",
          " --mem=200G --job-name=", paste0(soup.sample.names[i], "_split"),
          " -o ",script.filepath,".o%J -e ",script.filepath,".e%J\n",sep="")
        cmd.out <- paste(cmd.out,#"module load SAMtools\n",
#"samtools index ", soup.sample.dir.paths[i],"/aligned_fastq/Aligned.sortedByCoord.out.bam\n",
#"module unload SAMtools\n",
#"mkdir ", soup.sample.dir.paths[i], "/split_bam2\n",
"module load miniconda\n",
"conda activate SComatic\n",sep="")
        cmd.out <- paste(cmd.out, "python ", py.path, " ", 
                    " --bam ", soup.sample.dir.paths[i], "/aligned_fastq/Aligned.sortedByCoord.out.bam",
                    " --meta ", meta.path,
                    " --id ", soup.sample.names[i], " --max_nM 5 --max_NH 1",
                    " --outdir ", soup.sample.dir.paths[i], "/split_bam2", sep="") #NOTE THE 2
cat(cmd.out,file=script.filepath,append=F)
cat("sbatch ",script.filepath,"\n",sep="",file=jobsub.filepath, append=TRUE)
}
system(paste("chmod 700 ",jobsub.filepath,sep=""))

