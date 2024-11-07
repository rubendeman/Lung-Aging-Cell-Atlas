#####################################################################################################################################################
## argument 1: path to parent directory for the project
## eg: Rscript 10X_callable.R /gpfs/gibbs/pi/kaminski/public/Backup/Ruben/BaylorFASTQoutput
#####################################################################################################################################################

args <- commandArgs(trailingOnly=TRUE)
# args[1] is input parent directory

sample.parent.dir <- file.path(args[1], "sample_out")
scripts.parent.dir <- file.path(args[1], "scripts")
  callable.scripts.dir <- file.path(scripts.parent.dir, "callable")
    if(file.exists(callable.scripts.dir)==F){
                dir.create(callable.scripts.dir)
    }

soup.sample.names <- list.files(file.path(sample.parent.dir))
soup.sample.dir.paths <- file.path(sample.parent.dir, soup.sample.names)
jobsub.filepath <- file.path(callable.scripts.dir,"callable.jobsub.bat")

py.path <- "/home/rd796/project/SComatic-main/scripts/GetCallableSites/GetAllCallableSites.py"

## Wipe the old jobsub.bat file
cat("",sep="",file=jobsub.filepath, append=FALSE)
for(i in 1:length(soup.sample.names)){
        script.filepath <- file.path(callable.scripts.dir, paste(soup.sample.names[i], "callable.sh", sep = ""))
##### print the script ##  -N ",num.reads.per.cell,";
        cmd.out <- NULL
        cmd.out <- paste("#!/bin/bash\n")
        cmd.out <- paste(cmd.out,"#SBATCH --time=8:00:00 -p day --ntasks=1 --cpus-per-task=8",
          " --mem=60G --job-name=", paste0(soup.sample.names[i], "callable"),
          " -o ",script.filepath,".o%J -e ",script.filepath,".e%J\n",sep="")
        cmd.out <- paste(cmd.out,
        "mkdir ", soup.sample.dir.paths[i], "/callable2\n",
"module load miniconda\n",
"conda activate SComatic\n",sep="")
                    cmd.out <- paste(cmd.out, "python ", py.path, 
                    " --infile ", soup.sample.dir.paths[i], "/basecall2/",soup.sample.names[i],".calling.step1.tsv",
                    " --outfile ", soup.sample.dir.paths[i], "/callable2/",soup.sample.names[i],
                    " --max_cov 1\n",sep="")
cat(cmd.out,file=script.filepath,append=F)
cat("sbatch ",script.filepath,"\n",sep="",file=jobsub.filepath, append=TRUE)
}
system(paste("chmod 700 ",jobsub.filepath,sep=""))