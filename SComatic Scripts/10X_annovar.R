#####################################################################################################################################################
## argument 1: path to parent directory for the project
## eg: Rscript 10X_annovar.R /gpfs/gibbs/pi/kaminski/public/Backup/Ruben/BaylorFASTQoutput
#####################################################################################################################################################

args <- commandArgs(trailingOnly=TRUE)
# args[1] is input parent directory

sample.parent.dir <- file.path(args[1], "sample_out")
scripts.parent.dir <- file.path(args[1], "scripts")
  callable.scripts.dir <- file.path(scripts.parent.dir, "annovar")
    if(file.exists(callable.scripts.dir)==F){
                dir.create(callable.scripts.dir)
    }

soup.sample.names <- list.files(file.path(sample.parent.dir))
soup.sample.dir.paths <- file.path(sample.parent.dir, soup.sample.names)
jobsub.filepath <- file.path(callable.scripts.dir,"annovar.jobsub.bat")

## Wipe the old jobsub.bat file
cat("",sep="",file=jobsub.filepath, append=FALSE)
for(i in 1:length(soup.sample.names)){
        script.filepath <- file.path(callable.scripts.dir, paste(soup.sample.names[i], "annovar.sh", sep = ""))
##### print the script ##  -N ",num.reads.per.cell,";
        cmd.out <- NULL
        cmd.out <- paste("#!/bin/bash\n")
        cmd.out <- paste(cmd.out,"#SBATCH --time=8:00:00 -p day --ntasks=1 --cpus-per-task=8",
          " --mem=60G --job-name=", paste0(soup.sample.names[i], "annovar"),
          " -o ",script.filepath,".o%J -e ",script.filepath,".e%J\n",sep="")
        cmd.out <- paste(cmd.out,
        "mkdir ", soup.sample.dir.paths[i], "/annovar\n",
"module load miniconda\n",
"conda activate SComatic\n",sep="")
                    cmd.out <- paste(cmd.out, 
                    "ANNOVAR=/home/rd796/project/annovar\n",
"hummandb=/home/rd796/project/humandb\n",
"grep -v '#' ",soup.sample.dir.paths[i],"/basecall2/",soup.sample.names[i],".calling.step2.pass.tsv | tr '\\t' '-' | awk -F'-' -v OFS='\\t' '{print $1,$2,$3,$4,$5,$0}' > ",soup.sample.dir.paths[i],"/annovar/sample.variants.avinput\n",
"perl $ANNOVAR/table_annovar.pl ",soup.sample.dir.paths[i],"/annovar/sample.variants.avinput $hummandb -buildver hg38 -out ",soup.sample.dir.paths[i],"/annovar/sample.variants.annovar -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a,gnomad_genome -operation g,r,f,f,f,f -nastring . -csvout -polish --otherinfo",
sep="")
cat(cmd.out,file=script.filepath,append=F)
cat("sbatch ",script.filepath,"\n",sep="",file=jobsub.filepath, append=TRUE)
}
system(paste("chmod 700 ",jobsub.filepath,sep=""))