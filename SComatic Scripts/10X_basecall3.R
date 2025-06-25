#####################################################################################################################################################
## argument 1: path to parent directory for the project
## eg: Rscript 10X_basecall3.R /gpfs/gibbs/pi/kaminski/public/Backup/Ruben/BaylorFASTQoutput
#####################################################################################################################################################

args <- commandArgs(trailingOnly=TRUE)
# args[1] is input parent directory

sample.parent.dir <- file.path(args[1], "sample_out")
scripts.parent.dir <- file.path(args[1], "scripts")
  basecall.scripts.dir <- file.path(scripts.parent.dir, "basecall")
    if(file.exists(basecall.scripts.dir)==F){
                dir.create(basecall.scripts.dir)
    }

soup.sample.names <- list.files(file.path(sample.parent.dir))
soup.sample.dir.paths <- file.path(sample.parent.dir, soup.sample.names)
jobsub.filepath <- file.path(basecall.scripts.dir,"basecall3.jobsub.bat")

py.path1 <- "/home/rd796/project/SComatic-main/scripts/MergeCounts/MergeBaseCellCounts.py"
py.path2 <- "/home/rd796/project/SComatic-main/scripts/BaseCellCalling/BaseCellCalling.step1.py"
py.path3 <- "/home/rd796/project/SComatic-main/scripts/BaseCellCalling/BaseCellCalling.step2.py"
ref.path <- "/home/rd796/project/GRCh38.primary_assembly.genome.fa"
editing='/home/rd796/project/SComatic-main/RNAediting/AllEditingSites.hg38.txt'
PON='/home/rd796/project/SComatic-main/PoNs/PoN.scRNAseq.hg38.tsv'

## Wipe the old jobsub.bat file
cat("",sep="",file=jobsub.filepath, append=FALSE)
for(i in 1:length(soup.sample.names)){
        script.filepath <- file.path(basecall.scripts.dir, paste(soup.sample.names[i], "basecall3.sh", sep = ""))
##### print the script ##  -N ",num.reads.per.cell,";
        cmd.out <- NULL
        cmd.out <- paste("#!/bin/bash\n")
        cmd.out <- paste(cmd.out,"#SBATCH --time=8:00:00 -p day --ntasks=1 --cpus-per-task=8",
          " --mem=60G --job-name=", paste0(soup.sample.names[i], "basecall3"),
          " -o ",script.filepath,".o%J -e ",script.filepath,".e%J\n",sep="")
        cmd.out <- paste(cmd.out,
"module load miniconda\n",
"conda activate SComatic\n",sep="")
                    cmd.out <- paste(cmd.out, "python ", py.path3, 
                    " --infile ", soup.sample.dir.paths[i], "/basecall2/",soup.sample.names[i],".calling.step1.tsv",
                    " --outfile ", soup.sample.dir.paths[i], "/basecall2/", soup.sample.names[i],
                    " --editing ", editing,
                    " --pon ", PON, "\n",
                    "bedtools intersect -header -a ",soup.sample.dir.paths[i],"/basecall2/",soup.sample.names[i],".calling.step2.tsv -b /home/rd796/project/SComatic-main/bed_files_of_interest/UCSC.k100_umap.without.repeatmasker.bed | awk '$1 ~ /^#/ || $6 == \"PASS\"' > ",soup.sample.dir.paths[i],"/basecall2/",soup.sample.names[i],".calling.step2.pass.tsv",
                    sep="")
cat(cmd.out,file=script.filepath,append=F)
cat("sbatch ",script.filepath,"\n",sep="",file=jobsub.filepath, append=TRUE)
}
system(paste("chmod 700 ",jobsub.filepath,sep=""))
