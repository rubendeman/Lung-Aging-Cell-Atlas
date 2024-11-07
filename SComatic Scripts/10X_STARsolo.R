#####################################################################################################################################################
## argument 1: path to parent directory for the project
## eg: Rscript 10X_STARsolo.R /gpfs/gibbs/pi/kaminski/public/Backup/Ruben/BaylorFASTQoutput
#####################################################################################################################################################

args <- commandArgs(trailingOnly=TRUE)
# args[1] is input parent directory

sample.parent.dir <- file.path(args[1], "sample_out")
scripts.parent.dir <- file.path(args[1], "scripts")
  STARsolo.scripts.dir <- file.path(scripts.parent.dir, "STARsolo")
    if(file.exists(STARsolo.scripts.dir)==F){
                dir.create(STARsolo.scripts.dir)
    }

soup.sample.names <- list.files(file.path(sample.parent.dir))
soup.sample.dir.paths <- file.path(sample.parent.dir, soup.sample.names)
jobsub.filepath <- file.path(STARsolo.scripts.dir,"STARsolo.jobsub.bat")

STAR.path <- "/gpfs/gibbs/pi/kaminski/public/softwares/STAR-2.7.6a/bin/Linux_x86_64_static/STAR"
genome.dir <- "/gpfs/gibbs/pi/kaminski/public/Backup/Jonas/genome_index_human/GENCODE_release37_GRCh38.p13/STARindex_noGTF/"
GTF.path <-  "/gpfs/gibbs/pi/kaminski/public/Backup/Jonas/genome_index_human/GENCODE_release37_GRCh38.p13/GTF/gencode.v37.annotation.gtf"

####
CBwhitelist.path <- "/gpfs/gibbs/pi/kaminski/public/Backup/Jonas/softwares/10x_Whitelists/3M-february-2018.txt"
# CBwhitelist.path <- "/gpfs/gibbs/pi/kaminski/public/softwares/cellranger-7.1.0/lib/python/cellranger/barcodes/737K-arc-v1.txt"

STAR.options <- paste("--soloType Droplet --readFilesCommand zcat --soloFeatures Gene GeneFull SJ Velocyto --runThreadN 8 ", 
    "--soloBarcodeReadLength 151 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --twopassMode Basic --soloStrand Forward", sep="") # --soloUMIlen 10 for V2, --soloUMIlen 12 for V3, --soloBarcodeReadLength 151
    
BAM.options <-paste(" --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM",sep="")

## Wipe the old jobsub.bat file
cat("",sep="",file=jobsub.filepath, append=FALSE)
## Loop over samples to extract trimmed R1 and R2 1:length(soup.sample.dir.paths)
for(i in 1:length(soup.sample.names)){
      trimmed.folder.paths <- file.path(soup.sample.dir.paths[i], "trimmed_merged_fastq")
#### Distinguish Read1 and Read2
        read1.fastq.path <- file.path(trimmed.folder.paths , "trimmed.R1.fastq.gz")
        read2.fastq.path <- file.path(trimmed.folder.paths , "trimmed.R2.fastq.gz")
        script.filepath <- file.path(STARsolo.scripts.dir, paste(soup.sample.names[i], "_STARsolo.sh", sep = ""))
##### print the script ##  -N ",num.reads.per.cell,";
        cmd.out <- NULL
        cmd.out <- paste("#!/bin/bash\n")
        cmd.out <- paste(cmd.out,"#SBATCH --time=8:00:00 -p bigmem --ntasks=1 --cpus-per-task=8",
          " --mem=200G --job-name=", paste0(soup.sample.names[i], "_STARsolo"),
          " -o ",script.filepath,".o%J -e ",script.filepath,".e%J\n",sep="")
        cmd.out <- paste(cmd.out, STAR.path, " ", STAR.options, BAM.options," --soloCBwhitelist ", CBwhitelist.path, 
                    " --genomeDir ", genome.dir, " --sjdbGTFfile ", GTF.path,   
                    " --readFilesIn ", read2.fastq.path, " ", read1.fastq.path, 
                    " --outFileNamePrefix ", soup.sample.dir.paths[i], "/aligned_fastq/", sep="")
cat(cmd.out,file=script.filepath,append=F)
cat("sbatch ",script.filepath,"\n",sep="",file=jobsub.filepath, append=TRUE)
}
system(paste("chmod 700 ",jobsub.filepath,sep=""))

