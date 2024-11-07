#####################################################################################################################################################
## argument 1: parent directory for the relevant raw fastq files
## argument 2: parent directory for samples' respective output folders
## e.g. Rscript 10X_cutadapt.R /gpfs/gibbs/pi/kaminski/public/Backup/Ruben/BaylorFASTQ /gpfs/gibbs/pi/kaminski/public/Backup/Ruben/BaylorFASTQoutput
#####################################################################################################################################################

## Number of cores to multithread
num.cores <- 8
args <- commandArgs(trailingOnly=TRUE)
### Return error and cancel immediately if input directory doesn't exist
input.parent.dir <- args[1]
if(file.exists(input.parent.dir)==F){
    cat("Error: Input parent directory does not exist... Try again.\n")
    exit()
}
## Generate the output directory structure
output.parent.dir <- args[2]
if(file.exists(output.parent.dir)==F){
    cat("Generating top level directory for output...\n")
    dir.create(output.parent.dir)
}
output.dir <- file.path(output.parent.dir, "sample_out")
if(file.exists(output.dir)==F){
    cat("Generating parent directory for sample output...\n")
    dir.create(output.dir)
}
script.dir <- file.path(output.parent.dir, "scripts")
if(file.exists(script.dir)==F){
    cat("Generating script directory...\n")
    dir.create(script.dir)
}
  cut.script.dir <- file.path(script.dir, "cutadapt")
  if(file.exists(cut.script.dir)==F){
      cat("Generating cutadapt script directory...\n")
      dir.create(cut.script.dir)
  }

## wipe and remake the jobsub file
jobsub.filepath<-file.path(cut.script.dir,"cut.jobsub.bat")
if(file.exists(jobsub.filepath)){
    cat("Creating new jobsub.bat...\n")
    file.remove(jobsub.filepath)
}
file.create(jobsub.filepath)
##Extract the paths for the 10x CellSoup data parent directories (batches)
soup.batch.names <- list.files(file.path(input.parent.dir))
soup.batch.dir.paths <- file.path(input.parent.dir, soup.batch.names)
## make blank record for processed samples outside of the loop
completed.sample.names <- vector()
## Loop over batches
for(i in 1:length(soup.batch.dir.paths)){
      sample.names <- list.files(soup.batch.dir.paths[i])
      sample.paths <- file.path(soup.batch.dir.paths[i], sample.names)
## Loop over samples
    for(j in 1:length(sample.paths)){
#### return paths of all fastq files, recursively in case of subdirectories
#        fastq.paths <-list.files(file.path(sample.paths[j],raw.fastq.dir), full.names = TRUE, recursive = TRUE)
       fastq.paths <-list.files(sample.paths[j], recursive=TRUE, full.names=TRUE)
#Return warning and skip sample if there are non-fastq files included in the folder
#### Distinguish Read1 and Read2
        read1.fastq.paths <-fastq.paths[grep("_R1_", fastq.paths)]
        read2.fastq.paths <-fastq.paths[grep("_R2_", fastq.paths)]
################################################# generate zUMI script for each sample
            read1.fastqs <- paste(read1.fastq.paths, sep = "", collapse = " ")
            read2.fastqs <- paste(read2.fastq.paths, sep = "", collapse = " ")
            sample.output.dir <- file.path(output.dir, sample.names[j])
            if(file.exists(sample.output.dir)==F){
                dir.create(sample.output.dir)
                }
            script.filepath <- file.path(cut.script.dir, paste(sample.names[j], "_cutadapt.sh", sep = ""))
################## Handle R1 stuff
#### Make "R1 merged" directory
            merged.R1.dir <- file.path(sample.output.dir, "merged_r1")
              if(file.exists(merged.R1.dir)==F){
                    dir.create(merged.R1.dir)
              }
#### Make the a merged R1 file of R1.fastq.gz files from this sample in this batch
              if(length(read1.fastq.paths) > 1){
                  cat("Hold on... merging Read1 files for ", sample.names[j], "...\n", sep = "")
              } else {
                  cat("Hold on... transferring large Read1 file for ", sample.names[j], "...\n", sep = "")
              }
              merged.R1.file <- file.path(merged.R1.dir, paste(sample.names[j], "_R1_merged.fastq.gz", sep = ""))
######## Check if sample has replicate data already run w/in loop: if TRUE: append the merged fastq; if FALSE: overwrite
              if((sample.names[j]%in%completed.sample.names)==TRUE){
                  cat("Appending the reads from sequencing replicates.\n")
                  system(paste("cat ",read1.fastqs, " >> ", merged.R1.file, sep = ""))
              } else {
                  system(paste("cat ",read1.fastqs, " > ", merged.R1.file, sep = ""))
              }
################# Handle R2 stuff
#### Make "R2 merged" directory
            merged.R2.dir <- file.path(sample.output.dir, "merged_r2")
              if(file.exists(merged.R2.dir)==F){
                    dir.create(merged.R2.dir)
              }
#### Make the a merged R2 file of R2.fastq.gz files from this sample in this batch
              if(length(read2.fastq.paths) > 1){
                  cat("Hold on... merging Read2 files for ", sample.names[j], "...\n", sep = "")
              } else {
                  cat("Hold on... transferring large Read2 file for ", sample.names[j], "...\n", sep = "")
              }
              merged.R2.file <- file.path(merged.R2.dir, paste(sample.names[j], "_R2_merged.fastq.gz", sep = ""))
######## Check if sample has replicate data already run w/in loop: if TRUE: append the merged fastq; if FALSE: overwrite
              if((sample.names[j]%in%completed.sample.names)==TRUE){
                  cat("Appending the reads from sequencing replicates.\n")
                  system(paste("cat ",read2.fastqs, " >> ", merged.R2.file, sep = ""))
              } else {
              system(paste("cat ",read2.fastqs, " > ", merged.R2.file, sep = ""))
              }
################################### Make the output directory for merged and trimmed fastqs
            trimmed.merged.fastq.dir <- file.path(sample.output.dir,"trimmed_merged_fastq")
             if(file.exists(trimmed.merged.fastq.dir)==F){
                dir.create(trimmed.merged.fastq.dir)
                }
################################### Print the script
            cmd.out <- NULL
            cmd.out <- paste("#!/bin/bash\n\n")
            cmd.out <- paste(cmd.out,"#SBATCH --time=8:00:00 -p bigmem --ntasks=1 --cpus-per-task=8",
          " --mem=200G  --job-name=",paste(sample.names[j], "cut", sep = "."),
            " -o ",script.filepath,".o%J -e ",script.filepath,".e%J\n\n",sep="")
## enter environment w/ cutadapt
            cmd.out <- paste0(cmd.out, "module load miniconda\n", "conda init bash\n","source activate cutadapt_env\n\n")            
## Trim TSO reverse complement (SMART II A oligo); interleave the output and pipe it to another cutadapt run
            cmd.out<- paste0(cmd.out, "cutadapt -G AAGCAGTGGTATCAACGCAGAGTACATGGG --cores=",num.cores," --interleaved ",
## Input
            merged.R1.file," ",
            merged.R2.file, " | \n\n",
####### Run cutadapt again; trim poly A (up to 2 errors); discard reads with >20 N; discard reads w/ R2 length < 30bp
            "cutadapt -A \"A{20}\" --max-n 20 --pair-filter=both --minimum-length=30 --cores=",num.cores," --interleaved",
## Output
            " -o ",file.path(trimmed.merged.fastq.dir, "trimmed.R1.fastq.gz"),
            " -p ",file.path(trimmed.merged.fastq.dir, "trimmed.R2.fastq.gz -"))
        cat(cmd.out,file=script.filepath,append=F)
#### Make sure the sample wasn't already processed as a replicate before appending sbatch to jobsub.bat
        if((sample.names[j]%in%completed.sample.names)==FALSE){
            cat("sbatch ",script.filepath,"\n",sep="",file=jobsub.filepath, append=TRUE)
        }
## record what's so we can check for replicate sequencing data runs and append their reads together
   completed.sample.names <- append(completed.sample.names, sample.names[j])
    }
}
system(paste("chmod 700 ",file.path(cut.script.dir,"cut.jobsub.bat"),sep=""))
#################
#completed.sample.names
cat("Completed generating .sh files for ", length(unique(completed.sample.names)), " samples.\n", sep="")

