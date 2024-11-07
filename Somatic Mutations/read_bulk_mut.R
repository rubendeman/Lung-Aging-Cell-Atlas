setwd("/gpfs/gibbs/project/kaminski/rd796/ageproj")

colNames <- c("chr","pos","ref","alt","context","coverage","alt_count","tissue","sample_id","subject_id")
srr_index <- 9
tab <- read.table('13059_2019_1919_MOESM5_ESM.tsv', sep = "\t", stringsAsFactors = F, header = T)
# Create conversion vector fro sraIds
sraIds_all <- read.table('GtexSraRunTable.txt', sep = "\t", stringsAsFactors = F, header =T)[,16]
sraIds_uniq <- unique(sraIds_all)
sraIds <- setNames(add_lead_zeroes(1:length(sraIds_uniq)), sraIds_uniq)
# Create conversion data.frame for subject ids
getxIds_long <- sraToGtex(sraIds_uniq, formatOut = "long")
getxIds <- gsub("(GTEX-.+?)-.+", "\\1", getxIds_long)
uniq_gtexIds <- unique(getxIds)
uniq_gtexIds <- setNames(add_lead_zeroes(1:length(uniq_gtexIds)), uniq_gtexIds)
gtexIds <- data.frame(sraIds = names(getxIds), gtexIds = getxIds, gtexIds_samples = getxIds_long, gtexIds_de = uniq_gtexIds[getxIds], stringsAsFactors = F)
rownames(gtexIds) <- gtexIds$sraIds
# Convert
sraIds_original <- tab[,srr_index]
tab[,srr_index] <- sraIds[sraIds_original]
tab$subject <- gtexIds[sraIds_original,"gtexIds_de"]
if(colNames != "NULL"){colnames(tab) <- colNames}
tab_all <- tab
tab_all <- cbind(tab_all, gtexIds[sraIds_original,])

sraToGtex <- function(x, sraTable = read.table('GtexSraRunTable.txt', sep = "\t", stringsAsFactors = F, header =T), formatOut = "short") { 
    x <- as.character(x)
    if(!formatOut %in% c("short", "long"))
        stop("Format has to be 'short' or 'long'")
    if(formatOut == "short") {
        return(querySraRunTable(x, refCol = 16, queryCol = 32, sraTable = sraTable))
    } else if(formatOut == "long") {
        return(querySraRunTable(x, refCol = 16, queryCol = 18, sraTable = sraTable))
    }
}
    add_lead_zeroes <- function(x) {
    sprintf(paste0("%0", max(nchar(x)), "d"), x)
}
querySraRunTable <- function(x, refCol, queryCol, sraTable = read.table('GtexSraRunTable.txt', sep = "\t", stringsAsFactors = F, header =T), ordered = T){
    x_original <- x
    x <- as.character(unique(x))
    # Getting ids
    if (length(x) < 500) {
        cutFields <- paste0("cut -f ", refCol, ",", queryCol)
        pattern <- paste(x, collapse = "|")
        pattern <- paste0("'", pattern, "'")
        ids <- system(paste0(c("grep", "-E",  pattern, sraTable, "|", cutFields), collapse = " "), intern = T)
    } else {
        ids <- queryFile(x = x, refCol = refCol, queryCol = queryCol, filepath = sraTable)
    }
    ids <- read.table(textConnection(ids), sep = "\t", header = F, stringsAsFactors = F)
    if(refCol > queryCol)
        ids <- ids[,2:1]
    if(ordered) { 
        rownames(ids) <- ids[,1]  
        # Returning 
        result <- rep(NA, length(x_original))
        names(result) <- x_original
        present <-  names(result)[names(result) %in% rownames(ids)]
        result[present] <- ids[present, 2]
    } else {
        result <- ids
    }
    return (result)
}
queryFile <- function(x, refCol, queryCol, filepath, sep = "\t") {
    con <- file('GtexSraRunTable.txt', "r")
    result <- ""
    while(length(oneline <- readLines(con, n = 1, warn = FALSE)) > 0) {
        myvector <- unlist(strsplit(oneline, sep))
        if(myvector[refCol] %in% x) {
            current_cols <- paste(myvector[sort(c(refCol, queryCol))], collapse = sep)
            result <- paste0(c(result, current_cols), collapse = "\n")
            x <- x[x != myvector[refCol]]
        }
        if(length(x) == 0){break}
    }
    close(con)
    result <- sub("\\n","", result)
    return(result)
}