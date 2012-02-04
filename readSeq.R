##' <description>
##' Read sequence files from a directory.
##' <details>
##' @title
##' @param file_type The format of sequence file.
##' @param path Path of a directory.
##' @return A data frame of file names and sequences.
##' @author Ziqian Zhou
readSeq <- function(file_type='seq', path = "."){
    pathOld <- getwd()
    setwd(path)
    pat <- paste(".", file_type, "$", sep="")
    file_list <- dir(pattern=pat)
    seqs <- filenames <- vector()
    for(i in file_list){
        filenames <- c(filenames,substr(i, 1, regexpr(pat, i)-1))
        seqi <- scan(i, what="", sep="\n", quiet =TRUE)
        seqs <- c(seqs,seqi[length(seqi)])
    }
    print(paste(c("Minimum sequence length:", "Maximum sequence length:"), range(nchar(seqs))))
    setwd(pathOld)
    seqs <- tolower(seqs)
    return(data.frame(filenames, seqs, stringsAsFactors =FALSE))
}

##' <description>
##' Match the sequence names from .seq file with spreedsheet.
##' <details>
##' @title
##' @param seqDf
##' @param fileDf Isolates should not be factors
##' @return A data frame containing sequence and other informations and so other information of abnormality.
##' @author Ziqian Zhou
matchSeq <- function(seqDf, fileDf){
    matchedDf <- cbind(fileDf, NA)
    names(matchedDf) <- c(names(fileDf), "Seq")
    for(i in 1:dim(fileDf)[1]){
        isoPos <- grep(fileDf$Isolate[i], seqDf$filenames, value=FALSE)
        if(length(isoPos) > 1){
            matchedDf$Seq[i] <- "Multiple matches"
        }
        else{
            if(length(isoPos) == 0){
                matchedDf$Seq[i] <- "No match"
            }
            else matchedDf$Seq[i] <- seqDf$seqs[isoPos]
        }
    }
    mulMatch <- sum(matchedDf$Seq == "Multiple matches")
    noMatch <- sum(matchedDf$Seq == "No match")
    print(paste("The number of entries with multiple matched sequences:", mulMatch))
    print(paste("The number of entries with no matched sequence:", noMatch))
    return(matchedDf)
}

