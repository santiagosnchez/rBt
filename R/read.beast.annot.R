#' read.beast.annot
#'
#' This function reads a nexus formated MCC beast file and 
#' appends all meta-data and annotations. It uses a Perl 
#' script internally that quickly generates a table with 
#' annotations. Node ordering is the same as in the 
#' read.tree/read.nexus function from ape (cladewise).
#'
#' @param file The path to the MCC file from BEAST
#' @return annotated \code{phylo} object
#' @seealso \code{\link{read.nexus}}
#' @export
#' @examples
#' file <- system.file("data/trees/mcc.tre", package="rBt")
#' tr <- read.beast.annot(file)
#' class(tr)
#' # [1] "phylo"
#' names(tr)
#' # [1] "edge"                   "edge.length"            "Nnode"                 
#' # [4] "root.edge"              "tip.label"              "nodes"                 
#' # [7] "CAheight_mean"          "CAheight_median"        "height"                
#' #[10] "height_median"          "length"                 "length_median"         
#' #[13] "posterior"              "CAheight_95%_HPD_LOWER" "CAheight_95%_HPD_UPPER"
#' #[16] "height_95%_HPD_LOWER"   "height_95%_HPD_UPPER"   "length_95%_HPD_LOWER"  
#' #[19] "length_95%_HPD_UPPER"   "edge.ordered"  

read.beast.annot <- function(file){
    tr <- read.nexus(paste(file))
    perlscript <- system.file("data/Beast2AnnTabl.pl", package="rBt")
    out <- system(paste("perl", perlscript, file, collapse=" "), intern=TRUE)
    out <- strsplit(out, "\t")
    df <- do.call(rbind.data.frame, args=list(out, stringsAsFactors=FALSE))
    colnames(df) <- df[1,]
    rownames(df) <- df[,1]
    df <- as.data.frame(t(df[-1,-1]), stringsAsFactors=FALSE)
    rmrange <- grep("range", colnames(df))
    df <- df[,-rmrange]
    spltind <- grep("HPD",df[dim(df)[1],])
    dfsplt <- data.frame(rep(0,dim(df)[1]))
    for (i in spltind){
        tmp <- strsplit(df[,i], ",")
        for (j in grep("NA",tmp))
            tmp[[j]] <- c(NA,NA)
        tmp <- as.numeric(unlist(tmp))
        tmp <- data.frame(matrix(tmp, ncol=2, byrow=TRUE))
        dfsplt <- cbind(dfsplt,tmp)
    }
    dfsplt <- dfsplt[,-1]
    colnames(dfsplt) <- sort(c(paste(colnames(df)[spltind], "_LOWER", sep=""),paste(colnames(df)[spltind], "_UPPER", sep="")))
    df <- df[,-spltind]
    for (i in 1:dim(df)[2]){
        df[ df[,i] == "NA" ,i] <- NA
        df[ ,i] <- as.numeric(df[,i])
    }
    df <- cbind(df,dfsplt)
    rt <- length(tr$tip.label)+1
    nodes <- (length(tr$tip.label)+1):(length(tr$tip.label)+tr$Nnode)
    tr$nodes <- nodes
    for (i in colnames(df))
        tr[[paste(i)]] <- df[rt:dim(df)[1],paste(i)]
    pp <- prop.part(tr)
    edg <- sapply(pp, function(x,y=tr) min(which.edge(y,x))-1 )
    edg <- edg[-1]
    tr$edge.ordered <- edg
    return(tr)
}


