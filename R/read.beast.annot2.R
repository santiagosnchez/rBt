#' read.beast.annot2
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
#' file <- system.file("data/mcc.tre", package="rBt")
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
#' #[19] "length_95%_HPD_UPPER"   "edge.ordered"           "data"

read.annot.beast2 <- function(file){
	TREE <- scan(file=file, what=character(), sep="\n", quiet=T)
	TREE <- TREE[ grep("^tree", TREE) ]
	if (length(TREE) > 1){
		trs <- read.nexus(file)
		for (itr in 1:length(TREE)){
			tr <- TREE[[itr]]
			tr <- sub("tree.* ","", tr)
			edges <- list()
			annot <- list()
			for (i in strsplit(tr, ":")[[1]]){
				if (substr(i, nchar(i), nchar(i)) == "]"){
					start <- gregexpr(pattern ='\\[|\\]', i)[[1]][1]
					end <- gregexpr(pattern ='\\[|\\]', i)[[1]][2]
					annot <- c(annot, list(substr(i, start, end)))
					end <- start-1
					start <- 1
					edges <- c(edges, list(substr(i, start, end)))
				}
			}
			nodes <- lapply(edges, function(x) {
				if (substr(x, nchar(x), nchar(x)) == ")"){
					return("node")
				}
				else {
					start <- gregexpr(pattern ='[\\(|,][[:alnum:]]+', x)[[1]][1]
					return(substr(x, start+1, nchar(x)))
				}
				})
			backbone <- strsplit(tr, "")[[1]][ grep("\\(|\\)",strsplit(tr, "")[[1]]) ]
			tipsidx <- grep("node",nodes, invert=T)
			nodesidx <- grep("node",nodes)
			Ntips <- length(tipsidx)
			Nnodes <- length(nodesidx)
			lst <- list(opened=vector(), closed=vector(), currnode=Ntips)
			for (b in backbone){
				lst <- .cldws(b, lst)
			}
			nodes[nodesidx] <- as.character(lst$closed)
			nodes <- sapply(nodes, as.numeric)
			annot <- lapply(annot, .process_annot)
			annot_names <- sort(unlist(lapply(annot, function(x) names(x) )))
			annot_names <- annot_names[ !duplicated(annot_names) ]
			annot_mat <- data.frame(matrix(NA, ncol=length(annot_names)+1, nrow=length(nodes)))
			colnames(annot_mat) <- c('node', annot_names)
			annot_mat$node <- nodes
			for (i in annot_names){
				tmp <- lapply(annot, function(x, y=i) x[[paste(y)]] )
				tmp <- lapply(tmp, function(x) if (is.null(x)) return(NA) else return(x) )
				annot_mat[,paste(i)] <- unlist(tmp)
			}
			annot_mat <- annot_mat[order(annot_mat[,1]),]
			rownames(annot_mat) <- 1:dim(annot_mat)[1]
			trs[[itr]]$metadata <- annot_mat
		}
		message(paste("Read ",length(trs)," trees"))
		return(trs)
	} else {
		tr <- TREE[[1]]
		tr <- sub("tree.* ","", tr)
		edges <- list()
		annot <- list()
		for (i in strsplit(tr, ":")[[1]]){
			if (substr(i, nchar(i), nchar(i)) == "]"){
				start <- gregexpr(pattern ='\\[|\\]', i)[[1]][1]
				end <- gregexpr(pattern ='\\[|\\]', i)[[1]][2]
				annot <- c(annot, list(substr(i, start, end)))
				end <- start-1
				start <- 1
				edges <- c(edges, list(substr(i, start, end)))
			}
		}
		nodes <- lapply(edges, function(x) {
			if (substr(x, nchar(x), nchar(x)) == ")"){
				return("node")
			}
			else {
				start <- gregexpr(pattern ='[\\(|,][[:alnum:]]+', x)[[1]][1]
				return(substr(x, start+1, nchar(x)))
			}
			})
		backbone <- strsplit(tr, "")[[1]][ grep("\\(|\\)",strsplit(tr, "")[[1]]) ]
		tipsidx <- grep("node",nodes, invert=T)
		nodesidx <- grep("node",nodes)
		Ntips <- length(tipsidx)
		Nnodes <- length(nodesidx)
		lst <- list(opened=vector(), closed=vector(), currnode=Ntips)
		for (b in backbone){
			lst <- .cldws(b, lst)
		}
		nodes[nodesidx] <- as.character(lst$closed)
		nodes <- sapply(nodes, as.numeric)
		annot <- lapply(annot, .process_annot)
		annot_names <- sort(unlist(lapply(annot, function(x) names(x) )))
		annot_names <- annot_names[ !duplicated(annot_names) ]
		annot_mat <- data.frame(matrix(NA, ncol=length(annot_names)+1, nrow=length(nodes)))
		colnames(annot_mat) <- c('node', annot_names)
		annot_mat$node <- nodes
		for (i in annot_names){
			tmp <- lapply(annot, function(x, y=i) x[[paste(y)]] )
			tmp <- lapply(tmp, function(x) if (is.null(x)) return(NA) else return(x) )
			annot_mat[,paste(i)] <- unlist(tmp)
		}
		annot_mat <- annot_mat[order(annot_mat[,1]),]
		rownames(annot_mat) <- 1:dim(annot_mat)[1]
		tr <- read.nexus(file)
		tr$metadata <- annot_mat
		return(tr)
	}
}









