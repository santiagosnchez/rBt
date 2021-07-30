#' read.annot.beast
#'
#' @description
#' This function reads a nexus formated beast file and
#' appends all meta-data and annotations.
#'
#' @details
#' Node ordering is the same as in the \code{read.tree}/\code{read.nexus}
#' functions from ape (cladewise). The metadata is stored as a
#' data.frame named 'metadata'. Only node posterior
#' probabilites are passed as an additional numeric vector
#' named 'posterior'. Note that 'metadata' includes tip
#' annotations as well. The code might be a bit slow for large trees.
#'
#' @param file The path to the MCC file from BEAST
#' @return annotated \code{phylo} object
#' @seealso \code{\link{read.nexus}}
#' @export
#' @examples
#' file <- system.file("data/mcc.tre", package="rBt")
#' tr <- read.annot.beast(file)
#' class(tr)
#' # [1] "phylo"
#' names(tr)
#' #[1] "edge"        "edge.length" "Nnode"       "root.edge"   "tip.label"
#' #[6] "metadata"    "posterior"
#' # for add pp to nodes you could try:
#' plot(tr)
#' nodelabels(text=round(tr$posterior,2), cex=0.5)
#' # for edges, try:
#' tr$ordered.edges <- order.edges(tr)
#' plot(tr)
#' edgelabels(edge=tr$ordered.edges, text=round(tr$posterior,2), cex=0.5)

read.annot.beast <- function(file){
	TREE <- scan(file=file, what=character(), sep="\n", quiet=T)
	TREE <- TREE[ grep("^[[:space:]]*tree", TREE) ]
	if (length(TREE) > 1){
		trs <- read.nexus(file)
		for (itr in 1:length(TREE)){
			tr <- TREE[[itr]]
			tr <- sub("tree.* \\(","\\(", tr)
			edges <- list()
			annot <- list()
			for (i in strsplit(tr, ":")[[1]]){
				if (substr(i, nchar(i), nchar(i)) == "]" | substr(i, nchar(i)-1, nchar(i)-1) == "]"){
					brack <- gregexpr(pattern ='\\[|\\]', i)[[1]]
					start <- rev(brack)[2]
					end <- rev(brack)[1]
					annot <- c(annot, list(substr(i, start, end)))
					end <- start-1
					if (substr(i, 1, 1) == "["){
						start <- brack[2]+1
					} else {
						start <- 1
					}
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
				lst <- rBt:::cldws(b, lst)
			}
			tips = vector()
			tip_num = 1
			for (i in 1:length(nodes)){
			    if (nodes[[i]] != "node"){
					    tips = c(tips, nodes[[i]])
							nodes[[i]] = tip_num
							tip_num = tip_num + 1
					}
			}
			nodes[nodesidx] <- lst$closed
			nodes <- sapply(nodes, as.numeric)
			annot <- lapply(annot, rBt:::process_annot)
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
			trs[[itr]]$posterior <- as.numeric(annot_mat$posterior[lst$currnode:length(nodes)])
			cat("tree #",itr,"\r",sep='')
		}
		message(paste("\nRead ",length(trs)," trees"))
		return(trs)
	} else {
		tr <- TREE[[1]]
		tr <- sub("^[[:space:]]*tree.* \\(","\\(", tr)
		edges <- list()
		annot <- list()
		for (i in strsplit(tr, ":")[[1]]){
			if (substr(i, nchar(i), nchar(i)) == "]" | substr(i, nchar(i)-1, nchar(i)-1) == "]"){
				brack <- gregexpr(pattern ='\\[|\\]', i)[[1]]
				start <- rev(brack)[2]
				end <- rev(brack)[1]
				annot <- c(annot, list(substr(i, start, end)))
				end <- start-1
				if (substr(i, 1, 1) == "["){
					start <- brack[2]+1
				} else {
					start <- 1
				}
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
			lst <- rBt:::cldws(b, lst)
		}
		tips = vector()
		tip_num = 1
		for (i in 1:length(nodes)){
		    if (nodes[[i]] != "node"){
				    tips = c(tips, nodes[[i]])
						nodes[[i]] = tip_num
						tip_num = tip_num + 1
				}
		}
		nodes[nodesidx] <- lst$closed
		nodes <- sapply(nodes, as.numeric)
		annot <- lapply(annot, rBt:::process_annot)
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
		tr$posterior <- as.numeric(annot_mat$posterior[lst$currnode:length(nodes)])
		return(tr)
	}
}
