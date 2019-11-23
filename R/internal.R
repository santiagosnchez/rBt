# hidden functions used by read.annot.beast

cldws <- function(b, lst){
    if (b == "("){
        lst$currnode <- lst$currnode+1
        if (length(grep(paste(lst$currnode), lst$opened)) == 1 | length(grep(paste(lst$currnode), lst$closed)) == 1){
            return(cldws(b, lst))
        }
        else {
            lst$opened <- c(lst$opened, lst$currnode)
            return(lst)
        }
    }
    if (b == ")"){
        if (length(grep(paste(lst$currnode), lst$closed)) == 1){
            lst$currnode <- lst$currnode-1
            return(cldws(b, lst))
        }
        else {
            lst$closed <- c(lst$closed, lst$currnode)
            return(lst)
        }
    }
}

process_annot <- function(an){
    if (is.na(an)){
        return(NA)
    } else {
        res <- list()
        currlybr <- gregexpr(pattern ='\\{|\\}', an[[1]])[[1]]
        if (currlybr[1] == -1){
            sep_annot <- strsplit(an,",")[[1]]
        } else {
            commas <- vector()
            sep_annot <- vector()
            ctmp <- gregexpr(pattern =',', an[[1]])[[1]]
            currlybr <- data.frame(st=currlybr[seq(1,length(currlybr),2)],en=currlybr[seq(2,length(currlybr),2)])
            check_comma <- function(b, cm)
                return(any(apply(currlybr, 1, function(x, y=cm) return(x[1] < y & x[2] > y) )))
            for (cm in ctmp)
                if (!check_comma(currlybr, cm))
                    commas <- c(commas,cm)
            if (length(commas) == 0){
                sep_annot <- an
            }
            else if (length(commas) == 1){
                sep_annot <- c(sep_annot, substr(an, 1, commas-1))
                sep_annot <- c(sep_annot, substr(an, commas+1, nchar(an)))
            } else {
                for (i in 1:length(commas)){
                    if (i == 1){
                        sep_annot <- c(sep_annot, substr(an, 1, commas[i]-1))
                        sep_annot <- c(sep_annot, substr(an, commas[i]+1, commas[i+1]-1))
                    }
                    else if (i == length(commas)){
                        sep_annot <- c(sep_annot, substr(an, commas[i]+1, nchar(an)))
                    }
                    else {
                        sep_annot <- c(sep_annot, substr(an, commas[i]+1, commas[i+1]-1))
                    } 
                }
            }
        }
        sep_annot <- gsub("\\[|\\[\\&|\\]|\\{|\\}","", sep_annot)
        for (i in sep_annot){
            tmp <- strsplit(i, "=")[[1]]
            res[[paste(tmp[1])]] <- tmp[2]
        }
        return(res)
    }
}

# for genealogical sorting index
get_ancestors = function(phy, tips){
    if (is.character(tips) == TRUE)
        tips = sapply(tips, function(x) which(phy$tip.label == x))
    nodes = vector()
    mrca = which.node(phy, tips)
    for (tip in tips){
        x = 0
        e = tip
        while (x != mrca){
            x = edge[ which(phy$edge[,2] == e), 1]
            nodes = append(nodes, x)
            e = x
        }
    }
    nodes = unique(nodes)
    return(nodes)
}
