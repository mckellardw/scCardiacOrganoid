
library("oompaBase")

# Pull pathways out of a list of CellChat objects
getPathways <- function(
    chat.list,
    combo = c("intersection","union")
    # slot="netP" #TODO
){
  path.list <- lapply(
    chat.list,
    FUN = function(CHAT) CHAT@netP$pathways
  )
  
  if(combo[1]=="intersection"){
    return(
      sort(Reduce(intersect, path.list))
    )
  }else if(combo[1]=="union"){
    return(
      sort(unique(unlist(path.list)))
    )
  }
}