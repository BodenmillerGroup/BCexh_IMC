#' Define a cluster of cells by cellIDs
#' @param input_sce SingleCellExperiment object
#' @param IDs_of_interst vector containing the cellsID's of the cells of interest
#' @param cellID column name for object/cell IDs
#' @param X_coord column name for X coordinate
#' @param Y_coord column name for Y coordinate
#' @param ImageNumber column name for Image IDs
#' @param distance (numeric) maximal distance for a cell considered to be a neighbour
#' @param min_clust_size (numeric) min. number of cells for a community of cells to be considered as cluster
#' @param output_colname (optional) column name for the community ID column, defalut = "cell_cluster"
#' @return new colData entry called "cell_cluster" with unique ID for each cluster
#' @export

findClusters <- function(input_sce, 
                         IDs_of_interest, 
                         cellID, 
                         X_coord, Y_coord, 
                         ImageNumber, 
                         distance, 
                         min_clust_size, 
                         output_colname = "cell_cluster"){
  # start time
  start = Sys.time()

  # check if required colData is present
  .checkInputCluster(input_sce, IDs_of_interest, cellID, X_coord, Y_coord, ImageNumber, distance, min_clust_size, output_colname)

  # subset sce object
  input_sce_sub <- input_sce[,input_sce$cellID %in% IDs_of_interest]

  # create cur_df with
  input_df <- data.frame(cellID = input_sce_sub[[cellID]],
                       X = input_sce_sub[[X_coord]],
                       Y = input_sce_sub[[Y_coord]],
                       ImageNumber = input_sce_sub[[ImageNumber]])

  # add cluster column
  input_df$cluster = NA

  clusters = 1
  for(i in 1:nrow(input_df)){
    # check if there is another cell in the surrounding
    cur_neighbours = input_df[which(sqrt((input_df[i,]$X - input_df$X)^2 + (input_df[i,]$Y - input_df$Y)^2) <= distance &
                                input_df$ImageNumber == input_df[i,]$ImageNumber),]

    # no other cells surrounding the current cell - assign cluster 0
    if(nrow(cur_neighbours) == 1){
      input_df[i,]$cluster = clusters
      clusters = clusters + 1 
    }

    # multiple neighbours in a new cluster
    if(nrow(cur_neighbours) > 1 & all(is.na(cur_neighbours$cluster))){
      input_df[input_df$cellID %in% cur_neighbours$cellID,]$cluster = clusters
      clusters = clusters + 1
    }

    # multiple neighbours in an existing cluster
    if(nrow(cur_neighbours) > 1 & length(unique(sort(cur_neighbours$cluster))) == 1){
      input_df[input_df$cellID %in% cur_neighbours$cellID,]$cluster = unique(sort(cur_neighbours$cluster))
    }

    # multiple neighbours in different clusters
    if(nrow(cur_neighbours) > 1 & length(unique(sort(cur_neighbours$cluster))) > 1){
      existing_clusters = unique(sort(cur_neighbours$cluster))
      lowest_cluster_id = unique(sort(cur_neighbours$cluster))[1]
      input_df[input_df$cellID %in% cur_neighbours$cellID,]$cluster = lowest_cluster_id

      # overwrite cluster id of other cells which are not in cur_neighbours list but still part of this cluster
      input_df[input_df$cluster %in% existing_clusters,]$cluster = lowest_cluster_id
    }
  }

  # remove clusters with less than "min_clust_size number of cells
  for(i in unique(input_df$cluster)){
    if(nrow(input_df[input_df$cluster == i,]) < min_clust_size){
      input_df[input_df$cluster == i,]$cluster = 0
    }
  }

  # merge with complete data set
  IDs = as.data.frame(input_sce[[cellID]])
  colnames(IDs) = "cellID"
  IDs = left_join(IDs, input_df)
  IDs[is.na(IDs$cluster),]$cluster = 0

  # print computing time
  end = Sys.time()
  print(end-start)

  # check if cellIDs have same sorting in both data sets - if yes, return sce with new colData
  if(all(IDs$cellID == input_sce[[cellID]]) == TRUE){
    if("cell_cluster" %in% colnames(colData(input_sce))){
      colData(input_sce)$cell_cluster <- NULL
    }
    colData(input_sce)$cell_cluster <- IDs$cluster
    # overwrite column name for cluster ID
    names(colData(input_sce))[length(names(colData(input_sce)))] <- output_colname
    print("clusters successfully added to sce object")
    return(input_sce)
  }
  # check if cellIDs have same sorting in both data sets - if not, return input_df instead of sce
  if(!(all(IDs$cellID == input_sce[[cellID]]) == TRUE)){
    print("Output SCE does not contain same cells as input SCE. For safety reasons, data.frame instead of SCE is returned")
    return(input_df)
  }
}
