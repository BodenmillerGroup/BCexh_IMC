#' Measure the distance to the cluster center for each cluster cell
#' @param input_sce SingleCellExperiment object
#' @param cellID column name for object/cell IDs
#' @param X_coord column name for X coordinate
#' @param Y_coord column name for Y coordinate
#' @param ImageNumber column name for Image IDs
#' @param cluster column name for cell cluster IDs
#' @param output_colname (optional) column name for the output column, defalut = "distance_to_center"
#' @return new colData entries called "distance_to_center". If cell is not part of cluster, distance-to-center will be NA.
#' @export

calculateDistance <- function(input_sce,
                          cellID,
                          X_coord, Y_coord,
                          ImageNumber,
                          cluster,
                          output_colname = "distance_to_center"){
  # start time
  start = Sys.time()

  # check if required colData is present
  .checkInputDistance(input_sce, cellID, X_coord, Y_coord, ImageNumber, cluster, output_colname)

  # create cur_df with
  input_df <- data.frame(cellID = input_sce[[cellID]],
                         X = input_sce[[X_coord]],
                         Y = input_sce[[Y_coord]],
                         ImageNumber = input_sce[[ImageNumber]],
                         cluster = input_sce[[cluster]])

  # create data.frame for distance-to-center measure
  distance_to_center = data.frame(matrix(ncol=ncol(input_df)+1, nrow = 0))
  colnames(distance_to_center) = colnames(input_df)
  colnames(distance_to_center)[ncol(distance_to_center)] <- "distance_to_center"

  for(j in unique(input_df$ImageNumber)){

    # create multipoint object with all cells from one image
    cell_coord = as.matrix(input_df[input_df$ImageNumber == j,2:3])
    rownames(cell_coord) = input_df[input_df$ImageNumber == j,1]
    cells = sf::st_multipoint(cell_coord)
    cells_sfc = sf::st_cast(sf::st_sfc(cells), "POINT")

    # loop through all clusters and calculate distance to cluster center
    for(i in unique(input_df[input_df$ImageNumber == j & input_df$cluster != 0,]$cluster)){

      # remove cluster 0
      cur_df <- input_df[input_df$cluster == i & input_df$cluster != 0,]

      if(nrow(cur_df) == 0) {
        next
      }

      ## case where cluster is a single cell - no hull can be computed
      if (nrow(cur_df) == 1){
        polygon = sf::st_point(as.matrix(cur_df[,2:3]))
      }

      ## case where cluster includes two cells - no hull can be computed
      if (nrow(cur_df) == 2){
        polygon = sf::st_multipoint(as.matrix(cur_df[,2:3]))
      }

      ## more than 2 cells - hull can be computed
      if (nrow(cur_df) > 2){
        ## compute hull of cluster (default is concave)

          # compute concave hull
          coords <- as.matrix(cbind(cur_df$X, cur_df$Y))
          hull = data.table::as.data.table(concaveman::concaveman(coords, concavity = 1))
          colnames(hull) <- c("X", "Y")

          # create polygon object
          polygon = sf::st_polygon(list(as.matrix(hull)))
      }

      # center of polygon
      polygon_center <- st_centroid(polygon)

      # calculate distance
      cur_df$distance_to_center <- sqrt((cur_df[,2] - polygon_center[1])**2 + (cur_df[,3] - polygon_center[2])**2)

      # add to image-wide data.frame
      distance_to_center <- rbind(distance_to_center, cur_df)
    }
  }

  # add distance_to_center to input_df
  input_df <- left_join(input_df, distance_to_center[,c("cellID", "distance_to_center")])

  end = Sys.time()
  print(end-start)

  # check if cellIDs have same sorting in both data sets - if yes, return sce with new colData
  if(all(input_df$cellID == input_sce[[cellID]]) == TRUE){
    #colData(input_sce)$community <- input_df$community # only return the unique community
    if("distance_to_center" %in% colnames(colData(input_sce))){
      colData(input_sce)$distance_to_center <- NULL
    }
    colData(input_sce)$distance_to_center <- input_df$distance_to_center
    # overwrite column name of community ID column
    names(colData(input_sce))[length(names(colData(input_sce)))] <- output_colname
    print("distance-to-center successfully added to sce object")
    return(input_sce)
  }

  # check if cellIDs have same sorting in both data sets - if not, return input_df instead of sce
  if(!(all(IDs$cellID == input_sce[[cellID]]) == TRUE)){
    print("Output SCE does not contain same cells as input SCE. For safety reasons, data.frame instead of SCE is returned")
    return(input_df)
  }
}
