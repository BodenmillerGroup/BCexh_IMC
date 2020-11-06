#' Define a community within/surrounding a cluster of cells
#' @param input_sce SingleCellExperiment object
#' @param cellID column name for object/cell IDs
#' @param X_coord column name for X coordinate
#' @param Y_coord column name for Y coordinate
#' @param ImageNumber column name for Image IDs
#' @param cluster column name for cell cluster IDs
#' @param distance (numeric) distance to enlarge (buffer) the cluster
#' @param convex (optional) logical, should the cluster be enlarged based on its convex or concave hull? default = FALSE (concave)
#' @param plot logical, plot images showing cluster and communtiy borders using ggplot. default = FALSE
#' @param output_colname (optional) column name for the community ID column, defalut = "community"
#' @return new colData entries called "community" and "closest_community". Relevant is basically only the closest_community since an object can be part of multiple communities.
#' @export

findCommunity <- function(input_sce,
                          cellID,
                          X_coord, Y_coord,
                          ImageNumber,
                          cluster,
                          distance,
                          convex = F,
                          output_colname = "community",
                          plot = FALSE){
  # start time
  start = Sys.time()

  # check if required colData is present
  .checkInputCommunity(input_sce, cellID, X_coord, Y_coord, ImageNumber, cluster, distance, output_colname)

  # create cur_df with
  input_df <- data.frame(cellID = input_sce[[cellID]],
                         X = input_sce[[X_coord]],
                         Y = input_sce[[Y_coord]],
                         ImageNumber = input_sce[[ImageNumber]],
                         cluster = input_sce[[cluster]])

  # add community ID - empty vector
  input_df$community <- NA
  input_df$unique_community <- NA

  for(j in unique(input_df$ImageNumber)){

    # create multipoint object with all cells from one image
    cell_coord = as.matrix(input_df[input_df$ImageNumber == j,2:3])
    rownames(cell_coord) = input_df[input_df$ImageNumber == j,1]
    cells = sf::st_multipoint(cell_coord)
    cells_sfc = sf::st_cast(sf::st_sfc(cells), "POINT")

    # create data.frame for all border cells of all clusters (image-wise)
    border_cells_clusters = data.frame(matrix(ncol=ncol(input_df), nrow = 0))
    colnames(border_cells_clusters) = colnames(input_df)

    # create data.frame for all community cells
    community_cells = data.table::copy(border_cells_clusters)

    # create list for plots/polygon objects
    polygon_plots = ggplot()

    # loop through all clusters and identify community of each cluster
    for(i in unique(input_df[input_df$ImageNumber == j & input_df$cluster != 0,]$cluster)){
      # remove cluster 0
      cur_df <- input_df[input_df$cluster == i & input_df$cluster != 0,]
      if(nrow(cur_df) == 0) {
        next
      }

      ## case where cluster is a single cell - no hull can be computed
      if (nrow(cur_df) == 1){
        polygon = sf::st_point(as.matrix(cur_df[,2:3]))
        border_cells_clusters <- rbind(border_cells_clusters, cur_df)
      }

      ## case where cluster includes two cells - no hull can be computed
      if (nrow(cur_df) == 2){
        polygon = sf::st_multipoint(as.matrix(cur_df[,2:3]))
        border_cells_clusters <- rbind(border_cells_clusters, cur_df)
      }

      ## more than 2 cells - hull can be computed
      if (nrow(cur_df) > 2){
        ## compute hull of cluster (default is concave)
        if(convex == TRUE){
          # find convex hull
          hull = grDevices::chull(x = cur_df$X, y = cur_df$Y)

          # cells that build the border of a cluster
          border_cells = cur_df[hull,]

          # add cells to the list for the whole image
          border_cells_clusters = rbind(border_cells_clusters, border_cells)
          coordinates = as.matrix(border_cells[,2:3])

          # close the polygon (first row = last row)
          coordinates = rbind(coordinates, coordinates[1,])

          # create polygon object
          polygon = sf::st_polygon(list(coordinates))
        }

        if(convex == FALSE){
          # compute concave hull
          coords <- as.matrix(cbind(cur_df$X, cur_df$Y))
          hull = data.table::as.data.table(concaveman::concaveman(coords, concavity = 1))
          colnames(hull) <- c("X", "Y")

          # return common rows between input and the hull
          # first, we need to round the digits of cur_df since concaveman does not return full-length coordinates
          # (see github issue https://github.com/joelgombin/concaveman/issues/13)
          # is obsolete once the bug is fixed
          cur_df_2 <- data.table::copy(cur_df)
          cur_df_2$X <- round(cur_df_2$X, digits = 4)
          cur_df_2$Y <- round(cur_df_2$Y, digits = 4)

          # find common rows
          border_cells <- dplyr::inner_join(hull, cur_df_2, by = c("X", "Y"))

          # return original rows from data set with full-length cooridantes
          border_cells <- cur_df[cur_df$cellID %in% border_cells$cellID,]

          # add cells to the list for the whole image
          border_cells_clusters = rbind(border_cells_clusters, border_cells)

          # create polygon object
          polygon = sf::st_polygon(list(as.matrix(hull)))
        }
      }

      # buffer/enlarge polygon by "distance" pixels
      polygon_buff = sf::st_buffer(polygon, distance)
      polygon_buff_sfc = sf::st_sfc(polygon_buff)

      if(plot == TRUE){
      polygon_plots <- polygon_plots +
        geom_sf(data = polygon, fill = NA, size = 4, col = "deepskyblue2") +
        geom_sf(data = polygon_buff_sfc, fill = NA, size=4, col = "blue")
      }

      # return coordinates of the enlarged polygon
      enlarged_coordinates = as.data.frame(polygon_buff[1])

      # define cells within this enlarged polygon
      intersect = sf::st_intersects(polygon_buff_sfc, cells_sfc)
      intersect = rownames(as.data.frame(cells[intersect[[1]],]))
      community_cells = rbind(community_cells, input_df[input_df$cellID %in% intersect,])

      # assing cluster ID to cells - if a community already has been assigned
      if(nrow(input_df[input_df$cellID %in% intersect & !is.na(input_df$community),]) > 0){
        input_df[input_df$cellID %in% intersect & !is.na(input_df$community),]$community =
          paste(input_df[input_df$cellID %in% intersect & !is.na(input_df$community),]$community,i, sep="_")
      }
      # assing cluster ID to cells - if no community has been assinged yet
      if(nrow(input_df[input_df$cellID %in% intersect & is.na(input_df$community),]) > 0){
        input_df[input_df$cellID %in% intersect & is.na(input_df$community),]$community = i
      }
    }

    ## find the closest border_cell and assing same community (for visualization purposes)
    # cells with multiple communities
    multi_community <- input_df[grep("_", input_df$community),]

    # select current image
    multi_community <- multi_community[multi_community$ImageNumber == j, ]

    # only cells which are not part of a cluster - (non-chemokine producers)
    multi_community <- multi_community[multi_community$cluster == 0,]
    #print(paste("Number of cells with multiple communities: ", nrow(multi_community), sep = ""))

    # nearest neighbour search (for k = number of border cells)
    if(nrow(multi_community) > 0){
      nn <- RANN::nn2(data = border_cells_clusters[,2:3], query = multi_community[,2:3], k=nrow(border_cells_clusters))

      # first column of nn.idx is the row-index of the border_cell with the lowest distance
      input_df[input_df$cellID %in% multi_community$cellID,]$unique_community <- border_cells_clusters[nn$nn.idx[,1],]$cluster
    }

    if(plot == TRUE){
      # all cluster cells
      clust_cells <- input_df[input_df$ImageNumber == j & input_df$cluster != 0,]
      p <- polygon_plots +
             geom_sf(data = cells_sfc, color=alpha("black",0.3)) +
             geom_point(data = clust_cells, aes(x=X, y=Y), color="red") +
             ggtitle(paste(ImageNumber, j, sep = ": "))
      plot(p + theme_void() + theme(plot.title = element_text(hjust = 0.5)))
    }
  }

  # assign for each cell a unique community
  input_df[input_df$cluster == 0 & is.na(input_df$community),]$unique_community = 0
  input_df[input_df$cluster != 0,]$unique_community = input_df[input_df$cluster != 0,]$cluster
  input_df[is.na(input_df$unique_community),]$unique_community = input_df[is.na(input_df$unique_community),]$community

  end = Sys.time()
  print(end-start)

  # check if cellIDs have same sorting in both data sets - if yes, return sce with new colData
  if(all(input_df$cellID == input_sce[[cellID]]) == TRUE){
    #colData(input_sce)$community <- input_df$community # only return the unique community
    if("closest_community" %in% colnames(colData(input_sce))){
      colData(input_sce)$closest_community <- NULL
    }
    colData(input_sce)$closest_community <- input_df$unique_community
    # overwrite column name of community ID column
    names(colData(input_sce))[length(names(colData(input_sce)))] <- output_colname
    print("communities successfully added to sce object")
    return(input_sce)
  }

  # check if cellIDs have same sorting in both data sets - if not, return input_df instead of sce
  if(!(all(IDs$cellID == input_sce[[cellID]]) == TRUE)){
    print("Output SCE does not contain same cells as input SCE. For safety reasons, data.frame instead of SCE is returned")
    return(input_df)
  }
}
