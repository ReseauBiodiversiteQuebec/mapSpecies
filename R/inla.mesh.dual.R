#' @importFrom future.apply future_lapply
#' @importFrom future nbrOfWorkers

inla.mesh.dual <- function(mesh) {
### Function obtained from : http://www.math.ntnu.no/inla/ r-inla.org/tutorials/spde/R/spde-tutorial-functions.R
###
### This function compute a dual mesh diagram. A dual mesh diagram is
### constructed by joining the centroids of the a triangular mesh.
### The volumes of these dual cells define the weights of an integration scheme
### based at the nodes of the primal mesh.
###
  if (mesh$manifold=='R2') {
    ce <- t(sapply(1:nrow(mesh$graph$tv), function(i)
      colMeans(mesh$loc[mesh$graph$tv[i, ], 1:2])))
    cores<-nbrOfWorkers() # get nbr of workers from the chosen plan
    chunks <- split(1:mesh$n, rep(1:cores, each=ceiling(mesh$n/cores))[1:mesh$n]) # split task in nb cores chunks
    pls <- future_lapply(chunks, function(chunksi){
      plsin<-lapply(chunksi,function(i){
        p <- unique(do.call('rbind', lapply(1:3, function(k) {
          j <- which(mesh$graph$tv[, k] == i)
          if (length(j) > 0){
            return(rbind(ce[j, , drop = FALSE],
                         cbind(
                           mesh$loc[mesh$graph$tv[j, k], 1] +
                             mesh$loc[mesh$graph$tv[j, c(2:3, 1)[k]], 1],
                           mesh$loc[mesh$graph$tv[j, k], 2] +
                             mesh$loc[mesh$graph$tv[j, c(2:3, 1)[k]], 2]
                         ) / 2))
          }else{
            return(ce[j, , drop = FALSE])
          }
        })))
        j1 <- which(mesh$segm$bnd$idx[, 1] == i)
        j2 <- which(mesh$segm$bnd$idx[, 2] == i)
        if ((length(j1) > 0) | (length(j2) > 0)) {
          p <- unique(rbind(
            mesh$loc[i, 1:2],
            p,
            mesh$loc[mesh$segm$bnd$idx[j1, 1], 1:2] /
              2 +
              mesh$loc[mesh$segm$bnd$idx[j1, 2], 1:2] /
              2,
            mesh$loc[mesh$segm$bnd$idx[j2, 1], 1:2] /
              2 +
              mesh$loc[mesh$segm$bnd$idx[j2, 2], 1:2] /
              2
          ))
          yy <- p[, 2] - mean(p[, 2]) / 2 - mesh$loc[i, 2] / 2
          xx <- p[, 1] - mean(p[, 1]) / 2 - mesh$loc[i, 1] / 2
        }else{
          yy <- p[, 2] - mesh$loc[i, 2]
          xx <- p[, 1] - mesh$loc[i, 1]
        }
        pol <- cbind(p[order(atan2(yy, xx)), ])#,ID=i)
        st_polygon(list(pol[c(1:nrow(pol), 1), 1:2]))
      })
      st_as_sf(st_as_sfc(plsin),crs=mesh$crs)
    },future.seed=NULL) # disable warnings about rng
    dmesh<-do.call("rbind",pls)
    st_geometry(dmesh)<-"geometry"
    dmesh
  }else{ 
    stop("It only works for R2!")
  }
}
