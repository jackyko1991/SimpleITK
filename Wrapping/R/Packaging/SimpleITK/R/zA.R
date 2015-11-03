
enumToInteger <- function(name, type) {
    if (is.character(name)) {
        ans <- as.integer(get(paste(".__E__", type, sep = ""), pos='package:SimpleITK')[name])
        if (is.na(ans)) {
            warning("enum not found ", name, " ", type)
        }
        ans
    }
}
enumFromInteger = function(i, type) {
    itemlist <- get(paste(".__E__", type, sep = ""), pos='package:SimpleITK')
    names(itemlist)[match(i, itemlist)]
}


## Set up a map of functions for getting single voxels
createPixLookup <- function( where = topenv(parent.frame()))
{
    m <- get(".__E___itk__simple__PixelIDValueEnum")
    # get rid of unknown type - can't access anyway. There will be errors if
    # it happens.
    # Also need to map Float32 to Float and Float64 to Double
    m <- m[m>0]
    n <- names(m)
    # Turn the names into function names
    ff <- gsub("^sitk(.+)", "Image_GetPixelAs\\1", n)
    ff <- gsub("AsFloat32$", "AsFloat", ff)
    ff <- gsub("AsFloat64$", "AsDouble", ff)

    
    sitkPixelAccessMap <-  mget(ff, envir=where,
                                ifnotfound=rep(NA,length(ff)))
    names(sitkPixelAccessMap) <- n
    assign("sitkPixelAccessMap", sitkPixelAccessMap, envir=where)
}

                                        # experimental bracket operator for images
setMethod('[', "_p_itk__simple__Image",
          function(x,i, j, k, drop=TRUE) {
                                        # check to see whether this is returning a single number or an image
            m <- sys.call()

            imdim <- Image_GetDimension(x)
            if ((length(m)-2) < imdim)
              {
                stop("Image has more dimensions")
              }
            imsize <- rep(1, 5)
            imsize[1:imdim] <- Image_GetSize(x)

            if (missing(i)) {
              i <- 1:imsize[1]
            } else {
              i <- (1:imsize[1])[i]
            }

            if (missing(j)) {
              j <- 1:imsize[2]
            } else {
              j <- (1:imsize[2])[j]
            }
            if (missing(k)) {
              k <- 1:imsize[3]
            } else {
              k <- (1:imsize[3])[k]
            }


            if (any(is.na(c(i,j,k)))) {
              stop("Indexes out of range")
            }
            i <- i - 1
            j <- j - 1
            k <- k - 1
            if ((length(i) == 1) & (length(j) == 1) & (length(k) == 1) ) {
              ## return a single point
              pixtype <- x$GetPixelID()
              aF <- sitkPixelAccessMap[[pixtype]]
              if (inherits(aF, "SWIGFunction")) {
                ## need to check whether we are using R or C indexing.
                return(aF(x, c(i, j,k)))
              } else {
                  stop("Don't know how to return this pixel type\n")
              }
            } else {
              ## construct and return an image
              pixtype <- x$GetPixelIDValue()
              resIm <- SingleBracketOperator(i,j,k,x)
              return(resIm);
            }

          }

          )

setMethod('as.array', "_p_itk__simple__Image",
          function(x, drop=TRUE) {
            sz <- x$GetSize()
            if (.hasSlot(x, "ref")) x = slot(x,"ref")
            ans = .Call("R_swig_ImAsArray", x, FALSE, PACKAGE = "SimpleITK")
            dim(ans) <- sz
            if (drop)
              return(drop(ans))
            return(ans)

            }
          )

as.image <- function(arr, spacing=rep(1, length(dim(arr))),
                     origin=rep(0,length(dim(arr))))
  {
    size <- dim(arr)
    return(ArrayAsIm(arr, size, spacing,origin))
  }
