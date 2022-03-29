# Select Province and Field combination
provs    <- c("CK", "JB", "RI", "SS")
#provs    <- c("CK")
fields   <- paste0("F", seq(1,8))
combis   <- sort(apply(expand.grid(provs, fields), 1, paste, collapse="_"))

buffers  <- c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5)
#buffers  <- c(0.5, 1.0)

VIs <- c("NDRE", "NDVI")
