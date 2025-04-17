#' @title update_data_products
#' @description A function to ease the updating and management of NEON data products in a local database
#' @imports neonstore
#' @param products Character. String for class of products to download allowed values: 'all', 'resources', 'inverts','fish', 'temperature'.
#' @param set_env logical. Default is FALSE. If TRUE it sets the R environmental variable `NEONSTORE_HOME`.
#' @param create_dir logical. Default is FALSE. If TRUE the directory "./ignore/database-files" is created if it doesn't exist.
#'
#'

update_data_products = function(products = NULL, set_env = FALSE, create_dir = FALSE,...){
cat("This will take time for large initial downloads and updates.")

# check that R environment variables exist
if(is.null(Sys.getenv("NEONSTORE_HOME"))) stop("The environmental variable NEONSTORE_HOME does not exist. To create it, set `set_env` = TRUE")
if(!dir.exists("./ignore/database-files") & create_dir == FALSE)  stop("The database location './ignore/database-files' does not exist. To create it, set `create_dir` = TRUE")
## Set the local location of NEON download if `set_env` = TRUE
if(set_env)Sys.setenv("NEONSTORE_HOME" = paste0(getwd(),"/ignore/database-files"))
## Set the local location of NEON download if `create_dir` = TRUE
if(create_dir) if(dir.exists("./ignore/database-files")) cat("Database directory ('./ignore/database-files') already exists") else dir.create("./ignore/database-files")

# do the guts of the function
product_vars = c('all', 'inverts', 'fish', 'resources', 'temperature')
if(is.null(products)) stop(cat("`products` variable is NULL (Default). Must take one of values:",product_vars))  
if(!grep(products, product_vars)) stop(cat("`products` variable is",products,". Must take one of values:",product_vars))

dpList = list(temperature = c("DP1.20053.001","DP1.00002.001"), # surface water and air temperature, respectively
               inverts = c("DP1.20120.001"), # stream invert data product
               fish = c("DP1.20107.001"), # fish data product
               resources = c("DP1.20288.001", # water quality
                             "DP4.00130.001", # discharge
                             "DP1.20016.001", # depth/water surface height
                             "DP1.20053.001", # water temperature
                             # "DP1.20190.001", # rearation
                             "DP1.20166.001", # epilithon
                             "DP1.20191.001", # canopy cover
                             "DP1.00004.001", # air pressure
                             "DP1.20042.001") # PAR at surface
)
if(products == 'all'){
  dps = unname(unlist(dpList))
} else{
  dps = unname(unlist(dpList[grep(products, names(dpList))]))
}

neonstore::neon_download(product = dps)

}
