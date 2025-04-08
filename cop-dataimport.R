library(reticulate)
install_python() 
virtualenv_create(envname = "python_cop_env")
virtualenv_install("python_cop_env", packages = c("copernicusmarine"))
use_virtualenv("python_cop_env", required = TRUE)
py_install("copernicusmarine")
setwd("C:/Users/Theo/Desktop/R/whale")
use_virtualenv("python_cop_env", required = TRUE)
# Import the os module
os <- import("os")
cm <- import("copernicusmarine", convert = TRUE)

credentials_file <- os$path$expanduser("~/.copernicusmarine/.copernicusmarine-credentials")
if (os$path$exists(credentials_file)) {
  os$remove(credentials_file)
}

cm$login("username", "password")

cm <- import("copernicusmarine")


result <- cm$subset(
  dataset_id = "cmems_mod_glo_phy_my_0.083deg_P1D-m",
  start_datetime = "2021-01-01T00:00:00",
  end_datetime = "2021-06-30T00:00:00",
  variables = list("thetao"),
  minimum_longitude = -124.72,
  maximum_longitude = -119.76,
  minimum_latitude = 35.14,
  maximum_latitude = 38.26,
  #force_download = TRUE, - potentially unnecessary
)