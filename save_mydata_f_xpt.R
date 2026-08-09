library(haven)

# Save the filtered analytic sample (mydata_f) as a SAS transport (.xpt) file
write_xpt(mydata_f, "~/Desktop/mydata_f.xpt")

# Optional: for maximum SAS compatibility, use transport version 5 and set
# the internal SAS dataset name (must be <= 8 characters in v5)
# write_xpt(mydata_f, "~/Desktop/mydata_f.xpt", version = 5, name = "MYDATAF")

# Verify the file wrote correctly by reading it back
check <- read_xpt("~/Desktop/mydata_f.xpt")
dim(check)
