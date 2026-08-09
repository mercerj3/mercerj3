library(haven)

# Save the filtered analytic sample (mydata_f) as an SPSS (.sav) file
write_sav(mydata_f, "~/Desktop/mydata_f.sav")

# Optional: compress = "zsav" writes a smaller zlib-compressed file,
# but requires SPSS 21+ to open
# write_sav(mydata_f, "~/Desktop/mydata_f.zsav", compress = "zsav")

# Verify the file wrote correctly by reading it back
check <- read_sav("~/Desktop/mydata_f.sav")
dim(check)
