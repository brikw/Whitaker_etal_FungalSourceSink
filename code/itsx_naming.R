# write function to find start position of ITS1 from output names in ITSx
its1_start <- function(fname) {
 ITS1_positions <- strsplit(as.character(fname), ": ")[[1]][2]
 ITS1_start <- as.numeric(strsplit(ITS1_positions, "-")[[1]][1])
 print(ITS1_start)
}

# write function to find end position of 5.8S
x5.8s_end <- function(fname) {
 x5.8s_positions <- strsplit(as.character(fname), ": ")[[1]][2]
 x5.8s_END <- as.numeric(strsplit(x5.8s_positions, "-")[[1]][2])
 print(x5.8s_END)
}

# write function to get plain name from ITSx output
itsx_name <- function(fname) {
 name <- strsplit(fname, "\\|")[[1]][1]
 print(name)
}
