# write function to get out taxonomic levels information from UNITE headers
get.taxonomy <-function(fname) {
    info_kingdom <- strsplit(basename(fname), ";")[[1]][1] 
    phylum <- strsplit(basename(fname), ";")[[1]][2] 
        phylum2 <- strsplit(phylum, "__")[[1]][2]
    class <- strsplit(basename(fname), ";")[[1]][3]
        class2 <- strsplit(class, "__")[[1]][2]
    order <- strsplit(basename(fname), ";")[[1]][4] 
        order2 <- strsplit(order, "__")[[1]][2]
    family <- strsplit(basename(fname), ";")[[1]][5] 
        family2 <- strsplit(family, "__")[[1]][2]
    genus <- strsplit(basename(fname), ";")[[1]][6] 
        genus2 <- strsplit(genus, "__")[[1]][2]
    species <- strsplit(basename(fname), ";")[[1]][7] 
        species2 <- strsplit(species, "__")[[1]][2]
    combo <- rbind(info_kingdom, phylum2, class2, order2, family2, genus2, 
                   species2)
    colnames(combo) <- info_kingdom
    combo
}

# write function just to get out main name from UNITE sequences
spNameFunc <- function(fname) {
 inner <- strsplit(basename(fname), "\\|")[[1]][1] 
 print(inner)
 }
