setwd('~/virus-profile-mixture/vog237')
require(taxize)
# taxize_options(ncbi_sleep = 0.5)

species = read.csv('vogdb.species.txt', sep = '\t', header = FALSE)
species = species[2:length(species[, 1]), ]
colnames(species) = c('SpeciesName', 'TaxonID', 'Source', 'SourceVersion')

for (i in seq_along(species[, 1])) {
    id = species$TaxonID[i]
    # for (j in 1:4) { # suggestion from claude to store failed requests:
    #     taxon = tryCatch(ncbi_downstream(id = as.character(id),
    #                                    downto = 'family'),
    #                    error = function(e) e)
    #     if (!inherits(taxon, 'error')) break
    #     message(sprintf('ID %s, attempt %d: %s',
    #                     id, i, conditionMessage(res)))
    #     Sys.sleep(2^i)  # wait 2, 4, 8, 16 s
    # }
    taxon = ncbi_downstream(as.character(id), downto = 'family')
    # TODO not working right now, internal server errors
}
