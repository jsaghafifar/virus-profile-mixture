###### Gathering and filtering MSAs from VOGDB
###### Outputs:
       #### vog.select_names.csv (filenames and family names of filtered MSAs)
       #### vog.select_algs.tar.gz (archive of all filtered MSAs)
###### see sample.R for how MSAs were sub-sampled for training/hold-out

setwd('~/virus-profile-mixture/vog237')
require(taxonomizr)
require(dplyr)

dir.create('vogs_raw/') # where 49k unpacked MSAs go
dir.create('vogs_select/') # where 1k filtered MSAs go

# extract all MSAs. WARNING: 49k items, 800 MB
if (length(list.files('vogs_raw')) == 0) {
    untar('vog.raw_algs.tar.gz', exdir = 'vogs_raw')
}

###### IMPORT AND FILTER INTO TABLE

# create database of NCBI taxonomy. WARNING: 700 MB
if (!file.exists('taxonomy.sqlite')) { # one-time setup
   prepareDatabase('taxonomy.sqlite', getAccessions = FALSE)
}

# import VOGDB viral species names
species = read.csv('vogdb.species.txt', sep = '\t', header = FALSE)
species = species[2:length(species[, 1]), ]
colnames(species) = c('SpeciesName', 'TaxonID', 'Source', 'SourceVersion')

# match VOGDB entries to NCBI taxonomy, get families
taxa = getTaxonomy(species$TaxonID, 'taxonomy.sqlite', 'family')
species$Family = taxa[, 'family']

vogs = tibble()
x = 1

# filter MSAs in for-loop
for (i in seq_along(list.files('vogs_raw/msa'))) {
    filebase = list.files('vogs_raw/msa/')[i]
    msa = readLines(paste0('vogs_raw/msa/', filebase))
    seqs = grep('>', msa)

    if (length(seqs) >= 10) { # ≥10 taxa in each MSA
        id = sub('^>([0-9]+)\\..*', '\\1', msa[seqs[1]])

        # get family name for MSA
        family = species$Family[which(species$TaxonID == id)]
        family = family[1] # avoid error when >1 instance of ID matches
        if (is.na(family)) { # skip when family = 'NA'
            next
        }
        single_family = TRUE

        # check a single family is shared among all taxa in each MSA
        for (seq in seq_along(seqs)) {
            id = gsub('^>([0-9]+)\\..*', '\\1', msa[seqs[seq]])
            if (!(family %in%
                  species$Family[which(species$TaxonID == id)])) {
                single_family = FALSE
                next
            }
        }

        # store selected MSA in tibble and file
        if (single_family) {
            tibble = tibble(FileName = filebase,
                            Family = family,
                            NumberTaxa = length(seqs))
            vogs = bind_rows(vogs, tibble)
            # writeLines(msa, out)
            cat('Stored MSA', x, 'with', family,
                'family and', length(seqs), 'taxa.\n')
            x = x + 1
        }
    }
}

# store key info: filenames and family names. should be 1564 rows
write.csv(vogs, 'vog.select_names.csv', row.names = FALSE)

###### SAVE SELECT MSAs
#### Run the following code if prev. section was skipped:
# vogs = read.csv('vog.select_names.csv')

# only keep families with >10 MSAs. WARNING: 1k items, ?? MB
unique_families = unique(vogs$Family)
x = 1
for (i in seq_along(unique_families)) {
   family_msa = which(vogs$Family == unique_families[i])
   if (length(family_msa) > 10) {
      for (j in seq_along(family_msa)) {
         vog = vogs[family_msa[j], ]
         old = paste0('vogs_raw/msa/', vog$FileName)
         new = paste0('vogs_select/',
                      gsub('\\.msa', '', vog$FileName),
                      '_', vog$Family,
                      '_', vog$NumberTaxa, '.msa')
         file.copy(old, new)
      }
      cat(vogs$Family[family_msa[1]],
          'family copied. Family', x, 'with', length(family_msa), 'MSAs.\n')
      x = x + 1
   }
}

# compress selected MSAs
tar('vog.select_algs.tar.gz',
    paste0('vogs_select/', list.files('vogs_select/')),
    'gzip')
# items within the vogs_raw/ dir can now be deleted, not needed
