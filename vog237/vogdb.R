###### Gathering and filtering MSAs from VOGDB
###### Outputs:
       #### vog.select_names.csv (filenames and family names of filtered MSAs)
       #### vog.select_algs.tar.gz (archive of all filtered MSAs)
       #### vog.holdout_algs.tar.gz (archive of holdout data, for model testing)
       #### vog.train_algs.tar.gz (archive of model training data)

setwd('~/virus-profile-mixture/vog237')
require(taxonomizr)
require(dplyr)

dir.create('vogs_raw/') # where 49k unpacked MSAs go
dir.create('vogs_select/') # where 1k filtered MSAs go
dir.create('vogs_holdout/') # where 10% of filtered MSAs go
dir.create('vogs_training/') # where rest of filtered MSAs go



###### IMPORT AND FILTER INTO TABLE

# extract all MSAs. WARNING: 49k items, 800 MB
if (length(list.files('vogs_raw')) == 0) {
    untar('vog.raw_algs.tar.gz', exdir = 'vogs_raw')
}

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
# if (length(list.files('vogs_raw')) == 0) { # WARNING: 49k items, 800 MB
#     untar('vog.raw_algs.tar.gz', exdir = 'vogs_raw')
# }

# only keep families with >10 MSAs. WARNING: 1303 items, 6.4 MB
unique_families = unique(vogs$Family)
x = 1
vogs_10select = tibble()
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
         tibble = tibble(FileName = gsub('vogs_select\\/', '', new),
                         Family = vog$Family)
         vogs_10select = bind_rows(vogs_10select, tibble)
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

# store key info: filenames and family names. should be 1303 rows
write.csv(vogs_10select, 'vog.10select_names.csv', row.names = FALSE)
# items within the vogs_raw/ dir can now be deleted, not needed



###### SUB-SAMPLING INTO HOLD-OUT/TRAINING
#### Run the following code if prev. sections were skipped:
# vogs_10select = read.csv('vog.10select_names.csv')
# if (length(list.files('vogs_select')) == 0) { # WARNING: 1303 items, 6.4 MB
#     untar('vog.select_algs.tar.gz', exdir = 'vogs_select')
# }

families = unique(vogs_10select$Family)

if (length(list.files('vogs_holdout/')) != 0) {
    file.remove(list.files('vogs_holdout/', full.names = TRUE))
}

if (length(list.files('vogs_training/')) != 0) {
    file.remove(list.files('vogs_training/', full.names = TRUE))
}

for (i in seq_along(families)) {
   family_msa = which(vogs_10select$Family == families[i])
   holdout_data = sample(family_msa, ceiling(length(family_msa) / 10))
   training_data = family_msa[!family_msa %in% holdout_data]
   for (j in family_msa) {
      if (j %in% holdout_data) {
        if (j %in% training_data) {
            cat('Dev: found sequence in both holdout and training lists.\n')
            break
        }
        vogs_10select[j, 'Holdout'] = TRUE
        vog = vogs_10select[j, ]
        file.copy(paste0('vogs_select/', vog$FileName),
                  paste0('vogs_holdout/', vog$FileName))
      } else {
         if (j %in% training_data) {
            vogs_10select[j, 'Holdout'] = FALSE
            vog = vogs_10select[j, ]
            file.copy(paste0('vogs_select/', vog$FileName),
                  paste0('vogs_training/', vog$FileName))
         } else {
            cat('Dev: found sequence in neither holdout and training lists.\n')
            break
         }
      }
   }
}

# compress holdout and training MSAs
tar('vog.holdout_algs.tar.gz',
    paste0('vogs_holdout/', list.files('vogs_holdout/')),
    'gzip')
tar('vog.training_algs.tar.gz',
    paste0('vogs_training/', list.files('vogs_training/')),
    'gzip')

# rewrite stored key info to keep whether in training or holdout
write.csv(vogs_10select, 'vog.10select_names.csv', row.names = FALSE)
