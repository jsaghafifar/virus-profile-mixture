###### Processing results from models constructed by MixtureMaker
###### Outputs:
       #### results_${inputFileBase}.csv (logL, numParams, branchLength, AIC)

setwd('~/virus-profile-mixture/vog237')
require(dplyr)
require(stringr)
filebase = 'test2_vog20' # set filebase for input+output here

# extract .gz file
input_file = paste0(filebase, '.model')
if (!file.exists(paste0(filebase, '.model'))) {
    if (!file.exists(filebase, '.model.gz')) {
        cat('Input file not found!\n')
        break
    }
    untar(file.path(filebase, '.model.gz'))
}

# get indices of alignments and models in file
res = readLines(paste0(filebase, '.model'))
ind_aln = grep('msa\\:$', res) # finding alignments. check here if grep failing
aln = gsub('\\.msa\\:', '', res[ind_aln])
ind_mod = grep('C10\\:|R4\\:|G4\\:', res) # models. check here if grep failing
n_mod = length(ind_mod) / length(ind_aln)

x = 1
tibble = tibble()

# grab results for each alignment
for (i in seq_along(ind_aln)) {
   # grab for each model
   for (j in 1:n_mod) {
      if (ind_aln[i] + j != ind_mod[x]) {
        cat('grep failed?')
        break
      }
      mod = str_split_1(res[ind_mod[x]], ' ')
      val = as.numeric(mod[3:5])
      tib = tibble(Alignment = aln[i],
                   Model = gsub('\\:', '', mod[2]),
                   LogL = val[1],
                   NumParams = val[2],
                   BranchLength = val[3],
                   AIC = 2 * val[2] - 2 * val[1])
      tibble = bind_rows(tibble, tib)
      x = x + 1
   }
}

# read in CSV file if already processed:
# tibble = read.csv(paste0('results_', filebase, '.csv'))

# count of best model for each alignment, based on lowest AIC
tibble |>
    group_by(Alignment) |>
    slice_min(AIC) |>
    ungroup() |>
    count(Model)

tapply(tibble$AIC, tibble$Model, mean) |> sort() # avg AIC per model
tapply(tibble$BranchLength, tibble$Model, mean) |> sort() # avg branch length

# output file
write.csv(tibble, paste0('results_', filebase, '.csv'))
