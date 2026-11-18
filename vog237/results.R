###### Processing results from models constructed by MixtureMaker
###### Outputs:
       #### results_${inputFileBase}.csv (logL, numParams, branchLength, AIC)

setwd('~/virus-profile-mixture/vog237')
require(dplyr)
require(stringr)
filebase = 'vog' # set filebase for input+output here
file_iqtree = 'test_vog.iqtree'

# extract .gz file
input_file = paste0(filebase, '.model')
if (!file.exists(input_file)) {
    if (!file.exists(paste0(input_file, '.gz'))) {
        cat('Input file not found!\n')
        break
    }
    require(R.utils)
    gunzip(paste0(input_file, '.gz'), remove = FALSE)
}

# get indices of alignments and models in model file
res = readLines(input_file)
ind_aln = grep('msa\\:$', res) # finding alignments. check here if grep failing
aln = gsub('\\.msa\\:', '', res[ind_aln])
ind_mod = grep('C10\\:|R4\\:|G4\\:', res) # models. check here if grep failing
n_mod = length(ind_mod) / length(ind_aln)

# get indices of mixture weight info in iqtree file
res_mix = readLines(file_iqtree)
ind_head = grep('^SUBSTITUTION PROCESS$', res_mix)
ind_mix = c((ind_head + 6):(ind_head + 6 + length(ind_aln)))

x = 1
tibble = tibble()
mix_tibble = tibble()

# grab results for each alignment
for (i in seq_along(ind_aln)) {
   # grab best model
   best_mod = str_split_i(res_mix[ind_mix[i]], '  ', 3)

   # grab for each model
   for (j in 1:n_mod) {
      if (ind_aln[i] + j != ind_mod[x]) {
        cat('grep failed?')
        break
      }
      mod = str_split_1(res[ind_mod[x]], ' ')
      val = as.numeric(mod[3:5])
      mod = gsub('\\:', '', mod[2])
      if (best_mod == mod) {
         if (mod == 'MM10+R4' || mod == 'MM10+G4') {
            mix = str_split_1(res_mix[ind_mix[i]], 'MM10Q\\+FMM10pi')
            mix = gsub('\\d+\\:1:|\\,|\\}.+', '', mix[2:11])
            mix = format(as.numeric(mix), scientific = FALSE)
            tib = tibble(Alignment = aln[i],
                         BestModel = best_mod,
                         MixWeight1 = mix[1], # TODO this looks awful
                         MixWeight2 = mix[2], # TODO grab freqs for profiles too
                         MixWeight3 = mix[3],
                         MixWeight4 = mix[4],
                         MixWeight5 = mix[5],
                         MixWeight6 = mix[6],
                         MixWeight7 = mix[7],
                         MixWeight8 = mix[8],
                         MixWeight9 = mix[9],
                         MixWeight10 = mix[10])
            mix_tibble = bind_rows(mix_tibble, tib)
         }
      }

      tib = tibble(Alignment = aln[i],
                   Model = mod,
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

# output files
write.csv(tibble, paste0('results_', filebase, '.csv'))
write.csv(mix_tibble, paste0('mixweights_', filebase, '.csv'))
