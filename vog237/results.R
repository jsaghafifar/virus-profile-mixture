###### Processing results from models constructed by MixtureMaker
###### Outputs:
       #### results_${inputFileBase}.csv (logL, numParams, branchLength, AIC)

setwd('~/virus-profile-mixture/vog237')
require(dplyr)
require(stringr)
filebase = 'test_vog' # set filebase for input+output here

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
file = readLines(input_file)
# finding alignments. check here if grep failing
headers = grep('msa\\:$', file)
alignments = gsub('\\.msa\\:', '', file[headers])
# finding models. check here if grep failing
models = grep('LG\\:|I\\:|C10\\:|R\\d\\:|G\\d\\:', file)
profiles = grep('ModelMixtur.+\\d\\!Component', file)
weights = grep('ModelMixtur.+\\d\\!prop', file)

mixture_models = gsub('\\+.4.+| ', '',
                          file[grep('^ MM[2-9]|^ MM1\\d', file)]) |>
                        unique()
# amino acids, grouped by physicochemical properties
aa = c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
       'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')

x_model = 1
x_profile = 1
x_weight = 1
df_stats = tibble()
df_profiles = tibble()
df_weights = tibble()

# grab results for each alignment
for (i in seq_along(headers)) {
   # grab for each model
   stats_block = models[models >= headers[i] & models < headers[i + 1]]
   for (j in seq_along(stats_block)) {
      if (headers[i] + j != models[x_model]) {
        cat('grep failed?')
        break
      }
      stats = str_split_1(file[models[x_model]], ' ')
      model = gsub('\\:', '', stats[2])
      stats = as.numeric(stats[3:5])

      tib = tibble(Alignment = alignments[i],
                   Model = model,
                   LogL = stats[1],
                   NumParams = stats[2],
                   BranchLength = stats[3],
                   AIC = 2 * stats[2] - 2 * stats[1])
      df_stats = bind_rows(df_stats, tib)
      x_model = x_model + 1
   }

   for (j in seq_along(mixture_models)) {
      num_profiles = gsub('MM', '', mixture_models[j]) |> as.numeric()
      for (i2 in 1:num_profiles) {
         line_component = file[profiles[x_profile]]
         component = gsub('^.+Component|\\!ModelProtein.+',
                           '', line_component) |>
                           as.numeric()
         profile = gsub('^.+\\: ', '',
                                 line_component) |>
                           str_split_1(', ') |>
                           as.numeric() |>
                           format(scientific = FALSE)
         weight = gsub('^.+\\: ', '',
                                 file[weights[x_weight]]) |>
                           str_split_i(', ', i2) |>
                           as.numeric() |>
                           format(scientific = FALSE)
         tib = tibble(Alignment = alignments[i],
                      MixtureModel = mixture_models[j],
                      Component = component,
                      Weight = weight)
         # names(profile) = aa
         tib2 = tibble(profile) |> t()
         colnames(tib2) = aa
         tib = cbind(tib, tib2)
         df_profiles = rbind(df_profiles, tib)
         x_profile = x_profile + 1
      }

      x_weight = x_weight + 1
   }
}

df_profiles = df_profiles |>
                        arrange(Component) |>
                        arrange(MixtureModel) |>
                        arrange(Alignment)

# read in CSV files if already processed:
# df_stats = read.csv(paste0('results_', filebase, '.csv'))
# df_profiles = read.csv(paste0('profiles_', filebase, '.csv'))

# output files
write.csv(df_stats, paste0('results_', filebase, '.csv'))
write.csv(df_profiles, paste0('profiles_', filebase, '.csv'))
write.csv(df_weights, paste0('weights_', filebase, '.csv'))

### insights
# count of best model for each alignment, based on lowest AIC
df_stats |>
    group_by(Alignment) |>
    slice_min(AIC) |>
    ungroup() |>
    count(Model)

tapply(df_stats$AIC, df_stats$Model, mean) |> sort() # avg AIC per model
tapply(df_stats$BranchLength, df_stats$Model, mean) |> sort() # avg branch len

# mixture weights plot
mm10 = df_profiles[which(df_profiles$MixtureModel == 'MM10'), c(1, 3, 4)]
mm10 = mm10 |>
         mutate(Family = gsub('^.+\\d_|_\\d+', '', mm10$Alignment)) |>
         arrange(Family)

plot = ggplot(mm10, aes(Alignment, as.factor(Component))) +
   geom_tile(aes(fill = as.numeric(Weight))) +
   scale_y_discrete(limits = rev) +
   labs(y = 'Component', fill = 'Mixture Weight') +
   scale_x_discrete(label = NULL) +
   facet_grid(cols = vars(Family), scales = 'free_x',
              space = 'free_x', switch = 'x') +
   theme(strip.text.x = element_text(angle = 90),
         text = element_text(family = 'mono')) +
   scale_fill_gradient(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
                       limits = c(0, 1), low = '#f3f3f3', high = '#5e2b97')
ggsave('tmp.png', plot, units = 'px', width = 6000, height = 3000)
