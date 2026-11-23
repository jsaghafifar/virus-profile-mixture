###### Visualising substitution models by PCA plots
#### Input: NEXUS file of models
#### Output: PCA plot PNG files
####         (exchangeabilities, frequencies, frequencies w/ profile-mixture)

setwd('~/virus-profile-mixture/vog237')
require(tidyverse)
require(factoextra)
# filebase = 'allmodels' # set filebase for input+output here

input_file = 'allmodels.nex'
file = readLines(input_file) # read in file
headers = grep('model .+=$', file) # which lines precurse qmats

# amino acids, grouped by physicochemical properties
aa = c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
       'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')

q_indices = c() # to store lines with exchangeabilities
f_indices = c() # lines with frequencies
x = 1

for (model in seq_along(headers)) {
    name = str_split_i(file[headers[model]], ' ', 2)
    name = gsub('=', '', name)

    # build q-matrix
    q = c()
    q_indices = c(q_indices, (headers[model] + 1):(headers[model] + 19))
    for (i in 1:19) {
        line = gsub('^ ', '', file[q_indices[x]]) # trim space at start
        rates = str_split_1(line, ' ')
        x_rate = 1
        x = x + 1
        for (j in 1:i) {
            from = aa[i + 1]
            to = aa[j]
            q[paste0(from, '_', to)] = as.numeric(rates[x_rate])
            x_rate = x_rate + 1
        }
    }
    for (i in seq_along(q)) q[i] = q[i] / sum(q)

    # grab frequencies
    f_indices = c(f_indices, q_indices[x - 1] + 1) # freq in the line after q
    line = gsub('^ |;$', '', file[f_indices[model]]) # trim space and semicolon
    f = str_split_1(line, ' ') |> as.numeric()
    for (i in seq_along(f)) f[i] = f[i] / sum(f)

    # initialise tibbles
    if (model == 1) {
        q_df = t(tibble(q))
        colnames(q_df) = names(q)
        rownames(q_df) = name
        f_df = t(tibble(f))
        colnames(f_df) = aa
        rownames(f_df) = name
    } else {
        q_df = rbind(q_df, q)
        rownames(q_df)[model] = name
        f_df = rbind(f_df, f)
        rownames(f_df)[model] = name
    }
}

# model groups for PCA
groups = rownames(q_df)
groups[grepl('^Q\\.', groups)] = 'QMaker'
groups[grepl('DAYHOFF|.+CMUT|BLOSUM62|JTT|WAG|LG|PMB', groups)] = 'Classic'
groups[grepl('^MT', groups)] = 'mitochondrial'
groups[grepl('RTREV|HIV.|FLU|FLAVI|VOG', groups)] = 'Viral'

# square-root as scaling (TODO check if appropriate)
q_df2 = sqrt(q_df)

# PCA for exchangeabilities
q_pca = prcomp(q_df2, scale = FALSE, center = TRUE)
# fviz_eig(q_pca, addlabels = TRUE, barfill = 'steelblue',
#          barcolor = 'steelblue', linecolor = 'firebrick',
#          main = 'Scree Plot: Variance Explained by Component')
q_pca_plot = fviz_pca_biplot(q_pca, geom.var = FALSE, col.var = 'black',
                             col.ind = as.factor(groups),
                             repel = TRUE, legend.title = 'Model group') +
                             theme_minimal() +
    theme(text = element_text(family = 'mono'),
        plot.title = element_text(size = 16, hjust = 0.5)) +
        scale_colour_manual(values = c('#7e57c2', '#009688', '#c77025',
        '#000000')) +
                             labs(title = 'PCA of exchangeabilities')
ggsave('pca_q.png', q_pca_plot, units = 'px', width = 3000, height = 3000)

# grab frequency profiles from profile-mixture models
fmix_indices = grep('^frequency', file)
mix_names = c()
for (i in seq_along(fmix_indices)) {
   mixture = str_split_1(file[fmix_indices[i]], ' ')
   mix_names[i] = mixture[2]
   mixture[23] = gsub(';$', '', mixture[23])
   mixture = as.numeric(mixture[4:23])
   for (j in seq_along(mixture)) {
      mixture[j] = mixture[j] / sum(mixture)
   }
   f_df = rbind(f_df, mixture)
   rownames(f_df)[nrow(f_df)] = mix_names[i]
}

f_df2 = sqrt(f_df)
# colnames(f_df) = aa
# rownames(f_df) = c(rownames(q_data), mix_names)
freq_groups = c(groups,
                rep('viral-mixtures', times = nrow(f_df2) - length(groups)))

f_pca = prcomp(f_df2[seq_along(groups), ], scale = FALSE, center = TRUE)
# fviz_eig(f_pca, addlabels = TRUE, barfill = 'steelblue',
#          barcolor = 'steelblue', linecolor = 'firebrick',
#          main = 'Scree Plot: Variance Explained by Component')
f_pca_plot = fviz_pca_biplot(f_pca, geom.var = FALSE, col.var = 'black',
                             col.ind = as.factor(groups),
                             repel = TRUE, legend.title = 'Model group') +
                theme_minimal() +
    theme(text = element_text(family = 'mono'),
        plot.title = element_text(size = 16, hjust = 0.5)) +
        scale_colour_manual(values = c('#7e57c2', '#009688', '#c77025',
        '#000000')) +
                labs(title = 'PCA of frequencies')
ggsave('pca_f.png', f_pca_plot, units = 'px', width = 3000, height = 3000)

fmix_pca = prcomp(f_df, scale = FALSE, center = TRUE)
# fviz_eig(fmix_pca, addlabels = TRUE, barfill = 'steelblue',
#          barcolor = 'steelblue', linecolor = 'firebrick',
#          main = 'Scree Plot: Variance Explained by Component')
fmix_pca_plot = fviz_pca_biplot(fmix_pca, geom.var = FALSE, col.var = 'black',
                                col.ind = as.factor(freq_groups),
                                repel = TRUE, legend.title = 'Model group') +
                    theme_minimal() +
    theme(text = element_text(family = 'mono'),
        plot.title = element_text(size = 16, hjust = 0.5)) +
        scale_colour_manual(values = c('#7e57c2', '#009688', '#c77025',
        '#000000')) +
                    labs(title = 'PCA of frequencies (incl. profile-mixtures)')
ggsave('pca_fmix.png', fmix_pca_plot, units = 'px', width = 3000, height = 3000)
