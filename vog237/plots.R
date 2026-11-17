###### Visualising Q-matrices constructed by MixtureMaker

setwd('~/virus-profile-mixture/vog237')
require(stringr)
require(ggplot2)
require(tidyverse)
filebase = 'vog' # set filebase for input+output here

input_file = paste0(filebase, '.mm.nex')
file = readLines(input_file) # read in file
ind_head = grep('model .+=$', file) # which lines precurse qmats

# amino acids grouped by physicochemical properties
aa = c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
       'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')

ind_mod = c() # to store which lines from file
x_line = 1
q_list = list() # to store qmats
x_q = 1

for (mod in seq_along(ind_head)) {
    # build matrix
    q = matrix(NA, nrow = 20, ncol = 20, byrow = TRUE, dimnames = list(aa, aa))
    ind_mod = c(ind_mod, (ind_head[mod] + 1):(ind_head[mod] + 19))
    for (i in 2:20) {
        rates = str_split_1(file[ind_mod[x_line]], ' ')
        x_rate = 1
        x_line = x_line + 1
        for (j in 1:(i - 1)) {
            q[i, j] = rates[x_rate + 1]
            x_rate = x_rate + 1
        }
    }

    # convert matrix into tibble
    data = as_tibble(q, .name_repair = ~paste0('V', aa)) |>
        mutate(row = aa) |>
        pivot_longer(cols = starts_with('V'),
                     names_pattern = 'V(.+)', names_to = 'col')
    data$row = factor(data$row, levels = unique(data$row))
    data$col = factor(data$col, levels = unique(data$col))

    # plot
    plot = ggplot(data = data, aes(x = col, y = row,
                                            size = as.numeric(value))) +
        geom_point(shape = 21, colour = 'black', fill = 'white', stroke = 1) +
        scale_y_discrete(limits = rev) +
        labs(x = NULL, y = NULL) +
        theme(legend.position = 'none',
            panel.background = element_blank(),
            panel.grid = element_blank(),
            axis.ticks = element_blank())

    # save plot images
    name = str_split_i(file[ind_head[mod]], ' ', 2)
    pngtitle = paste0(filebase, '_', name, '.png')
    ggsave(pngtitle, plot, units = 'px', width = 1440, height = 1440)

    q_list[[x_q]] = q
    x_q = x_q + 1
}
