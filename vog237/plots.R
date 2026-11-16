###### Visualising Q-matrices constructed by MixtureMaker

setwd('~/virus-profile-mixture/vog237')
require(stringr)
filebase = 'vog' # set filebase for input+output here

input_file = paste0(filebase, '.mm.nex')

# read in file
file = readLines(input_file)
ind_head = grep('model .+=$', file)
ind_mod = c()
aa = c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
       'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
q_list = list()
x_line = 1
x_q = 1

for (mod in seq_along(ind_head)) {
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
    q_list[[x_q]] = q
    x_q = x_q + 1
}

### TODO: plot each matrix as bubble plots. maybe include pi
