library(FCBF)
library(CEMiTool)
data(expr0)
data("sample_annot")
target <- as.factor(sample_annot$Class)
source("./R/fcoex.R")

fc <- new_fcoex(expr0, target)

fc <- discretize(fc)

fc <- find_cbf_modules(fc,FCBF_threshold = 0.15)

fc_man <- fc
fc <- mod_colors(fc)
fc@mod_colors
