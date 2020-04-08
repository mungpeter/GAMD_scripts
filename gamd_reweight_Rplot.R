
library(tidyr)
library(ggplot2)
library(reshape2)
library(graphics)

#jet.colors = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
jet.colors = colorRampPalette(c("#00007F", "cyan", "#7FFF7F", "yellow", "red","red","red","red","red", "#7F0000"))


setwd('~/xxx_data/5_klb/5_md/6_fgf_site1/2_md/fgf19-wt')

in1 = read.csv(file = 'pmf-fgf19-wt.run_all.rmsd.gyration_max.2D-CE2.xvg', sep = '\t',
                stringsAsFactors = F, header = T, comment.char = '@'); head(in1); str(in1)

# in2 = dcast(in1, sc_4D.7S ~ RMSD); head(in2)

p = ggplot(in1, aes(x = X.RMSD, y = Gyration_max, z = PMF.kcal.mol.)) +
      geom_raster(aes(fill = PMF.kcal.mol.), interpolate = T) +
      scale_fill_gradientn(colours = jet.colors(9)) +
      geom_contour(binwidth = 1, col = 'black', size = .2)
p
