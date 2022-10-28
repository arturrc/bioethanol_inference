theme_set(theme_bw())
theme_update(plot.background = element_rect(fill = "transparent", color = "transparent"),
             panel.grid = element_blank(),
             panel.border = element_rect(color = "black", fill = NA, size = 0.5),
             text = element_text(family = "Helvetica", size = 6),
             axis.text = element_text(color = "black"),
             line = element_line(size = 0.17, color = "black"),
             axis.ticks.length = unit(0.75, units = "mm"),
             axis.ticks = element_line(color = "black", size = 0.5))

theme_inset <- theme(text = element_text(family = "Helvetica", size = 5),
                     panel.border = element_rect(color = "black", fill = NA, size = 0.25),
                     axis.ticks.length = unit(0.3, units = "mm"),
                     axis.ticks = element_line(color = "black", size = 0.15))
