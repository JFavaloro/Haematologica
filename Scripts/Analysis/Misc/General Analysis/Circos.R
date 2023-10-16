#Circos plots

# Create an empty list to store results
circos_ready <- list()

#1.1 Extract V & J name data & Count from each sample
for(i in 1:length(immdata$data)) {
  sample = immdata$data[[i]] %>%
    group_by(V.name, J.name) %>%
    summarise (count = sum(Clones))
  circos_ready[[i]] = sample[order(-sample$count),]
}

# Fix the names
names(circos_ready) <- names(immdata$data)

# Split the circos_ready list into BM and PB lists
circos_BM = circos_ready[1:4]
circos_PB = circos_ready[5:8]

# Assign colour palette for top 10 clones
vcol <- c("#CC6600", "#6600CC", "#00FFFF", "#CC0000", "#CC0066", 
          "#66CC00", "#00CC66", "#0066CC", "#00CC00", "#CC00CC")
jcol <- c("#FF56C6", "#AA56FF", "#5672FF", "#56BBFF", "#56FFDD", 
          "#56FF7E", "#E3FF56", "#FFC156", "#FF7856", "#0415AC")
colours <- c(vcol, jcol)

# Clean the environment
rm(vcol, jcol, sample)

#1.3 Select the top n pairs
x = circos_BM$BM43 %>% head(n = 100)
y = circos_PB$PB43 %>% head(n = 100)

#Circos plots
#Give colours to common genes across both BM and PB
x_v = intersect(x$V.name, y$V.name)
y_j = intersect(x$J.name, y$J.name)

names(colours) = c(x_v[c(1,2,3,4,5,6,7,8,9,10)], y_j[c(1,2,3,4,5,6,7,8,9,10)])

chordDiagram(x, grid.col = colours, annotationTrack = "grid", preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(x))))))
circos.track(track.index = 1, panel.fun = function(x,y){
  xplot = get.cell.meta.data("xplot")
  if(abs(xplot[2] - xplot[1]) > 3.5)
  {circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
               facing = "clockwise", niceFacing = TRUE, adj = c(0,0.5), cex = 1)}
  else if (abs(xplot[2] - xplot[1]) <= 3.5 &abs(xplot[2] - xplot[1]) >  1)
  {circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
               facing = "clockwise", niceFacing = TRUE, adj = c(0,0.5), cex = 0.75)}
  else
  {circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
               facing = "clockwise", niceFacing = TRUE, adj = c(0,0.5), col = "#606060", cex = 0.5)}
},bg.border = NA)


