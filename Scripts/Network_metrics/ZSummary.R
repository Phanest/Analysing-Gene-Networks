ref = 1
test = 2

name = 'Cancer'

statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);

modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];

plotMods = !(modColors %in% c("grey", "gold"));
text = modColors[plotMods];

plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
mains = c("Preservation Median rank", "Preservation Zsummary");

Zsummary.values = mp$preservation$Z[[ref]][[test]][, 2]
Zsummary.values = Zsummary.values[plotMods]

names(Zsummary.values) = text
Zsummary.values = sort(Zsummary.values)

#################################

dev.new()

sizeGrWindow(12,9)

par(mfrow = c(1,2))
par(mar = c(3,3,2,1))
par(mgp = c(1.6, 0.4, 0));

for (p in 1:2){
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  
  if (p==2){
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  }
  else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p],
       xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  
  #labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
  
  # For Zsummary, add threshold lines
  if (p==2){
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  }
}

######################################
# Re-initialize module color labels and sizes
modColors = rownames(statsZ)
moduleSizes = mp$quality$Z[[ref]][[test]][, 1];

plotMods = !(modColors %in% c("grey", "gold"));
labs = match(modColors[plotMods], standardColors(62));

sizeGrWindow(12,9)

par(mfrow = c(4,4))
par(mar = c(3,3,2,1))
par(mgp = c(1.6, 0.4, 0));

for (s in 1:ncol(statsZ)){
  min = min(statsZ[plotMods, s], na.rm = TRUE);
  max = max(statsZ[plotMods, s], na.rm = TRUE);
  
  if (min > -max/12) min = -max/12
  
  plot(moduleSizes[plotMods], statsZ[plotMods, s], col = 1, bg = modColors[plotMods], pch = 21,
       main = colnames(statsZ)[s],
       cex = 2.2,ylab = colnames(statsZ)[s],
       xlab = "Module size",
       log = "x",
       ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min)),
       xlim = c(30, 800),
       cex.lab = 1.2,
       cex.axis = 1.2)
  
  labelPoints(moduleSizes[plotMods], statsZ[plotMods, s], labs, cex = 0.8, offs = 0.06);
  #text(moduleSizes[-1], statsZ[-c(1:2), s], labels = letter[-c(1:2)], col = "black"); #modColors[-2]);
  abline(h=0)
  abline(h=2, col = "blue", lty = 2)
  abline(h=10, col = "darkgreen", lty = 2)
}

moduleColor = data.frame(color = modColors[plotMods], label = labs)

write.csv(moduleColor, sprintf('./moduleColor2Number%s.csv', name))