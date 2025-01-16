## ----setup, include = FALSE---------------------------------------------------
oldopt <- options(width = 80)
## Use pngquant to optimize PNG images
library(knitr)
knit_hooks$set(pngquant = hook_pngquant)
pngquant <- "--speed=1 --quality=0-25"
if (!nzchar(Sys.which("pngquant"))) pngquant <- NULL

## ----HTML, include = FALSE----------------------------------------------------
Zc <- "<i>Z</i><sub>C</sub>"
nH2O <- "<i>n</i>H<sub>2</sub>O"

## ----load_packages, message = FALSE-------------------------------------------
library(chem16S)

## ----taxon_AA-----------------------------------------------------------------
taxon_AA <- read.csv(system.file("RefDB/RefSeq_206/taxon_AA.csv.xz", package = "chem16S"))
ranks <- taxon_AA$protein
table(ranks)[unique(ranks)]

## ----Zc_boxplot, fig.width = 5, fig.height = 5, fig.align = "center", pngquant = pngquant----
taxon_Zc <- canprot::calc_metrics(taxon_AA, "Zc")[, 1]
Zc_list <- sapply( unique(ranks), function(rank) taxon_Zc[ranks == rank] )
opar <- par(mar = c(4, 7, 1, 1))
boxplot(Zc_list, horizontal = TRUE, las = 1, xlab = chemlab("Zc"))
par(opar)

## ----phylum_to_genus----------------------------------------------------------
taxnames <- read.csv(system.file("RefDB/RefSeq_206/taxonomy.csv.xz", package = "chem16S"))
phylum_to_genus <- function(phylum) na.omit(unique(taxnames$genus[taxnames$phylum == phylum]))
get_Zc <- function(genera) na.omit(taxon_Zc[match(genera, taxon_AA$organism)])
sapply(sapply(sapply(c("Crenarchaeota", "Euryarchaeota"), phylum_to_genus), get_Zc), mean)

## ----class_to_genus-----------------------------------------------------------
class_to_genus <- function(class) na.omit(unique(taxnames$genus[taxnames$class == class]))
sapply(sapply(sapply(c("Methanococci", "Halobacteria"), class_to_genus), get_Zc), mean)

## ----phylum_Zc, fig.width = 7, fig.height = 6, fig.align = "center", pngquant = pngquant----
taxnames2 <- taxnames[taxnames$superkingdom != "Viruses", ]
taxnames3 <- taxnames2[!duplicated(taxnames2$genus), ]
(top20_phyla <- head(sort(table(taxnames3$phylum), decreasing = TRUE), 20))
Zc_list <- sapply(sapply(names(top20_phyla), phylum_to_genus), get_Zc)
order_Zc <- order(sapply(Zc_list, mean))
Zc_list <- Zc_list[order_Zc]
opar <- par(mar = c(4, 13, 1, 1))
boxplot(Zc_list, horizontal = TRUE, las = 1, xlab = chemlab("Zc"))
par(opar)

## ----class_Zc, fig.width = 10, fig.height = 5, fig.align = "center", pngquant = pngquant----
opar <- par(mfrow = c(1, 2), mar = c(4, 10, 1, 1))
for(phylum in c("Euryarchaeota", "Proteobacteria")) {
  taxnames4 <- taxnames3[taxnames3$phylum == phylum, ]
  classes <- na.omit(unique(taxnames4$class))
  Zc_list <- sapply(sapply(classes, class_to_genus), get_Zc)
  order_Zc <- order(sapply(Zc_list, mean))
  Zc_list <- Zc_list[order_Zc]
  boxplot(Zc_list, horizontal = TRUE, las = 1, xlab = chemlab("Zc"))
}
par(opar)

## ----class_nH2O, fig.width = 10, fig.height = 5, fig.align = "center", pngquant = pngquant----
taxon_nH2O <- canprot::calc_metrics(taxon_AA, "nH2O")[, 1]
get_nH2O <- function(genera) na.omit(taxon_nH2O[match(genera, taxon_AA$organism)])

opar <- par(mfrow = c(1, 2), mar = c(4, 10, 1, 1))
for(phylum in c("Euryarchaeota", "Proteobacteria")) {
  taxnames4 <- taxnames3[taxnames3$phylum == phylum, ]
  classes <- na.omit(unique(taxnames4$class))
  nH2O_list <- sapply(sapply(classes, class_to_genus), get_nH2O)
  order_nH2O <- order(sapply(nH2O_list, mean))
  nH2O_list <- nH2O_list[order_nH2O]
  boxplot(nH2O_list, horizontal = TRUE, las = 1, xlab = chemlab("nH2O"))
}
par(opar)

## ----viruses_top50------------------------------------------------------------
(top50_phyla <- head(sort(table(taxnames$phylum), decreasing = TRUE), 50))

## ----viruses_plot, fig.width = 7, fig.height = 5, fig.align = "center", pngquant = pngquant----
Zc_mean <- sapply(sapply(sapply(names(top50_phyla), phylum_to_genus), get_Zc), mean)
nH2O_mean <- sapply(sapply(sapply(names(top50_phyla), phylum_to_genus), get_nH2O), mean)
domain <- taxnames$superkingdom[match(names(top50_phyla), taxnames$phylum)]
pchs <- c(24, 21, 23)
pch <- sapply(domain, switch, Archaea = pchs[1], Bacteria = pchs[2], Viruses = pchs[3])
bgs <- topo.colors(3, alpha = 0.5)
bg <- sapply(domain, switch, Archaea = bgs[1], Bacteria = bgs[2], Viruses = bgs[3])
opar <- par(mar = c(4, 4, 1, 1))
plot(Zc_mean, nH2O_mean, xlab = chemlab("Zc"), ylab = chemlab("nH2O"), pch = pch, bg = bg)
ilow <- nH2O_mean < -0.77 & domain == "Bacteria"
xadj <- c(-0.9, -0.8, 0.8, 1, -0.8)
yadj <- c(0, 1, 1, -1, -1)
text(Zc_mean[ilow] + 0.02 * xadj, nH2O_mean[ilow] + 0.005 * yadj, names(top50_phyla[ilow]), cex = 0.9)
legend("bottomleft", c("Archaea", "Bacteria", "Viruses"), pch = pchs, pt.bg = bgs)
par(opar)

## ----other_metrics, fig.width = 8, fig.height = 5, fig.align = "center", pngquant = pngquant----
AAcomp <- taxon_AA[match(classes, taxon_AA$organism), ]
metrics <- canprot::calc_metrics(AAcomp, c("HC", "OC", "NC", "SC", "GRAVY", "pI", "MW", "plength"))
layout(rbind(c(1, 2, 5), c(3, 4, 5)), widths = c(2, 2, 1.5))
opar <- par(mar = c(4.5, 4, 1, 1), cex = 1)
plot(metrics$OC, metrics$HC, col = 1:10, pch = 1:10, xlab = "O/C", ylab = "H/C")
plot(metrics$NC, metrics$SC, col = 1:10, pch = 1:10, xlab = "N/C", ylab = "S/C")
plot(metrics$pI, metrics$GRAVY, col = 1:10, pch = 1:10, xlab = "pI", ylab = "GRAVY")
plot(metrics$plength, metrics$MW, col = 1:10, pch = 1:10, xlab = "Length", ylab = "MW")
plot.new()
legend("right", classes, col = 1:10, pch = 1:10, bty = "n", xpd = NA)
par(opar)

## ----cleanup, include = FALSE-------------------------------------------------
options(oldopt)

