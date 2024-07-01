## ----HTML, include = FALSE----------------------------------------------------
Zc <- "<i>Z</i><sub>C</sub>"
nH2O <- "<i>n</i><sub>H<sub>2</sub>O</sub>"
nO2 <- "<i>n</i><sub>O<sub>2</sub></sub>"
H2O <- "H<sub>2</sub>O"
O2 <- "O<sub>2</sub>"

## ----setup, include = FALSE---------------------------------------------------
oldopt <- options(width = 80)
# Use pngquant to optimize PNG images
library(knitr)
knit_hooks$set(pngquant = hook_pngquant)
pngquant <- "--speed=1 --quality=0-25"
if (!nzchar(Sys.which("pngquant"))) pngquant <- NULL
# To make warnings appear in text box 20230619
# https://selbydavid.com/2017/06/18/rmarkdown-alerts/
knitr::knit_hooks$set(
   error = function(x, options) {
     paste('\n\n<div class="alert alert-danger">',
           gsub('##', '\n', gsub('^##\ Error', '**Error:**', x)),
           '</div>', sep = '\n')
   },
   warning = function(x, options) {
     paste('\n\n<div class="alert alert-warning">',
           gsub('##', '\n', gsub('^##\ Warning:', '**Warning:**', x)),
           '</div>', sep = '\n')
   },
   message = function(x, options) {
     paste('\n\n<div class="alert alert-info">',
           gsub('##', '\n', x),
           '</div>', sep = '\n')
   }
)
## Don't evaluate chunks if phyloseq is not available 20230619
#if(!requireNamespace("phyloseq", quietly = TRUE)) {
#  knitr::opts_chunk$set(eval = FALSE)
#  day <- imin <- AA.RDP <- map.RDP <- map.silva <- NULL
#  warning("The **phyloseq** package is not available, so this vignette shows only the code without the results.")
#}

## ----load_packages------------------------------------------------------------
library(chem16S)
library(phyloseq)
library(ggplot2)
theme_set(theme_bw())

## ----load_mouse.silva---------------------------------------------------------
data(mouse.silva)
mouse.silva

## ----ps_dada2_barplot, fig.width = 7, fig.height = 5, fig.align = "center", pngquant = pngquant----
top20.silva <- names(sort(taxa_sums(mouse.silva), decreasing = TRUE))[1:20]
mouse.silva.top20 <- transform_sample_counts(mouse.silva, function(OTU) OTU/sum(OTU))
mouse.silva.top20 <- prune_taxa(top20.silva, mouse.silva.top20)
plot_bar(mouse.silva.top20, x = "Day", fill = "Family") + facet_wrap(~When, scales = "free_x")

## ----mouse.silva_taxacounts---------------------------------------------------
tc.silva <- ps_taxacounts(mouse.silva)
head(tc.silva[, -1])

## ----mouse.silva_taxid--------------------------------------------------------
head(tc.silva$taxid)

## ----mouse.silva_levels, collapse = TRUE--------------------------------------
table(tc.silva$rank)

## ----map.silva, collapse = TRUE-----------------------------------------------
map.silva <- map_taxa(tc.silva, refdb = "RefSeq_206")

## ----map.RDP, collapse = TRUE-------------------------------------------------
data(mouse.RDP)
tc.RDP <- ps_taxacounts(mouse.RDP)
map.RDP <- map_taxa(tc.RDP, refdb = "RefSeq_206")

## ----AA.RDP, collapse = TRUE--------------------------------------------------
AA.RDP <- ps_metrics(mouse.RDP, refdb = "RefSeq_206", quiet = TRUE, return_AA = TRUE)
head(AA.RDP)

## ----AA.RDP_length, fig.width = 7, fig.height = 5, fig.align = "center", pngquant = pngquant----
lengths <- canprot::calc_metrics(AA.RDP, "length")[, 1]
hist(lengths)
imin <- which.min(lengths)
text(lengths[imin], 1.5, AA.RDP$Run[imin], adj = 1)

## ----metrics.RDP, collapse = TRUE---------------------------------------------
metrics.RDP <- ps_metrics(mouse.RDP, refdb = "RefSeq_206", quiet = TRUE)
head(metrics.RDP)

## ----cor.RDP, collapse = TRUE-------------------------------------------------
cor(metrics.RDP$Zc, metrics.RDP$nO2)
cor(metrics.RDP$Zc, metrics.RDP$nH2O)

## ----plot_metrics.RDP, fig.width = 7, fig.height = 5, fig.align = "center", pngquant = pngquant----
plot_ps_metrics(mouse.RDP, x = "Day", color = "When", shape = "When", refdb = "RefSeq_206", quiet = TRUE) +
  geom_point(size = 3)

## ----plot_metrics2.RDP, fig.width = 6, fig.height = 4, fig.align = "center", pngquant = pngquant----
plot_ps_metrics2(mouse.RDP, color = "When", shape = "When", refdb = "RefSeq_206", quiet = TRUE) +
  geom_point(size = 3)

## ----data.GTDB, collapse = TRUE, fig.width = 6, fig.height = 4, fig.align = "center", pngquant = pngquant----
data(mouse.GTDB_214)
plot_ps_metrics2(mouse.GTDB_214, refdb = "GTDB_214", color = "When", shape = "When") + geom_point(size = 3)

## ----early.GTDB, collapse = TRUE, fig.width = 6, fig.height = 4, fig.align = "center", pngquant = pngquant----
metrics.GTDB <- ps_metrics(mouse.GTDB_214)
is.early <- sample_data(mouse.GTDB_214)$When == "Early"
iout <- which.min(metrics.GTDB[is.early, ]$nH2O)
(day <- sample_data(mouse.GTDB_214)[is.early, ]$Day[iout])

## ----barplot.GTDB, fig.width = 7, fig.height = 5, fig.align = "center", pngquant = pngquant----
top20.GTDB <- names(sort(taxa_sums(mouse.GTDB_214), decreasing = TRUE))[1:20]
mouse.GTDB_214.top20 <- transform_sample_counts(mouse.GTDB_214, function(OTU) OTU/sum(OTU))
mouse.GTDB_214.top20 <- prune_taxa(top20.GTDB, mouse.GTDB_214.top20)
plot_bar(mouse.GTDB_214.top20, x = "Day", fill = "Phylum") + facet_wrap(~When, scales = "free_x")

## ----cleanup, include = FALSE-------------------------------------------------
options(oldopt)

