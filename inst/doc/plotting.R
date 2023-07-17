## ----HTML, include = FALSE----------------------------------------------------
Zc <- "<i>Z</i><sub>C</sub>"
nH2O <- "<i>n</i><sub>H<sub>2</sub>O</sub>"

## ----setup, include = FALSE---------------------------------------------------
oldopt <- options(width = 80)
## Use pngquant to optimize PNG images
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
#  warning("The **phyloseq** package is not available, so this vignette shows only the code without the results.")
#}

## ----load_packages, message = FALSE-------------------------------------------
require(chem16S)
require(phyloseq)
require(ggplot2)
theme_set(theme_bw())

# For composing plots and making a common legend (plot_layout())
require(patchwork)

# For annotating plots with regression coefficients (stat_poly_line())
require(ggpmisc)

## ----load_GlobalPatterns------------------------------------------------------
data(GlobalPatterns)
Human = get_variable(GlobalPatterns, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue")
sample_data(GlobalPatterns)$Human <- factor(Human)

## ----plot_metrics2.GlobalPatterns, collapse = TRUE----------------------------
p2 <- plot_ps_metrics2(GlobalPatterns, color = "SampleType", shape = "Human", refdb = "RefSeq")

## ----plot_metrics2.GlobalPatterns_geom, fig.width = 6, fig.height = 4, fig.align = "center", pngquant = pngquant----
p2 + geom_polygon(aes(fill = SampleType), alpha = 0.5) + geom_point(size = 3) +
  guides(colour = guide_legend(override.aes = list(shape = c(17, 19, 19, 17, 19, 19, 17, 19, 17))))

## ----read_Humboldt------------------------------------------------------------
psfile <- system.file("extdata/DADA2/FEN+22/ps_FEN+22.rds", package = "chem16S")
ps <- readRDS(psfile)
ps <- prune_samples(sample_names(ps) != "SRR1346095", ps)
ps

## ----plot_metrics2.Humboldt, fig.width = 6, fig.height = 4, fig.align = "center", pngquant = pngquant----
plot_ps_metrics2(ps, color = "Location") +
  geom_polygon(aes(fill = Location), alpha = 0.5) + geom_point(size = 3)

## ----check_patchwork, echo = FALSE--------------------------------------------
# Don't evaluate remaining chunks if patchwork is not available 20230627
if(!requireNamespace("patchwork", quietly = TRUE)) {
  knitr::opts_chunk$set(eval = FALSE)
  warning("The **patchwork** package is not available, so the remaining plots are not shown.")
}

## ----plot_metrics2.Humboldt_redox_OM, fig.width = 6, fig.height = 5, fig.align = "center", pngquant = pngquant----
p2 <- plot_ps_metrics2(ps, color = "Sediment_redox") +
  geom_point(size = 4) + labs(color = "Sediment redox (Eh)")
p3 <- plot_ps_metrics2(ps, color = "Organic_matter") +
  geom_point(size = 4) + labs(color = "Organic matter (%)")
p2 / p3

## ----sample.data.and.chemical.metrics.for.communities-------------------------
sample.data.and.chemical.metrics.for.communities <- cbind(sample_data(ps), ps_metrics(ps))

## ----check_ggpmisc, echo = FALSE----------------------------------------------
# Don't evaluate remaining chunks if ggpmisc is not available 20230627
if(!requireNamespace("ggpmisc", quietly = TRUE)) {
  knitr::opts_chunk$set(eval = FALSE)
  warning("The **ggpmisc** package is not available, so the remaining plots are not shown.")
}

## ----scatter_plot, eval = FALSE-----------------------------------------------
#  # Defuse (enquo) and Inject (!!) from https://www.tidyverse.org/blog/2018/07/ggplot2-tidy-evaluation/
#  # Regression line and equation from https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph
#  scatter_plot <- function(data = sample.data.and.chemical.metrics.for.communities, x, y, xlab, ylab) {
#    x <- enquo(x)
#    y <- enquo(y)
#    ggplot(data, aes(x = !!x, y = !!y, color = .data[["Location"]])) +
#      geom_point() + xlab(xlab) + ylab(ylab) +
#      # Override aes to plot one regression line for samples from all locations
#      stat_poly_line(aes(x = !!x, y = !!y), inherit.aes = FALSE) +
#      stat_poly_eq(aes(x = !!x, y = !!y), inherit.aes = FALSE, label.x = "center")
#  }
#  
#  sp1 <- scatter_plot(x = Sediment_redox, y = Zc, xlab = "Sediment redox (mV)", ylab = chemlab("Zc"))
#  sp2 <- scatter_plot(x = Sediment_redox, y = nH2O, xlab = "Sediment redox (mV)", ylab = chemlab("nH2O"))
#  sp3 <- scatter_plot(x = Organic_matter, y = Zc, xlab = "Organic matter (%)", ylab = chemlab("Zc"))
#  sp4 <- scatter_plot(x = Organic_matter, y = nH2O, xlab = "Organic matter (%)", ylab = chemlab("nH2O"))
#  sp1 + sp2 + sp3 + sp4 + plot_layout(guides = "collect")

## ----scatter_plot, echo = FALSE, fig.width = 7, fig.height = 6, fig.align = "center", pngquant = pngquant----
# Defuse (enquo) and Inject (!!) from https://www.tidyverse.org/blog/2018/07/ggplot2-tidy-evaluation/
# Regression line and equation from https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph
scatter_plot <- function(data = sample.data.and.chemical.metrics.for.communities, x, y, xlab, ylab) {
  x <- enquo(x)
  y <- enquo(y)
  ggplot(data, aes(x = !!x, y = !!y, color = .data[["Location"]])) +
    geom_point() + xlab(xlab) + ylab(ylab) +
    # Override aes to plot one regression line for samples from all locations
    stat_poly_line(aes(x = !!x, y = !!y), inherit.aes = FALSE) +
    stat_poly_eq(aes(x = !!x, y = !!y), inherit.aes = FALSE, label.x = "center")
}

sp1 <- scatter_plot(x = Sediment_redox, y = Zc, xlab = "Sediment redox (mV)", ylab = chemlab("Zc"))
sp2 <- scatter_plot(x = Sediment_redox, y = nH2O, xlab = "Sediment redox (mV)", ylab = chemlab("nH2O"))
sp3 <- scatter_plot(x = Organic_matter, y = Zc, xlab = "Organic matter (%)", ylab = chemlab("Zc"))
sp4 <- scatter_plot(x = Organic_matter, y = nH2O, xlab = "Organic matter (%)", ylab = chemlab("nH2O"))
sp1 + sp2 + sp3 + sp4 + plot_layout(guides = "collect")

## ----load_Qinghai.Tibet-------------------------------------------------------
psfile2 <- system.file("extdata/DADA2/ZFZ+23/ps_ZFZ+23.rds", package = "chem16S")
ps2 <- readRDS(psfile2)
data.and.metrics <- cbind(sample_data(ps2), ps_metrics(ps2))
ps2

## ----scatter_plot_2, eval = FALSE---------------------------------------------
#  scatter_plot_2 <- function(data = data.and.metrics, x, y, xlab, ylab) {
#    x <- enquo(x)
#    y <- enquo(y)
#    ggplot(data, aes(x = !!x, y = !!y)) +
#      geom_point() + xlab(xlab) + ylab(ylab) +
#      stat_poly_line() +
#      stat_poly_eq(label.x = "center")
#  }
#  sp1 <- scatter_plot_2(x = ORP, y = Zc, xlab = "ORP (mV)", ylab = chemlab("Zc"))
#  sp2 <- scatter_plot_2(x = ORP, y = nH2O, xlab = "ORP (mV)", ylab = chemlab("nH2O"))
#  sp3 <- scatter_plot_2(x = T, y = Zc, xlab = "T (°C)", ylab = chemlab("Zc"))
#  sp4 <- scatter_plot_2(x = T, y = nH2O, xlab = "T (°C)", ylab = chemlab("nH2O"))
#  sp1 + sp2 + sp3 + sp4

## ----scatter_plot_2, echo = FALSE, fig.width = 7, fig.height = 6, fig.align = "center", pngquant = pngquant----
scatter_plot_2 <- function(data = data.and.metrics, x, y, xlab, ylab) {
  x <- enquo(x)
  y <- enquo(y)
  ggplot(data, aes(x = !!x, y = !!y)) +
    geom_point() + xlab(xlab) + ylab(ylab) +
    stat_poly_line() +
    stat_poly_eq(label.x = "center")
}
sp1 <- scatter_plot_2(x = ORP, y = Zc, xlab = "ORP (mV)", ylab = chemlab("Zc"))
sp2 <- scatter_plot_2(x = ORP, y = nH2O, xlab = "ORP (mV)", ylab = chemlab("nH2O"))
sp3 <- scatter_plot_2(x = T, y = Zc, xlab = "T (°C)", ylab = chemlab("Zc"))
sp4 <- scatter_plot_2(x = T, y = nH2O, xlab = "T (°C)", ylab = chemlab("nH2O"))
sp1 + sp2 + sp3 + sp4

## ----merge_datasets-----------------------------------------------------------
taxa_names(ps2) <- paste0(taxa_names(ps2), "b")
ps_merged <- merge_phyloseq(ps, ps2)
ps_merged

## ----plot_metrics2.merged, fig.width = 6, fig.height = 4, fig.align = "center", pngquant = pngquant----
sample_data(ps_merged)$Environment <-
  ifelse(is.na(sample_data(ps_merged)$Depth), "Hot spring", "Marine sediment")
plot_ps_metrics2(ps_merged, color = "Environment", shape = "Environment") + geom_point(size = 3)

## ----cleanup, include = FALSE-------------------------------------------------
options(oldopt)

