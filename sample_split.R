# 
# Example run
# Rscript sample_split.R -i filtered/ -m filtered/samples.tsv --pdf plots.pdf -o filtered/per_sample


#######################################
# Inputs
#######################################

library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(optparse)


option_list = list(
  make_option(c("-i", "--indir"), type="character", default=NULL, 
              help="Input directory", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL,
              help="Output directory", metavar="character"),
  make_option(c("--pdf"), type="character", default=NULL,
              help="Path to output PDF", metavar="character"),
  make_option(c("-m", "--metadata"), type="character", default=NULL,
              help="Tsv file with barcode to sample information", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
params = parse_args(opt_parser);


if(!dir.exists(params$outdir)){
  dir.create(params$outdir)
}

plotfile = params$pdf

#######################################
# Read data
#######################################

meta = read.table(params$metadata, sep = "\t", header = T)
meta$SampleRep = paste(meta$Sample, meta$Replicate, sep = "_rep")
cat("Read metadata with:", dim(meta), "entries\n")

mtx = ReadMtx(mtx = file.path(params$indir, "matrix.mtx"), 
                cells = file.path(params$indir, "barcodes.tsv"), 
                features = file.path(params$indir, "features.tsv"),
                cell.column = 1, 
                feature.column = 2)
cat("Read matrix with:", dim(mtx), "entries\n")

cell2bc = substr(colnames(mtx),1,10)

stopifnot(!is.na(match(cell2bc, meta$Barcode)))

names(cell2bc) = meta$Sample[match(cell2bc, meta$Barcode)]

meta.cell = data.frame(FullBC = colnames(mtx), CellBC = cell2bc,
                       Sample = meta$Sample[match(cell2bc, meta$Barcode)],
                       SampleReplicate = meta$SampleRep[match(cell2bc, meta$Barcode)], 
                       Replicate = factor(meta$Replicate[match(cell2bc, meta$Barcode)]),
                       row.names = colnames(mtx))   


sdata = CreateSeuratObject(mtx, meta.data = meta.cell)
sdata = SetIdent(sdata, value ="Sample") 

saveRDS(sdata, file.path(params$outdir, paste0("allsamples.rds")))

t = table(sdata$Sample)
cat("Cells per sample:\n", paste(names(t),": ",t,"\n"))

#######################################
# Plot
#######################################

pdf(plotfile,10,6)

# stats on number of cells per sample/replicate
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
replicate_colors = sample(col_vector, length(unique(meta$Replicate)))

ggplot(meta.cell, aes(x=Sample,fill=Replicate)) + geom_bar( stat = "count")  + 
  RotatedAxis() + theme_classic() + 
  scale_fill_manual(values = replicate_colors) + 
  ggtitle("Number of cells per sample and replicate") + RotatedAxis()

ggplot(meta.cell, aes(x=Sample,fill=Replicate)) + geom_bar( stat = "count", position = "fill")  + 
  RotatedAxis() + theme_classic() + 
  scale_fill_manual(values = replicate_colors) + 
  ggtitle("Proportion of replicate per sample") + RotatedAxis()


# QC stats
VlnPlot(sdata, "nFeature_RNA", group.by = "Sample", pt.size = .1) + NoLegend()
VlnPlot(sdata, "nCount_RNA", group.by = "Sample", pt.size = .1) + NoLegend()

cat("Range for nFeatures:", range(sdata$nFeature_RNA), "\n")
cat("Range for nCounts:", range(sdata$nCount_RNA), "\n")


# Gene QC stats

# Boxplot with top expressed genes from a Seurat object.

top_expressed_genes = function(sobject, nPlot=20, title="Top expressed genes", assay = "RNA"){
  C = sobject@assays[[assay]]@layers$counts
  C@x = C@x/rep.int(colSums(C), diff(C@p))
  most_expressed <- order(Matrix::rowSums(C), decreasing = T)[nPlot:1]
  boxplot(t(as.matrix(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "Fraction counts per cell",
    col = (scales::hue_pal())(20)[20:1], horizontal = TRUE, main=title)
  invisible(rownames(C)[most_expressed])
}

par(mar=c(5,6,2,2))
top_expressed_genes(sdata)

dev.off()

#######################################
# Write as seurat objects per sample
#######################################

sobjects = SplitObject(sdata, split.by = "Sample")

for(so in names(sobjects)){
  outfile = file.path(params$outdir, paste0(so, ".rds"))
  saveRDS(sobjects[so], outfile)
}
