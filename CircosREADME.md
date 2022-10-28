# Circos is a visualization tool to facilitate the identification and analysis of similarities and differences arising from comparisons of genomes. This tool is effective in displaying variation in genome structure and, generally, any other kind of positional relationships between genomic intervals. Such data are routinely produced by sequence alignments, hybridization arrays, genome mapping, and genotyping studies. Circos uses a circular ideogram layout to facilitate the display of relationships between pairs of positions by the use of ribbons, which encode the position, size, and orientation of related genomic elements. Circos is capable of displaying data as scatter, line, and histogram plots, heat maps, tiles, connectors, and text.
# Initializing the library for Circos
library(RCircos)
data("UCSC.HG38.Human.CytoBandIdeogram")
cyto.info = UCSC.HG38.Human.CytoBandIdeogram
RCircos.Set.Core.Components(cyto.info, 
                            chr.exclude=NULL, 
                            tracks.inside=10, 
                            tracks.outside=0)
                            
# Modify the core components                            
RCircos.Set.Plot.Area() 
RCircos.Chromosome.Ideogram.Plot()

# Customising the outer band
cyto.info = UCSC.HG38.Human.CytoBandIdeogram
cyto.info$Name = NA
cyto.info$Stain = NA
RCircos.Set.Core.Components(cyto.info, 
                            chr.exclude=NULL, 
                            tracks.inside=10, 
                            tracks.outside=0)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
ideo = RCircos.Get.Plot.Ideogram()
ideo$BandColor = 'purple'
num = which(ideo$Chromosome == 'chrX')
ideo[num, 'BandColor'] = 'chartreuse'

num = which(ideo$Chromosome == 'chrY')
ideo[num, 'BandColor'] = 'orange'


RCircos.Reset.Plot.Ideogram(ideo)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

num = which(ideo$Chromosome == 'chr1')
ideo[num, 'ChrColor'] = 'goldenrod2'

RCircos.Reset.Plot.Ideogram(ideo)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

library(biomaRt)

library(org.Hs.eg.db)
matC= readRDS("mat.RDS")
View(matC)
m = useMart('ensembl', dataset='hsapiens_gene_ensembl')
coords = getBM(attributes=c('chromosome_name', 'start_position', 
                            'end_position', 'hgnc_symbol'),
               filters = c('hgnc_symbol'),
               values = list(rownames(mat)),
               mart = m)
View(coords)
write.csv(coords, file = 'coords1.csv')
coords$chromosome_name = paste0('chr', coords$chromosome_name)
chr_order= unique(cyto.info$Chromosome)
coords$chromosome_name = factor(coords$chromosome_name, levels = chr_order)
num = which(is.na(coords$chromosome_name))
coords = coords[-num, ]

up = which((matC$pval < 0.01) &
             (matC$log2FC > 1))
upmat = mat[up, ]
num = which(coords$hgnc_symbol %in% rownames(upmat))
coords1 = coords[num, ]

RCircos.Gene.Name.Plot(coords1, name.col=4, track.num = 2, side = "in",
                       is.sorted = F)
# To see the number of labels for each chromosome.
genes = intersect(rownames(mat), coords$hgnc_symbol)

mat1 = matC[genes, ]
df = cbind.data.frame(rownames(mat1), mat1[, c(1,2,4)])
colnames(df)[1] = 'hgnc_symbol'

data = merge(coords, df, by = 'hgnc_symbol', all.x = T)

data = data[, c('chromosome_name', 'start_position',
                'end_position', 'hgnc_symbol',
                'test', 'control', 'pval','log2FC')]
RCircos.Heatmap.Plot(data, data.col = 7, track.num = 6, side = "in",
                     min.value = -3, max.value = 6, 
                    is.sorted = F)    
 RC.param = RCircos.Get.Plot.Parameters()
RC.param['heatmap.color'] = "GreenWhiteRed"
RCircos.Reset.Plot.Parameters(RC.param)

RCircos.Heatmap.Plot(data, data.col = 7, track.num = 10, side = "in",
                     min.value = -2, max.value = 2,
                     is.sorted = F)                   
