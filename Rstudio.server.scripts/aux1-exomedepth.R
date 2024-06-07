require("ExomeDepth")
require("RIdeogram")

exomeDepth <- function(my.test.num=NULL, outdir=NULL, sample=NULL, ExomeCount, my.ref.samples, exons.mm10.GRanges){
  # perform the ExomeDepth calculations to generate all.exons
  ExomeCount.dafr <- as(ExomeCount, 'data.frame')
  
  # reformat bam file names to match R naming conversion
  my.ref.samples <- gsub('-','.',my.ref.samples)
  my.test.num <- gsub('-','.',my.test.num)
  my.reference.set <- as.matrix(ExomeCount.dafr[, basename(my.ref.samples)])
  my.test <- ExomeCount.dafr[,my.test.num]
  my.choice <- select.reference.set(test.counts = my.test,
                                     reference.counts = my.reference.set,
                                     bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start)/1000,
                                     n.bins.reduced = 10000
                                    )
  
  my.matrix <- as.matrix( ExomeCount.dafr[, my.choice$reference.choice, drop = FALSE])
  
  my.reference.selected <- apply(X = my.matrix,
                                 MAR = 1,
                                 FUN = sum
                                )
  
  all.exons <- new('ExomeDepth',
                   test = my.test,
                   reference = my.reference.selected,
                   formula = 'cbind(test, reference) ~ 1'
  )
  
  #ExomeCount.dafr$chromosome <- gsub(as.character(ExomeCount.dafr$chromosome),pattern = 'chr',replacement = '')
  
  
  all.exons <- CallCNVs(x = all.exons,
                        transition.probability = 10^-4,
                        chromosome = ExomeCount.dafr$chromosome,
                        start = ExomeCount.dafr$start,
                        end = ExomeCount.dafr$end,
                        name = as.character(exons.mm10.GRanges$exon_id)
                        #name = exons.mm10.GRanges$names
                        )
  

  all.exons <- AnnotateExtra(x = all.exons,
                             reference.annotation = exons.mm10.GRanges,
                             min.overlap = 0.0001,
                             column.name = 'genes.mm10'
  )

  
  # a %>% stringi::stri_split(regex = ",") %>% unlist() %>% unique() %>% stringr::str_flatten_comma()
  collapsed_gene_names <- sapply(all.exons@CNV.calls[,"genes.mm10"], function(x) {
                                          stringi::stri_split(str = x,regex = ",") %>%
                                          unlist() %>% 
                                          unique() %>% 
                                          stringr::str_flatten_comma() 
                                        }
                                      )
  names(collapsed_gene_names) <- NULL
  
  all.exons@CNV.calls[,"genes.mm10"] <- collapsed_gene_names
  
  # output the data
  write.csv(file = paste0(outdir,"/exomeDepth_",sample,".csv"),
            x = all.exons@CNV.calls,
            row.names = FALSE
  )
  
  assign("all.exons", all.exons, .GlobalEnv)
  
  # make the ideogram displaying only those with BF values above the mean
  
  mm10data <- data.frame(V1=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY"), 
                         V2=1, 
                         V3=c(195471971, 182113224,160039680,156508116,151834684,149736546,145441459,129401213,124595110,130694993,122082543,120129022,120421639,124902244,104043685,98207768,94987271,90702639,61431566,171031299,91744698)
                          )
  
  colnames(mm10data) <- c("Chr","Start","End")
  mm10data$Chr <- as.character(mm10data$Chr)
  
  chrdata <- data.frame(V1=all.exons@CNV.calls$chromosome,
                        V2=all.exons@CNV.calls$start,
                        V3=all.exons@CNV.calls$end,
                        V4=all.exons@CNV.calls$BF
  )
  
  colnames(chrdata) <- c("Chr","Start","End","Value")
  overmean <- subset(all.exons@CNV.calls, BF >= mean(BF))
  overmean <- overmean[,c(3,5,6,7)]
  colnames(overmean) <- c("Type","Start","End","Chr")
  overmean$Type <- ifelse(overmean$Type == "duplication", "Gain", "Loss")
  overmean$Shape <- "triangle"
  overmean$color <- ifelse(overmean$Type == "Gain", "ff0000", "008000")
  overmean <- overmean[,c(1,5,4,2,3,6)]
  ideogram(karyotype=mm10data,overlaid=chrdata,label=overmean,label_type="marker")
  convertSVG("chromosome.svg", device = "png")
  file.copy(from = "./chromosome.svg", to = paste0("./",outdir,"/", "chromosome_",sample,".svg"), overwrite = TRUE)
  file.copy(from = "./chromosome.png", to = paste0("./",outdir,"/", "chromosome_",sample,".png"), overwrite = TRUE)
}
