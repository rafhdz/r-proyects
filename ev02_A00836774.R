#librerias a descargar
library(BiocManager)

library(Biostrings)

library(ggtree)

library(ggmsa)

library(DECIPHER)

library(BiocParallel)

library(seqinr)

library(adegenet)

library(ape)

library(viridis)

library(RSQLite)

library(viridisLite)

library(ggplot2)

#Código del programa
virus <- c(  "JX869059", "AY508724", "MN908947", "AY390556", "AY278489", "MN985325","AY485277","MT292571")
virus_sequences <- read.GenBank(virus)
write.dna(virus_sequences,  file ="virus_seqs.fasta", format = "fasta", append = FALSE, nbcol = 6, colsep = " ", colw = 10)
virus_seq_not_align <- readDNAStringSet("virus_seqs.fasta", format = "fasta")
virus_seq_not_align <- OrientNucleotides(virus_seq_not_align)
virus_seq_align <- AlignSeqs(virus_seq_not_align)
BrowseSeqs(virus_seq_align, highlight=0)

# Conclusión e interpretación de las graficas
# En base a la pregunta ¿Son muy diferentes las variantes entre cada país?, se puede referenciar a la grafica
# anteriormente mencionada que entre las secuencias de código geneticos de las variantes del virus de comparten 
# 28683 bases, en promedio las 8 sequencias tienen 29835 bases entonces tomando en cuenta ese promedio, se 
# comparten aproximadamente 96% del código genetico de las varientes del virus. Llevando estos datos al usos o
# utilidad de la bioinformatica en relación al sector salud (tocando un virus que causa una enfermedad en las personas).
# En relación a la pregunta, a las varientes del virus no son muy diferentes entre sí aunque es importante mencionar
# que los cambios en el código (aunque sean pequeños) pueden tener un impacto significativo en las caracteristicas que
# contiene cada variente como tener un impacto en la transmisibilidad del virus mediante la modificación de las funciones
# biologicas. De la misma manera, entre más mutaciones se tiene apartir del virus original o las demás varientes de este 
# es más probable que este tenga inmunidad de las vacunas.
# Durante la elaboración de las vacunas para el Covid-19 (de las que se incluyen: de virus inactivo, de acido nucleico, de 
# subunidades de proteinas y de vector adenoviral), las cuales se llevaron a cabo dentro de un tiempo record mediante
# la colaboración entre la comunidad cientifica y el campo de investigación del biocomputo (en el analisis del código 
# genetico del virus), se demostro que los continuos cambios en el código generados por las mutaciones en las varrientes
# dificultaron el desarrollo de estas. Las más significativas basada en el ARN mensajero que emplean el código genetico del
# virus para su funcionamiento puede ser beneficiadas directamente con la biocomputo ya que tiene campo permite el analisis de
# las mutaciones de una manera rápida para crear adaptaciones continuas en base a las mutaciones de las varientes más 
# predominantes demanera global. Permitiendo un seguimiento continuo del virus mediante la adaptación de la vacunación.
# Referencias:
# Zhou, W., & Wang, W. (2021). Fast-spreading SARS-CoV-2 variants: challenges to and new design strategies of COVID-19 vaccines. 
# Signal Transduction and Targeted Therapy, 6(1). 
# https://doi.org/10.1038/s41392-021-00644-x