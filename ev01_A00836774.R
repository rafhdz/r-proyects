#Rafael Lorenzo Hernández Zambrano A00836774
# librerias a usar 
library(seqinr)
library(ape)
# leer fastas
var_B117 = read.fasta("B.1.1.7.fasta")
var_B1351 = read.fasta("B.1.351.fasta")
var_B1525 = read.fasta("B.1.525.fasta")
var_P1 = read.fasta("P.1.fasta")
var_P2 = read.fasta("P.2.fasta")
# calcular longitudes 
length(var_B117[["HG999939.1"]])
length(var_B1351[["FR988030.1"]])
length(var_B1525[["MW809054.1"]])
length(var_P1[["MW892185.1"]])
length(var_P2[["MW966398.1"]])
# grafica bases ADN 
var1_b = print(count(var_B117[["HG999939.1"]],1)) 
var2_b = print(count(var_B1351[["FR988030.1"]],1))
var3_b = print(count(var_B1525[["MW809054.1"]],1))
var4_b = print(count(var_P1[["MW892185.1"]],1))
var5_b = print(count(var_P2[["MW966398.1"]],1))
var_n = c("A","C","G","T")
var1_bc = c(8375/length(var_B117[["HG999939.1"]])*100,5115/length(var_B117[["HG999939.1"]])*100,5501/length(var_B117[["HG999939.1"]])*100,8974/length(var_B117[["HG999939.1"]])*100)
var2_bc = c(8735/length(var_B1351[["FR988030.1"]])*100, 5350/length(var_B1351[["FR988030.1"]])*100, 5750/length(var_B1351[["FR988030.1"]])*100, 9363/length(var_B1351[["FR988030.1"]])*100)
var3_bc = c(8871/length(var_B1525[["MW809054.1"]])*100, 5458/length(var_B1525[["MW809054.1"]])*100, 5837/length(var_B1525[["MW809054.1"]])*100, 9559/length(var_B1525[["MW809054.1"]])*100)
var4_bc = c(8897/length(var_P1[["MW892185.1"]])*100, 5464/length(var_P1[["MW892185.1"]])*100, 5839/length(var_P1[["MW892185.1"]])*100, 9577/length(var_P1[["MW892185.1"]])*100)
var5_bc = c(8953/length(var_P2[["MW966398.1"]])*100, 5483/length(var_P2[["MW966398.1"]])*100, 5856/length(var_P2[["MW966398.1"]])*100, 9606/length(var_P2[["MW966398.1"]])*100)
barplot(var1_bc, names.arg = var_n, xlab = "Bases de ADN", ylab = "Numero de base en (%)", main = "Comparación de bases de ADN B117 ")
barplot(var2_bc, names.arg = var_n, xlab = "Bases de ADN", ylab = "Numero de base en (%)", main = "Comparación de bases de ADN B1351 ")
barplot(var3_bc, names.arg = var_n, xlab = "Bases de ADN", ylab = "Numero de base en (%)", main = "Comparación de bases de ADN B1525 ")
barplot(var4_bc, names.arg = var_n, xlab = "Bases de ADN", ylab = "Numero de base en (%)", main = "Comparación de bases de ADN P1 ")
barplot(var5_bc, names.arg = var_n, xlab = "Bases de ADN", ylab = "Numero de base en (%)", main = "Comparación de bases de ADN P2 ")
# interpretación de las graficas: debido a que los numeros de las longitudes 
# y las proporciónes de las bases del ADN son bastantes similares, de la misma
# manera son las graficas creadas en relación a los porcentages y cantidades 
# determinadas
# en relación a los porcentages que se presentan la base A y T tienen aproximadamente
# 30 % de presencia (teniendo la T más que la A) y por el otro lado las base
# C y G tienen un porcentage de 16 % a 17 % respectivamente por lo que la distribucción
# de la cantidad de las bases no es igualitaria en todas las variantes
# GC de cada variante 
var1_c = print(count(var_B117[["HG999939.1"]],2)) 
var2_c = print(count(var_B1351[["FR988030.1"]],2))
var3_c = print(count(var_B1525[["MW809054.1"]],2))
var4_c = print(count(var_P1[["MW892185.1"]],2))
var5_c = print(count(var_P2[["MW966398.1"]],2))
cant_gc = c(1089,1147,1162,1166,1167)
var_gcn = c("B117","B1351","B1525","P1","P2")
barplot(cant_gc, names.arg = var_gcn, xlab = "Variantes", ylab = "Cantidad de GC", main = "Cantidad de GC")
# interpretación de la grafica: todas las variantes tienen una cantidad similar aproximadamente 1150 GC de su código 
# genetico, teniendo la variante B117 una cantidad notablemente menor de GC en comparación al resto que presentan una menor 
# diferencia
# secuencias contrasentido
hebrac_invertir = function(hebra){
  sv <- unlist(strsplit(hebra, ""))
  si <- rev(sv)
  si <- paste(si, collapse = "")
  return(si)
}
var1_cs = print(hebrac_invertir(var_B117[["HG999939.1"]])) 
var2_cs = print(hebrac_invertir(var_B1351[["FR988030.1"]]))
var3_cs = print(hebrac_invertir(var_B1525[["MW809054.1"]]))
var4_cs = print(hebrac_invertir(var_P1[["MW892185.1"]]))
var5_cs = print(hebrac_invertir(var_P2[["MW966398.1"]]))