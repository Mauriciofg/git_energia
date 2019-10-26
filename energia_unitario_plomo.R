
#CONTERO DE ATP POR SECUENCIA
#Es un conteo de genes totales 


setwd("~/Escritorio/DE/cd/DE_cd/ENERGIA_FINAL_CD/energia_unitaria/plomo")
conteo<-read.table("salida_conteo_todos_aa.txt", header = F)
#El vector conteo contiene el conteo de cada uno de los aminácidos para todos los
#genes, no separado por sobre y sub.
#Esto fue obtenido con el srcipt de perl "contador_aa.pl" de la misma carpeta


#View(conteo)

C=c(conteo[2:29128,2])
#View(C)

A=c(conteo[29130:58256,2])

R=c(conteo[58258:87384,2])

N=c(conteo[87386:116512,2])

D=c(conteo[116514:145640,2])

Q=c(conteo[145642:174768,2])

E=c(conteo[174770:203896,2])

G=c(conteo[203898:233024,2])

H=c(conteo[233026:262152,2])

I=c(conteo[262154:291280,2])

L=c(conteo[291282:320408,2])

K=c(conteo[320410:349536,2])

M=c(conteo[349538:378664,2])

f=c(conteo[378666:407792,2])

P=c(conteo[407794:436920,2])

S=c(conteo[436922:466048,2])

t=c(conteo[466050:495176,2])

W=c(conteo[495178:524304,2])

Y=c(conteo[524306:553432,2])

V=c(conteo[553434:582560,2])

##Calculo gasto por AA
#df=(cbind(C,A,R,N,D,Q,E,G,H,I,L,K,M,f,P,S,t,W,Y,V))
#Los números corresponden al número de moléculas necesarias para sintetizar el Aa
df2=(data.frame(C*24,A*12.5,R*18.5,N*4,D*1,Q*8.5,E*9.5,G*14.5,H*33,I*20,L*33,K*18.5,M*18.5,f*63,P*12.5,S*15,t*6,W*78.5,Y*56.5,V*25))



#ACA OBTENTO LA TABLA CON LOS VALORES TOTALES DE costo de todos los AA POR GEN, 
#SIN LA EXPRESION.
sss<-rowSums(df2)
View(sss)

conteoscd1<-read.csv("Pb.csv", header = T)
View(conteoscd1)
#Lo convierto a data frame ya que sss no tiene header, por lo que no hace el cbind
#ya que no estarían valanceados. Por eso el dataframe con el header =T para sss.
sss1=as.data.frame(sss, header = T)
p=cbind(conteoscd1[,1],sss1)
View(p)
write.table(p, file="unitario_pb_energia_de.txt", quote = T, sep="\t", row.names = F)

#Ahora busco los datos de expresión normalizada



library("DESeq2")
# Leemos conteos 
conteoscd<-read.table("Pb.csv", header = T, row.names = "GEN", sep = ",")


#View(conteoscd)
dim(conteoscd)# dimensiones de la tabla

# La serie "X1000_follow_up" son los No Meditadores, y 
# la serie "X2000_follow_up" son los Meditadores regulares.
colnames(conteoscd[1])  #nombre de columna 1
colnames(conteoscd[2]) #nombre de columna 56
colnames(conteoscd[7]) #nombre de columna 57
#col 2 a 56 son los no meditadores 1 mil ...
#57 adelante meditadores 2 mil..


#Creaci?n de objeto countdata, a partir del Subseting de los datos "ConteosMeditacion". 
# Se selecciona Todas las filas y desde la columna 2 hasta la 7:
countdata<-conteoscd[, c(1:6)]
#filtrado de conteos
#me quedo solo con los conteos, no con los nombres
dim(countdata) #dimensiones del objeto "countdata"


#coldata: tabla con informacion de las muestras
coldata = data.frame(
  #crear data frame con dos columnas
  row.names = colnames( countdata ),
  #nombres con lo anterior
  tipo = c( rep("Control", 3), rep("CD", 3)))
# a los datos de countdata, escribir No meditador 55 veces, etc


# Construimos el objeto DESeqDataSet a partir de la matriz de conteos y la informaci?n de las muestras
dds<- DESeqDataSetFromMatrix(countData = countdata,
                             #Entrega conteos, metadata, y permite comparar entre tipo Nomed y medit (que es desing=tipo)
                             colData = coldata,
                             design = ~ tipo)
t <- proc.time() # Inicia el cron?metro
dds <- DESeq(dds)
proc.time()-t    # Detiene el cron?metro. 
#En mi compu demora alrededor de 74.42 s aprox.



##########
#OBTENCION DE DATOS NORMALIZADOS - OPCIONAL
##########



# Obtener conteos normalizados (opcional, si se necesita el dato de normalizacion)
sf <- estimateSizeFactors(dds)
normctn<-counts(sf, normalized=TRUE)
#Si imprimo aca conserva los nombres.
write.table(normctn, file="tratamiento_normalizados_PB.txt", quote = T, sep="\t", row.names = T)


#####MULTIPLICACION
#LEO EL ARCHIVO ANTERIOR

norm_pb=read.table("tratamiento_normalizados_PB.txt", header = T)
View(norm_pb)

#multiplicamos

energia_de= norm_pb*sss
View(energia_de)

#Para tener solo los nombres de que genes están solo en control y solo en Pb debo:
#Calcular logaritmo
#Restarle -1
#Calcular logaritmo
#eliminar todas las filas que tengan NA y INF (0 no)

log_energia_de=log(energia_de)
uno_log_energia_de=log_energia_de-1
na_uno_log_energia_de=log(uno_log_energia_de)
View(na_uno_log_energia_de)
#Ahora los separo en dos grupos para que al eliminar -Inf y Na por separado me queden
#Control y tratamiento con sus correspondientes genes expr. que pueden no estar en
#El otro.
control=as.data.frame(na_uno_log_energia_de[,1:3])
View(control)
write.table(control, file="control_na_uno_log_energia_de.txt", quote = T, sep="\t", row.names = T)
#Ahora en linux elimino las filas que contengan NA y -Inf
#sed '/-Inf/d; /NA/d' control_na_uno_log_energia_de.txt > control_sin_na_uno_log_energia_de.txt

pb=as.data.frame(na_uno_log_energia_de[,4:6])
View(pb)
write.table(pb, file="pb_na_uno_log_energia_de.txt", quote = T, sep="\t", row.names = T)
#Ahora en linux elimino las filas que contengan NA y -Inf
#sed '/-Inf/d; /NA/d' pb_na_uno_log_energia_de.txt > pb_sin_na_uno_log_energia_de.txt



#Hago un grep entre solo los nombres de los archivos obtenidos y su valor normalizado.
#awk '{print $1}' control_sin_na_uno_log_energia_de.txt > nombres_solo_control
#awk '{print $1}' pb_sin_na_uno_log_energia_de.txt > nombres_solo_pb

#grep -wf nombres_solo_control tratamiento_normalizados_PB.txt > solo_control_norm.csv
#grep -wf nombres_solo_pb tratamiento_normalizados_PB.txt > solo_pb_norm.csv


#Ahora leo lo anterior y multiplico por el valor de energía de 1 gen

solo_control_norm=read.csv("solo_control_norm.csv", header = T,sep = "\t")
#View(solo_control_norm)
solo_pb_norm=read.table("solo_pb_norm.csv", header = T,sep = "\t")
#View(solo_pb_norm)


###LOS ARCHIVOS ANTERIORES TIENEN QUE MULTIPLICARSE POR EL VALOR UNITARIO DE ENERGIA

#hacemos un grep con los nombres de cada trat contra los valores de energia uniq

#grep -wf nombres_solo_control unitario_pb_energia_de.txt > control_de_en.csv
#grep -wf nombres_solo_pb unitario_pb_energia_de.txt > pb_de_en.csv

#Leo los archivos anteriores del grep

control_en_unit=read.csv("control_de_en.csv", header = F, sep = "\t")
View(control_en_unit)
pb_en_unit=read.csv("pb_de_en.csv", header = F, sep = "\t")
View(pb_en_unit)

control1=control_en_unit$V2*solo_control_norm[,1:3]
View(control1)
pb_1=pb_en_unit$V2*solo_pb_norm[,4:6]
View(pb_1)


vector_de_pb=c(pb_1$XI_PB_r1_07,pb_1$XI_PB_r2_08,pb_1$XI_PB_r3_09)
View(vector_de_pb)
vector_de_co=c(control1$XI_CO_r1_04,control1$XI_CO_r2_05,control1$XI_CO_r3_06)
View(vector_de_co)

vector_de11=c(vector_de_co,vector_de_pb)
View(vector_de11)

#Ahora crearemos los factores  vector_factores
#Creamos una columna para hacer la tabla de factores
v <- c("control")
w <- rep(v, times = 65124)
View(w)
ww=cbind(vector_de_co,w)
View(ww)

v1 <- c("Pb")
w1 <- rep(v1, times = 63426)
#View(w1)
ww1=cbind(vector_de_pb,w1)

View(ww1)

el_bind0_factores=c(w,w1)
el_bind0_en=c(vector_de_co,vector_de_pb)
el_bind1_en=log(el_bind0_en)
el_bind22=cbind(el_bind0_factores,el_bind1_en)


el_bind2=as.data.frame(el_bind22)


library(Rcmdr)

#'2.1' %in% factores4$Control
data= el_bind2


plotMeans(as.numeric(el_bind2$el_bind0_factores),el_bind2$el_bind0_en,error.bars="conf.int", level=0.95)

round_de=as.character(el_bind2$el_bind1_en)
View(round_de)
round_de1=(round((as.numeric(round_de)), digits = 3))
View(round_de1)
###Con esto elimino los puntos flotantes producidos por mucho decimal
###Y no se generan vectores de doble posicion que matan el eje Y

library(plotrix) #Para axis.break
plotMeans(round_de1,el_bind2$el_bind0_factores, ylim=c(16.8, 17.2),error.bars="conf.int", level=0.95, main = NULL, ylab="Numbers of ATP molecules", xlab= "Treatments")
axis.break(axis=2,breakpos=16.8,style="slash")


#Mas iguales, mas juntos
plotMeans(round_de1,el_bind2$el_bind0_factores, ylim=c(0, 30),error.bars="conf.int", level=0.95, main = NULL, ylab="Numbers of ATP molecules", xlab= "Treatments")
axis.break(5,10, axis=2,style="zigzag")


#### ACA VOY


#axis 2 es el eje y
#breakpos=Donde quiere el corte
#style =estilo, en este caso dos barras paralelas
#Ojo que al graficar ordena por roden alfabético, por eso puede ser que quede cadmio antes que control

#otro corte distinto

View(el_bind2$vector_factores)bla4=cbind(cad,w22)

# Distribuyen normales?
#Esto se hace NO con los valores log?
#Kolmogorov-Smirnov conocida como test Lilliefors que es para sobre 50 observaciones.
#Interpretación:
#Si el valor del p-value es menor a 0.05, la normal y mis datos son distintos. 

vector_de_co
View(vector_de_co)
vector_de_pb

#Como Kruskall-Wallis no acepta datos duplicados
#Elimino datos duplicados
vector_de_co_norep <- unique(vector_de_co)
View(vector_de_co_norep)

vector_de_pb_norep <- unique(vector_de_pb)


ks.test(x = vector_de_co_norep,pnorm, mean(vector_de_co_norep), sd(vector_de_co_norep))
#One-sample Kolmogorov-Smirnov test

#data:  vector_de_co_norep
#D = 0.45104, p-value < 2.2e-16
#alternative hypothesis: two-sided
# Resultado: Distintos, no normal
ks.test(x = vector_de_pb_norep,"pnorm", mean(vector_de_pb_norep), sd(vector_de_pb_norep))
#One-sample Kolmogorov-Smirnov test

#data:  vector_de_pb_norep
#D = 0.45988, p-value < 2.2e-16
#alternative hypothesis: two-sided
#Son iguales.


#Otra forma de calcularlo
#kruskal.test(vector_de11 ~ vector_factores22, data = el_bind222)





######
#Grafico de ven para genes compartidos y no

#hice en la terminal
#grep -wf solo_nombres_control_en_unit.csv solo_nombres_cadmium_en_unit.csv > compartidos_co_ca.csv
#20112 compartidos.
#Total control = 21102
#Total plomo = 20554

solo_control= 21102 - 20112
solo_plomo=20554 - 20112

#Solo control=990
#Solo plomo= 442

# Calculo para estimar la pbb de que los resultados anteriores se den por azar.

library(gmp)

####Correr por default
#####Crear la función que calcula p = Σ (m,i)(N-m,n-i)/(N,n)
enrich_pvalue <- function(N, A, B, k)
{
  m <- A + k
  n <- B + k
  i <- k:min(m,n)
  
  as.numeric( sum(chooseZ(m,i)*chooseZ(N-m,n-i))/chooseZ(N,n) )
}

#Ingresar los valores
#N=14800, A=Total gupo a (menos los comunes), B = Total gupo B (menoslos comunes)
#K=comunes

#From 15220 genes, set A is 1850+195 genes, 
#set B is 195+596 genes, overlap is 195 genes. Their p value is 2e-26. 
#enrich_pvalue(15220, 1850, 596, 195)

valorp=enrich_pvalue(15220, 1850, 596, 195)
#todos los genes involucrados  # genes solo co #genes solo trat #genes comunes

#Para los datos de Pb
valorp=enrich_pvalue(29127, 990, 442, 20112)
From 15220 genes, set A is 1850+195 genes, set B is 195+596 genes, overlap is 195 genes. Their p value is 2e-26. 
valorp1=enrich_pvalue(N=15220, A=1850, B=596, k=195)
valorp1=enrich_pvalue(N=29127, A=990, B=442, k=20112)
#Resultado = 0

summary(valorp1)

#RESULTADO: 1.91221e-18

#Numero esperado por azar

#Solo co(solo_cd/total)=
#15


#### INTERVALOS DE CONFIANZA

vector_de=c(con,cad)
View(vector_de)

n<-41656

#CALCULAR EL INTERVALO DE CONFIANZA PARA I
#Obtener la columna promedio para I

media <- mean(vector_de)
desv <- sd(vector_de)
nivelconfianza = 0.95



error.est <- desv/sqrt(n) # Calculamos el error estándar
margen.error <- 1.960 * error.est  # nivel de confianza de 95% 
#si fuera del 90% =1.644854


lim.inf <- media - margen.error # Límite inferior del intervalo
lim.inf
#17.11395

lim.sup <- media + margen.error # Límite superior del intervalo
lim.sup
#17.15393

sd(vector_de)
#2.081326


#### OTRA FORMA DE CALCULAR IC
# Pero tiene que ser por separado, ya que es por tratamiento.

wilcox.test(con,
            alternative="two.sided",
            correct=TRUE,
            conf.int=TRUE,
            conf.level=0.95)

#Wilcoxon signed rank test with continuity correction

#data:  con
#V = 222657753, p-value < 2.2e-16
#alternative hypothesis: true location is not equal to 0
#95 percent confidence interval:
#  17.13084 17.19370
#sample estimates:
# (pseudo)median 
#17.16228 

wilcox.test(log(vector_de_co),
            alternative="two.sided",
            correct=TRUE,
            conf.int=TRUE,
            conf.level=0.95)