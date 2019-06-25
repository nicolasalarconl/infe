install.packages("futile.matrix")
library("futile.matrix")

setwd("D:\\OneDrive\\Magister\\Optimizacion\\Metaheuristica\\T1")
 

leeArchivo <- function(filename){
  READONLY = "r"
  
  datafile = file(filename, READONLY)
  
  # lee la primera línea (número de filas/columnas)
  nn <- as.integer(readLines(datafile, n=1))
  n <- nn
  
  m <- read.table(datafile, header = FALSE)
  
  # lee la primera matriz
  m1 <- m[1:n,]
  #m1 <- t(m1)
  
  mm <- 2*n
  n <- n+1
  
  # lee la segunda matriz
  m2 <- m[n:mm,]
  
  
  return (list(N=nn, MF=m1, MK=m2))
}

evaluarQAP<-function(sol, f, d, n){
  
  acum<-0
  for(i in 1:n){
    for(j in 1:n){
      acum = acum + f[i,j]*d[sol[i],sol[j]]   
    }
  }
  return(acum)
}

# función: construyeVecindad1 
# Mueve n posiciones (pasos) los valores del arreglo 
# ej: 
#   array1: (1, 2, 3, 4, 5, 6, 7, 8, 9, 0)
#   con n = 6 quedaría
#   arrayvecino: (5, 6, 7, 8, 9, 0, 1, 2, 3, 4)
# Parámetros: 
#   arreglo: arreglo solución inicial
#   paso: cantidad de indices a mover
construyeVecindad1 <- function(arreglo, paso){
  largo <- length(arreglo)
  arrayvecino <- c(1:largo)
  
  for(i in 1:largo){
    a <- i + paso
    
    if(a>largo)
      a <- a - largo
    
    arrayvecino[a] <- arreglo[i]
  }
  
  return(arrayvecino) 
  
}

construyeVecindad2 <- function(arreglo){
  largo <- as.integer(length(arreglo))
  arrayvecino <- arreglo
 
   for(i in 1:largo){
    random_index <- as.integer(runif(1,min=1,max=largo))
    temporal <- arrayvecino[i]
    arrayvecino[i] <- arrayvecino[random_index]
    arrayvecino[random_index] <- temporal
  }
  return(arrayvecino)
}

simulatedAnnealing<-function(temperatura, mx){
  #sol_actual <- matrix()
  fitness <- vector()
  #s<-c(7,5,12,2,1,3,9,11,10,6,8,4) #Solución Óptima conocida para chr12a.dat
  s<-c(1:mx$N) 
  temperatura<-temperatura
  t_min<-temperatura*0.01
  
  while(t_min<temperatura){
    iteracion <- as.integer(0.25*mx$N) #25% de la vecindad   ???
    while(iteracion > 0){
        
        #s1 <- construyeVecindad1(s, as.integer(runif(1,min=1,max=mx$N-1)))
        s1 <- construyeVecindad2(s) 
        f_s1 <- evaluarQAP(s1, mx$MF, mx$MK, mx$N)
        f_s  <- evaluarQAP(s,  mx$MF, mx$MK, mx$N)
        dE <- f_s1 - f_s # si s1 es mayor que s, dE es positivo y significa que s1 es peor solucion
                         # que s 
        
        if(dE<=0){ 
          s <- s1 #Acepta la solución vecina
          f_s <- f_s1
        }else{
          prob <- exp(-dE/temperatura)
          random <- runif(1)
          if(prob>random){
            s <- s1
            f_s <- f_s1
          } 
        } 
        iteracion <- iteracion - 1
    }
    
    fitness <- c(fitness, f_s)
    #temperatura <- temperatura - 0.05*temperatura #enfriamiento lineal
    temperatura <- 0.96*temperatura #enfriamiento geométrico
  }
  
  
  return(list(FIT=fitness))
  #return(sol_actual[which.min(fitness)]) 
  
}
 

archivo = "chr12a.dat" #"bur26a.dat" # 
mx = leeArchivo(filename = archivo)

best <- simulatedAnnealing(300, mx)

print(best$FIT)
plot(best$FIT)
