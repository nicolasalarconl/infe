library(qap)
library(arrangements)

evaluarQAP<-function(sol, f, d){
  acum<-0
  n = length(sol)
  for(i in 1:n){
    for(j in 1:n){
        acum = acum + f[i,j]*d[sol[i],sol[j]]   
    }
  }
  return(acum)
}

evaluar_poblacion <- function(poblacion, n, p){
  return (
    t(apply(poblacion, 1, function(x){
      x[n+1]<-evaluarQAP(x[1:n],p$A,p$B)
      return(x)
    }))
  )
}

# funcion para hacer ranking
ordenar_por_valor <- function(poblacion){
  poblacion[order(poblacion[,ncol(poblacion)]),]
}

# swapea dos posiciones
mutar <- function(padres, repeticiones, n){
  mutaciones <- matrix(ncol = 27, nrow = 0)
  
  for(i in 1:repeticiones){
    mutaciones <- rbind(mutaciones, t(apply(padres, 1, function(padre){
      swap <- sample(1:n,size=2, replace = TRUE)
      swapped <- padre
      swapped[swap] <- swapped[rev(swap)]
      return(swapped)
    })))
  }
  return(mutaciones)
}

# toma dos padres, del primero saca los primeros p elementos, 
# elimina las caracterisicas ya seleccionadas del segundo padre, 
# y concatena los restantes del segundo padre respetando el orden que tenian
cruzar <- function(padres, repeticiones, n){
  cruzamientos <- matrix(ncol = 27, nrow = 0)
  n_padres = nrow(padres)
  
  for(i in 1:repeticiones){
    cruzamientos <- rbind(cruzamientos, t(apply(padres, 1, function(padre){
      index_pareja <- sample(1:n_padres,size=1)
      pareja <- padres[index_pareja,]
      pivote <- sample(1:n,size=1)
      
      hijo <- padre[0:pivote]
      hijo <- c(hijo, pareja[! pareja %in% hijo])

      return(hijo)
    })))
  }
  return(cruzamientos)
}

snapshot <- function(poblacion, generation){
  generacion <- replicate(nrow(poblacion),generation)
  valores <- poblacion[,ncol(poblacion)]
  return(data.frame(generacion,valores))
}


#### PARAMETROS ####

n <- 26
n_poblacion <- 100
generaciones <- 100
mutaciones_por_padre <- 2
cruzamientos_por_padre <- 2

#### PARAMETROS ####


p <- read_qaplib(system.file("qaplib","bur26a.dat",package = "qap"))

poblacion_inicial <- permutations(1:n,nsample = n_poblacion)
poblacion <- cbind(poblacion_inicial,replicate(n_poblacion,0)) # agrega una columna adicional de ceros

mejores <- c() # mejores segun iteracion
history <- data.frame()
for(i in 1:generaciones){
  # Evaluar y ordenar
  poblacion <- evaluar_poblacion(poblacion, n, p)
  poblacion <- ordenar_por_valor(poblacion)
  mejores <- c(mejores, poblacion[1,n+1])
  history <- rbind(history, snapshot(poblacion,i))
  
  # Seleccion
  padres <- poblacion[1:(n_poblacion/5),] # 20% mejor
  
  # Mutacion
  mutaciones <- mutar(padres,mutaciones_por_padre,n)
  
  # Cruzamiento
  cruzamiento <- cruzar(padres,cruzamientos_por_padre,n)
  
  # Siguiente poblacion
  poblacion <- rbind(padres,mutaciones)
  poblacion <- rbind(poblacion,cruzamiento)
}


# plot(history$generacion,history$valores)
plot(1:generaciones,mejores)
mejor <- mejores[length(mejores)]
optimo <- p$opt
print((mejor-optimo)*100/optimo)

