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
  mutaciones <- matrix(ncol = n+1, nrow = 0)
  
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
  cruzamientos <- matrix(ncol = n+1, nrow = 0)
  n_padres = nrow(padres)
  
  for(i in 1:(repeticiones*n_padres)){
    cruzamientos <- rbind(cruzamientos, t(apply(padres, 1, function(padre){
      index_pareja <- ceiling((abs(rnorm(1))*n_padres)%%n_padres)
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

genetic <- function(n_poblacion, n_generaciones, n_mutaciones, n_cruzamientos,p,n){
  adaptativo <- FALSE
  if(n_generaciones == 0){ # adaptativo
    adaptativo <- TRUE
    n_generaciones <- 9999999
  }
  
  start.time <- Sys.time()
  
  poblacion_inicial <- permutations(1:n,nsample = n_poblacion)
  poblacion <- cbind(poblacion_inicial,replicate(n_poblacion,0)) # agrega una columna adicional de ceros
  
  mejores <- c() # mejores segun iteracion
  mejor <- c()
  sin_cambios <- 0
  history <- data.frame()
  
  for(i in 1:n_generaciones){
    poblacion <- evaluar_poblacion(poblacion, n, p)
    poblacion <- ordenar_por_valor(poblacion)
    mejores <- c(mejores, poblacion[1,n+1])
    mejor <- c(mejor, poblacion[1,1])
    history <- rbind(history, snapshot(poblacion,i))
    
    if(adaptativo && i>1 && mejor[i] == mejor[i-1]){
      sin_cambios <- sin_cambios + 1
      if(sin_cambios == 10){ # parada adaptativa
        break
      }
    }else{
      sin_cambios <- 0
    }
    
    # Seleccion
    padres <- poblacion[1:(n_poblacion/10),] # 10% mejor
    random <- poblacion[sample(nrow(poblacion),size=n_poblacion/10,replace=FALSE),]
    padres <- rbind(padres,random)
    
    # Mutacion
    mutaciones <- mutar(padres,n_mutaciones,n)
    
    # Cruzamiento
    cruzamiento <- cruzar(padres,n_cruzamientos,n)
    
    # Siguiente poblacion
    poblacion <- rbind(padres,mutaciones)
    poblacion <- rbind(poblacion,cruzamiento)
    
    if(n_poblacion > nrow(poblacion)){
      filler <- poblacion[sample(nrow(poblacion),size=(n_poblacion-nrow(poblacion)),replace=FALSE),]
      poblacion <- rbind(poblacion,filler)
    }
  }
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  
  
  #title <- sprintf("AG P=%i G=%i C=%i M=%i", n_poblacion, generaciones, cruzamientos_por_padre*20, mutaciones_por_padre*20)
  ## 1. Open jpeg file
  #filename <-sprintf("poblacion %s.png",title)
  #png(filename=filename, width = 800, height = 600)
  ## 2. Create the plot
  #plot(history$generacion,history$valores, main=title, xlab = "Generación", ylab = "Valor F.O.")
  ## 3. Close the file
  #dev.off()
  
  ## 1. Open jpeg file
  #filename <-sprintf("mejor %s.png",title)
  #png(filename=filename, width = 800, height = 600)
  ## 2. Create the plot
  #plot(1:generaciones,mejores, main=title, xlab = "Generación", ylab = "Valor F.O.")
  ## 3. Close the file
  #dev.off()
  
  mejor <- mejores[length(mejores)]
  optimo <- p$opt
  
  if(adaptativo){
    n_generaciones <- -1
  }
  
  return(data.frame(
    p=c(n_poblacion),
    g=c(n_generaciones),
    m=c(n_mutaciones*20),
    c=c(n_cruzamientos*20),
    mejor=c(mejor),
    error=c((mejor-optimo)*100/optimo),
    tiempo=c(as.integer(time.taken, units="secs"))
  ))
}

#### PARAMETROS ####

n <- 26
n_poblacion <- 50
generaciones <- 50
mutaciones_por_padre <- 2
cruzamientos_por_padre <- 2
p <- read_qaplib(system.file("qaplib","bur26a.dat",package = "qap"))

#### PARAMETROS ####

N <- 10 # iteraciones por configuracion
result = data.frame(
  p=c(),
  g=c(),
  m=c(),
  c=c(),
  mejor=c(),
  error=c(),
  tiempo=c()
)


# generaciones
for(generaciones in c(0, 25, 50, 100)){
  for(i in 1:N) {
    result <- rbind(result,genetic(n_poblacion,generaciones,mutaciones_por_padre, cruzamientos_por_padre,p,n))
  }
}


boxplot(error~g,data=result, main="Error Porcentual con p=50, g variable", ylab="% Error")
boxplot(tiempo~g,data=result, main="Tiempo para p=50, g variable", ylab="Tiempo [s]")


