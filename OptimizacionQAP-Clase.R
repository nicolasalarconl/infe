

library(ez) # ANOVA
library(ggpubr) # GRAFICOS
library(WRS2)

readQAP<-function(name){ 
a <- read.delim(name,header=FALSE, sep ="")
n<-as.integer(a[1,1])
fl<-a[2:(n+1),1:n]
dis<-a[(n+2):(n+n+1),1:n]
d <- list(n=n, f= fl, d = dis)
return(d)
}

evaluarQAP<-function(sol, f, d){
  
  acum<-0
  for(i in 1:n){
    for(j in 1:n){
      acum = acum + f[i,j]*d[sol[i],sol[j]]   
    }
  }
  return(acum)
}

intercambio<- function(entrada,cambio1,cambio2){
  aux =entrada[cambio1]
  entrada[cambio1]=entrada[cambio2]
  entrada[cambio2]=aux
  return(entrada)
}

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

insercion  <- function(entrada,Mover,Hacia) {
  aux <-entrada[Hacia]
  entrada[Hacia]<- entrada[Mover]
  if (Mover < Hacia){
    #realiza solo un cambio
    entrada[Mover]<-aux
  }
  
  if  (Mover >Hacia){
    Hacia<- Hacia +1 
    while (Hacia <= Mover){
      aux2 = entrada[Hacia]
      entrada[Hacia]=aux
      aux = aux2
      Hacia <-Hacia +1
    }
  }
  
  return (entrada)
}
#swap<-function(sol,i,j){
#  piv<-sol[i]
#  sol[i]<-sol[j]
#  sol[j]<-piv
#  return(sol)
#}


#leer instancia, crear y evaluar una solucion inicial
instancia<-readQAP("bur26a.dat")
#s<-c(1:instancia$n) #solucion inicial 
#fitness<-evaluarQAP(sol,instancia$f,instancia$d)
#fitness


actualizar<-function(temperatura){
  temperatura <- temperatura-1
  return(temperatura)
}

#generar soluci?n aleatoria y evaluar


simulatedVecindadBloque<- function(TemperaturaInicial,
                     condicionTemperatura,condicionEquilibrio){
  #random 
  s <- sample(1:instancia$n,instancia$n,replace=F)
  evaluarS<-evaluarQAP(s,instancia$f,instancia$d)
  
  
  resultados <-evaluarS

  Temperatura<- TemperaturaInicial # temperatura inicial
  condicion <- condicionEquilibrio
  
  while(Temperatura > condicionTemperatura){
    while(condicion > 0){
 
    
      s2<-insercion(s,10,2)
      evaluarS2<-evaluarQAP(s2,instancia$f,instancia$d)
      
      if (condicion ==1){
        s2 <- sample(1:instancia$n,instancia$n,replace=F)
        evaluarS2<-evaluarQAP(s,instancia$f,instancia$d)
      }
      
      delta <- evaluarS2 -evaluarS
      
      print(c(evaluarS2,evaluarS,evaluarS2-evaluarS,evaluarS-evaluarS2,delta))
 
      if (delta <= 0){
         s<-s2
         evaluarS<-evaluarS2
      }else{
         probabilidad <- exp(-delta/Temperatura)
         if(runif(1,0,1)<probabilidad){
           s<-s2
           evaluarS<-evaluarS2
         }
      }
      resultados <-c(resultados,evaluarS)
      condicion <- condicion -1 
    }
    condicion <-condicionEquilibrio
    Temperatura <-actualizar(Temperatura)
  }
  return (resultados)
}

prueba1<-c(simulatedVecindadBloque(10000,0,10)) 
plot(prueba1)

  
#########################################

simulatedInsert<- function(TemperaturaInicial,
                                   condicionTemperatura,condicionEquilibrio){
  s <- sample(1:instancia$n,instancia$n,replace=F)
  evaluarS<-evaluarQAP(s,instancia$f,instancia$d)
  resultados <-evaluarS
  Temperatura<- TemperaturaInicial # temperatura inicial
  condicion <- condicionEquilibrio
  while(Temperatura > condicionTemperatura){
    while(condicion > 0){
      # generar vecino con insercion
      # s2 <-insercion(s,10,2)
      if (condicion ==1){
        s <- sample(1:instancia$n,instancia$n,replace=F)
        evaluarS<-evaluarQAP(s,instancia$f,instancia$d)
      }
      s2<-insercion(s,10,2)
      evaluarS2<-evaluarQAP(s2,instancia$f,instancia$d)
      delta <- evaluarS2-evaluarS
      if (delta<= 0){
        s<-s2
      }else{
        probabilidad <- exp(1)^(-delta/Temperatura)
        aceptar<-sample(c(0,1), size=1,prob=c(1-probabilidad,probabilidad))
        if (aceptar ==1){
          s<-s2
        }
      }
      resultados <-c(resultados,evaluarS)
      condicion <- condicion -1 
    }
    condicion <-condicionEquilibrio
    Temperatura <-actualizar(Temperatura)
  }
  return (resultados)
}


#########

simulatedIntercambio<- function(TemperaturaInicial,
                                   condicionTemperatura,condicionEquilibrio){
  s <- sample(1:instancia$n,instancia$n,replace=F)
  evaluarS<-evaluarQAP(s,instancia$f,instancia$d)
  resultados <-evaluarS
  Temperatura<- TemperaturaInicial # temperatura inicial
  condicion <- condicionEquilibrio
  while(Temperatura > condicionTemperatura){
    while(condicion > 0){
      # generar vecino con insercion
      # s2 <-insercion(s,10,2)
      if (condicion ==1){
        s <- sample(1:instancia$n,instancia$n,replace=F)
        evaluarS<-evaluarQAP(s,instancia$f,instancia$d)
      }
      s2<-intercambio(s,10,1)
      evaluarS2<-evaluarQAP(s2,instancia$f,instancia$d)
      delta <- evaluarS2-evaluarS
      if (delta<= 0){
        s<-s2
      }else{
        probabilidad <- exp(1)^(-delta/Temperatura)
        aceptar<-sample(c(0,1), size=1,prob=c(1-probabilidad,probabilidad))
        if (aceptar ==1){
          s<-s2
        }
      }
      resultados <-c(resultados,evaluarS)
      condicion <- condicion -1 
    }
    condicion <-condicionEquilibrio
    Temperatura <-actualizar(Temperatura)
  }
  return (resultados)
}



prueba1min <-c(min(prueba1))
for (k in 1:14){
  nuevaprueba <- simulatedVecindadBloque(100,0,10)
  prueba1<-c(prueba1,nuevaprueba)
  prueba1min <-c(prueba1min,min(nuevaprueba))
}

plot(prueba1,main="14 ejecuciones de vecindad por bloque")
plot(prueba1min,main="14 minimos de 14 ejecuciones de vecindad por bloque")

prueba2<-c(simulatedInsert(100,0,10))
prueba2min<-c(min(prueba2))
for (k in 1:14){
  nuevaprueba2 <- simulatedInsert(100,0,10)
  prueba2<-c(prueba2,nuevaprueba2)
  prueba2min <-c(prueba2min,min(nuevaprueba2))
}

#plot(prueba2,main="14 ejecuciones de vecindad por Inserción")
#plot(prueba2min,main="14 minimos de 14 ejecuciones de vecindad por Inserción")


prueba3<-c(simulatedIntercambio(100,0,10))
prueba3min<-c(min(prueba3))
for (k in 1:14){
  nuevaprueba3 <- simulatedIntercambio(100,0,10)
  prueba3<-c(prueba3,nuevaprueba3)
  prueba3min <-c(prueba3min,min(nuevaprueba3))
}


id<- 1:15
datos.wide <- data.frame(id,prueba1min,prueba2min,prueba3min)

# Pero los procedimientos para hacer ANOVA, y muchas rutinas para
# graficar, requieren los datos en formalo largo (long).
dl <- gather(
  data = datos.wide,
  key = "Algoritmo",
  value = "Mejor",
  -id
)


p1 <- ggboxplot(
  dl,
  x="Algoritmo",y="Mejor",
  color="Algoritmo"
)
print(p1)


p2 <- ggboxplot(
  prueba2,
  add = "jitter",
)
print(p2)

apply(datos.wide, 2, shapiro.test)
rmanova(dl$Mejor, dl$Algoritmo,dl$id, tr = 0.05)








