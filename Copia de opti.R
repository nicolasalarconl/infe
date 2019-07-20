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
instancia<-readQAP("bur26a.dat")

intercambio<- function(entrada,cambio1,cambio2){
  aux =entrada[cambio1]
  entrada[cambio1]=entrada[cambio2]
  entrada[cambio2]=aux
  return(entrada)
}


actualizar<-function(temperatura){
 # temperatura <- 0.85*temperatura
  temperatura <- 0.95*temperatura
  return(temperatura)
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
simulatedVecindadBloque<- function(TemperaturaInicial,condicionTemperatura,condicionEquilibrio){
  #random 
  s <- sample(1:instancia$n,instancia$n,replace=F)
  evaluarS<-evaluarQAP(s,instancia$f,instancia$d)
  resultados <-evaluarS
  Temperatura<- TemperaturaInicial
  condicion <- condicionEquilibrio

  while(Temperatura > condicionTemperatura){
    while(condicion > 0){
      s2<-insercion(s,10,2)
      evaluarS2<-evaluarQAP(s2,instancia$f,instancia$d)
      condicion <- condicion -1  
      if (condicion ==1){
            s2 <- sample(1:instancia$n,instancia$n,replace=F)
            evaluarS2<-evaluarQAP(s,instancia$f,instancia$d)
      }
      delta <- evaluarS2 -evaluarS
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
    condicion<-condicionEquilibrio
    Temperatura <-actualizar(Temperatura) 

  }
  
  
  return (resultados)
}  
  




#prueba1<-c(simulatedVecindadBloque(100000000,10,100)) 
prueba1<-c(simulatedVecindadBloque(100000000000000,1,10)) 
plot(prueba1)




