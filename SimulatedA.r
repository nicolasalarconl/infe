    #install.packages("futile.matrix")
    library("futile.matrix")
    library(qap)
    library(arrangements)
    library(ez)
    library(ggpubr)
    library(tidyr)
    library(WRS2)
    
    leeArchivo <- function(filename){
      READONLY = "r"
      
      datafile = file(filename, READONLY)
      
      # lee la primera l?nea (n?mero de filas/columnas)
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
    
    evaluarQAP<-function(sol, f, d,n){
    
      acum<-0
      for(i in 1:n){
        for(j in 1:n){
          acum = acum + f[i,j]*d[sol[i],sol[j]]   
        }
      }
      return(acum)
    }
    
    # funci?n: construyeVecindad1 
    # Mueve n posiciones (pasos) los valores del arreglo 
    # ej: 
    #   array1: (1, 2, 3, 4, 5, 6, 7, 8, 9, 0)
    #   con n = 6 quedar?a
    #   arrayvecino: (5, 6, 7, 8, 9, 0, 1, 2, 3, 4)
    # Par?metros: 
    #   arreglo: arreglo soluci?n inicial
    #   paso: cantidad de indices a mover
    construyeVecindad1 <- function(arreglo){
      largo <- as.integer(length(arreglo))
      paso <- as.integer(runif(1,min=1,max=largo))
    
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
    
    
    insercion  <- function(entrada){
      largo <- as.integer(length(entrada))
      Mover <- as.integer(runif(1,min=1,max=largo))
      Hacia <- as.integer(runif(1,min=1,max=largo))
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
    
    simulatedAnnealing1<-function(temperaturaInicial, mx){
    
      #sol_actual <- matrix()
      fitness <- vector()
      #s<-c(7,5,12,2,1,3,9,11,10,6,8,4) #Soluci?n ?ptima conocida para chr12a.dat
      s <- sample(1:mx$N,replace=FALSE)
      f_s  <- evaluarQAP(s,  mx$MF, mx$MK, mx$N)
    
      temperatura<-temperaturaInicial
      t_min<-temperaturaInicial*0.00001
      
      while(t_min<temperatura){
        iteracion <- 6 #as.integer(0.25*mx$N) #25% de la vecindad   ???
        while(iteracion > 0){
            if (iteracion ==1 ){
              s1 <- sample(1:mx$N,replace=FALSE)
              f_s1  <- evaluarQAP(s1, mx$MF, mx$MK, mx$N)
            }
            else{
              
              s1 <- construyeVecindad1(s) 
              f_s1 <- evaluarQAP(s1, mx$MF, mx$MK, mx$N)  
            }
            
            dE <- f_s1 - f_s 
            if(dE<=0){ 
              s <- s1 #Acepta la soluci?n vecina
              f_s <- f_s1
            }else{
              prob <- exp(-(dE/temperatura))
              random <- runif(1, 0, 1) #runif(1)
              if(prob>random){
                s <- s1
                f_s <- f_s1
              } 
            } 
            iteracion <- iteracion - 1
            fitness <- c(fitness, f_s)
        }
        
        
        #temperatura <- temperatura - 0.05*temperatura #enfriamiento lineal
        temperatura <- 0.90*temperatura #enfriamiento geom?trico
      }

      return(fitness) 
      
    }
    
    simulatedAnnealing2<-function(temperaturaInicial, mx){
      
      #sol_actual <- matrix()
      fitness <- vector()
      #s<-c(7,5,12,2,1,3,9,11,10,6,8,4) #Soluci?n ?ptima conocida para chr12a.dat
      s <- sample(1:mx$N,replace=FALSE)
      f_s  <- evaluarQAP(s,  mx$MF, mx$MK, mx$N)
      
      temperatura<-temperaturaInicial
      t_min<-temperaturaInicial*0.00001
      
      while(t_min<temperatura){
        iteracion <- 6 #as.integer(0.25*mx$N) #25% de la vecindad   ???
        while(iteracion > 0){
          if (iteracion ==1 ){
            s1 <- sample(1:mx$N,replace=FALSE)
            f_s1  <- evaluarQAP(s1, mx$MF, mx$MK, mx$N)
          }
          else{
            
            s1 <- construyeVecindad2(s) 
            f_s1 <- evaluarQAP(s1, mx$MF, mx$MK, mx$N)  
          }
          
          dE <- f_s1 - f_s 
          if(dE<=0){ 
            s <- s1 #Acepta la soluci?n vecina
            f_s <- f_s1
          }else{
            prob <- exp(-(dE/temperatura))
            random <- runif(1, 0, 1) #runif(1)
            if(prob>random){
              s <- s1
              f_s <- f_s1
            } 
          } 
          iteracion <- iteracion - 1
          fitness <- c(fitness, f_s)
        }
        
        
        #temperatura <- temperatura - 0.05*temperatura #enfriamiento lineal
        temperatura <- 0.90*temperatura #enfriamiento geom?trico
      }
      resp <- min(fitness)
      return(resp)
      #return(sol_actual[which.min(fitness)]) 
      
    }
    
    simulatedAnnealing3<-function(temperaturaInicial, mx){
      
      #sol_actual <- matrix()
      fitness <- vector()
      #s<-c(7,5,12,2,1,3,9,11,10,6,8,4) #Soluci?n ?ptima conocida para chr12a.dat
      s <- sample(1:mx$N,replace=FALSE)
      f_s  <- evaluarQAP(s,  mx$MF, mx$MK, mx$N)
      
      temperatura<-temperaturaInicial
      t_min<-temperaturaInicial*0.00001
      
      while(t_min<temperatura){
        iteracion <- 6 #as.integer(0.25*mx$N) #25% de la vecindad   ???
        while(iteracion > 0){
          if (iteracion ==1 ){
            s1 <- sample(1:mx$N,replace=FALSE)
            f_s1  <- evaluarQAP(s1, mx$MF, mx$MK, mx$N)
          }
          else{
            
            s1 <- insercion(s) 
            f_s1 <- evaluarQAP(s1, mx$MF, mx$MK, mx$N)  
          }
          
          dE <- f_s1 - f_s 
          if(dE<=0){ 
            s <- s1 #Acepta la soluci?n vecina
            f_s <- f_s1
          }else{
            prob <- exp(-(dE/temperatura))
            random <- runif(1, 0, 1) #runif(1)
            if(prob>random){
              s <- s1
              f_s <- f_s1
            } 
          } 
          iteracion <- iteracion - 1
          fitness <- c(fitness, f_s)
        }
        
        
        #temperatura <- temperatura - 0.05*temperatura #enfriamiento lineal
        temperatura <- 0.90*temperatura #enfriamiento geom?trico
        
      }s
      resp <- min(fitness)
      return(resp)
      #return(sol_actual[which.min(fitness)]) 
      
    }
    
    archivo = "bur26a.dat"
    mx = leeArchivo(filename = archivo)
    start.time <- Sys.time()

  
    best1 <- simulatedAnnealing1(3000000, mx)
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    plot(best1)
    s1<-c(best1)
    for (i in 1:15){
        best1 <- simulatedAnnealing1(3000000, mx)
      s1<-c(s1,best1)
    }
    plot(s1)
    start.time2 <- Sys.time()
    best2 <- simulatedAnnealing2(3000000, mx)
    end.time2 <- Sys.time()
    time.taken2 <- end.time2 - start.time2
    s2<-c(best2)
    for (i in 1:15){
      best2 <- simulatedAnnealing2(3000000, mx)
      s2<-c(s2,best2)
    }
    plot(s2)
    start.time3 <- Sys.time()
    best3 <- simulatedAnnealing3(3000000, mx)
    end.time3 <- Sys.time()
    time.taken3 <- end.time3 - start.time3
    s3<-c(best2)
    for (i in 1:15){
      best3 <- simulatedAnnealing3(3000000, mx)
      s3<-c(s3,best3)
    }
    
    plot(s3)

   
id<- 1:16
    datos.wide <- data.frame(id,s1,s2,s3)
    dl <- gather(
      data = datos.wide,
      key = "Algoritmo",
      value = "Fitness",
      -id
    )
    p1 <- ggboxplot(
      dl,
      x="Algoritmo",y="Fitness",
      color="Algoritmo"
    )
    
    plot(p1)
    
    test <- rmanova(dl$Fitness, dl$Algoritmo,dl$id, tr = 0.01)
    ph2 <- pairwise.wilcox.test(dl[["Fitness"]], dl[["Algoritmo"]],
                                paired = TRUE,
                                p.adjust.method = "BH")
    
    cat("\n\n")
    cat("Comparaciones entre pares de algoritmos (post-hoc)\n")
    cat("--------------------------------------------------\n")
    print(ph2)
    
    
  bestSoluction <-   5426670
  err1<- (100*(s1[1]-bestSoluction))/bestSoluction
  error1 <- (err)
  for (i in 2:16){
    err1<- (100*(s1[i]-bestSoluction))/bestSoluction
    error1 <- c(err1,error1)
  }
  errorPromedio1<-mean(error1)
  
  err2<- (100*(s2[1]-bestSoluction))/bestSoluction
  error2 <- (err2)
  for (i in 2:16){
    err2<- (100*(s2[i]-bestSoluction))/bestSoluction
    error2 <- c(err2,error2)
  }
  errorPromedio2<-mean(error2)
  err3<- (100*(s3[1]-bestSoluction))/bestSoluction
  error3 <- (err3)
  for (i in 2:16){
    err3<- (100*(s3[i]-bestSoluction))/bestSoluction
    error3 <- c(err3,error3)
  }
  errorPromedio3<-mean(error3)
  errorPromedio1
  errorPromedio2
  errorPromedio3
  
  
  datos2.wide <- data.frame(id,error1,error2,error3)
  dl2 <- gather(
    data = datos2.wide,
    key = "Algoritmo",
    value = "Error",
    -id
  )
  p2 <- ggboxplot(
    dl2,
    x="Algoritmo",y="Error",
    color="Algoritmo"
  )
  
  plot(p2)
  
  