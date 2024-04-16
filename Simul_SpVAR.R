#### SIMULACIONES ####
####    S-VAR     ####

#PAQUETES
require(abind)
require(sandwich)
require(AER)

set.seed(123)
#---------------------x-------------------------#
#FUNCIONES
#---------------------x-------------------------#
generar_datos <- function(t, K, N) {
  # Crear un array vacío para almacenar los datos
  datos <- array(NA, dim = c(t, K, N))
  
  # Etiquetas para las filas (dimensión t)
  row_labels <- paste0("t", 1:t)
  
  # Etiquetas para las columnas (dimensión n)
  col_labels <- paste0("Var", 1:K)
  
  # Etiquetas para las capas (dimensión k)
  layer_labels <- paste0("Reg", 1:N)
  
  # Asignar nombres a las dimensiones
  dimnames(datos) <- list(row_labels, col_labels, layer_labels)
  
  # Simular los datos y almacenarlos en el array
  for (i in 1:t) {
    for (j in 1:K) {
      for (m in 1:N) {
        datos[i, j, m] <- runif(1) # Simulación de datos aleatorios, ajusta según tu necesidad
      }
    }
  }
  
  return(datos)
} #Función generadora de datos
generar_matriz_contiguidad <- function(num_regiones) {
  # Crear una matriz de contigüidad inicializada con ceros
  matriz_contiguidad <- matrix(0, nrow = num_regiones, ncol = num_regiones)
  
  # Generar distancias aleatorias para las regiones
  distancias <- runif(n = num_regiones^2, min = 0, max = 2)
  
  # Asegurar que la distancia entre una región y la misma región sea 0
  #distancias[seq(1, num_regiones^2, by = num_regiones + 1)] <- 0
  
  # Reflejar las distancias para que la matriz sea simétrica
  distancias <- matrix(distancias, nrow = num_regiones)
  distancias <- (distancias + t(distancias)) / 2
  
  # Calcular la matriz de contigüidad inversa
  matriz_contiguidad <- 1 / distancias
  
  # Asegurar que la diagonal de la matriz de contigüidad sea 0
  diag(matriz_contiguidad) <- 0
  
  return(matriz_contiguidad)
} #Función generadora de la matriz de contiguidad


# PARÁMETROS
N  = 4  # Regiones
K  = 2  # Variables
T  = 20 #Tiempo
P  = 2  #Rezagos
J  = 10 #Horizontes adelante


# Generar las series de tiempo
Y = generar_datos(T, K, N)

# Generar la matriz de contigüidad
W = generar_matriz_contiguidad(N)

#Generar los vectores Nx1
Y.ast = array(NA, dim = c(T, K, N),dimnames=list(paste0("t", 1:T),
                                                 paste0("Var.ast", 1:K),
                                                 paste0("Reg", 1:N)))
for (t in 1:T){
  for (k in 1:K)
    Y.ast[t,k,] = W%*%Y[t,k,]
}
data <- abind(Y, Y.ast, along = 2)

# Bucle para crear matrices retrasadas y agregarlas a <data>
lab.vars.lag = NULL; lab.vars.ast.lag = NULL; lab.var.hat.ast.lag = NULL
for(p in 1:P){
  lab.vars.lag        = c(lab.vars.lag, paste0("Var", 1:K,".lag",p))
  lab.vars.ast.lag    = c(lab.vars.ast.lag, paste0("Var.ast", 1:K,".lag",p))
  lab.var.hat.ast.lag = c(lab.var.hat.ast.lag, paste0("Var.ast.hat", 1:K,".lag",p))
}

Y.Lag <- array(NA, dim = c(T, K*P, N),dimnames=list(paste0("t", 1:T),
                                                    lab.vars.lag,
                                                    paste0("Reg", 1:N)))
Y.ast.Lag <- array(NA, dim = c(T, K*P, N),dimnames=list(paste0("t", 1:T),
                                                        lab.vars.ast.lag,
                                                        paste0("Reg", 1:N)))

# Crear un nuevo array para almacenar las variables retrasadas
for (p in 1:P) {  
  for (t in (p+1):T){ 
    Y.Lag[t,(K*(p-1)+1):(K*p) , ]     = Y[t-p, , ]
    Y.ast.Lag[t,(K*(p-1)+1):(K*p) , ] = Y.ast[t-p, , ]
  }    
}
data <- abind(data, Y.Lag, Y.ast.Lag, along = 2)

# Ecuación 30
# Y.ast.hat.temp    = list()
Y.ast.hat         = array(NA, dim = c(T, K, N),dimnames=list(paste0("t", 1:T),
                                                             paste0("Var.ast.hat", 1:K),
                                                             paste0("Reg", 1:N)))

for (k in 1:K) {
  for (n in 1:N) {
    modelo = lm(data[, paste0('Var',k), n] ~ data[, lab.vars.lag, n] + data[,lab.vars.ast.lag , n])
    # Almacenar los coeficientes del modelo en la lista
    Y.ast.hat[c(-1:-P),k,n] = predict(modelo)
  }
}

for (t in 1:T){
  for (k in 1:K)
    Y.ast.hat[t,k,] = W%*%Y.ast.hat[t,k,]
}

Y.ast.hat.lag <- array(NA, dim = c(T, K*P, N),dimnames=list(paste0("t", 1:T),
                                                            lab.var.hat.ast.lag,
                                                            paste0("Reg", 1:N)))
for (p in 1:P) {  
  for (t in (p+1):T)
    Y.ast.hat.lag[t,(K*(p-1)+1):(K*p) , ] = Y.ast.hat[t-p, , ]
}
data <- abind(data, Y.ast.hat, Y.ast.hat.lag, along = 2)

#Ecuacion 16a
#coef.modelo = list()
num.coef = 1 + K + (K*P) + (K*P)
coef.modelo = array(NA, dim = c(K*N,num.coef))
list.cov = list()
#LM
i=0
for (k in 1:K) {
  for (n in 1:N) {
    modelo = lm(data[, paste0('Var',k), n] ~ data[, lab.vars.lag, n] +
                  data[, paste0("Var.ast.hat",1:K), n] +
                  data[,lab.var.hat.ast.lag , n])
    
    # Almacenar los coeficientes del modelo en la lista
    # coef.modelo[[paste("k", k, "n", n, sep = "_")]] <- coef(modelo)
    i=i+1
    coef.modelo[i,] <- coef(modelo)
    # list.cov[[paste("k", k, "n", n, sep = "_")]] <- sqrt(diag(vcovHC(modelo)))
  }
}


#IVREG
# for (k in 1:K) {
#   for (n in 1:N) {
#     # Modelo con variables instrumentales
#     data_df <- as.data.frame(data)
#     modelo <- ivreg(data[, paste0('Var', k), n] ~ data[, lab.vars.lag, n]  
#                     | data[, paste0("Var.ast.hat",1:K), n] +
#                       data[, lab.var.hat.ast.lag , n], data = data_df)
#     
#     # Almacenar los coeficientes del modelo en la lista
#     coef.modelo[k*n,] <- coef(modelo)
#     
#     # Calcular intervalos de confianza para los coeficientes
#     intervalos_confianza <- confint(modelo)
#     print(intervalos_confianza)  # Esto imprimirá los intervalos de confianza en la consola
#   }
# }
#FUNCION DE IMPULSO RESPUESTA
wkron = kronecker(diag(K),W)
coef(modelo)
theta  = coef.modelo[,(1+(P*K)+1):(1+(P*K)+K)]
beta   = coef.modelo[,(1+1):(1+(P*K))]
lambda = coef.modelo[,(1+(P*K)+K+1):(1+(P*K)+K+(P*K))]
theta  = kronecker(diag(N),theta)
beta   = kronecker(diag(N),beta)
lambda = kronecker(diag(N),lambda)

theta       = coef.modelo[,(1+(P*K)+1)]
thetak      = matrix(0, nrow = length(theta), ncol = length(theta))
diag(thetak) = theta


# for (i in (1+(P*K)+1):(1+(P*K)+K)){
#   theta       = coef.modelo[,i]
#   paste0("thetak",i) = matrix(0, nrow = length(theta), ncol = length(theta))
#   diag(paste0("thetak",i)) = theta
# }

NK = N*K
c  = 0
matrices_theta <- array(0, dim = c(NK,NK, K))
for (i in (1+(P*K)+1):(1+(P*K)+K)) {
  theta <- coef.modelo[,i]
  thetak <- matrix(0, nrow = NK, ncol = NK)
  diag(thetak) <- theta
  c = c+1
  matrices_theta[,,c] <- thetak
}
NK = N*K
c  = 0
matrices_theta <- array(0, dim = c(NK,NK, K))
for (i in (1+(P*K)+1):(1+(P*K)+K)) {
  theta <- coef.modelo[,i]
  thetak <- matrix(0, nrow = NK, ncol = NK)
  diag(thetak) <- theta
  c = c+1
  matrices_theta[,,c] <- thetak
}

resul  = array(NA, dim = J+1)
for (j in 1:J){
  resul.j <- ((solve(diag(N*K)-(theta*wkron)))*(beta+(lambda*wkron)))^{j}
  resul[j] <- resul.j
}



