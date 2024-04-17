#### SIMULACIONES ####
####    S-VAR     ####

#PAQUETES
require(abind)
require(sandwich)
require(AER)

set.seed(123)
#---------------------x-------------------------
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
coef.Mu     = array(NA, dim=c(K,N), dimnames=list(paste0("Var",1:K),paste0("Reg",1:N)))
coef.Beta   = array(NA, dim=c(K,K*P,N), dimnames=list(paste0("Var.",1:K), lab.vars.lag, paste0("Reg",1:N)))
coef.Theta  = array(NA, dim=c(K,K,N),dimnames=list(paste0("Var.ast.hat",1:K),paste0("Var.ast.hat",1:K), paste0("Reg", 1:N)))
coef.Lambda = array(NA, dim=c(K,K*P,N), dimnames=list(paste0("Var.",1:K), lab.var.hat.ast.lag, paste0("Reg",1:N)))
num.coef = 1 + K + (K*P) + (K*P)
#----- Estimacion del modelo SP-VAR obtenido de la segunda etapa -----#
for (k in 1:K) {
  for (n in 1:N) {
    modelo = lm(data[, paste0('Var',k), n] ~ data[, lab.vars.lag, n] +
                  data[, paste0("Var.ast.hat",1:K), n] +
                  data[,lab.var.hat.ast.lag , n])
    # Almacenar los coeficientes del modelo en el array
    coef.Mu[k,n]      = coef(modelo)[1]
    coef.Beta[k,,n]   = coef(modelo)[(1+1):(1+K*P)]
    coef.Theta[k,,n]  = coef(modelo)[(1+K*P +1):(1+K*P +K)]
    coef.Lambda[k,,n] = coef(modelo)[(1+K*P+K+1):(1+K*P+K+K*P)]
  }
}

if(0){
 print(coef.Mu)
 print(coef.Beta)
 print(coef.Theta)
 print(coef.Lambda)
}

# Construccion de las matrices tilde usadas para calcular la IRF 
coef.Theta.tilde = matrix(0, nrow = K*N, ncol = K*N)
for (n in 1:N)
  coef.Theta.tilde[((n-1)*K+1):(n*K), ((n-1)*K+1):(n*K)] = coef.Theta[,,n]



#FUNCION DE IMPULSO RESPUESTA
wkron              = kronecker(diag(K),W)
theta.tilde.array  = kronecker(diag(N),coef.Theta)
beta.tilde.array   = kronecker(diag(N),coef.Beta)
lambda.tilde.array = kronecker(diag(N),coef.Lambda)

resul  = array(NA, dim = J+1)
for (j in 1:J){
  resul.j <- ((solve(diag(N*K)-(theta*wkron)))*(beta+(lambda*wkron)))^{j}
  resul[j] <- resul.j
}



