#### SIMULACIONES ####
####    S-VAR     ####

#PAQUETES
#require()

set.seed(123)

#FUNCIONES
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
}
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
}



# PARÁMETROS
N  = 4  # Regiones
K  = 2  # Variables
T = 10 #Tiempo
p  = 1  #Rezagos


# Generar las series de tiempo
Y = generar_datos(T, K, N)

# Generar la matriz de contigüidad
W   = generar_matriz_contiguidad(N)

#Generar los vectores Nx1
Y.ast = array(NA, dim = c(T, K, N),dimnames=list(paste0("t", 1:T),
                                                 paste0("Var", 1:K),
                                                 paste0("Reg", 1:N)))
for (t in 1:T){
  for (k in 1:K)
    Y.ast[t,k,] = W%*%series_tiempo[t,k,]
}

###Loop t-1

