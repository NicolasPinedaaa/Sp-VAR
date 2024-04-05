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
T = 20 #Tiempo
P  = 1  #Rezagos


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
    Y.ast[t,k,] = W%*%Y[t,k,]
}

###Loop t-1
for (p in 1:P) {
  # Crear un nuevo array para almacenar las variables retrasadas
  Y_Lag <- array(NA, dim = dim(Y),dimnames=list(paste0("t", 1:T),
                                                paste0("Var", 1:K),
                                                paste0("Reg", 1:N)))
  
  # Aplicar el retraso en cada dimensión
  for (t in 1:(T-p)) {
    Y_Lag[t + p, , ] <- Y[t, , ]
  }
  assign(paste0("Ylag_", p), Y_Lag)
}
for (p in 1:P) {
  # Crear un nuevo array para almacenar las variables retrasadas
  Yast_Lag <- array(NA, dim = dim(Y),dimnames=list(paste0("t", 1:T),
                                                paste0("Var", 1:K),
                                                paste0("Reg", 1:N)))
  
  # Aplicar el retraso en cada dimensión
  for (t in 1:(T-p)) {
    Yast_Lag[t + p, , ] <- Y.ast[t, , ]
  }
  assign(paste0("Yastlag_", p), Yast_Lag)
}

# Crear una lista para almacenar los dataframes por cada [,k,n]
df_list <- list()


for (k in 1:K) {
  for (n in 1:N) {
    # Crear un dataframe combinando los arrays para [,k,n]
    df <- data.frame(Y = Y[, k, n], 
                     Y.ast = Y.ast[, k, n])
    
    # Agregar las columnas para los lags de Y
    for (p in 1:P) {
      df[paste0("Ylag_", p)] <- get(paste0("Ylag_", p))[, k, n]  # Acceder directamente a los valores de Y
    }
    
    # Agregar las columnas para los lags de Y.ast
    for (p in 1:P) {
      df[paste0("Yastlag_", p)] <- get(paste0("Yastlag_", p))[, k, n]  # Acceder a cada Yastlag_p
    }
    
    # Eliminar filas con NA
    df <- na.omit(df)
    
    # Guardar el dataframe en la lista
    df_list[[paste("k", k, "n", n, sep = "_")]] <- df
  }
}

modelo_list=list()
for (k in 1:K) {
  for (n in 1:N) {
    # Seleccionar el dataframe correspondiente a [,k,n] de la lista df_list
    df <- df_list[[paste("k", k, "n", n, sep = "_")]]
    
    # Crear una fórmula dinámica para incluir las variables lag
    formula <- as.formula(paste("Y ~", 
                                paste0("Ylag_", 1:P, collapse = "+"), "+", 
                                paste0("Y.ast +", 
                                       paste0("Yastlag_", 1:P, collapse = " + "), 
                                       collapse = " + ")))
    
    # Ajustar el modelo lineal utilizando el dataframe df
    modelo <- lm(formula, data = df)
    
    # Almacenar los coeficientes del modelo en la lista
    modelo_list[[paste("k", k, "n", n, sep = "_")]] <- coef(modelo)
  }
}
