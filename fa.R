# analisis faktor

# 1. cari matriks var kovar
# 2. eigen value and eigen vector
# 3. menentukan jumlah faktor
# 4. menaksir L
# 5. menentukan factor score
data = matrix(c(1,6,9,4,12,10,3,12,15,4,10,12), nrow = 4, ncol = 3, byrow = 3)
data

# function fa stands for factor analysis
factorAnalysis = function(data, standardize = FALSE){
  # process to standardize the data
  if (standardize == TRUE){
    data = scale(data)
  }
  #----------------------------------------------------------------------------------------
  # finding covariance matrix of the data
  S = cov(data)
  #----------------------------------------------------------------------------------------
  # finding the eigen value and eigen vector of S
  eigen_value = eigen(S)$values
  eigen_vector = eigen(S)$vector
  #----------------------------------------------------------------------------------------
  # estimating the loading score (L) WITH PRINCIPAL COMPONENT METHOD
  # lambda is an diagonal matrix with eigen values from covariance matrix as the element of 
  # the main diagonal
  lambda = diag(eigen_value, nrow = length(eigen_value), ncol = length(eigen_value))
  L = eigen_vector%*%sqrt(lambda)
  #----------------------------------------------------------------------------------------
  # estimating F score
  # 1. first, calculate the centralized X
  mean_vectors = colMeans(data)
  central = c()
  for (i in (1:nrow(data))){
    for (j in (1:ncol(data))){
      central = c(c(central), c(data[i,j]-mean_vectors[j]))
    }
  }
  Xc = matrix(central, nrow = nrow(data), ncol = ncol(data), byrow = TRUE)
  # next, estimate the F score with the equation F = Xc * S inverse * L
  F = Xc%*%solve(S)%*%L
  biplot(princomp(data))
  hasil = list(covariance = S, eigen_value = eigen_value, eigen_vector = eigen_vector, 
               lambda = lambda, loading_score = L, centralized_data = Xc, F_score = F)
  return (hasil)
}
#------------------------------------------------------------------------------------------
S = cov(data)
eigen_vector = eigen(S)$vector
eigen_value = eigen(S)$values
eigen_vector

V = eigen_vector[,1:2]
V


cor = matrix(c(1, 0.55, 0.58, 0.42, 0.77, 0.29, 0.55, 1, 0.49, 0.48, 0.63, 0.13
               , 0.58, 0.49, 1, 0.29, 0.81, 0.22, 0.42, 0.48, 0.29, 1, 0.29,
               -0.28, 0.77, 0.63, 0.81, 0.29, 1, 0.24, 0.29, 0.13, 0.22, -0.28,
               0.24, 1), nrow = 6, ncol = 6, byrow = TRUE)
cor
e_val_cor = eigen(cor)$values
e_vec_cor = eigen(cor)$vector


lambda = diag(e_val_cor, nrow = 6, ncol = 6, names = TRUE)
lambda
lambda1 = diag(c(3.22,1.3,0.58,0.43,0.33,0.1), nrow = 6, ncol = 6, names = TRUE)
lambda1

L = e_vec_cor%*%sqrt(lambda)
L


L1 = e_vec_cor%*%sqrt(lambda1)

hasil = c()
for (i in (1:nrow(data))){
  for (j in (1:ncol(data))){
    hasil = c(c(hasil), c(data[i,j]-mean(data[,j])))
  }
  mat = matrix(hasil, nrow = 4, ncol = 3, byrow = TRUE)
  print(mat)
}
pca(data)$mat
L1

sqrt(3.22*0.48)

data1 = read.csv('F:/data1.csv')
data1 = data1[-9]
data1 = as.matrix(data1)

Xc = pca(data1, standardize = TRUE)$mat
S1 = pca(data1, standardize = TRUE)$covariance
S1
e_val_S1 = eigen(S1)$values
e_val_S1
e_vec_S1 = eigen(S1)$vector
e_vec_S1
lambda_S1 = diag(e_val_S1, nrow = 11, ncol = 11, names = TRUE)
L2 = e_vec_S1%*%sqrt(lambda_S1)
L2
F2 = Xc%*%solve(S1)%*%L2
F2

Sdb = cov(data)
eigen(Sdb)
# -------------------------------------------------------------------------
# test pake data mr ireland
dataTest = read.csv("F:/pcairlan.csv")
dataTest = data.frame(dataTest[,-1], row.names = dataTest[,1])
dataTest

fa(dataTest, standardize = TRUE)
vektorEigen = matrix(c(0.48, 0.43, 0.46, 0.29, 0.51, 0.15, 0.06, -0.15, 0.11, -0.63, 0.11, 0.74)
                     , nrow = 2, ncol = 6, byrow = TRUE)
vektorEigen
nilaiEigen = diag(c(3.22,1.3,0.59,0.44,0.34,0.11))
nilaiEigen
LTest = vektorEigen%*%sqrt(nilaiEigen)
Ftest = 