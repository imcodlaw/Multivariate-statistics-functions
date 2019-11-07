# function PCA, dengan parameter data dan standardize
# parameter data adalah input dari function, input dapat berupa matriks dan dataframe
# parameter standardize adalah untuk menentukan apakah data perlu distandarisasi atau tidak
# jika tidak memerlukan standarisasi, isi parameter standardize == FALSE, atau standardize == TRUE jika
# memerlukan proses standarisas
pca = function(data, standardize = FALSE){
  # untuk menstandarisasi data
  if(standardize == TRUE){
    data = scale(data)
    # outputnya adalah data yang sudah distandarisasi
  }
  # mencari nilai rata-rata dari tiap kolom
  mean = c() # vektor rata-rata
  vektor = c() # vektor nilai x dikurang dengan rata-rata tiap kolomnya
  n_row = nrow(data) # banyak baris
  n_col = ncol(data) # banyak kolom
  for (i in (1:n_col)){
    mean = c(c(mean), c((sum(data[,i]))/n_row))
  }
  # hasil dari looping di atas adalah vektor rata-rata untuk tiap variabel
  #-------------------------------------------------------------------------
  # PROSES PERHITUNGAN SEMI MANUAL (KARENA MENGGUNAKAN RUMUS TAPI OTOMATIS, JADINYA SEMI MANUAL HEHEHEHEHE)
  # 1. melakukan operasi pengurangan antara matriks dengan matriks rata-ratanya
  for (i in (1:n_row)){
    for(j in (1:n_col)){
      vektor = c(c(vektor), c(data[i,j] - mean[j]))
    }
  }
  #-------------------------------------------------------------------------
  # 2. melakukan perkalian matriks antara transpose matriks hasil perhitungan di atas
  #    dengan matriks tersebut
  mat = matrix(vektor, nrow = n_row, ncol = n_col, byrow = TRUE)
  S = (t(mat)%*%mat)/(n_row-1)
  # untuk data yang sudah distandarisasi, hasil matriks varians kovarians adalah berupa matriks korelasi
  # sedangkan untuk data yang belum distandarisasi, hasilnya berupa matriks varians kovarians
  #-------------------------------------------------------------------------
  # memperoleh nilai eigen dari matriks varians kovarians
  eigen_value = eigen(S)$values
  eigen_vector = eigen(S)$vector
  #-------------------------------------------------------------------------
  # mencari proporsi dari total populasi varians
  nilai_proporsi = c()
  for (i in (1:length(eigen_value))){
    nilai_proporsi = c(c(nilai_proporsi), c(eigen_value[i]/sum(eigen_value)))
  }
  #-------------------------------------------------------------------------
  # mencari matriks korelasi variabel lama dengan variabel baru
  korelasibaru = c()
  for (i in (1:nrow(eigen_vector))){
    for (j in (1:nrow(eigen_vector))){
      korelasibaru = c(c(korelasibaru), c(eigen_vector[i,j]*sqrt(eigen_value[j])/S[j,j]))
    }
  }
  # ini intinya masukin nilai korelasi tersebut ke dalam matriks
  mat_korelasi_baru = matrix(korelasibaru, nrow = nrow(eigen_vector), ncol = ncol(eigen_vector)
                             , byrow = TRUE)
  #-------------------------------------------------------------------------
  # menentukan persamaan variabel baru
  varbaru = matrix(eigen_vector, nrow = nrow(eigen_vector), ncol = ncol(eigen_vector), byrow =TRUE)
  hasil = list(mat = mat, covariance = S, eigen_value = eigen_value, eigen_vector = eigen_vector, persamaan_variabel_baru = varbaru, 
               proporsi = nilai_proporsi, new_correlation = mat_korelasi_baru)
  # return(mat_korelasi_baru)
  return(hasil)
}


# cov(c)

# data = read.csv("D:/tugas3.csv")