
mdsnon = function(data, k){
  # masukin stress di parameter kalo udah beres coba-coba
  row1 = rownames(data)
  stress_hat = 1
  iterasi = 1
  D0 = c()
  for (i in (1:nrow(data))){
    for (j in (1:nrow(data))){
      hasil = sqrt(sum((data[i,] - data[j,])^2))
      D0 = c(c(D0), c(hasil))
    }
  } 
  D0 = matrix(D0, nrow = nrow(data), ncol = nrow(data), byrow = FALSE)
  dimnames(D0) = list(c(row1), c(row1))
  while(k > 0){
    # syarat while diganti sama yang di bawah, kalo udah beres:
    # stress_hat >= stress
    
    # untuk memperoleh matriks distance (D) dengan euclidean
    # --------------------------------------------------------------------
    if(iterasi == 1){
      D = D0
    }
    # else{
    #   D = c()
    #   for (i in (1:nrow(data))){
    #     for (j in (1:nrow(data))){
    #       hasil = sqrt(sum((data[i,] - data[j,])^2))
    #       D = c(c(D), c(hasil))
    #     }
    #   }
    #   D = matrix(D, nrow = nrow(data), ncol = nrow(data), byrow = FALSE)
    #   dimnames(D) = list(c(row1), c(row1))
    # }
    
    # --------------------------------------------------------------------
    # D = data
    # dimnames(D) = list(c(row1), c(row1)) 
    
    # CLEAR!
    
    # --------------------------------------------------------------------
    # untuk memperoleh sum of square
    A = -0.5*(D^2)
    ai. = c()
    # a.j = c()
    for(i in (1 : nrow(D))){
      di. = (sum(A[i,]))/(ncol(A))
      ai. = c(c(ai.), c(di.))
      # --------------------------------------------------------------------
    }
    # untuk menghasilkan matriks B
    dsquare = sum(A)/(length(A))
    elemen = c()
    for(i in (1:(ncol(A)))){
      for(j in (1:(nrow(A)))){
        # hasil = -0.5*((D^2)[j,i]) - h.j[j] - hi.[i] + dsquare
        hasil = A[j,i]-ai.[j]-ai.[i]+dsquare
        # -4 + 7.2 + 3.2 - 6.4
        # A[2,1] - ai.[2] - ai.[1] + dsquare
        elemen = c(c(elemen), c(hasil))
        # print(A[j,i]-ai.[j]-ai.[i]+dsquare)
      }
    }
    B = matrix(c(elemen), nrow = nrow(D), ncol = ncol(D), byrow = FALSE)
    rownames(B) = row1
    
    # CLEAR!
    
    # --------------------------------------------------------------------
    
    # mencari eigen value dan vector untuk matriks koordinat
    eigen_vector = eigen(B)$vectors
    eigen_value = eigen(B)$values
    
    F_hasil = matrix(c(c(sqrt(eigen_value[1]) * eigen_vector[,1]), c(sqrt(eigen_value[2]) 
    * eigen_vector[,2])), nrow = nrow(B), ncol = 2, byrow = FALSE)
    rownames(F_hasil) = row1
    
    # CLEAR!
    
    # --------------------------------------------------------------------
    
    # untuk memperoleh matriks D(F)
    D_hat = c()
    for (i in (1:nrow(F_hasil))){
      for (j in (1:nrow(F_hasil))){
        hasil = sqrt(sum((F_hasil[i,] - F_hasil[j,])^2))
        D_hat = c(c(D_hat), c(hasil))
      }
    }
    D_hat = matrix(D_hat, nrow = nrow(data), ncol = nrow(data), byrow = FALSE)
    dimnames(D_hat) = list(c(row1), c(row1))
    
    # CLEAR!
    
    # --------------------------------------------------------------------
    
    stress_hat = sqrt(sum((D - D_hat)^2)/(sum(D0^2)))
    
    
    # CLEAR !
    # --------------------------------------------------------------------
    D = D_hat
    
    print(paste(c(c("iterasi ke-"), c(iterasi)), collapse = ""))
    print("ini nilai A")
    print(A)
    print("--------------------------------------------------------------------")
    print("ini nilai B")
    print(B)
    print("--------------------------------------------------------------------")
    print("ini nilai F")
    print(F_hasil)
    print("--------------------------------------------------------------------")
    print("ini nilai stress hat")
    print(stress_hat)
    print("--------------------------------------------------------------------")
    iterasi = iterasi + 1
    
    k = k - 1
  
  }
  plot(F_hasil)
  text(F_hasil, labels = row1, cex = 0.6, pos = 4)
  return(F_hasil)
}
# 0.0015
kota = read.csv("kota.csv")
rownames(kota)
kota = kota[,-1]
kota
kota = matrix(as.numeric(unlist(kota)),nrow=nrow(kota))
rownames(kota) = nama

mdsnon(kota, 5) 

plot(F_hasil)
text(F_hasil, labels = nama, cex = 0.6, pos = 3)
