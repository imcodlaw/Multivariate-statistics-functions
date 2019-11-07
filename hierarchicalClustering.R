hc = function(data, method){
  row1 = rownames(data)
  D = c()
  for (i in (1:nrow(data))){
     for (j in (1:nrow(data))){
       hasil = sqrt(sum((data[i,] - data[j,])^2))
       D = c(c(D), c(hasil)) 
     }
   }
  D = matrix(D, nrow = nrow(data), ncol = nrow(data), byrow = FALSE)
  dimnames(D) = list(c(row1), c(row1))
  
  step = 1
  if(method == "single"){
    n = length(D)
    while(n > 1){
      # untuk mencari nilai yang terkecil dari matriks D
      result = which(D == min(D[D > 0]), arr.ind = TRUE)
      result = result[2,]
      # -------------------------------------------------------------------------
      
      # looping buat mencari nilai jarak antara cluster dengan pengamatan di
      # luar cluster
      hasil = c()
      for(i in (1:(ncol(D)))){
        if((D[result[1],i]) == 0 || (D[result[2], i]) == 0){
          next
        }
        else{
          sel = min((D[result[1], i]), (D[result[2], i]))
          hasil = c(c(hasil), c(sel))
        }
      }
      varsisa = c(1:ncol(D))
      varsisa = setdiff(varsisa, result)
      # -------------------------------------------------------------------------
      
      # looping untuk mencari jarak antar observasi yang di luar cluster
      for(i in (1:length(varsisa))){
        for(j in (1:length(varsisa))){
          if(i == j || i >= j){
            next
          }
          else{
            nilai = D[varsisa[i],varsisa[j]]
            hasil = c(c(hasil), c(nilai))
          }
        }
      }
      # -------------------------------------------------------------------------
      
      # proses pembentukan dan penamaan pada matriks
      kolom = colnames(D) # untuk memperoleh nama kolom
      
      # -------------------------------------------------------------------------
      
      # pembentukan matriks 
      D = diag(0, (nrow(D)-1), (ncol(D)-1))
      D[lower.tri(D)] = hasil
      D[upper.tri(D)] = t(D)[upper.tri(D)]
      # -------------------------------------------------------------------------
      
      # penamaan matriks 
      dimnames(D) = list(c(paste(kolom[result], collapse = " ") , kolom[-result]), 
                         c(paste(kolom[result], collapse = " ") , kolom[-result]))
      n = length(D)
      print(paste(c("iterasi ke-", step), collapse = ""))
      print(D)
      
      step = step + 1
      print("-------------------------------------------------------------------------")
      # -------------------------------------------------------------------------
    } # akhir while
    # return(D)
    dendogram = dist(data, method = "euclidean")
    dendogram = hclust(dendogram, method = "single")
    plot(dendogram)
  } # ini akhir method
  else if(method == "complete"){
    n = length(D)
    while(n > 1){
      # untuk mencari nilai yang terkecil dari matriks D
      result = which(D == min(D[D > 0]), arr.ind = TRUE)
      result = result[2,]
      # -------------------------------------------------------------------------
      
      # looping untuk mencari nilai jarak antara cluster dengan pengamatan di
      # luar cluster
      hasil = c()
      for(i in (1:(ncol(D)))){
        if((D[result[1],i]) == 0 || (D[result[2], i]) == 0){
          next
        }
        else{
          sel = max((D[result[1], i]), (D[result[2], i]))
          hasil = c(c(hasil), c(sel))
        }
      }
      varsisa = c(1:ncol(D))
      varsisa = setdiff(varsisa, result)
      # -------------------------------------------------------------------------
      
      # looping untuk mencari jarak antar observasi di luar cluster
      for(i in (1:length(varsisa))){
        for(j in (1:length(varsisa))){
          if(i == j || i >= j){
            next
          }
          else{
            nilai = D[varsisa[i],varsisa[j]]
            hasil = c(c(hasil), c(nilai))
          }
        }
      }
      # -------------------------------------------------------------------------
      
      # pembentukan matriks 
      kolom = colnames(D)
      D = diag(0, (nrow(D)-1), (ncol(D)-1))
      D[lower.tri(D)] = hasil
      D[upper.tri(D)] = t(D)[upper.tri(D)]
      # -------------------------------------------------------------------------
      
      # penamaan matriks 
      dimnames(D) = list(c(paste(kolom[result], collapse = "") , kolom[-result]), 
                         c(paste(kolom[result], collapse = "") , kolom[-result]))
      # -------------------------------------------------------------------------

      n = length(D)
      print(paste(c("iterasi ke-", step), collapse = ""))
      print(D)
      
      step = step + 1
      print("-------------------------------------------------------------------------")
    } # akhir while
    # return(D)
    dendogram = dist(data, method = "euclidean")
    dendogram = hclust(dendogram, method = "complete")
    plot(dendogram)
  } # akhir method complete
  else if(method == "average"){
    n = length(D)
    while(n > 1){
      # untuk mencari nilai yang terkecil dari matriks D
      result = which(D == min(D[D > 0]), arr.ind = TRUE)
      result = result[2,]
      # -------------------------------------------------------------------------
      
      # looping buat mencari nilai jarak antara cluster dengan pengamatan di
      # luar cluster
      hasil = c()
      for(i in (1:(ncol(D)))){
        if((D[result[1],i]) == 0 || (D[result[2], i]) == 0){
          next
        }
        else{
          sel = (D[result[1], i] + D[result[2], i])/2
          hasil = c(c(hasil), c(sel))
        }
      }
      varsisa = c(1:ncol(D))
      varsisa = setdiff(varsisa, result)
      # -------------------------------------------------------------------------
      
      # looping untuk mencari jarak antar observasi yang di luar cluster
      for(i in (1:length(varsisa))){
        for(j in (1:length(varsisa))){
          if(i == j || i >= j){
            next
          }
          else{
            nilai = D[varsisa[i],varsisa[j]]
            hasil = c(c(hasil), c(nilai))
          }
        }
      }
      # -------------------------------------------------------------------------
      
      # proses pembentukan dan penamaan pada matriks
      kolom = colnames(D) # untuk memperoleh nama kolom
      
      # -------------------------------------------------------------------------
      
      # pembentukan matriks 
      D = diag(0, (nrow(D)-1), (ncol(D)-1))
      D[lower.tri(D)] = hasil
      D[upper.tri(D)] = t(D)[upper.tri(D)]
      # -------------------------------------------------------------------------
      
      # penamaan matriks 
      dimnames(D) = list(c(paste(kolom[result], collapse = "") , kolom[-result]), 
                         c(paste(kolom[result], collapse = "") , kolom[-result]))
      n = length(D)
      print(paste(c("iterasi ke-", step), collapse = ""))
      print(D)
      
      step = step + 1
      print("-------------------------------------------------------------------------")
      # -------------------------------------------------------------------------
    } # akhir while
    # return(D)
    dendogram = dist(data, method = "euclidean")
    dendogram = hclust(dendogram, method = "average")
    plot(dendogram)
  } # ini akhir method average
  else if(method == "centroid"){
    D = D^2
    n = length(D)
    while(n > 1){
      # untuk mencari nilai yang terkecil dari matriks D
      result = which(D == min(D[D > 0]), arr.ind = TRUE)
      result = result[2,]
      # -------------------------------------------------------------------------
      
      kolom = colnames(D) # untuk memperoleh nama kolom
      
      # looping buat mencari nilai jarak antara cluster dengan pengamatan di
      # luar cluster
      hasil = c()
      for(i in (1:(ncol(D)))){
        if((D[result[1],i]) == 0 || (D[result[2], i]) == 0){
          next
        }
        else{
          sel = (length(result[1])*D[result[1],i] + length(result[2])*D[result[2],i] - length(result[1])*length(result[2])*D[result]/sum(nchar(kolom[result])))/(sum(nchar(kolom[result])))
          hasil = c(c(hasil), c(sel))
        }
      }
      
      varsisa = c(1:ncol(D))
      varsisa = setdiff(varsisa, result)
      # -------------------------------------------------------------------------
      
      # looping untuk mencari jarak antar observasi yang di luar cluster
      for(i in (1:length(varsisa))){
        for(j in (1:length(varsisa))){
          if(i == j || i >= j){
            next
          }
          else{
            nilai = D[varsisa[i],varsisa[j]]
            hasil = c(c(hasil), c(nilai))
          }
        }
      }
      # -------------------------------------------------------------------------
      
      # proses pembentukan dan penamaan pada matriks
      # kolom = colnames(D) # untuk memperoleh nama kolom
      
      # -------------------------------------------------------------------------
      
      # pembentukan matriks 
      D = diag(0, (nrow(D)-1), (ncol(D)-1))
      D[lower.tri(D)] = hasil
      D[upper.tri(D)] = t(D)[upper.tri(D)]
      # -------------------------------------------------------------------------
      
      # penamaan matriks 
      dimnames(D) = list(c(paste(kolom[result], collapse = "") , kolom[-result]), 
                         c(paste(kolom[result], collapse = "") , kolom[-result]))
      n = length(D)
      print(paste(c("iterasi ke-", step), collapse = ""))
      print(D)
      
      step = step + 1
      print("-------------------------------------------------------------------------")
      # -------------------------------------------------------------------------
    } # akhir while
    dendogram = dist(data, method = "euclidean")
    dendogram = hclust(dendogram, method = "centroid")
    plot(dendogram)
  } # end of method
} # end of code
