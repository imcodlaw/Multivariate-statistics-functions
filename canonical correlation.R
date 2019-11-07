kk = function(data, p, q, standardize = FALSE){
  # parameter data untuk menyatakan data input
  # parameter p untuk menyatakan banyak variabel 1
  # parameter q untuk menyatakan banyak variabel 2
  # dengan asumsi, matriks data p berurutan kemudian matriks q, p <= q
  if(standardize == TRUE){
    data = scale(data)
  }
  # menghitung matriks S dan partisinya
  S = cov(data)
  s11 = S[1:p, 1:p]
  s12 = S[1:p, (p+1):(p+q)]
  s21 = S[(p+1):(p+q), 1:p]
  s22 = S[(p+1):(p+q), (p+1):(p+q)]
  
  # mencari nilai s11^(-1/2)
  s11_12 = (1/sqrt(eigen(s11)$values[1]))*eigen(s11)$vector[,1]%*%t(eigen(s11)$vector[,1])+
    (1/sqrt(eigen(s11)$values[2]))*eigen(s11)$vector[,2]%*%t(eigen(s11)$vector[,2])
  
  # rhofix merupakan perkalian antara rho11, rho12, rho21, rho22
  rhofix = s11_12%*%s12%*%(solve(s22))%*%s21%*%s11_12
  
  # eigen value dari perkalian nilai rho
  # eigen value ini berguna untuk melakukan signifikansi parameter?
  # detail tentang ujinya ada di pdf
  eigen_value = eigen(rhofix)$values
  
  # nilai a1 merupakan nilai koefisien variabel pada u1
  hasil_a = c()
  hasil_b = c()
  for(i in (1:min(p,q))){
    a = s11_12%*%eigen(rhofix)$vector[,i]
    hasil_a = c(c(hasil_a), c(a))
    
    b1 = (solve(s22))%*%s21%*%a
    b2 = t(b1)%*%s22%*%b1
    b11 = 1/sqrt(b2[1])*b1
    hasil_b = c(c(hasil_b), c(b11)) 
    # b1 = (solve(s22))%*%s21%*%a1
  }
  A = matrix(c(hasil_a), nrow = min(p,q), ncol = p, byrow = FALSE)
  B = matrix(c(hasil_b), nrow = min(p,q), ncol = q, byrow = FALSE)
  
  # r = t(A)%*%s12%*%B/(sqrt(t(A)%*%s11%*%A * sqrt(t(B)%*%s22%*%B)))
  # hasil_r = c()
  # for(i in (1:min(p,q))){
  #   vu = t(A[,i])%*%s11%*%A[,i]
  #   vv = t(B[,i])%*%s22%*%B[,i]
  #   cov_uv = t(A[,i])%*%s12%*%B[,i]
  #   r = cov_uv/(sqrt(vu)*sqrt(vv))
  #   hasil_r = c(c(hasil_r), c(r))
  # }
  # mat_r = matrix(c(hasil_r), nrow = min(p,q), ncol = 1, byrow = FALSE)
  
  # # vu1 merupakan nilai varians dari variabel u1
  # vu1 = t(a1)%*%s11%*%a1
  # 
  # # vv1 merupakan nilai varians dari variabel v1
  # vv1 = t(b11)%*%s22%*%b11
  # 
  # # cov_u1v1 merupakan nilai kovarians dari variabel u1 dengan v1
  # cov_u1v1 = t(a1)%*%s12%*%b11
  # 
  # # r1 merupakan nilai korelasi kanonik antara u1 dengan v1
  # r1 = cov_u1v1/(sqrt(vu1)*sqrt(vv1))
  # 
  # # vu2 merupakan nilai varians variabel u2
  # vu2 = t(a2)%*%s11%*%a2
  # 
  # # vv2 merupakan nilai varians variabel v2
  # vv2 = t(b21)%*%s22%*%b21
  # 
  # # cov_u2v2 merupakan nilai kovarians dari variabel u2 dengan v2
  # cov_u2v2 = t(a2)%*%s12%*%b21
  # 
  # # r2 merupakan nilai korelasi kanonik antara u2 dengan v2
  # r2 = cov_u2v2/(sqrt(vu2)*sqrt(vv2))
  # 
  # # r1. dan r2. merupakan nilai untuk mencari seberapa (isi titik-titik)
  # r1. = r1/(r1+r2)
  # r2. = r2/(r1+r2)
  
  # intinya .........
  # cor_ux1 = t(a1)%*%s11
  # cor_ux2 = t(a2)%*%s12
  # cor_vx1 = t(b11)%*%s22
  # cor_vx2 = t(b21)%*%s21
  
  # return(list("A" = A, "B" = B, "R1" = r1, 
  #             "R2" = r2, eigen_value, 'cor_ux1' = cor_ux1, 'cor_ux2' = cor_ux2, 'cor_vx1'
  #             = cor_vx1, 'cor_vx2' = cor_vx2, "r1." = r1., "r2." = r2.))
  return(list("A" = A, "B" = B, "eigen_value" = eigen_value))
}