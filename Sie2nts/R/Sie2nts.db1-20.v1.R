# Daubechies Wavelet 1 and some general functions to be used for Daubechies, the basis of wavelet is formulated by the method of Meyer (S.2 in the paper)
# However, we also try the S.1 for estimating coefficients and choose the J_0 = 0.
library(ggplot2)


# find the minimal value correspond, inter is the interval created, valtable is the db2 value table which is generated in advance.upper is the upper bound of basis. val express the input x.
w_find = function(valtable, inter, upper, val){
  if(val < 0 | val >= upper){
    return(0)
  } else{
    return(valtable[which(abs(inter-val) == min(abs(inter-val)))])
  }
}

db1_f = function(t){
  return(ifelse(t<1 & t>=0, 1,0))
}


dbplot = function(w, ops, title){ # c indicate the order of db and w indicate the number of basis(the true basis is 2^w).
  if(ops == "db1"){
    point = 10000
    dbt = valdb(w,psi.f = 0, 1, point)
    n = dim(dbt)[1]
    x = seq(0, 1, length=n)

    df = data.frame()
    for (i in 2:(2^w + 1)){
      res = as.data.frame(dbt[, i]) #unlist(poly_val(coeffi, i))
      df = rbind(df, res)
    }
    f.df = df
    colnames(f.df) = "new"
    f.df$x = rep(x, 2^w)
    f.df$order = as.factor(rep(1:(2^w), each = n))
    theme_update(plot.title = element_text(hjust = 0.5))
    p1 = ggplot(f.df, aes(x=x, y=new, group=order, colour = order))+ geom_line()  + ggtitle(title) +
      xlab("") + ylab("") + scale_colour_discrete(name  ="order")+theme(plot.title = element_text(size=18, face="bold"),
                                                                        legend.text=element_text(size=24, face = "bold"),
                                                                        axis.text.x = element_text(face="bold", color="#993333",                                                                                                                              size=22, angle=0),
                                                                        axis.text.y = element_text(face="bold", color="#993333",size=22, angle=0),
                                                                        axis.title.x=element_text(size=22,face='bold'),
                                                                        axis.title.y=element_text(angle=90, face='bold', size=22),
                                                                        legend.title = element_text(face = "bold"))
    return(p1)
  } else{
    library(stringr)
    filename = paste(ops,"_fa_table", sep = "")
    aux_str = str_split(ops, "")[[1]]
    if(aux_str[1] == "d"){
      if(length(aux_str) == 3){
        c = as.numeric(aux_str[3])
      }else
      {
        c = as.numeric(paste(aux_str[3], aux_str[length(aux_str)], sep = ""))
      }

    } else{
      c = as.numeric(aux_str[3])*3
    }

    library(RCurl)

    x <- getURL(paste("https://raw.githubusercontent.com/xcding1212/Sie2nts/main/db_table/", filename, sep = ""))
    dbt = read.csv(text = x)

    dbt = dbt[,2]

    df = data.frame()
    t=seq(0, 2*c - 1, length.out = length(dbt))
    n = length(t)
    x =seq(0,1, length.out = n)
    for(h in 0:(2^w-1)){
      p=rep(0, n)
      for (k in 1:n){
        for (l in -70:70){
          p[k]=p[k]+ w_find(dbt, t, 2*c - 1, (2^(w)*(x[k]+l)-h))
        }
        p[k]=2^(w/2)*p[k]
      }

      df=rbind(df, data.frame(value = p))
      p=rep(0, n)
    }
    df$x = rep(x, dim(df)[2])
    df$class = as.factor(rep(1:(2^w), each = n))
    theme_update(plot.title = element_text(hjust = 0.5))
    p1 = ggplot(df, aes(x=x, y=value, group=class, colour = class))+ geom_line()  + ggtitle(title) +
      xlab("") + ylab("") + scale_colour_discrete(name  ="order")+theme(plot.title = element_text(size=18, face="bold"),
                                                                        legend.text=element_text(size=24, face = "bold"),
                                                                        axis.text.x = element_text(face="bold", color="#993333",                                                                                                                              size=22, angle=0),
                                                                        axis.text.y = element_text(face="bold", color="#993333",size=22, angle=0),
                                                                        axis.title.x=element_text(size=22,face='bold'),
                                                                        axis.title.y=element_text(angle=90, face='bold', size=22),
                                                                        legend.title = element_text(face = "bold"))
    return(p1)
  }
}


wavelet_kth_b = function(k, ops){ # the k-th basis in 2^k total basis
  w = k
  if(ops == "db1"){
    point = 10000
    dbt = valdb(w,psi.f = 0, 1, point)
    n = dim(dbt)[1]
    x = seq(0, 1, length=n)

    df = data.frame(x)
    for (i in 2:(2^w + 1)){
      res = as.data.frame(dbt[, i]) #unlist(poly_val(coeffi, i))
      df = cbind(df, res)
    }
    df = data.frame(basis_value = df[,k+1])
    return(df)
  } else{
    library(stringr)
    filename = paste(ops,"_fa_table", sep = "")
    aux_str = str_split(ops, "")[[1]]
    if(aux_str[1] == "d"){
      if(length(aux_str) == 3){
        c = as.numeric(aux_str[3])
      }else
      {
        c = as.numeric(paste(aux_str[3], aux_str[length(aux_str)], sep = ""))
      }

    } else{
      c = as.numeric(aux_str[3])*3
    }

    library(RCurl)
    x <- getURL(paste("https://raw.githubusercontent.com/xcding1212/Sie2nts/main/db_table/", filename, sep = ""))
    dbt = read.csv(text = x)
    dbt = dbt[,2]
    t=seq(0, 2*c - 1, length.out = length(dbt))
    n = length(t)
    x =seq(0,1, length.out = n)
    df = data.frame(x)
    for(h in 0:(2^w-1)){
      p=rep(0, n)
      for (k in 1:n){
        for (l in -70:70){
          p[k]=p[k]+ w_find(dbt, t, 2*c - 1, (2^(w)*(x[k]+l)-h))
        }
        p[k]=2^(w/2)*p[k]
      }

      df=cbind(df, data.frame(value = p))
      p=rep(0, n)
    }

    # db.leng = seq(0,1, length.out = point)
    df = data.frame(basis_value = df[,w+1])
    return(df)
  }

}



# This basis is formulated by Meyer for Daubechies1
valdb1 = function(w, n){
  x =seq(0,1, length.out = n)
  df.db1 = data.frame(x = x)

  for(h in 0:(2^w-1)){
    p=rep(0, n);
    for (k in 1:n){
      for (l in 0:70){
        p[k]=p[k]+ db1_f((2^(w)*(x[k]+l)-h))
      }
      p[k]=2^(w/2)*p[k]
    }
    df.db1[,h+2] = p
  }
  return(df.db1)
}


# get the res with different w.

# db_number represent which order of Daubechies to be used. 1-20 right now could be chosen. w indicates the number of basis functions.

valdb = function(w, psi.f=0, db_number, len.n){ # len.n indicates the point chosen
  if(db_number == 1){
    return(valdb1(w, len.n))
  } else{
    n = length(psi.f)
    x =seq(0,1, length.out = n)
    df.db = data.frame(x)
    t=seq(0, 2*db_number - 1, length.out = n)

    for(k in 0:(2^w-1)){
      p=rep(0, n);
      for (ind in 1:n){
        for (l in 0:70){
          p[ind]=p[ind]+ w_find(psi.f, t, 2*db_number - 1, (2^(w)*(x[ind]+l)-k))
        }
        p[ind]=(2^(w/2))*p[ind]
      }
      df.db[,k+2] = p
    }
    return(df.db)
  }

}



phi_db = function(bspline, beta, b){
  c = length(bspline)
  b_res = list()
  for(i in 0:b){
    B.aux = matrix(c(rep(0, c*i), bspline, rep(0, c*(b-i))), ncol = 1)
    b_res[[i+1]] = as.numeric(t(beta)%*%B.aux)
  }

  return(b_res)
}

# For resolving the issue system is computationally singular, we need to put options tol in the function solve(). One issue, why the last one coefficient inflated.

simu_db = function(ts, df.db, b, m = 500){
  l.alpha = list()
  aux.alpha = c()

  n <- length(ts)
  x <- seq(0, 1, length=n)

  Base = df.db[which(abs(df.db$x - x[1]) == min(abs(df.db$x - x[1]))), -1]
  for(i in 2:length(x)){
    Base = rbind(Base, df.db[which(abs(df.db$x - x[i]) == min(abs(df.db$x - x[i]))), -1])
  }
  B = as.matrix(Base)

  aux_B = B[(b+1):n,]
  ind = b
  aux.Y = matrix(rep(ts[ind:(n-b+ind-1)], dim(aux_B)[2]), ncol = dim(aux_B)[2])
  Y = cbind(aux_B, aux_B*aux.Y)
  ind = ind - 1
  while(ind >= 1){
    aux.Y = matrix(rep(ts[ind:(n-b+ind-1)], dim(aux_B)[2]), ncol = dim(aux_B)[2])
    Y = cbind(Y, aux_B*aux.Y)
    ind = ind - 1
  }

  X = matrix(ts[(b+1):n], ncol = 1)
  beta = solve(t(Y)%*%Y, tol = 1e-20)%*%t(Y)%*%X

  x <- seq(0, 1, length = m)

  Base = df.db[which(abs(df.db$x - x[1]) == min(abs(df.db$x - x[1]))), -1]
  for(i in 2:length(x)){
    Base = rbind(Base, df.db[which(abs(df.db$x - x[i]) == min(abs(df.db$x - x[i]))), -1])
  }
  B = as.matrix(Base)

  for(j in 1:(b+1)){
    for(i in 1:m){
      bspline = B[i,]
      aux.alpha[i] = phi_db(bspline, beta, b)[[j]]
    }
    l.alpha[[j]] = aux.alpha
    aux.alpha = c()
  }

  return(list(l.alpha, beta, Y))
}



db.loocv = function(ts, df.db, b){
  c = dim(df.db)[2] - 1
  n = length(ts)
  aux.true = ts[(b+1):n]
  aux.esti = c()
  leve.i = c()

  beta.es = simu_db(ts, df.db, b, n)
  hat = beta.es[[3]]%*%solve(t(beta.es[[3]])%*%beta.es[[3]])%*%t(beta.es[[3]])

  for(i in (b+1):n){
    aux.esti[i-b]  = matrix(beta.es[[3]][i-b,], nrow = 1)%*%beta.es[[2]]
    leve.i[i-b] = as.numeric(hat[i-b,i-b]) # hii is the diagonal of the hat matrix
  }
  error = (aux.true - aux.esti)^2
  lever = sum(error/((1-leve.i)^2))/(n-b)
  return(c(log(c,2),b,lever))
}


db.cv = function(ts, c, b, ops){
  n = length(ts)
  l = floor(3*log2(n))
  aux.train = ts[1:(n-l)]
  aux.vali = ts[(n-l+1):n]
  tt = fix.fit.wavelet(aux.train, c, b, length(aux.train), ops)
  pre = predict.wavelet(aux.train, tt, length(aux.vali))
  error = sum((aux.vali - pre)^2)/l
  return(c(c,b,error))
}


# estimate contains the coefficients, beta and the design matrix.

# prediction
predict.wavelet = function(ts, esti.li, k){ # k indicates the number of predictions
  ts.pre = c()
  phi.h = esti.li[[2]]
  n = length(ts)
  b = length(phi.h)

  for(h in 1:k){
    aux.pre = phi.h[[1]][n]
    for(j in 2:b){
      aux.pre = aux.pre + phi.h[[j]][n]*ts[n-h-j]
    }
    ts.pre[h] = aux.pre
  }
  return(ts.pre)
}

fix.fit.wavelet = function(ts, k, b, m = 500, ops){ # this k indicates that the number of basis is 2^k

  if(ops == "db1"){
    basis_db1 = valdb(k, psi.f=0, db_number=1,len.n=10000)
    ols = simu_db(ts, basis_db1, b, m)[[2]]
    aux.ts = simu_db(ts, basis_db1, b, m)[[1]]

    error.s = c()
    n = length(ts)
    es.alpha = simu_db(ts, basis_db1, b, n)[[1]]
    aux.len = length(es.alpha)
    for(i in (b+1):n){
      val.aux = es.alpha[[1]][i]
      for(j in 2:aux.len){
        val.aux = val.aux + es.alpha[[j]][i]*ts[i-j+1]
      }
      error.s[i-b] = ts[i] - val.aux
    }

    return(list(ols.coef = ols, ts.coef = aux.ts, Residuals = error.s))

  } else{
    library(stringr)
    filename = paste(ops,"_fa_table", sep = "")
    aux_str = str_split(ops, "")[[1]]
    if(aux_str[1] == "d"){
      if(length(aux_str) == 3){
        c = as.numeric(aux_str[3])
      }else
      {
        c = as.numeric(paste(aux_str[3], aux_str[length(aux_str)], sep = ""))
      }

    } else{
      c = as.numeric(aux_str[3])*3
    }

    library(RCurl)
    x <- getURL(paste("https://raw.githubusercontent.com/xcding1212/Sie2nts/main/db_table/", filename, sep = ""))
    dbt = read.csv(text = x)
    dbt = dbt[,2]
    db.basis = valdb(k, dbt, db_number = c, len.n = 0)

    ols = simu_db(ts, db.basis, b, m)[[2]] # first k indicates the k of 2^k number of basis, second c indicates db number, b indicates b.
    ts.c = simu_db(ts, db.basis, b, m)[[1]]

    error.s = c()
    n = length(ts)
    es.alpha = simu_db(ts, db.basis, b, n)[[1]]
    aux.len = length(es.alpha)
    for(i in (b+1):n){
      val.aux = es.alpha[[1]][i]
      for(j in 2:aux.len){
        val.aux = val.aux + es.alpha[[j]][i]*ts[i-j+1]
      }
      error.s[i-b] = ts[i] - val.aux
    }

    return(list(ols.coef = ols, ts.coef = ts.c, Residuals = error.s))

    }

}

auto.fit.wavelet = function(ts, c=3, b = 2, m=500, ops, method = "CV", threshold = 0){ # c indicates 2^c number of basis function

  res.bc = matrix(ncol = 3, nrow = c*b)
  ind = 1
  for(i in 1:c){
    for(j in 1:b){
      if(method == "CV"){
        res.bc[ind, ] = db.cv(ts, i, j, ops)
      } else{
        if(ops == "db1"){
          res.bc[ind, ] = db.loocv(ts, valdb(i, psi.f=0, db_number=1,len.n=10000), j)
        }else{
          library(stringr)
          filename = paste(ops,"_fa_table", sep = "")
          aux_str = str_split(ops, "")[[1]]
          if(aux_str[1] == "d"){
            if(length(aux_str) == 3){
              c = as.numeric(aux_str[3])
            }else
            {
              c = as.numeric(paste(aux_str[3], aux_str[length(aux_str)], sep = ""))
            }

          } else{
            c = as.numeric(aux_str[3])*3
          }
          library(RCurl)
          x <- getURL(paste("https://raw.githubusercontent.com/xcding1212/Sie2nts/main/db_table/", filename, sep = ""))
          dbt = read.csv(text = x)
          dbt = dbt[,2]
          res.bc[ind, ] = db.loocv(ts, valdb(i, dbt, db_number = c, len.n = 0), j)
          }
      }

      ind = ind + 1
    }
  }
  colnames(res.bc) = c("c", "b", "cv")
  if(method == "Elbow"){

    b.s = res.bc[which(res.bc[,3] == min(res.bc[, 3])),2]
    res.bc = res.bc[which(res.bc[,2] == b.s), ]

    if(threshold == 0){
      c.s = 1 + which(abs(res.bc[1:(length(res.bc[,3])-1),3]/res.bc[-1,3] - 1) == max(abs(res.bc[1:(length(res.bc[,3])-1),3]/res.bc[-1,3] - 1)))
    } else{
      c.s = max(which(abs(res.bc[1:(length(res.bc[,3])-1),3]/res.bc[-1,3] - 1) >= threshold)) + 1
    }
    estimate = fix.fit.wavelet(ts, c.s, b.s, m, ops)

  }else{
    b.s = res.bc[which(res.bc[,3] == min(res.bc[, 3])),2]
    c.s = res.bc[which(res.bc[,3] == min(res.bc[, 3])),1]
    estimate = fix.fit.wavelet(ts, c.s, b.s, m, ops)
  }

  return(list(Estimate = estimate[[2]], CV = res.bc, Coefficients = estimate[[1]], BC = c(c.s, b.s)))
}


db_basis.f = function(b.table, x){
  return(as.numeric(b.table[which(abs(b.table[,1]-x) == min(abs(b.table[,1]-x)))[1], -1]))
}

# Testing
mv_method.wav = function(timese, db.basis, k, b, ops){ # c means c+2 basis line
  library(splines)
  h.0 = 3
  m.li = c(1:25)
  library(Matrix)
  # Design matrix

  Y = simu_db(timese, db.basis, b, 10000)[[3]]

  n = length(timese)
  # li.res = list()
  # m = 6

  # Error, i = b* + 1... n
  error.s = c()
  es.alpha = fix.fit.wavelet(timese, k, b, n, ops)[[2]]
  aux.len = length(es.alpha)
  for(i in (b+1):n){
    val.aux = es.alpha[[1]][i]
    for(j in 2:aux.len){
      val.aux = val.aux + es.alpha[[j]][i]*timese[i-j+1]
    }
    error.s[i-b] = timese[i] - val.aux
  }

  Phi.li = list()
  for (m in m.li){
    aux_Phi=0
    Phi = 0
    for(i in (b+1):(n-m)){
      h = 0
      for(j in i:(i+m)){
        aux.h = matrix(rev(c(timese[(j- b):(j - 1)],1)), ncol = 1)*error.s[j-b]
        h = h + aux.h
      }
      B = matrix(db_basis.f(db.basis, i/n), ncol = 1)
      Phi = Phi + kronecker(h, B)
      aux_Phi = aux_Phi + Phi%*%t(Phi)
    }
    Phi.li[[m]] = 1/((n-m-b+1)*m)*aux_Phi
  }

  se.li = list()
  for(mj in (min(m.li)+h.0):(max(m.li)-h.0)){
    av.Phi = 0
    se = 0
    for (k in -3:3){
      av.Phi = av.Phi + Phi.li[[mj + k]]
    }
    av.Phi = av.Phi/7

    for(k in -3:3){
      se = se + norm(av.Phi - Phi.li[[mj + k]], "2")^2
    }
    se.li[[mj-3]] = sqrt(se/6)
  }
  return(m.op = which(unlist(se.li) == min(unlist(se.li))) + 3)
  # return(unlist(se.li))
}


fix.test.wavelet = function(timese, k, b, ops, B.s, m){
  library(Matrix)
  # Design matrix
  if(ops == "db1"){
    db.basis = valdb(k, psi.f=0, db_number=1,len.n=10000)
  }else{
    library(stringr)
    filename = paste(ops,"_fa_table", sep = "")
    aux_str = str_split(ops, "")[[1]]
    if(aux_str[1] == "d"){
      if(length(aux_str) == 3){
        c = as.numeric(aux_str[3])
      }else
      {
        c = as.numeric(paste(aux_str[3], aux_str[length(aux_str)], sep = ""))
      }

    } else{
      c = as.numeric(aux_str[3])*3
    }
    library(RCurl)
    x <- getURL(paste("https://raw.githubusercontent.com/xcding1212/Sie2nts/main/db_table/", filename, sep = ""))
    dbt = read.csv(text = x)
    dbt = dbt[,2]
    db.basis = valdb(k, dbt, db_number = c, len.n = 0)
  }

  c = dim(db.basis)[2] - 1
  Y = simu_db(timese, db.basis, b, 10000)[[3]]
  n = length(timese)
  # li.res = list()
  if(m == 0){
    m = mv_method.wav(timese, db.basis, k, b, ops)  #floor(n^(1/3))
  }
  # m = 6

  esti = simu_db(timese, db.basis, b, 10000)[[1]] # the estimate of coefficients

  # Error, i = b* + 1... n
  error.s = c()
  es.alpha = simu_db(timese, db.basis, b, n)[[1]]
  aux.len = length(es.alpha)
  for(i in (b+1):n){
    val.aux = es.alpha[[1]][i]
    for(j in 2:aux.len){
      val.aux = val.aux + es.alpha[[j]][i]*timese[i-j+1]
    }
    error.s[i-b] = timese[i] - val.aux
  }

  # B
  inte = db_basis.f(db.basis, 1/10000)*(1/10000)
  for(i in 2:10000){
    inte = inte + db_basis.f(db.basis, i/10000)*(1/10000)
  }
  r.c = c # 2*(c-1)+1 only for tri basis.



  # I.bc
  I = matrix(rep(0, ((b+1)*r.c)^2), ncol = (b+1)*r.c)
  for(ind in 0:(b*r.c-1)){
    I[dim(I)[1] - ind, dim(I)[2] -ind] = 1
  }

  nT = 0
  for (i in 2:length(esti)){
    nT = nT + sum(((esti[[i]] - sum(esti[[i]]/10000))^2)/10000)
  }
  nT = n*nT

  Sigma = n*solve(t(Y)%*%Y)
  inte = matrix(inte, ncol =1)
  W = diag(r.c) - inte%*%t(inte)
  W = matrix(bdiag(replicate(b+1,W,simplify=FALSE)), ncol = (b+1)*r.c)

  Tao = Sigma%*%I%*%W%*%Sigma

  # hist(unlist(Sta))
  # print(ite)

  Sta = list()
  Phi.li = list()

  for(k in 1:B.s){
    R = rnorm(n-m-b, 0, 1)
    Phi = 0
    for(i in (b+1):(n-m)){
      h = 0
      for(j in i:(i+m)){
        aux.h = matrix(rev(c(timese[(j- b):(j - 1)],1)), ncol = 1)*error.s[j-b]
        h = h + aux.h
      }
      B = matrix(db_basis.f(db.basis, i/n), ncol = 1)
      Phi = Phi + kronecker(h, B)*R[i-b]
    }

    Phi = (1/sqrt((n-m-b+1)*m))*Phi
    Phi.li[[k]] = Phi
  }
  # image(W)
  # W[(c+1):dim(W)[1], (c+1):dim(W)[2]] = 0
  for(k in 1:B.s){
    Sta[[k]] = t(Phi.li[[k]])%*%Tao%*%Phi.li[[k]]
  }

  # nT > sort(unlist(Sta))[950] if TRUE reject the null
  return(1 - sum(unlist(Sta) <= nT)/B.s)

}



# testing b function(timese, k, b, ops)
fit.testing.b.wavelet = function(timese, k, b.0 = 3,  ops, b = 8, B.s, m){
  if(b.0 >= b){return(FALSE)}
  library(Matrix)
  if(ops == "db1"){
    db.basis = valdb(k, psi.f=0, db_number=1,len.n=10000)
  }else{
    library(stringr)
    filename = paste(ops,"_fa_table", sep = "")
    aux_str = str_split(ops, "")[[1]]
    if(aux_str[1] == "d"){
      if(length(aux_str) == 3){
        c = as.numeric(aux_str[3])
      }else
      {
        c = as.numeric(paste(aux_str[3], aux_str[length(aux_str)], sep = ""))
      }

    }  else{
      c = as.numeric(aux_str[3])*3
    }
    library(RCurl)
    x <- getURL(paste("https://raw.githubusercontent.com/xcding1212/Sie2nts/main/db_table/", filename, sep = ""))
    dbt = read.csv(text = x)
    dbt = dbt[,2]
    db.basis = valdb(k, dbt, db_number = c, len.n = 0)
  }
  c = dim(db.basis)[2] - 1
  Y = simu_db(timese, db.basis, b, 10000)[[3]]
  n = length(timese)
  # li.res = list()
  if(m == 0){
    m = mv_method.wav(timese, db.basis, k, b, ops)  #floor(n^(1/3))
  }
  # m = 6

  esti = simu_db(timese, db.basis, b, 10000)[[1]] # the estimate of coefficients

  # Error, i = b* + 1... n
  error.s = c()
  es.alpha = simu_db(timese, db.basis, b, n)[[1]]
  aux.len = length(es.alpha)
  for(i in (b+1):n){
    val.aux = es.alpha[[1]][i]
    for(j in 2:aux.len){
      val.aux = val.aux + es.alpha[[j]][i]*timese[i-j+1]
    }
    error.s[i-b] = timese[i] - val.aux
  }

  r.c = c # 2*(c-1)+1 only for tri basis.
  aux.pval = list()


  # B
  for(k.aux in 0:(b.0-1)){ # 0 ---> 1-15
    # 1 ---> 2-15
    nT = 0
    for (i in (2+k.aux):length(esti)){
      nT = nT + sum((esti[[i]]^2)/10000)
    }
    nT = n*nT

    # I.bc
    I = matrix(rep(0, ((b+1)*r.c)^2), ncol = (b+1)*r.c)
    for(ind in 0:((b-k.aux)*r.c-1)){
      I[dim(I)[1] - ind, dim(I)[2] -ind] = 1
    }

    Sigma = n*solve(t(Y)%*%Y)
    Tao = Sigma%*%I%*%Sigma

    Sta = list()
    Phi.li = list()

    for(k in 1:B.s){
      R = rnorm(n-m-b, 0, 1)
      Phi = 0
      for(i in (b+1):(n-m)){
        h = 0
        for(j in i:(i+m)){
          aux.h = matrix(rev(c(timese[(j- b):(j - 1)],1)), ncol = 1)*error.s[j-b]
          h = h + aux.h
        }
        B = matrix(db_basis.f(db.basis, i/n), ncol = 1)
        Phi = Phi + kronecker(h, B)*R[i-b]
      }

      Phi = (1/sqrt((n-m-b+1)*m))*Phi
      Phi.li[[k]] = Phi
    }
    for(k.aux2 in 1:B.s){
      Sta[[k.aux2]] = t(Phi.li[[k.aux2]])%*%Tao%*%Phi.li[[k.aux2]]

    }
    aux.pval[[k.aux+1]] = 1 - sum(unlist(Sta) <= nT)/B.s
  }
  # nT > sort(unlist(Sta))[950] if TRUE reject the null
  return(aux.pval)
}




