# c indicate the number of tri basis function. The true number of basis functions are 2*c - 1. For the rest basis, c expresses the number of basis function.
# Compared with version1, F_general() is removed

F_tri = function(c, x){
  aux = c()
  if(c == 1){
    aux[c] = 1
    return(aux)
  }
  ind = 1
  aux[ind] = 1
  ind = ind + 1
  for(i in 1:(c-1)){
    aux[ind] = sqrt(2)*sin(2*i*pi*x)
    aux[ind+1] = sqrt(2)*cos(2*i*pi*x)
    ind = ind + 2
  }
  return(aux)

}

F_cosPol = function(c, x){
  aux = c()
  if(c == 1){
    aux[c] = 1
    return(aux)
  }
  for(i in 1:c){
    if(i == 1){
      aux[i] = 1
    } else{
      aux[i] = sqrt(2)*cos((i-1)*pi*x)
    }
  }
  return(aux)

}

F_sinPol = function(c, x){
  aux = c()
  for(i in 1:c){
    aux[i] = sqrt(2)*sin(i*pi*x)

  }
  return(aux)
}

# F_general = function(c, x){
#   aux = c()
#   for(i in 1:c){
#     if (i == 1){
#       aux[i] = 1
#     } else if (i %% 2 == 0){
#       aux[i] = sqrt(2)*sin(i*pi*x)
#     } else{
#       aux[i] = sqrt(2)*cos((i-1)*pi*x)
#     }
#   }
#   return(aux)
# }

# If the option parameter equals tri, it means we choose trigometric basis, cos means cospol, sin means sinpol.
select.basis = function(c,x, ops = "tri"){
  if(ops == "tri"){
    return(F_tri(c, x))
  } else if (ops == "cos"){
    return(F_cosPol(c,x))
  } else{
    return(F_sinPol(c,x))
  }
}


Fourier_kth_b = function(k, n, ops){
  df = data.frame()
  aux_x = seq(0,1, length.out = n)
  # coeffi = legendre_coeff(k)[k]
  for (i in aux_x){
    res = select.basis(k, i, ops)
    df = rbind(df, res)
  }
  df = data.frame(basis_value = df[,k])
  return(df)
}


fourier_plot = function(c, ops = "tri", title){

  df = data.frame()
  aux_x = seq(0,1,0.005)

  for (i in aux_x){
    res = select.basis(c, i, ops)
    df = rbind(df, res)
  }

  new = c()
  for(i in 1:dim(df)[2]){
    new = c(new, df[, i])
  }
  f.df = as.data.frame(new)
  f.df$x = rep(aux_x, dim(df)[2])
  f.df$order = as.factor(rep(0:(dim(df)[2]-1), each = 201))
  theme_update(plot.title = element_text(hjust = 0.5))
  p1 <- ggplot(f.df, aes(x=x, y=new, group=order, colour = order))+ geom_line()  + ggtitle(title) +
    xlab("") + ylab("") + scale_colour_discrete(name  ="order")+theme(plot.title = element_text(size=18, face="bold"),
                                                                      legend.text=element_text(size=24, face = "bold"),
                                                                      axis.text.x = element_text(face="bold", color="#993333",                                                                                                                              size=22, angle=0),
                                                                      axis.text.y = element_text(face="bold", color="#993333",size=22, angle=0),
                                                                      axis.title.x=element_text(size=22,face='bold'),
                                                                      axis.title.y=element_text(angle=90, face='bold', size=22),
                                                                      legend.title = element_text(face = "bold"))
  return(p1)
}

# fourier_plot(2)

# In the function beta_f(), alpha.fou(), basis depends on the options, we have 4 choice for fourier basis.
# beta_f returns the list which contains beta and design matrix Y

beta_f = function(ts, c, b, ops = "tri"){
  n = length(ts)
  X = matrix(ts[(b+1):n], ncol = 1)
  aux = c(select.basis(c, (b+1)/n, ops))
  for(j in 1:b){
    aux = c(aux, select.basis(c, ((b+1)/n), ops)*ts[b+1-j])
  }
  Y = matrix(aux, nrow = 1) # i =2 j >= 1

  for(i in (b+2):n){
    aux = c(select.basis(c, i/n, ops))
    for(j in 1:b){
      aux = c(aux, select.basis(c, (i/n),ops)*ts[i-j])
    }

    aux_Y = matrix(aux, nrow = 1)
    Y = rbind(Y, aux_Y)
  }
  beta = solve(t(Y)%*%Y, tol = 1e-40)%*%t(Y)%*%X
  return(list(beta, Y))
}

phi_f = function(fourier, beta, b){
  c = length(fourier)
  b_res = list()
  for(i in 0:b){
    B.aux = matrix(c(rep(0, c*i), fourier, rep(0, c*(b-i))), ncol = 1)
    b_res[[i+1]] = as.numeric(t(beta)%*%B.aux)
  }
  return(b_res)

}


# we want to generate m points of coefficients of time series and the default number is 500.
alpha.fou = function(ts, c, b, m=500, ops){

  l.alpha = list()
  aux.alpha = c()
  beta.es = beta_f(ts, c, b, ops)

  for(j in 1:(b+1)){
    for(i in 1:m){
      aux.alpha[i] = phi_f(select.basis(c, i/m, ops), beta.es[[1]], b)[[j]]
    }
    l.alpha[[j]] = aux.alpha
    aux.alpha = c()
  }
  return(l.alpha)
}

alpha.loocv.f = function(ts, c, b, ops){
  n = length(ts)
  aux.true = ts[(b+1):n]
  aux.esti = c()
  leve.i = c()

  beta.es = beta_f(ts, c, b, ops)
  hat = beta.es[[2]]%*%solve(t(beta.es[[2]])%*%beta.es[[2]])%*%t(beta.es[[2]])

  for(i in (b+1):n){
    aux.esti[i-b]  = matrix(beta.es[[2]][i-b,], nrow = 1)%*%beta.es[[1]]
    leve.i[i-b] = as.numeric(hat[i-b,i-b]) # hii is the diagonal of the hat matrix
  }
  error = (aux.true - aux.esti)^2
  lever = sum(error/((1-leve.i)^2))/(n-b)
  return(c(c,b,lever))
}

#  CV in the paper
alpha.cv.f = function(ts, c, b, ops){
  n = length(ts)
  l = floor(3*log2(n))
  aux.train = ts[1:(n-l)]
  aux.vali = ts[(n-l+1):n]
  tt = fix.fit.four(aux.train, c, b, length(aux.train), ops)
  pre = predict.Four(aux.train, tt, length(aux.vali))
  error = sum((aux.vali - pre)^2)/l
  return(c(c,b,error))
}


# prediction
predict.Four = function(ts, esti.li, k){ # k indicates the number of predictions
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


# The return of fit.ts.f() is the list contains 4 parts named Estimate, cv, coefficients and bc. Rather, Estimate is the estimate of coefficient for the time series
# cv is the cross validation matrix, Coefficients is the estimate of coefficient for each basis function, BC contains the best number of b and c basis on LOOCV method.
# fit.ts.f() automatically choose the best b and c for the time series and get the estimate basis on that.


fix.fit.four = function(ts, c, b, m, ops){

  error.s = c()
  n = length(ts)
  es.alpha = alpha.fou(ts, c, b, n, ops)
  aux.len = length(es.alpha)
  for(i in (b+1):n){
    val.aux = es.alpha[[1]][i]
    for(j in 2:aux.len){
      val.aux = val.aux + es.alpha[[j]][i]*ts[i-j+1]
    }
    error.s[i-b] = ts[i] - val.aux
  }

  return(list(ols.coef = beta_f(ts, c, b, ops)[[1]], ts.coef = alpha.fou(ts, c, b, m, ops), Residuals = error.s))

}

auto.fit.four = function(ts, c = 10, b = 3, m = 500, ops, method = "LOOCV", threshold = 0){
  res.bc = matrix(ncol = 3, nrow = c*b)
  ind = 1
  for(i in 1:c){
    for(j in 1:b){
      if(method == "CV"){
        res.bc[ind, ] = alpha.cv.f(ts, i, j, ops)
      } else{
        res.bc[ind, ] = alpha.loocv.f(ts, i, j, ops)
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
    estimate = alpha.fou(ts, c.s, b.s, m, ops)

  }else{
    b.s = res.bc[which(res.bc[,3] == min(res.bc[, 3])),2]
    c.s = res.bc[which(res.bc[,3] == min(res.bc[, 3])),1]
    estimate = alpha.fou(ts, c.s, b.s, m, ops)
  }


  return(list(Estimate = estimate, CV = res.bc, Coefficients = beta_f(ts, c.s, b.s)[[1]], BC = c(c.s, b.s)))
}

# Testing
mv_method.four = function(timese, c, b, ops){
  h.0 = 3
  m.li = c(1:25)
  library(Matrix)
  # Design matrix
  Y = beta_l(timese, c, b)[[2]]
  n = length(timese)
  # li.res = list()
  # m = 6

  # Error, i = b* + 1... n
  error.s = c()
  es.alpha = alpha.fou(timese, c, b, n, ops)
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
      B = matrix(select.basis(c, i/n, ops), ncol = 1)
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


fix.test.four = function(timese, c, b, ops, B.s, m){
  library(Matrix)
  # Design matrix
  Y = beta_f(timese, c, b, ops)[[2]]
  n = length(timese)
  # li.res = list()
  # m = 6
  if(m == 0){
    m = mv_method.four(timese, c, b, ops) #mv_method(timese, c, b)  #floor(n^(1/3))
  }

  esti = alpha.fou(timese, c, b, 10000, ops) # the estimate of coefficients

  # Error, i = b* + 1... n
  error.s = c()
  es.alpha = alpha.fou(timese, c, b, n, ops)
  aux.len = length(es.alpha)
  for(i in (b+1):n){
    val.aux = es.alpha[[1]][i]
    for(j in 2:aux.len){
      val.aux = val.aux + es.alpha[[j]][i]*timese[i-j+1]
    }
    error.s[i-b] = timese[i] - val.aux
  }
  length(error.s)


  # B
  inte = select.basis(c, 1/10000, ops)*(1/10000)
  for(i in 2:10000){
    inte = inte + select.basis(c, i/10000, ops)*(1/10000)
  }
  if(ops == "tri"){
    r.c = 2*(c-1)+1

  }else{
    r.c = c # 2*(c-1)+1 only for tri basis.
  }

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

  Sigma = n*solve(t(Y)%*%Y, tol = 1e-40)
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
      B = matrix(select.basis(c, i/n, ops), ncol = 1)
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
  return(1 - sum(unlist(Sta) <= nT)/B.s) # P value

}



# testing b


fit.testing.b.four = function(timese, c, b.0 = 3, ops, b = 8, B.s, m){
  if(b.0 >= b){return(FALSE)}
  library(Matrix)
  # Design matrix
  Y = beta_f(timese, c, b, ops)[[2]]
  n = length(timese)
  # li.res = list()
  if(m == 0){
    m = mv_method.four(timese, c, b, ops) # does m influenced by b ??? mv_method(timese, c, b)  #floor(n^(1/3))
  }

  esti = alpha.fou(timese, c, b, 10000, ops) # the estimate of coefficients
  # Error, i = b* + 1... n
  error.s = c()
  es.alpha = alpha.fou(timese, c, b, n, ops)
  aux.len = length(es.alpha)
  for(i in (b+1):n){
    val.aux = es.alpha[[1]][i]
    for(j in 2:aux.len){
      val.aux = val.aux + es.alpha[[j]][i]*timese[i-j+1]
    }
    error.s[i-b] = timese[i] - val.aux
  }
  if(ops == "tri"){
    r.c = 2*(c-1)+1

  }else{
    r.c = c # 2*(c-1)+1 only for tri basis.
  }

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

    Sigma = n*solve(t(Y)%*%Y, tol = 1e-40)
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
        B = matrix(select.basis(c, i/n, ops), ncol = 1)
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




