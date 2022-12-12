# Legendre basis function
# input is the coefficients of polynomial order by (1, x, x^2, x^3,...) until the highest order. The out put is this polynomial multiple by x.
move_order = function(poly){ # get the coefficient when polynomial multiple 1 times of x
  return(c(0,poly))
}

# we can get the plot from the coefficients
# input is n (the number of basis function)
# output is the coefficients of polynomial order by (1, x, x^2, x^3,...) until the highest order.
# always normalized

legendre_coeff = function(n){
  p_c = list()
  p_c[1] = c(1)
  if(n == 1){
    return(p_c[[1]])
  } else if(n==2){
    p_c[[2]] = c(0,sqrt(3))
    return(p_c)

  } else{
    for(i in 3:n){
      aux.i = i - 1
      p_c[[2]] = c(0, 1)
      p_c[[i]] = ((2*aux.i-1)/aux.i)*move_order(p_c[[i-1]]) - ((aux.i-1)/aux.i)*c(p_c[[i-2]], 0, 0)
    }
    p_n = list(p_c[[1]])
    for (i in 2:length(p_c)){
      p_n[[i]] = sqrt((2*i-1))*p_c[[i]]
    }
    # p_n[[1]] = 1/sqrt(2)
    return(p_n)
  }
}


# x \in [0,1]
# input is n (the number of basis function). coeffi is the list and contain the coefficients of polynomial order by (1, x, x^2, x^3,...) until the highest order.
poly_val = function(coeffi, x){
  for(i in 1:length(coeffi)){
    aux = 0
    for(j in 1:length(coeffi[[i]])){
      aux = aux + coeffi[[i]][j]*(2*x-1)^(j - 1)
    }
    coeffi[[i]] = aux
  }
  return(coeffi)
}

# c is number of basis function, x is the inputs value.
Legendre_basis = function(c, x){
  aux_li = legendre_coeff(c)
  return(unlist(poly_val(aux_li, x)))
}


Legendre_kth_b = function(k, n){
  df = data.frame()
  aux_x = seq(0,1, length.out = n)
  coeffi = legendre_coeff(k)[k]
  for (i in aux_x){
    res = unlist(poly_val(coeffi, i))
    df = rbind(df, res)
  }
  colnames(df) = c("basis_value")
  return(df)
}


legendre_plot = function(n, title){
  coeffi = legendre_coeff(n)
  df = data.frame()
  aux_x = seq(0,1,0.005)

  for (i in aux_x){
    res = unlist(poly_val(coeffi, i))
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



# ggtitle(TeX("L_1"))

# beta is the fitted estimators and Y is the design matrix. This Y is for cv.

beta_l = function(ts, c, b){
  n = length(ts)
  X = matrix(ts[(b+1):n], ncol = 1)
  aux = c(Legendre_basis(c, (b+1)/n))
  for(j in 1:b){
    aux = c(aux, Legendre_basis(c, (b+1)/n)*ts[b+1-j])
  }
  Y = matrix(aux, nrow = 1) # i =2 j >= 1

  for(i in (b+2):n){
    aux = c(Legendre_basis(c, i/n))
    for(j in 1:b){
      aux = c(aux, Legendre_basis(c, i/n)*ts[i-j])
    }

    aux_Y = matrix(aux, nrow = 1)
    Y = rbind(Y, aux_Y)
  }
  beta = solve(t(Y)%*%Y, tol = 1e-40)%*%t(Y)%*%X
  return(list(beta, Y))
}



# for each i/n, we can get the coefficient estimated.
phi_l = function(legendre, beta, b){
  c = length(legendre)
  b_res = list()
  for(i in 0:b){
    B.aux = matrix(c(rep(0, c*i), legendre, rep(0, c*(b-i))), ncol = 1)
    b_res[[i+1]] = as.numeric(t(beta)%*%B.aux)
  }

  return(b_res)
}


# ts is the time series data, c is the number of basis, b indicates AR(b_i)

alpha.legen = function(ts, c, b, m=500){

  l.alpha = list()
  aux.alpha = c()
  beta.es = beta_l(ts, c, b)

  for(j in 1:(b+1)){
    for(i in 1:m){
      aux.alpha[i] = phi_l(Legendre_basis(c, i/m), beta.es[[1]], b)[[j]]
    }
    l.alpha[[j]] = aux.alpha
    aux.alpha = c()
  }
  return(l.alpha)
}

# Leave one out
alpha.loocv.l = function(ts, c, b){
  n = length(ts)
  aux.true = ts[(b+1):n]
  aux.esti = c()
  leve.i = c()

  beta.es = beta_l(ts, c, b)
  hat = beta.es[[2]]%*%solve(t(beta.es[[2]])%*%beta.es[[2]])%*%t(beta.es[[2]])

  for(i in (b+1):n){
    aux.esti[i-b]  = matrix(beta.es[[2]][i-b,], nrow = 1)%*%beta.es[[1]]
    leve.i[i-b] = as.numeric(hat[i-b,i-b]) # hii is the diagonal of the hat matrix
  }
  error = (aux.true - aux.esti)^2
  lever = sum(error/((1-leve.i)^2))/(n-b)
  return(c(c,b,lever))
}

#  CV in the paper, the b incorrect when the data set is small.

alpha.cv.l = function(ts, c, b){
  n = length(ts)
  l = floor(3*log2(n))
  aux.train = ts[1:(n-l)]
  aux.vali = ts[(n-l+1):n]
  tt = fix.fit.legen(aux.train, c, b, length(aux.train))
  pre = predict.legen(aux.train, tt, length(aux.vali))
  error = sum((aux.vali - pre)^2)/l
  return(c(c,b,error))
}



# prediction

predict.legen = function(ts, esti.li, k){ # k indicates the number of predictions
  ts.pre = c()
  if(length(esti.li) == 3){
    phi.h = esti.li[[2]]
  }else{
    phi.h = esti.li[[1]]
  }

  n = length(phi.h[[1]])
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


# the default number of c and b indicate the maximum number, fit.ts function will get the best c* and b* automatically as the returns.
# The return is the list contains estimate coefficients, the cross validation matrix, the estimate coefficients for basis functions, and the best c* and b* in order.
# ops = CV indicates cross validation and method = COOV means leave one out cross validation.(Elbow)

fix.fit.legen = function(ts, c, b, m){

  error.s = c()
  n = length(ts)
  es.alpha = alpha.legen(ts, c, b, n)
  aux.len = length(es.alpha)
  for(i in (b+1):n){
    val.aux = es.alpha[[1]][i]
    for(j in 2:aux.len){
      val.aux = val.aux + es.alpha[[j]][i]*ts[i-j+1]
    }
    error.s[i-b] = ts[i] - val.aux
  }

  return(list(ols.coef = beta_l(ts, c, b)[[1]], ts.coef = alpha.legen(ts, c, b, m), Residuals = error.s))
}
#timese = generate_AR1(786, 20)
#pp = fix.fit.legen(timese, 5, 1, 500)
auto.fit.legen = function(ts, c = 10, b = 3, m=500, method = "CV", threshold = 0){
  res.bc = matrix(ncol = 3, nrow = c*b)
  ind = 1
  for(i in 1:c){
    for(j in 1:b){
      if(method == "CV"){
        res.bc[ind, ] = alpha.cv.l(ts, i, j)
      } else{
        res.bc[ind, ] = alpha.loocv.l(ts, i, j)
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
    estimate = alpha.legen(ts, c.s, b.s, m)

  }else{
    b.s = res.bc[which(res.bc[,3] == min(res.bc[, 3])),2]
    c.s = res.bc[which(res.bc[,3] == min(res.bc[, 3])),1]
    estimate = alpha.legen(ts, c.s, b.s, m)
  }

  return(list(Estimate = estimate, CV = res.bc, Coefficients = beta_l(ts, c.s, b.s)[[1]], BC = c(c.s, b.s)))
}

# Testing
mv_method.legen = function(timese, c, b){
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
  es.alpha = alpha.legen(timese, c, b, n)
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
      B = matrix(Legendre_basis(c, i/n), ncol = 1)
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


fix.test.legen = function(timese, c, b, B.s, m){
  library(Matrix)
  # Design matrix
  Y = beta_l(timese, c, b)[[2]]
  n = length(timese)
  # li.res = list()
  if(m == 0){
    m = mv_method.legen(timese, c, b) #mv_method(timese, c, b)  #floor(n^(1/3))
  }

  # m = 6
  esti = alpha.legen(timese, c, b, 10000) # the estimate of coefficients

  # Error, i = b* + 1... n
  error.s = c()
  es.alpha = alpha.legen(timese, c, b, n)
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
  inte = Legendre_basis(c, 1/10000)*(1/10000)
  for(i in 2:10000){
    inte = inte + Legendre_basis(c, i/10000)*(1/10000)
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
      B = matrix(Legendre_basis(c, i/n), ncol = 1)
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




# Testing b

fit.testing.b.legen = function(timese, c, b.0 = 3, b = 8, B.s, m){
  # lag = 7, b = 10, we calculate, 1-10(reject), 2-10(not reject, statistics is small), 3-10 until 7-10 (not reject, statistics is small),
  if(b.0 >= b){return(FALSE)}
  library(Matrix)
  # Design matrix
  Y = beta_l(timese, c, b)[[2]]
  n = length(timese)
  # li.res = list()
  if(m == 0){
    m = mv_method.legen(timese, c, b) #  mv_method(timese, c, b)  #floor(n^(1/3))
  }

  esti = alpha.legen(timese, c, b, 10000) # the estimate of coefficients

  # Error, i = b* + 1... n
  error.s = c()
  es.alpha = alpha.legen(timese, c, b, n)
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

  # testing checking.
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
        B = matrix(Legendre_basis(c, i/n), ncol = 1)
        Phi = Phi + kronecker(h, B)*R[i-b]
      }

      Phi = (1/sqrt((n-m-b+1)*m))*Phi
      Phi.li[[k]] = Phi
    }
    for(k.aux2 in 1:B.s){
      Sta[[k.aux2]] = t(Phi.li[[k.aux2]])%*%Tao%*%Phi.li[[k.aux2]]

    }
    #print(sort(unlist(Sta)))
    #print(nT)
    aux.pval[[k.aux+1]] = 1 - sum(unlist(Sta) <= nT)/B.s
  }

  # nT > sort(unlist(Sta))[950] if TRUE reject the null

  return(aux.pval)
}

