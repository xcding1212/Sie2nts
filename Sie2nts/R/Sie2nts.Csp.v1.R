# Cubic-spline Demo
# Write function for C-spline, # c is the number of the cubic spline.In detail, splines library is called.

## ord options default 4   number of basis c-2+r
# ord = 5 we need add 4 additional points

# c is at least 2, which is the number of the basis function - 2, the true number of basis is c+2

phi_csp = function(cspline, beta, b){
  c = length(cspline)
  b_res = list()
  for(i in 0:b){
    B.aux = matrix(c(rep(0, c*i), cspline, rep(0, c*(b-i))), ncol = 1)
    b_res[[i+1]] = as.numeric(t(beta)%*%B.aux)
  }

  return(b_res)
}

Cspline_kth_b = function(k, n, c, or){
  #library(splines)
  x = seq(0, 1, length=n)
  knots=c(rep(0, or - 1), seq(0,1, length.out = c), rep(1, or - 1)) ##need to add three additional knots at the two ends;
  B=splineDesign(knots, x, ord = or)                ##default order: ord =4, corresponds to cubic splines
  csptable = cbind(as.matrix(x, ncol = 1), B)
  inte = csp_basis.f(csptable, 1/10000)^2*(1/10000)
  for(i in 2:10000){
    inte = inte + csp_basis.f(csptable, i/10000)^2*(1/10000)
    }

  for(i in 2:dim(csptable)[2]){
    csptable[, i] = sqrt(1/inte)[i-1]*csptable[, i]
    }
  B = csptable[, -1]
  df = data.frame(basis_value = B[, k])
  return(df)
}



cspline_plot = function(c, or, title){
  #library(splines)
  n = 1024
  x = seq(0, 1, length=n)
  knots=c(rep(0, or - 1), seq(0,1, length.out = c), rep(1, or - 1)) ##need to add three additional knots at the two ends;
  B=splineDesign(knots, x, ord = or)                ##default order: ord =4, corresponds to cubic splines
  csptable = cbind(as.matrix(x, ncol = 1), B)
  inte = csp_basis.f(csptable, 1/10000)^2*(1/10000)
  for(i in 2:10000){
    inte = inte + csp_basis.f(csptable, i/10000)^2*(1/10000)
    }

  for(i in 2:dim(csptable)[2]){
    csptable[, i] = sqrt(1/inte)[i-1]*csptable[, i]
    }
  B = csptable[, -1]
  df = data.frame()
  for (i in 1:(c-2+or)){
    res = as.data.frame(B[, i]) #unlist(poly_val(coeffi, i))
    df = rbind(df, res)
  }
  f.df = df
  colnames(f.df) = "new"
  f.df$x = rep(x, c-2+or)
  f.df$order = as.factor(rep(1:(c-2+or), each = 1024))
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


alpha.Cspline = function(ts, c, b, or, m=500){
  #library(splines)
  l.alpha = list()
  aux.alpha = c()
  n = length(ts)
  x = seq(0, 1, length=n)
  knots=c(rep(0, or - 1), seq(0,1, length.out = c), rep(1, or - 1)) ##need to add three additional knots at the two ends;
  B=splineDesign(knots, x, ord = or)
  csptable = cbind(as.matrix(x, ncol = 1), B)
  inte = csp_basis.f(csptable, 1/10000)^2*(1/10000)
  for(i in 2:10000){
    inte = inte + csp_basis.f(csptable, i/10000)^2*(1/10000)
    }

  for(i in 2:dim(csptable)[2]){
    csptable[, i] = sqrt(1/inte)[i-1]*csptable[, i]
    }
  B = csptable[, -1]
  aux_B = B[(b+1):n,]
  ind = b
  aux.Y = matrix(rep(ts[ind:(n-b+ind-1)],dim(aux_B)[2]), ncol = dim(aux_B)[2])
  Y = cbind(aux_B, aux_B*aux.Y)
  ind = ind - 1
  while(ind >= 1){
    aux.Y = matrix(rep(ts[ind:(n-b+ind-1)],dim(aux_B)[2]), ncol = dim(aux_B)[2])
    Y = cbind(Y, aux_B*aux.Y)
    ind = ind - 1
  }
  X = matrix(ts[(b+1):n], ncol = 1)
  beta = solve(t(Y)%*%Y, tol = 1e-40)%*%t(Y)%*%X

  x = seq(0, 1, length=m)
  knots=c(rep(0, or - 1), seq(0,1, length.out = c), rep(1, or - 1)) ##need to add three additional knots at the two ends;
  B=splineDesign(knots, x, ord = or)
  csptable = cbind(as.matrix(x, ncol = 1), B)
  inte = csp_basis.f(csptable, 1/10000)^2*(1/10000)
  for(i in 2:10000){
    inte = inte + csp_basis.f(csptable, i/10000)^2*(1/10000)
  }

  for(i in 2:dim(csptable)[2]){
    csptable[, i] = sqrt(1/inte)[i-1]*csptable[, i]
  }
  B = csptable[, -1]
  for(j in 1:(b+1)){
    for(i in 1:m){
      cspline = B[i,]
      aux.alpha[i] = phi_csp(cspline, beta, b)[[j]]
    }
    l.alpha[[j]] = aux.alpha
    aux.alpha = c()
  }

  return(list(l.alpha, beta, Y))
}


alpha.loocv.mc = function(ts, c, b, or){
  #library(splines)
  n = length(ts)
  aux.true = ts[(b+1):n]
  aux.esti = c()
  leve.i = c()

  beta.es = alpha.Cspline(ts, c, b, or=or, n)
  hat = beta.es[[3]]%*%solve(t(beta.es[[3]])%*%beta.es[[3]])%*%t(beta.es[[3]])

  for(i in (b+1):n){
    aux.esti[i-b]  = matrix(beta.es[[3]][i-b,], nrow = 1)%*%beta.es[[2]]
    leve.i[i-b] = as.numeric(hat[i-b,i-b]) # hii is the diagonal of the hat matrix
  }
  error = (aux.true - aux.esti)^2
  lever = sum(error/((1-leve.i)^2))/(n-b)
  return(c(c,b,lever))
}

alpha.cv.mc = function(ts, c, b, or){
  #library(splines)
  n = length(ts)
  l = floor(3*log2(n))
  aux.train = ts[1:(n-l)]
  aux.vali = ts[(n-l+1):n]
  tt = fix.fit.cspline(aux.train, c, b,or=or, length(aux.train))
  pre = predict.cspline(aux.train, tt, length(aux.vali))
  error = sum((aux.vali - pre)^2)/l
  return(c(c,b,error))
}

# prediction
predict.cspline = function(ts, esti.li, k){
  #library(splines) # k indicates the number of predictions
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

fix.fit.cspline = function(ts, c, b, or, m){
  #library(splines)
  error.s = c()
  n = length(ts)
  es.alpha = alpha.Cspline(ts, c, b,or=or, n)[[1]]
  aux.len = length(es.alpha)
  for(i in (b+1):n){
    val.aux = es.alpha[[1]][i]
    for(j in 2:aux.len){
      val.aux = val.aux + es.alpha[[j]][i]*ts[i-j+1]
    }
    error.s[i-b] = ts[i] - val.aux
  }

  return(list(ols.coef = alpha.Cspline(ts, c, b, or=or, m)[[2]], ts.coef = alpha.Cspline(ts, c, b, or=or, m)[[1]], Residuals = error.s))
}

auto.fit.cspline = function(ts, c = 10, b = 3, or, m=500 ,method = "CV", threshold = 0){
  #library(splines)
  res.bc = matrix(ncol = 3, nrow = (c-1)*b)
  ind = 1
  for(i in 2:c){
    for(j in 1:b){
      if(method == "CV"){
        res.bc[ind, ] = alpha.cv.mc(ts, i, j, or = or)
      } else{
        res.bc[ind, ] = alpha.loocv.mc(ts, i, j, or = or)
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
    estimate = alpha.Cspline(ts, c.s, b.s, or = or, m)

  }else{
    b.s = res.bc[which(res.bc[,3] == min(res.bc[, 3])),2]
    c.s = res.bc[which(res.bc[,3] == min(res.bc[, 3])),1]
    estimate = alpha.Cspline(ts, c.s, b.s, or = or, m)
  }

  return(list(Estimate = estimate[[1]], CV = res.bc, Coefficients = estimate[[2]], BC = c(c.s, b.s)))
}



csp_basis.f = function(b.table, x){
  return(as.numeric(b.table[which(abs(b.table[,1]-x) == min(abs(b.table[,1]-x)))[1], -1]))
}

# Testing
mv_method.cspline = function(timese, c, b, or){ # c means c+2 basis line
  #library(splines)
  #library(Matrix)
  h.0 = 3
  m.li = c(1:25)
  # Design matrix
  Y = alpha.Cspline(timese, c, b, or = or, 10000)[[3]]
  n = length(timese)
  # li.res = list()
  # m = 6
  x = seq(0, 1, length=5000)
  knots=c(rep(0, or - 1), seq(0,1, length.out = c), rep(1, or - 1)) ##need to add three additional knots at the two ends;
  csptable=splineDesign(knots, x, ord = or)
  csptable = cbind(as.matrix(x, ncol = 1), csptable)
  inte = csp_basis.f(csptable, 1/10000)^2*(1/10000)
  for(i in 2:10000){
    inte = inte + csp_basis.f(csptable, i/10000)^2*(1/10000)
  }

  for(i in 2:dim(csptable)[2]){
    csptable[, i] = sqrt(1/inte)[i-1]*csptable[, i]
  }

  # Error, i = b* + 1... n
  error.s = c()
  es.alpha = alpha.Cspline(timese, c, b,or =or, n)[[1]]
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
      B = matrix(csp_basis.f(csptable, i/n), ncol = 1)
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

fix.test.cspline = function(timese, c, b, or, B.s, m){
  #library(splines)
  #library(Matrix)
  # Design matrix
  Y = alpha.Cspline(timese, c, b,or = or, 10000)[[3]]
  n = length(timese)
  # li.res = list()
  if(m == 0){
    m = mv_method.cspline(timese, c, b, or = or)  #floor(n^(1/3))
  }
  # m = 6

  x = seq(0, 1, length=5000)
  knots=c(rep(0, or - 1), seq(0,1, length.out = c), rep(1, or - 1)) ##need to add three additional knots at the two ends;
  B=splineDesign(knots, x, ord = or)

  csptable = cbind(as.matrix(x, ncol = 1), B)
  inte = csp_basis.f(csptable, 1/10000)^2*(1/10000)
  for(i in 2:10000){
    inte = inte + csp_basis.f(csptable, i/10000)^2*(1/10000)
    }

  for(i in 2:dim(csptable)[2]){
    csptable[, i] = sqrt(1/inte)[i-1]*csptable[, i]
    }


  esti = alpha.Cspline(timese, c, b, or = or, 10000)[[1]] # the estimate of coefficients

  # Error, i = b* + 1... n
  error.s = c()
  es.alpha = alpha.Cspline(timese, c, b, or = or, n)[[1]]
  aux.len = length(es.alpha)
  for(i in (b+1):n){
    val.aux = es.alpha[[1]][i]
    for(j in 2:aux.len){
      val.aux = val.aux + es.alpha[[j]][i]*timese[i-j+1]
    }
    error.s[i-b] = timese[i] - val.aux
  }

  # B
  inte = csp_basis.f(csptable, 1/10000)*(1/10000)
  for(i in 2:10000){
    inte = inte + csp_basis.f(csptable, i/10000)*(1/10000)
  }
  r.c = c - 2 + or  # 2*(c-1)+1 only for tri basis. c-2+or only for cspline

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
      B = matrix(csp_basis.f(csptable, i/n), ncol = 1)
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



fit.testing.b.cspline = function(timese, c, or, b.0 = 3, b = 8, B.s, m){
  #library(splines)
  if(b.0 >= b){return(FALSE)}
  #library(Matrix)
  # Design matrix
  Y = alpha.Cspline(timese, c, b,or=or, 10000)[[3]]   # can not fit very well with b = 8
  n = length(timese)
  # li.res = list()
  if(m == 0){
    m = mv_method.cspline(timese, c, b, or = or) # does m influenced by b ??? mv_method(timese, c, b)  #floor(n^(1/3))
  }

  esti = alpha.Cspline(timese, c, b, or=or, 10000)[[1]]  # the estimate of coefficients

  x = seq(0, 1, length= 5000)
  knots=c(rep(0, or - 1), seq(0,1, length.out = c), rep(1, or - 1)) ##need to add three additional knots at the two ends;
  B=splineDesign(knots, x, ord = or)

  csptable = cbind(as.matrix(x, ncol = 1), B)
  inte = csp_basis.f(csptable, 1/10000)^2*(1/10000)
  for(i in 2:10000){
    inte = inte + csp_basis.f(csptable, i/10000)^2*(1/10000)
    }

  for(i in 2:dim(csptable)[2]){
    csptable[, i] = sqrt(1/inte)[i-1]*csptable[, i]
    }

  # Error, i = b* + 1... n
  error.s = c()
  es.alpha = alpha.Cspline(timese, c, b,or=or, n)[[1]]
  aux.len = length(es.alpha)
  for(i in (b+1):n){
    val.aux = es.alpha[[1]][i]
    for(j in 2:aux.len){
      val.aux = val.aux + es.alpha[[j]][i]*timese[i-j+1]
    }
    error.s[i-b] = timese[i] - val.aux
  }

  r.c = c - 2 + or # 2*(c-1)+1 only for tri basis.
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
        B = matrix(csp_basis.f(csptable, i/n), ncol = 1)
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





