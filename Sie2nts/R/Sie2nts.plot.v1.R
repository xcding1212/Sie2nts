
# This file contains several functions used in plot

fit.plot.estimate = function(alpha.e, c, b, basis, title){
  library(ggplot2)
  n = length(alpha.e[[2]][[1]])
  res = list()
  df = data.frame(phi = alpha.e[[2]][[1]])

  for(i in 2:(b+1)){
    aux_alpha.g = as.data.frame(alpha.e[[2]][[i]])
    colnames(aux_alpha.g) = "phi"
    df = rbind(df, aux_alpha.g)
  }

  df$t = rep(seq(0, 1, length.out = n), b+1)
  df$Phi = as.factor(rep(0:(length(alpha.e[[2]])-1), each = n))
  theme_update(plot.title = element_text(hjust = 0.5))

  res = ggplot(df, aes(x=t, y=phi, group=Phi, color=Phi))+ geom_line()  + ggtitle(title) +
                xlab("t") + ylab("phi") + scale_colour_discrete(name  ="phi")+theme(plot.title = element_text(size=18, face="bold"),
                                                                      legend.text=element_text(size=24, face = "bold"),
                                                                      axis.text.x = element_text(face="bold", color="#993333",size=22, angle=0),
                                                                      axis.text.y = element_text(face="bold", color="#993333",size=22, angle=0),
                                                                      axis.title.x=element_text(size=22,face='bold'),
                                                                      axis.title.y=element_text(angle=90, face='bold', size=22),
                                                                      legend.title = element_text(face = "bold"))


  return(res)
}


fit.plot.elbow = function(alpha.e, c, b, basis, title){
  library(ggplot2)
  df.cv = data.frame(alpha.e$CV)
  theme_update(plot.title = element_text(hjust = 0.5))
  res = ggplot(df.cv[df.cv$b == alpha.e$BC[2], ], aes(x=c, y=cv)) + geom_point() +geom_line(color='darkblue') + scale_x_continuous(limits=c(1,10), breaks=seq(1, 10, 1)) + ggtitle(title) +
                xlab("c") + ylab("cv") + scale_colour_discrete(name  ="phi")+theme(plot.title = element_text(size=18, face="bold"),
                                                                        legend.text=element_text(size=24, face = "bold"),
                                                                        axis.text.x = element_text(face="bold", color="#993333",size=22, angle=0),
                                                                        axis.text.y = element_text(face="bold", color="#993333",size=22, angle=0),
                                                                        axis.title.x=element_text(size=22,face='bold'),
                                                                        axis.title.y=element_text(angle=90, face='bold', size=22),
                                                                        legend.title = element_text(face = "bold"))
  return(res)
}


fit.plot.cvm = function(alpha.e, basis){
  library(plotly)
  cv_m = alpha.e$CV
  df.cv = data.frame(c = cv_m[,1], b = as.factor(cv_m[,2]), cv = cv_m[,3])
  fig <- plot_ly(df.cv, x = ~c, y = ~b, z = ~cv,
                 marker = list(color = ~cv, colorscale = 'Viridis', showscale = TRUE),
                 text = ~paste('c:', c, '<br>b:', b, '<br>cv:', cv))
  fig <- fig %>% add_markers()
  fig <- fig %>% layout(title = basis, scene = list(camera = list(eye = list(x = -1.68, y = 1.68, z = 1.3)), xaxis = list(title = 'c'),
                                                    yaxis = list(title = 'b', tickvals = list(1,2)),
                                                    zaxis = list(title = 'cv')),
                        annotations = list(
                          x = 1.13,
                          y = 1.05,
                          text = 'cv',
                          xref = 'paper',
                          yref = 'paper',
                          showarrow = FALSE
                        ))
  return(fig)

}





fit.plot.estimate.aux = function(alpha.e, c, b, basis, title){
  library(ggplot2)
  n = length(alpha.e[[1]])

  res = list()
  df = data.frame(phi = alpha.e[[1]])

  for(i in 2:(b+1)){
    aux_alpha.g = as.data.frame(alpha.e[[i]])
    colnames(aux_alpha.g) = "phi"
    df = rbind(df, aux_alpha.g)
  }
  df$t = rep(seq(0, 1, length.out = n), b+1)
  df$Phi = as.factor(rep(0:b, each = n))
  theme_update(plot.title = element_text(hjust = 0.5))
  res = ggplot(df, aes(x=t, y=phi, group=Phi, color=Phi)) + geom_line() + ggtitle(title) +
              xlab("t") + ylab("phi") + scale_colour_discrete(name  ="phi")+theme(plot.title = element_text(size=18, face="bold"),
                                                                        legend.text=element_text(size=24, face = "bold"),
                                                                        axis.text.x = element_text(face="bold", color="#993333",size=22, angle=0),
                                                                        axis.text.y = element_text(face="bold", color="#993333",size=22, angle=0),
                                                                        axis.title.x=element_text(size=22,face='bold'),
                                                                        axis.title.y=element_text(angle=90, face='bold', size=22),
                                                                        legend.title = element_text(face = "bold"))

  return(res)
}



pacf_3dplot = function(ff){
  fig <- plot_ly(ff, x = ~t, y = ~class, z = ~pacf, type="mesh3d",intensity= ~pacf)

  fig <- fig %>% layout(title = "legen", scene = list(camera = list(eye = list(x = 2.6, y = 0.15, z = 1)),xaxis = list(title = 't'),
                                                      yaxis = list(title = 'lag'),
                                                      zaxis = list(title = 'PACF')))

  fig
}

