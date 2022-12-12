# plot for pacf

source(paste(getwd(), "/R/Sie2nts.Legen.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.Cheby.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.Four.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.Csp.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.db1-20.v1.R", sep = ""))
source(paste(getwd(), "/R/sie.auto.pacf.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.plot.v1.R", sep = ""))

#' Plot Pacf
#' @description sie.plot.pacf() shows the PACF with different lag.
#'
#' @param ts ts is the data set which is a time series data typically
#' @param c c indicates the number of basis used to estimate (For wavelet, the number of basis is 2^c. If
#'     Cspli is chosen, the real number of basis is c-2+or)
#' @param lag lag b is the lag for auto-regressive model
#' @param type type indicates which type of basis is used (There are 32 types in this package)
#' @param m m indicates the number of points of coefficients to estimate
#' @param ops choose 2D plot or 3D plot ("2d" inicates 2D plot and "3d" indicates 3D plot)
#' @param title give the title for the pacf plot
#' @param or or indicates the order of spline, default is 4 which indicates cubic spline
#'
#' @return The plot of pacf basis on the time series data
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(137)
#' time.series = c()
#' n = 1024
#' v = 25
#' w = rnorm(n, 0, 1) / v
#' x_ini = runif(1,0,1)

#' for(i in 1:n){
#'   if(i == 1){
#'     time.series[i] = 0.2 + 0.6*cos(2*pi*(i/n))*x_ini  + w[i] #
#'   } else{
#'     time.series[i] = 0.2 + 0.6*cos(2*pi*(i/n))*time.series[i-1] + w[i]
#'   }
#' }
#'
#' sie.plot.pacf(time.series, 5, 10, "Legen")
#' sie.plot.pacf(time.series, 5, 10, "Legen", "3d")
#'}


sie.plot.pacf = function(ts, c, lag, type, ops = "2d", title = "", m=500, or =4){
  library(ggplot2)
  aux = list()
  wavelet_basis = c("db1", "db2", "db3", "db4", "db5",
                    "db6", "db7", "db8", "db9", "db10",
                    "db11", "db12", "db13", "db14", "db15",
                    "db16", "db17", "db18", "db19", "db20",
                    "cf1", "cf2", "cf3", "cf4", "cf5"
  )

  if(ops == "2d"){
    if(type == "Legen"){
      if(lag == 1){
        res = fix.fit.legen(ts, c, 1, m)
        val = list()
        val[[1]] =  res$ts.coef[[2]]

        ff = data.frame(pacf = unlist(val))
        ff$t = rep(seq(0,1,length.out = m), dim(ff)[1]/m)
        ff$class = as.factor(rep(c(1:lag), each = m))

        theme_update(plot.title = element_text(hjust = 0.5))
        aux[[1]] = ggplot(ff, aes(x=t, y=pacf, group=class, colour = class))+ geom_line()  + ylim(-1,1) + ggtitle(title) +
          geom_segment(aes(x=0,xend=1, y=1.96*(length(ts)^(-1/2)), yend=1.96*(length(ts)^(-1/2))), linetype="dashed", color = "black", size=0.4) +
          geom_segment(aes(x=0,xend=1, y=-1.96*(length(ts)^(-1/2)), yend=-1.96*(length(ts)^(-1/2))), linetype="dashed", color = "black", size=0.4) +
          xlab("t") + ylab("pacf") + scale_colour_discrete(name  ="lag")+theme(plot.title = element_text(size=18, face="bold"),
                                                                               legend.text=element_text(size=24, face = "bold"),
                                                                               axis.text.x = element_text(face="bold", color="#993333",                                                                                                                              size=22, angle=0),
                                                                               axis.text.y = element_text(face="bold", color="#993333",size=22, angle=0),
                                                                               axis.title.x=element_text(size=22,face='bold'),
                                                                               axis.title.y=element_text(angle=90, face='bold', size=22),
                                                                               legend.title = element_text(face = "bold"))

        aux[[2]] = ff
        return(aux[[1]])

      } else{
        res = fix.fit.legen(ts, c, 1, m)
        val = list()
        val[[1]] =  res$ts.coef[[1+1]]
        for(b in 2:lag){
          res = fix.fit.legen(ts, c, b, m)
          val[[b]] = res$ts.coef[[b+1]]
        }
        ff = data.frame(pacf = unlist(val))
        ff$t = rep(seq(0,1,length.out = m), dim(ff)[1]/m)
        ff$class = as.factor(rep(c(1:lag), each = m))

        theme_update(plot.title = element_text(hjust = 0.5))
        aux[[1]] = ggplot(ff, aes(x=t, y=pacf, group=class, colour = class))+ geom_line()  + ylim(-1,1) + ggtitle(title) +
          geom_segment(aes(x=0,xend=1, y=1.96*(length(ts)^(-1/2)), yend=1.96*(length(ts)^(-1/2))), linetype="dashed", color = "black", size=0.4) +
          geom_segment(aes(x=0,xend=1, y=-1.96*(length(ts)^(-1/2)), yend=-1.96*(length(ts)^(-1/2))), linetype="dashed", color = "black", size=0.4) +
          xlab("t") + ylab("pacf") + scale_colour_discrete(name  ="lag")+theme(plot.title = element_text(size=18, face="bold"),
                                                                               legend.text=element_text(size=24, face = "bold"),
                                                                               axis.text.x = element_text(face="bold", color="#993333",                                                                                                                              size=22, angle=0),
                                                                               axis.text.y = element_text(face="bold", color="#993333",size=22, angle=0),
                                                                               axis.title.x=element_text(size=22,face='bold'),
                                                                               axis.title.y=element_text(angle=90, face='bold', size=22),
                                                                               legend.title = element_text(face = "bold"))
        aux[[2]] = ff
        return(aux[[1]])
      }
    } else if (type == "Cheby"){
      if(lag == 1){
        res = fix.fit.cheby(ts, c, 1, m)
        val = list()
        val[[1]] =  res$ts.coef[[2]]

        ff = data.frame(pacf = unlist(val))
        ff$t = rep(seq(0,1,length.out = m), dim(ff)[1]/m)
        ff$class = as.factor(rep(c(1:lag), each = m))
        theme_update(plot.title = element_text(hjust = 0.5))
        aux[[1]] = ggplot(ff, aes(x=t, y=pacf, group=class, colour = class))+ geom_line()  + ylim(-1,1) + ggtitle(title) +
          geom_segment(aes(x=0,xend=1, y=1.96*(length(ts)^(-1/2)), yend=1.96*(length(ts)^(-1/2))), linetype="dashed", color = "black", size=0.4) +
          geom_segment(aes(x=0,xend=1, y=-1.96*(length(ts)^(-1/2)), yend=-1.96*(length(ts)^(-1/2))), linetype="dashed", color = "black", size=0.4) +
          xlab("t") + ylab("pacf") + scale_colour_discrete(name  ="lag")+theme(plot.title = element_text(size=18, face="bold"),
                                                                               legend.text=element_text(size=24, face = "bold"),
                                                                               axis.text.x = element_text(face="bold", color="#993333",                                                                                                                              size=22, angle=0),
                                                                               axis.text.y = element_text(face="bold", color="#993333",size=22, angle=0),
                                                                               axis.title.x=element_text(size=22,face='bold'),
                                                                               axis.title.y=element_text(angle=90, face='bold', size=22),
                                                                               legend.title = element_text(face = "bold"))

        aux[[2]] = ff
        return(aux[[1]])

      } else{
        res = fix.fit.cheby(ts, c, 1, m)
        val = list()
        val[[1]] =  res$ts.coef[[1+1]]
        for(b in 2:lag){
          res = fix.fit.cheby(ts, c, b, m)
          val[[b]] = res$ts.coef[[b+1]]
        }
        ff = data.frame(pacf = unlist(val))
        ff$t = rep(seq(0,1,length.out = m), dim(ff)[1]/m)
        ff$class = as.factor(rep(c(1:lag), each = m))
        theme_update(plot.title = element_text(hjust = 0.5))
        aux[[1]] = ggplot(ff, aes(x=t, y=pacf, group=class, colour = class))+ geom_line()  + ylim(-1,1) + ggtitle(title) +
          geom_segment(aes(x=0,xend=1, y=1.96*(length(ts)^(-1/2)), yend=1.96*(length(ts)^(-1/2))), linetype="dashed", color = "black", size=0.4) +
          geom_segment(aes(x=0,xend=1, y=-1.96*(length(ts)^(-1/2)), yend=-1.96*(length(ts)^(-1/2))), linetype="dashed", color = "black", size=0.4) +
          xlab("t") + ylab("pacf") + scale_colour_discrete(name  ="lag")+theme(plot.title = element_text(size=18, face="bold"),
                                                                               legend.text=element_text(size=24, face = "bold"),
                                                                               axis.text.x = element_text(face="bold", color="#993333",                                                                                                                              size=22, angle=0),
                                                                               axis.text.y = element_text(face="bold", color="#993333",size=22, angle=0),
                                                                               axis.title.x=element_text(size=22,face='bold'),
                                                                               axis.title.y=element_text(angle=90, face='bold', size=22),
                                                                               legend.title = element_text(face = "bold"))

        aux[[2]] = ff
        return(aux[[1]])
      }
    } else if (type %in% c("tri", "cos", "sin")){
      if(lag == 1){
        res = fix.fit.four(ts, c, 1, m, ops = type)
        val = list()
        val[[1]] =  res$ts.coef[[2]]

        ff = data.frame(pacf = unlist(val))
        ff$t = rep(seq(0,1,length.out = m), dim(ff)[1]/m)
        ff$class = as.factor(rep(c(1:lag), each = m))
        theme_update(plot.title = element_text(hjust = 0.5))
        aux[[1]] = ggplot(ff, aes(x=t, y=pacf, group=class, colour = class))+ geom_line()  + ylim(-1,1) + ggtitle(title) +
          geom_segment(aes(x=0,xend=1, y=1.96*(length(ts)^(-1/2)), yend=1.96*(length(ts)^(-1/2))), linetype="dashed", color = "black", size=0.4) +
          geom_segment(aes(x=0,xend=1, y=-1.96*(length(ts)^(-1/2)), yend=-1.96*(length(ts)^(-1/2))), linetype="dashed", color = "black", size=0.4) +
          xlab("t") + ylab("pacf") + scale_colour_discrete(name  ="lag")+theme(plot.title = element_text(size=18, face="bold"),
                                                                               legend.text=element_text(size=24, face = "bold"),
                                                                               axis.text.x = element_text(face="bold", color="#993333",                                                                                                                              size=22, angle=0),
                                                                               axis.text.y = element_text(face="bold", color="#993333",size=22, angle=0),
                                                                               axis.title.x=element_text(size=22,face='bold'),
                                                                               axis.title.y=element_text(angle=90, face='bold', size=22),
                                                                               legend.title = element_text(face = "bold"))

        aux[[2]] = ff

        return(aux[[1]])
      } else{
        res = fix.fit.four(ts, c, 1, m, ops = type)
        val = list()
        val[[1]] =  res$ts.coef[[1+1]]
        for(b in 2:lag){
          res = fix.fit.four(ts, c, b, m, ops = type)
          val[[b]] = res$ts.coef[[b+1]]
        }

        ff = data.frame(pacf = unlist(val))
        ff$t = rep(seq(0,1,length.out = m), dim(ff)[1]/m)
        ff$class = as.factor(rep(c(1:lag), each = m))
        theme_update(plot.title = element_text(hjust = 0.5))
        aux[[1]] = ggplot(ff, aes(x=t, y=pacf, group=class, colour = class))+ geom_line()  + ylim(-1,1) + ggtitle(title) +
          geom_segment(aes(x=0,xend=1, y=1.96*(length(ts)^(-1/2)), yend=1.96*(length(ts)^(-1/2))), linetype="dashed", color = "black", size=0.4) +
          geom_segment(aes(x=0,xend=1, y=-1.96*(length(ts)^(-1/2)), yend=-1.96*(length(ts)^(-1/2))), linetype="dashed", color = "black", size=0.4) +
          xlab("t") + ylab("pacf") + scale_colour_discrete(name  ="lag")+theme(plot.title = element_text(size=18, face="bold"),
                                                                               legend.text=element_text(size=24, face = "bold"),
                                                                               axis.text.x = element_text(face="bold", color="#993333",                                                                                                                              size=22, angle=0),
                                                                               axis.text.y = element_text(face="bold", color="#993333",size=22, angle=0),
                                                                               axis.title.x=element_text(size=22,face='bold'),
                                                                               axis.title.y=element_text(angle=90, face='bold', size=22),
                                                                               legend.title = element_text(face = "bold"))

        aux[[2]] = ff
        return(aux[[1]])
      }
    } else if (type == "Cspli"){
      if(lag == 1){
        res = fix.fit.cspline(ts, c, 1,or =or, m)
        val = list()
        val[[1]] =  res$ts.coef[[2]]

        ff = data.frame(pacf = unlist(val))
        ff$t = rep(seq(0,1,length.out = m), dim(ff)[1]/m)
        ff$class = as.factor(rep(c(1:lag), each = m))
        theme_update(plot.title = element_text(hjust = 0.5))
        aux[[1]] = ggplot(ff, aes(x=t, y=pacf, group=class, colour = class))+ geom_line()  + ylim(-1,1) + ggtitle(title) +
          geom_segment(aes(x=0,xend=1, y=1.96*(length(ts)^(-1/2)), yend=1.96*(length(ts)^(-1/2))), linetype="dashed", color = "black", size=0.4) +
          geom_segment(aes(x=0,xend=1, y=-1.96*(length(ts)^(-1/2)), yend=-1.96*(length(ts)^(-1/2))), linetype="dashed", color = "black", size=0.4) +
          xlab("t") + ylab("pacf") + scale_colour_discrete(name  ="lag")+theme(plot.title = element_text(size=18, face="bold"),
                                                                               legend.text=element_text(size=24, face = "bold"),
                                                                               axis.text.x = element_text(face="bold", color="#993333",                                                                                                                              size=22, angle=0),
                                                                               axis.text.y = element_text(face="bold", color="#993333",size=22, angle=0),
                                                                               axis.title.x=element_text(size=22,face='bold'),
                                                                               axis.title.y=element_text(angle=90, face='bold', size=22),
                                                                               legend.title = element_text(face = "bold"))

        aux[[2]] = ff

        return(aux[[1]])
      } else{
        res = fix.fit.cspline(ts, c, 1,or=or, m)
        val = list()
        val[[1]] =  res$ts.coef[[1+1]]
        for(b in 2:lag){
          res = fix.fit.cspline(ts, c, b,or=or, m)
          val[[b]] = res$ts.coef[[b+1]]
        }

        ff = data.frame(pacf = unlist(val))
        ff$t = rep(seq(0,1,length.out = m), dim(ff)[1]/m)
        ff$class = as.factor(rep(c(1:lag), each = m))
        theme_update(plot.title = element_text(hjust = 0.5))
        aux[[1]] = ggplot(ff, aes(x=t, y=pacf, group=class, colour = class))+ geom_line()  + ylim(-1,1) + ggtitle(title) +
          geom_segment(aes(x=0,xend=1, y=1.96*(length(ts)^(-1/2)), yend=1.96*(length(ts)^(-1/2))), linetype="dashed", color = "black", size=0.4) +
          geom_segment(aes(x=0,xend=1, y=-1.96*(length(ts)^(-1/2)), yend=-1.96*(length(ts)^(-1/2))), linetype="dashed", color = "black", size=0.4) +
          xlab("t") + ylab("pacf") + scale_colour_discrete(name  ="lag")+theme(plot.title = element_text(size=18, face="bold"),
                                                                               legend.text=element_text(size=24, face = "bold"),
                                                                               axis.text.x = element_text(face="bold", color="#993333",                                                                                                                              size=22, angle=0),
                                                                               axis.text.y = element_text(face="bold", color="#993333",size=22, angle=0),
                                                                               axis.title.x=element_text(size=22,face='bold'),
                                                                               axis.title.y=element_text(angle=90, face='bold', size=22),
                                                                               legend.title = element_text(face = "bold"))

        aux[[2]] = ff
        return(aux[[1]])
      }
    } else if (type %in% wavelet_basis){
      if(lag == 1){
        res = fix.fit.wavelet(ts, c, 1, m, ops = type)
        val = list()
        val[[1]] =  res$ts.coef[[2]]

        ff = data.frame(pacf = unlist(val))
        ff$t = rep(seq(0,1,length.out = m), dim(ff)[1]/m)
        ff$class = as.factor(rep(c(1:lag), each = m))
        theme_update(plot.title = element_text(hjust = 0.5))
        aux[[1]] = ggplot(ff, aes(x=t, y=pacf, group=class, colour = class))+ geom_line()  + ylim(-1,1) + ggtitle(title) +
          geom_segment(aes(x=0,xend=1, y=1.96*(length(ts)^(-1/2)), yend=1.96*(length(ts)^(-1/2))), linetype="dashed", color = "black", size=0.4) +
          geom_segment(aes(x=0,xend=1, y=-1.96*(length(ts)^(-1/2)), yend=-1.96*(length(ts)^(-1/2))), linetype="dashed", color = "black", size=0.4) +
          xlab("t") + ylab("pacf") + scale_colour_discrete(name  ="lag")+theme(plot.title = element_text(size=18, face="bold"),
                                                                               legend.text=element_text(size=24, face = "bold"),
                                                                               axis.text.x = element_text(face="bold", color="#993333",                                                                                                                              size=22, angle=0),
                                                                               axis.text.y = element_text(face="bold", color="#993333",size=22, angle=0),
                                                                               axis.title.x=element_text(size=22,face='bold'),
                                                                               axis.title.y=element_text(angle=90, face='bold', size=22),
                                                                               legend.title = element_text(face = "bold"))

        aux[[2]] = ff

        return(aux[[1]])
      } else{
        res = fix.fit.wavelet(ts, c, 1, m, ops = type)
        val = list()
        val[[1]] =  res$ts.coef[[1+1]]
        for(b in 2:lag){
          res = fix.fit.wavelet(ts, c, b, m, ops = type)
          val[[b]] = res$ts.coef[[b+1]]
        }
        ff = data.frame(pacf = unlist(val))
        ff$t = rep(seq(0,1,length.out = m), dim(ff)[1]/m)
        ff$class = as.factor(rep(c(1:lag), each = m))
        theme_update(plot.title = element_text(hjust = 0.5))
        aux[[1]] = ggplot(ff, aes(x=t, y=pacf, group=class, colour = class))+ geom_line()  + ylim(-1,1) + ggtitle(title) +
          geom_segment(aes(x=0,xend=1, y=1.96*(length(ts)^(-1/2)), yend=1.96*(length(ts)^(-1/2))), linetype="dashed", color = "black", size=0.4) +
          geom_segment(aes(x=0,xend=1, y=-1.96*(length(ts)^(-1/2)), yend=-1.96*(length(ts)^(-1/2))), linetype="dashed", color = "black", size=0.4) +
          xlab("t") + ylab("pacf") + scale_colour_discrete(name  ="lag")+theme(plot.title = element_text(size=18, face="bold"),
                                                                               legend.text=element_text(size=24, face = "bold"),
                                                                               axis.text.x = element_text(face="bold", color="#993333",                                                                                                                              size=22, angle=0),
                                                                               axis.text.y = element_text(face="bold", color="#993333",size=22, angle=0),
                                                                               axis.title.x=element_text(size=22,face='bold'),
                                                                               axis.title.y=element_text(angle=90, face='bold', size=22),
                                                                               legend.title = element_text(face = "bold"))

        aux[[2]] = ff
        return(aux[[1]])
      }
    } else{
      return(cat("Invalid option!"))
    }
  }else if (ops == "3d"){
    library(plotly)
    if(type == "Legen"){
      if(lag == 1){
        res = fix.fit.legen(ts, c, 1, m)
        val = list()
        val[[1]] =  res$ts.coef[[2]]

        ff = data.frame(pacf = unlist(val))
        ff$t = rep(seq(0,1,length.out = m), dim(ff)[1]/m)
        ff$class = as.factor(rep(c(1:lag), each = m))

        return(pacf_3dplot(ff))

      } else{
        res = fix.fit.legen(ts, c, 1, m)
        val = list()
        val[[1]] =  res$ts.coef[[1+1]]
        for(b in 2:lag){
          res = fix.fit.legen(ts, c, b, m)
          val[[b]] = res$ts.coef[[b+1]]
        }
        ff = data.frame(pacf = unlist(val))
        ff$t = rep(seq(0,1,length.out = m), dim(ff)[1]/m)
        ff$class = as.factor(rep(c(1:lag), each = m))
        return(pacf_3dplot(ff))
      }
    } else if (type == "Cheby"){
      if(lag == 1){
        res = fix.fit.cheby(ts, c, 1, m)
        val = list()
        val[[1]] =  res$ts.coef[[2]]

        ff = data.frame(pacf = unlist(val))
        ff$t = rep(seq(0,1,length.out = m), dim(ff)[1]/m)
        ff$class = as.factor(rep(c(1:lag), each = m))

        return(pacf_3dplot(ff))

      } else{
        res = fix.fit.cheby(ts, c, 1, m)
        val = list()
        val[[1]] =  res$ts.coef[[1+1]]
        for(b in 2:lag){
          res = fix.fit.cheby(ts, c, b, m)
          val[[b]] = res$ts.coef[[b+1]]
        }
        ff = data.frame(pacf = unlist(val))
        ff$t = rep(seq(0,1,length.out = m), dim(ff)[1]/m)
        ff$class = as.factor(rep(c(1:lag), each = m))
        return(pacf_3dplot(ff))
      }
    } else if (type %in% c("tri", "cos", "sin")){
      if(lag == 1){
        res = fix.fit.four(ts, c, 1, m, ops = type)
        val = list()
        val[[1]] =  res$ts.coef[[2]]

        ff = data.frame(pacf = unlist(val))
        ff$t = rep(seq(0,1,length.out = m), dim(ff)[1]/m)
        ff$class = as.factor(rep(c(1:lag), each = m))
        return(pacf_3dplot(ff))
      } else{
        res = fix.fit.four(ts, c, 1, m, ops = type)
        val = list()
        val[[1]] =  res$ts.coef[[1+1]]
        for(b in 2:lag){
          res = fix.fit.four(ts, c, b, m, ops = type)
          val[[b]] = res$ts.coef[[b+1]]
        }

        ff = data.frame(pacf = unlist(val))
        ff$t = rep(seq(0,1,length.out = m), dim(ff)[1]/m)
        ff$class = as.factor(rep(c(1:lag), each = m))
        return(pacf_3dplot(ff))
      }
    } else if (type == "Cspli"){
      if(lag == 1){
        res = fix.fit.cspline(ts, c, 1,or=or, m)
        val = list()
        val[[1]] =  res$ts.coef[[2]]

        ff = data.frame(pacf = unlist(val))
        ff$t = rep(seq(0,1,length.out = m), dim(ff)[1]/m)
        ff$class = as.factor(rep(c(1:lag), each = m))
        return(pacf_3dplot(ff))
      } else{
        res = fix.fit.cspline(ts, c, 1,or=or, m)
        val = list()
        val[[1]] =  res$ts.coef[[1+1]]
        for(b in 2:lag){
          res = fix.fit.cspline(ts, c, b,or=or, m)
          val[[b]] = res$ts.coef[[b+1]]
        }

        ff = data.frame(pacf = unlist(val))
        ff$t = rep(seq(0,1,length.out = m), dim(ff)[1]/m)
        ff$class = as.factor(rep(c(1:lag), each = m))
        return(pacf_3dplot(ff))
      }
    } else if (type %in% wavelet_basis){
      if(lag == 1){
        res = fix.fit.wavelet(ts, c, 1, m, ops = type)
        val = list()
        val[[1]] =  res$ts.coef[[2]]

        ff = data.frame(pacf = unlist(val))
        ff$t = rep(seq(0,1,length.out = m), dim(ff)[1]/m)
        ff$class = as.factor(rep(c(1:lag), each = m))
        return(pacf_3dplot(ff))
      } else{
        res = fix.fit.wavelet(ts, c, 1, m, ops = type)
        val = list()
        val[[1]] =  res$ts.coef[[1+1]]
        for(b in 2:lag){
          res = fix.fit.wavelet(ts, c, b, m, ops = type)
          val[[b]] = res$ts.coef[[b+1]]
        }
        ff = data.frame(pacf = unlist(val))
        ff$t = rep(seq(0,1,length.out = m), dim(ff)[1]/m)
        ff$class = as.factor(rep(c(1:lag), each = m))
        return(pacf_3dplot(ff))
      }
    } else{
      return(cat("Invalid option!"))
    }
    # 3d plot
  }else{
    return(cat("Invalid option!"))
  }
}
