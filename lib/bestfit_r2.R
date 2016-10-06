bestfit_r2 = function(m,multiline=NULL) {
  yobs=m$model[[1]]
  yhat=fitted(m)

  l <- list(r2 = sprintf(cor(yobs,yhat)^2, fmt='%.2f'),
            pval = sprintf(anova(m,test='F')$P[2],fmt='%.2f'))

  if(!isTRUE(multiline)){
    eq <- substitute(italic(r)^2~"="~r2*","~~italic(p)~"="~pval,l)
  } else {
    ## multiline ----
    eq <- substitute(
      atop(italic(r)^2~"="~r2,
           italic(p)~"="~pval
      )
      ,l)
  }
  ##
  as.character(as.expression(eq))
}
