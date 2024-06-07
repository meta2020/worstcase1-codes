library(meta)

ml.lnDOR <- metabin(TP, with(example2.1, TP+FN), FP, with(example2.1, TN+FP), data = example2.1, sm ="OR", fixed = FALSE)
tf.lndor <- trimfill(ml.lnDOR, fixed = FALSE)
funnel(tf.lndor, xlab = "lnDOR", 
       #pch = c(shapes, rep(1, sum(tf.lndor[["trimfill"]]))),
       # bg = example2.1$cutoff.grp+1,
       # col = c(example2.1$cutoff.grp+1, rep(1, sum(tf.lndor[["trimfill"]]))),
       cex = 1.5, studlab = FALSE, backtransf = FALSE,
       # fixed = FALSE, random = TRUE,
       ref = exp(ml.lnDOR$TE.random)
       # level = 0.95,
       # contour = c(0.9, 0.95, 0.99)
)

mtext("Standar Error", side = 2, line = 2, at = 1, cex = 0.7)

## Sens
ml.se <- metabin(TP, with(example2.1, TP+FN), FN, with(example2.1, TP+FN), data = example2.1, sm ="OR", fixed = FALSE)
tf.se <- trimfill(ml.se, fixed = FALSE)
tf.se[["TE"]] <- tf.se[["TE"]]/2
tf.se[["seTE"]] <- tf.se[["seTE"]]/sqrt(2) 
funnel(tf.se, xlab = "logit-sensitivity",
       #pch = c(shapes, rep(1, sum(tf.se[["trimfill"]]))),
       yaxt = FALSE,
       # bg = example2.1$cutoff.grp+1,
       cex = 1.5, studlab = FALSE, backtransf = FALSE, 
       ref = exp(ml.se$TE.random)
       # level = 0.95,
       # contour = c(0.9, 0.95, 0.99)
)


## Spec
ml.sp <- metabin(TN, with(example2.1, TN+FP), FP, with(example2.1, TN+FP), data = example2.1, sm ="OR", fixed = FALSE)
tf.sp <- trimfill(ml.sp, fixed = FALSE)
tf.sp[["TE"]] <- tf.sp[["TE"]]/2
tf.sp[["seTE"]] <- tf.sp[["seTE"]]/sqrt(2) 
funnel(tf.sp, xlab = "logit-specificity", fixed = FALSE,
       #pch = c(shapes, rep(1, sum(tf.sp[["trimfill"]]))),
       # bg = example2.1$cutoff.grp+1,
       cex = 1.5, studlab = FALSE, backtransf = FALSE,
       ref = exp(ml.sp$TE.random)
       # level = 0.95,
       # contour = c(0.9, 0.95, 0.99)
)
