TCGA_PathObj = NetFeatureCalculate(TCGA_PathObj)

library(poweRlaw)
x = x[x != 0]
m_pl = poweRlaw::displ$new(x)
est = poweRlaw::estimate_xmin(m_pl)
poweRlaw::m_pl$setXmin(est)
bs_p = poweRlaw::bootstrap_p(m_pl, threads = nthreads)$p
