# R CMD SHLIB fit_nb_sf.c lbfgsb_dr.c

dyn.load( "fit_nb_sf.so")
.Call( "fit_nb_sf", 
	c( 0, 0, 0, 0, 0, 0, 1, 2, 6, 6 ), 
	c(0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.30, 1.50, 1.90, 2.00) )