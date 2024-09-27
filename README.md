# nonequilibriumMatrixAnalysis

function to analyze a two-times correlation matrix from an out-of equilibrium system.
Currently performs a KWW fit on the effective g2(tw,t)-1 extracted from a correlation matrix C(t1,t2).

supported "cuts":
fast =  tw=t1; t = |t2-t1|. This method corresponds to the correlation at a given age with all the future configurations
ACS  =  tw=(t1+t2)/2; t = |t2-t1|. This method considers the effective g2 as the curve composed by the elements at same average age 
classic = calculates a g2 function from a diagonal submatrix 
