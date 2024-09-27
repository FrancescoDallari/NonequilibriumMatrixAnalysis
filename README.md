# NonequilibriumMatrixAnalysis

function to analyze a two-times correlation matrix from an out-of equilibrium system.
Currently performs a KWW fit on the effective g2(tw,t)-1 extracted from a correlation matrix C(t1,t2).

supported "cuts":
fast =  tw=t1; t = |t2-t1|. This method corresponds to the correlation at a given age with all the future configurations
ACS  =  tw=(t1+t2)/2; t = |t2-t1|. This method considers the effective g2 as the curve composed by the elements at same average age 
classic = calculates a g2 function from a diagonal submatrix 


    ATTENTION! ALL VALUES CALCULATED IN FRAMES
    
    input:
        TT_in = Autocorellation matrix; it should be the g2(t1,t2)-1 !
        tStep = step between cuts
        lastTw= maximum age to be used. If the mode is not 'classic' it's better to keep id a bit smaller than TT_in.shape[0]
     optional:
         firstTw = first frame to be considered
         usefulBaseline = can be None, 'auto' or a number; if None the fit is done normally, if 'Auto' a baseline value is calculated taking the average of the corner (last 20 lags) of the TT_in; if it's a number the baseline in the g2 willbe fixed at that value 
         method = way in which the cut is performed:
                'fast' -> t_w = t1 ; t = |t2-t1|
                'ACS' -> t_w=(t1+t2)/2 ; t = |t2-t1|
                'classic' -> calculate the g2 over a submatrix between lagmin and lagmax (requires a rather large tStep size)
         userbounds = boundaries for the fitting function, if it's not None it has to be a touple with the bounds for: contrast, relax. rate, exponent, baseline. 
    
    output:
        outDset =  a dataset containing the results of a KWW fit and the age axis




example of usage: 

let's say that we have a corraltion matrix called TT_all 
then we ca run 
outDset = smf.analyzeAgedepTT(TT_all,5,200)
to perform the fits in "fast" mode from tw=0 to tw=200 with steps of 5 frames

If we want to restrain the fit within certain boundaries and show in a single plot the g2 and their respective fits we can use:

user_bounds = ([0.01, 0, 0.5, -0.01],[0.44, np.Inf, 2, 0.8*0.044])

outDset = smf.analyzeAgedepTT(TT_all,5,200,
                              firstTw=0,
                              usefulBaseline=None,
                              method='fast',
                              showPlot=True,
                              userbounds=user_bounds)
