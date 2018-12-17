function [yy, Q] = fixed_mass_deconvolve(W, xx)
    
    global fmax Term1boot fctargs  penalite Termequality;
    
    n = length(W);
    varW = var(W);
    dx=xx(2)-xx(1);
    
    % FIX THIS TO CHOOSE GOOD VALUES AS IN DECONVOLVE PACKAGE
    longt = 50;
    a = -8;
    b = 8;
    tt = a:((b - a) / longt):b;
    
    % Compute phiW and psiW
    [tt1, tt2, rehatphiW, imhatphiW, normhatphiW] = computephiW(tt,longt,W,n);
    tt = tt1:(tt2-tt1)/longt:tt2;
    [~, ~, sqabshatpsi] = computepsiW(tt,W,n);

    %----------------------------------------------------------------------
    % Optimization
    %----------------------------------------------------------------------
    
    %Compute a small grid and the longer grid of x_j-values on which we want to find our solution

    m=round(2*sqrt(n));
    mlong=round(5*sqrt(n));
    xgridlongW=sort(unifrnd(min(W),max(W),[1,mlong]));

    deltam=round((mlong-1)/m);
    indm=1:deltam:mlong;
    m=length(indm);
    xgridW=xgridlongW(indm);
    xgrid=xgridW;

    %-----------------------------------------------------------
    % First step: compute solution on a small grid of x_j values
    %-----------------------------------------------------------

    %Using the small grid of x_j-values, find the minimum value of the integral (our objective function to minimise) without any constraint on the variance of X
    %Solution (p_j's) found by minimising the objective function without imposing any variance constraint. Hope is to find the minimum value that the objective function can take

    fctargs=struct('xgrid',xgrid,'n',n,'tt',tt,'rehatphiW',rehatphiW,'imhatphiW',imhatphiW,'normhatphiW',normhatphiW,'sqabshatpsi',sqabshatpsi,'varW',varW);

    display('Doing First Minimization')
    [psolA,fvalA]=findpsolBoot2(m,'fobjUnconst');
    display('Finished First Minimization')
    Termequality

    %Compute the corresponding value of the objective function (fvalA does NOT provide this as it contains some penalties we do not want to take into account here)
    fmax=fobjUnconstB(psolA(1:(m-1)));
    %display(['fmax = ',fmax])

    %Now recompute the estimator but impose a minimum variance constraint. 
    %That is, find p_j's that minimise the variance of X under the constraint that the integral (main objective function) is no larger than fmax computed above
    fctargs=struct('xgrid',xgrid,'n',n,'tt',tt,'rehatphiW',rehatphiW,'imhatphiW',imhatphiW,'normhatphiW',normhatphiW,'sqabshatpsi',sqabshatpsi,'varW',varW);
    pstartB=psolA;

    display('Doing Second Minimization')
    [psolB,fvalB]=findpsolBoot(m,'fobjBoot',pstartB);
    display('Finished Second Minimization')

    %--------------------------------------------------------------------------------------------------------------------------
    % Second step: recompute solution on a larger grid of x_j values in a neighborhood of what we just found on the small grid
    %--------------------------------------------------------------------------------------------------------------------------

    %Third compute psol on a longer grid taking this psol as a starting point with some adjustments to account for the increase of grid size:

    m=mlong;
    xgridW=xgridlongW;
    xgrid=xgridW;
    pstart=0*xgrid;
    pstartB=pstart;

    pstart(indm)=psolA;
    pstartB(indm)=psolB;

    for k=1:deltam
        indm2=indm+k;
        indm2(indm2>mlong)=[];

        pstart(indm2)=psolA(1:length(indm2));
        pstartB(indm2)=psolB(1:length(indm2));
    end

    pstart=pstart/sum(pstart);
    pstartB=pstartB/sum(pstartB);



    fctargs=struct('xgrid',xgrid,'n',n,'tt',tt,'rehatphiW',rehatphiW,'imhatphiW',imhatphiW,'normhatphiW',normhatphiW,'sqabshatpsi',sqabshatpsi,'varW',varW);

    %As above, first compute p_j's without any constraint on variance of X:
    [psolA,fvalA]=findpsolWithStart(pstart,m,'fobjUnconst');
    %Compute the corresponding value of the objective function (fvalA does NOT provide this as it contains some penalties we do not want to take into account here)
    fmax=fobjUnconstB(psolA(1:(m-1)));


    %Now recompute the estimator but impose a minimum variance constraint. 
    %That is, find p_j's that minimise the variance of X under the constraint that the integral (main objective function) is no larger than fmax computed above
    [psolB,fvalB]=findpsolBootWithStart(m,'fobjBoot',pstartB);
    psol=psolB;
    Termequality



    %---------------------------------------------------
    %---------------------------------------------------
    % Now turn the discrete distribution into a density
    %---------------------------------------------------
    %---------------------------------------------------
    % [pj_new,xj_new] = simplify_masses(psol,xgrid);
    % psol = pj_new;
    % xgrid = xj_new;
    % m = length(psol);
    %----------------------------
    % Estimate phi_X and phi_U
    %----------------------------

    [rephip,imphip,normphip]=computephiX(tt,xgrid,psol,n);
    hatphiU=normhatphiW./normphip;

    %----------------------------------------------------------------------------------------
    %Estimate var(U): approximate phi_U by poly of degree 2, and estimate varU by -2 phi_U''
    %----------------------------------------------------------------------------------------

    %For this we use a finer grid than the grid tt
    ttBB=tt1:(tt2-tt1)/200:tt2;
    OO=outerop(ttBB,W,'*');

    %Estimate empirical characersitic function of W
    rehatphiWBB=sum(cos(OO),2)/n;
    imhatphiWBB=sum(sin(OO),2)/n;

    clear OO;
    normhatphiWBB=sqrt(rehatphiWBB.^2+imhatphiWBB.^2);
    clear rehatphiWBB,imhatphiWBB;

    [rephipBB,imphipBB,normphipBB]=computephiX(ttBB,xgrid,psol,n);
    hatphiUBB=normhatphiWBB./normphipBB;

    tvec=ttBB(hatphiUBB'>=1-0.05);
    phiUtvec=hatphiUBB(hatphiUBB'>=1-0.05);
    pp=polyfit(tvec,phiUtvec',2);
    hatvarU=-2*pp(1);



    %Then compute density estimator as indicated in the paper

    tlim=[min(tt),max(tt)];

    %Adjust estimator of phi_U as recommended in the paper
    ppphiU=spline(tt,hatphiU);

    %Compute bandwidth as recommended in the paper

    hPIc=PI_deconvUestth4(W,tlim,ppphiU,hatvarU);

    %Compute density estmator
    fXdecc=fXKernDec2(xx,hPIc,W,tlim,ppphiU,hatvarU);

    %Remove negative parts and resclae to integrate to 1
    fXdecc(fXdecc<0)=0*fXdecc(fXdecc<0);
    fXdecc=fXdecc/sum(fXdecc)/dx;
    fXtest=fXdecc;
    
    % ---------------------------------------------------------------------
    % Return values
    % ---------------------------------------------------------------------
    
    yy = fXtest;
    Q.Support = xgrid;
    Q.ProbWeights = psol;

end