function XAMil=schemeAMilCPU(parms, XAMil)
% % To generate the terminal value of the proposed Adaptive Milstein method
% %             without a reference solution
t=0;
while abs(t-parms.T)>parms.Tol && t<=parms.T
    HcoreFunc=AMilCoreH(XAMil,parms);
    huse=max(parms.hmin, min(parms.hmax, HcoreFunc)); % choice of h
    if t+huse > parms.T  % last tep
        huse=parms.T-t;
    end
    WincAMil=sqrt(huse)*randn;  % Wiener increment
    if huse <= parms.hmin  % backstop case
        XAMil=schemePMil(parms, XAMil, WincAMil, huse);
    else  % explicit Milstein
        [dft,Diff,Ddiff]=coefficients(XAMil, parms);
        XAMil=XAMil+dft*huse+diffTerms(WincAMil,Diff,Ddiff,huse);
    end
    t=t+huse;
end