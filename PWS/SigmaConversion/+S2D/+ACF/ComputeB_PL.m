function bnr_theory    = ComputeB_PL(D,lmin,lmax,r,bnr)
% Whoever wrote this function didn't bother to add any comments. :(
 [~,is] = min(abs(r-lmin));
    [~,ie] = min(abs(r-lmax));
    if(ie<is)
        temp = is; 
        is = ie; 
        ie = temp;         
    end
bnr_theory = (r/lmin).^(D-3);
bnr_theory = bnr_theory * median(bnr(is:ie)'./bnr_theory(is:ie));    
