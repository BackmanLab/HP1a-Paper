function fun=gamma_incomplete(x,a)
% GAMMA_INCOMPLETE evaluates the upper incomplete gamma function 
% (incomplete gamma function of the second kind) $\Gamma(a,x)$ 
% at non-negative values of the argument. This function extends the 
% MATLAB function gammainc to negative values of the parameter a.
% 
%  SYNOPSIS:  fun=gamma_incomplete(x,a)
%            
% 
%  INPUT  x       : function argument
%         a       : parameter
% 
%  OUTPUT fun     : a vector of the same length as x; it contains NaN values 
%  at places where elements of x are negative.
%  
%  REMARK This function extends the MATLAB function gammainc 
%  to negative values of the parameter a.
%  
%  EXAMPLES
%
%     x=0.01:0.01:8;
%     f=gamma_incomplete(x,1);
%     plot(x,f);
%  
% 
%     x=0.001:0.001:.1;
%     f=gamma_incomplete(x,-2.3);
%     plot(x,f);
% 
% 
%     x=0.001:0.001:.2;
%     f=gamma_incomplete(x,-1);
%     plot(x,f);
    if length(a) > 1
        error('`a` must be a scalar, not an array');
    end

    p=find(x>=0);
    q=find(x<0);
    f = gamma_inc(a,x(p));
    fun(p)=f;
    fun(q)=NaN;
    
    function res=gamma_inc(a,x)
        if a < -10
            res = NaN;  % Due to the recursive handling of a<0 we get a recursion error if a is too far below 0. This data isn't valid anyway.
        elseif a==0
            res= expint(x);    
        elseif a>0
            res=gamma(a).*gammainc(x,a,'upper');
        elseif a<0
            res=(gamma_inc(a+1,x)-x.^a.*exp(-x))./a; 
        else
            error('No case matched conditions.')
        end
    end
end