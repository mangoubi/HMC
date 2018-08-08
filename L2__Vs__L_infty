clear
tic

for dd=1:14
   dd
   toc
    
clear global U

global U; d=2^dd; r=d;
U=synthetic_data_generate(d,r);
U_collection{dd}=U;

options = optimoptions(@fminunc,'Algorithm','quasi-newton','MaxIter',100);
[z,fval,exitflag,output] = fminunc(@objfun,randn(2*d,1),options);

%%%%%compute $\sqrt{L_2}\|p\|_2$%%%%%
L2(dd) = sqrt(-fval)*norm(randn(d,1));

%%%%%compute $\sqrt{L_\infty}r^{1/4} \|p\|_{u,\infty}$ $%%%%%
[z,fval,exitflag,output] = fminunc(@objfun2,randn(3*d,1),options);
L_infinity(dd) = sqrt(-fval)*(r^(0.25))*max(abs(U'*(randn(d,1))));


end
