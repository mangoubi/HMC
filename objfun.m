function f = objfun(z)
global U
[d,r]=size(U);

x=z(1:d);
y=z((d+1):(2*d));

f = -svds((Hessian(y,eye(d),U)-Hessian(x,eye(d),U)),1)/norm(y-x);

end
