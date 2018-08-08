function f = objfun2(z)
global U
[d,r]=size(U);


x=z(1:d);
y=z((d+1):(2*d));
v=z((2*d+1):(3*d));

f = -norm((Hessian(y,eye(d),U)*v-Hessian(x,eye(d),U)*v))/(sqrt(r)*max(abs(U'*(y-x)))*max(abs(U'*(v))));

end
