function [U] = synthetic_data_generate(d,r)
%%%% "d" is dimension and "r" is number of data vectors


X = randn(d,r);
X=X.*(vecnorm(X).^(-1));

beta = randn(d,1);beta = beta/norm(beta);


h=X'*beta;

p = (1+exp(h)).^(-1);

u = rand(r,1);

%%%generate binary dependent variables%%%%
Y = (sign(p - u)');

U = X.*Y;
end
