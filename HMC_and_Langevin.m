%%%%%%%HMC and Langevin for Bayesian Logistic regression%%%%%%%%%




%%%%%loading the data%%%%%%%%

%load('musk_dataset.mat')

%X_data=table2array(musk(:,1:166));
%Y_data=table2array(musk(:,167));
%U = ((Y_data-0.5)*2).*X_data;



%%%%%Specify which algorithm to use%%%%%%%%

%Algorithm ("1" = Langevin, "0" = HMC)
Langevin = 0

%%%%%Specify whether to use Metropolis "accept/reject" filter%%%%%%%%

%Metropolis filter ("1" = Metropolis filter, "0" = no filter)
Metropolis =0


%%%%specify number of Markov chain steps%%%
k_max = 10000

%%%%specify discretization level%%%
eta = 0.1



[d,r]=size(U);




%Starting point
X(:,1) = randn(d,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%constants needed to compute the gradient
%2xBayesian prior inverse covariance matrix
%A = diag([ones(d/2,1);5*ones(d/2,1)]);%1*eye(d);
A = eye(d);%diag(1:d);


%%%%%%%%%%%%%%%%%%%%%%%Due to numerical stability issues, we scale down U
%%%%%%%%%%%%%%%%%%%%%%%and A, then scale X back up again at the end.

%U = U/norm(U);

scale=1;%norm(A);

%inverse Mass matrix
M = scale*eye(d);%A; %inv(eye(d));
MM = sqrtm(inv(M));

%Integration time (we set T= eta for Langevin)
if Langevin==1
    T=eta
else
T = (pi/3)*1/sqrt(scale*norm(M));
end

%Number of numerical steps
i_max = floor(T/eta)+1;






for k = 1:k_max
    if mod(k,1000)==1
        k
        toc
        X(:,k) = randn(d,1);
    end
    q(:,1) = X(:,k);
    p(:,1) = MM*randn(d,1);
    
    
    
    for i = 1:(i_max -1)
        
        qq = q(:,i);
        pp = p(:,i);
        
        
        g = grad(qq,A,U);

        
        %%%
        
        
        
    q(:,i+1)  = qq + eta*M*pp - 0.5*(eta^2)*M*g;
    
                    g_2 = grad(q(:,i+1),A,U);

    p(:,i+1) = pp - eta*g -0.5*eta*(g_2-g);
    end
    
    
    %Metropolis step
 
    
    H1 = Ham(q(:,1),p(:,1),A,U,M);
    
    H2 = Ham(q(:,i_max),p(:,i_max),A,U,M);
    

    alpha = rand;
           
    
    if Metropolis == 1
    if alpha < exp(H1-H2)
    X(:,k+1) = q(:,i_max);
    accept(k) = 1;
    else
    X(:,k+1) = X(:,k);
    accept(k) = 0;
    end
    else
        X(:,k+1) = q(:,i_max);
        accept(k) = 1;
    end
    
end


%print average acceptance probability
if Metropolis == 1
mean(accept);
end

%c = linspace(1,10,length(X(1,:)));
%scatter(X(1,:), X(2,:),[],c)



%Autocorellation of l1 test function   (not for Octave)
v=ones(d,1);
v=v/norm(v);
xx=X(:,k_max/2:k_max)'*v;z=autocorr(repelem(xx,floor(T/eta)),100);autocorr(repelem(xx,floor(T/eta)),100)
mean(xx);
std(xx);


%%%%%%define gradient of potential
function G = grad(x,AA,UU)

  G = AA*x;
       
        Q = x'*UU;
        G = G + sum(((1+(exp(1)).^(-Q)).^(-1)).*UU,2);
end


%%%%%define the Hamiltonian (only needed for Metropolis step)
function HH = Ham(x,v,AA,UU,IMass)
        Q = x'*UU;
        
  HH = 0.5*x'*AA*x - sum(log(((1+(exp(1)).^(Q)).^(-1))),2)   +  0.5*v'*IMass*v;
end
