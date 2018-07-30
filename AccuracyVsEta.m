tic;
clear
Auto=[];

Eta = 0.1:0.05:0.6

for J = 1:length(Eta)
clearvars -except J  Auto Eta X_eta;

eta = Eta(J)

load('synthetic_data.mat');

Langevin=1;Metropolis=0;
if Langevin==1
    T=eta
else
T = (pi/3);
end
i_max = floor(T/eta)+1;


k_max=floor(100000/i_max)



HMC_synthetic;toc

Auto = [Auto, sum(z)]
X_eta{J} = X;
end

plot(Eta,Auto)
