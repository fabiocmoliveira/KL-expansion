
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% CAMPINAS 2022 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Karhunen-Loève Enpansion %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Fábio César Miranda de Oliveira %%%%%%%%%%%%%%%%%%%%

%Este código realiza e expansão de Karhunen-Loève para uma resposta em 
%frequência de uma barra excitada em uma extremidade




clear all
close all
clc

load(['XKL1.mat']) %FRF das n observações
load(['H.mat']) %FRF das n observações
%load(['viga8_1.mat']) %Dados de uma viga no tempo utilizado para obter o vetor do tempo

Nobs = 512;
Nt = 32;

% dt=c1x(2)-c1x(1);   %incremento no tempo
% Df=1/(2*Nobs*dt);
% Freq=((1:Nobs)-1)*Df;

load(['Freq.mat'])
Df = f(2)-f(1);
writematrix(f,'Frequencia.csv');
writematrix(res,'FRF_KL1.csv');
writematrix(H,'FRF_H.csv');


X=(res);
X=X';

t = f;

for i=1:1:Nt
    plot(t,20*log10(X(:,i)),'.');
    hold on
end
title('Amostras','fontsize',22)
xlabel('Freq. [Hz]','fontsize',22)
ylabel('FRF [dB]','fontsize',22)


uns = ones(Nt,1);

Mu = [];

for i = 1:1:Nt
Mu(:,i) = X*uns/(Nt);
end

X0 = X-Mu;

C_hat = (1/(Nt))*X0*X0'*Df;

[V, D, W] = eig(C_hat);

D=diag(sort(diag(D),'descend')); % make diagonal matrix out of sorted diagonal values of input D
[c, ind]=sort(diag(D),'descend'); % store the indices of which columns the sorted eigenvalues come from
V=V(:,ind); % arrange the columns in this order

D = abs(diag(D));

figure
subplot(2,1,1)
plot(1:Nobs,D)

subplot(2,1,2)
plot(1:20,D(1:20))

Z = [];

for i=1:1:Nt
Z = [Z; 1/D(i)*X0(:,i)'.*V(i,:)];
end


figure
plot(t,(20*log10(Mu(:,1))),'LineWidth',2);
hold on
plot(t,(20*log10(H(1:512))),'LineWidth',2);
legend('Estimation of H','H','fontsize',22)
xlabel('Freq. [Hz]','fontsize',22)
ylabel('FRF [dB]','fontsize',22)
title('Estimation of H (discrete KL mean)','fontsize',22)


XKL = zeros(Nt,Nobs);

Mu=Mu';
XKL = Mu;
for i = 1:1:Nt
for j = 1:1:Nobs
for k = 1:1:1
XKL(i,j) = XKL(i,j) + sqrt(D(i))*V(j,k)'.*Z(k,j);
end
end
end

figure
subplot(2,2,1)
ksdensity(abs(20*log10(XKL(:,1))))
title('fdp f=0')

subplot(2,2,2)
ksdensity(abs(20*log10(XKL(:,220))))
title('fdp f=6.45')

subplot(2,2,3)
ksdensity(abs(20*log10(XKL(:,295))))
title('fdp t=8.65')

subplot(2,2,4)
ksdensity(abs(20*log10(XKL(:,480))))
title('fdp t=14.05')

XKL_max = [];
XKL_min = [];

for i = 1:1:Nobs
    XKL_max = [XKL_max max((XKL(:,i)))];
    XKL_min = [XKL_min min((XKL(:,i)))];
end



figure
plot(t,20*log10(XKL_min),'LineWidth',2)
hold on
plot(t,20*log10(XKL_max),'LineWidth',2)
plot(t,20*log10(H(1:1:512)),'LineWidth',2)
title('X_{KL} with k = 32','fontsize',22)
legend({'X_{KL,max,tot}','X_{KL,min,tot}','H'},'fontsize',12)
xlabel('Freq. [Hz]','fontsize',22)
ylabel('FRF [dB]','fontsize',22)
title('Estimation of H (discrete KL mean)','fontsize',22)



es = mahal(20*log10(abs(Mu(1,:)')),20*log10(abs(XKL_max))');
ei = mahal(20*log10(abs(Mu(1,:)')),20*log10(abs(XKL_min))');

esm = sum(es)/length(es)
eim = sum(ei)/length(ei)

e = sum((20*log10(abs(Mu(1,:)'))-20*log10(abs(XKL_max))').^2);
e = e+sum((20*log10(abs(Mu(1,:)'))-20*log10(abs(XKL_min))').^2);
e = sqrt(e)/(2*length(XKL_max))

XKL_min = smooth(XKL_min);
XKL_max = smooth(XKL_max);

figure
plot(t,20*log10(XKL_min),'LineWidth',2)
hold on
plot(t,20*log10(XKL_max),'LineWidth',2)
title('X_{KL} with k = 1024','fontsize',22)
xlabel('Freq. [Hz]','fontsize',22)
ylabel('FRF [dB]','fontsize',22)
legend({'X_{KL,max,sys}','X_{KL,min,sys}'},'fontsize',12)
set(gca,'FontSize',12)

x = [1 2 3 5 10 20 31];
y = [0.0332 0.1410 0.2459 0.3882 0.5807 0.7620 0.8953];


figure
plot(x,y)
title('Estimated error vs number of applied eigenvectors','fontsize',22)
xlabel('Number of applied eigenvectors','fontsize',22)
ylabel('Estimated error [dB]','fontsize',22)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%           KL por autofunções        %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return
load('raizes.mat')
load('lamb.mat')


x=0:1:1024;
y=zeros(1,length(x));

eps = normrnd(0,1,[1025,length(x)]);

for i=1:2:50
    y=y+sqrt((D(i)))*cos(raizes(i)*(x-1/2))./(sqrt(1/2+sin(raizes(i))./(2*raizes(i))));
    hold on
end

for i=2:2:50
    y=y+sqrt((D(i)))*sin(raizes(i)*(x-1/2))./(sqrt(1/2-sin(raizes(i))./(2*raizes(i))));
    hold on
end

figure
plot(x,20*log10(y))


x=0:0.01:1;
y=zeros(1,length(x));

eps = normrnd(0,1,[1025,length(x)]);

for i=1:2:50
    y=y+sqrt((D(i)))*cos(raizes(i)*(x-1/2))./(sqrt(1/2+sin(raizes(i))./(2*raizes(i))));
    hold on
end

for i=2:2:50
    y=y+sqrt((D(i)))*sin(raizes(i)*(x-1/2))./(sqrt(1/2-sin(raizes(i))./(2*raizes(i))));
    hold on
end

figure
plot(x,20*log10(y))


x=1:1:1024;
y=zeros(1,length(x));

for i=1:2:50
    y=y+sqrt((D(i)))*cos(raizes(i)*(x-1/2))./(sqrt(1/2+sin(raizes(i))./(2*raizes(i))));
    hold on
end

for i=2:2:50
    y=y+sqrt((D(i)))*sin(raizes(i)*(x-1/2))./(sqrt(1/2-sin(raizes(i))./(2*raizes(i))));
    hold on
end


y2 = [];

for i = 1:1:length(y)
    y2(:,i)=y(i)*eps(:,1);
end


figure
for i = 1:1:1000 
    plot(x,20*log10(y2(i,:))),'.';
    hold on
end


x = 0:0.01:1;

y=zeros(1,length(x));

for i=1:2:50
    y=y+sqrt((D(i)))*cos(raizes(i)*(x-1/2))./(sqrt(1/2+sin(raizes(i))./(2*raizes(i))));
    hold on
end

for i=2:2:50
    y=y+sqrt((D(i)))*sin(raizes(i)*(x-1/2))./(sqrt(1/2-sin(raizes(i))./(2*raizes(i))));
    hold on
end

eps = normrnd(0,1,[1000,1]);

y2 = [];

for i = 1:1:length(y)
    y2(:,i)=y(i)*eps(:,1);
end


figure
for i = 1:1:1000 
    plot(x,20*log10(y2(i,:))),'.';
    hold on
end



