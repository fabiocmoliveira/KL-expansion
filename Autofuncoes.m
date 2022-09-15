clear all
close all
clc

b=1;
w=[];

w=[];

f= @(w)((1-b*w*tan(w/2))*(b*w+tan(w/2)));
w0=[1 2];
sol=fzero(f,w0);
w=[sol];

for i=1:1:50
f= @(w)((1-b*w*tan(w/2))*(b*w+tan(w/2)));
w0=[3+3.1*(i-1) 6.2+3.14*(i-1)];
sol=fzero(f,w0);
w=[w;sol];
end

raizes = w;

w = -100:0.1:100;
y=(1-b*w.*tan(w/2)).*(b*w+tan(w/2));

figure
plot(w,y)
ylim([-1 1])


lamb = [];

sigma2u = 1;

for i=1:1:50
    lamb = [lamb;sigma2u*2*b/(1+(raizes(i)*b)^2)];
end

x=0:0.01:1;
y2 = zeros(1,length(x));

figure
subplot(3,1,2)
for i=1:2:10
    y=sqrt(lamb(i))*cos(raizes(i)*(x-1/2))./(sqrt(1/2+sin(raizes(i))./(2*raizes(i))));
    y2=y2+y;
    plot(x,y)
    hold on
end

for i=2:2:10
    y=sqrt(lamb(i))*sin(raizes(i)*(x-1/2))./(sqrt(1/2-sin(raizes(i))./(2*raizes(i))));
    plot(x,y)
    y2=y2+y;
    hold on
end

subplot(3,1,1)
plot(1:1:10,lamb(1:10),'+','LineWidth',1.5,'color','r')

subplot(3,1,3)
plot(x,y2,'+','LineWidth',1.5)

save('raizes.mat','raizes');
save('lamb.mat','lamb');