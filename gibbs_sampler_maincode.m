% %% Gibbs sampler code.
% % 
clear all
clc
close all
for ni=1:10
Tbig=600;
theta_org=zeros(1,200);
delta_org=zeros(1,200);
y_org=zeros(1,200);

par.W_theta.var=20;
par.W_theta.mu=0;
par.W_theta.sigma=sqrt(par.W_theta.var);

par.W_delta.var=0.2;
par.W_delta.mu=0;
par.W_delta.sigma=sqrt(par.W_delta.var);

par.error.var=260;
par.error.mu=0;
par.error.sigma=sqrt(par.error.var);

par.covariate.mean=5;
par.covariate.var=0.3;
par.covariate.sigma=sqrt(par.covariate.var);
par.XB=100+normrnd(0,10,[1 Tbig]);


theta_org(1,1)=normrnd(par.W_theta.mu,par.W_theta.sigma);
delta_org(1,1)=normrnd(par.W_delta.mu,par.W_delta.sigma);
beta_org(1,1)=normrnd(par.covariate.mean,par.covariate.sigma);
y_org(1,1)=theta_org(1,1)+beta_org(1,1)*par.XB(1,1)+normrnd(par.error.mu,par.error.sigma);

for t=2:Tbig
    theta_org(1,t)=theta_org(1,t-1)+delta_org(1,t-1)+normrnd(0,par.W_theta.sigma);
    delta_org(1,t)=delta_org(1,t-1)+normrnd(0,par.W_delta.sigma);
    beta_org(1,t)=beta_org(1,t-1)+normrnd(0,par.covariate.sigma);
    y_org(1,t)=theta_org(1,t)+beta_org(1,t)*par.XB(1,t)+normrnd(0,par.error.sigma);
end
figure;plot(y_org)
prompt='Do you want to continue, press 1 to continue or zero to break';
user_ip = input(prompt);

if (user_ip==1)
data=xlsread('Book1.xlsx');
y_org=log(data(:,1)');
V=0.1;
Wmat.Wtrend=0.002;
Wmat.Wdelta=0.0002;
Wmat.Wseas=0.0002;
Wmat.Wreg=1;
nsimu=40000;
timeperiod=12;
maxterms=6;
[G,F,W]=cons_mat(timeperiod,maxterms,2,Wmat);
% W=[W1 0 0;
%     0 W2 0
%     0 0 W3
%     ];
sysmat.G=G;
sysmat.W=W;
sysmat.F=F;
sysmat.V=V;

[gibbs]=gibbs_sampler_2ndorder_withcovariates(y_org,nsimu,par,sysmat);
figure;
subplot(221)
plot(gibbs.V)
subplot(223)
histogram(gibbs.V)
vline(par.error.var);

subplot(222)
plot(gibbs.W_theta)
subplot(224)
histogram(gibbs.W_theta)
vline(par.W_theta.var);

figure;
boundedline(1:length(gibbs.theta1mean),gibbs.theta1mean,gibbs.theta1var,'alpha')
hold on
plot(y_org,'o','MarkerEdgeColor',[0.4940, 0.1840, 0.5560],'MarkerSize',2);
% hold on
% plot(theta_org,'Color',[0.8500, 0.3250, 0.0980],'LineWidth',1.5,'LineStyle','--')
legend('Var-Trend','Trend mean','Org measurement');

% figure;


break
else 
    continue
end
end
