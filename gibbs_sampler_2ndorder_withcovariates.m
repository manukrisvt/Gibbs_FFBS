function [gibbs]=gibbs_sampler_2ndorder_withcovariates(y_org,V,W,nsimu,par)

% V=100;
% W1=100;
% W2=2;
% nsimu=5000;
Tbig=length(y_org);
gibbs.theta1=zeros(nsimu,Tbig);
gibbs.theta2=zeros(nsimu,Tbig);
gibbs.V=zeros(nsimu,1);
gibbs.W_theta=zeros(nsimu,1);
gibbs.W_delta=zeros(nsimu,1);

G=[1 1 0
    0 1 0
    0 0 1];
B=[0;0];
F=[1 0 1];
% W=[W1 0;
%     0 W2];
% R=V;
z=y_org;
disp('%%%%%%%%%%%%%%%%%GIBBS RUNNING%%%%%%%%%%%%%%%%%%%%%%%');
flg=0;
for ng=1:nsimu
    ng;
    if (rem(ng*100/nsimu,20)==0)
    flg=flg+1;
    fprintf('Gibbs sampling: %d percent done\n',20*flg);
    end
    cov_X=par.XB;
    [theta]=FFBS_mk(z,G,B,F,W,V,cov_X);
    vt=y_org(1,:)-theta(1,:)-theta(3,:).*par.XB;
    wt1=theta(1,2:Tbig)-theta(1,1:Tbig-1)-theta(2,1:Tbig-1);
    wt2=theta(2,2:Tbig)-theta(2,1:Tbig-1);
    wt3=theta(3,2:Tbig)-theta(3,1:Tbig-1);
    
    
    gibbs.V(ng,1)=sample_fromIG(Tbig,vt);
    gibbs.W_theta(ng,1)=sample_fromIG(Tbig,wt1);
    gibbs.W_delta(ng,1)=sample_fromIG(Tbig,wt2);
    gibbs.W_beta(ng,1)=sample_fromIG(Tbig,wt3);
    gibbs.theta1(ng,:)=theta(1,:);
    gibbs.theta2(ng,:)=theta(2,:);
    gibbs.theta3(ng,:)=theta(3,:);
    
    W=[gibbs.W_theta(ng,1) 0 0;
        0 gibbs.W_delta(ng,1) 0
        0 0 gibbs.W_beta(ng,1)];
    V=gibbs.V(ng,1);
    
end
%%
burnin=2000;
gibbs.theta1mean=mean(gibbs.theta1(burnin:end,:),1);
gibbs.theta1var=var(gibbs.theta1(burnin:end,:),1);
gibbs.theta2mean=mean(gibbs.theta2(burnin:end,:),1);
gibbs.theta2var=var(gibbs.theta2(burnin:end,:),1);
gibbs.theta3mean=mean(gibbs.theta3(burnin:end,:),1);
gibbs.theta3var=var(gibbs.theta3(burnin:end,:),1);
% figure;
% plot(1:Tbig,y_org);
% hold on
% boundedline(1:Tbig,theta1mean,theta1var,'alpha');
% xlabel('Time');
% ylabel('Magnitude');
% legend('Org meas','Trend');
%%
end
function [otpt]=sample_fromIG(Tbig,dx)
a_hat=0.001;b_hat=0.001;
a=a_hat+Tbig/2;
b=b_hat+0.5*sum(dx.^2);
x=gamrnd(a,1/b);
otpt=1/x;
end
