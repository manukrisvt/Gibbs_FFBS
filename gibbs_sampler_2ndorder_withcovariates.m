function [gibbs]=gibbs_sampler_2ndorder_withcovariates(y_org,nsimu,par,sysmat);

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
gibbs.W_seas=zeros(nsimu,1);


B=[0;0];
F=sysmat.F;
G=sysmat.G;
W=sysmat.W;
V=sysmat.V;
% W=[W1 0;
%     0 W2];
% R=V;
z=y_org;
disp('%%%%%%%%%%%%%%%%%GIBBS RUNNING%%%%%%%%%%%%%%%%%%%%%%%');
flg=0;
W_dm=zeros(length(G),length(G));
for ng=1:nsimu
    ng;
    if (rem(ng*100/nsimu,20)==0)
    flg=flg+1;
    fprintf('Gibbs sampling: %d percent done\n',20*flg);
    end
%     cov_X=par.XB;
    [theta]=FFBS_mk(z,G,B,F,W,V);
    vt=y_org(1,:)-F*theta;
    wt1=theta(1,2:Tbig)-theta(1,1:Tbig-1)-theta(2,1:Tbig-1);
    wt2=theta(2,2:Tbig)-theta(2,1:Tbig-1);
    wt3_dm=[];
    for dt=1:length(G)-2
        wt3_dm=[wt3_dm;theta(dt+2,2:Tbig)-theta(dt+2,1:Tbig-1)];
    end
    wt3=mean(wt3_dm,1);
    
    
    gibbs.V(ng,1)=sample_fromIG(Tbig,vt);
    gibbs.W_theta(ng,1)=sample_fromIG(Tbig,wt1);
    gibbs.W_delta(ng,1)=sample_fromIG(Tbig,wt2);
    gibbs.W_seas(ng,1)=sample_fromIG(Tbig,wt3);
    gibbs.theta1(ng,:)=theta(1,:);
    gibbs.theta2(ng,:)=theta(2,:);
    gibbs.theta3(ng,:)=theta(3,:);
    gibbs.theta4(ng,:)=theta(4,:);
    gibbs.theta5(ng,:)=theta(5,:);
    gibbs.theta6(ng,:)=theta(6,:);
    gibbs.theta7(ng,:)=theta(7,:);
    gibbs.theta8(ng,:)=theta(8,:);
    gibbs.theta9(ng,:)=theta(9,:);
    gibbs.theta10(ng,:)=theta(10,:);
    gibbs.theta11(ng,:)=theta(11,:);
    gibbs.theta12(ng,:)=theta(12,:);
    gibbs.theta13(ng,:)=theta(13,:);
    
    
    W_dm(1,1)=gibbs.W_theta(ng,1);
    W_dm(2,2)= gibbs.W_delta(ng,1);
    W_dm(3:length(G),3:length(G))= diag(gibbs.W_seas(ng,1)*ones(length(G)-2,1));
     
  
    
end
%%
burnin=500;
gibbs.theta1mean=mean(gibbs.theta1(burnin:end,:),1);
gibbs.theta1var=var(gibbs.theta1(burnin:end,:),1);
gibbs.theta2mean=mean(gibbs.theta2(burnin:end,:),1);
gibbs.theta2var=var(gibbs.theta2(burnin:end,:),1);
gibbs.theta3mean=mean(gibbs.theta3(burnin:end,:),1);
gibbs.theta3var=var(gibbs.theta3(burnin:end,:),1);
gibbs.theta5mean=mean(gibbs.theta5(burnin:end,:),1);
gibbs.theta5var=var(gibbs.theta5(burnin:end,:),1);
gibbs.theta7mean=mean(gibbs.theta7(burnin:end,:),1);
gibbs.theta7var=var(gibbs.theta7(burnin:end,:),1);
gibbs.theta9mean=mean(gibbs.theta9(burnin:end,:),1);
gibbs.theta9var=var(gibbs.theta9(burnin:end,:),1);
gibbs.theta11mean=mean(gibbs.theta11(burnin:end,:),1);
gibbs.theta11var=var(gibbs.theta11(burnin:end,:),1);
gibbs.theta13mean=mean(gibbs.theta13(burnin:end,:),1);
gibbs.theta13var=var(gibbs.theta13(burnin:end,:),1);

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
