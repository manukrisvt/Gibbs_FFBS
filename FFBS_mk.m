function [lat_state]=FFBS_mk(z,G,B,F,W,V,cov_X)
% F is the measurement equation (equivalent to our F in time series class)
% G is our process matrix (equivalenet to our G in time series class)
% B is the loading matrix.
% z is the measurements in row matrix
% W is covariance matrix
% V is the measurement covariance matrix

lsim=length(z);
no_measurements=size(z,1);

%Kalman filter
% F=[1 0 cov_X(1,1)];
a0=zeros(length(G),1);
R0=diag(10^3*ones(length(G),1));
f0=F*a0;
Q0=F*R0*F'+V;
A0=R0*F'*pinv(Q0);
e0=z(1,1)-f0;
C0=R0-A0*A0'*e0;
m0=a0+A0*e0;

C=zeros(length(G),length(G),lsim);
m=zeros(length(G),lsim);

m(:,1)=m0;
C(:,:,1)=C0;

for it = 2:lsim
%     F=[1 0 cov_X(1,it)];
    a=G*m(:,it-1);
    R=G*C(:,:,it-1)*G'+W;
    f=F*a;
    Q=F*R*F'+V;
    K=R*F'*pinv(Q);
    e=z(1,it)-f;
    
    m(:,it)=a+K*e;
    C(:,:,it)=R-K*Q*K';
   
end
% figure(fig)
% subplot(2,1,1);
% plot(t,X_(1,:));
% subplot(2,1,2)
% plot(t,X_(2,:));
% refresh(fig);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Kalman smooth - Rauch?ung?triebel
theta=zeros(length(G),lsim);
C(:,:,lsim) = (C(:,:,lsim) +C(:,:,lsim).') / 2;
theta(:,lsim)=mvnrnd(m(:,lsim),C(:,:,lsim));
for it = lsim-1:-1:1
    dmvar=pinv(C(:,:,it))+(G'*pinv(W)*G);
    Ht=pinv(dmvar);
    ht=Ht*(pinv(C(:,:,it))*m(:,it)+G'*pinv(W)*theta(:,it+1));
    Ht = (Ht +Ht.') / 2;
    theta(:,it)=mvnrnd(ht,Ht);
end

lat_state=theta;
end