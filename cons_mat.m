function [G,F,W]=cons_mat(timeperiod,nt_u,order,Wmat)

if rem(timeperiod,2)==0
    nt=nt_u-1;
    d=2*nt_u-1+order;
    G=zeros(d,d);
    G(end,end)=-1;
    idx=order+1:2:d;
    F=zeros(1,d);
    F(1,1:2)=[1 0];
    F(idx)=1;
else
    nt=nt_u;
    d=2*nt_u+order;
    G=zeros(d,d);
    idx=order+1:2:d;
    F=zeros(1,d);
    F(1,1:2)=[1 0];
    F(idx)=1;
end
    W=zeros(d,d);
if order==2
    G(1:2,1:2)=[1 1;
        0 1];
    W(1:2,1:2)=[Wmat.Wtrend 0;
                0 Wmat.Wdelta];
    F(1:2)=[1 0];
end

for ni=1:nt
    str_id=order+(ni-1)*2+1;
    end_id=order+(ni-1)*2+2;
    
    G(str_id:end_id,str_id:end_id)=[cos(2*pi*ni/timeperiod) sin(2*pi*ni/timeperiod);
                                -sin(2*pi*ni/timeperiod) cos(2*pi*ni/timeperiod)];
    W(str_id,str_id)=Wmat.Wseas;
    

end

