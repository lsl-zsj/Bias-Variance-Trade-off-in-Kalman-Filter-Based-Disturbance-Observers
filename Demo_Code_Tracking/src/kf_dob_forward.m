function kff=kf_dob_forward(kf,u,z)

% kf : the kalman fitler instance 
Q=kf.Q;
R=kf.R;
F=kf.F;
H=kf.H;
G=kf.G;
Ga=kf.Ga;
len=kf.len;
n=kf.n;
m=kf.m;
% stored state
statef_=zeros(n,len);
statef=zeros(n,len);
covf_=zeros(n,n,len);
covf=zeros(n,n,len);
tic
for i=1:len
    if(i==1)
        x=kf.x0;
        P=kf.P0;
    end
    % prediction 
    x_=F*x+Ga*u(:,i);
    P_=F*P*F'+Q;
    % update
    K=P_*H'/(H*P_*H'+R);
    P=(eye(n)-K*H)*P_;
    %P=(eye(2)-K*H)*P_*(eye(2)-K*H)'+K*R*K';
    x=x_+K*(z(:,i)-H*x_);
    % store the data
    statef_(:,i)=x_';
    statef(:,i)=x';
    covf_(:,:,i)=P_;
    covf(:,:,i)=P;
    kd(:,i)=K(1,1:2)';
    pd(i)=covf(1,1,i);
end
tcost=toc;

kff.tcost=tcost;
kff.statef_=statef_;
kff.statef=statef;
kff.covf_=covf_;
kff.covf=covf;
kff.kd=kd;
kff.pd=pd;
end