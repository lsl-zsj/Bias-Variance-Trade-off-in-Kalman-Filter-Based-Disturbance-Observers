function kff=kf_forward(kf,u,z)

% kf : the kalman fitler instance 
Q=kf.Q;
R=kf.R;
F=kf.F;
H=kf.H;
G=kf.G;
len=kf.len;
n=kf.n;
m=kf.m;
D=kf.D;
% stored state
statef_=zeros(n,len);
statef=zeros(n,len);
covf_=zeros(n,n,len);
covf=zeros(n,n,len);

for i=1:len
    if(i==1)
        x=kf.x0;
        P=kf.P0;
        d=0;
    end
    % prediction 
    x_=F*x+G*u(:,i);
    P_=F*P*F'+Q+G*kf.D*G';
    % update
    K=P_*H'/(H*P_*H'+R);
    P=(eye(n)-K*H)*P_;
    %P=(eye(2)-K*H)*P_*(eye(2)-K*H)'+K*R*K';
    x=x_+K*(z(:,i)-H*x_);
    % dist estimate
    %d=inv(G'*G)*(G'*(x-F*x_last));
    R_t=H*P_*H'+R;
    M=(G'*H'*inv(R_t)*H*G)\G'*H'*inv(R_t);
    d= M*(z(:,i)-H*x_);
    % store the data
    kgain(:,:,i)=K;
    statef_(:,i)=x_';
    statef(:,i)=x';
    covf_(:,:,i)=P_;
    covf(:,:,i)=P;
    dist(:,i)=d;
end

kff.statef_=statef_;
kff.statef=statef;
kff.covf_=covf_;
kff.covf=covf;
kff.kgain=kgain;
kff.dist=dist;
end