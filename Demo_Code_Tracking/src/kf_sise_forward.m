function kff=kf_sise_forward(kf,un,z)

% kf : the kalman fitler instance 
Q=kf.Q;
R=kf.R;
F=kf.F;
H=kf.H;
G=kf.G;
len=kf.len;
n=kf.n;
m=kf.m;
% stored state
statef_=zeros(n,len);
statef=zeros(n,len);
covf_=zeros(n,n,len);
covf=zeros(n,n,len);

for i=1:len
    if(i==1)
        x=kf.x0;
        P=kf.P0;
    end
    % prediction 
    x_=F*x+G*un(:,i); % note that un is the noise here
    P_=F*P*F'+Q;
    % unknown input estimation
    R_t=H*P_*H'+R;
    M=(G'*H'*inv(R_t)*H*G)\(G'*H'*inv(R_t));
    d=M*(z(:,i)-H*x_);
    % update
    x__=x_+G*d;
    K=P_*H'/(H*P_*H'+R);
    x=x__+K*(z(:,i)-H*x__);

    P=(eye(n)-K*H)*((eye(n)-G*M*H)*P_*(eye(n)-G*M*H)'+G*M*R*M'*G')+K*R*M'*G';
    %P=(eye(2)-K*H)*P_*(eye(2)-K*H)'+K*R*K';

    Pd=inv(G'*H'*inv(R_t)*H*G);

    % store the data
    statef_(:,i)=x_';
    statef(:,i)=x';
    covf_(:,:,i)=P_;
    covf(:,:,i)=P;
    dist(:,i)=d;
    covd(:,i)=Pd;
    kd(:,i)=M';
end

kff.statef_=statef_;
kff.statef=statef;
kff.covf_=covf_;
kff.covf=covf;
kff.dist=dist;
kff.covd=covd;
kff.kd=kd;
end