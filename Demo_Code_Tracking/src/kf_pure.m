function kff=kf_pure(kf,z)

% kf_pure : traditional Kalman filter 
Q=kf.Q;
R=kf.R;
F=kf.F;
H=kf.H;
len=kf.len;
n=kf.n;
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
    x_=F*x;
    P_=F*P*F'+Q;
    % test
    %com_test=P_*H'*inv(R)*H-H'*inv(R)*H*P_
    % update
    K=P_*H'/(H*P_*H'+R);
    P=(eye(n)-K*H)*P_;
    %P=(eye(2)-K*H)*P_*(eye(2)-K*H)'+K*R*K';
    x=x_+K*(z(:,i)-H*x_);
    % store the data
    kgain(:,:,i)=K;
    statef_(:,i)=x_';
    statef(:,i)=x';
    covf_(:,:,i)=P_;
    covf(:,:,i)=P;
    traceP_(i)=trace(P_);
    traceP(i)=trace(P);
    ident(:,:,i)=(eye(n)-K*H)-inv(eye(n)+P_*H'*inv(R)*H); % yes and happy
end

kff.statef_=statef_;
kff.statef=statef;
kff.covf_=covf_;
kff.covf=covf;
kff.kgain=kgain;
kff.pinf_=traceP_;
kff.pinf=traceP;
kff.ident=ident;
end