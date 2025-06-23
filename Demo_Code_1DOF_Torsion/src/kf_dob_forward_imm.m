function kff=kf_dob_forward_imm(kf,u,z,T)

% kf : the kalman fitler instance 

Q1=kf.Q;
Q2=kf.Qsise;
R=kf.R;
F=kf.F;
H=kf.H;
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
        x1=kf.x0;
        P1=kf.P0;
        x2=kf.x0;
        P2=kf.P0; 
        mu1=0.5;
        mu2=0.5;
    end
    % 1: input interaction: c=[mu1 mu2]*T
    c1=T(1,1)*mu1+T(2,1)*mu2;
    c2=T(1,2)*mu1+T(2,2)*mu2;
    % mumu=diag([mu1,mu1])*T   mu=mumu./sum(mumu); % colum normalize
    mu11=T(1,1)*mu1/c1;
    mu12=T(1,2)*mu1/c2;
    mu21=T(2,1)*mu2/c1;
    mu22=T(2,2)*mu2/c2;
    % 2: calculate the initial fusion state vector of model m
    % [x_pre1;x_pre2]=[x1, x2] [mu11 mu12; mu21; mu22]
    x_pre1=x1*mu11+x2*mu21;
    x_pre2=x1*mu12+x2*mu22;
    % 3: calculate the mixed covariance
    P_pre1=mu11*(P1+(x1-x_pre1)*(x1-x_pre1)')+mu21*(P2+(x2-x_pre1)*(x2-x_pre1)');
    P_pre2=mu12*(P1+(x1-x_pre2)*(x1-x_pre2)')+mu22*(P2+(x2-x_pre2)*(x2-x_pre2)');
    % 4: the kalman filter
    % update 1 
    x1_=F*x_pre1+Ga*u(:,i);
    P1_=F*P_pre1*F'+Q1;
    v1=z(:,i)-H*x1_;
    S1=H*P1_*H'+R;
    K1=P1_*H'/S1;
    x1= x1_+  K1 *v1; % posteriori
    % update 2
    x2_=F*x_pre2+Ga*u(:,i);
    P2_=F*P_pre2*F'+Q2;
    v2=z(:,i)-H*x2_;
    S2=H*P2_*H'+R;
    K2=P2_*H'/S2;
    x2= x2_+ K2*v2; % posteriori
    % 5: observe the likelihood
    Gamma1=1/sqrt(2*pi*det(S1))*exp(-0.5*v1'*inv(S1)*v1);
    Gamma2=1/sqrt(2*pi*det(S2))*exp(-0.5*v2'*inv(S2)*v2);
    % 6: calculate the filtering covariance matrix
    P1=(eye(n)-K1*H)*P1_*(eye(n)-K1*H)'+K1*R*K1';
    P2=(eye(n)-K2*H)*P2_*(eye(n)-K2*H)'+K2*R*K2';
    % 7: update the model probability
    c= Gamma1*c1+Gamma2*c2;
    mu1= Gamma1*c1/c;
    mu2= Gamma2*c2/c;
    % 8: output interaction
    x=mu1*x1+mu2*x2;
    % 9: interact the error covariance matrix
    P=mu1*(P1+(x1-x)*(x1-x)')+mu2*(P2+(x2-x)*(x2-x)');

    muv=[mu1;mu2];
    % store the data
    statef1_(:,i)=x1_';
    statef2_(:,i)=x2_';
    statef(:,i)=x';
    covf(:,:,i)=P;
    MU(:,i)=muv';
end
tcost=toc;

kff.tcost=tcost;
kff.statef1_=statef1_;
kff.statef2_=statef2_;
kff.statef=statef;
kff.covf=covf;
kff.MU=MU;


end