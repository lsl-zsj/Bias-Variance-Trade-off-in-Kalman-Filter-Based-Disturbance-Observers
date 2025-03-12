function kff=mkc_dob_forward(kf,u,z)

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
sigma_p=kf.sigma_p;
sigma_r=kf.sigma_r;
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
    %
    bp = chol(P_,'lower') ;
    %br = chol(R,'lower') ;
    cnt=3;
    num=3;
    while(num>0)
        %  
        if(num==cnt)
          x_tlast=x_; 
        else  
          x_tlast=x_t; 
        end
        num=num-1;
        dp= bp\x_;
        z_=z(:,i);
        %dr= br\z_;
        wp= bp\x_tlast;
        %wr= inv(br)*H*x_tlast;
        ep=dp-wp;
        %er=dr-wr;
        %  P_ and R
        Cx=diag(exp(-ep.*ep./(2*sigma_p.*sigma_p)));
        %Cy=diag(exp(-er.*er./(2*sigma_r.*sigma_r)));
        for kk=1:n
            if(Cx(kk,kk)<0.0001)
                Cx(kk,kk)=0.0001;
            end   
        end
        % obtain the D 
        % Pe=bp/Cx*bp'-P_;
        % Dtemp=Pe(1);
        % P_1=bp/Cx*bp';
        % P_1(1,2:3)=P_1(1,2:3)+Dtemp*G';
        % P_1(2:3,1)=P_1(2:3,1)+G*Dtemp;
        P_1=bp/Cx*bp';
        R_1=R;
        K_1=P_1*H'/(H*P_1*H'+R_1);
        % obtain the new measurement
        %X_t=X_+K_1*(z(:,i)-H*X_); 
        x_t=x_+K_1*(z_-H*x_);
        xe(cnt-num,i)=norm(x_t-x_tlast)/(norm(x_tlast)+0.001);
        % stored data for inspectation
        if(num==cnt-1)
%        er_matrix(:,i)=er;
        end
        % if(num==cnt-2)
        % ep_matrix(:,i)=ep;
        % end
        if(xe(cnt-num,i)<0.01)
            break
        end
        threshold(cnt-num,i)=xe(cnt-num,i);
    end 

    % update
    x=x_t;
    P=(eye(n)-K_1*H)*P_*(eye(n)-K_1*H)'+K_1*R*K_1';
    % store the data
    statef_(:,i)=x_';
    statef(:,i)=x';
    covf_(:,:,i)=P_;
    covf(:,:,i)=P;
    covf_1(:,:,i)=P_1;
end
tcost=toc;

kff.tcost=tcost;
kff.statef_=statef_;
kff.statef=statef;
kff.covf_=covf_;
kff.covf=covf;
kff.covf_1=covf_1;
%kff.ep_matrix=ep_matrix;
%kff.er_matrix=er_matrix;
kff.threshold=threshold;

end