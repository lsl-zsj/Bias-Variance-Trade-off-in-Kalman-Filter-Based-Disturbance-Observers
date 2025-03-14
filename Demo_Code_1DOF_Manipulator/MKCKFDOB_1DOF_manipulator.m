function MKCKFDOB_1DOF_manipulator()

% Brief: One line description of what the function or class performs
% Details:
%    None
% 
% Syntax:  
%     MKCKFDOB_1DOF_manipulator()
% 
% Inputs:
%    None
% 
% Outputs:
%    None
% 
% Example: 
%    None
% 
% See also: None

% Author:                          ShileiLi
% Email:                           slidk@connect.ust.hk
% Created:                         14-Mar-2025 10:52:08
% Version history revision notes:
%                                  None
% Implementation In Matlab R2023a
% Copyright © 2025 ShileiLi.All Rights Reserved.
%

clear all
addpath(genpath(pwd));  %
%
dt=0.01;
m=0.1; k=0.1; b=1; T=0.01;
% kd-dob 
F=[1,0,0; T/m, 1- b/m*T, -k/m*T; 0, T, 1];
G=[0;T/m;0];              % input matrix
H=[0,0,1];     % observation matrix 

x0=[0;0;0];        % initial guess
P0=eye(3);             % initial covariance

tsim = 10;
len=tsim/dt;            % step number
t    = 0:dt:tsim-dt;          % time

Vwd=0.01;
Vwn=0.0001;
Vvn=0.0001;

Qx=eye(2);
Q=diag([Vwd;Vwn*diag(Qx)]);
R=Vvn*eye(1);                % measurement covariance

%% start simulation
% command angle
theta=15*sin(2*pi*0.2*t);
theta_d=15*2*pi*0.2*cos(2*pi*0.2*t);
theta_dd=-15*2*pi*0.2*2*pi*0.2*sin(2*pi*0.2*t);
theta=theta';theta_d=theta_d';theta_dd=theta_dd';

% initial state
X=[0;0;0];
X_last=[0;0;0];
S=zeros(len,3);
U=zeros(len,5);
X_R=zeros(len,3);
Y=zeros(len,1);

d=zeros(len,1);
d(400:600)=50;

% real system
X_r=[0;0;0];
u=0;
rng(42)

sigma_p=[1.5 100000000 100000000]';

tic
for i=1:len
    %% real system 
    X_rlast=X_r;
    %% generate state
    [X_r,y]=systemDynamic(m,k,b,T,X_rlast,u,d(i)); 
    %% estimate
    if(i==1)
        X=x0;
        P=P0;
    end
    % state prediction 
    % prediction 
    X_=F*X+G*u;
    % prediction covariance
    P_=F*P*F'+Q;
    %%
    bp = chol(P_,'lower') ;
    cnt=5;
    num=5;
    x_=X_;
    while(num>0)
        %  
        if(num==cnt)
          x_tlast=X_; 
        else  
          x_tlast=x_t; 
        end
        num=num-1;
        dp= bp\x_;
        wp= bp\x_tlast;
        ep=dp-wp;
        %  P_ and R
        Cx=diag(exp(-ep.*ep./(2*sigma_p.*sigma_p)));
        for kk=1:3
            if(Cx(kk,kk)<0.0001)
                Cx(kk,kk)=0.0001;
            end   
        end
        P_1=bp/Cx*bp';
        R_1=R;
        K_1=P_1*H'/(H*P_1*H'+R_1);
        x_t=x_+K_1*(y-H*x_);
        xe(cnt-num,i)=norm(x_t-x_tlast)/(norm(x_tlast)+0.001);
        % stored data for inspectation
        if(xe(cnt-num,i)<0.01)
            break
        end
        threshold(cnt-num,i)=xe(cnt-num,i);
    end 
        %P=(eye(3)-K_1*H)*P_1;
        temp=eye(3)-K_1*H;
        P= temp*P_*temp'+K_1*R*K_1';
        X=x_t;
        S(i,:)=X';
    
    %% control
    u1=m*theta_dd(i)+b*theta_d(i)+k*theta(i);
    u2=-X(1);
    kp=100;kd=10;
    u3=kp*(theta(i)-X(3));
    u4=kd*(theta_d(i)-X(2));
    u=u1+u2+u3+u4;
    
    %% data
     S(i,:)=X';
     U(i,:)=[u1,u2,u3,u4,u];
     X_R(i,:)=X_r;
     Y(i)=y;
end
toc
error.d=S(:,1)-X_R(:,1);
error.v=S(:,2)-X_R(:,2);
error.p=S(:,3)-X_R(:,3);
error.p_d=theta;
error.p_a=X_R(:,3);
error.p_err=theta-X_R(:,3);
error.u=U(:,5);

figure
subplot(2,1,1)
plot(t,theta,t,S(:,3))
subplot(2,1,2)
plot(t,theta_d,t,S(:,2))

save MKC_error.mat error

figure
subplot(3,1,1)
plot(t,S(:,1),t,X_R(:,1))
hold on
plot(t,S(:,1)-X_R(:,1))
subplot(3,1,2)
plot(t,S(:,2),t,X_R(:,2))
hold on
plot(t,S(:,2)-X_R(:,2))
subplot(3,1,3)
plot(t,S(:,3),t,X_R(:,3))
hold on
plot(t,S(:,3)-X_R(:,3))


end