function IMMKFDOB_1DOF_manipulator()

% Brief: One line description of what the function or class performs
% Details:
%    None
% 
% Syntax:  
%     IMMKFDOB_1DOF_manipulator()
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
% Created:                         14-Mar-2025 11:06:05
% Version history revision notes:
%                                  None
% Implementation In Matlab R2023a
% Copyright © 2025 ShileiLi.All Rights Reserved.
%


clear all
addpath(genpath(pwd));  % 将当前目录加入路径

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

tic
for i=1:len
    %% real system 
    X_rlast=X_r;
    %% generate state
    [X_r,y]=systemDynamic(m,k,b,dt,X_rlast,u,d(i)); 
    %% estimate
    %
    Q1=diag([Vwd*exp(0);Vwn*diag(Qx)]);
    Q2=diag([Vwd*exp(4);Vwn*diag(Qx)]);
    T=[0.98,0.02;0.5,0.5];
    n=3;
    %% IMMKFDOB
    if(i==1)
        x1=x0;
        P1=P0;
        x2=x0;
        P2=P0; 
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
    x1_=F*x_pre1+G*u;
    P1_=F*P_pre1*F'+Q1;
    v1=y-H*x1_;
    S1=H*P1_*H'+R;
    K1=P1_*H'/S1;
    x1= x1_+  K1 *v1; % posteriori
    % update 2
    x2_=F*x_pre2+G*u;
    P2_=F*P_pre2*F'+Q2;
    v2=y-H*x2_;
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
    
    X=x;
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

 save IMM_error.mat error

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

figure
hold on
plot(t,MU(1,:))
plot(t,MU(2,:))


end