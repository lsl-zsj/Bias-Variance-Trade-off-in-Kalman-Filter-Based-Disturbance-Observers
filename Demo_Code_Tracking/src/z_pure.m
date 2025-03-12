function [state,z]=z_pure(kf,len)
%% state and measurement generation
Q=kf.Q;
R=kf.R;
F=kf.F;
H=kf.H;
x0=kf.x0;
G=kf.G;

for i = 1:len
    %
    u_n(:,i) = randn(1); % this actually is the sequence of 0 to N-1
    w_xn= randn(2,1);
    % real state generation
    if(i==1)
    x(:,i) =F*x0+G*u_n(:,i)+w_xn;
    else
    %x(:,i) = F * x(:,i-1) + G*u(:,i)+  G* u_n(:,i);   % Generate truth  u(:,i) actually is u(:,i-1)
    x(:,i) = F * x(:,i-1) + G*u_n(:,i)+w_xn;   % Generate truth  u(:,i) actually is u(:,i-1)
    end
    % measurement
    z(:,i) = H*x(:,i) + sqrt(R) * randn(2,1);      
end
state=x; % state is the ground truth state
end