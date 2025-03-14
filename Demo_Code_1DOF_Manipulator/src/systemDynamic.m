function [X_,y]=systemDynamic(m,k,b,T,X,u,d)

% X is the state vector 

% X(1) is the external disturbance 
% X(2) is the angular velocity
% X(3) is the angle 

% X_(1)=X(1)+randn(1,1);
X_(1)=d+0.1*randn(1,1);
X_(2)=T/m*X(1)+(1-b*T/m)*X(2)+(-k*T/m)*X(3)+T/m*u+0.01*randn(1,1);
X_(3)=T*X(2)+X(3)+0.01*randn(1,1);

y=X_(3)+0.01*randn(1,1);

end