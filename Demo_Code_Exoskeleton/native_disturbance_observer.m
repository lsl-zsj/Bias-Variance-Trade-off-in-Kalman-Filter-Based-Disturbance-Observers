function D=native_disturbance_observer(q1,dq1,ddq1,cur1)

load('optimresultsLeg.mat')
para=optimresultsLeg.x;
l1=0.42;
tau_c2=para(1);
tau_s2=para(2);
dqs2=para(3);
eta2=para(4);
X2=para(5);
Y2=para(6);
J2=para(7);

tau_c1=para(8);
tau_s1=para(9);
dqs1=para(10);
eta1=para(11);
X1=para(12);
Y1=para(13);
J1=para(14);


len=length(q1);
D=zeros(len,2);
for i=1:len
q=q1(i,:);
dq=dq1(i,:);
ddq=ddq1(i,:);
cur=cur1(i,:);
%% M matrix
M11=J1+2*l1*(X2*cos(q(2))-Y2*sin(q(2)));
M12=J2+l1*(X2*cos(q(2))-Y2*sin(q(2)));
M21=J2+l1*(X2*cos(q(2))-Y2*sin(q(2)));
M22=J2;
M=[M11,M12;M21,M22];
%% C matrix
C11=-2*l1*(X2*sin(q(2))+Y2*cos(q(2))).*dq(2);
C12=-l1*(X2*sin(q(2))+Y2*cos(q(2))).*dq(2);
C21=l1*(X2*sin(q(2))+Y2*cos(q(2))).*dq(1);
C22=0;
C=[C11,C12;C21,C22];
%% G matrix
G1=9.7925*(X1*sin(q(1))+Y1*cos(q(1)+X2*sin(q(1)+q(2))+Y2*cos(q(1)+q(2))));
G2=9.7925*(X2*sin(q(1)+q(2))+Y2*cos(q(1)+q(2)));
G=[G1;G2];
%% F matrix：fraction
% F1=(tau_c1+(tau_s1-tau_c1)*exp(-dq(1)/dqs1)).*customSIGN(dq(1))+eta1*dq(1);
% F2=(tau_c2+(tau_s2-tau_c2)*exp(-dq(2)/dqs2)).*customSIGN(dq(2))+eta2*dq(2);
% F=[F1;F2];
% F=[0;0];
%% native disturbance
dist= -cur' + (M*ddq'+C*dq'+G);

D(i,:)=dist';
end

end



function y=customSIGN(u)
threshold=0.05;
if(u>threshold)
    y=1;
elseif(u<-threshold)
    y=-1;
else
    y=u/threshold;
end
end



