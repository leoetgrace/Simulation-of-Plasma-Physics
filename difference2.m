%%This function is used to solve the difference-differential equation given
%%a certain set of parameter. The result will be shown by eigenvalues and
%%eigenvactors.


function [A,B]=difference2(b,ita,step)
%parameter definition
%step=0.1;
%ita=-20:step:20;
q=1;
delta=1;
s=1;
tau=1;
epls=0.2;
Ti=1;

%difference matrix
l=length(ita);
dia=ones(1,l)*(-2);
D1=diag(dia,0);
dia=ones(1,l-1);
D2=diag(dia,1);
D3=diag(dia,-1);
D=D1+D2+D3;
D=D/step^2;

%denifition of Eigen Matrix & Identity Matrix
I=ones(1,l);
I=diag(I,0);
Z=zeros(l);
%definition of M0
M0=-D;
%definition of M1
v1=cos(ita)+s*ita.*sin(ita);
V1=diag(v1,0);
M1=-(tau*sqrt(epls)*D+2*q^2*b*V1);
%denifition of M2
v21=1+s^2*ita.^2;
V21=diag(v21,0);
v22=cos(ita)+s*ita.*sin(ita);
V22=diag(v22,0);
M2=-(q*q*b*b*V21+2*q*q*b*tau*sqrt(epls)*V22);
%definition of M3
v3=1+s^2*ita.^2;
V3=diag(v3,0);
M3=(q^2*b*I+q^2*b^2*V3*tau*epls^(1/2));
%definition of companion matrix as well as the primitive one
Er=[Z,I,Z;Z,Z,I;M0,M1,M2];
El=[I,Z,Z;Z,I,Z;Z,Z,M3];

%find the eigenvalue using QR
[Ql,Rl]=qr(El);
[Qr,Rr]=qr(Er);
E=inv(Rl)*inv(Ql)*(Qr*Rr);
[A,B]=eigs(E,3*l,'largestimag');
end


