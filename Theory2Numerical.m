%This program is used to compare the theoretical solution to numerical
%solution to verify the accurate of theoretical solution

%Parameter define
b=1;
step=0.1;
ita=-20:step:20;
q=1;
delta=1;
s=1;
tao=1;
epls=0.2;
Ti=1;
%
l=length(ita);
R1=zeros(1,3*l);
R2=zeros(1,3*l);
[A B]=difference2(b,ita,step);
for i=1:3*l
   R1(1,i)=real(B(i,i))*sqrt(0.5*b/sqrt(epls));%UNITY
   R2(1,i)=imag(B(i,i))*sqrt(0.5*b/sqrt(epls));
end
plot(R1,R2,'o')
xlabel('real(¦¸)','FontSize',12);
ylabel('imag(¦¸)','FontSize',12);
%title('The struture of \Omega when b=?','FontSize',12); %This is a plot of
%numerical simulation
hold on


%%%Next we solve theoretical solution numerically
syms omega
eqn=(tao*omega/(1+tao*omega*sqrt(epls))+b+2/omega)^2+1/q^2/omega^2/b*(b*s^2+(2*s-1)/omega)==0;
sol=solve(eqn,omega);
sol=double(sol);
m=max(imag(sol));
[x0,y0]=find(imag(sol)==m);
plot(real(sol(5,1))*sqrt(0.5*b/sqrt(epls)),imag(sol(5,1))*sqrt(0.5*b/sqrt(epls)),'*') %This is the solution of theoretical result

