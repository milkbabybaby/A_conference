clear
clc
i=1
B=[ 1 0 0; 0 0 1]';
N=3
tf=1
x0=[];
y0=[];
for x=0.01:0.01:2
for y=0.01:0.01:2
     x0=[x0,x];
    y0=[y0,y];
    
A=[0 0 z; x 0 0; 0 y 0];

 H=zeros(N,N);

 ot=0.01; 
for k=1:tf/ot
    H=H+expm(-A*(ot*k))*B*B'*expm(-A'*(ot*k))*ot;
end
 C=pinv(H);
 z(i)=trace(C);
 i=i+1;

end 
end 
 
%z=trace(inv(double(int(expm(-A*t)*B*B'*expm(-A'*t),t,0,tf))));

% 
% x0=[];
% y0=[];
% for i=0.01:0.01:1
%     for j=0.01:0.01:1
%    
%     end
% end
z0=log10(z);
[X,Y,Z]=griddata(x0,y0,z0,linspace(0,2,1000)',linspace(0,2,1000),'linear');
%pcolor(X,Y,Z);shading interp
contourf(X,Y,Z)
xlabel('a_1')
ylabel('b_1')
colorbar