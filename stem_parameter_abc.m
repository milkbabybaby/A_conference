clear;
clc;

tf=1

syms t s

M=1
N=5
A_diag=sym('a',[1 N-1])



A_diag=sym(A_diag,'positive');




I=eye(N);



B=zeros(N,M);
B(1,1)=1
%B(4,2)=1
%A=[zeros(1,N);diag(A_diag),zeros(N-1,1)]; %stem
A=[zeros(1,N);diag(A_diag(1:N-1)),zeros(N-1,1)];%circle 
%A(2,4)=A_diag(N);
W=int(expm(-A*t)*B*B.'*expm(-A.'*t),t,0,s);%

H=inv(W);

result0=subs(trace(H),s,tf);
result=subs(result0,{A_diag},{A_diag.^0.5});

% ff=a+b+c+d+e+f+g
% 
% d_f=[diff(ff,a,1),diff(ff,b,1),diff(ff,c,1),diff(ff,d,1),diff(ff,e,1),diff(ff,f,1),diff(ff,g,1)]

d_f2=[];


for ii=1:N-1
d_f2=[d_f2,diff(result,A_diag(ii),1)];
end 

    
eq=[];
for i=1:N-2
    eq=[eq,d_f2(N-1)==d_f2(i)];
 end 

soll=struct2cell(solve(eq, sum(A_diag)==N));

sol = sqrt(double(cat(1,soll{:})))

% Function F=myfun(x)
% a=x(1)
% b=x(2)
% c=x(3)
% d=x(4);
% e=x(5);
% f=x(6);
% g=x(7);
% F=[d_f2(7)-d_f2(1),d_f2(7)-d_f2(2),d_f2(7)-d_f2(3),d_f2(7)-d_f2(4),d_f2(7)-d_f2(5),d_f2(7)-d_f2(6), a+b+c+d+e+f+g-8]
% 
% sol=fsolve('myfun',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]',optimset('Display','off'))