clear;
clc;
 load('D:\code\A_conference\A_20_4_4_SF.mat');
 load('D:\code\A_conference\Bb_20_4_4.mat');
 A=A';
% N=3
% A=[zeros(1,N);diag(ones(N-1,1)),zeros(N-1,1)]
tf=1
A_topo=sym(A);
N=length(A_topo);
syms t s;

A0=sym('A_%d_%d',[N N]);
A0=sym(A0,'positive');

A=A0.*A_topo;
Index = find(A);

A_var=A(Index);

I=eye(N);



B=I(:,1);

trace_AA=trace(A*A.');

W=int(expm(-A*t)*B*B.'*expm(-A.'*t),t,0,s);

H=inv(W);

result0=subs(trace(H),s,tf);
result=subs(result0,{A_var},{A_var.^0.5});

% ff=a+b+c+d+e+f+g
% 
% d_f=[diff(ff,a,1),diff(ff,b,1),diff(ff,c,1),diff(ff,d,1),diff(ff,e,1),diff(ff,f,1),diff(ff,g,1)]

d_f2=[];


for ii=1:N-1
d_f2=[d_f2,diff(result,A_var(ii),1)];
end 

    
eq=[];
for i=1:N-2
    eq=[eq,d_f2(N-1)==d_f2(i)];
 end 

soll=struct2cell(solve(eq, sum(A_var)==N));

sol = sqrt(double(cat(1,soll{:})))