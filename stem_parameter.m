clear;
clc;
N=5

tf=1

A=[zeros(1,N);eye(N-1),zeros(N-1,1)]; 
I=eye(N);



B=I(:,1)
%B=I
syms t s


W=int(expm(-A*t)*B*B'*expm(-A'*t),t,0,s);

H=inv(W)

c=subs(trace(H),s,tf)

[v e]=eig(subs(W,s,tf))