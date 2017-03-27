%%simple gradient

close all
clear all
 clc
%format long 
%N=100

  
 


%A=[0,0,0,0,0,0,0,0,0,0;1,0,0,0,0,0,0,0,0,0;0,1,0,0,0,0,0,0,0,0;0,0,1,0,0,0,0,0,0,0;1,0,0,0,0,0,0,0,0,0;0,0,0,0,1,0,0,0,0,0;0,0,0,0,0,1,0,0,0,0;0,0,0,0,0,1,0,0,0,1;0,0,0,0,0,0,0,1,0,0;0,0,0,0,0,0,0,0,1,0;]
%A=[0,0,0,0,0,0,0,0;1,0,0,0,0,0,0,0;0,1,0,0,0,0,0,0;0,0,1,0,0,0,0,0;0,0,0,1,0,0,0,0;0,1,0,0,0,0,0,1;0,0,0,0,0,1,0,0;0,0,0,0,0,0,1,0;];
%A=[zeros(1,N);diag(ones(1,N-1)),zeros(N-1,1)];
% A(1,6)=1;
%  load('D:\code\NJP\NJP\dataset\circuit-s208.mat');
%  [A,~] = Net_Consecutive(double(A));
load('D:\code\dataset\Foodweb\Rhode.mat')
[A,~] = Net_Consecutive(double(A));
Topo=A;

 load('D:\code\A_conference\Bb.mat')
B=Bb;
%%random  initilization
%A=rand(N,N);
%A=Topo.*A;

 N=size(A,1);  
 M=size(B,2);
A=(sqrt((N)/trace(A'*A)))*A; 

I=eye(N);



%B=I(:,1:M)

 %B=zeros(N,M);
 %B(1,1)=1
%B(3,2)=1

%B(8,3)=1
%B=I;

I1=eye(M);



tf=1









v=0.0004;


iterations=500
Cost1=zeros(iterations,1);
% Cost2=zeros(iterations,1);
% Cost3=zeros(iterations,1);
% Cost4=zeros(iterations,1);  
Cost_cos=zeros(iterations,1);  
weight_changes=[];


 
for ii=1:iterations
   
 
 H=zeros(N,N);

 ot=0.01; 
for k=1:tf/ot
    H=H+expm(-A*(ot*k))*B*B'*expm(-A'*(ot*k))*ot;
end
 C=pinv(H); 

   
Cost1(ii)=trace(C);

AA(:,:,ii)=A;
   
V1=eig(A);

V=esort(V1);

 II=eye(N);
for  iii=1:N-1
    
    if V(iii)==conj(V(iii+1)) & V(iii)~=V(iii+1)
        
        II(iii,iii)=0.5;
        
        II(iii,iii+1)=0.5;
           II(iii+1,iii)=0.5*i;
        
        II(iii+1,iii+1)=-0.5*i;
        
  
        
    end 
end




 chongshu=zeros(N,1);
ik=1;
while  ik>=1 && ik <=N
    index=find(V(ik:N)==V(ik));
    chongshu(ik)=length(index);
    ik=ik+length(index);
end 
 
 


BB=fliplr(vander(V));
CC1=exp(-V*tf);
for ik=1:N-1
    if chongshu(ik)>1
syms s
f=[];
for ppp=0:N-1
    f=[f,s^ppp];
end 
    

ff=f;
for pp=1:chongshu(ik)-1
    ff=[ff;diff(f,'s',pp)]; %%%求导
    CC1(ik+pp)=(-tf)^pp*exp(-V(ik)*tf);%%%求导
end 

fff=subs(ff,s,V(ik));      
BB(ik:ik+chongshu(ik)-1,:)=fff;%%%
 end 
end 



    

BBB=II*BB;
CCC=II*CC1;    


%alpha_tf=inv(BB)*CC1;
alpha_tf=BBB\CCC;

%alpha_tf =TVreg(BBB,CCC2,100)



F111=zeros(N,N);


for kk=1:tf/ot
    
   CC2=exp(-V*ot*kk); 
    
     for ik=1:N-1
    if chongshu(ik)>1


for pp=1:chongshu(ik)-1
     CC2(ik+pp)=(-ot*kk)^pp*exp(-V(ik)*(ot*kk));
end
    end
     end


CCC2=II*CC2;    

alpha_t=BBB\CCC2;


 pf=expm(-A*ot*kk); 

CB=C^2*pf*B*B';


F11=zeros(N,N);
for iiii =1:N-1
    
F1=zeros(N,N);
    for k=1:iiii
        
        F1=A'^(k-1)*CB*A'^(iiii-k)+F1;

    end

 F11=-2*alpha_t (iiii+1)*F1+F11;
end

F111=F11*ot+F111;


  
end 

    
    
  
      
 F=Topo.*F111/norm(Topo.*F111);
 
 F3=2*A;%(trace(A'*A)-N)*4*A; % gradient of norm function
 
 F3=F3/norm(F3);


F0=reshape((eye(N*N)-F3(:)*pinv(F3(:)))*F(:),N,N); 

Cost_cos(ii)=F(:)'*F3(:)/(norm(F(:))*norm(F3(:)));
A=A-v*F0/norm(F0);
%   F=F'+F-diag(diag(F));
%  for kkk=1:N
%     F(kkk,kkk)=0;
% end

%A=A-v*Topo.*F111/norm(F111);
A=(sqrt((N)/trace(A'*A)))*A; 



     weight=diag(A(2:N,1:N-1));
     weight(N,1)=A(1,N);
weight_changes=[weight_changes,weight];

ii  

 
end




figure(1)
  plot((Cost1),'r-*')
  
  

  
  figure(3)
plot(log10(Cost1),'r-*','MarkerSize',5)
   xlabel('Iterations','FontName','Times New Roman','FontWeight','bold');
 ylabel('lg E(A) ','FontName','Times New Roman','FontWeight','bold');
    set(gca, 'LineWidth', 1.5);
      set(gca,'FontName','Times New Roman','FontWeight','bold')
  export_fig stem_converge.eps -painters -transparent 
  
    figure(4)
  plot(Cost_cos,'r-*','MarkerSize',5)
 xlabel('Iterations','FontName','Times New Roman','FontWeight','bold');
 ylabel('cos(\theta)','FontName','Times New Roman','FontWeight','bold');
    set(gca, 'LineWidth', 1.5);
      set(gca,'FontName','Times New Roman','FontWeight','bold')
  export_fig stem_cos.eps -painters -transparent 
  
  

  figure(5)
  plot(weight_changes(:,iterations),'r-*');
  
%   figure (6);
%   plot(weight_changes(1,:),'r-*','MarkerSize',2,'LineWidth', 2)
%   hold on;
%   plot(weight_changes(2,:),'b--','MarkerSize',2,'LineWidth', 2)
%   hold on;
%   plot(weight_changes(3,:),'g-','MarkerSize',2,'LineWidth', 2)
%   hold on;
%   plot(weight_changes(4,:),'m-.','MarkerSize',2,'LineWidth', 2)
%   hold on;
%    plot(weight_changes(5,:),'k:','MarkerSize',2,'LineWidth',2)
%   hold on;
%    plot(weight_changes(6,:),'c:','MarkerSize',2,'LineWidth',2)
%  ylim([0 1.6]) 
%   legend('A_{21}','A_{32}','A_{43}','A_{54}','A_{65}','A_{16}','Location','East');   
%   set(gca, 'LineWidth', 1.5);
%   xlabel('Iterations','FontName','Times New Roman','FontWeight','bold');
%  ylabel('Interaction strength ','FontName','Times New Roman','FontWeight','bold');
%    set(gca,'FontName','Times New Roman','FontWeight','bold')
%   export_fig weight_changes_13_circle.eps -painters -transparent 
  
  figure (6);
  plot(weight_changes(1,:),'r-*','MarkerSize',2,'LineWidth', 2)
  hold on;
  plot(weight_changes(2,:),'b--','MarkerSize',2,'LineWidth', 2)
  hold on;
  plot(weight_changes(3,:),'g-','MarkerSize',2,'LineWidth', 2)
  hold on;
  plot(weight_changes(4,:),'m-.','MarkerSize',2,'LineWidth', 2)
  hold on;
   plot(weight_changes(5,:),'k:','MarkerSize',2,'LineWidth',2)

 ylim([0 1.6]) 
  legend('A_{21}','A_{32}','A_{43}','A_{54}','A_{65}','Location','East');   
  set(gca, 'LineWidth', 1.5);
  xlabel('Iterations','FontName','Times New Roman','FontWeight','bold');
 ylabel('Interaction strength ','FontName','Times New Roman','FontWeight','bold');
   set(gca,'FontName','Times New Roman','FontWeight','bold')
  export_fig weight_changes_13.eps -painters -transparent 


  
    AL=zeros(N,N);
for kkkk=1:N
    
   AL= alpha_tf(kkkk)*A^(kkkk-1)+AL;
end

AAAA=expm(-A*tf);
AAAA-AL


 H1=zeros(N,N);
 
for k=1:tf/ot
    H1=H1+expm(-A*(ot*k))*B(:,1)*B(:,1)'*expm(-A'*(ot*k))*ot;
end
 C1=pinv(H1); 

  trace(C1)
  
   H2=zeros(N,N);
 
for k=1:tf/ot
    H2=H2+expm(-A*(ot*k))*B(:,2)*B(:,2)'*expm(-A'*(ot*k))*ot;
end
 C2=pinv(H2); 

  trace(C2)

  