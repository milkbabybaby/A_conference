%%simple gradient

close all
clear all
 clc
format long 
%N=100


  
   
load('A_20_6_5_SF.mat');
load('Bb_20_6_5.mat')
B=Bb;
[N,M]=size(B)
A=A';

Topo=A;

A=(sqrt((N)/trace(A'*A)))*A; 

I=eye(N);



%B=I(:,1:M)

%  B=zeros(N,M);
%  B(1,1)=1
% B(3,2)=1

%B=I;

I1=eye(M);



tf=1









v=0.0025;


iterations=100
Cost1=zeros(iterations,1);
Cost2=zeros(iterations,1);
Cost3=zeros(iterations,1);
Cost4=zeros(iterations,1);  
Cost_cos=zeros(iterations,1);  
weight_changes=[];

AA=zeros(N,N);
for iiiii=1:1
 
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

%V2=roundn(V1,-10);

V1=esort(V1);
for iii=1:N
%     if real(V(iii))-real(V(iii+1))<0.02
%         V(iii)=real(V(iii+1))+i*imag(V(iii))
%     end 
if abs(imag(V1(iii)))<0.0001
    V1(iii)=real(V1(iii));
end 
if abs(V1(iii))<0.05
    V1(iii)=0;
end
end 


V=esort(V1);

for iii=1:N-1
    if imag(V(iii))==0 && imag(V(iii+1))==0 && V(iii)-(V(iii+1))<0.05
        V(iii+1)=V(iii);
    end 
end 


   
V=esort(V);
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

%BBB=sym(BBB);
%invBBB=pinv(BBB);
%invBBB=double(invBBB);

%alpha_tf=inv(BB)*CC1;
alpha_tf=BBB\CCC;

%alpha_tf =TVreg(BBB,CCC2,100)

%alpha_tf=invBBB*CCC

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
%alpha_t=invBBB*CCC2;

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

Cost2(ii)=F(:)'*F3(:)/(norm(F(:))*norm(F3(:)));
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

AA=AA+A;
end 


figure(1)
  plot((Cost1),'r-*')
  
  

  
  
  figure(2)
plot(log10(Cost1),'r-*','MarkerSize',5)
   
    set(gca, 'LineWidth', 1.5);
      set(gca,'FontName','Times New Roman','FontWeight','bold')
 xlabel('Iterations','FontName','Times New Roman','FontWeight','bold');
 ylabel('lg E(A) ','FontName','Times New Roman','FontWeight','bold');
  export_fig stem_converge.eps -painters -transparent 
  
    figure(3)
  plot(Cost2,'r-*','MarkerSize',5)
   
    set(gca, 'LineWidth', 1.5);
      set(gca,'FontName','Times New Roman','FontWeight','bold')
 xlabel('Iterations','FontName','Times New Roman','FontWeight','bold');
 ylabel('cos ','FontName','Times New Roman','FontWeight','bold');
  export_fig stem_cos.eps -painters -transparent 
  
  

%   figure(5)
%   plot(weight_changes(:,iterations),'r-*');
%   
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
%     hold on;
%    plot(weight_changes(6,:),'c:','MarkerSize',2,'LineWidth',2)
%  ylim([0 1.6]) 
%   legend('A_{12}','A_{23}','A_{34}','A_{45}','A_{56}','A_{61}','Location','East');   
%   set(gca, 'LineWidth', 1.5);
%       set(gca,'FontName','Times New Roman','FontWeight','bold')
%  xlabel('Iterations','FontName','Times New Roman','FontWeight','bold');
%  ylabel('Edge weight ','FontName','Times New Roman','FontWeight','bold');
%   export_fig weight_changes_14.eps -painters -transparent 
  
  
  
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

  
AAA=AA/max(max(AA));  
  
  
  
  
 k_out= sum(Topo,1)';
 k_in=sum(Topo,2);
 k_total=k_in+k_out;
 
 w=[];
 k_outmax=[];
 k_outmin=[];
  k_inmax=[];
 k_inmin=[];
  k_totalmax=[];
 k_totalmin=[];
 k_product=[];
for i=1:N
  for   j=1:N
    if AAA(i,j)~=0
        w=[w,AAA(i,j)]
        k_outmax=[k_outmax,max(k_out(i),k_out(j))];
        k_outmin=[k_outmin,min(k_out(i),k_out(j))];
         k_inmax=[k_inmax,max(k_in(i),k_in(j))];
        k_inmin=[k_inmin,min(k_in(i),k_in(j))];
         k_totalmax=[k_totalmax,max(k_total(i),k_total(j))];
        k_totalmin=[k_totalmin,min(k_total(i),k_total(j))];
        k_product=[k_product,k_out(j)*k_in(i)];
     end 
   end 
end 
 w_n=w/max(w);   
figure(1)
  plot(k_outmax,w_n,'r*')
  hold on
   plot(k_outmin,w_n,'ro')
   hold on
   plot(k_inmax,w_n,'c*')
      hold on
   plot(k_inmin,w_n,'co')
   hold on
   plot(k_totalmax,w_n,'b*')
   hold on
   plot(k_totalmin,w_n,'bo')
         hold on
   plot(k_product,w_n,'rx')
   legend('k_outmax','koutmin','k_inmax','k_inmin','k_totalmax','k_totalmin','k_product')
   plot(w)
 
w_out=sum(AAA,1)';
w_in=sum(AAA,2);
w_nodes=w_out+w_in;
plot(k_total,w_nodes,'r')
  