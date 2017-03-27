clear all;
close all;
clc;


%{
%%%%%%%%%%%%%%%real-life网络直接载入Edge List%%%%%%%%%%%%%
dataset_name = 'Maspalomas'
load_path = strcat('./dataset/',dataset_name,'.mat');
load(load_path);
L = size(A,1)%L为边的个数
N2 = length(unique(A))
[A,~] = Net_Consecutive(double(A)); %调用汤沛的代码防止A编号不连续，反馈的A为邻接矩阵
%}


%%%%%%%%%%%%%%%ER-SF网络直接载入邻接矩阵%%%%%%%%%%%%%%%%%%%
A=ErSfNetGen(2,100,3);
dataset_name = 'A100ER'
L=sum(sum(A))


%{
%产生stem的邻接矩阵
A=[0 0 0 0 0 1 0 1 0 0;
   0 0 0 0 0 0 0 0 0 0;
   0 1 0 0 1 1 0 1 0 0;
   1 0 0 0 0 0 0 0 0 0;
   0 1 0 0 0 0 0 0 0 0;
   0 0 0 0 0 0 0 0 1 0;
   0 1 0 0 0 0 0 1 0 0;
   0 0 0 0 0 0 1 0 0 0;
   0 0 0 0 1 1 0 1 0 0;
   0 0 0 1 0 0 0 1 1 0;];
%}
N=size(A,1)
M=50;

sum1=zeros(N,M);
v=0.05; %迭代步长
ot=0.025; %积分时间步
iterations=1;  %总的循环次数
Iterations_number=140; 
tf = 1;



Cost11=zeros(N,M);
for ii=1:iterations
    
    B=randn(N,M); %B为控制权重矩阵，随机初始值
    MM=M; %规定了一个norm值
    
    B=(sqrt(1*(MM+0.1)/trace(B'*B)))*B; %initiate the norm to MM-1

    M=min(size((B)));  %M这里一般情况下就是M
    WB0=zeros(N,N); 
    F0=zeros(N,M);
    Xf=expm(A*tf)*expm(A'*tf); %对矩阵求指数，非每个元素求指数
    
    Cost=zeros(Iterations_number,1);
    Cost3=zeros(Iterations_number,1);
    Cost5=zeros(Iterations_number,1);

    
    %1200次梯度下降迭代
    for k=1:Iterations_number
        fprintf('No. %d Epoch, No. %d Descent \n',ii,k);
        
        WB0=zeros(N,N);
        F0=zeros(N,M);
        %先求出文章中的WB，即为这里的WB0
        for k1=1:tf/ot
            WB0 = WB0+expm(A*(ot*k1))*B*B'*expm(A'*(ot*k1))*ot;
        end
        C = pinv(WB0);
        
        Cost(k)=trace(C*Xf);
        Cost3(k)=trace(B'*B);
        
        %再求出文章中的能量函数对B矩阵的梯度。
        for k1=1:tf/ot
            F0=F0-(expm(A'*(ot*k1))*C*Xf*C*expm(A*(ot*k1)))*B*ot;
        end
        F00=F0/norm(F0); 
        F22=F0; %F22为本次迭代的梯度
        
        %F1为Norm函数对B的梯度
        F1=(1/ot)*(2*trace(B'*B)-2*MM)*2*B; % gradient of norm function
        F11=F1/norm(F1); %归一化Norm函数梯度
        F11=-(sqrt(1*(MM+0.001)/trace(F11'*F11)))*F11;    %Normailize gradient of norm function
        
        Cost5(k)=-sum(F00(:).*F11(:))/(norm(F00(:))*norm(F11(:)));
        
        %如果cost5接近-1便跳出梯度下降迭代循环
        if  Cost5(k)+1<0.000005 
            break;
        end
        
        if  Cost3(k)>=MM  %MM is the fixed energy bound
            %F0=reshape((eye(N*M)-F1(:)*pinv(F1(:)))*F0(:),N,M);  %the direction of gradient energy function  (projection)
            
            %这一部分完全为了节省空间代替上面这句代码
            F = F1(:); %F为列向量NM*1
            P = F'*F;
            P = 1./P;
            Fp = P*F'; %Fp为行向量1*NM
            
            FF = F0(:); %FF为列向量NM*1
            Fv = zeros(N*M,1);
            C = Fp*FF;
            Fv = FF - F.*C;
            F0 = reshape(Fv,N,M);%the direction of gradient energy function  (projection)
     
            F0=F0/norm(F0); %unit vector of gradient energy function         
            if sum(sum(F0.*F22))<=0
                F0=-F0;
            end
  
            B = B-v*F0/norm(F0);  % F0 is the projected direction
            B = sqrt(1*(MM+0.001)/trace(B'*B))*B; %initiate the norm to MM+1
        else
            B = B-v*F0/norm(F0);
        end
   
    end
    
end

%输出能量的log图
figure(1)
plot(log10(Cost),'r-*')

figure(2)
plot(Cost3,'r-*')
legend('trace BB')

figure(3)
plot(Cost5,'r-*')


Init = Cost(1)
PGM = Cost(end)


%下面求PGME
%按照Science paper的方法求Bs
[N,M]=size(B);
Bb = zeros(N,M);
sum_B=sum(abs(B),2);
[max_v order]=sort(sum_B);%max_v为值，order为在sum_B中的序号
for k=1:M
    Bb(order(end-k+1),k)=1;
end

WB0=zeros(N,N);
F0=zeros(N,M);

%先求出文章中的WB，即为这里的WB0
for k1=1:tf/ot
    WB0=WB0+expm(A*(ot*k1))*Bb*Bb'*expm(A'*(ot*k1))*ot;
end
C=pinv(WB0);
PGME=trace(C*Xf)

save_path = strcat('PGM-M',num2str(M),'tf1_',dataset_name);
save(save_path,'Cost','B','Bb','Init','PGM', 'PGME','Iterations_number');
%save(save_path,'Cost','B','Bs1','Bs2','Init','PGM', 'PGME1','PGME2','Iterations_number');
