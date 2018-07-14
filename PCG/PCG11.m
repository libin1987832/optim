% N=10;
% p=rand([3,N])*100;
% p=[aa(:,3:5);ones(1,N)];
% 
% yaw_i = 0.5854; pitch_i = 0.3; roll_i = 0.2; tran_i=[1,3,4];
% dcm_i = angle2dcm( yaw_i, pitch_i, roll_i );
% T_i =  [dcm_i,tran_i';0 0 0 1];
% T_i_0 = [angle2dcm( yaw_i+0.01, pitch_i+0.02, roll_i+0.01 ),tran_i'+[0.01;0.02;0.01];0 0 0 1];
% T_i_inv=[dcm_i',-dcm_i'*tran_i';0 0 0 1];
% p_i = T_i_inv*p+[rand([3,N])*0.01;zeros(1,N)];
% 
% yaw_j = 0.4854; pitch_j = 0.5; roll_j = 0.41; tran_j=[1,4,6];
% dcm_j = angle2dcm( yaw_j, pitch_j, roll_j );
% T_j = [dcm_j,tran_j';0 0 0 1];
% T_j_0 = [angle2dcm( yaw_j+0.01, pitch_j+0.02, roll_j+0.01 ),tran_j'+[0.01;0.02;0.01];0 0 0 1];
% T_j_inv=[dcm_j',-dcm_j'*tran_j';0 0 0 1];
% p_j = T_j_inv*p+[rand([3,N])*0.01;zeros(1,N)];
% 
% 
% yaw_k = 0.3854; pitch_k = 0.1; roll_k = 0.3; tran_k=[1,40,60];
% dcm_k = angle2dcm( yaw_k, pitch_k, roll_k );
% T_k = [dcm_k,tran_k';0 0 0 1];
% T_k_0 = [angle2dcm( yaw_k+0.01, pitch_k+0.02, roll_k+0.01 ),tran_k'+[0.01;0.02;0.01];0 0 0 1];
% T_k_inv=[dcm_k',-dcm_k'*tran_k';0 0 0 1];
% p_k = T_k_inv*p+[rand([3,N])*0.01;zeros(1,N)];
load data.mat

% load test1_3.mat
% load test2_3.mat
% load tt_3.mat
% data=[ss1;ss2];
% Transform=[
%     [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1];
% %     [0.796264498872805,0.527854767232124,-0.295520206661340,1;
% %     -0.492584019961313,0.849316453694221,0.189796060978687,3;
% %     0.351174929506071,-0.00555933400618633,0.936293363584199,4;
% %     0,0,0,1];
% % [0.785895001742623,0.532369162497330,-0.314566560616118,1.01000000000000;-0.494228099222101,0.846512299346681,0.197877520183836,3.02000000000000;0.371628352222044,-4.33207375641065e-05,0.928381584235729,4.01000000000000;0,0,0,1];
% %  [0.787126063354628,0.530135956919702,-0.315259302112845,1.39136040130988;-0.492383253214346,0.847912457216081,0.196476962649172,2.85079478764810;0.371471792121425,0.000576262630489977,0.928444061631867,4.12646486264528;0,0,0,1];
% %  tt1(5:12,:)tt1(5:12,:)
% % tt1(1:8,:);
% % [0.922018753951339,0.374051742787868,-0.0998334166468282,1;-0.331800773085830,0.896351755264367,0.294043836551856,40;0.199473467763676,-0.237989126961396,0.950563785922063,60;0,0,0,1]
% %  
% tt1
% ];
N=max(data(:,2))+1;
M=size(data,1);
% N=4;
% M=20;

f=zeros(M*3,1);
tq_i=zeros(3,M);
tq_j=zeros(3,M);
TS=zeros(4*N,4);
NTS=zeros(4*N,4);
TTransform=zeros(1,N);
TTransform(1)=1;
TTransform(2)=0;
TTransform(3)=1;
TTransform(4)=0;
real_N=N-sum(TTransform);
map_T_realT=zeros(1,N);
J=zeros(M*3,6*real_N);
for i=1:N
    if TTransform(i)==0
%         TS(i*4-3:i*4,:)=[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1];
        TS(i*4-3:i*4,:)=Transform(i*4-3:i*4,:);
        map_T_realT(i)=i-sum(TTransform(1:i));
    else
        TS(i*4-3:i*4,:)=Transform(i*4-3:i*4,:);
    end
    NTS(i*4-3:i*4,:)=[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1];
end
% for i=1:N-1
%     TS(i*4+1:i*4+4,:)=tt1(4*i-3:4*i,:);
% end
for i=1:M
    T_i_index=data(i,1);
    T_j_index=data(i,2);
    cor_i=TS(4*T_i_index+1:4*T_i_index+4,:)*[data(i,3:5),1]';
    cor_j=TS(4*T_j_index+1:4*T_j_index+4,:)*[data(i,6:8),1]';
    tq_i(:,i)=[cor_i(1),cor_i(2),cor_i(3)]';
    tq_j(:,i)=[cor_j(1),cor_j(2),cor_j(3)]';
end
real_M=0;
for i=1:M
    T_i_index=data(i,1);
    T_j_index=data(i,2);
%   T_T=data(i,9);
    real_M=real_M+1;
    real_T_i=map_T_realT(T_i_index+1);
    real_T_j=map_T_realT(T_j_index+1);
    if real_T_i == 0 && real_T_j ==0 
        disp 'two fix';
        real_M=real_M-1;
    elseif real_T_i == 0
        J(3*real_M-2:3*real_M,6*real_T_j-5:6*real_T_j)= -[eye(3,3),-1*[0 -tq_j(3,i) tq_j(2,i);tq_j(3,i) 0 -tq_j(1,i);-tq_j(2,i) tq_j(1,i) 0]];
        f(3*real_M-2:3*real_M,1)=tq_i(:,i)-tq_j(:,i);
    elseif real_T_j == 0
        J(3*real_M-2:3*real_M,6*real_T_i-5:6*real_T_i)= [eye(3,3),-1*[0 -tq_i(3,i) tq_i(2,i);tq_i(3,i) 0 -tq_i(1,i);-tq_i(2,i) tq_i(1,i) 0]];
        f(3*real_M-2:3*real_M,1)=tq_i(:,i)-tq_j(:,i);
    else
        J(3*real_M-2:3*real_M,6*real_T_i-5:6*real_T_i)=[eye(3,3),-1*[0 -tq_i(3,i) tq_i(2,i);tq_i(3,i) 0 -tq_i(1,i);-tq_i(2,i) tq_i(1,i) 0]];
        J(3*real_M-2:3*real_M,6*real_T_j-5:6*real_T_j)=-[eye(3,3),-1*[0 -tq_j(3,i) tq_j(2,i);tq_j(3,i) 0 -tq_j(1,i);-tq_j(2,i) tq_j(1,i) 0]];
        f(3*real_M-2:3*real_M,1)=tq_i(:,i)-tq_j(:,i);
    end
    
%     if T_i == 0
%         T_T = 1;
%     else
%         T_T = 0;
%     end
%     if T_T==1
%        J(3*i-2:3*i,6*T_j-5:6*T_j)=[eye(3,3),-1*[0 -tq_j(3,i) tq_j(2,i);tq_j(3,i) 0 -tq_j(1,i);-tq_j(2,i) tq_j(1,i) 0]];
%        f(3*i-2:3*i,1)=tq_j(:,i)-tq_i(:,i);
%         J(3*i-2:3*i,6*T_j-5:6*T_j)= -[eye(3,3),-1*[0 -tq_j(3,i) tq_j(2,i);tq_j(3,i) 0 -tq_j(1,i);-tq_j(2,i) tq_j(1,i) 0]];
%     elseif T_T==2
%         J(3*i-2:3*i,6*T_i-5:6*T_i)=  [eye(3,3),-1*[0 -tq_i(3,i) tq_i(2,i);tq_j(3,i) 0 -tq_i(1,i);-tq_i(2,i) tq_i(1,i) 0]];
%        f(3*i-2:3*i,1)=tq_i(:,i)-tq_j(:,i);
%     else
%         J(3*i-2:3*i,6*T_i-5:6*T_i)=[eye(3,3),-1*[0 -tq_i(3,i) tq_i(2,i);tq_i(3,i) 0 -tq_i(1,i);-tq_i(2,i) tq_i(1,i) 0]];
%         J(3*i-2:3*i,6*T_j-5:6*T_j)=-[eye(3,3),-1*[0 -tq_j(3,i) tq_j(2,i);tq_j(3,i) 0 -tq_j(1,i);-tq_j(2,i) tq_j(1,i) 0]];
%         f(3*i-2:3*i,1)=tq_i(:,i)-tq_j(:,i);
%     end
%     f(3*i-2:3*i,1)=tq_i(:,i)-tq_j(:,i);
end
AJ=J(1:real_M*3,:);
Af=f(1:real_M*3,:);
   t=tic;
% detX=-(AJ'*AJ+0.001*diag(ones(1,(real_N)*6)))\(AJ'*Af);
    toc(t);
       t=tic;
 J = sparse(J);
%  L = ichol(J'*J);
 L = ichol(J'*J,struct('michol','on'));
[detX,fl2,rr2,it2,rv2]= pcg(J'*J,-J'*f,1e-10,200,L,L');
   toc(t);
% detX=dDetX;
% detX=-(J'*J+0.001*diag(ones(1,(N-1)*6)))\(J'*f);
for i=1:N
    real_T_i=map_T_realT(i);
    if real_T_i>0
        tran3f_i=detX(6*real_T_i-5:6*real_T_i-3);
        rot3f_i=detX(6*real_T_i-2:6*real_T_i);
        [rot3f_i_M,tran3f_i_M]=computerTran(rot3f_i,tran3f_i);
        NTS(i*4-3:i*4,:)=[rot3f_i_M,tran3f_i_M;0 0 0 1]*TS(i*4-3:i*4,:);
    else
        NTS(i*4-3:i*4,:)=TS(i*4-3:i*4,:);
    end
end
meam_init=f'*f/M;
max_init=max(f);
res=zeros(M*3,1);
for i=1:M/2
    T_i_index=data(i,1)+1;
    T_j_index=data(i,2)+1;
    pp=NTS(4*T_i_index-3:4*T_i_index,:)*[data(i,3:5),1]'-NTS(4*T_j_index-3:4*T_j_index,:)*[data(i,6:8),1]';
    res(3*i-2:3*i,1)=pp(1:3);
end
mean_opt=res'*res/M;
max_opt=max(res);
disp(['init error:',num2str(meam_init),',optimization error:',num2str(mean_opt),',init maxRes:',num2str(max_init),',opt maxRes:',num2str(max_opt)]);
% fprintf('初始化误差 %f，优化后误差 %f,初始最大反投影误差 %f，优化后反投影误差 %f\n',meam_init,mean_opt,max_init,max_opt)
% tran3f_i=detX(1:3);
% rot3f_i=detX(4:6);
% tran3f_j=detX(7:9);
% rot3f_j=detX(10:12);
% tran3f_k=detX(13:15);
% rot3f_k=detX(16:18);
%  rot3f_i_M=rotationVectorToMatrix(rot3f_i)';
%  rot3f_j_M=rotationVectorToMatrix(rot3f_j)';

% [rot3f_i_M,tran3f_i_M]=computerTran(rot3f_i,tran3f_i);
% [rot3f_j_M,tran3f_j_M]=computerTran(rot3f_j,tran3f_j);
% [rot3f_k_M,tran3f_k_M]=computerTran(rot3f_k,tran3f_k);
% T_i_N=[rot3f_i_M,tran3f_i_M;0 0 0 1]*T_i_0;
% T_j_N=[rot3f_j_M,tran3f_j_M;0 0 0 1]*T_j_0;
% T_k_N=[rot3f_k_M,tran3f_k_M;0 0 0 1]*T_k_0;
% norm(T_i*T_i_inv*p-T_j*T_j_inv*p)+norm(T_i*T_i_inv*p-T_k*T_k_inv*p)
% norm(T_i*p_i-T_j*p_j)+norm(T_i*p_i-T_k*p_k)
% norm(T_i_0*p_i-T_j_0*p_j)+norm(T_i_0*p_i-T_k_0*p_k)
% norm(T_i_N*p_i-T_j_N*p_j)+norm(T_i_N*p_i-T_k_N*p_k)