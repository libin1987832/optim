N=200000;
p=rand([3,N])*100;
p=[p;ones(1,N)];
% p=[76.1970638864839,5.98170332245082,3.19635005438385,1.20592196038944,26.6226961629021,77.2538811819699,54.2985240569339,17.3408933462850,66.7090315156029,68.9027527983544;11.2879308313449,72.5013455549013,61.7829662063213,84.7864536455246,79.5587688244189,47.2764802896122,6.58595255023609,13.5600359400232,23.6448663451618,58.3249631802497;69.2189789903680,25.1168265388705,51.8971165452512,37.2179849793529,49.1392749898114,4.78548473023002,35.7950408601128,47.1028852661705,19.4276944620092,77.7128433485589;1,1,1,1,1,1,1,1,1,1];
% load p_3.mat;
yaw_i = 0.5854; pitch_i = 0.3; roll_i = 0.2; tran_i=[1,3,4];
dcm_i = angle2dcm( yaw_i, pitch_i, roll_i );
T_i =  [dcm_i,tran_i';0 0 0 1];
T_i_0 = [angle2dcm( yaw_i+0.01, pitch_i+0.02, roll_i+0.01 ),tran_i'+[0.01;0.02;0.01];0 0 0 1];
T_i_inv=[dcm_i',-dcm_i'*tran_i';0 0 0 1];
% p_i = T_i_inv*p+[rand([3,N])*0.01;zeros(1,N)];
p_i = T_i_inv*p+[ones(3,N)*0.01;zeros(1,N)];

yaw_j = 0.4854; pitch_j = 0.5; roll_j = 0.41; tran_j=[1,4,6];
dcm_j = angle2dcm( yaw_j, pitch_j, roll_j );
T_j = [dcm_j,tran_j';0 0 0 1];
T_j_0 = [angle2dcm( yaw_j+0.01, pitch_j+0.02, roll_j+0.01 ),tran_j'+[0.01;0.02;0.01];0 0 0 1];
T_j_inv=[dcm_j',-dcm_j'*tran_j';0 0 0 1];
% p_j = T_j_inv*p+[rand([3,N])*0.01;zeros(1,N)];
p_j = T_j_inv*p+[ones(3,N)*0.01;zeros(1,N)];

yaw_k = 0.3854; pitch_k = 0.1; roll_k = 0.3; tran_k=[1,40,60];
dcm_k = angle2dcm( yaw_k, pitch_k, roll_k );
T_k = [dcm_k,tran_k';0 0 0 1];
% T_k_0 = [angle2dcm( yaw_k+0.01, pitch_k+0.02, roll_k+0.01 ),tran_k'+[0.01;0.02;0.01];0 0 0 1];
T_k_0=T_k ;
T_k_inv=[dcm_k',-dcm_k'*tran_k';0 0 0 1];
% p_k = T_k_inv*p+[rand([3,N])*0.01;zeros(1,N)];
 p_k = T_k_inv*p+[ones(3,N)*0.01;zeros(1,N)];
ss1=[ones(N,1),ones(N,1)+1,p_i(1:3,:)',p_j(1:3,:)'];
ss2=[ones(N,1),ones(N,1)+2,p_i(1:3,:)',p_k(1:3,:)'];
tt1=[T_i_0;T_j_0;T_k_0];
J=[];
f=[];
tq_i=T_i_0*p_i;
tq_j=T_j_0*p_j;
tq_k=T_k_0*p_k;
f=zeros(2*N*3,1);
J=zeros(2*N*3,12);
for i=1:N
    Jtq_i=[eye(3,3),-1*[0 -tq_i(3,i) tq_i(2,i);tq_i(3,i) 0 -tq_i(1,i);-tq_i(2,i) tq_i(1,i) 0]];
    Jtq_j=[eye(3,3),-1*[0 -tq_j(3,i) tq_j(2,i);tq_j(3,i) 0 -tq_j(1,i);-tq_j(2,i) tq_j(1,i) 0]];
     Jtq_k=[eye(3,3),-1*[0 -tq_k(3,i) tq_k(2,i);tq_k(3,i) 0 -tq_k(1,i);-tq_k(2,i) tq_k(1,i) 0]];
    J(3*i-2:3*i,:)=[Jtq_i -Jtq_j];
     J(3*i-2+3*N:3*i+3*N,:)=[Jtq_i zeros(3,6) -Jtq_k];
    J(3*i-2+3*N:3*i+3*N,:)=[Jtq_i zeros(3,6)];
    pp=T_i_0*p_i(:,i)-T_j_0*p_j(:,i);
    pp2=T_i_0*p_i(:,i)-T_k_0*p_k(:,i);
    f(3*i-2:3*i,:)=[pp(1);pp(2);pp(3)];
    f(3*i-2+3*N:3*i+N*3,:)=[pp2(1);pp2(2);pp2(3)];
end

% r0=-J'*f-J'*J*[T_i_0;T_j_0]
JTT=J;
fffTT=f;
J';
J'*J;
detX=-(J'*J+0.001*diag(ones(1,12)))\(J'*f);
fxx=f;
detxx=detX;
tran3f_i=detX(1:3);
rot3f_i=detX(4:6);
tran3f_j=detX(7:9);
rot3f_j=detX(10:12);
% tran3f_k=detX(13:15);
% rot3f_k=detX(16:18);
%  rot3f_i_M=rotationVectorToMatrix(rot3f_i)';
%  rot3f_j_M=rotationVectorToMatrix(rot3f_j)';

[rot3f_i_M,tran3f_i_M]=computerTran(rot3f_i,tran3f_i);
[rot3f_j_M,tran3f_j_M]=computerTran(rot3f_j,tran3f_j);
% [rot3f_k_M,tran3f_k_M]=computerTran(rot3f_k,tran3f_k);

T_i_N=[rot3f_i_M,tran3f_i_M;0 0 0 1]*T_i_0;
T_j_N=[rot3f_j_M,tran3f_j_M;0 0 0 1]*T_j_0;
% T_k_N=[rot3f_k_M,tran3f_k_M;0 0 0 1]*T_k_0;
T_k_N=T_k_0;

norm(T_i*T_i_inv*p-T_j*T_j_inv*p)+norm(T_i*T_i_inv*p-T_k*T_k_inv*p);
norm(T_i*p_i-T_j*p_j)+norm(T_i*p_i-T_k*p_k);
norm(T_i_0*p_i-T_j_0*p_j)^2+norm(T_i_0*p_i-T_k_0*p_k)^2;
norm(T_i_N*p_i-T_j_N*p_j)+norm(T_i_N*p_i-T_k_N*p_k);

f_N=zeros(2*N*3,1);
for i=1:N
    pp=T_i_N*p_i(:,i)-T_j_N*p_j(:,i);
    pp2=T_i_N*p_i(:,i)-T_k_N*p_k(:,i);
    f_N(3*i-2:3*i,:)=[pp(1);pp(2);pp(3)];
    f_N(3*i-2+N*3:3*i+N*3,:)=[pp2(1);pp2(2);pp2(3)];
end
f'*f
f_N'*f_N
max(max(f_N))