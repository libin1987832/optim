N=2000;
p=rand([3,N])*100;
p=[p;ones(1,N)];

yaw_i = 0.5854; pitch_i = 0.3; roll_i = 0.2; tran_i=[1,3,4];
dcm_i = angle2dcm( yaw_i, pitch_i, roll_i );
T_i =  [dcm_i,tran_i';0 0 0 1];
T_i_0 = [angle2dcm( yaw_i+0.01, pitch_i+0.02, roll_i+0.01 ),tran_i'+[0.01;0.02;0.01];0 0 0 1];
T_i_inv=[dcm_i',-dcm_i'*tran_i';0 0 0 1];
p_i = T_i_inv*p+[rand([3,N])*0.01;zeros(1,N)];

yaw_j = 0.4854; pitch_j = 0.5; roll_j = 0.41; tran_j=[1,4,6];
dcm_j = angle2dcm( yaw_j, pitch_j, roll_j );
T_j = [dcm_j,tran_j';0 0 0 1];
T_j_0 = [angle2dcm( yaw_j+0.01, pitch_j+0.02, roll_j+0.01 ),tran_j'+[0.01;0.02;0.01];0 0 0 1];
T_j_inv=[dcm_j',-dcm_j'*tran_j';0 0 0 1];
p_j = T_j_inv*p+[rand([3,N])*0.01;zeros(1,N)];

J=[];
f=[];
tq_i=T_i_0*p_i;
tq_j=T_j_0*p_j;
for i=1:N
    Jtq_i=[eye(3,3),-1*[0 -tq_i(3,i) tq_i(2,i);tq_i(3,i) 0 -tq_i(1,i);-tq_i(2,i) tq_i(1,i) 0]];
    Jtq_j=[-eye(3,3),  [0 -tq_j(3,i) tq_j(2,i);tq_j(3,i) 0 -tq_j(1,i);-tq_j(2,i) tq_j(1,i) 0]];
    J=[J ; Jtq_i Jtq_j];
    pp=T_i_0*p_i(:,i)-T_j_0*p_j(:,i);
    f=[f;pp(1);pp(2);pp(3)];
end

% r0=-J'*f-J'*J*[T_i_0;T_j_0]

detX=-(J'*J)\(J'*f);
tran3f_i=detX(1:3);
rot3f_i=detX(4:6);
tran3f_j=detX(7:9);
rot3f_j=detX(10:12);

%  rot3f_i_M=rotationVectorToMatrix(rot3f_i)';
%  rot3f_j_M=rotationVectorToMatrix(rot3f_j)';

[rot3f_i_M,tran3f_i_M]=computerTran(rot3f_i,tran3f_i);
[rot3f_j_M,tran3f_j_M]=computerTran(rot3f_j,tran3f_j);

T_i_N=[rot3f_i_M,tran3f_i_M;0 0 0 1]*T_i_0;
T_j_N=[rot3f_j_M,tran3f_j_M;0 0 0 1]*T_j_0;

norm(T_i*T_i_inv*p-T_j*T_j_inv*p)
norm(T_i*p_i-T_j*p_j)
norm(T_i_0*p_i-T_j_0*p_j)
norm(T_i_N*p_i-T_j_N*p_j)