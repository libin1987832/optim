function [rotdetX,trandetX]=computerTran(rot3f,tran3f)
theta_sq=rot3f'*rot3f;
cr=cross(rot3f,tran3f);
theta=norm(rot3f);

if theta_sq < 1e-8
    disp 'theta_sq smallest'
    A=1-theta_sq/6;
    B=0.5;
    trandetX=tran3f+0.5*cr;
elseif theta_sq < 1e-6
    disp 'theta_sq small'
    C=(1-theta_sq/20)/6;
    A=1-theta_sq*C;
    B=0.5-0.25*theta_sq/6;
    w_cr=cross(rot3f,cr);
    trandetX=tran3f+B*cr+C*w_cr;
else
    itheta=1/theta;
    A=sin(theta)*itheta;
    B=(1-cos(theta))*(itheta*itheta);
    C=(1-A)*(itheta*itheta);
    w_cr=cross(rot3f,cr);
    trandetX=tran3f+B*cr+C*w_cr;
end
  rotdetX=rotationVectorToMatrix(rot3f)';
% wx=rot3f(1);
% wy=rot3f(2);
% wz=rot3f(3);
% wx2=wx^2;
% wy2=wy^2;
% wz2=wz^2;
% rotdetX=[1-B*(wy2+wz2),B*(wx*wy)-A*wz,B*(wx*wz)+A*wy;
%          B*(wx*wy)+A*wz,1-B*(wx2+wz2),B*(wy*wz)-A*wx;
%          B*(wx*wz)-A*wy,B*(wy*wz)+A*wx,1-B*(wx2+wy2)];