function[error, R] = update(R,x)
%
% performs single row update of R
%
%   Beata Winnicka
%   Dept of Computer Science
%   Indiana University - Bloomington
%   1 October 1993 
%
%

    [m,p]=size(R);
%if(m<p)
%    error =1; 
%    return
%    R=[R; zeros(p-m,m),eye(p-m)];
%    disp(['had to resize R in update'])
%end

    for j=1:min(m,p)
        [G, Y]=planerot([R(j,j) x(j)]');
        t=G*[R(j,j:p);x(j:p)];
        x(j:p)=t(2,1:(p-j+1));
        R(j,j:p)=t(1,1:(p-j+1));
    end
    R=[R;x];
    error=0;

