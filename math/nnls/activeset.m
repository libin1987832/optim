function [x,y] = activeset(A,b)
%[x,y] = activeset(A,b)
%solves the linear least squares problem with nonnegative variables using the active set algorithm in [1].
%Input:
%   A:      [MxN] matrix 
%   b:      [Mx1] vector
%Output
%   x:      solution
%   y:      complementary solution
%
% [1] Portugal, Judice and Vicente, A comparison of block pivoting and
% interior point algorithms for linear least squares problems with
% nonnegative variables, Mathematics of Computation, 63(1994), pp. 625-643
%
%Uriel Roque
%28.4.2006
[m,n] = size(A);
F = [];
G = 1:n;
x = zeros(n,1);
Atb = A'*b;
y = -Atb;
noready = 1;
while noready
    
    %step 1    
    yG = y(G);
    if isempty(yG)
        break;
    end
    r = G(yG == min(yG));
    r = r(1);
    
    if (y(r) < 0) 
        H1 = [];
        H2 = r;
        F = union(setdiff(F,H1),H2);
        G = union(setdiff(G,H2),H1);
    else
        noready = 0;
        break;
    end
    
    noready2 = 1;
    while noready2;
    
        %step 2
        AF = A(:,F);
        xF = AF\b;
        nF = length(xF);
        if all(xF >= 0)
        
            x = [xF; zeros(n-nF,1)];
            %goto 3
            noready2 = 0;
            break;
       
        else
          
            index = find(xF<0);
            t = -x(index)./(xF(index) - x(index));
            tetha = min(t);
            r = F(index(t == tetha));
            r = r(1);
            x = [(1-tetha)*x(F) + tetha*xF ; zeros(n-nF,1)];
            H1 = r;
            H2 = [];
            F = union(setdiff(F,H1),H2);
            G = union(setdiff(G,H2),H1);
            %goto 2
            
        end
    
    end %while noready2
    
    %step 3
    AG = A(:,G);
    yG = AG' * (AF * xF - b);
    y(G) = yG;
        
end
%put results into their corresponding positions
a = x;
b = y;
x = zeros(n,1);
y = zeros(n,1);
x(F) = a(1:length(F));
y(G) = b(1:length(G));
