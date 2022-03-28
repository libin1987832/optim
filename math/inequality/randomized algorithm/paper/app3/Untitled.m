str=['P','R','P','C','N'];
for i = 1:5
r = -u + A1 * xs(:,i);
rnum1=sum(r>0);
r = l - A2 * xs(:,i);
rnum2=sum(r>0);
r=b-H*xs(:,i);
r(r<0)=0;
r_GS = norm(r);
rum3=sum(x<0);
fprintf('& %sCD &%d &%d &%d & %g & %g \\\\\n', str(i), rnum1,rnum2,rum3,r_GS, xst(i));
end