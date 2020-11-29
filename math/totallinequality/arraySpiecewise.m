function a = arraySpiecewise(A,b,x0,u)
nonzerou = u>0;
arange = -x0(nonzerou)./u(nonzerou);
aranges = sort(arange(arange > 0 & arange < 1));
aranges = [0,aranges];
aranges = unique(aranges);
all = [];
for i = 2:length(aranges)
    x = x0 + aranges( i - 1 )*u;
    Iu = x < 0;
    x(Iu) = 0;
    r = b - A * x;
    u(Iu) = 0;
    ap = A*u;
    ai=r./ap;
    alength = aranges( i ) - aranges( i -1);
    as=sort(ai(ai > 0 && alength >ai));
    ass = alength(i-1) + as;
    ass = [aranges( i-1 ) ass aranges(i)]; 
    all = [all,ass];
end


