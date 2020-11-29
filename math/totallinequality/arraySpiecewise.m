function a = arraySpiecewise(A,b,x,u)
nonzerou = u>0;
arange = -x(nonzerou)./u(nonzerou);
aranges = sort(arange(arange > 0 & arange < 1));

