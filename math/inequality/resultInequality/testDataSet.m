function output = testDataSet(A1,At1,At2,fm1,fm2,xkh,gamm)
for k=1:2
    if k==2
    tAt1 = At1;
    tAt2 = At2;
    At1 = A1(1:fm1,:);
    At2 = A1(fm1+1:fm1+fm2,:);
    end
    for i = 1:2
    count1 = sum(At1*xkh(:,i) < gamm(1,i));
    count2 = size(At1,1);
%     percent1 = count1/count2;
%     [-count1+count2 percent1]
    count21 = sum(At2*xkh(:,i) > gamm(1,i));
    count22 = size(At2,1);
%     percent2 = count21/count22;
%     [-count21+count22 percent2]
    sumerror = -count21+count22-count1+count2;
    sumcount = count2+count22;
    output(1,3*(i+2*k-2)-2:3*(i+2 *k-2)) = [sumerror, sumcount, 1-sumerror/sumcount];
end
end