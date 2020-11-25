function [A1, b1, A2, b2, At1, At2, fm1, fm2, AL1, bL1, AL2, bL2] = readHeart(gamm)
filename = 'processed.cleveland.data';
format = '%n%n%n%n%n%n%n%n%n%n%n%n%n%n';
[data1,data2,data3,data4,data5,data6,data7,data8,data9,data10,data11,data12,data13,data14] = textread(filename , format , 'delimiter', ',');
AB = [data1,data2,data3,data4,data5,data6,data7,data8,data9,data10,data11,data12,data13,data14];
% if there is the nonnegative value, -1 so sum great 0 
nonnegative = not(sum(AB<0,2));
AB=AB(nonnegative, : );
class2 = (AB( : , end ) == 0);
class4 = (AB( : , end ) > 0);
A = AB(class2,1:end-1);% postive
B = AB(class4,1:end-1);% non positive
[m1,n1] = size(A);
[m2,n2] = size(B);
rm1 = randperm(m1);
rm2 = randperm(m2);
fm1 = floor(m1/3*2);
fm2 = floor(m2/3*2);
A1 = [-A(rm1(1:fm1),:);B(rm2(1:fm2),:)];% train data
b1 = [(1 - gamm)*ones(fm1,1);(1+gamm)*ones(fm2,1)];% 
A2 = [-A(rm1(1:fm1),:), ones(fm1,1);B(rm2(1:fm2),:), -ones(fm2,1)];% gamm also is parametra
b2 = ones(fm1+fm2, 1);
% test data
At1 = A(rm1(fm1+1:end),:);
At2 = B(rm2(fm2+1:end),:);

AL1 = -[B(rm2(1:fm2),:);-A(rm1(1:fm1),:)];
AL1 = [AL1, -eye(size(AL1,1))];
bL1 = [-(1+gamm)*ones(fm2,1);-(1 - gamm)*ones(fm1,1)];

AL2 = -[B(rm2(1:fm2),:) -ones(fm2,1);-A(rm1(1:fm1),:) ones(fm1,1)];
AL2 = [AL2, -eye(size(AL2,1))];
bL2 = [-ones(fm2,1);-ones(fm1,1)];


