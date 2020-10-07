filename = 'breast-cancer-wisconsin.data';
format = '%n%n%n%n%n%n%n%n%n%n%n';

[data1,data2,data3,data4,data5,data6,data7,data8,data9,data10,data11] = textread(filename , format , 'delimiter', ',');
class2 = (data11 == 2);
class4 = (data11 == 4);
AB = [data2,data3,data4,data5,data6,data7,data8,data9,data10];
A = AB(class2);
B = AB(class4);