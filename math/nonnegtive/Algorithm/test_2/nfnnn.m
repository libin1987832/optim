
g1000=[0.0457999999999996;0.0474000000000003;0.0471999999999995;0.0462999999999989;0.0475000000000005;0.0475000000000001;0.0473000000000011;0.0477000000000002;0.0476000000000011;0.0484000000000009;0.0494999999999996;0.0482000000000001;0.0471999999999997;0.0464999999999998;0.0473999999999997;0.0498999999999991;0.0466000000000003;0.0483000000000005;0.0481999999999997;0.0479000000000003;0.0484000000000000];
g500=[6.22320000000000;6.20790000000000;6.19490000000000;6.21320000000000;6.01930000000000;6.06250000000000;6.07490000000000;5.92710000000000;5.99220000000000;5.84390000000000;6.03860000000000;5.91530000000000;5.95950000000000;5.94740000000000;5.80260000000000;5.90090000000000;6.03760000000000;5.91870000000000;5.89810000000000;5.94410000000000;5.88250000000000];
g100=[0.0937000000000001;0.0951999999999991;0.0980999999999995;0.0936000000000007;0.104599999999999;0.113600000000000;0.0923000000000002;0.0881999999999991;0.0791000000000004;0.0896000000000008;0.112400000000000;0.112500000000001;0.102800000000000;0.104900000000000;0.120499999999998;0.0989000000000000;0.103700000000000;0.103200000000001;0.124300000000001;0.126600000000001;0.119000000000000];

c1000=[0.0457000000000006;0.0479999999999999;0.0470999999999989;0.0487999999999993;0.0484999999999997;0.0470999999999994;0.0483999999999996;0.0494000000000002;0.0458000000000000;0.0487999999999993;0.0492000000000004;0.0481000000000003;0.0498999999999983;0.0484999999999998;0.0480999999999983;0.0468999999999995;0.0480999999999999;0.0491000000000005;0.0470999999999991;0.0493000000000006;0.0475999999999989];
c500=[5.71360000000000;5.67990000000000;5.71470000000000;5.69170000000000;5.67630000000000;5.67710000000000;5.70790000000000;5.67190000000000;5.68050000000000;5.55540000000000;5.65590000000000;5.60680000000000;5.60390000000000;5.70530000000000;5.63670000000000;5.69100000000000;5.70250000000000;5.68810000000000;5.77050000000000;5.65230000000000;5.70110000000000];
c100=[0.0916000000000008;0.0961000000000009;0.0924999999999990;0.0906999999999996;0.108600000000000;0.0961999999999996;0.0912999999999997;0.0901999999999987;0.0809000000000012;0.0905999999999985;0.110900000000000;0.109000000000000;0.0993999999999996;0.108400000000000;0.124499999999999;0.0980999999999996;0.106400000000000;0.105500000000001;0.125200000000001;0.131000000000001;0.116900000000000];

 p1000=[0.0549000000000004;0.0554999999999996;0.0537000000000005;0.0550000000000003;0.0553000000000002;0.0535999999999994;0.0545999999999997;0.0551999999999993;0.0541999999999994;0.0542000000000004;0.0547999999999999;0.0540999999999991;0.0553999999999990;0.0537000000000002;0.0548000000000012;0.0564999999999995;0.0546000000000001;0.0561999999999995;0.0536999999999998;0.0573000000000006;0.0526000000000009];
 p500=[6.43740000000000;6.21690000000000;6.30430000000000;6.29950000000000;6.15490000000000;6.16170000000000;6.06240000000000;6.12240000000000;6.07150000000000;6.20780000000000;6.09360000000000;6.10650000000000;6.03810000000000;5.97490000000000;6.03400000000000;6.02800000000000;5.96670000000000;6.02900000000000;6.00610000000000;5.94440000000000;5.86210000000000;0.0549000000000004];
 p100=[0.0964;0.0937000000000000;0.0906999999999996;0.0930000000000007;0.111300000000000;0.0953000000000003;0.0941999999999997;0.0891000000000002;0.0819000000000003;0.0914999999999992;0.111500000000000;0.115300000000000;0.0981999999999996;0.112700000000002;0.121699999999998;0.103300000000000;0.106400000000000;0.111700000000001;0.125600000000001;0.129500000000000;0.115000000000000;6.43740000000000];
 
 nfn=[10:30];
%  for i=1:21
%     fprintf('%d & %4.5f & %4.5f & %4.5f  \\\\ \n',nfn(i),c100(i),g100(i),p100(i));
%  end

 %for i=1:21
 %   fprintf('%d & %4.5f & %4.5f & %4.5f  \\\\ \n',nfn(i),c500(i),g500(i),p500(i));
 %end
  for i=1:21
    fprintf('%d & %4.5f & %4.5f & %4.5f  \\\\ \n',nfn(i),c1000(i),g1000(i),p1000(i));
 end