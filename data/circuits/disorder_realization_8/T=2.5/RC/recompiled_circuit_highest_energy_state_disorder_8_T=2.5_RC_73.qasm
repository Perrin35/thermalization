OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3266069) q[0];
sx q[0];
rz(-2.9171483) q[0];
sx q[0];
rz(1.8087968) q[0];
rz(2.3092071) q[1];
sx q[1];
rz(-1.7008984) q[1];
sx q[1];
rz(2.031215) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8827986) q[0];
sx q[0];
rz(-1.7340394) q[0];
sx q[0];
rz(1.4397653) q[0];
rz(-1.215685) q[2];
sx q[2];
rz(-0.74374226) q[2];
sx q[2];
rz(-1.6031934) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.346525) q[1];
sx q[1];
rz(-1.5988886) q[1];
sx q[1];
rz(-3.0884519) q[1];
rz(-pi) q[2];
rz(0.40544647) q[3];
sx q[3];
rz(-0.67908248) q[3];
sx q[3];
rz(-2.9630141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9139468) q[2];
sx q[2];
rz(-0.0089184428) q[2];
sx q[2];
rz(1.439636) q[2];
rz(1.7281744) q[3];
sx q[3];
rz(-0.011912502) q[3];
sx q[3];
rz(-0.2695151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6972167) q[0];
sx q[0];
rz(-1.5378636) q[0];
sx q[0];
rz(0.61475301) q[0];
rz(-2.6055824) q[1];
sx q[1];
rz(-3.1159846) q[1];
sx q[1];
rz(-2.8047309) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5081072) q[0];
sx q[0];
rz(-2.9683873) q[0];
sx q[0];
rz(2.2499535) q[0];
x q[1];
rz(2.6652671) q[2];
sx q[2];
rz(-1.5850955) q[2];
sx q[2];
rz(2.4313025) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.58648169) q[1];
sx q[1];
rz(-0.046849061) q[1];
sx q[1];
rz(1.2270801) q[1];
rz(-pi) q[2];
rz(0.19014374) q[3];
sx q[3];
rz(-1.2570253) q[3];
sx q[3];
rz(-2.6699175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.88379318) q[2];
sx q[2];
rz(-3.1287441) q[2];
sx q[2];
rz(-0.69965714) q[2];
rz(2.9720225) q[3];
sx q[3];
rz(-3.1278059) q[3];
sx q[3];
rz(2.5794896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70497847) q[0];
sx q[0];
rz(-2.5956557) q[0];
sx q[0];
rz(-0.28526947) q[0];
rz(2.8778695) q[1];
sx q[1];
rz(-0.00058760651) q[1];
sx q[1];
rz(-2.347351) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7274311) q[0];
sx q[0];
rz(-0.88849706) q[0];
sx q[0];
rz(-1.4201565) q[0];
rz(1.2474485) q[2];
sx q[2];
rz(-1.4686606) q[2];
sx q[2];
rz(2.7030526) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2546187) q[1];
sx q[1];
rz(-1.5363664) q[1];
sx q[1];
rz(1.5772343) q[1];
rz(0.66624347) q[3];
sx q[3];
rz(-2.1123721) q[3];
sx q[3];
rz(0.2592087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0245725) q[2];
sx q[2];
rz(-0.046155013) q[2];
sx q[2];
rz(-2.2771007) q[2];
rz(2.3965059) q[3];
sx q[3];
rz(-2.2741208) q[3];
sx q[3];
rz(0.043070506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2660148) q[0];
sx q[0];
rz(-0.046253007) q[0];
sx q[0];
rz(0.86638802) q[0];
rz(0.59151793) q[1];
sx q[1];
rz(-2.1951127) q[1];
sx q[1];
rz(-2.0902858) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6503904) q[0];
sx q[0];
rz(-1.5728381) q[0];
sx q[0];
rz(1.5709086) q[0];
x q[1];
rz(-1.573083) q[2];
sx q[2];
rz(-1.5703452) q[2];
sx q[2];
rz(1.5754896) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.72877872) q[1];
sx q[1];
rz(-0.83910131) q[1];
sx q[1];
rz(1.0207291) q[1];
rz(0.71612181) q[3];
sx q[3];
rz(-0.72343003) q[3];
sx q[3];
rz(1.0264068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6136578) q[2];
sx q[2];
rz(-0.064585678) q[2];
sx q[2];
rz(0.85395542) q[2];
rz(-0.69771403) q[3];
sx q[3];
rz(-1.8055975) q[3];
sx q[3];
rz(2.8310988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39603221) q[0];
sx q[0];
rz(-0.10265352) q[0];
sx q[0];
rz(2.7295617) q[0];
rz(1.283006) q[1];
sx q[1];
rz(-2.367815) q[1];
sx q[1];
rz(-0.7512908) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64941943) q[0];
sx q[0];
rz(-0.26628263) q[0];
sx q[0];
rz(2.6800781) q[0];
x q[1];
rz(0.89116781) q[2];
sx q[2];
rz(-1.749265) q[2];
sx q[2];
rz(-2.8986069) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4163782) q[1];
sx q[1];
rz(-2.0592923) q[1];
sx q[1];
rz(-2.4644771) q[1];
rz(-pi) q[2];
rz(-0.12215944) q[3];
sx q[3];
rz(-1.7137104) q[3];
sx q[3];
rz(1.1003332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6777307) q[2];
sx q[2];
rz(-0.025721392) q[2];
sx q[2];
rz(-0.99979293) q[2];
rz(0.60100466) q[3];
sx q[3];
rz(-0.060636245) q[3];
sx q[3];
rz(-2.291688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26337013) q[0];
sx q[0];
rz(-0.090494089) q[0];
sx q[0];
rz(-2.7409842) q[0];
rz(-1.406631) q[1];
sx q[1];
rz(-0.35146439) q[1];
sx q[1];
rz(-1.3946784) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18086223) q[0];
sx q[0];
rz(-0.020190857) q[0];
sx q[0];
rz(-1.3391311) q[0];
rz(-0.80446316) q[2];
sx q[2];
rz(-1.7906252) q[2];
sx q[2];
rz(1.1457844) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.38401022) q[1];
sx q[1];
rz(-1.3802802) q[1];
sx q[1];
rz(-1.5662929) q[1];
rz(-pi) q[2];
rz(-1.7441148) q[3];
sx q[3];
rz(-2.299274) q[3];
sx q[3];
rz(-2.4761379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3957735) q[2];
sx q[2];
rz(-2.5534111) q[2];
sx q[2];
rz(-1.0751209) q[2];
rz(-0.51946688) q[3];
sx q[3];
rz(-2.963701) q[3];
sx q[3];
rz(-1.5943257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9290685) q[0];
sx q[0];
rz(-1.2146177) q[0];
sx q[0];
rz(-1.529083) q[0];
rz(0.56135881) q[1];
sx q[1];
rz(-3.1293271) q[1];
sx q[1];
rz(-0.57048172) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9425526) q[0];
sx q[0];
rz(-1.6919475) q[0];
sx q[0];
rz(-0.051661927) q[0];
rz(-pi) q[1];
rz(-1.8884185) q[2];
sx q[2];
rz(-1.3537242) q[2];
sx q[2];
rz(1.0509911) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7509527) q[1];
sx q[1];
rz(-1.5696073) q[1];
sx q[1];
rz(-1.5863064) q[1];
x q[2];
rz(-0.53855207) q[3];
sx q[3];
rz(-0.90044124) q[3];
sx q[3];
rz(0.30645257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.059171112) q[2];
sx q[2];
rz(-2.8736281) q[2];
sx q[2];
rz(-1.1632261) q[2];
rz(1.6022812) q[3];
sx q[3];
rz(-3.0395165) q[3];
sx q[3];
rz(2.1479837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38692835) q[0];
sx q[0];
rz(-2.8480777) q[0];
sx q[0];
rz(2.5013404) q[0];
rz(-2.2786268) q[1];
sx q[1];
rz(-2.9997365) q[1];
sx q[1];
rz(-0.77756768) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3024738) q[0];
sx q[0];
rz(-0.57838744) q[0];
sx q[0];
rz(-2.3443787) q[0];
rz(-2.8993494) q[2];
sx q[2];
rz(-2.3628334) q[2];
sx q[2];
rz(2.3295516) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7217162) q[1];
sx q[1];
rz(-1.6437746) q[1];
sx q[1];
rz(3.130129) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6396995) q[3];
sx q[3];
rz(-1.1560901) q[3];
sx q[3];
rz(0.98687526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5295279) q[2];
sx q[2];
rz(-2.7175588) q[2];
sx q[2];
rz(0.5303793) q[2];
rz(-0.027675962) q[3];
sx q[3];
rz(-0.045462463) q[3];
sx q[3];
rz(1.3191222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6767122) q[0];
sx q[0];
rz(-0.16093971) q[0];
sx q[0];
rz(-0.93223923) q[0];
rz(1.5981916) q[1];
sx q[1];
rz(-2.3379969) q[1];
sx q[1];
rz(2.7995321) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0823183) q[0];
sx q[0];
rz(-3.0558028) q[0];
sx q[0];
rz(0.17753521) q[0];
rz(-pi) q[1];
rz(3.0717631) q[2];
sx q[2];
rz(-2.402555) q[2];
sx q[2];
rz(-2.7122517) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.910557) q[1];
sx q[1];
rz(-1.0536095) q[1];
sx q[1];
rz(1.3523577) q[1];
rz(-1.6430278) q[3];
sx q[3];
rz(-2.9119303) q[3];
sx q[3];
rz(-3.0185452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0590234) q[2];
sx q[2];
rz(-3.1408568) q[2];
sx q[2];
rz(1.9334582) q[2];
rz(2.1900603) q[3];
sx q[3];
rz(-0.0080684302) q[3];
sx q[3];
rz(-2.4212196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35219881) q[0];
sx q[0];
rz(-2.3645526) q[0];
sx q[0];
rz(0.081789516) q[0];
rz(0.31546047) q[1];
sx q[1];
rz(-3.0888562) q[1];
sx q[1];
rz(-1.2299406) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8393676) q[0];
sx q[0];
rz(-1.7291843) q[0];
sx q[0];
rz(1.7296289) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66218485) q[2];
sx q[2];
rz(-2.9176468) q[2];
sx q[2];
rz(2.0020773) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0198063) q[1];
sx q[1];
rz(-2.9933194) q[1];
sx q[1];
rz(-0.91839494) q[1];
rz(-3.0026912) q[3];
sx q[3];
rz(-1.0437696) q[3];
sx q[3];
rz(-2.4580818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6007467) q[2];
sx q[2];
rz(-0.019286152) q[2];
sx q[2];
rz(-1.1136327) q[2];
rz(0.039637808) q[3];
sx q[3];
rz(-3.1318635) q[3];
sx q[3];
rz(2.5447194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46350805) q[0];
sx q[0];
rz(-1.2152553) q[0];
sx q[0];
rz(-1.8334462) q[0];
rz(2.5254163) q[1];
sx q[1];
rz(-2.0694852) q[1];
sx q[1];
rz(-2.9328666) q[1];
rz(-1.3327333) q[2];
sx q[2];
rz(-1.3849266) q[2];
sx q[2];
rz(1.0746335) q[2];
rz(-0.039327903) q[3];
sx q[3];
rz(-1.7915184) q[3];
sx q[3];
rz(0.030727006) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
