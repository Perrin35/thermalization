OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.964736) q[0];
sx q[0];
rz(-2.2778947) q[0];
sx q[0];
rz(2.9404844) q[0];
rz(-1.3970628) q[1];
sx q[1];
rz(-1.8048598) q[1];
sx q[1];
rz(2.5147658) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7456822) q[0];
sx q[0];
rz(-2.9183309) q[0];
sx q[0];
rz(-1.0049055) q[0];
rz(-pi) q[1];
rz(-2.7896499) q[2];
sx q[2];
rz(-1.8004187) q[2];
sx q[2];
rz(-0.72598347) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.37571733) q[1];
sx q[1];
rz(-2.1581576) q[1];
sx q[1];
rz(-2.767763) q[1];
x q[2];
rz(-1.2960035) q[3];
sx q[3];
rz(-1.7405675) q[3];
sx q[3];
rz(-0.78026375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8864002) q[2];
sx q[2];
rz(-2.2108086) q[2];
sx q[2];
rz(1.8480776) q[2];
rz(2.0375997) q[3];
sx q[3];
rz(-2.2938426) q[3];
sx q[3];
rz(-2.3867992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8936546) q[0];
sx q[0];
rz(-2.0825443) q[0];
sx q[0];
rz(-0.98518103) q[0];
rz(0.39918104) q[1];
sx q[1];
rz(-1.2626516) q[1];
sx q[1];
rz(-0.10736297) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3283008) q[0];
sx q[0];
rz(-2.11587) q[0];
sx q[0];
rz(-2.6585447) q[0];
x q[1];
rz(1.2843578) q[2];
sx q[2];
rz(-1.0019433) q[2];
sx q[2];
rz(-1.4494277) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0641891) q[1];
sx q[1];
rz(-1.940155) q[1];
sx q[1];
rz(2.2662152) q[1];
x q[2];
rz(1.739005) q[3];
sx q[3];
rz(-1.1871927) q[3];
sx q[3];
rz(-0.83354359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.62464109) q[2];
sx q[2];
rz(-1.7525643) q[2];
sx q[2];
rz(-0.66169935) q[2];
rz(-3.0858357) q[3];
sx q[3];
rz(-2.7825833) q[3];
sx q[3];
rz(-2.8957446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4645638) q[0];
sx q[0];
rz(-0.55389535) q[0];
sx q[0];
rz(0.04034986) q[0];
rz(-0.57998747) q[1];
sx q[1];
rz(-2.2433387) q[1];
sx q[1];
rz(2.0592164) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45535116) q[0];
sx q[0];
rz(-0.65128122) q[0];
sx q[0];
rz(2.5097646) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9627377) q[2];
sx q[2];
rz(-2.760663) q[2];
sx q[2];
rz(2.4286963) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1632417) q[1];
sx q[1];
rz(-0.32838531) q[1];
sx q[1];
rz(0.82113012) q[1];
rz(-pi) q[2];
rz(-1.247418) q[3];
sx q[3];
rz(-1.3591027) q[3];
sx q[3];
rz(-2.4993103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0064156) q[2];
sx q[2];
rz(-1.2536896) q[2];
sx q[2];
rz(-2.5877) q[2];
rz(1.4931549) q[3];
sx q[3];
rz(-2.0886383) q[3];
sx q[3];
rz(0.55251399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64788139) q[0];
sx q[0];
rz(-0.58409062) q[0];
sx q[0];
rz(-0.035382263) q[0];
rz(2.1454051) q[1];
sx q[1];
rz(-1.532998) q[1];
sx q[1];
rz(2.6534973) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29823829) q[0];
sx q[0];
rz(-2.6697864) q[0];
sx q[0];
rz(1.7422373) q[0];
rz(-pi) q[1];
x q[1];
rz(0.79779063) q[2];
sx q[2];
rz(-0.75227037) q[2];
sx q[2];
rz(-2.4385902) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9490337) q[1];
sx q[1];
rz(-2.9395736) q[1];
sx q[1];
rz(-0.87534027) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.048269) q[3];
sx q[3];
rz(-2.100088) q[3];
sx q[3];
rz(-2.570591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7867243) q[2];
sx q[2];
rz(-1.0151261) q[2];
sx q[2];
rz(2.6341237) q[2];
rz(-1.9594225) q[3];
sx q[3];
rz(-1.2927262) q[3];
sx q[3];
rz(1.8866084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59947157) q[0];
sx q[0];
rz(-0.71389714) q[0];
sx q[0];
rz(-2.7713293) q[0];
rz(0.26369357) q[1];
sx q[1];
rz(-1.5779243) q[1];
sx q[1];
rz(2.945074) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4465966) q[0];
sx q[0];
rz(-1.2906404) q[0];
sx q[0];
rz(2.8269935) q[0];
rz(-2.0318803) q[2];
sx q[2];
rz(-1.5421252) q[2];
sx q[2];
rz(1.2210786) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.12381183) q[1];
sx q[1];
rz(-1.31243) q[1];
sx q[1];
rz(-1.3892795) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7408095) q[3];
sx q[3];
rz(-1.8432872) q[3];
sx q[3];
rz(1.221399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1281517) q[2];
sx q[2];
rz(-2.1898966) q[2];
sx q[2];
rz(2.2244942) q[2];
rz(-2.3069978) q[3];
sx q[3];
rz(-1.83225) q[3];
sx q[3];
rz(-2.807907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31345263) q[0];
sx q[0];
rz(-1.7001292) q[0];
sx q[0];
rz(0.20203461) q[0];
rz(1.0728041) q[1];
sx q[1];
rz(-1.6348811) q[1];
sx q[1];
rz(-1.3938168) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4112339) q[0];
sx q[0];
rz(-0.3404091) q[0];
sx q[0];
rz(-0.22986869) q[0];
rz(-pi) q[1];
rz(1.2932111) q[2];
sx q[2];
rz(-0.38799122) q[2];
sx q[2];
rz(0.1462305) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2653633) q[1];
sx q[1];
rz(-0.47788844) q[1];
sx q[1];
rz(2.1761314) q[1];
rz(2.4059541) q[3];
sx q[3];
rz(-1.9197575) q[3];
sx q[3];
rz(0.58128202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1918734) q[2];
sx q[2];
rz(-1.6615901) q[2];
sx q[2];
rz(-0.8423155) q[2];
rz(1.0874282) q[3];
sx q[3];
rz(-2.4098586) q[3];
sx q[3];
rz(1.4830164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6120646) q[0];
sx q[0];
rz(-1.2696711) q[0];
sx q[0];
rz(2.1202309) q[0];
rz(2.0102603) q[1];
sx q[1];
rz(-2.1919057) q[1];
sx q[1];
rz(-2.9313415) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70049858) q[0];
sx q[0];
rz(-2.9919618) q[0];
sx q[0];
rz(2.0307226) q[0];
rz(-pi) q[1];
rz(-2.7619751) q[2];
sx q[2];
rz(-2.9656177) q[2];
sx q[2];
rz(-1.1773674) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5308295) q[1];
sx q[1];
rz(-1.3794583) q[1];
sx q[1];
rz(2.8192239) q[1];
x q[2];
rz(-0.85150163) q[3];
sx q[3];
rz(-1.1151033) q[3];
sx q[3];
rz(-1.15629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.82688275) q[2];
sx q[2];
rz(-1.3898536) q[2];
sx q[2];
rz(1.7604527) q[2];
rz(2.3794877) q[3];
sx q[3];
rz(-2.6657181) q[3];
sx q[3];
rz(2.7281318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.64637) q[0];
sx q[0];
rz(-0.78283993) q[0];
sx q[0];
rz(-0.58724171) q[0];
rz(-2.6953221) q[1];
sx q[1];
rz(-1.279436) q[1];
sx q[1];
rz(-0.97672021) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5644972) q[0];
sx q[0];
rz(-1.7194887) q[0];
sx q[0];
rz(0.45678267) q[0];
x q[1];
rz(3.061053) q[2];
sx q[2];
rz(-1.7320247) q[2];
sx q[2];
rz(1.5921519) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.8969438) q[1];
sx q[1];
rz(-1.9172501) q[1];
sx q[1];
rz(-1.1036554) q[1];
x q[2];
rz(0.55240734) q[3];
sx q[3];
rz(-2.2033785) q[3];
sx q[3];
rz(1.6747024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7405159) q[2];
sx q[2];
rz(-0.70827168) q[2];
sx q[2];
rz(-0.55244279) q[2];
rz(-2.6756514) q[3];
sx q[3];
rz(-1.7539141) q[3];
sx q[3];
rz(2.1140816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4733646) q[0];
sx q[0];
rz(-2.3478274) q[0];
sx q[0];
rz(-2.2209432) q[0];
rz(-3.0905511) q[1];
sx q[1];
rz(-1.5565245) q[1];
sx q[1];
rz(2.2424973) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2313401) q[0];
sx q[0];
rz(-1.9950584) q[0];
sx q[0];
rz(0.92696855) q[0];
rz(-0.33234889) q[2];
sx q[2];
rz(-2.1045448) q[2];
sx q[2];
rz(-0.74408434) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5324459) q[1];
sx q[1];
rz(-1.4557314) q[1];
sx q[1];
rz(-1.902641) q[1];
x q[2];
rz(-1.2006239) q[3];
sx q[3];
rz(-1.5705974) q[3];
sx q[3];
rz(-1.9716594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.94545323) q[2];
sx q[2];
rz(-0.070572704) q[2];
sx q[2];
rz(-0.10458874) q[2];
rz(2.3032522) q[3];
sx q[3];
rz(-2.392277) q[3];
sx q[3];
rz(-2.2588363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33912441) q[0];
sx q[0];
rz(-0.81766468) q[0];
sx q[0];
rz(1.1178281) q[0];
rz(-0.71240187) q[1];
sx q[1];
rz(-2.5103705) q[1];
sx q[1];
rz(0.39696473) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8660276) q[0];
sx q[0];
rz(-2.4336928) q[0];
sx q[0];
rz(0.78535725) q[0];
x q[1];
rz(-2.2592779) q[2];
sx q[2];
rz(-0.46122641) q[2];
sx q[2];
rz(-3.0416833) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9576479) q[1];
sx q[1];
rz(-0.34788222) q[1];
sx q[1];
rz(-1.8326879) q[1];
rz(-pi) q[2];
rz(2.2393353) q[3];
sx q[3];
rz(-2.5359631) q[3];
sx q[3];
rz(2.827364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.94840702) q[2];
sx q[2];
rz(-0.90017515) q[2];
sx q[2];
rz(-3.0449384) q[2];
rz(-1.4139253) q[3];
sx q[3];
rz(-2.0035412) q[3];
sx q[3];
rz(1.1269425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5554572) q[0];
sx q[0];
rz(-2.0270892) q[0];
sx q[0];
rz(-2.8591697) q[0];
rz(-1.8718406) q[1];
sx q[1];
rz(-1.7813663) q[1];
sx q[1];
rz(-2.355994) q[1];
rz(-0.10791049) q[2];
sx q[2];
rz(-1.6549329) q[2];
sx q[2];
rz(-1.8555117) q[2];
rz(0.11937033) q[3];
sx q[3];
rz(-1.3138249) q[3];
sx q[3];
rz(-1.8520595) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];