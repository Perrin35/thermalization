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
rz(-0.62682682) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1683567) q[0];
sx q[0];
rz(-1.3827948) q[0];
sx q[0];
rz(-3.0204535) q[0];
x q[1];
rz(-1.8148412) q[2];
sx q[2];
rz(-1.913117) q[2];
sx q[2];
rz(-2.2133881) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.37571733) q[1];
sx q[1];
rz(-0.98343508) q[1];
sx q[1];
rz(0.37382965) q[1];
rz(0.17625531) q[3];
sx q[3];
rz(-1.8415383) q[3];
sx q[3];
rz(-0.83812974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8864002) q[2];
sx q[2];
rz(-0.93078405) q[2];
sx q[2];
rz(1.2935151) q[2];
rz(1.1039929) q[3];
sx q[3];
rz(-2.2938426) q[3];
sx q[3];
rz(-0.75479341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8936546) q[0];
sx q[0];
rz(-2.0825443) q[0];
sx q[0];
rz(0.98518103) q[0];
rz(-0.39918104) q[1];
sx q[1];
rz(-1.2626516) q[1];
sx q[1];
rz(-3.0342297) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1185703) q[0];
sx q[0];
rz(-1.1623628) q[0];
sx q[0];
rz(-0.9704216) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58789247) q[2];
sx q[2];
rz(-1.8111472) q[2];
sx q[2];
rz(0.27871486) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0774035) q[1];
sx q[1];
rz(-1.940155) q[1];
sx q[1];
rz(-2.2662152) q[1];
x q[2];
rz(0.38856216) q[3];
sx q[3];
rz(-1.7266759) q[3];
sx q[3];
rz(-2.4678096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.62464109) q[2];
sx q[2];
rz(-1.3890283) q[2];
sx q[2];
rz(0.66169935) q[2];
rz(-3.0858357) q[3];
sx q[3];
rz(-2.7825833) q[3];
sx q[3];
rz(-2.8957446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4645638) q[0];
sx q[0];
rz(-0.55389535) q[0];
sx q[0];
rz(3.1012428) q[0];
rz(2.5616052) q[1];
sx q[1];
rz(-0.89825392) q[1];
sx q[1];
rz(-2.0592164) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4989935) q[0];
sx q[0];
rz(-1.9369619) q[0];
sx q[0];
rz(0.55143349) q[0];
rz(-pi) q[1];
rz(1.9627377) q[2];
sx q[2];
rz(-0.38092962) q[2];
sx q[2];
rz(-0.71289635) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3861474) q[1];
sx q[1];
rz(-1.8090994) q[1];
sx q[1];
rz(2.9134681) q[1];
x q[2];
rz(1.8941746) q[3];
sx q[3];
rz(-1.3591027) q[3];
sx q[3];
rz(-2.4993103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.13517705) q[2];
sx q[2];
rz(-1.8879031) q[2];
sx q[2];
rz(0.55389261) q[2];
rz(1.4931549) q[3];
sx q[3];
rz(-1.0529543) q[3];
sx q[3];
rz(2.5890787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4937113) q[0];
sx q[0];
rz(-2.557502) q[0];
sx q[0];
rz(-3.1062104) q[0];
rz(-2.1454051) q[1];
sx q[1];
rz(-1.6085947) q[1];
sx q[1];
rz(-0.48809537) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6513718) q[0];
sx q[0];
rz(-2.0351366) q[0];
sx q[0];
rz(0.086829348) q[0];
rz(-pi) q[1];
x q[1];
rz(0.79779063) q[2];
sx q[2];
rz(-0.75227037) q[2];
sx q[2];
rz(0.70300245) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.19255895) q[1];
sx q[1];
rz(-2.9395736) q[1];
sx q[1];
rz(2.2662524) q[1];
rz(-0.66588464) q[3];
sx q[3];
rz(-2.4443691) q[3];
sx q[3];
rz(-2.9149885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.35486832) q[2];
sx q[2];
rz(-1.0151261) q[2];
sx q[2];
rz(-0.50746894) q[2];
rz(-1.1821702) q[3];
sx q[3];
rz(-1.8488665) q[3];
sx q[3];
rz(1.8866084) q[3];
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
rz(0.59947157) q[0];
sx q[0];
rz(-2.4276955) q[0];
sx q[0];
rz(2.7713293) q[0];
rz(-0.26369357) q[1];
sx q[1];
rz(-1.5636684) q[1];
sx q[1];
rz(2.945074) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6949961) q[0];
sx q[0];
rz(-1.2906404) q[0];
sx q[0];
rz(-2.8269935) q[0];
rz(-2.0318803) q[2];
sx q[2];
rz(-1.5994674) q[2];
sx q[2];
rz(1.920514) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6424751) q[1];
sx q[1];
rz(-0.31458464) q[1];
sx q[1];
rz(0.59928526) q[1];
rz(-pi) q[2];
rz(1.7408095) q[3];
sx q[3];
rz(-1.2983054) q[3];
sx q[3];
rz(1.221399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1281517) q[2];
sx q[2];
rz(-2.1898966) q[2];
sx q[2];
rz(-0.91709843) q[2];
rz(-2.3069978) q[3];
sx q[3];
rz(-1.3093427) q[3];
sx q[3];
rz(2.807907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.82814) q[0];
sx q[0];
rz(-1.4414635) q[0];
sx q[0];
rz(2.939558) q[0];
rz(1.0728041) q[1];
sx q[1];
rz(-1.6348811) q[1];
sx q[1];
rz(1.7477759) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0840543) q[0];
sx q[0];
rz(-1.49465) q[0];
sx q[0];
rz(-2.8094859) q[0];
rz(-3.0300573) q[2];
sx q[2];
rz(-1.9431912) q[2];
sx q[2];
rz(-0.44484777) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.85642568) q[1];
sx q[1];
rz(-1.3060096) q[1];
sx q[1];
rz(-1.9734288) q[1];
x q[2];
rz(-1.1145669) q[3];
sx q[3];
rz(-0.8884512) q[3];
sx q[3];
rz(1.2896468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.94971925) q[2];
sx q[2];
rz(-1.6615901) q[2];
sx q[2];
rz(-2.2992772) q[2];
rz(-1.0874282) q[3];
sx q[3];
rz(-0.73173404) q[3];
sx q[3];
rz(1.4830164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52952805) q[0];
sx q[0];
rz(-1.8719215) q[0];
sx q[0];
rz(1.0213617) q[0];
rz(1.1313324) q[1];
sx q[1];
rz(-2.1919057) q[1];
sx q[1];
rz(2.9313415) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7267669) q[0];
sx q[0];
rz(-1.504577) q[0];
sx q[0];
rz(-1.7050752) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5049997) q[2];
sx q[2];
rz(-1.4074667) q[2];
sx q[2];
rz(2.3492299) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.47760689) q[1];
sx q[1];
rz(-2.7684282) q[1];
sx q[1];
rz(0.54877703) q[1];
rz(-pi) q[2];
rz(-2.290091) q[3];
sx q[3];
rz(-2.0264894) q[3];
sx q[3];
rz(-1.15629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.82688275) q[2];
sx q[2];
rz(-1.751739) q[2];
sx q[2];
rz(1.7604527) q[2];
rz(2.3794877) q[3];
sx q[3];
rz(-0.47587454) q[3];
sx q[3];
rz(-2.7281318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49522266) q[0];
sx q[0];
rz(-2.3587527) q[0];
sx q[0];
rz(-2.5543509) q[0];
rz(2.6953221) q[1];
sx q[1];
rz(-1.8621567) q[1];
sx q[1];
rz(-0.97672021) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5644972) q[0];
sx q[0];
rz(-1.7194887) q[0];
sx q[0];
rz(2.68481) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0301684) q[2];
sx q[2];
rz(-2.9615235) q[2];
sx q[2];
rz(-2.0153231) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8752012) q[1];
sx q[1];
rz(-2.5677263) q[1];
sx q[1];
rz(-0.89504524) q[1];
rz(-pi) q[2];
rz(0.85985698) q[3];
sx q[3];
rz(-1.1338187) q[3];
sx q[3];
rz(-0.45339282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4010767) q[2];
sx q[2];
rz(-2.433321) q[2];
sx q[2];
rz(0.55244279) q[2];
rz(2.6756514) q[3];
sx q[3];
rz(-1.7539141) q[3];
sx q[3];
rz(1.027511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4733646) q[0];
sx q[0];
rz(-0.79376525) q[0];
sx q[0];
rz(2.2209432) q[0];
rz(-0.051041516) q[1];
sx q[1];
rz(-1.5850681) q[1];
sx q[1];
rz(2.2424973) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2313401) q[0];
sx q[0];
rz(-1.1465342) q[0];
sx q[0];
rz(0.92696855) q[0];
x q[1];
rz(1.0663509) q[2];
sx q[2];
rz(-0.62014183) q[2];
sx q[2];
rz(1.8014185) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.60914674) q[1];
sx q[1];
rz(-1.6858613) q[1];
sx q[1];
rz(-1.902641) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9409688) q[3];
sx q[3];
rz(-1.5709953) q[3];
sx q[3];
rz(1.1699333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.94545323) q[2];
sx q[2];
rz(-3.0710199) q[2];
sx q[2];
rz(0.10458874) q[2];
rz(2.3032522) q[3];
sx q[3];
rz(-0.74931562) q[3];
sx q[3];
rz(-0.88275638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8024682) q[0];
sx q[0];
rz(-2.323928) q[0];
sx q[0];
rz(-1.1178281) q[0];
rz(-0.71240187) q[1];
sx q[1];
rz(-2.5103705) q[1];
sx q[1];
rz(0.39696473) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8660276) q[0];
sx q[0];
rz(-0.70789982) q[0];
sx q[0];
rz(-2.3562354) q[0];
rz(-pi) q[1];
rz(-0.88231477) q[2];
sx q[2];
rz(-2.6803662) q[2];
sx q[2];
rz(-3.0416833) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.18394477) q[1];
sx q[1];
rz(-2.7937104) q[1];
sx q[1];
rz(-1.8326879) q[1];
x q[2];
rz(-1.0730562) q[3];
sx q[3];
rz(-1.9314249) q[3];
sx q[3];
rz(-1.3090759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1931856) q[2];
sx q[2];
rz(-2.2414175) q[2];
sx q[2];
rz(0.096654264) q[2];
rz(1.7276673) q[3];
sx q[3];
rz(-2.0035412) q[3];
sx q[3];
rz(1.1269425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5861355) q[0];
sx q[0];
rz(-2.0270892) q[0];
sx q[0];
rz(-2.8591697) q[0];
rz(1.8718406) q[1];
sx q[1];
rz(-1.3602263) q[1];
sx q[1];
rz(0.7855986) q[1];
rz(-3.0336822) q[2];
sx q[2];
rz(-1.4866598) q[2];
sx q[2];
rz(1.286081) q[2];
rz(3.0222223) q[3];
sx q[3];
rz(-1.8277677) q[3];
sx q[3];
rz(1.2895332) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
