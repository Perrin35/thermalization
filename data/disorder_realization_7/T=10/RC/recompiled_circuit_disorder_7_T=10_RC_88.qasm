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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62030828) q[0];
sx q[0];
rz(-1.6897908) q[0];
sx q[0];
rz(1.760153) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5457382) q[2];
sx q[2];
rz(-2.7239954) q[2];
sx q[2];
rz(2.8516304) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.37571733) q[1];
sx q[1];
rz(-0.98343508) q[1];
sx q[1];
rz(2.767763) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1342282) q[3];
sx q[3];
rz(-0.3218739) q[3];
sx q[3];
rz(-0.25062996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.25519249) q[2];
sx q[2];
rz(-2.2108086) q[2];
sx q[2];
rz(-1.2935151) q[2];
rz(-2.0375997) q[3];
sx q[3];
rz(-0.84775001) q[3];
sx q[3];
rz(-2.3867992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8936546) q[0];
sx q[0];
rz(-2.0825443) q[0];
sx q[0];
rz(-0.98518103) q[0];
rz(-2.7424116) q[1];
sx q[1];
rz(-1.2626516) q[1];
sx q[1];
rz(3.0342297) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8132919) q[0];
sx q[0];
rz(-1.0257226) q[0];
sx q[0];
rz(-0.48304793) q[0];
x q[1];
rz(2.5537002) q[2];
sx q[2];
rz(-1.3304454) q[2];
sx q[2];
rz(-0.27871486) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.7992236) q[1];
sx q[1];
rz(-2.2111726) q[1];
sx q[1];
rz(0.46701042) q[1];
rz(-pi) q[2];
rz(-1.739005) q[3];
sx q[3];
rz(-1.1871927) q[3];
sx q[3];
rz(-2.3080491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.62464109) q[2];
sx q[2];
rz(-1.7525643) q[2];
sx q[2];
rz(0.66169935) q[2];
rz(-0.055756904) q[3];
sx q[3];
rz(-0.35900933) q[3];
sx q[3];
rz(0.24584809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67702883) q[0];
sx q[0];
rz(-2.5876973) q[0];
sx q[0];
rz(-0.04034986) q[0];
rz(-2.5616052) q[1];
sx q[1];
rz(-2.2433387) q[1];
sx q[1];
rz(1.0823762) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6425991) q[0];
sx q[0];
rz(-1.2046308) q[0];
sx q[0];
rz(2.5901592) q[0];
rz(-0.15180397) q[2];
sx q[2];
rz(-1.9215343) q[2];
sx q[2];
rz(1.1317859) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1632417) q[1];
sx q[1];
rz(-0.32838531) q[1];
sx q[1];
rz(2.3204625) q[1];
rz(-pi) q[2];
rz(1.8941746) q[3];
sx q[3];
rz(-1.3591027) q[3];
sx q[3];
rz(0.64228234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.13517705) q[2];
sx q[2];
rz(-1.2536896) q[2];
sx q[2];
rz(2.5877) q[2];
rz(1.4931549) q[3];
sx q[3];
rz(-1.0529543) q[3];
sx q[3];
rz(2.5890787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64788139) q[0];
sx q[0];
rz(-2.557502) q[0];
sx q[0];
rz(-3.1062104) q[0];
rz(0.9961876) q[1];
sx q[1];
rz(-1.532998) q[1];
sx q[1];
rz(0.48809537) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29823829) q[0];
sx q[0];
rz(-0.47180628) q[0];
sx q[0];
rz(1.3993553) q[0];
rz(0.98056356) q[2];
sx q[2];
rz(-2.0681941) q[2];
sx q[2];
rz(-2.892708) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0779873) q[1];
sx q[1];
rz(-1.4418774) q[1];
sx q[1];
rz(-1.414826) q[1];
x q[2];
rz(-2.5591764) q[3];
sx q[3];
rz(-1.1629259) q[3];
sx q[3];
rz(1.2553314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.35486832) q[2];
sx q[2];
rz(-2.1264666) q[2];
sx q[2];
rz(0.50746894) q[2];
rz(1.9594225) q[3];
sx q[3];
rz(-1.2927262) q[3];
sx q[3];
rz(-1.8866084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59947157) q[0];
sx q[0];
rz(-0.71389714) q[0];
sx q[0];
rz(2.7713293) q[0];
rz(-0.26369357) q[1];
sx q[1];
rz(-1.5636684) q[1];
sx q[1];
rz(2.945074) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21393468) q[0];
sx q[0];
rz(-1.2688583) q[0];
sx q[0];
rz(1.2769804) q[0];
rz(-pi) q[1];
rz(-1.1097124) q[2];
sx q[2];
rz(-1.5421252) q[2];
sx q[2];
rz(-1.2210786) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6477485) q[1];
sx q[1];
rz(-1.7462247) q[1];
sx q[1];
rz(-0.26248787) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8653141) q[3];
sx q[3];
rz(-1.7344788) q[3];
sx q[3];
rz(-0.39556634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1281517) q[2];
sx q[2];
rz(-0.95169607) q[2];
sx q[2];
rz(-0.91709843) q[2];
rz(-2.3069978) q[3];
sx q[3];
rz(-1.83225) q[3];
sx q[3];
rz(-2.807907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.82814) q[0];
sx q[0];
rz(-1.7001292) q[0];
sx q[0];
rz(0.20203461) q[0];
rz(-2.0687885) q[1];
sx q[1];
rz(-1.5067116) q[1];
sx q[1];
rz(1.3938168) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48702792) q[0];
sx q[0];
rz(-1.9019039) q[0];
sx q[0];
rz(1.6513255) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2932111) q[2];
sx q[2];
rz(-2.7536014) q[2];
sx q[2];
rz(-0.1462305) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.60336941) q[1];
sx q[1];
rz(-1.1829611) q[1];
sx q[1];
rz(0.28660725) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4059541) q[3];
sx q[3];
rz(-1.9197575) q[3];
sx q[3];
rz(-2.5603106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1918734) q[2];
sx q[2];
rz(-1.4800025) q[2];
sx q[2];
rz(0.8423155) q[2];
rz(-2.0541644) q[3];
sx q[3];
rz(-0.73173404) q[3];
sx q[3];
rz(1.6585763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52952805) q[0];
sx q[0];
rz(-1.8719215) q[0];
sx q[0];
rz(-1.0213617) q[0];
rz(-1.1313324) q[1];
sx q[1];
rz(-2.1919057) q[1];
sx q[1];
rz(-2.9313415) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70049858) q[0];
sx q[0];
rz(-2.9919618) q[0];
sx q[0];
rz(-1.11087) q[0];
x q[1];
rz(0.16367754) q[2];
sx q[2];
rz(-1.5058766) q[2];
sx q[2];
rz(2.373873) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.47760689) q[1];
sx q[1];
rz(-0.37316445) q[1];
sx q[1];
rz(2.5928156) q[1];
rz(-2.2103569) q[3];
sx q[3];
rz(-0.82914549) q[3];
sx q[3];
rz(-0.88013807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3147099) q[2];
sx q[2];
rz(-1.3898536) q[2];
sx q[2];
rz(1.38114) q[2];
rz(-2.3794877) q[3];
sx q[3];
rz(-2.6657181) q[3];
sx q[3];
rz(0.41346082) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49522266) q[0];
sx q[0];
rz(-0.78283993) q[0];
sx q[0];
rz(0.58724171) q[0];
rz(-2.6953221) q[1];
sx q[1];
rz(-1.279436) q[1];
sx q[1];
rz(2.1648724) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92102345) q[0];
sx q[0];
rz(-2.0221634) q[0];
sx q[0];
rz(-1.4054106) q[0];
rz(-pi) q[1];
rz(1.4090528) q[2];
sx q[2];
rz(-1.4913034) q[2];
sx q[2];
rz(0.0083991945) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8752012) q[1];
sx q[1];
rz(-0.57386639) q[1];
sx q[1];
rz(2.2465474) q[1];
rz(-pi) q[2];
rz(-0.85985698) q[3];
sx q[3];
rz(-1.1338187) q[3];
sx q[3];
rz(0.45339282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7405159) q[2];
sx q[2];
rz(-0.70827168) q[2];
sx q[2];
rz(-2.5891499) q[2];
rz(-2.6756514) q[3];
sx q[3];
rz(-1.3876785) q[3];
sx q[3];
rz(-2.1140816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4733646) q[0];
sx q[0];
rz(-0.79376525) q[0];
sx q[0];
rz(2.2209432) q[0];
rz(3.0905511) q[1];
sx q[1];
rz(-1.5565245) q[1];
sx q[1];
rz(0.89909536) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9793648) q[0];
sx q[0];
rz(-2.3875036) q[0];
sx q[0];
rz(0.92569949) q[0];
rz(-0.33234889) q[2];
sx q[2];
rz(-2.1045448) q[2];
sx q[2];
rz(2.3975083) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2832665) q[1];
sx q[1];
rz(-0.35052931) q[1];
sx q[1];
rz(-1.2298898) q[1];
x q[2];
rz(1.5702463) q[3];
sx q[3];
rz(-0.37017248) q[3];
sx q[3];
rz(0.40035029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1961394) q[2];
sx q[2];
rz(-3.0710199) q[2];
sx q[2];
rz(3.0370039) q[2];
rz(0.83834046) q[3];
sx q[3];
rz(-2.392277) q[3];
sx q[3];
rz(2.2588363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8024682) q[0];
sx q[0];
rz(-0.81766468) q[0];
sx q[0];
rz(1.1178281) q[0];
rz(-2.4291908) q[1];
sx q[1];
rz(-0.63122216) q[1];
sx q[1];
rz(0.39696473) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3544281) q[0];
sx q[0];
rz(-2.0485326) q[0];
sx q[0];
rz(0.54425311) q[0];
rz(-pi) q[1];
rz(-2.8357387) q[2];
sx q[2];
rz(-1.9216188) q[2];
sx q[2];
rz(2.2985814) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9576479) q[1];
sx q[1];
rz(-0.34788222) q[1];
sx q[1];
rz(1.3089048) q[1];
x q[2];
rz(1.0730562) q[3];
sx q[3];
rz(-1.9314249) q[3];
sx q[3];
rz(-1.8325168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.94840702) q[2];
sx q[2];
rz(-2.2414175) q[2];
sx q[2];
rz(-3.0449384) q[2];
rz(-1.7276673) q[3];
sx q[3];
rz(-2.0035412) q[3];
sx q[3];
rz(-1.1269425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5861355) q[0];
sx q[0];
rz(-2.0270892) q[0];
sx q[0];
rz(-2.8591697) q[0];
rz(-1.2697521) q[1];
sx q[1];
rz(-1.3602263) q[1];
sx q[1];
rz(0.7855986) q[1];
rz(1.6554228) q[2];
sx q[2];
rz(-1.6783236) q[2];
sx q[2];
rz(2.8477737) q[2];
rz(1.9962911) q[3];
sx q[3];
rz(-2.8588061) q[3];
sx q[3];
rz(-2.2929946) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];