OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.1768567) q[0];
sx q[0];
rz(5.4194874) q[0];
sx q[0];
rz(9.6258862) q[0];
rz(1.7445298) q[1];
sx q[1];
rz(-1.3367329) q[1];
sx q[1];
rz(0.62682682) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5212844) q[0];
sx q[0];
rz(-1.4518019) q[0];
sx q[0];
rz(-1.760153) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7896499) q[2];
sx q[2];
rz(-1.3411739) q[2];
sx q[2];
rz(2.4156092) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.99170291) q[1];
sx q[1];
rz(-2.457379) q[1];
sx q[1];
rz(1.0690772) q[1];
x q[2];
rz(1.8455891) q[3];
sx q[3];
rz(-1.7405675) q[3];
sx q[3];
rz(2.3613289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.25519249) q[2];
sx q[2];
rz(-2.2108086) q[2];
sx q[2];
rz(1.8480776) q[2];
rz(-1.1039929) q[3];
sx q[3];
rz(-0.84775001) q[3];
sx q[3];
rz(2.3867992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24793808) q[0];
sx q[0];
rz(-2.0825443) q[0];
sx q[0];
rz(-2.1564116) q[0];
rz(-2.7424116) q[1];
sx q[1];
rz(-1.8789411) q[1];
sx q[1];
rz(-3.0342297) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1185703) q[0];
sx q[0];
rz(-1.9792299) q[0];
sx q[0];
rz(0.9704216) q[0];
rz(-pi) q[1];
x q[1];
rz(0.58789247) q[2];
sx q[2];
rz(-1.3304454) q[2];
sx q[2];
rz(0.27871486) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3423691) q[1];
sx q[1];
rz(-2.2111726) q[1];
sx q[1];
rz(0.46701042) q[1];
rz(-1.4025877) q[3];
sx q[3];
rz(-1.1871927) q[3];
sx q[3];
rz(-0.83354359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5169516) q[2];
sx q[2];
rz(-1.7525643) q[2];
sx q[2];
rz(0.66169935) q[2];
rz(0.055756904) q[3];
sx q[3];
rz(-0.35900933) q[3];
sx q[3];
rz(-0.24584809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67702883) q[0];
sx q[0];
rz(-2.5876973) q[0];
sx q[0];
rz(3.1012428) q[0];
rz(-0.57998747) q[1];
sx q[1];
rz(-2.2433387) q[1];
sx q[1];
rz(-1.0823762) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28856453) q[0];
sx q[0];
rz(-2.0819426) q[0];
sx q[0];
rz(1.9938064) q[0];
rz(-1.9627377) q[2];
sx q[2];
rz(-0.38092962) q[2];
sx q[2];
rz(0.71289635) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7554453) q[1];
sx q[1];
rz(-1.8090994) q[1];
sx q[1];
rz(2.9134681) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9186967) q[3];
sx q[3];
rz(-1.8867023) q[3];
sx q[3];
rz(0.99881682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.13517705) q[2];
sx q[2];
rz(-1.2536896) q[2];
sx q[2];
rz(2.5877) q[2];
rz(-1.4931549) q[3];
sx q[3];
rz(-2.0886383) q[3];
sx q[3];
rz(-0.55251399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4937113) q[0];
sx q[0];
rz(-2.557502) q[0];
sx q[0];
rz(-3.1062104) q[0];
rz(2.1454051) q[1];
sx q[1];
rz(-1.532998) q[1];
sx q[1];
rz(2.6534973) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49022084) q[0];
sx q[0];
rz(-1.106456) q[0];
sx q[0];
rz(-0.086829348) q[0];
rz(0.98056356) q[2];
sx q[2];
rz(-2.0681941) q[2];
sx q[2];
rz(-2.892708) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6285703) q[1];
sx q[1];
rz(-1.4161308) q[1];
sx q[1];
rz(-3.0111074) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.475708) q[3];
sx q[3];
rz(-0.69722359) q[3];
sx q[3];
rz(-2.9149885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7867243) q[2];
sx q[2];
rz(-2.1264666) q[2];
sx q[2];
rz(-2.6341237) q[2];
rz(-1.1821702) q[3];
sx q[3];
rz(-1.2927262) q[3];
sx q[3];
rz(-1.8866084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59947157) q[0];
sx q[0];
rz(-0.71389714) q[0];
sx q[0];
rz(2.7713293) q[0];
rz(2.8778991) q[1];
sx q[1];
rz(-1.5636684) q[1];
sx q[1];
rz(2.945074) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4465966) q[0];
sx q[0];
rz(-1.2906404) q[0];
sx q[0];
rz(-2.8269935) q[0];
rz(1.1097124) q[2];
sx q[2];
rz(-1.5421252) q[2];
sx q[2];
rz(1.2210786) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0177808) q[1];
sx q[1];
rz(-1.8291626) q[1];
sx q[1];
rz(-1.7523132) q[1];
rz(-pi) q[2];
rz(-0.27627857) q[3];
sx q[3];
rz(-1.7344788) q[3];
sx q[3];
rz(0.39556634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.013441) q[2];
sx q[2];
rz(-0.95169607) q[2];
sx q[2];
rz(0.91709843) q[2];
rz(-0.83459485) q[3];
sx q[3];
rz(-1.3093427) q[3];
sx q[3];
rz(-2.807907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.82814) q[0];
sx q[0];
rz(-1.7001292) q[0];
sx q[0];
rz(-2.939558) q[0];
rz(-1.0728041) q[1];
sx q[1];
rz(-1.5067116) q[1];
sx q[1];
rz(-1.3938168) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0840543) q[0];
sx q[0];
rz(-1.6469427) q[0];
sx q[0];
rz(0.3321068) q[0];
x q[1];
rz(1.9453085) q[2];
sx q[2];
rz(-1.4669344) q[2];
sx q[2];
rz(1.166677) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8762293) q[1];
sx q[1];
rz(-0.47788844) q[1];
sx q[1];
rz(-0.96546121) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6447547) q[3];
sx q[3];
rz(-2.3416069) q[3];
sx q[3];
rz(0.62832384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.94971925) q[2];
sx q[2];
rz(-1.6615901) q[2];
sx q[2];
rz(0.8423155) q[2];
rz(2.0541644) q[3];
sx q[3];
rz(-2.4098586) q[3];
sx q[3];
rz(1.6585763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52952805) q[0];
sx q[0];
rz(-1.2696711) q[0];
sx q[0];
rz(1.0213617) q[0];
rz(2.0102603) q[1];
sx q[1];
rz(-2.1919057) q[1];
sx q[1];
rz(-2.9313415) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7267669) q[0];
sx q[0];
rz(-1.504577) q[0];
sx q[0];
rz(1.7050752) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5049997) q[2];
sx q[2];
rz(-1.4074667) q[2];
sx q[2];
rz(-0.79236275) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6639858) q[1];
sx q[1];
rz(-2.7684282) q[1];
sx q[1];
rz(-0.54877703) q[1];
rz(-pi) q[2];
rz(-0.85150163) q[3];
sx q[3];
rz(-2.0264894) q[3];
sx q[3];
rz(1.15629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.82688275) q[2];
sx q[2];
rz(-1.751739) q[2];
sx q[2];
rz(1.38114) q[2];
rz(0.76210493) q[3];
sx q[3];
rz(-0.47587454) q[3];
sx q[3];
rz(-0.41346082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49522266) q[0];
sx q[0];
rz(-2.3587527) q[0];
sx q[0];
rz(0.58724171) q[0];
rz(0.44627055) q[1];
sx q[1];
rz(-1.8621567) q[1];
sx q[1];
rz(0.97672021) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57709549) q[0];
sx q[0];
rz(-1.4221039) q[0];
sx q[0];
rz(-0.45678267) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7325399) q[2];
sx q[2];
rz(-1.4913034) q[2];
sx q[2];
rz(3.1331935) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2446489) q[1];
sx q[1];
rz(-1.9172501) q[1];
sx q[1];
rz(2.0379373) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2817357) q[3];
sx q[3];
rz(-2.007774) q[3];
sx q[3];
rz(-0.45339282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4010767) q[2];
sx q[2];
rz(-0.70827168) q[2];
sx q[2];
rz(-0.55244279) q[2];
rz(2.6756514) q[3];
sx q[3];
rz(-1.7539141) q[3];
sx q[3];
rz(1.027511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4733646) q[0];
sx q[0];
rz(-2.3478274) q[0];
sx q[0];
rz(2.2209432) q[0];
rz(3.0905511) q[1];
sx q[1];
rz(-1.5850681) q[1];
sx q[1];
rz(2.2424973) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9793648) q[0];
sx q[0];
rz(-2.3875036) q[0];
sx q[0];
rz(2.2158932) q[0];
rz(-pi) q[1];
rz(-2.0752418) q[2];
sx q[2];
rz(-0.62014183) q[2];
sx q[2];
rz(1.8014185) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.92210772) q[1];
sx q[1];
rz(-1.9003632) q[1];
sx q[1];
rz(-0.12164128) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5702463) q[3];
sx q[3];
rz(-0.37017248) q[3];
sx q[3];
rz(2.7412424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1961394) q[2];
sx q[2];
rz(-0.070572704) q[2];
sx q[2];
rz(-0.10458874) q[2];
rz(-0.83834046) q[3];
sx q[3];
rz(-0.74931562) q[3];
sx q[3];
rz(-0.88275638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8024682) q[0];
sx q[0];
rz(-0.81766468) q[0];
sx q[0];
rz(2.0237645) q[0];
rz(-0.71240187) q[1];
sx q[1];
rz(-2.5103705) q[1];
sx q[1];
rz(0.39696473) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3544281) q[0];
sx q[0];
rz(-1.0930601) q[0];
sx q[0];
rz(0.54425311) q[0];
rz(-pi) q[1];
x q[1];
rz(0.88231477) q[2];
sx q[2];
rz(-2.6803662) q[2];
sx q[2];
rz(-0.099909401) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.46170235) q[1];
sx q[1];
rz(-1.2352518) q[1];
sx q[1];
rz(-0.093613503) q[1];
rz(-pi) q[2];
rz(-2.7361717) q[3];
sx q[3];
rz(-1.1076895) q[3];
sx q[3];
rz(0.4511569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1931856) q[2];
sx q[2];
rz(-0.90017515) q[2];
sx q[2];
rz(-0.096654264) q[2];
rz(-1.4139253) q[3];
sx q[3];
rz(-1.1380514) q[3];
sx q[3];
rz(-1.1269425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5554572) q[0];
sx q[0];
rz(-1.1145034) q[0];
sx q[0];
rz(0.282423) q[0];
rz(-1.8718406) q[1];
sx q[1];
rz(-1.7813663) q[1];
sx q[1];
rz(-2.355994) q[1];
rz(0.66432129) q[2];
sx q[2];
rz(-0.13673377) q[2];
sx q[2];
rz(-2.7665334) q[2];
rz(-1.8295287) q[3];
sx q[3];
rz(-1.4553634) q[3];
sx q[3];
rz(2.8298557) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
