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
rz(-0.86369792) q[0];
sx q[0];
rz(0.20110826) q[0];
rz(1.7445298) q[1];
sx q[1];
rz(4.9464524) q[1];
sx q[1];
rz(10.051605) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62030828) q[0];
sx q[0];
rz(-1.4518019) q[0];
sx q[0];
rz(-1.760153) q[0];
rz(1.8148412) q[2];
sx q[2];
rz(-1.2284757) q[2];
sx q[2];
rz(0.92820456) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7658753) q[1];
sx q[1];
rz(-2.1581576) q[1];
sx q[1];
rz(-0.37382965) q[1];
rz(-pi) q[2];
rz(-1.0073644) q[3];
sx q[3];
rz(-2.8197188) q[3];
sx q[3];
rz(0.25062996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8864002) q[2];
sx q[2];
rz(-0.93078405) q[2];
sx q[2];
rz(1.8480776) q[2];
rz(2.0375997) q[3];
sx q[3];
rz(-0.84775001) q[3];
sx q[3];
rz(2.3867992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.8936546) q[0];
sx q[0];
rz(-1.0590483) q[0];
sx q[0];
rz(-2.1564116) q[0];
rz(-2.7424116) q[1];
sx q[1];
rz(-1.8789411) q[1];
sx q[1];
rz(0.10736297) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8132919) q[0];
sx q[0];
rz(-2.11587) q[0];
sx q[0];
rz(-2.6585447) q[0];
rz(-pi) q[1];
rz(0.58789247) q[2];
sx q[2];
rz(-1.8111472) q[2];
sx q[2];
rz(-0.27871486) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0641891) q[1];
sx q[1];
rz(-1.2014376) q[1];
sx q[1];
rz(-0.87537745) q[1];
x q[2];
rz(-1.739005) q[3];
sx q[3];
rz(-1.1871927) q[3];
sx q[3];
rz(-2.3080491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.62464109) q[2];
sx q[2];
rz(-1.7525643) q[2];
sx q[2];
rz(2.4798933) q[2];
rz(-0.055756904) q[3];
sx q[3];
rz(-2.7825833) q[3];
sx q[3];
rz(2.8957446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4645638) q[0];
sx q[0];
rz(-0.55389535) q[0];
sx q[0];
rz(-0.04034986) q[0];
rz(0.57998747) q[1];
sx q[1];
rz(-2.2433387) q[1];
sx q[1];
rz(-2.0592164) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28856453) q[0];
sx q[0];
rz(-1.0596501) q[0];
sx q[0];
rz(1.9938064) q[0];
x q[1];
rz(2.9897887) q[2];
sx q[2];
rz(-1.9215343) q[2];
sx q[2];
rz(1.1317859) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3861474) q[1];
sx q[1];
rz(-1.3324932) q[1];
sx q[1];
rz(-2.9134681) q[1];
x q[2];
rz(-2.1654487) q[3];
sx q[3];
rz(-2.757132) q[3];
sx q[3];
rz(2.773075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.13517705) q[2];
sx q[2];
rz(-1.2536896) q[2];
sx q[2];
rz(-2.5877) q[2];
rz(1.4931549) q[3];
sx q[3];
rz(-1.0529543) q[3];
sx q[3];
rz(-0.55251399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64788139) q[0];
sx q[0];
rz(-0.58409062) q[0];
sx q[0];
rz(0.035382263) q[0];
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
rz(-1.1195388) q[0];
sx q[0];
rz(-1.6484123) q[0];
sx q[0];
rz(2.0366497) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.79779063) q[2];
sx q[2];
rz(-0.75227037) q[2];
sx q[2];
rz(2.4385902) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0636053) q[1];
sx q[1];
rz(-1.6997153) q[1];
sx q[1];
rz(-1.414826) q[1];
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
rz(pi/2) q[1];
sx q[1];
rz(-0.35486832) q[2];
sx q[2];
rz(-1.0151261) q[2];
sx q[2];
rz(0.50746894) q[2];
rz(-1.9594225) q[3];
sx q[3];
rz(-1.2927262) q[3];
sx q[3];
rz(1.8866084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59947157) q[0];
sx q[0];
rz(-2.4276955) q[0];
sx q[0];
rz(-0.3702634) q[0];
rz(-0.26369357) q[1];
sx q[1];
rz(-1.5779243) q[1];
sx q[1];
rz(-2.945074) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6949961) q[0];
sx q[0];
rz(-1.2906404) q[0];
sx q[0];
rz(2.8269935) q[0];
x q[1];
rz(-3.1095805) q[2];
sx q[2];
rz(-1.1099166) q[2];
sx q[2];
rz(-0.33547685) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.12381183) q[1];
sx q[1];
rz(-1.31243) q[1];
sx q[1];
rz(-1.3892795) q[1];
x q[2];
rz(0.27627857) q[3];
sx q[3];
rz(-1.7344788) q[3];
sx q[3];
rz(-0.39556634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.013441) q[2];
sx q[2];
rz(-0.95169607) q[2];
sx q[2];
rz(0.91709843) q[2];
rz(0.83459485) q[3];
sx q[3];
rz(-1.3093427) q[3];
sx q[3];
rz(-0.33368567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.82814) q[0];
sx q[0];
rz(-1.7001292) q[0];
sx q[0];
rz(-0.20203461) q[0];
rz(2.0687885) q[1];
sx q[1];
rz(-1.5067116) q[1];
sx q[1];
rz(1.7477759) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6545647) q[0];
sx q[0];
rz(-1.9019039) q[0];
sx q[0];
rz(-1.4902671) q[0];
rz(-pi) q[1];
rz(-1.2932111) q[2];
sx q[2];
rz(-0.38799122) q[2];
sx q[2];
rz(-0.1462305) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5382232) q[1];
sx q[1];
rz(-1.1829611) q[1];
sx q[1];
rz(0.28660725) q[1];
rz(-2.4059541) q[3];
sx q[3];
rz(-1.9197575) q[3];
sx q[3];
rz(2.5603106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1918734) q[2];
sx q[2];
rz(-1.6615901) q[2];
sx q[2];
rz(-0.8423155) q[2];
rz(-1.0874282) q[3];
sx q[3];
rz(-0.73173404) q[3];
sx q[3];
rz(1.4830164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
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
rz(0.21025118) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70049858) q[0];
sx q[0];
rz(-0.14963089) q[0];
sx q[0];
rz(2.0307226) q[0];
x q[1];
rz(0.16367754) q[2];
sx q[2];
rz(-1.6357161) q[2];
sx q[2];
rz(-2.373873) q[2];
rz(-pi) q[3];
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
rz(2.5928156) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.85150163) q[3];
sx q[3];
rz(-1.1151033) q[3];
sx q[3];
rz(1.9853026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.82688275) q[2];
sx q[2];
rz(-1.3898536) q[2];
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
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(0.49522266) q[0];
sx q[0];
rz(-0.78283993) q[0];
sx q[0];
rz(0.58724171) q[0];
rz(-0.44627055) q[1];
sx q[1];
rz(-1.8621567) q[1];
sx q[1];
rz(-0.97672021) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8551089) q[0];
sx q[0];
rz(-2.6628462) q[0];
sx q[0];
rz(2.8141862) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1114242) q[2];
sx q[2];
rz(-0.18006912) q[2];
sx q[2];
rz(2.0153231) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2446489) q[1];
sx q[1];
rz(-1.2243425) q[1];
sx q[1];
rz(-1.1036554) q[1];
rz(-pi) q[2];
rz(2.1920491) q[3];
sx q[3];
rz(-2.3275259) q[3];
sx q[3];
rz(-2.4809588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4010767) q[2];
sx q[2];
rz(-0.70827168) q[2];
sx q[2];
rz(-2.5891499) q[2];
rz(2.6756514) q[3];
sx q[3];
rz(-1.3876785) q[3];
sx q[3];
rz(-1.027511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66822806) q[0];
sx q[0];
rz(-0.79376525) q[0];
sx q[0];
rz(0.92064944) q[0];
rz(-0.051041516) q[1];
sx q[1];
rz(-1.5565245) q[1];
sx q[1];
rz(-2.2424973) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5024922) q[0];
sx q[0];
rz(-2.149625) q[0];
sx q[0];
rz(-0.51410189) q[0];
rz(-2.0752418) q[2];
sx q[2];
rz(-0.62014183) q[2];
sx q[2];
rz(-1.3401741) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2832665) q[1];
sx q[1];
rz(-2.7910633) q[1];
sx q[1];
rz(-1.9117029) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2006239) q[3];
sx q[3];
rz(-1.5705974) q[3];
sx q[3];
rz(1.1699333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.94545323) q[2];
sx q[2];
rz(-0.070572704) q[2];
sx q[2];
rz(-0.10458874) q[2];
rz(-0.83834046) q[3];
sx q[3];
rz(-2.392277) q[3];
sx q[3];
rz(-2.2588363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8024682) q[0];
sx q[0];
rz(-0.81766468) q[0];
sx q[0];
rz(-1.1178281) q[0];
rz(-2.4291908) q[1];
sx q[1];
rz(-0.63122216) q[1];
sx q[1];
rz(0.39696473) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8660276) q[0];
sx q[0];
rz(-2.4336928) q[0];
sx q[0];
rz(-0.78535725) q[0];
rz(-0.30585395) q[2];
sx q[2];
rz(-1.2199739) q[2];
sx q[2];
rz(-0.84301126) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.18394477) q[1];
sx q[1];
rz(-2.7937104) q[1];
sx q[1];
rz(-1.3089048) q[1];
x q[2];
rz(1.0730562) q[3];
sx q[3];
rz(-1.9314249) q[3];
sx q[3];
rz(1.3090759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1931856) q[2];
sx q[2];
rz(-0.90017515) q[2];
sx q[2];
rz(-0.096654264) q[2];
rz(1.4139253) q[3];
sx q[3];
rz(-2.0035412) q[3];
sx q[3];
rz(-1.1269425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5861355) q[0];
sx q[0];
rz(-2.0270892) q[0];
sx q[0];
rz(-2.8591697) q[0];
rz(-1.8718406) q[1];
sx q[1];
rz(-1.7813663) q[1];
sx q[1];
rz(-2.355994) q[1];
rz(3.0336822) q[2];
sx q[2];
rz(-1.6549329) q[2];
sx q[2];
rz(-1.8555117) q[2];
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
