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
rz(-0.20110826) q[0];
rz(-1.3970628) q[1];
sx q[1];
rz(-1.8048598) q[1];
sx q[1];
rz(2.5147658) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1683567) q[0];
sx q[0];
rz(-1.7587979) q[0];
sx q[0];
rz(3.0204535) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5457382) q[2];
sx q[2];
rz(-2.7239954) q[2];
sx q[2];
rz(2.8516304) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1498897) q[1];
sx q[1];
rz(-0.68421364) q[1];
sx q[1];
rz(-1.0690772) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9653373) q[3];
sx q[3];
rz(-1.3000543) q[3];
sx q[3];
rz(0.83812974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.25519249) q[2];
sx q[2];
rz(-0.93078405) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8936546) q[0];
sx q[0];
rz(-1.0590483) q[0];
sx q[0];
rz(0.98518103) q[0];
rz(-0.39918104) q[1];
sx q[1];
rz(-1.8789411) q[1];
sx q[1];
rz(3.0342297) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1185703) q[0];
sx q[0];
rz(-1.1623628) q[0];
sx q[0];
rz(-2.1711711) q[0];
x q[1];
rz(1.2843578) q[2];
sx q[2];
rz(-2.1396494) q[2];
sx q[2];
rz(-1.6921649) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0641891) q[1];
sx q[1];
rz(-1.2014376) q[1];
sx q[1];
rz(-0.87537745) q[1];
rz(-pi) q[2];
rz(2.7483838) q[3];
sx q[3];
rz(-2.7244096) q[3];
sx q[3];
rz(1.259491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.62464109) q[2];
sx q[2];
rz(-1.3890283) q[2];
sx q[2];
rz(-0.66169935) q[2];
rz(3.0858357) q[3];
sx q[3];
rz(-2.7825833) q[3];
sx q[3];
rz(2.8957446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67702883) q[0];
sx q[0];
rz(-2.5876973) q[0];
sx q[0];
rz(-3.1012428) q[0];
rz(2.5616052) q[1];
sx q[1];
rz(-2.2433387) q[1];
sx q[1];
rz(-1.0823762) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6425991) q[0];
sx q[0];
rz(-1.2046308) q[0];
sx q[0];
rz(-0.55143349) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9252831) q[2];
sx q[2];
rz(-1.4282994) q[2];
sx q[2];
rz(-2.6500677) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.12990002) q[1];
sx q[1];
rz(-1.7923647) q[1];
sx q[1];
rz(1.8151912) q[1];
rz(-pi) q[2];
rz(-0.2228959) q[3];
sx q[3];
rz(-1.2548903) q[3];
sx q[3];
rz(-0.99881682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0064156) q[2];
sx q[2];
rz(-1.8879031) q[2];
sx q[2];
rz(2.5877) q[2];
rz(1.4931549) q[3];
sx q[3];
rz(-2.0886383) q[3];
sx q[3];
rz(-2.5890787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4937113) q[0];
sx q[0];
rz(-2.557502) q[0];
sx q[0];
rz(-0.035382263) q[0];
rz(2.1454051) q[1];
sx q[1];
rz(-1.6085947) q[1];
sx q[1];
rz(-2.6534973) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6513718) q[0];
sx q[0];
rz(-2.0351366) q[0];
sx q[0];
rz(-3.0547633) q[0];
x q[1];
rz(-2.1610291) q[2];
sx q[2];
rz(-2.0681941) q[2];
sx q[2];
rz(0.24888466) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.19255895) q[1];
sx q[1];
rz(-0.20201905) q[1];
sx q[1];
rz(0.87534027) q[1];
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
rz(-pi) q[1];
rz(-2.7867243) q[2];
sx q[2];
rz(-2.1264666) q[2];
sx q[2];
rz(0.50746894) q[2];
rz(-1.9594225) q[3];
sx q[3];
rz(-1.8488665) q[3];
sx q[3];
rz(1.2549843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5421211) q[0];
sx q[0];
rz(-0.71389714) q[0];
sx q[0];
rz(-2.7713293) q[0];
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
rz(1.6949961) q[0];
sx q[0];
rz(-1.8509522) q[0];
sx q[0];
rz(-2.8269935) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1095805) q[2];
sx q[2];
rz(-1.1099166) q[2];
sx q[2];
rz(0.33547685) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6424751) q[1];
sx q[1];
rz(-0.31458464) q[1];
sx q[1];
rz(-0.59928526) q[1];
x q[2];
rz(2.8653141) q[3];
sx q[3];
rz(-1.4071138) q[3];
sx q[3];
rz(-0.39556634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1281517) q[2];
sx q[2];
rz(-2.1898966) q[2];
sx q[2];
rz(0.91709843) q[2];
rz(0.83459485) q[3];
sx q[3];
rz(-1.3093427) q[3];
sx q[3];
rz(2.807907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.82814) q[0];
sx q[0];
rz(-1.4414635) q[0];
sx q[0];
rz(-0.20203461) q[0];
rz(-1.0728041) q[1];
sx q[1];
rz(-1.5067116) q[1];
sx q[1];
rz(1.7477759) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73035875) q[0];
sx q[0];
rz(-0.3404091) q[0];
sx q[0];
rz(0.22986869) q[0];
rz(-pi) q[1];
rz(0.11153531) q[2];
sx q[2];
rz(-1.9431912) q[2];
sx q[2];
rz(-0.44484777) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8762293) q[1];
sx q[1];
rz(-0.47788844) q[1];
sx q[1];
rz(0.96546121) q[1];
x q[2];
rz(2.0270258) q[3];
sx q[3];
rz(-2.2531415) q[3];
sx q[3];
rz(-1.2896468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1918734) q[2];
sx q[2];
rz(-1.4800025) q[2];
sx q[2];
rz(0.8423155) q[2];
rz(-2.0541644) q[3];
sx q[3];
rz(-2.4098586) q[3];
sx q[3];
rz(1.4830164) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6120646) q[0];
sx q[0];
rz(-1.8719215) q[0];
sx q[0];
rz(1.0213617) q[0];
rz(-1.1313324) q[1];
sx q[1];
rz(-0.94968692) q[1];
sx q[1];
rz(2.9313415) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7267669) q[0];
sx q[0];
rz(-1.504577) q[0];
sx q[0];
rz(1.4365175) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9779151) q[2];
sx q[2];
rz(-1.6357161) q[2];
sx q[2];
rz(2.373873) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6107632) q[1];
sx q[1];
rz(-1.7621343) q[1];
sx q[1];
rz(2.8192239) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93123575) q[3];
sx q[3];
rz(-0.82914549) q[3];
sx q[3];
rz(0.88013807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3147099) q[2];
sx q[2];
rz(-1.3898536) q[2];
sx q[2];
rz(-1.7604527) q[2];
rz(-0.76210493) q[3];
sx q[3];
rz(-0.47587454) q[3];
sx q[3];
rz(0.41346082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.64637) q[0];
sx q[0];
rz(-2.3587527) q[0];
sx q[0];
rz(2.5543509) q[0];
rz(-0.44627055) q[1];
sx q[1];
rz(-1.279436) q[1];
sx q[1];
rz(0.97672021) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8551089) q[0];
sx q[0];
rz(-0.47874641) q[0];
sx q[0];
rz(0.32740645) q[0];
rz(-1.7325399) q[2];
sx q[2];
rz(-1.4913034) q[2];
sx q[2];
rz(0.0083991945) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2663914) q[1];
sx q[1];
rz(-0.57386639) q[1];
sx q[1];
rz(0.89504524) q[1];
x q[2];
rz(-2.1920491) q[3];
sx q[3];
rz(-2.3275259) q[3];
sx q[3];
rz(2.4809588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7405159) q[2];
sx q[2];
rz(-0.70827168) q[2];
sx q[2];
rz(-0.55244279) q[2];
rz(0.46594122) q[3];
sx q[3];
rz(-1.7539141) q[3];
sx q[3];
rz(-1.027511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4733646) q[0];
sx q[0];
rz(-0.79376525) q[0];
sx q[0];
rz(2.2209432) q[0];
rz(-3.0905511) q[1];
sx q[1];
rz(-1.5565245) q[1];
sx q[1];
rz(-0.89909536) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16222787) q[0];
sx q[0];
rz(-2.3875036) q[0];
sx q[0];
rz(2.2158932) q[0];
rz(1.012072) q[2];
sx q[2];
rz(-1.2860824) q[2];
sx q[2];
rz(-2.4887091) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8583262) q[1];
sx q[1];
rz(-2.7910633) q[1];
sx q[1];
rz(1.9117029) q[1];
rz(-3.1413792) q[3];
sx q[3];
rz(-1.2006239) q[3];
sx q[3];
rz(-2.7406524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.94545323) q[2];
sx q[2];
rz(-0.070572704) q[2];
sx q[2];
rz(-0.10458874) q[2];
rz(0.83834046) q[3];
sx q[3];
rz(-0.74931562) q[3];
sx q[3];
rz(0.88275638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8024682) q[0];
sx q[0];
rz(-2.323928) q[0];
sx q[0];
rz(-1.1178281) q[0];
rz(0.71240187) q[1];
sx q[1];
rz(-0.63122216) q[1];
sx q[1];
rz(-2.7446279) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94496942) q[0];
sx q[0];
rz(-2.0485749) q[0];
sx q[0];
rz(2.1150132) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.88231477) q[2];
sx q[2];
rz(-2.6803662) q[2];
sx q[2];
rz(0.099909401) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6798903) q[1];
sx q[1];
rz(-1.2352518) q[1];
sx q[1];
rz(-3.0479792) q[1];
rz(-pi) q[2];
rz(2.0685365) q[3];
sx q[3];
rz(-1.9314249) q[3];
sx q[3];
rz(1.8325168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.94840702) q[2];
sx q[2];
rz(-2.2414175) q[2];
sx q[2];
rz(-0.096654264) q[2];
rz(-1.4139253) q[3];
sx q[3];
rz(-2.0035412) q[3];
sx q[3];
rz(1.1269425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5861355) q[0];
sx q[0];
rz(-1.1145034) q[0];
sx q[0];
rz(0.282423) q[0];
rz(-1.2697521) q[1];
sx q[1];
rz(-1.3602263) q[1];
sx q[1];
rz(0.7855986) q[1];
rz(-3.0336822) q[2];
sx q[2];
rz(-1.4866598) q[2];
sx q[2];
rz(1.286081) q[2];
rz(-0.11937033) q[3];
sx q[3];
rz(-1.8277677) q[3];
sx q[3];
rz(1.2895332) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
