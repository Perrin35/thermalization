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
rz(2.3020701) q[0];
sx q[0];
rz(4.7293651) q[0];
sx q[0];
rz(11.601996) q[0];
rz(0.40274611) q[1];
sx q[1];
rz(-1.6778814) q[1];
sx q[1];
rz(-0.95199624) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.096199252) q[0];
sx q[0];
rz(-2.0598167) q[0];
sx q[0];
rz(3.0351991) q[0];
x q[1];
rz(1.9929041) q[2];
sx q[2];
rz(-1.4808345) q[2];
sx q[2];
rz(0.46101704) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.18033548) q[1];
sx q[1];
rz(-1.9587659) q[1];
sx q[1];
rz(-1.9684102) q[1];
rz(0.49523109) q[3];
sx q[3];
rz(-0.81988813) q[3];
sx q[3];
rz(-0.26671975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6160148) q[2];
sx q[2];
rz(-1.3495477) q[2];
sx q[2];
rz(0.40526059) q[2];
rz(2.0112093) q[3];
sx q[3];
rz(-0.78961343) q[3];
sx q[3];
rz(1.9270886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90527117) q[0];
sx q[0];
rz(-3.1089678) q[0];
sx q[0];
rz(0.76500934) q[0];
rz(2.8969823) q[1];
sx q[1];
rz(-1.2449934) q[1];
sx q[1];
rz(2.5027221) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34341223) q[0];
sx q[0];
rz(-1.596178) q[0];
sx q[0];
rz(-2.4973386) q[0];
x q[1];
rz(-2.5193754) q[2];
sx q[2];
rz(-1.8832616) q[2];
sx q[2];
rz(2.0553596) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.42145211) q[1];
sx q[1];
rz(-1.1126592) q[1];
sx q[1];
rz(-1.1790183) q[1];
x q[2];
rz(1.8539393) q[3];
sx q[3];
rz(-0.75732175) q[3];
sx q[3];
rz(0.18268302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2834393) q[2];
sx q[2];
rz(-0.44698134) q[2];
sx q[2];
rz(-2.8972674) q[2];
rz(1.1089995) q[3];
sx q[3];
rz(-1.2764443) q[3];
sx q[3];
rz(0.23787704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7636488) q[0];
sx q[0];
rz(-2.0574594) q[0];
sx q[0];
rz(-1.9042683) q[0];
rz(1.7297277) q[1];
sx q[1];
rz(-0.4487764) q[1];
sx q[1];
rz(-1.809874) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6520572) q[0];
sx q[0];
rz(-2.3416828) q[0];
sx q[0];
rz(1.4481988) q[0];
rz(-pi) q[1];
rz(-1.9218512) q[2];
sx q[2];
rz(-1.7760385) q[2];
sx q[2];
rz(0.83265162) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6159042) q[1];
sx q[1];
rz(-1.5153023) q[1];
sx q[1];
rz(-1.5374119) q[1];
rz(-pi) q[2];
rz(0.94966763) q[3];
sx q[3];
rz(-2.0916208) q[3];
sx q[3];
rz(0.68062147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1639013) q[2];
sx q[2];
rz(-2.0745664) q[2];
sx q[2];
rz(-0.65995556) q[2];
rz(-1.6948949) q[3];
sx q[3];
rz(-1.4283254) q[3];
sx q[3];
rz(-2.0383294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8834943) q[0];
sx q[0];
rz(-0.75089184) q[0];
sx q[0];
rz(1.1335491) q[0];
rz(-2.6253888) q[1];
sx q[1];
rz(-1.8323003) q[1];
sx q[1];
rz(-0.54905999) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20401317) q[0];
sx q[0];
rz(-1.4268865) q[0];
sx q[0];
rz(-1.4043441) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1604161) q[2];
sx q[2];
rz(-1.7038888) q[2];
sx q[2];
rz(-2.1504064) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.258865) q[1];
sx q[1];
rz(-0.2828063) q[1];
sx q[1];
rz(-2.905719) q[1];
rz(0.28259773) q[3];
sx q[3];
rz(-1.9535086) q[3];
sx q[3];
rz(0.0569009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.05222008) q[2];
sx q[2];
rz(-1.0129656) q[2];
sx q[2];
rz(-0.29903665) q[2];
rz(-0.61797577) q[3];
sx q[3];
rz(-1.0733913) q[3];
sx q[3];
rz(-1.7834024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1842136) q[0];
sx q[0];
rz(-3.1035677) q[0];
sx q[0];
rz(2.6755565) q[0];
rz(-2.274463) q[1];
sx q[1];
rz(-0.51141089) q[1];
sx q[1];
rz(0.83494157) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4049618) q[0];
sx q[0];
rz(-1.5190655) q[0];
sx q[0];
rz(-1.542227) q[0];
rz(-1.6500128) q[2];
sx q[2];
rz(-1.4822373) q[2];
sx q[2];
rz(1.7065602) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9795436) q[1];
sx q[1];
rz(-1.0912553) q[1];
sx q[1];
rz(0.060258106) q[1];
x q[2];
rz(-0.68170516) q[3];
sx q[3];
rz(-2.2870664) q[3];
sx q[3];
rz(2.0264421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.20778188) q[2];
sx q[2];
rz(-0.89263478) q[2];
sx q[2];
rz(-3.1100698) q[2];
rz(2.3837714) q[3];
sx q[3];
rz(-1.1008681) q[3];
sx q[3];
rz(-2.5789564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9071478) q[0];
sx q[0];
rz(-1.7666768) q[0];
sx q[0];
rz(-3.1100519) q[0];
rz(3.0517598) q[1];
sx q[1];
rz(-2.641771) q[1];
sx q[1];
rz(2.8057742) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9851097) q[0];
sx q[0];
rz(-1.2439737) q[0];
sx q[0];
rz(0.37460297) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0903942) q[2];
sx q[2];
rz(-2.5246387) q[2];
sx q[2];
rz(-1.761957) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0586186) q[1];
sx q[1];
rz(-1.0320283) q[1];
sx q[1];
rz(-2.8584216) q[1];
rz(2.8434136) q[3];
sx q[3];
rz(-1.1345053) q[3];
sx q[3];
rz(2.0661708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1699367) q[2];
sx q[2];
rz(-1.1145096) q[2];
sx q[2];
rz(-1.2706832) q[2];
rz(3.0830248) q[3];
sx q[3];
rz(-1.4437557) q[3];
sx q[3];
rz(-0.34846714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5295277) q[0];
sx q[0];
rz(-0.53851524) q[0];
sx q[0];
rz(-2.8756323) q[0];
rz(2.3174441) q[1];
sx q[1];
rz(-1.8310841) q[1];
sx q[1];
rz(-1.5305653) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1743721) q[0];
sx q[0];
rz(-1.8215113) q[0];
sx q[0];
rz(0.27085568) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5179727) q[2];
sx q[2];
rz(-2.3011156) q[2];
sx q[2];
rz(-2.2234302) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2648523) q[1];
sx q[1];
rz(-1.0050259) q[1];
sx q[1];
rz(1.8840204) q[1];
rz(-0.77218382) q[3];
sx q[3];
rz(-1.8265542) q[3];
sx q[3];
rz(2.135475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1223569) q[2];
sx q[2];
rz(-1.1160206) q[2];
sx q[2];
rz(-2.1717211) q[2];
rz(0.14361778) q[3];
sx q[3];
rz(-0.93846455) q[3];
sx q[3];
rz(2.5041049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.403991) q[0];
sx q[0];
rz(-0.25594512) q[0];
sx q[0];
rz(-3.0373489) q[0];
rz(-0.741611) q[1];
sx q[1];
rz(-2.9545018) q[1];
sx q[1];
rz(2.7033477) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5271664) q[0];
sx q[0];
rz(-0.90133716) q[0];
sx q[0];
rz(-1.7015083) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.97317071) q[2];
sx q[2];
rz(-1.7665553) q[2];
sx q[2];
rz(0.51908609) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.91854535) q[1];
sx q[1];
rz(-2.1761736) q[1];
sx q[1];
rz(-0.041408509) q[1];
rz(1.7167818) q[3];
sx q[3];
rz(-2.101482) q[3];
sx q[3];
rz(-2.3455181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.40615842) q[2];
sx q[2];
rz(-0.27190748) q[2];
sx q[2];
rz(-2.820365) q[2];
rz(-0.99212232) q[3];
sx q[3];
rz(-2.0480053) q[3];
sx q[3];
rz(-1.7295624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.853249) q[0];
sx q[0];
rz(-0.11985954) q[0];
sx q[0];
rz(-2.993592) q[0];
rz(2.0026228) q[1];
sx q[1];
rz(-2.482246) q[1];
sx q[1];
rz(0.99050561) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97738701) q[0];
sx q[0];
rz(-1.6456104) q[0];
sx q[0];
rz(-2.6590024) q[0];
x q[1];
rz(2.4100346) q[2];
sx q[2];
rz(-1.9266945) q[2];
sx q[2];
rz(-1.6936091) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.77978342) q[1];
sx q[1];
rz(-1.6020007) q[1];
sx q[1];
rz(-0.0051910944) q[1];
rz(-pi) q[2];
x q[2];
rz(0.93659525) q[3];
sx q[3];
rz(-0.59118012) q[3];
sx q[3];
rz(-1.9612966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1691576) q[2];
sx q[2];
rz(-0.7908228) q[2];
sx q[2];
rz(2.5992744) q[2];
rz(-1.4646685) q[3];
sx q[3];
rz(-1.2277536) q[3];
sx q[3];
rz(0.98384583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-1.9024314) q[0];
sx q[0];
rz(-0.79817525) q[0];
sx q[0];
rz(-0.27957988) q[0];
rz(-2.3565049) q[1];
sx q[1];
rz(-0.79477349) q[1];
sx q[1];
rz(-0.40207544) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3668985) q[0];
sx q[0];
rz(-2.0089564) q[0];
sx q[0];
rz(1.9209981) q[0];
rz(-pi) q[1];
rz(-0.76519544) q[2];
sx q[2];
rz(-0.94586654) q[2];
sx q[2];
rz(-0.13002061) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6732503) q[1];
sx q[1];
rz(-1.3788917) q[1];
sx q[1];
rz(-0.90190355) q[1];
rz(-pi) q[2];
rz(-2.8567186) q[3];
sx q[3];
rz(-1.9683629) q[3];
sx q[3];
rz(-1.3314825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1606007) q[2];
sx q[2];
rz(-2.4805785) q[2];
sx q[2];
rz(0.92657363) q[2];
rz(2.5911234) q[3];
sx q[3];
rz(-0.87369839) q[3];
sx q[3];
rz(-1.8615104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4704623) q[0];
sx q[0];
rz(-1.067938) q[0];
sx q[0];
rz(-2.9965042) q[0];
rz(3.0628224) q[1];
sx q[1];
rz(-1.4046184) q[1];
sx q[1];
rz(1.2974993) q[1];
rz(1.5552022) q[2];
sx q[2];
rz(-0.77257665) q[2];
sx q[2];
rz(2.0333124) q[2];
rz(-0.13194247) q[3];
sx q[3];
rz(-2.0880825) q[3];
sx q[3];
rz(-1.1121398) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
