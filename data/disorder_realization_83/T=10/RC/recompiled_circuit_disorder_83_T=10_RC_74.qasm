OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5542334) q[0];
sx q[0];
rz(-2.1589307) q[0];
sx q[0];
rz(0.76173705) q[0];
rz(-2.3770483) q[1];
sx q[1];
rz(-1.0772871) q[1];
sx q[1];
rz(2.397937) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0775454) q[0];
sx q[0];
rz(-0.24928688) q[0];
sx q[0];
rz(-1.3011978) q[0];
x q[1];
rz(1.450199) q[2];
sx q[2];
rz(-2.2199124) q[2];
sx q[2];
rz(-2.803034) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.57063369) q[1];
sx q[1];
rz(-1.966241) q[1];
sx q[1];
rz(-0.30084893) q[1];
x q[2];
rz(-0.071804382) q[3];
sx q[3];
rz(-1.2076488) q[3];
sx q[3];
rz(0.24955173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7068229) q[2];
sx q[2];
rz(-0.40199026) q[2];
sx q[2];
rz(-3.0337231) q[2];
rz(0.14262959) q[3];
sx q[3];
rz(-1.4167891) q[3];
sx q[3];
rz(-0.67255783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0649439) q[0];
sx q[0];
rz(-2.3363484) q[0];
sx q[0];
rz(2.9108677) q[0];
rz(1.327286) q[1];
sx q[1];
rz(-0.67115152) q[1];
sx q[1];
rz(3.1006295) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36386585) q[0];
sx q[0];
rz(-1.7330609) q[0];
sx q[0];
rz(-1.8198387) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3054753) q[2];
sx q[2];
rz(-2.7364519) q[2];
sx q[2];
rz(-0.9678313) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4119271) q[1];
sx q[1];
rz(-1.1150868) q[1];
sx q[1];
rz(2.5018442) q[1];
rz(-pi) q[2];
rz(-2.326194) q[3];
sx q[3];
rz(-2.6147463) q[3];
sx q[3];
rz(2.6032053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.162398) q[2];
sx q[2];
rz(-1.4731864) q[2];
sx q[2];
rz(0.56817788) q[2];
rz(-0.5125106) q[3];
sx q[3];
rz(-0.5957225) q[3];
sx q[3];
rz(-2.5797243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.7115962) q[0];
sx q[0];
rz(-1.6226409) q[0];
sx q[0];
rz(-0.45021737) q[0];
rz(-1.2954767) q[1];
sx q[1];
rz(-1.9772915) q[1];
sx q[1];
rz(-0.67726642) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60904658) q[0];
sx q[0];
rz(-0.17621528) q[0];
sx q[0];
rz(0.30058582) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0063653221) q[2];
sx q[2];
rz(-1.5718939) q[2];
sx q[2];
rz(-1.416625) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9159569) q[1];
sx q[1];
rz(-1.3924053) q[1];
sx q[1];
rz(-0.13326463) q[1];
rz(-pi) q[2];
rz(-0.47795313) q[3];
sx q[3];
rz(-1.9216929) q[3];
sx q[3];
rz(1.705845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9230817) q[2];
sx q[2];
rz(-1.1542902) q[2];
sx q[2];
rz(0.046860524) q[2];
rz(0.81165195) q[3];
sx q[3];
rz(-2.4620158) q[3];
sx q[3];
rz(2.9714382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1470404) q[0];
sx q[0];
rz(-3.0299598) q[0];
sx q[0];
rz(3.0901093) q[0];
rz(0.4908081) q[1];
sx q[1];
rz(-2.1665116) q[1];
sx q[1];
rz(1.1725918) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9358312) q[0];
sx q[0];
rz(-1.4371705) q[0];
sx q[0];
rz(2.4806116) q[0];
rz(-1.5863717) q[2];
sx q[2];
rz(-2.1767463) q[2];
sx q[2];
rz(1.7473999) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6231411) q[1];
sx q[1];
rz(-1.8703096) q[1];
sx q[1];
rz(1.2007984) q[1];
rz(-pi) q[2];
rz(3.067279) q[3];
sx q[3];
rz(-1.1399817) q[3];
sx q[3];
rz(-2.5527692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2525758) q[2];
sx q[2];
rz(-2.2106407) q[2];
sx q[2];
rz(0.48689294) q[2];
rz(1.9481109) q[3];
sx q[3];
rz(-0.32961696) q[3];
sx q[3];
rz(-0.025432767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.837773) q[0];
sx q[0];
rz(-1.0205512) q[0];
sx q[0];
rz(1.3254962) q[0];
rz(2.5505113) q[1];
sx q[1];
rz(-2.4609844) q[1];
sx q[1];
rz(-0.65471929) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7514873) q[0];
sx q[0];
rz(-2.562398) q[0];
sx q[0];
rz(0.67436995) q[0];
rz(-pi) q[1];
rz(-2.7227313) q[2];
sx q[2];
rz(-1.3992116) q[2];
sx q[2];
rz(0.78531839) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3760738) q[1];
sx q[1];
rz(-0.79080938) q[1];
sx q[1];
rz(-2.5942624) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4773024) q[3];
sx q[3];
rz(-1.1870541) q[3];
sx q[3];
rz(1.547471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4513662) q[2];
sx q[2];
rz(-1.9407242) q[2];
sx q[2];
rz(-0.5711242) q[2];
rz(-0.58978224) q[3];
sx q[3];
rz(-0.46001205) q[3];
sx q[3];
rz(3.0736198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-0.94474435) q[0];
sx q[0];
rz(-1.1905043) q[0];
sx q[0];
rz(0.29904547) q[0];
rz(-1.3202745) q[1];
sx q[1];
rz(-2.8838005) q[1];
sx q[1];
rz(-1.6437795) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2346674) q[0];
sx q[0];
rz(-1.5186131) q[0];
sx q[0];
rz(0.90953565) q[0];
x q[1];
rz(0.80070337) q[2];
sx q[2];
rz(-0.75468894) q[2];
sx q[2];
rz(-0.18037361) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1491579) q[1];
sx q[1];
rz(-0.5545534) q[1];
sx q[1];
rz(-1.0013594) q[1];
rz(-1.8201581) q[3];
sx q[3];
rz(-1.8421838) q[3];
sx q[3];
rz(-2.4621778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.51612878) q[2];
sx q[2];
rz(-1.7136145) q[2];
sx q[2];
rz(0.027475474) q[2];
rz(-2.6190858) q[3];
sx q[3];
rz(-2.3404739) q[3];
sx q[3];
rz(0.81645042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7291173) q[0];
sx q[0];
rz(-1.1183879) q[0];
sx q[0];
rz(3.0274042) q[0];
rz(0.97822899) q[1];
sx q[1];
rz(-0.54662919) q[1];
sx q[1];
rz(2.3506929) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7366911) q[0];
sx q[0];
rz(-1.4762523) q[0];
sx q[0];
rz(-0.94876429) q[0];
x q[1];
rz(2.3979264) q[2];
sx q[2];
rz(-1.3197834) q[2];
sx q[2];
rz(-2.2199092) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6939047) q[1];
sx q[1];
rz(-1.4920007) q[1];
sx q[1];
rz(1.9359971) q[1];
rz(-pi) q[2];
rz(2.0049719) q[3];
sx q[3];
rz(-0.73528157) q[3];
sx q[3];
rz(-2.6649464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2699282) q[2];
sx q[2];
rz(-2.0285138) q[2];
sx q[2];
rz(0.86501914) q[2];
rz(-0.67251742) q[3];
sx q[3];
rz(-1.1285684) q[3];
sx q[3];
rz(-3.045936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4928116) q[0];
sx q[0];
rz(-0.5195986) q[0];
sx q[0];
rz(-0.46736026) q[0];
rz(2.6043747) q[1];
sx q[1];
rz(-2.1610114) q[1];
sx q[1];
rz(-2.8875202) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5131322) q[0];
sx q[0];
rz(-1.5043133) q[0];
sx q[0];
rz(-2.9456338) q[0];
rz(2.8587564) q[2];
sx q[2];
rz(-3.1361702) q[2];
sx q[2];
rz(-1.1495513) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.2368187) q[1];
sx q[1];
rz(-0.33729759) q[1];
sx q[1];
rz(0.39223139) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0244175) q[3];
sx q[3];
rz(-2.5391038) q[3];
sx q[3];
rz(0.52091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.49999985) q[2];
sx q[2];
rz(-0.77074146) q[2];
sx q[2];
rz(2.8059778) q[2];
rz(-2.8149758) q[3];
sx q[3];
rz(-2.2461522) q[3];
sx q[3];
rz(-1.7232822) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.083387233) q[0];
sx q[0];
rz(-2.931262) q[0];
sx q[0];
rz(-1.0539508) q[0];
rz(0.15696934) q[1];
sx q[1];
rz(-1.418768) q[1];
sx q[1];
rz(-0.98186791) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40030038) q[0];
sx q[0];
rz(-1.9809082) q[0];
sx q[0];
rz(-1.6923231) q[0];
x q[1];
rz(-1.4023196) q[2];
sx q[2];
rz(-0.57541621) q[2];
sx q[2];
rz(1.1489431) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.30584221) q[1];
sx q[1];
rz(-1.3853449) q[1];
sx q[1];
rz(3.0830543) q[1];
rz(-2.8544159) q[3];
sx q[3];
rz(-2.4099726) q[3];
sx q[3];
rz(-1.8457796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2173569) q[2];
sx q[2];
rz(-1.057426) q[2];
sx q[2];
rz(-2.1441933) q[2];
rz(2.635397) q[3];
sx q[3];
rz(-0.96118569) q[3];
sx q[3];
rz(-3.1072646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2845594) q[0];
sx q[0];
rz(-0.71145809) q[0];
sx q[0];
rz(-0.51666623) q[0];
rz(-0.38756469) q[1];
sx q[1];
rz(-2.0811847) q[1];
sx q[1];
rz(0.26836747) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70290138) q[0];
sx q[0];
rz(-0.99150204) q[0];
sx q[0];
rz(2.2772574) q[0];
rz(-2.1568314) q[2];
sx q[2];
rz(-1.1981989) q[2];
sx q[2];
rz(-0.44024703) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.2942218) q[1];
sx q[1];
rz(-0.98681824) q[1];
sx q[1];
rz(1.8933312) q[1];
rz(-pi) q[2];
rz(1.9645421) q[3];
sx q[3];
rz(-2.5265794) q[3];
sx q[3];
rz(-2.7494489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9364075) q[2];
sx q[2];
rz(-0.78339094) q[2];
sx q[2];
rz(-0.56383413) q[2];
rz(1.1901723) q[3];
sx q[3];
rz(-2.1824013) q[3];
sx q[3];
rz(-0.64175516) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.664809) q[0];
sx q[0];
rz(-1.8905147) q[0];
sx q[0];
rz(2.074194) q[0];
rz(1.339636) q[1];
sx q[1];
rz(-1.4383153) q[1];
sx q[1];
rz(-1.7972606) q[1];
rz(2.6396991) q[2];
sx q[2];
rz(-0.5770275) q[2];
sx q[2];
rz(1.8736476) q[2];
rz(2.8292538) q[3];
sx q[3];
rz(-1.3795508) q[3];
sx q[3];
rz(3.0637904) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
