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
rz(4.1242546) q[0];
sx q[0];
rz(10.186515) q[0];
rz(-2.3770483) q[1];
sx q[1];
rz(-1.0772871) q[1];
sx q[1];
rz(-0.74365562) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7997768) q[0];
sx q[0];
rz(-1.3306949) q[0];
sx q[0];
rz(-0.067702985) q[0];
x q[1];
rz(-0.6526297) q[2];
sx q[2];
rz(-1.666781) q[2];
sx q[2];
rz(-1.8362311) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.248677) q[1];
sx q[1];
rz(-2.6495653) q[1];
sx q[1];
rz(0.9534652) q[1];
x q[2];
rz(0.071804382) q[3];
sx q[3];
rz(-1.9339438) q[3];
sx q[3];
rz(0.24955173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4347697) q[2];
sx q[2];
rz(-2.7396024) q[2];
sx q[2];
rz(0.10786954) q[2];
rz(0.14262959) q[3];
sx q[3];
rz(-1.4167891) q[3];
sx q[3];
rz(-0.67255783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.07664872) q[0];
sx q[0];
rz(-2.3363484) q[0];
sx q[0];
rz(-0.23072492) q[0];
rz(1.327286) q[1];
sx q[1];
rz(-2.4704411) q[1];
sx q[1];
rz(0.040963106) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.165867) q[0];
sx q[0];
rz(-1.3250933) q[0];
sx q[0];
rz(-2.9742572) q[0];
rz(-1.8789005) q[2];
sx q[2];
rz(-1.3034046) q[2];
sx q[2];
rz(1.2958796) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7663824) q[1];
sx q[1];
rz(-2.3751405) q[1];
sx q[1];
rz(-2.4541928) q[1];
rz(-pi) q[2];
rz(2.7621272) q[3];
sx q[3];
rz(-1.1960408) q[3];
sx q[3];
rz(-1.3665762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9791947) q[2];
sx q[2];
rz(-1.6684063) q[2];
sx q[2];
rz(-0.56817788) q[2];
rz(-0.5125106) q[3];
sx q[3];
rz(-2.5458702) q[3];
sx q[3];
rz(2.5797243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7115962) q[0];
sx q[0];
rz(-1.5189518) q[0];
sx q[0];
rz(2.6913753) q[0];
rz(1.2954767) q[1];
sx q[1];
rz(-1.9772915) q[1];
sx q[1];
rz(-2.4643262) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5325461) q[0];
sx q[0];
rz(-0.17621528) q[0];
sx q[0];
rz(-0.30058582) q[0];
rz(-pi) q[1];
x q[1];
rz(0.1707503) q[2];
sx q[2];
rz(-3.1351334) q[2];
sx q[2];
rz(-0.016575459) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.87264204) q[1];
sx q[1];
rz(-2.9193455) q[1];
sx q[1];
rz(2.2058669) q[1];
rz(-pi) q[2];
rz(0.47795313) q[3];
sx q[3];
rz(-1.9216929) q[3];
sx q[3];
rz(1.4357476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9230817) q[2];
sx q[2];
rz(-1.1542902) q[2];
sx q[2];
rz(3.0947321) q[2];
rz(-0.81165195) q[3];
sx q[3];
rz(-2.4620158) q[3];
sx q[3];
rz(-2.9714382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(2.1470404) q[0];
sx q[0];
rz(-0.11163286) q[0];
sx q[0];
rz(-0.051483367) q[0];
rz(-2.6507846) q[1];
sx q[1];
rz(-2.1665116) q[1];
sx q[1];
rz(1.1725918) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9358312) q[0];
sx q[0];
rz(-1.4371705) q[0];
sx q[0];
rz(2.4806116) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1191191) q[2];
sx q[2];
rz(-0.60612504) q[2];
sx q[2];
rz(1.7200574) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6231411) q[1];
sx q[1];
rz(-1.8703096) q[1];
sx q[1];
rz(1.9407942) q[1];
rz(-pi) q[2];
rz(-3.067279) q[3];
sx q[3];
rz(-2.001611) q[3];
sx q[3];
rz(0.58882344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.88901687) q[2];
sx q[2];
rz(-2.2106407) q[2];
sx q[2];
rz(-2.6546997) q[2];
rz(-1.1934818) q[3];
sx q[3];
rz(-0.32961696) q[3];
sx q[3];
rz(-0.025432767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30381969) q[0];
sx q[0];
rz(-1.0205512) q[0];
sx q[0];
rz(1.8160965) q[0];
rz(-2.5505113) q[1];
sx q[1];
rz(-0.68060827) q[1];
sx q[1];
rz(-0.65471929) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7692208) q[0];
sx q[0];
rz(-2.0125611) q[0];
sx q[0];
rz(1.9584993) q[0];
x q[1];
rz(-1.3833369) q[2];
sx q[2];
rz(-1.9831295) q[2];
sx q[2];
rz(2.2802441) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3760738) q[1];
sx q[1];
rz(-2.3507833) q[1];
sx q[1];
rz(0.54733025) q[1];
x q[2];
rz(0.38526411) q[3];
sx q[3];
rz(-1.6574727) q[3];
sx q[3];
rz(-0.011766089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4513662) q[2];
sx q[2];
rz(-1.9407242) q[2];
sx q[2];
rz(-2.5704685) q[2];
rz(-2.5518104) q[3];
sx q[3];
rz(-0.46001205) q[3];
sx q[3];
rz(0.067972876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1968483) q[0];
sx q[0];
rz(-1.9510883) q[0];
sx q[0];
rz(2.8425472) q[0];
rz(1.3202745) q[1];
sx q[1];
rz(-2.8838005) q[1];
sx q[1];
rz(-1.4978131) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5446536) q[0];
sx q[0];
rz(-0.66300809) q[0];
sx q[0];
rz(1.4859499) q[0];
rz(-pi) q[1];
rz(0.97700714) q[2];
sx q[2];
rz(-2.0679681) q[2];
sx q[2];
rz(-2.0055111) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.076768) q[1];
sx q[1];
rz(-1.2829363) q[1];
sx q[1];
rz(1.0899781) q[1];
rz(-pi) q[2];
x q[2];
rz(0.72553708) q[3];
sx q[3];
rz(-0.36645884) q[3];
sx q[3];
rz(-3.0612502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.51612878) q[2];
sx q[2];
rz(-1.7136145) q[2];
sx q[2];
rz(3.1141172) q[2];
rz(2.6190858) q[3];
sx q[3];
rz(-0.80111879) q[3];
sx q[3];
rz(0.81645042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7291173) q[0];
sx q[0];
rz(-2.0232047) q[0];
sx q[0];
rz(3.0274042) q[0];
rz(2.1633637) q[1];
sx q[1];
rz(-0.54662919) q[1];
sx q[1];
rz(0.79089975) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7366911) q[0];
sx q[0];
rz(-1.4762523) q[0];
sx q[0];
rz(-2.1928284) q[0];
rz(-2.3979264) q[2];
sx q[2];
rz(-1.3197834) q[2];
sx q[2];
rz(2.2199092) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6939047) q[1];
sx q[1];
rz(-1.4920007) q[1];
sx q[1];
rz(-1.9359971) q[1];
x q[2];
rz(1.1366208) q[3];
sx q[3];
rz(-0.73528157) q[3];
sx q[3];
rz(2.6649464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.87166446) q[2];
sx q[2];
rz(-1.1130788) q[2];
sx q[2];
rz(0.86501914) q[2];
rz(-2.4690752) q[3];
sx q[3];
rz(-1.1285684) q[3];
sx q[3];
rz(3.045936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6487811) q[0];
sx q[0];
rz(-2.6219941) q[0];
sx q[0];
rz(-2.6742324) q[0];
rz(0.53721792) q[1];
sx q[1];
rz(-2.1610114) q[1];
sx q[1];
rz(-0.25407243) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5131322) q[0];
sx q[0];
rz(-1.6372794) q[0];
sx q[0];
rz(2.9456338) q[0];
x q[1];
rz(0.28283624) q[2];
sx q[2];
rz(-3.1361702) q[2];
sx q[2];
rz(1.1495513) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.2368187) q[1];
sx q[1];
rz(-0.33729759) q[1];
sx q[1];
rz(-0.39223139) q[1];
rz(1.4905606) q[3];
sx q[3];
rz(-2.1685765) q[3];
sx q[3];
rz(2.762592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.49999985) q[2];
sx q[2];
rz(-0.77074146) q[2];
sx q[2];
rz(0.33561486) q[2];
rz(2.8149758) q[3];
sx q[3];
rz(-0.89544046) q[3];
sx q[3];
rz(1.4183104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083387233) q[0];
sx q[0];
rz(-2.931262) q[0];
sx q[0];
rz(1.0539508) q[0];
rz(-0.15696934) q[1];
sx q[1];
rz(-1.7228246) q[1];
sx q[1];
rz(-0.98186791) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4440585) q[0];
sx q[0];
rz(-2.7148348) q[0];
sx q[0];
rz(2.8696637) q[0];
rz(-0.1083381) q[2];
sx q[2];
rz(-1.0045369) q[2];
sx q[2];
rz(0.94891753) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.30584221) q[1];
sx q[1];
rz(-1.7562477) q[1];
sx q[1];
rz(-0.058538392) q[1];
x q[2];
rz(-1.8198265) q[3];
sx q[3];
rz(-2.2662275) q[3];
sx q[3];
rz(1.673656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.92423576) q[2];
sx q[2];
rz(-2.0841667) q[2];
sx q[2];
rz(0.99739933) q[2];
rz(2.635397) q[3];
sx q[3];
rz(-2.180407) q[3];
sx q[3];
rz(-0.034328073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2845594) q[0];
sx q[0];
rz(-2.4301346) q[0];
sx q[0];
rz(-2.6249264) q[0];
rz(-0.38756469) q[1];
sx q[1];
rz(-2.0811847) q[1];
sx q[1];
rz(0.26836747) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4386913) q[0];
sx q[0];
rz(-2.1500906) q[0];
sx q[0];
rz(2.2772574) q[0];
rz(-pi) q[1];
rz(2.7029413) q[2];
sx q[2];
rz(-1.029656) q[2];
sx q[2];
rz(1.7739319) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8473709) q[1];
sx q[1];
rz(-0.98681824) q[1];
sx q[1];
rz(-1.8933312) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1488232) q[3];
sx q[3];
rz(-1.3475932) q[3];
sx q[3];
rz(-2.2900667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2051852) q[2];
sx q[2];
rz(-0.78339094) q[2];
sx q[2];
rz(-2.5777585) q[2];
rz(1.9514203) q[3];
sx q[3];
rz(-0.95919132) q[3];
sx q[3];
rz(-0.64175516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.664809) q[0];
sx q[0];
rz(-1.8905147) q[0];
sx q[0];
rz(2.074194) q[0];
rz(-1.8019567) q[1];
sx q[1];
rz(-1.4383153) q[1];
sx q[1];
rz(-1.7972606) q[1];
rz(2.6230326) q[2];
sx q[2];
rz(-1.3052365) q[2];
sx q[2];
rz(-2.4076751) q[2];
rz(0.56223829) q[3];
sx q[3];
rz(-2.7769965) q[3];
sx q[3];
rz(-2.1806352) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
