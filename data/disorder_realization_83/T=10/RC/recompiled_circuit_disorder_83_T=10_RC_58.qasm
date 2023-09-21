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
rz(0.76454437) q[1];
sx q[1];
rz(4.2188797) q[1];
sx q[1];
rz(10.168434) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0640472) q[0];
sx q[0];
rz(-2.8923058) q[0];
sx q[0];
rz(-1.3011978) q[0];
x q[1];
rz(-2.488963) q[2];
sx q[2];
rz(-1.4748117) q[2];
sx q[2];
rz(-1.8362311) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.570959) q[1];
sx q[1];
rz(-1.966241) q[1];
sx q[1];
rz(2.8407437) q[1];
rz(-pi) q[2];
rz(-3.0697883) q[3];
sx q[3];
rz(-1.9339438) q[3];
sx q[3];
rz(0.24955173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7068229) q[2];
sx q[2];
rz(-0.40199026) q[2];
sx q[2];
rz(-3.0337231) q[2];
rz(0.14262959) q[3];
sx q[3];
rz(-1.7248036) q[3];
sx q[3];
rz(-2.4690348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07664872) q[0];
sx q[0];
rz(-0.80524421) q[0];
sx q[0];
rz(-0.23072492) q[0];
rz(-1.8143066) q[1];
sx q[1];
rz(-0.67115152) q[1];
sx q[1];
rz(-0.040963106) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.165867) q[0];
sx q[0];
rz(-1.3250933) q[0];
sx q[0];
rz(2.9742572) q[0];
rz(-pi) q[1];
rz(-1.2626921) q[2];
sx q[2];
rz(-1.8381881) q[2];
sx q[2];
rz(-1.845713) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6662308) q[1];
sx q[1];
rz(-2.1365709) q[1];
sx q[1];
rz(2.1192141) q[1];
x q[2];
rz(1.1702234) q[3];
sx q[3];
rz(-1.922732) q[3];
sx q[3];
rz(0.34917253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9791947) q[2];
sx q[2];
rz(-1.6684063) q[2];
sx q[2];
rz(0.56817788) q[2];
rz(-0.5125106) q[3];
sx q[3];
rz(-2.5458702) q[3];
sx q[3];
rz(2.5797243) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7115962) q[0];
sx q[0];
rz(-1.6226409) q[0];
sx q[0];
rz(0.45021737) q[0];
rz(1.2954767) q[1];
sx q[1];
rz(-1.9772915) q[1];
sx q[1];
rz(-2.4643262) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8836425) q[0];
sx q[0];
rz(-1.5188688) q[0];
sx q[0];
rz(0.16846637) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1352273) q[2];
sx q[2];
rz(-1.5718939) q[2];
sx q[2];
rz(-1.7249677) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3213768) q[1];
sx q[1];
rz(-1.7019338) q[1];
sx q[1];
rz(1.3908435) q[1];
rz(-pi) q[2];
rz(-1.1797818) q[3];
sx q[3];
rz(-2.0174332) q[3];
sx q[3];
rz(3.1004268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.21851097) q[2];
sx q[2];
rz(-1.9873025) q[2];
sx q[2];
rz(-3.0947321) q[2];
rz(-2.3299407) q[3];
sx q[3];
rz(-0.67957687) q[3];
sx q[3];
rz(0.17015447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99455225) q[0];
sx q[0];
rz(-3.0299598) q[0];
sx q[0];
rz(0.051483367) q[0];
rz(-0.4908081) q[1];
sx q[1];
rz(-0.97508109) q[1];
sx q[1];
rz(-1.9690008) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8797982) q[0];
sx q[0];
rz(-0.9167295) q[0];
sx q[0];
rz(-1.7394702) q[0];
rz(-3.1191191) q[2];
sx q[2];
rz(-2.5354676) q[2];
sx q[2];
rz(-1.4215353) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9753032) q[1];
sx q[1];
rz(-1.9235833) q[1];
sx q[1];
rz(0.3198448) q[1];
rz(-pi) q[2];
rz(-1.4106393) q[3];
sx q[3];
rz(-2.7048116) q[3];
sx q[3];
rz(2.3763451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.88901687) q[2];
sx q[2];
rz(-2.2106407) q[2];
sx q[2];
rz(-0.48689294) q[2];
rz(1.1934818) q[3];
sx q[3];
rz(-2.8119757) q[3];
sx q[3];
rz(3.1161599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.837773) q[0];
sx q[0];
rz(-1.0205512) q[0];
sx q[0];
rz(-1.8160965) q[0];
rz(-0.59108132) q[1];
sx q[1];
rz(-0.68060827) q[1];
sx q[1];
rz(0.65471929) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3723719) q[0];
sx q[0];
rz(-2.0125611) q[0];
sx q[0];
rz(1.9584993) q[0];
rz(2.7388217) q[2];
sx q[2];
rz(-0.45071128) q[2];
sx q[2];
rz(0.41926256) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4795585) q[1];
sx q[1];
rz(-0.91842945) q[1];
sx q[1];
rz(1.0865092) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6642903) q[3];
sx q[3];
rz(-1.1870541) q[3];
sx q[3];
rz(-1.547471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4513662) q[2];
sx q[2];
rz(-1.9407242) q[2];
sx q[2];
rz(2.5704685) q[2];
rz(-2.5518104) q[3];
sx q[3];
rz(-2.6815806) q[3];
sx q[3];
rz(3.0736198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(0.94474435) q[0];
sx q[0];
rz(-1.1905043) q[0];
sx q[0];
rz(2.8425472) q[0];
rz(1.3202745) q[1];
sx q[1];
rz(-2.8838005) q[1];
sx q[1];
rz(-1.4978131) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70443557) q[0];
sx q[0];
rz(-2.2309982) q[0];
sx q[0];
rz(-0.066083834) q[0];
rz(-pi) q[1];
rz(-0.97700714) q[2];
sx q[2];
rz(-1.0736246) q[2];
sx q[2];
rz(1.1360816) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.99243473) q[1];
sx q[1];
rz(-2.5870393) q[1];
sx q[1];
rz(-1.0013594) q[1];
x q[2];
rz(-1.3214345) q[3];
sx q[3];
rz(-1.2994088) q[3];
sx q[3];
rz(0.67941487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6254639) q[2];
sx q[2];
rz(-1.7136145) q[2];
sx q[2];
rz(0.027475474) q[2];
rz(0.52250683) q[3];
sx q[3];
rz(-2.3404739) q[3];
sx q[3];
rz(0.81645042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1066125) q[0];
sx q[0];
rz(-2.5133586) q[0];
sx q[0];
rz(1.4094704) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2355455) q[2];
sx q[2];
rz(-0.8555879) q[2];
sx q[2];
rz(2.7170979) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.447688) q[1];
sx q[1];
rz(-1.4920007) q[1];
sx q[1];
rz(-1.2055956) q[1];
rz(-pi) q[2];
rz(2.2579455) q[3];
sx q[3];
rz(-1.2847319) q[3];
sx q[3];
rz(-1.4253695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.87166446) q[2];
sx q[2];
rz(-1.1130788) q[2];
sx q[2];
rz(-0.86501914) q[2];
rz(-0.67251742) q[3];
sx q[3];
rz(-2.0130242) q[3];
sx q[3];
rz(-0.095656693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4928116) q[0];
sx q[0];
rz(-2.6219941) q[0];
sx q[0];
rz(2.6742324) q[0];
rz(2.6043747) q[1];
sx q[1];
rz(-0.98058128) q[1];
sx q[1];
rz(-0.25407243) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8763037) q[0];
sx q[0];
rz(-2.9348001) q[0];
sx q[0];
rz(0.32949038) q[0];
rz(-pi) q[1];
rz(2.8587564) q[2];
sx q[2];
rz(-3.1361702) q[2];
sx q[2];
rz(1.9920414) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.961812) q[1];
sx q[1];
rz(-1.4439549) q[1];
sx q[1];
rz(-0.3133874) q[1];
x q[2];
rz(1.4905606) q[3];
sx q[3];
rz(-2.1685765) q[3];
sx q[3];
rz(-0.3790006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6415928) q[2];
sx q[2];
rz(-0.77074146) q[2];
sx q[2];
rz(-2.8059778) q[2];
rz(2.8149758) q[3];
sx q[3];
rz(-0.89544046) q[3];
sx q[3];
rz(-1.7232822) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.083387233) q[0];
sx q[0];
rz(-0.21033062) q[0];
sx q[0];
rz(1.0539508) q[0];
rz(-0.15696934) q[1];
sx q[1];
rz(-1.418768) q[1];
sx q[1];
rz(0.98186791) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4440585) q[0];
sx q[0];
rz(-2.7148348) q[0];
sx q[0];
rz(0.27192893) q[0];
x q[1];
rz(-1.0018714) q[2];
sx q[2];
rz(-1.47942) q[2];
sx q[2];
rz(2.5779974) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.30584221) q[1];
sx q[1];
rz(-1.3853449) q[1];
sx q[1];
rz(-3.0830543) q[1];
x q[2];
rz(0.28717678) q[3];
sx q[3];
rz(-0.73162006) q[3];
sx q[3];
rz(-1.2958131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.92423576) q[2];
sx q[2];
rz(-2.0841667) q[2];
sx q[2];
rz(2.1441933) q[2];
rz(0.50619566) q[3];
sx q[3];
rz(-2.180407) q[3];
sx q[3];
rz(-3.1072646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.2845594) q[0];
sx q[0];
rz(-0.71145809) q[0];
sx q[0];
rz(-0.51666623) q[0];
rz(-0.38756469) q[1];
sx q[1];
rz(-1.060408) q[1];
sx q[1];
rz(-0.26836747) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29753387) q[0];
sx q[0];
rz(-0.8807655) q[0];
sx q[0];
rz(0.78155078) q[0];
x q[1];
rz(0.98476121) q[2];
sx q[2];
rz(-1.1981989) q[2];
sx q[2];
rz(2.7013456) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.2942218) q[1];
sx q[1];
rz(-0.98681824) q[1];
sx q[1];
rz(-1.2482615) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8769365) q[3];
sx q[3];
rz(-1.0088682) q[3];
sx q[3];
rz(0.86268007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2051852) q[2];
sx q[2];
rz(-2.3582017) q[2];
sx q[2];
rz(-0.56383413) q[2];
rz(-1.9514203) q[3];
sx q[3];
rz(-0.95919132) q[3];
sx q[3];
rz(0.64175516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47678369) q[0];
sx q[0];
rz(-1.2510779) q[0];
sx q[0];
rz(-1.0673987) q[0];
rz(-1.339636) q[1];
sx q[1];
rz(-1.7032774) q[1];
sx q[1];
rz(1.3443321) q[1];
rz(-0.50189353) q[2];
sx q[2];
rz(-0.5770275) q[2];
sx q[2];
rz(1.8736476) q[2];
rz(1.3700804) q[3];
sx q[3];
rz(-1.8772535) q[3];
sx q[3];
rz(1.5542961) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];