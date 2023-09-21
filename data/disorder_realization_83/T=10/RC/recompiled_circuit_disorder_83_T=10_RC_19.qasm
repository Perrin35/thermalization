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
rz(0.76454437) q[1];
sx q[1];
rz(-2.0643056) q[1];
sx q[1];
rz(0.74365562) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3418158) q[0];
sx q[0];
rz(-1.8108978) q[0];
sx q[0];
rz(-0.067702985) q[0];
x q[1];
rz(-1.450199) q[2];
sx q[2];
rz(-2.2199124) q[2];
sx q[2];
rz(-0.33855864) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.570959) q[1];
sx q[1];
rz(-1.966241) q[1];
sx q[1];
rz(-0.30084893) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2067912) q[3];
sx q[3];
rz(-1.6379106) q[3];
sx q[3];
rz(-1.2957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4347697) q[2];
sx q[2];
rz(-0.40199026) q[2];
sx q[2];
rz(3.0337231) q[2];
rz(0.14262959) q[3];
sx q[3];
rz(-1.4167891) q[3];
sx q[3];
rz(-0.67255783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.07664872) q[0];
sx q[0];
rz(-2.3363484) q[0];
sx q[0];
rz(2.9108677) q[0];
rz(-1.8143066) q[1];
sx q[1];
rz(-2.4704411) q[1];
sx q[1];
rz(0.040963106) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9757257) q[0];
sx q[0];
rz(-1.3250933) q[0];
sx q[0];
rz(0.16733549) q[0];
rz(-pi) q[1];
rz(-2.8616521) q[2];
sx q[2];
rz(-1.8676114) q[2];
sx q[2];
rz(0.1910301) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.47536182) q[1];
sx q[1];
rz(-1.0050217) q[1];
sx q[1];
rz(1.0223785) q[1];
rz(0.81539865) q[3];
sx q[3];
rz(-0.5268464) q[3];
sx q[3];
rz(-2.6032053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.162398) q[2];
sx q[2];
rz(-1.4731864) q[2];
sx q[2];
rz(0.56817788) q[2];
rz(0.5125106) q[3];
sx q[3];
rz(-2.5458702) q[3];
sx q[3];
rz(-2.5797243) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42999643) q[0];
sx q[0];
rz(-1.5189518) q[0];
sx q[0];
rz(-0.45021737) q[0];
rz(-1.2954767) q[1];
sx q[1];
rz(-1.9772915) q[1];
sx q[1];
rz(2.4643262) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2579502) q[0];
sx q[0];
rz(-1.6227239) q[0];
sx q[0];
rz(2.9731263) q[0];
x q[1];
rz(1.5696987) q[2];
sx q[2];
rz(-1.5771616) q[2];
sx q[2];
rz(-2.9874143) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.87264204) q[1];
sx q[1];
rz(-2.9193455) q[1];
sx q[1];
rz(0.93572576) q[1];
rz(-pi) q[2];
rz(2.4694091) q[3];
sx q[3];
rz(-0.58478343) q[3];
sx q[3];
rz(2.4206846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.21851097) q[2];
sx q[2];
rz(-1.1542902) q[2];
sx q[2];
rz(3.0947321) q[2];
rz(-0.81165195) q[3];
sx q[3];
rz(-0.67957687) q[3];
sx q[3];
rz(-0.17015447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99455225) q[0];
sx q[0];
rz(-3.0299598) q[0];
sx q[0];
rz(0.051483367) q[0];
rz(2.6507846) q[1];
sx q[1];
rz(-0.97508109) q[1];
sx q[1];
rz(1.1725918) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2617944) q[0];
sx q[0];
rz(-2.2248631) q[0];
sx q[0];
rz(-1.4021224) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1191191) q[2];
sx q[2];
rz(-0.60612504) q[2];
sx q[2];
rz(-1.4215353) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7395775) q[1];
sx q[1];
rz(-2.6699454) q[1];
sx q[1];
rz(2.2775843) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.067279) q[3];
sx q[3];
rz(-1.1399817) q[3];
sx q[3];
rz(2.5527692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.88901687) q[2];
sx q[2];
rz(-2.2106407) q[2];
sx q[2];
rz(-2.6546997) q[2];
rz(1.1934818) q[3];
sx q[3];
rz(-0.32961696) q[3];
sx q[3];
rz(-3.1161599) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30381969) q[0];
sx q[0];
rz(-2.1210414) q[0];
sx q[0];
rz(1.3254962) q[0];
rz(-0.59108132) q[1];
sx q[1];
rz(-2.4609844) q[1];
sx q[1];
rz(2.4868734) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7692208) q[0];
sx q[0];
rz(-2.0125611) q[0];
sx q[0];
rz(-1.9584993) q[0];
rz(-pi) q[1];
rz(-1.3833369) q[2];
sx q[2];
rz(-1.9831295) q[2];
sx q[2];
rz(2.2802441) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.66203413) q[1];
sx q[1];
rz(-2.2231632) q[1];
sx q[1];
rz(-1.0865092) q[1];
rz(-pi) q[2];
x q[2];
rz(2.914364) q[3];
sx q[3];
rz(-0.39441808) q[3];
sx q[3];
rz(1.3487181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4513662) q[2];
sx q[2];
rz(-1.2008685) q[2];
sx q[2];
rz(-0.5711242) q[2];
rz(-2.5518104) q[3];
sx q[3];
rz(-2.6815806) q[3];
sx q[3];
rz(-0.067972876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1968483) q[0];
sx q[0];
rz(-1.1905043) q[0];
sx q[0];
rz(0.29904547) q[0];
rz(1.8213182) q[1];
sx q[1];
rz(-0.25779217) q[1];
sx q[1];
rz(1.6437795) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70443557) q[0];
sx q[0];
rz(-0.91059443) q[0];
sx q[0];
rz(-0.066083834) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3408893) q[2];
sx q[2];
rz(-0.75468894) q[2];
sx q[2];
rz(0.18037361) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1491579) q[1];
sx q[1];
rz(-2.5870393) q[1];
sx q[1];
rz(1.0013594) q[1];
x q[2];
rz(2.8619814) q[3];
sx q[3];
rz(-1.8108484) q[3];
sx q[3];
rz(0.9595426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6254639) q[2];
sx q[2];
rz(-1.7136145) q[2];
sx q[2];
rz(-0.027475474) q[2];
rz(-0.52250683) q[3];
sx q[3];
rz(-2.3404739) q[3];
sx q[3];
rz(2.3251422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7291173) q[0];
sx q[0];
rz(-2.0232047) q[0];
sx q[0];
rz(0.11418848) q[0];
rz(-2.1633637) q[1];
sx q[1];
rz(-0.54662919) q[1];
sx q[1];
rz(-0.79089975) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7366911) q[0];
sx q[0];
rz(-1.4762523) q[0];
sx q[0];
rz(0.94876429) q[0];
rz(-pi) q[1];
rz(2.7795243) q[2];
sx q[2];
rz(-0.77713359) q[2];
sx q[2];
rz(2.2287378) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.32614) q[1];
sx q[1];
rz(-2.7683612) q[1];
sx q[1];
rz(1.7883854) q[1];
rz(-pi) q[2];
rz(0.36356504) q[3];
sx q[3];
rz(-2.225038) q[3];
sx q[3];
rz(-0.082106575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2699282) q[2];
sx q[2];
rz(-1.1130788) q[2];
sx q[2];
rz(2.2765735) q[2];
rz(2.4690752) q[3];
sx q[3];
rz(-2.0130242) q[3];
sx q[3];
rz(3.045936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6487811) q[0];
sx q[0];
rz(-0.5195986) q[0];
sx q[0];
rz(-2.6742324) q[0];
rz(0.53721792) q[1];
sx q[1];
rz(-0.98058128) q[1];
sx q[1];
rz(0.25407243) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5131322) q[0];
sx q[0];
rz(-1.6372794) q[0];
sx q[0];
rz(0.19595887) q[0];
x q[1];
rz(-2.8587564) q[2];
sx q[2];
rz(-0.0054224646) q[2];
sx q[2];
rz(-1.1495513) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1797807) q[1];
sx q[1];
rz(-1.6976377) q[1];
sx q[1];
rz(-2.8282053) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5423126) q[3];
sx q[3];
rz(-1.5044971) q[3];
sx q[3];
rz(-1.9950206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.49999985) q[2];
sx q[2];
rz(-0.77074146) q[2];
sx q[2];
rz(-2.8059778) q[2];
rz(-0.32661682) q[3];
sx q[3];
rz(-2.2461522) q[3];
sx q[3];
rz(1.7232822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0582054) q[0];
sx q[0];
rz(-2.931262) q[0];
sx q[0];
rz(2.0876419) q[0];
rz(2.9846233) q[1];
sx q[1];
rz(-1.7228246) q[1];
sx q[1];
rz(-0.98186791) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0197524) q[0];
sx q[0];
rz(-1.4593908) q[0];
sx q[0];
rz(-0.41282546) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1397212) q[2];
sx q[2];
rz(-1.47942) q[2];
sx q[2];
rz(2.5779974) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2541483) q[1];
sx q[1];
rz(-1.5132628) q[1];
sx q[1];
rz(-1.385034) q[1];
rz(-pi) q[2];
rz(-2.8544159) q[3];
sx q[3];
rz(-0.73162006) q[3];
sx q[3];
rz(1.8457796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.92423576) q[2];
sx q[2];
rz(-2.0841667) q[2];
sx q[2];
rz(0.99739933) q[2];
rz(-2.635397) q[3];
sx q[3];
rz(-0.96118569) q[3];
sx q[3];
rz(3.1072646) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2845594) q[0];
sx q[0];
rz(-0.71145809) q[0];
sx q[0];
rz(-2.6249264) q[0];
rz(0.38756469) q[1];
sx q[1];
rz(-1.060408) q[1];
sx q[1];
rz(0.26836747) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29753387) q[0];
sx q[0];
rz(-0.8807655) q[0];
sx q[0];
rz(-2.3600419) q[0];
rz(2.1568314) q[2];
sx q[2];
rz(-1.1981989) q[2];
sx q[2];
rz(0.44024703) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8473709) q[1];
sx q[1];
rz(-2.1547744) q[1];
sx q[1];
rz(1.8933312) q[1];
rz(-pi) q[2];
rz(-2.1488232) q[3];
sx q[3];
rz(-1.3475932) q[3];
sx q[3];
rz(2.2900667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2051852) q[2];
sx q[2];
rz(-0.78339094) q[2];
sx q[2];
rz(0.56383413) q[2];
rz(1.1901723) q[3];
sx q[3];
rz(-2.1824013) q[3];
sx q[3];
rz(-0.64175516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47678369) q[0];
sx q[0];
rz(-1.8905147) q[0];
sx q[0];
rz(2.074194) q[0];
rz(1.8019567) q[1];
sx q[1];
rz(-1.7032774) q[1];
sx q[1];
rz(1.3443321) q[1];
rz(-2.6230326) q[2];
sx q[2];
rz(-1.8363562) q[2];
sx q[2];
rz(0.73391757) q[2];
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
