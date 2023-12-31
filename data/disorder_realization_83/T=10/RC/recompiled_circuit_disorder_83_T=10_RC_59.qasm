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
rz(0.24510358) q[0];
sx q[0];
rz(-1.6365543) q[0];
sx q[0];
rz(-1.8114281) q[0];
rz(2.488963) q[2];
sx q[2];
rz(-1.4748117) q[2];
sx q[2];
rz(-1.3053615) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.88120645) q[1];
sx q[1];
rz(-1.8477866) q[1];
sx q[1];
rz(1.1587515) q[1];
x q[2];
rz(3.0697883) q[3];
sx q[3];
rz(-1.2076488) q[3];
sx q[3];
rz(0.24955173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7068229) q[2];
sx q[2];
rz(-0.40199026) q[2];
sx q[2];
rz(0.10786954) q[2];
rz(0.14262959) q[3];
sx q[3];
rz(-1.4167891) q[3];
sx q[3];
rz(2.4690348) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07664872) q[0];
sx q[0];
rz(-0.80524421) q[0];
sx q[0];
rz(2.9108677) q[0];
rz(-1.327286) q[1];
sx q[1];
rz(-2.4704411) q[1];
sx q[1];
rz(-0.040963106) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7777268) q[0];
sx q[0];
rz(-1.7330609) q[0];
sx q[0];
rz(-1.8198387) q[0];
rz(-pi) q[1];
x q[1];
rz(0.83611739) q[2];
sx q[2];
rz(-0.40514075) q[2];
sx q[2];
rz(-0.9678313) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7296655) q[1];
sx q[1];
rz(-1.1150868) q[1];
sx q[1];
rz(-0.63974849) q[1];
rz(-1.9713692) q[3];
sx q[3];
rz(-1.2188606) q[3];
sx q[3];
rz(-0.34917253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.162398) q[2];
sx q[2];
rz(-1.4731864) q[2];
sx q[2];
rz(-0.56817788) q[2];
rz(-2.6290821) q[3];
sx q[3];
rz(-2.5458702) q[3];
sx q[3];
rz(-2.5797243) q[3];
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
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8836425) q[0];
sx q[0];
rz(-1.6227239) q[0];
sx q[0];
rz(-0.16846637) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9708423) q[2];
sx q[2];
rz(-0.0064592529) q[2];
sx q[2];
rz(3.1250172) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2689506) q[1];
sx q[1];
rz(-2.9193455) q[1];
sx q[1];
rz(-2.2058669) q[1];
rz(-pi) q[2];
rz(-1.9618109) q[3];
sx q[3];
rz(-1.1241594) q[3];
sx q[3];
rz(-0.041165813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.21851097) q[2];
sx q[2];
rz(-1.1542902) q[2];
sx q[2];
rz(0.046860524) q[2];
rz(0.81165195) q[3];
sx q[3];
rz(-2.4620158) q[3];
sx q[3];
rz(-0.17015447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99455225) q[0];
sx q[0];
rz(-0.11163286) q[0];
sx q[0];
rz(3.0901093) q[0];
rz(-0.4908081) q[1];
sx q[1];
rz(-2.1665116) q[1];
sx q[1];
rz(-1.1725918) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6068891) q[0];
sx q[0];
rz(-0.67236116) q[0];
sx q[0];
rz(-2.9260203) q[0];
x q[1];
rz(3.1191191) q[2];
sx q[2];
rz(-2.5354676) q[2];
sx q[2];
rz(1.4215353) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.51845156) q[1];
sx q[1];
rz(-1.8703096) q[1];
sx q[1];
rz(-1.2007984) q[1];
rz(0.074313642) q[3];
sx q[3];
rz(-1.1399817) q[3];
sx q[3];
rz(-0.58882344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2525758) q[2];
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
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30381969) q[0];
sx q[0];
rz(-1.0205512) q[0];
sx q[0];
rz(1.8160965) q[0];
rz(0.59108132) q[1];
sx q[1];
rz(-0.68060827) q[1];
sx q[1];
rz(-0.65471929) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3712758) q[0];
sx q[0];
rz(-1.2219984) q[0];
sx q[0];
rz(-2.6692997) q[0];
rz(-1.7582558) q[2];
sx q[2];
rz(-1.1584632) q[2];
sx q[2];
rz(-0.86134855) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.66203413) q[1];
sx q[1];
rz(-2.2231632) q[1];
sx q[1];
rz(2.0550834) q[1];
rz(-pi) q[2];
rz(0.22722865) q[3];
sx q[3];
rz(-0.39441808) q[3];
sx q[3];
rz(1.7928746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4513662) q[2];
sx q[2];
rz(-1.2008685) q[2];
sx q[2];
rz(-0.5711242) q[2];
rz(-0.58978224) q[3];
sx q[3];
rz(-2.6815806) q[3];
sx q[3];
rz(-3.0736198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1968483) q[0];
sx q[0];
rz(-1.9510883) q[0];
sx q[0];
rz(0.29904547) q[0];
rz(1.3202745) q[1];
sx q[1];
rz(-0.25779217) q[1];
sx q[1];
rz(-1.6437795) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90692524) q[0];
sx q[0];
rz(-1.5186131) q[0];
sx q[0];
rz(-2.232057) q[0];
rz(-pi) q[1];
x q[1];
rz(0.80070337) q[2];
sx q[2];
rz(-0.75468894) q[2];
sx q[2];
rz(-0.18037361) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.99243473) q[1];
sx q[1];
rz(-0.5545534) q[1];
sx q[1];
rz(2.1402332) q[1];
x q[2];
rz(-2.8619814) q[3];
sx q[3];
rz(-1.3307443) q[3];
sx q[3];
rz(0.9595426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6254639) q[2];
sx q[2];
rz(-1.7136145) q[2];
sx q[2];
rz(-3.1141172) q[2];
rz(2.6190858) q[3];
sx q[3];
rz(-0.80111879) q[3];
sx q[3];
rz(-2.3251422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41247535) q[0];
sx q[0];
rz(-1.1183879) q[0];
sx q[0];
rz(3.0274042) q[0];
rz(-0.97822899) q[1];
sx q[1];
rz(-2.5949635) q[1];
sx q[1];
rz(-0.79089975) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1066125) q[0];
sx q[0];
rz(-2.5133586) q[0];
sx q[0];
rz(1.4094704) q[0];
rz(-1.2355455) q[2];
sx q[2];
rz(-0.8555879) q[2];
sx q[2];
rz(2.7170979) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6939047) q[1];
sx q[1];
rz(-1.649592) q[1];
sx q[1];
rz(-1.2055956) q[1];
rz(-0.36356504) q[3];
sx q[3];
rz(-2.225038) q[3];
sx q[3];
rz(-3.0594861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2699282) q[2];
sx q[2];
rz(-1.1130788) q[2];
sx q[2];
rz(-2.2765735) q[2];
rz(-2.4690752) q[3];
sx q[3];
rz(-2.0130242) q[3];
sx q[3];
rz(-3.045936) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6487811) q[0];
sx q[0];
rz(-2.6219941) q[0];
sx q[0];
rz(2.6742324) q[0];
rz(2.6043747) q[1];
sx q[1];
rz(-0.98058128) q[1];
sx q[1];
rz(2.8875202) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8763037) q[0];
sx q[0];
rz(-2.9348001) q[0];
sx q[0];
rz(-2.8121023) q[0];
x q[1];
rz(2.8587564) q[2];
sx q[2];
rz(-0.0054224646) q[2];
sx q[2];
rz(-1.9920414) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1797807) q[1];
sx q[1];
rz(-1.6976377) q[1];
sx q[1];
rz(0.3133874) q[1];
rz(-pi) q[2];
rz(-2.5423126) q[3];
sx q[3];
rz(-1.6370956) q[3];
sx q[3];
rz(-1.9950206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.49999985) q[2];
sx q[2];
rz(-0.77074146) q[2];
sx q[2];
rz(-2.8059778) q[2];
rz(2.8149758) q[3];
sx q[3];
rz(-2.2461522) q[3];
sx q[3];
rz(1.7232822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0582054) q[0];
sx q[0];
rz(-2.931262) q[0];
sx q[0];
rz(-1.0539508) q[0];
rz(-2.9846233) q[1];
sx q[1];
rz(-1.7228246) q[1];
sx q[1];
rz(0.98186791) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40030038) q[0];
sx q[0];
rz(-1.1606845) q[0];
sx q[0];
rz(-1.6923231) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0018714) q[2];
sx q[2];
rz(-1.47942) q[2];
sx q[2];
rz(2.5779974) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8874444) q[1];
sx q[1];
rz(-1.5132628) q[1];
sx q[1];
rz(-1.7565586) q[1];
rz(-pi) q[2];
rz(0.710886) q[3];
sx q[3];
rz(-1.7611739) q[3];
sx q[3];
rz(-0.058660942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.92423576) q[2];
sx q[2];
rz(-1.057426) q[2];
sx q[2];
rz(0.99739933) q[2];
rz(2.635397) q[3];
sx q[3];
rz(-2.180407) q[3];
sx q[3];
rz(3.1072646) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85703325) q[0];
sx q[0];
rz(-2.4301346) q[0];
sx q[0];
rz(0.51666623) q[0];
rz(2.754028) q[1];
sx q[1];
rz(-2.0811847) q[1];
sx q[1];
rz(0.26836747) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8440588) q[0];
sx q[0];
rz(-2.2608272) q[0];
sx q[0];
rz(-0.78155078) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1860113) q[2];
sx q[2];
rz(-0.68253839) q[2];
sx q[2];
rz(2.5126484) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.2942218) q[1];
sx q[1];
rz(-2.1547744) q[1];
sx q[1];
rz(-1.8933312) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1770505) q[3];
sx q[3];
rz(-2.5265794) q[3];
sx q[3];
rz(0.39214373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9364075) q[2];
sx q[2];
rz(-0.78339094) q[2];
sx q[2];
rz(2.5777585) q[2];
rz(-1.9514203) q[3];
sx q[3];
rz(-2.1824013) q[3];
sx q[3];
rz(-0.64175516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.664809) q[0];
sx q[0];
rz(-1.2510779) q[0];
sx q[0];
rz(-1.0673987) q[0];
rz(-1.8019567) q[1];
sx q[1];
rz(-1.4383153) q[1];
sx q[1];
rz(-1.7972606) q[1];
rz(-2.6230326) q[2];
sx q[2];
rz(-1.8363562) q[2];
sx q[2];
rz(0.73391757) q[2];
rz(1.7715122) q[3];
sx q[3];
rz(-1.2643391) q[3];
sx q[3];
rz(-1.5872965) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
