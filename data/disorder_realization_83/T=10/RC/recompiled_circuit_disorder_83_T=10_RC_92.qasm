OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5873592) q[0];
sx q[0];
rz(-0.98266196) q[0];
sx q[0];
rz(2.3798556) q[0];
rz(-2.3770483) q[1];
sx q[1];
rz(-1.0772871) q[1];
sx q[1];
rz(2.397937) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7997768) q[0];
sx q[0];
rz(-1.8108978) q[0];
sx q[0];
rz(3.0738897) q[0];
rz(-pi) q[1];
rz(-2.9843569) q[2];
sx q[2];
rz(-0.6586282) q[2];
sx q[2];
rz(3.0008891) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.570959) q[1];
sx q[1];
rz(-1.1753517) q[1];
sx q[1];
rz(0.30084893) q[1];
x q[2];
rz(-0.071804382) q[3];
sx q[3];
rz(-1.2076488) q[3];
sx q[3];
rz(-2.8920409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7068229) q[2];
sx q[2];
rz(-0.40199026) q[2];
sx q[2];
rz(-0.10786954) q[2];
rz(-2.9989631) q[3];
sx q[3];
rz(-1.7248036) q[3];
sx q[3];
rz(0.67255783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07664872) q[0];
sx q[0];
rz(-0.80524421) q[0];
sx q[0];
rz(-0.23072492) q[0];
rz(1.327286) q[1];
sx q[1];
rz(-2.4704411) q[1];
sx q[1];
rz(-3.1006295) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3687392) q[0];
sx q[0];
rz(-2.8452747) q[0];
sx q[0];
rz(-0.98451891) q[0];
rz(-1.2626921) q[2];
sx q[2];
rz(-1.8381881) q[2];
sx q[2];
rz(1.2958796) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.47536182) q[1];
sx q[1];
rz(-2.1365709) q[1];
sx q[1];
rz(1.0223785) q[1];
rz(-pi) q[2];
rz(2.326194) q[3];
sx q[3];
rz(-0.5268464) q[3];
sx q[3];
rz(-0.53838733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.162398) q[2];
sx q[2];
rz(-1.6684063) q[2];
sx q[2];
rz(-0.56817788) q[2];
rz(2.6290821) q[3];
sx q[3];
rz(-0.5957225) q[3];
sx q[3];
rz(-2.5797243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42999643) q[0];
sx q[0];
rz(-1.6226409) q[0];
sx q[0];
rz(-2.6913753) q[0];
rz(1.2954767) q[1];
sx q[1];
rz(-1.1643012) q[1];
sx q[1];
rz(2.4643262) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60904658) q[0];
sx q[0];
rz(-2.9653774) q[0];
sx q[0];
rz(0.30058582) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9708423) q[2];
sx q[2];
rz(-0.0064592529) q[2];
sx q[2];
rz(-3.1250172) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8202159) q[1];
sx q[1];
rz(-1.4396588) q[1];
sx q[1];
rz(-1.7507491) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.47795313) q[3];
sx q[3];
rz(-1.2198997) q[3];
sx q[3];
rz(-1.705845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9230817) q[2];
sx q[2];
rz(-1.9873025) q[2];
sx q[2];
rz(0.046860524) q[2];
rz(2.3299407) q[3];
sx q[3];
rz(-2.4620158) q[3];
sx q[3];
rz(0.17015447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99455225) q[0];
sx q[0];
rz(-3.0299598) q[0];
sx q[0];
rz(3.0901093) q[0];
rz(-0.4908081) q[1];
sx q[1];
rz(-2.1665116) q[1];
sx q[1];
rz(-1.1725918) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5347036) q[0];
sx q[0];
rz(-0.67236116) q[0];
sx q[0];
rz(0.21557233) q[0];
rz(-3.1191191) q[2];
sx q[2];
rz(-2.5354676) q[2];
sx q[2];
rz(-1.4215353) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9753032) q[1];
sx q[1];
rz(-1.9235833) q[1];
sx q[1];
rz(-0.3198448) q[1];
x q[2];
rz(2.0026607) q[3];
sx q[3];
rz(-1.5032839) q[3];
sx q[3];
rz(2.1907012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2525758) q[2];
sx q[2];
rz(-0.93095195) q[2];
sx q[2];
rz(0.48689294) q[2];
rz(-1.1934818) q[3];
sx q[3];
rz(-0.32961696) q[3];
sx q[3];
rz(3.1161599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30381969) q[0];
sx q[0];
rz(-2.1210414) q[0];
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
rz(-1.7703169) q[0];
sx q[0];
rz(-1.9195942) q[0];
sx q[0];
rz(2.6692997) q[0];
rz(-pi) q[1];
rz(-2.7388217) q[2];
sx q[2];
rz(-2.6908814) q[2];
sx q[2];
rz(0.41926256) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4795585) q[1];
sx q[1];
rz(-0.91842945) q[1];
sx q[1];
rz(1.0865092) q[1];
rz(-pi) q[2];
rz(-1.6642903) q[3];
sx q[3];
rz(-1.9545385) q[3];
sx q[3];
rz(-1.547471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6902265) q[2];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94474435) q[0];
sx q[0];
rz(-1.9510883) q[0];
sx q[0];
rz(0.29904547) q[0];
rz(-1.8213182) q[1];
sx q[1];
rz(-2.8838005) q[1];
sx q[1];
rz(-1.4978131) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5446536) q[0];
sx q[0];
rz(-2.4785846) q[0];
sx q[0];
rz(-1.6556428) q[0];
rz(-pi) q[1];
rz(2.1645855) q[2];
sx q[2];
rz(-2.0679681) q[2];
sx q[2];
rz(-1.1360816) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.3469997) q[1];
sx q[1];
rz(-1.1113249) q[1];
sx q[1];
rz(-0.32230349) q[1];
x q[2];
rz(-1.3214345) q[3];
sx q[3];
rz(-1.2994088) q[3];
sx q[3];
rz(0.67941487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.51612878) q[2];
sx q[2];
rz(-1.4279782) q[2];
sx q[2];
rz(-3.1141172) q[2];
rz(0.52250683) q[3];
sx q[3];
rz(-0.80111879) q[3];
sx q[3];
rz(-0.81645042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7291173) q[0];
sx q[0];
rz(-1.1183879) q[0];
sx q[0];
rz(-3.0274042) q[0];
rz(-0.97822899) q[1];
sx q[1];
rz(-0.54662919) q[1];
sx q[1];
rz(0.79089975) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1066125) q[0];
sx q[0];
rz(-0.62823409) q[0];
sx q[0];
rz(-1.4094704) q[0];
rz(-pi) q[1];
rz(2.7795243) q[2];
sx q[2];
rz(-0.77713359) q[2];
sx q[2];
rz(2.2287378) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6939047) q[1];
sx q[1];
rz(-1.4920007) q[1];
sx q[1];
rz(-1.2055956) q[1];
x q[2];
rz(-2.0049719) q[3];
sx q[3];
rz(-0.73528157) q[3];
sx q[3];
rz(2.6649464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2699282) q[2];
sx q[2];
rz(-1.1130788) q[2];
sx q[2];
rz(2.2765735) q[2];
rz(0.67251742) q[3];
sx q[3];
rz(-2.0130242) q[3];
sx q[3];
rz(0.095656693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4928116) q[0];
sx q[0];
rz(-2.6219941) q[0];
sx q[0];
rz(0.46736026) q[0];
rz(-0.53721792) q[1];
sx q[1];
rz(-2.1610114) q[1];
sx q[1];
rz(-2.8875202) q[1];
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
rz(-3.1363856) q[2];
sx q[2];
rz(-1.5723096) q[2];
sx q[2];
rz(3.0031799) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.2368187) q[1];
sx q[1];
rz(-2.8042951) q[1];
sx q[1];
rz(-2.7493613) q[1];
rz(-pi) q[2];
rz(1.4905606) q[3];
sx q[3];
rz(-2.1685765) q[3];
sx q[3];
rz(-0.3790006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6415928) q[2];
sx q[2];
rz(-2.3708512) q[2];
sx q[2];
rz(2.8059778) q[2];
rz(-0.32661682) q[3];
sx q[3];
rz(-0.89544046) q[3];
sx q[3];
rz(1.4183104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0582054) q[0];
sx q[0];
rz(-2.931262) q[0];
sx q[0];
rz(1.0539508) q[0];
rz(0.15696934) q[1];
sx q[1];
rz(-1.418768) q[1];
sx q[1];
rz(-0.98186791) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40030038) q[0];
sx q[0];
rz(-1.9809082) q[0];
sx q[0];
rz(1.6923231) q[0];
rz(-1.7392731) q[2];
sx q[2];
rz(-0.57541621) q[2];
sx q[2];
rz(-1.1489431) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5280142) q[1];
sx q[1];
rz(-2.947223) q[1];
sx q[1];
rz(-1.8730875) q[1];
rz(2.8544159) q[3];
sx q[3];
rz(-2.4099726) q[3];
sx q[3];
rz(-1.2958131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2173569) q[2];
sx q[2];
rz(-1.057426) q[2];
sx q[2];
rz(0.99739933) q[2];
rz(-2.635397) q[3];
sx q[3];
rz(-2.180407) q[3];
sx q[3];
rz(-3.1072646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
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
rz(2.754028) q[1];
sx q[1];
rz(-2.0811847) q[1];
sx q[1];
rz(0.26836747) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70290138) q[0];
sx q[0];
rz(-2.1500906) q[0];
sx q[0];
rz(2.2772574) q[0];
rz(0.95558138) q[2];
sx q[2];
rz(-0.68253839) q[2];
sx q[2];
rz(-2.5126484) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.2942218) q[1];
sx q[1];
rz(-0.98681824) q[1];
sx q[1];
rz(-1.2482615) q[1];
x q[2];
rz(1.1770505) q[3];
sx q[3];
rz(-0.6150133) q[3];
sx q[3];
rz(-2.7494489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9364075) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.51856002) q[2];
sx q[2];
rz(-1.8363562) q[2];
sx q[2];
rz(0.73391757) q[2];
rz(-2.8292538) q[3];
sx q[3];
rz(-1.7620419) q[3];
sx q[3];
rz(-0.077802303) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
