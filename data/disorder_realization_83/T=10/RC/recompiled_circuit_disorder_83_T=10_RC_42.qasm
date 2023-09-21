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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3418158) q[0];
sx q[0];
rz(-1.8108978) q[0];
sx q[0];
rz(-0.067702985) q[0];
rz(-pi) q[1];
rz(-2.488963) q[2];
sx q[2];
rz(-1.4748117) q[2];
sx q[2];
rz(-1.8362311) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.570959) q[1];
sx q[1];
rz(-1.966241) q[1];
sx q[1];
rz(0.30084893) q[1];
x q[2];
rz(1.7573962) q[3];
sx q[3];
rz(-2.7717234) q[3];
sx q[3];
rz(2.6922525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4347697) q[2];
sx q[2];
rz(-2.7396024) q[2];
sx q[2];
rz(3.0337231) q[2];
rz(2.9989631) q[3];
sx q[3];
rz(-1.4167891) q[3];
sx q[3];
rz(0.67255783) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07664872) q[0];
sx q[0];
rz(-0.80524421) q[0];
sx q[0];
rz(2.9108677) q[0];
rz(1.8143066) q[1];
sx q[1];
rz(-2.4704411) q[1];
sx q[1];
rz(-0.040963106) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9757257) q[0];
sx q[0];
rz(-1.3250933) q[0];
sx q[0];
rz(-0.16733549) q[0];
rz(-pi) q[1];
rz(2.3054753) q[2];
sx q[2];
rz(-2.7364519) q[2];
sx q[2];
rz(-0.9678313) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7296655) q[1];
sx q[1];
rz(-2.0265059) q[1];
sx q[1];
rz(2.5018442) q[1];
rz(-pi) q[2];
x q[2];
rz(0.37946545) q[3];
sx q[3];
rz(-1.9455519) q[3];
sx q[3];
rz(1.7750164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.162398) q[2];
sx q[2];
rz(-1.4731864) q[2];
sx q[2];
rz(-0.56817788) q[2];
rz(-0.5125106) q[3];
sx q[3];
rz(-0.5957225) q[3];
sx q[3];
rz(0.56186831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7115962) q[0];
sx q[0];
rz(-1.6226409) q[0];
sx q[0];
rz(-2.6913753) q[0];
rz(1.2954767) q[1];
sx q[1];
rz(-1.1643012) q[1];
sx q[1];
rz(-0.67726642) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2579502) q[0];
sx q[0];
rz(-1.6227239) q[0];
sx q[0];
rz(-2.9731263) q[0];
x q[1];
rz(2.9708423) q[2];
sx q[2];
rz(-0.0064592529) q[2];
sx q[2];
rz(-0.016575459) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.22563572) q[1];
sx q[1];
rz(-1.7491873) q[1];
sx q[1];
rz(0.13326463) q[1];
rz(-2.6636395) q[3];
sx q[3];
rz(-1.9216929) q[3];
sx q[3];
rz(-1.705845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.21851097) q[2];
sx q[2];
rz(-1.9873025) q[2];
sx q[2];
rz(0.046860524) q[2];
rz(2.3299407) q[3];
sx q[3];
rz(-2.4620158) q[3];
sx q[3];
rz(-2.9714382) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1470404) q[0];
sx q[0];
rz(-3.0299598) q[0];
sx q[0];
rz(3.0901093) q[0];
rz(2.6507846) q[1];
sx q[1];
rz(-0.97508109) q[1];
sx q[1];
rz(1.1725918) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20576142) q[0];
sx q[0];
rz(-1.7044221) q[0];
sx q[0];
rz(0.66098102) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1191191) q[2];
sx q[2];
rz(-0.60612504) q[2];
sx q[2];
rz(1.4215353) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7395775) q[1];
sx q[1];
rz(-0.47164729) q[1];
sx q[1];
rz(2.2775843) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7309534) q[3];
sx q[3];
rz(-2.7048116) q[3];
sx q[3];
rz(-2.3763451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2525758) q[2];
sx q[2];
rz(-2.2106407) q[2];
sx q[2];
rz(0.48689294) q[2];
rz(-1.1934818) q[3];
sx q[3];
rz(-2.8119757) q[3];
sx q[3];
rz(0.025432767) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.837773) q[0];
sx q[0];
rz(-2.1210414) q[0];
sx q[0];
rz(1.3254962) q[0];
rz(-0.59108132) q[1];
sx q[1];
rz(-0.68060827) q[1];
sx q[1];
rz(0.65471929) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7692208) q[0];
sx q[0];
rz(-1.1290316) q[0];
sx q[0];
rz(1.1830933) q[0];
x q[1];
rz(2.7388217) q[2];
sx q[2];
rz(-2.6908814) q[2];
sx q[2];
rz(2.7223301) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3760738) q[1];
sx q[1];
rz(-0.79080938) q[1];
sx q[1];
rz(0.54733025) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6642903) q[3];
sx q[3];
rz(-1.1870541) q[3];
sx q[3];
rz(1.547471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4513662) q[2];
sx q[2];
rz(-1.2008685) q[2];
sx q[2];
rz(0.5711242) q[2];
rz(2.5518104) q[3];
sx q[3];
rz(-0.46001205) q[3];
sx q[3];
rz(-0.067972876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(0.94474435) q[0];
sx q[0];
rz(-1.9510883) q[0];
sx q[0];
rz(-2.8425472) q[0];
rz(1.8213182) q[1];
sx q[1];
rz(-0.25779217) q[1];
sx q[1];
rz(1.6437795) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4371571) q[0];
sx q[0];
rz(-2.2309982) q[0];
sx q[0];
rz(0.066083834) q[0];
x q[1];
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
rz(0.3469997) q[1];
sx q[1];
rz(-1.1113249) q[1];
sx q[1];
rz(-0.32230349) q[1];
x q[2];
rz(-2.4160556) q[3];
sx q[3];
rz(-2.7751338) q[3];
sx q[3];
rz(3.0612502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6254639) q[2];
sx q[2];
rz(-1.4279782) q[2];
sx q[2];
rz(3.1141172) q[2];
rz(0.52250683) q[3];
sx q[3];
rz(-2.3404739) q[3];
sx q[3];
rz(0.81645042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-2.7291173) q[0];
sx q[0];
rz(-2.0232047) q[0];
sx q[0];
rz(3.0274042) q[0];
rz(0.97822899) q[1];
sx q[1];
rz(-2.5949635) q[1];
sx q[1];
rz(0.79089975) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1066125) q[0];
sx q[0];
rz(-2.5133586) q[0];
sx q[0];
rz(-1.4094704) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3979264) q[2];
sx q[2];
rz(-1.3197834) q[2];
sx q[2];
rz(-0.92168346) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.32614) q[1];
sx q[1];
rz(-0.37323144) q[1];
sx q[1];
rz(-1.3532072) q[1];
x q[2];
rz(-0.88364717) q[3];
sx q[3];
rz(-1.8568608) q[3];
sx q[3];
rz(1.4253695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.87166446) q[2];
sx q[2];
rz(-2.0285138) q[2];
sx q[2];
rz(2.2765735) q[2];
rz(-2.4690752) q[3];
sx q[3];
rz(-2.0130242) q[3];
sx q[3];
rz(-3.045936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4928116) q[0];
sx q[0];
rz(-0.5195986) q[0];
sx q[0];
rz(-2.6742324) q[0];
rz(-0.53721792) q[1];
sx q[1];
rz(-2.1610114) q[1];
sx q[1];
rz(0.25407243) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92914903) q[0];
sx q[0];
rz(-1.7663167) q[0];
sx q[0];
rz(-1.6385727) q[0];
x q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(2.4916363) q[1];
sx q[1];
rz(-1.8815814) q[1];
sx q[1];
rz(-1.4375356) q[1];
rz(-pi) q[2];
rz(0.11717511) q[3];
sx q[3];
rz(-0.60248884) q[3];
sx q[3];
rz(-2.6206827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.49999985) q[2];
sx q[2];
rz(-0.77074146) q[2];
sx q[2];
rz(-0.33561486) q[2];
rz(0.32661682) q[3];
sx q[3];
rz(-0.89544046) q[3];
sx q[3];
rz(1.7232822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0582054) q[0];
sx q[0];
rz(-0.21033062) q[0];
sx q[0];
rz(-1.0539508) q[0];
rz(-0.15696934) q[1];
sx q[1];
rz(-1.7228246) q[1];
sx q[1];
rz(-0.98186791) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4440585) q[0];
sx q[0];
rz(-0.42675787) q[0];
sx q[0];
rz(2.8696637) q[0];
rz(-pi) q[1];
rz(2.1397212) q[2];
sx q[2];
rz(-1.47942) q[2];
sx q[2];
rz(-0.56359529) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8357504) q[1];
sx q[1];
rz(-1.7562477) q[1];
sx q[1];
rz(-0.058538392) q[1];
rz(-1.8198265) q[3];
sx q[3];
rz(-0.8753652) q[3];
sx q[3];
rz(1.4679366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2173569) q[2];
sx q[2];
rz(-1.057426) q[2];
sx q[2];
rz(-0.99739933) q[2];
rz(-2.635397) q[3];
sx q[3];
rz(-2.180407) q[3];
sx q[3];
rz(-3.1072646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85703325) q[0];
sx q[0];
rz(-2.4301346) q[0];
sx q[0];
rz(-2.6249264) q[0];
rz(2.754028) q[1];
sx q[1];
rz(-2.0811847) q[1];
sx q[1];
rz(0.26836747) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29753387) q[0];
sx q[0];
rz(-2.2608272) q[0];
sx q[0];
rz(-0.78155078) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.43865135) q[2];
sx q[2];
rz(-1.029656) q[2];
sx q[2];
rz(-1.3676608) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.25071535) q[1];
sx q[1];
rz(-2.4836575) q[1];
sx q[1];
rz(0.4472181) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1770505) q[3];
sx q[3];
rz(-0.6150133) q[3];
sx q[3];
rz(-0.39214373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2051852) q[2];
sx q[2];
rz(-2.3582017) q[2];
sx q[2];
rz(2.5777585) q[2];
rz(1.9514203) q[3];
sx q[3];
rz(-2.1824013) q[3];
sx q[3];
rz(-2.4998375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47678369) q[0];
sx q[0];
rz(-1.8905147) q[0];
sx q[0];
rz(2.074194) q[0];
rz(-1.8019567) q[1];
sx q[1];
rz(-1.4383153) q[1];
sx q[1];
rz(-1.7972606) q[1];
rz(-2.6396991) q[2];
sx q[2];
rz(-2.5645651) q[2];
sx q[2];
rz(-1.2679451) q[2];
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