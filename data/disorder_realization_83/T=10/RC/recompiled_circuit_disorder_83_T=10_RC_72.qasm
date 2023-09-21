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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24510358) q[0];
sx q[0];
rz(-1.5050383) q[0];
sx q[0];
rz(-1.8114281) q[0];
x q[1];
rz(-1.450199) q[2];
sx q[2];
rz(-0.92168027) q[2];
sx q[2];
rz(0.33855864) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.570959) q[1];
sx q[1];
rz(-1.966241) q[1];
sx q[1];
rz(-0.30084893) q[1];
rz(-1.9348014) q[3];
sx q[3];
rz(-1.5036821) q[3];
sx q[3];
rz(1.2957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7068229) q[2];
sx q[2];
rz(-2.7396024) q[2];
sx q[2];
rz(3.0337231) q[2];
rz(2.9989631) q[3];
sx q[3];
rz(-1.7248036) q[3];
sx q[3];
rz(2.4690348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0649439) q[0];
sx q[0];
rz(-0.80524421) q[0];
sx q[0];
rz(2.9108677) q[0];
rz(1.327286) q[1];
sx q[1];
rz(-0.67115152) q[1];
sx q[1];
rz(-0.040963106) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3687392) q[0];
sx q[0];
rz(-2.8452747) q[0];
sx q[0];
rz(0.98451891) q[0];
x q[1];
rz(2.3054753) q[2];
sx q[2];
rz(-2.7364519) q[2];
sx q[2];
rz(2.1737614) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6662308) q[1];
sx q[1];
rz(-1.0050217) q[1];
sx q[1];
rz(2.1192141) q[1];
rz(-pi) q[2];
rz(2.326194) q[3];
sx q[3];
rz(-2.6147463) q[3];
sx q[3];
rz(-2.6032053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.162398) q[2];
sx q[2];
rz(-1.6684063) q[2];
sx q[2];
rz(0.56817788) q[2];
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
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7115962) q[0];
sx q[0];
rz(-1.5189518) q[0];
sx q[0];
rz(-2.6913753) q[0];
rz(-1.2954767) q[1];
sx q[1];
rz(-1.1643012) q[1];
sx q[1];
rz(-2.4643262) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60904658) q[0];
sx q[0];
rz(-2.9653774) q[0];
sx q[0];
rz(2.8410068) q[0];
x q[1];
rz(0.1707503) q[2];
sx q[2];
rz(-0.0064592529) q[2];
sx q[2];
rz(0.016575459) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8202159) q[1];
sx q[1];
rz(-1.4396588) q[1];
sx q[1];
rz(1.7507491) q[1];
x q[2];
rz(-2.4694091) q[3];
sx q[3];
rz(-2.5568092) q[3];
sx q[3];
rz(2.4206846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.21851097) q[2];
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
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99455225) q[0];
sx q[0];
rz(-3.0299598) q[0];
sx q[0];
rz(-0.051483367) q[0];
rz(2.6507846) q[1];
sx q[1];
rz(-0.97508109) q[1];
sx q[1];
rz(1.1725918) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6068891) q[0];
sx q[0];
rz(-2.4692315) q[0];
sx q[0];
rz(0.21557233) q[0];
rz(-pi) q[1];
rz(-0.6060068) q[2];
sx q[2];
rz(-1.5835985) q[2];
sx q[2];
rz(-0.16773227) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7395775) q[1];
sx q[1];
rz(-0.47164729) q[1];
sx q[1];
rz(0.86400835) q[1];
rz(-pi) q[2];
rz(1.138932) q[3];
sx q[3];
rz(-1.5032839) q[3];
sx q[3];
rz(0.95089144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.88901687) q[2];
sx q[2];
rz(-0.93095195) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30381969) q[0];
sx q[0];
rz(-1.0205512) q[0];
sx q[0];
rz(-1.3254962) q[0];
rz(2.5505113) q[1];
sx q[1];
rz(-2.4609844) q[1];
sx q[1];
rz(2.4868734) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3723719) q[0];
sx q[0];
rz(-2.0125611) q[0];
sx q[0];
rz(-1.1830933) q[0];
rz(-pi) q[1];
rz(-0.40277092) q[2];
sx q[2];
rz(-0.45071128) q[2];
sx q[2];
rz(-2.7223301) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.66203413) q[1];
sx q[1];
rz(-0.91842945) q[1];
sx q[1];
rz(2.0550834) q[1];
rz(-pi) q[2];
rz(2.7563285) q[3];
sx q[3];
rz(-1.6574727) q[3];
sx q[3];
rz(0.011766089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4513662) q[2];
sx q[2];
rz(-1.2008685) q[2];
sx q[2];
rz(-2.5704685) q[2];
rz(2.5518104) q[3];
sx q[3];
rz(-0.46001205) q[3];
sx q[3];
rz(-0.067972876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94474435) q[0];
sx q[0];
rz(-1.1905043) q[0];
sx q[0];
rz(-2.8425472) q[0];
rz(-1.8213182) q[1];
sx q[1];
rz(-2.8838005) q[1];
sx q[1];
rz(1.6437795) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70443557) q[0];
sx q[0];
rz(-2.2309982) q[0];
sx q[0];
rz(0.066083834) q[0];
rz(0.80070337) q[2];
sx q[2];
rz(-0.75468894) q[2];
sx q[2];
rz(-0.18037361) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7945929) q[1];
sx q[1];
rz(-2.0302677) q[1];
sx q[1];
rz(0.32230349) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.27961126) q[3];
sx q[3];
rz(-1.8108484) q[3];
sx q[3];
rz(0.9595426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.51612878) q[2];
sx q[2];
rz(-1.4279782) q[2];
sx q[2];
rz(3.1141172) q[2];
rz(-0.52250683) q[3];
sx q[3];
rz(-0.80111879) q[3];
sx q[3];
rz(0.81645042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41247535) q[0];
sx q[0];
rz(-1.1183879) q[0];
sx q[0];
rz(-3.0274042) q[0];
rz(-2.1633637) q[1];
sx q[1];
rz(-2.5949635) q[1];
sx q[1];
rz(-2.3506929) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0349802) q[0];
sx q[0];
rz(-2.5133586) q[0];
sx q[0];
rz(-1.7321222) q[0];
rz(-0.74366624) q[2];
sx q[2];
rz(-1.3197834) q[2];
sx q[2];
rz(0.92168346) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6939047) q[1];
sx q[1];
rz(-1.649592) q[1];
sx q[1];
rz(1.2055956) q[1];
x q[2];
rz(0.88364717) q[3];
sx q[3];
rz(-1.2847319) q[3];
sx q[3];
rz(1.4253695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.87166446) q[2];
sx q[2];
rz(-2.0285138) q[2];
sx q[2];
rz(-2.2765735) q[2];
rz(0.67251742) q[3];
sx q[3];
rz(-1.1285684) q[3];
sx q[3];
rz(3.045936) q[3];
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
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6487811) q[0];
sx q[0];
rz(-0.5195986) q[0];
sx q[0];
rz(-2.6742324) q[0];
rz(2.6043747) q[1];
sx q[1];
rz(-0.98058128) q[1];
sx q[1];
rz(2.8875202) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2652889) q[0];
sx q[0];
rz(-2.9348001) q[0];
sx q[0];
rz(-0.32949038) q[0];
rz(-pi) q[1];
rz(3.1363856) q[2];
sx q[2];
rz(-1.569283) q[2];
sx q[2];
rz(3.0031799) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.64995631) q[1];
sx q[1];
rz(-1.2600113) q[1];
sx q[1];
rz(1.4375356) q[1];
rz(-pi) q[2];
rz(-1.4905606) q[3];
sx q[3];
rz(-0.9730162) q[3];
sx q[3];
rz(2.762592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-1.7232822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0582054) q[0];
sx q[0];
rz(-2.931262) q[0];
sx q[0];
rz(2.0876419) q[0];
rz(-2.9846233) q[1];
sx q[1];
rz(-1.418768) q[1];
sx q[1];
rz(2.1597247) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1218402) q[0];
sx q[0];
rz(-1.6822018) q[0];
sx q[0];
rz(-2.7287672) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0018714) q[2];
sx q[2];
rz(-1.6621727) q[2];
sx q[2];
rz(-0.56359529) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8874444) q[1];
sx q[1];
rz(-1.6283298) q[1];
sx q[1];
rz(1.385034) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.710886) q[3];
sx q[3];
rz(-1.3804187) q[3];
sx q[3];
rz(3.0829317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.92423576) q[2];
sx q[2];
rz(-2.0841667) q[2];
sx q[2];
rz(0.99739933) q[2];
rz(-0.50619566) q[3];
sx q[3];
rz(-0.96118569) q[3];
sx q[3];
rz(-3.1072646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.2845594) q[0];
sx q[0];
rz(-2.4301346) q[0];
sx q[0];
rz(-0.51666623) q[0];
rz(0.38756469) q[1];
sx q[1];
rz(-1.060408) q[1];
sx q[1];
rz(-2.8732252) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3049406) q[0];
sx q[0];
rz(-2.1450844) q[0];
sx q[0];
rz(0.71026295) q[0];
x q[1];
rz(2.1568314) q[2];
sx q[2];
rz(-1.9433937) q[2];
sx q[2];
rz(-0.44024703) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8908773) q[1];
sx q[1];
rz(-0.6579352) q[1];
sx q[1];
rz(-2.6943745) q[1];
rz(2.8769365) q[3];
sx q[3];
rz(-2.1327244) q[3];
sx q[3];
rz(-2.2789126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9364075) q[2];
sx q[2];
rz(-2.3582017) q[2];
sx q[2];
rz(0.56383413) q[2];
rz(-1.9514203) q[3];
sx q[3];
rz(-0.95919132) q[3];
sx q[3];
rz(-2.4998375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47678369) q[0];
sx q[0];
rz(-1.2510779) q[0];
sx q[0];
rz(-1.0673987) q[0];
rz(1.339636) q[1];
sx q[1];
rz(-1.4383153) q[1];
sx q[1];
rz(-1.7972606) q[1];
rz(-2.6230326) q[2];
sx q[2];
rz(-1.8363562) q[2];
sx q[2];
rz(0.73391757) q[2];
rz(0.31233882) q[3];
sx q[3];
rz(-1.7620419) q[3];
sx q[3];
rz(-0.077802303) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];