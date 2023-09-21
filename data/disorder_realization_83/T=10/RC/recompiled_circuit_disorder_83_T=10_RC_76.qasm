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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0640472) q[0];
sx q[0];
rz(-0.24928688) q[0];
sx q[0];
rz(-1.3011978) q[0];
x q[1];
rz(-0.15723575) q[2];
sx q[2];
rz(-0.6586282) q[2];
sx q[2];
rz(0.14070357) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.57063369) q[1];
sx q[1];
rz(-1.966241) q[1];
sx q[1];
rz(2.8407437) q[1];
rz(-pi) q[2];
rz(-1.7573962) q[3];
sx q[3];
rz(-0.36986923) q[3];
sx q[3];
rz(2.6922525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7068229) q[2];
sx q[2];
rz(-0.40199026) q[2];
sx q[2];
rz(0.10786954) q[2];
rz(-0.14262959) q[3];
sx q[3];
rz(-1.4167891) q[3];
sx q[3];
rz(-2.4690348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07664872) q[0];
sx q[0];
rz(-2.3363484) q[0];
sx q[0];
rz(-2.9108677) q[0];
rz(1.8143066) q[1];
sx q[1];
rz(-0.67115152) q[1];
sx q[1];
rz(-3.1006295) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7777268) q[0];
sx q[0];
rz(-1.7330609) q[0];
sx q[0];
rz(-1.8198387) q[0];
rz(-pi) q[1];
rz(2.3054753) q[2];
sx q[2];
rz(-0.40514075) q[2];
sx q[2];
rz(0.9678313) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6662308) q[1];
sx q[1];
rz(-1.0050217) q[1];
sx q[1];
rz(-2.1192141) q[1];
rz(1.9713692) q[3];
sx q[3];
rz(-1.922732) q[3];
sx q[3];
rz(-0.34917253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.162398) q[2];
sx q[2];
rz(-1.6684063) q[2];
sx q[2];
rz(-0.56817788) q[2];
rz(0.5125106) q[3];
sx q[3];
rz(-2.5458702) q[3];
sx q[3];
rz(-2.5797243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5325461) q[0];
sx q[0];
rz(-0.17621528) q[0];
sx q[0];
rz(0.30058582) q[0];
rz(0.1707503) q[2];
sx q[2];
rz(-3.1351334) q[2];
sx q[2];
rz(3.1250172) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3213768) q[1];
sx q[1];
rz(-1.4396588) q[1];
sx q[1];
rz(-1.3908435) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.67218353) q[3];
sx q[3];
rz(-0.58478343) q[3];
sx q[3];
rz(-0.72090805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.21851097) q[2];
sx q[2];
rz(-1.9873025) q[2];
sx q[2];
rz(-0.046860524) q[2];
rz(2.3299407) q[3];
sx q[3];
rz(-0.67957687) q[3];
sx q[3];
rz(2.9714382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1470404) q[0];
sx q[0];
rz(-0.11163286) q[0];
sx q[0];
rz(3.0901093) q[0];
rz(-0.4908081) q[1];
sx q[1];
rz(-2.1665116) q[1];
sx q[1];
rz(-1.1725918) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5347036) q[0];
sx q[0];
rz(-2.4692315) q[0];
sx q[0];
rz(-0.21557233) q[0];
rz(-pi) q[1];
rz(-2.5355859) q[2];
sx q[2];
rz(-1.5835985) q[2];
sx q[2];
rz(0.16773227) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6231411) q[1];
sx q[1];
rz(-1.8703096) q[1];
sx q[1];
rz(-1.2007984) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.138932) q[3];
sx q[3];
rz(-1.5032839) q[3];
sx q[3];
rz(2.1907012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.88901687) q[2];
sx q[2];
rz(-0.93095195) q[2];
sx q[2];
rz(0.48689294) q[2];
rz(-1.1934818) q[3];
sx q[3];
rz(-2.8119757) q[3];
sx q[3];
rz(-3.1161599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30381969) q[0];
sx q[0];
rz(-1.0205512) q[0];
sx q[0];
rz(1.3254962) q[0];
rz(0.59108132) q[1];
sx q[1];
rz(-2.4609844) q[1];
sx q[1];
rz(-2.4868734) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3723719) q[0];
sx q[0];
rz(-1.1290316) q[0];
sx q[0];
rz(-1.9584993) q[0];
x q[1];
rz(-1.3833369) q[2];
sx q[2];
rz(-1.1584632) q[2];
sx q[2];
rz(0.86134855) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4795585) q[1];
sx q[1];
rz(-0.91842945) q[1];
sx q[1];
rz(2.0550834) q[1];
rz(-1.4773024) q[3];
sx q[3];
rz(-1.1870541) q[3];
sx q[3];
rz(-1.547471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6902265) q[2];
sx q[2];
rz(-1.2008685) q[2];
sx q[2];
rz(-0.5711242) q[2];
rz(0.58978224) q[3];
sx q[3];
rz(-2.6815806) q[3];
sx q[3];
rz(3.0736198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1968483) q[0];
sx q[0];
rz(-1.1905043) q[0];
sx q[0];
rz(-0.29904547) q[0];
rz(-1.8213182) q[1];
sx q[1];
rz(-0.25779217) q[1];
sx q[1];
rz(-1.6437795) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90692524) q[0];
sx q[0];
rz(-1.6229796) q[0];
sx q[0];
rz(-2.232057) q[0];
rz(-2.3408893) q[2];
sx q[2];
rz(-2.3869037) q[2];
sx q[2];
rz(-2.961219) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7945929) q[1];
sx q[1];
rz(-2.0302677) q[1];
sx q[1];
rz(0.32230349) q[1];
rz(1.8201581) q[3];
sx q[3];
rz(-1.8421838) q[3];
sx q[3];
rz(2.4621778) q[3];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41247535) q[0];
sx q[0];
rz(-1.1183879) q[0];
sx q[0];
rz(-0.11418848) q[0];
rz(-0.97822899) q[1];
sx q[1];
rz(-0.54662919) q[1];
sx q[1];
rz(-2.3506929) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9081156) q[0];
sx q[0];
rz(-0.95196264) q[0];
sx q[0];
rz(-3.0254362) q[0];
x q[1];
rz(-1.2355455) q[2];
sx q[2];
rz(-2.2860048) q[2];
sx q[2];
rz(-2.7170979) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0930209) q[1];
sx q[1];
rz(-1.2067817) q[1];
sx q[1];
rz(3.0572592) q[1];
rz(0.88364717) q[3];
sx q[3];
rz(-1.8568608) q[3];
sx q[3];
rz(1.7162232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2699282) q[2];
sx q[2];
rz(-1.1130788) q[2];
sx q[2];
rz(0.86501914) q[2];
rz(2.4690752) q[3];
sx q[3];
rz(-1.1285684) q[3];
sx q[3];
rz(0.095656693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4928116) q[0];
sx q[0];
rz(-0.5195986) q[0];
sx q[0];
rz(-0.46736026) q[0];
rz(-0.53721792) q[1];
sx q[1];
rz(-2.1610114) q[1];
sx q[1];
rz(0.25407243) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2652889) q[0];
sx q[0];
rz(-2.9348001) q[0];
sx q[0];
rz(2.8121023) q[0];
rz(-1.569283) q[2];
sx q[2];
rz(-1.5655893) q[2];
sx q[2];
rz(-1.7092012) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1797807) q[1];
sx q[1];
rz(-1.4439549) q[1];
sx q[1];
rz(2.8282053) q[1];
rz(-2.5423126) q[3];
sx q[3];
rz(-1.6370956) q[3];
sx q[3];
rz(-1.9950206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6415928) q[2];
sx q[2];
rz(-0.77074146) q[2];
sx q[2];
rz(0.33561486) q[2];
rz(-0.32661682) q[3];
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
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.083387233) q[0];
sx q[0];
rz(-2.931262) q[0];
sx q[0];
rz(-1.0539508) q[0];
rz(-2.9846233) q[1];
sx q[1];
rz(-1.418768) q[1];
sx q[1];
rz(2.1597247) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4440585) q[0];
sx q[0];
rz(-2.7148348) q[0];
sx q[0];
rz(0.27192893) q[0];
x q[1];
rz(-2.1397212) q[2];
sx q[2];
rz(-1.6621727) q[2];
sx q[2];
rz(-0.56359529) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5280142) q[1];
sx q[1];
rz(-2.947223) q[1];
sx q[1];
rz(-1.8730875) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4307067) q[3];
sx q[3];
rz(-1.7611739) q[3];
sx q[3];
rz(-3.0829317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2173569) q[2];
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
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85703325) q[0];
sx q[0];
rz(-0.71145809) q[0];
sx q[0];
rz(2.6249264) q[0];
rz(-0.38756469) q[1];
sx q[1];
rz(-2.0811847) q[1];
sx q[1];
rz(-2.8732252) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.836652) q[0];
sx q[0];
rz(-0.9965082) q[0];
sx q[0];
rz(-0.71026295) q[0];
rz(-pi) q[1];
x q[1];
rz(0.43865135) q[2];
sx q[2];
rz(-1.029656) q[2];
sx q[2];
rz(-1.7739319) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.25071535) q[1];
sx q[1];
rz(-2.4836575) q[1];
sx q[1];
rz(2.6943745) q[1];
rz(-0.99276944) q[3];
sx q[3];
rz(-1.3475932) q[3];
sx q[3];
rz(0.85152599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2051852) q[2];
sx q[2];
rz(-2.3582017) q[2];
sx q[2];
rz(2.5777585) q[2];
rz(1.9514203) q[3];
sx q[3];
rz(-0.95919132) q[3];
sx q[3];
rz(-0.64175516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
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
rz(1.2673169) q[2];
sx q[2];
rz(-2.0694642) q[2];
sx q[2];
rz(2.4533761) q[2];
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