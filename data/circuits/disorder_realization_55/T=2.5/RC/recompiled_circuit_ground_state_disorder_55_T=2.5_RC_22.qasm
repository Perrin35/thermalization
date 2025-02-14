OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3352873) q[0];
sx q[0];
rz(-1.0399613) q[0];
sx q[0];
rz(-0.34997532) q[0];
rz(1.5179874) q[1];
sx q[1];
rz(4.0859695) q[1];
sx q[1];
rz(10.610217) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1040504) q[0];
sx q[0];
rz(-0.54706956) q[0];
sx q[0];
rz(2.4614835) q[0];
rz(0.41857206) q[2];
sx q[2];
rz(-1.3579733) q[2];
sx q[2];
rz(-0.71453071) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6373316) q[1];
sx q[1];
rz(-2.0954014) q[1];
sx q[1];
rz(0.27224147) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1864263) q[3];
sx q[3];
rz(-1.7209779) q[3];
sx q[3];
rz(2.3038441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8650633) q[2];
sx q[2];
rz(-0.41215602) q[2];
sx q[2];
rz(2.7318447) q[2];
rz(-2.3661738) q[3];
sx q[3];
rz(-1.6131468) q[3];
sx q[3];
rz(-2.8342136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2839117) q[0];
sx q[0];
rz(-2.6429521) q[0];
sx q[0];
rz(2.6890802) q[0];
rz(-0.2287989) q[1];
sx q[1];
rz(-2.6401873) q[1];
sx q[1];
rz(2.9948044) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2499224) q[0];
sx q[0];
rz(-1.2091214) q[0];
sx q[0];
rz(-0.2043496) q[0];
x q[1];
rz(1.8092621) q[2];
sx q[2];
rz(-1.9386282) q[2];
sx q[2];
rz(0.17300082) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2816946) q[1];
sx q[1];
rz(-1.4364151) q[1];
sx q[1];
rz(-1.6107537) q[1];
rz(-pi) q[2];
rz(2.4815219) q[3];
sx q[3];
rz(-1.1997031) q[3];
sx q[3];
rz(-2.0160272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.16810922) q[2];
sx q[2];
rz(-2.8066469) q[2];
sx q[2];
rz(-0.87146634) q[2];
rz(0.35428366) q[3];
sx q[3];
rz(-0.93577093) q[3];
sx q[3];
rz(1.3080477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63534921) q[0];
sx q[0];
rz(-2.5214218) q[0];
sx q[0];
rz(1.1059906) q[0];
rz(0.94331074) q[1];
sx q[1];
rz(-2.3425075) q[1];
sx q[1];
rz(-1.6518637) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7794253) q[0];
sx q[0];
rz(-0.12813103) q[0];
sx q[0];
rz(-1.4904899) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9331384) q[2];
sx q[2];
rz(-1.4770419) q[2];
sx q[2];
rz(1.1445647) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5691487) q[1];
sx q[1];
rz(-0.97578543) q[1];
sx q[1];
rz(1.5250564) q[1];
rz(-1.4964281) q[3];
sx q[3];
rz(-1.8083124) q[3];
sx q[3];
rz(0.46922153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6279383) q[2];
sx q[2];
rz(-1.1955806) q[2];
sx q[2];
rz(2.1453843) q[2];
rz(0.57322383) q[3];
sx q[3];
rz(-2.6300391) q[3];
sx q[3];
rz(0.63428026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2117677) q[0];
sx q[0];
rz(-1.8695762) q[0];
sx q[0];
rz(1.5856532) q[0];
rz(-2.8171483) q[1];
sx q[1];
rz(-0.8322081) q[1];
sx q[1];
rz(-1.9437887) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75173855) q[0];
sx q[0];
rz(-0.9969939) q[0];
sx q[0];
rz(1.874254) q[0];
x q[1];
rz(2.156331) q[2];
sx q[2];
rz(-1.6949495) q[2];
sx q[2];
rz(-2.2145766) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8363322) q[1];
sx q[1];
rz(-1.700749) q[1];
sx q[1];
rz(1.0856245) q[1];
rz(-pi) q[2];
rz(-2.1230202) q[3];
sx q[3];
rz(-0.30507824) q[3];
sx q[3];
rz(-2.6570002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.045362554) q[2];
sx q[2];
rz(-2.5924293) q[2];
sx q[2];
rz(1.9722923) q[2];
rz(-2.5176288) q[3];
sx q[3];
rz(-1.6922765) q[3];
sx q[3];
rz(0.72871488) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28155926) q[0];
sx q[0];
rz(-0.39880729) q[0];
sx q[0];
rz(-1.1997724) q[0];
rz(-1.8343605) q[1];
sx q[1];
rz(-2.2016826) q[1];
sx q[1];
rz(-1.7050381) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67302895) q[0];
sx q[0];
rz(-1.1832245) q[0];
sx q[0];
rz(-3.0762071) q[0];
rz(-pi) q[1];
x q[1];
rz(0.13135037) q[2];
sx q[2];
rz(-1.7997912) q[2];
sx q[2];
rz(-0.74818767) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5561394) q[1];
sx q[1];
rz(-1.2188497) q[1];
sx q[1];
rz(-1.2519678) q[1];
x q[2];
rz(0.9229696) q[3];
sx q[3];
rz(-1.2296835) q[3];
sx q[3];
rz(0.34527147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3522219) q[2];
sx q[2];
rz(-1.4381961) q[2];
sx q[2];
rz(2.5686725) q[2];
rz(0.88417435) q[3];
sx q[3];
rz(-2.1638162) q[3];
sx q[3];
rz(0.49527112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5837412) q[0];
sx q[0];
rz(-1.252625) q[0];
sx q[0];
rz(-2.0649233) q[0];
rz(2.63511) q[1];
sx q[1];
rz(-2.5264637) q[1];
sx q[1];
rz(-1.3004251) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5426555) q[0];
sx q[0];
rz(-3.0910203) q[0];
sx q[0];
rz(-1.845832) q[0];
rz(-pi) q[1];
rz(1.7188273) q[2];
sx q[2];
rz(-0.99642053) q[2];
sx q[2];
rz(0.21166052) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.36339736) q[1];
sx q[1];
rz(-1.7018688) q[1];
sx q[1];
rz(1.78472) q[1];
x q[2];
rz(-0.20008101) q[3];
sx q[3];
rz(-1.225138) q[3];
sx q[3];
rz(2.8576771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7252245) q[2];
sx q[2];
rz(-0.84364265) q[2];
sx q[2];
rz(-1.9492524) q[2];
rz(3.1277749) q[3];
sx q[3];
rz(-0.71607029) q[3];
sx q[3];
rz(-2.1684087) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.248812) q[0];
sx q[0];
rz(-2.2154494) q[0];
sx q[0];
rz(-2.7698351) q[0];
rz(0.88428503) q[1];
sx q[1];
rz(-2.6339032) q[1];
sx q[1];
rz(2.5977792) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.666558) q[0];
sx q[0];
rz(-2.130593) q[0];
sx q[0];
rz(-2.0535357) q[0];
rz(0.45525785) q[2];
sx q[2];
rz(-2.2374059) q[2];
sx q[2];
rz(-2.8316759) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4298599) q[1];
sx q[1];
rz(-1.5774343) q[1];
sx q[1];
rz(-0.37488519) q[1];
x q[2];
rz(3.032148) q[3];
sx q[3];
rz(-1.4196803) q[3];
sx q[3];
rz(1.5990638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1112993) q[2];
sx q[2];
rz(-1.0617504) q[2];
sx q[2];
rz(-0.57478762) q[2];
rz(2.287367) q[3];
sx q[3];
rz(-1.6653019) q[3];
sx q[3];
rz(-2.0265719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2146571) q[0];
sx q[0];
rz(-0.68100005) q[0];
sx q[0];
rz(-2.8171203) q[0];
rz(-0.60943162) q[1];
sx q[1];
rz(-0.84171265) q[1];
sx q[1];
rz(0.51469523) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1457532) q[0];
sx q[0];
rz(-1.1356562) q[0];
sx q[0];
rz(0.010220411) q[0];
x q[1];
rz(1.1616522) q[2];
sx q[2];
rz(-0.7557286) q[2];
sx q[2];
rz(-2.8314396) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0275484) q[1];
sx q[1];
rz(-2.1308772) q[1];
sx q[1];
rz(-1.8482988) q[1];
rz(-1.9402766) q[3];
sx q[3];
rz(-2.8469128) q[3];
sx q[3];
rz(-1.0436383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9784166) q[2];
sx q[2];
rz(-1.5055483) q[2];
sx q[2];
rz(-0.27463883) q[2];
rz(-0.93242532) q[3];
sx q[3];
rz(-0.43007389) q[3];
sx q[3];
rz(-2.8453804) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4579929) q[0];
sx q[0];
rz(-1.4608811) q[0];
sx q[0];
rz(-3.006111) q[0];
rz(2.1645891) q[1];
sx q[1];
rz(-2.3833387) q[1];
sx q[1];
rz(-0.0010842222) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8086804) q[0];
sx q[0];
rz(-1.4152193) q[0];
sx q[0];
rz(-0.498144) q[0];
x q[1];
rz(2.094467) q[2];
sx q[2];
rz(-1.4677248) q[2];
sx q[2];
rz(1.6873941) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.407775) q[1];
sx q[1];
rz(-2.1548591) q[1];
sx q[1];
rz(-2.255846) q[1];
rz(-pi) q[2];
rz(-2.3850609) q[3];
sx q[3];
rz(-2.0755141) q[3];
sx q[3];
rz(2.0677419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.5929426) q[2];
sx q[2];
rz(-2.103001) q[2];
sx q[2];
rz(-2.9366117) q[2];
rz(0.18260469) q[3];
sx q[3];
rz(-2.9349116) q[3];
sx q[3];
rz(-0.73777795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36271998) q[0];
sx q[0];
rz(-0.39411476) q[0];
sx q[0];
rz(2.6642098) q[0];
rz(-0.1167156) q[1];
sx q[1];
rz(-1.621403) q[1];
sx q[1];
rz(2.3027072) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4981094) q[0];
sx q[0];
rz(-2.2543397) q[0];
sx q[0];
rz(-0.34434005) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7493189) q[2];
sx q[2];
rz(-0.76425154) q[2];
sx q[2];
rz(1.4832254) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2005208) q[1];
sx q[1];
rz(-1.3844191) q[1];
sx q[1];
rz(2.446387) q[1];
rz(-pi) q[2];
rz(-1.8096883) q[3];
sx q[3];
rz(-1.5419067) q[3];
sx q[3];
rz(-0.98945557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.78187752) q[2];
sx q[2];
rz(-0.30656591) q[2];
sx q[2];
rz(0.26620418) q[2];
rz(-2.6340458) q[3];
sx q[3];
rz(-2.4188953) q[3];
sx q[3];
rz(2.7005196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97892852) q[0];
sx q[0];
rz(-1.4586108) q[0];
sx q[0];
rz(2.3391531) q[0];
rz(2.2294527) q[1];
sx q[1];
rz(-1.1911387) q[1];
sx q[1];
rz(-2.3008507) q[1];
rz(1.4100762) q[2];
sx q[2];
rz(-1.3257324) q[2];
sx q[2];
rz(-0.64468862) q[2];
rz(-1.9515362) q[3];
sx q[3];
rz(-1.1968812) q[3];
sx q[3];
rz(1.4637917) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
