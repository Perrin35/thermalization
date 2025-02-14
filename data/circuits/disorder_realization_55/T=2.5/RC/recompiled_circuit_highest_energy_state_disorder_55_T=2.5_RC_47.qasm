OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.363938) q[0];
sx q[0];
rz(3.5600297) q[0];
sx q[0];
rz(10.726396) q[0];
rz(-2.9786181) q[1];
sx q[1];
rz(1.4459223) q[1];
sx q[1];
rz(12.583153) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34384511) q[0];
sx q[0];
rz(-0.36917403) q[0];
sx q[0];
rz(-1.3626211) q[0];
rz(1.7503337) q[2];
sx q[2];
rz(-1.554536) q[2];
sx q[2];
rz(-0.017841466) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.7960733) q[1];
sx q[1];
rz(-1.5757478) q[1];
sx q[1];
rz(1.5657182) q[1];
x q[2];
rz(2.7990325) q[3];
sx q[3];
rz(-0.16157074) q[3];
sx q[3];
rz(-1.8752961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5825384) q[2];
sx q[2];
rz(-0.6607008) q[2];
sx q[2];
rz(1.5793229) q[2];
rz(2.9301379) q[3];
sx q[3];
rz(-0.00051694218) q[3];
sx q[3];
rz(2.9805984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.529539) q[0];
sx q[0];
rz(-2.8654629) q[0];
sx q[0];
rz(-1.8147234) q[0];
rz(-0.57948411) q[1];
sx q[1];
rz(-0.0038298413) q[1];
sx q[1];
rz(2.5025867) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76692894) q[0];
sx q[0];
rz(-2.2013454) q[0];
sx q[0];
rz(2.2847644) q[0];
rz(3.0045549) q[2];
sx q[2];
rz(-3.0195022) q[2];
sx q[2];
rz(-3.0220495) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9663869) q[1];
sx q[1];
rz(-1.5537973) q[1];
sx q[1];
rz(-1.5592038) q[1];
rz(3.0538179) q[3];
sx q[3];
rz(-0.78237659) q[3];
sx q[3];
rz(-2.6915472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4352033) q[2];
sx q[2];
rz(-3.0053164) q[2];
sx q[2];
rz(1.603568) q[2];
rz(-1.5806574) q[3];
sx q[3];
rz(-3.1272562) q[3];
sx q[3];
rz(3.1106136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.2494217) q[0];
sx q[0];
rz(-0.51512655) q[0];
sx q[0];
rz(2.7601335) q[0];
rz(-0.70746607) q[1];
sx q[1];
rz(-0.019376945) q[1];
sx q[1];
rz(-2.0170508) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1046421) q[0];
sx q[0];
rz(-0.35250394) q[0];
sx q[0];
rz(-0.80926766) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3493645) q[2];
sx q[2];
rz(-0.11781684) q[2];
sx q[2];
rz(0.066514579) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2644653) q[1];
sx q[1];
rz(-1.6337602) q[1];
sx q[1];
rz(1.5829045) q[1];
rz(1.0550523) q[3];
sx q[3];
rz(-2.0655144) q[3];
sx q[3];
rz(-0.62913857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6996998) q[2];
sx q[2];
rz(-0.012233891) q[2];
sx q[2];
rz(0.049467889) q[2];
rz(2.5345645) q[3];
sx q[3];
rz(-0.0012461239) q[3];
sx q[3];
rz(1.9342669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98168755) q[0];
sx q[0];
rz(-2.9732381) q[0];
sx q[0];
rz(-0.017177563) q[0];
rz(2.848564) q[1];
sx q[1];
rz(-0.79048645) q[1];
sx q[1];
rz(-1.5471829) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38774224) q[0];
sx q[0];
rz(-2.1374636) q[0];
sx q[0];
rz(2.500534) q[0];
rz(1.3457005) q[2];
sx q[2];
rz(-1.0402816) q[2];
sx q[2];
rz(-0.38166416) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3354825) q[1];
sx q[1];
rz(-1.6948189) q[1];
sx q[1];
rz(1.5107442) q[1];
rz(-pi) q[2];
rz(-1.5805575) q[3];
sx q[3];
rz(-0.92484821) q[3];
sx q[3];
rz(0.6432337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8881417) q[2];
sx q[2];
rz(-2.6710822) q[2];
sx q[2];
rz(-0.68224254) q[2];
rz(0.055179723) q[3];
sx q[3];
rz(-3.1339055) q[3];
sx q[3];
rz(-1.8365708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-2.3799915) q[0];
sx q[0];
rz(-2.70607) q[0];
sx q[0];
rz(-0.4250266) q[0];
rz(1.60166) q[1];
sx q[1];
rz(-0.48307499) q[1];
sx q[1];
rz(-0.81533283) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7777625) q[0];
sx q[0];
rz(-1.6319931) q[0];
sx q[0];
rz(-1.5817002) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0032665423) q[2];
sx q[2];
rz(-1.5940394) q[2];
sx q[2];
rz(-0.88752623) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1734577) q[1];
sx q[1];
rz(-1.6575529) q[1];
sx q[1];
rz(-0.11958338) q[1];
rz(-pi) q[2];
rz(2.7975122) q[3];
sx q[3];
rz(-1.6852143) q[3];
sx q[3];
rz(-2.2215869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1347947) q[2];
sx q[2];
rz(-0.012601348) q[2];
sx q[2];
rz(1.4726144) q[2];
rz(-0.93004477) q[3];
sx q[3];
rz(-3.127122) q[3];
sx q[3];
rz(-2.3013733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8188266) q[0];
sx q[0];
rz(-3.1184986) q[0];
sx q[0];
rz(1.7148788) q[0];
rz(2.4115883) q[1];
sx q[1];
rz(-0.58861029) q[1];
sx q[1];
rz(2.0402562) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5752714) q[0];
sx q[0];
rz(-0.76714095) q[0];
sx q[0];
rz(-1.1772035) q[0];
rz(1.3404714) q[2];
sx q[2];
rz(-1.7353976) q[2];
sx q[2];
rz(2.4245461) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7633865) q[1];
sx q[1];
rz(-2.9855033) q[1];
sx q[1];
rz(0.915145) q[1];
rz(-pi) q[2];
rz(-1.3188478) q[3];
sx q[3];
rz(-0.12853208) q[3];
sx q[3];
rz(2.4737308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4250028) q[2];
sx q[2];
rz(-3.0807107) q[2];
sx q[2];
rz(1.8343743) q[2];
rz(2.763125) q[3];
sx q[3];
rz(-3.1186447) q[3];
sx q[3];
rz(2.4881261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3017479) q[0];
sx q[0];
rz(-1.2675588) q[0];
sx q[0];
rz(2.2140455) q[0];
rz(1.357366) q[1];
sx q[1];
rz(-0.83186847) q[1];
sx q[1];
rz(1.6102788) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0831628) q[0];
sx q[0];
rz(-1.6091804) q[0];
sx q[0];
rz(-3.1192664) q[0];
rz(-pi) q[1];
x q[1];
rz(0.39674098) q[2];
sx q[2];
rz(-1.7974241) q[2];
sx q[2];
rz(0.47725783) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5824052) q[1];
sx q[1];
rz(-1.4417164) q[1];
sx q[1];
rz(0.0027133769) q[1];
rz(-pi) q[2];
rz(1.0590943) q[3];
sx q[3];
rz(-1.8476433) q[3];
sx q[3];
rz(0.83864318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.36074582) q[2];
sx q[2];
rz(-0.0043892269) q[2];
sx q[2];
rz(-1.2537664) q[2];
rz(0.69416657) q[3];
sx q[3];
rz(-0.73918754) q[3];
sx q[3];
rz(0.23301253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53741443) q[0];
sx q[0];
rz(-2.1359213) q[0];
sx q[0];
rz(-1.0373212) q[0];
rz(1.6088156) q[1];
sx q[1];
rz(-2.9210126) q[1];
sx q[1];
rz(1.6745837) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4983385) q[0];
sx q[0];
rz(-2.9349929) q[0];
sx q[0];
rz(-1.905196) q[0];
rz(-2.8634818) q[2];
sx q[2];
rz(-1.5486954) q[2];
sx q[2];
rz(1.8430361) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8407708) q[1];
sx q[1];
rz(-1.5696708) q[1];
sx q[1];
rz(0.00043934396) q[1];
x q[2];
rz(-2.6862304) q[3];
sx q[3];
rz(-1.6557367) q[3];
sx q[3];
rz(-2.7698539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.011270114) q[2];
sx q[2];
rz(-2.932817) q[2];
sx q[2];
rz(-0.043896349) q[2];
rz(-0.52055001) q[3];
sx q[3];
rz(-0.0046516727) q[3];
sx q[3];
rz(-1.4588149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3540102) q[0];
sx q[0];
rz(-3.1391322) q[0];
sx q[0];
rz(1.3186697) q[0];
rz(-1.4175381) q[1];
sx q[1];
rz(-0.28957614) q[1];
sx q[1];
rz(1.5971378) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8108687) q[0];
sx q[0];
rz(-1.6858584) q[0];
sx q[0];
rz(2.6601391) q[0];
rz(-pi) q[1];
rz(-0.33716069) q[2];
sx q[2];
rz(-2.4558407) q[2];
sx q[2];
rz(-2.8711988) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.456157) q[1];
sx q[1];
rz(-0.36000571) q[1];
sx q[1];
rz(-2.092157) q[1];
rz(-pi) q[2];
x q[2];
rz(0.61267743) q[3];
sx q[3];
rz(-0.78734382) q[3];
sx q[3];
rz(-1.2932216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3370207) q[2];
sx q[2];
rz(-1.8420409) q[2];
sx q[2];
rz(-0.19680944) q[2];
rz(-1.9491516) q[3];
sx q[3];
rz(-2.9350023) q[3];
sx q[3];
rz(2.9406252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30109677) q[0];
sx q[0];
rz(-1.3571285) q[0];
sx q[0];
rz(1.1916196) q[0];
rz(1.6169029) q[1];
sx q[1];
rz(-2.4947417) q[1];
sx q[1];
rz(1.5764538) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8218481) q[0];
sx q[0];
rz(-0.34043202) q[0];
sx q[0];
rz(2.6186187) q[0];
rz(-0.2290957) q[2];
sx q[2];
rz(-2.014403) q[2];
sx q[2];
rz(-2.6227621) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8818156) q[1];
sx q[1];
rz(-1.5702973) q[1];
sx q[1];
rz(-3.1406443) q[1];
rz(-pi) q[2];
rz(1.7129536) q[3];
sx q[3];
rz(-1.6831493) q[3];
sx q[3];
rz(-1.6522828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9578751) q[2];
sx q[2];
rz(-2.5536733) q[2];
sx q[2];
rz(1.4584165) q[2];
rz(0.030473907) q[3];
sx q[3];
rz(-3.1320429) q[3];
sx q[3];
rz(-2.9392346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.6644345) q[0];
sx q[0];
rz(-1.3344593) q[0];
sx q[0];
rz(1.6819171) q[0];
rz(1.5741813) q[1];
sx q[1];
rz(-1.8125143) q[1];
sx q[1];
rz(0.090851091) q[1];
rz(1.5051928) q[2];
sx q[2];
rz(-3.0658683) q[2];
sx q[2];
rz(0.22584596) q[2];
rz(0.81672698) q[3];
sx q[3];
rz(-1.1882267) q[3];
sx q[3];
rz(-1.4878275) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
