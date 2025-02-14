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
rz(-0.37201878) q[0];
sx q[0];
rz(-2.7899185) q[0];
sx q[0];
rz(-0.055211842) q[0];
rz(1.4456324) q[1];
sx q[1];
rz(-0.90336019) q[1];
sx q[1];
rz(3.0076495) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75823821) q[0];
sx q[0];
rz(-1.7899492) q[0];
sx q[0];
rz(-1.4683649) q[0];
rz(-pi) q[1];
rz(-0.29965286) q[2];
sx q[2];
rz(-2.5404394) q[2];
sx q[2];
rz(0.047182949) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6559927) q[1];
sx q[1];
rz(-2.6645086) q[1];
sx q[1];
rz(2.361195) q[1];
rz(-3.1024801) q[3];
sx q[3];
rz(-1.0457504) q[3];
sx q[3];
rz(2.0994179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.17406164) q[2];
sx q[2];
rz(-1.920819) q[2];
sx q[2];
rz(-0.81452149) q[2];
rz(-0.19541611) q[3];
sx q[3];
rz(-0.82868367) q[3];
sx q[3];
rz(1.0403847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.279351) q[0];
sx q[0];
rz(-0.26222721) q[0];
sx q[0];
rz(2.1694515) q[0];
rz(-2.5402918) q[1];
sx q[1];
rz(-2.1763132) q[1];
sx q[1];
rz(-1.345529) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1336254) q[0];
sx q[0];
rz(-2.1012126) q[0];
sx q[0];
rz(-2.6625457) q[0];
rz(-pi) q[1];
rz(-2.7072979) q[2];
sx q[2];
rz(-2.5611097) q[2];
sx q[2];
rz(0.24809294) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7899155) q[1];
sx q[1];
rz(-0.43444217) q[1];
sx q[1];
rz(1.2277568) q[1];
rz(-pi) q[2];
rz(2.9604994) q[3];
sx q[3];
rz(-1.5646184) q[3];
sx q[3];
rz(2.1712077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.95416516) q[2];
sx q[2];
rz(-1.8623872) q[2];
sx q[2];
rz(0.97174755) q[2];
rz(1.0176954) q[3];
sx q[3];
rz(-1.965799) q[3];
sx q[3];
rz(-0.051232256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.510842) q[0];
sx q[0];
rz(-2.172281) q[0];
sx q[0];
rz(0.72072679) q[0];
rz(3.0156056) q[1];
sx q[1];
rz(-2.4469913) q[1];
sx q[1];
rz(-0.40649498) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46716786) q[0];
sx q[0];
rz(-2.1851288) q[0];
sx q[0];
rz(-2.7697639) q[0];
x q[1];
rz(0.5936908) q[2];
sx q[2];
rz(-0.83752692) q[2];
sx q[2];
rz(1.6551203) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5169591) q[1];
sx q[1];
rz(-2.0849094) q[1];
sx q[1];
rz(2.2463283) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8536989) q[3];
sx q[3];
rz(-1.9215309) q[3];
sx q[3];
rz(-0.76606295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5410109) q[2];
sx q[2];
rz(-1.7966248) q[2];
sx q[2];
rz(0.37008944) q[2];
rz(-0.078941405) q[3];
sx q[3];
rz(-2.168096) q[3];
sx q[3];
rz(2.6551042) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9321891) q[0];
sx q[0];
rz(-1.7904733) q[0];
sx q[0];
rz(-0.47602794) q[0];
rz(1.1141106) q[1];
sx q[1];
rz(-0.94534355) q[1];
sx q[1];
rz(1.6260737) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1080804) q[0];
sx q[0];
rz(-0.42244222) q[0];
sx q[0];
rz(2.6080934) q[0];
rz(0.8117746) q[2];
sx q[2];
rz(-2.0305131) q[2];
sx q[2];
rz(0.22961337) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.042121902) q[1];
sx q[1];
rz(-1.9212706) q[1];
sx q[1];
rz(1.480669) q[1];
rz(-pi) q[2];
rz(-1.2021121) q[3];
sx q[3];
rz(-2.3895309) q[3];
sx q[3];
rz(-0.78898174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.49426207) q[2];
sx q[2];
rz(-1.6202972) q[2];
sx q[2];
rz(-0.94804135) q[2];
rz(2.7028132) q[3];
sx q[3];
rz(-0.56883562) q[3];
sx q[3];
rz(-0.89890629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
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
rz(0.62094837) q[0];
sx q[0];
rz(-2.181894) q[0];
sx q[0];
rz(0.4278675) q[0];
rz(2.1828792) q[1];
sx q[1];
rz(-1.9045279) q[1];
sx q[1];
rz(-1.0520891) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3593345) q[0];
sx q[0];
rz(-2.3365031) q[0];
sx q[0];
rz(0.87053086) q[0];
x q[1];
rz(2.1053931) q[2];
sx q[2];
rz(-1.9048759) q[2];
sx q[2];
rz(0.47185791) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3173988) q[1];
sx q[1];
rz(-1.6643859) q[1];
sx q[1];
rz(0.99653523) q[1];
rz(-pi) q[2];
rz(0.73486272) q[3];
sx q[3];
rz(-0.70847337) q[3];
sx q[3];
rz(2.3264309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.170257) q[2];
sx q[2];
rz(-1.3054138) q[2];
sx q[2];
rz(-0.35856426) q[2];
rz(1.1809008) q[3];
sx q[3];
rz(-2.2660393) q[3];
sx q[3];
rz(2.6023279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81826687) q[0];
sx q[0];
rz(-1.8805255) q[0];
sx q[0];
rz(-2.4928424) q[0];
rz(2.5541041) q[1];
sx q[1];
rz(-1.8193388) q[1];
sx q[1];
rz(2.9041451) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5013468) q[0];
sx q[0];
rz(-3.0398439) q[0];
sx q[0];
rz(-0.32008381) q[0];
rz(-pi) q[1];
x q[1];
rz(0.35570972) q[2];
sx q[2];
rz(-0.28537073) q[2];
sx q[2];
rz(-1.4089546) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9623649) q[1];
sx q[1];
rz(-0.74848524) q[1];
sx q[1];
rz(-0.039123936) q[1];
rz(-pi) q[2];
rz(-2.1319207) q[3];
sx q[3];
rz(-2.3497407) q[3];
sx q[3];
rz(-0.71969024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4805523) q[2];
sx q[2];
rz(-0.66142267) q[2];
sx q[2];
rz(-2.1858369) q[2];
rz(1.3794948) q[3];
sx q[3];
rz(-0.62164128) q[3];
sx q[3];
rz(-0.1563589) q[3];
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
rz(1.5697923) q[0];
sx q[0];
rz(-0.0026230165) q[0];
sx q[0];
rz(3.0963335) q[0];
rz(2.3472002) q[1];
sx q[1];
rz(-2.2519799) q[1];
sx q[1];
rz(-0.20283595) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5595374) q[0];
sx q[0];
rz(-1.7235316) q[0];
sx q[0];
rz(-0.57270292) q[0];
rz(-2.4100609) q[2];
sx q[2];
rz(-1.712633) q[2];
sx q[2];
rz(2.518613) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6306298) q[1];
sx q[1];
rz(-1.65224) q[1];
sx q[1];
rz(2.5982473) q[1];
rz(1.6129812) q[3];
sx q[3];
rz(-0.72244553) q[3];
sx q[3];
rz(2.9056321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5016735) q[2];
sx q[2];
rz(-1.7379652) q[2];
sx q[2];
rz(0.61140927) q[2];
rz(-1.1008788) q[3];
sx q[3];
rz(-0.63488638) q[3];
sx q[3];
rz(-2.8618405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6708577) q[0];
sx q[0];
rz(-1.9453456) q[0];
sx q[0];
rz(-2.0580976) q[0];
rz(0.96393839) q[1];
sx q[1];
rz(-2.0699392) q[1];
sx q[1];
rz(-2.6458157) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4300604) q[0];
sx q[0];
rz(-1.4976131) q[0];
sx q[0];
rz(-3.1203695) q[0];
x q[1];
rz(-2.9588885) q[2];
sx q[2];
rz(-1.2983381) q[2];
sx q[2];
rz(-2.8962816) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.116175) q[1];
sx q[1];
rz(-2.1350645) q[1];
sx q[1];
rz(-0.49617824) q[1];
rz(-pi) q[2];
rz(-3.1115948) q[3];
sx q[3];
rz(-0.98326937) q[3];
sx q[3];
rz(0.8647747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.44020161) q[2];
sx q[2];
rz(-2.3728366) q[2];
sx q[2];
rz(2.4526147) q[2];
rz(-1.6280599) q[3];
sx q[3];
rz(-0.83008927) q[3];
sx q[3];
rz(-1.0015063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5643519) q[0];
sx q[0];
rz(-1.5654726) q[0];
sx q[0];
rz(-2.054457) q[0];
rz(2.3785036) q[1];
sx q[1];
rz(-1.7440081) q[1];
sx q[1];
rz(-0.22689247) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.034982) q[0];
sx q[0];
rz(-0.7312432) q[0];
sx q[0];
rz(1.909622) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.82590286) q[2];
sx q[2];
rz(-2.7159198) q[2];
sx q[2];
rz(0.22944006) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.58851885) q[1];
sx q[1];
rz(-1.5344193) q[1];
sx q[1];
rz(-0.22471551) q[1];
rz(-2.2506892) q[3];
sx q[3];
rz(-2.4197787) q[3];
sx q[3];
rz(2.8310404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6202966) q[2];
sx q[2];
rz(-2.2038348) q[2];
sx q[2];
rz(2.1916981) q[2];
rz(-2.1096443) q[3];
sx q[3];
rz(-2.2622006) q[3];
sx q[3];
rz(0.43752813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6501605) q[0];
sx q[0];
rz(-3.037945) q[0];
sx q[0];
rz(-1.1520804) q[0];
rz(0.092747124) q[1];
sx q[1];
rz(-0.68789613) q[1];
sx q[1];
rz(-0.50382096) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9418056) q[0];
sx q[0];
rz(-1.5212987) q[0];
sx q[0];
rz(-1.0440396) q[0];
rz(-0.70561346) q[2];
sx q[2];
rz(-2.793987) q[2];
sx q[2];
rz(-0.068989601) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.42119869) q[1];
sx q[1];
rz(-0.47537801) q[1];
sx q[1];
rz(0.092751547) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9484048) q[3];
sx q[3];
rz(-1.2758299) q[3];
sx q[3];
rz(2.0279391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4656226) q[2];
sx q[2];
rz(-2.0268107) q[2];
sx q[2];
rz(-0.62391227) q[2];
rz(0.55656773) q[3];
sx q[3];
rz(-2.4534295) q[3];
sx q[3];
rz(-2.9437959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.9482166) q[0];
sx q[0];
rz(-2.2103136) q[0];
sx q[0];
rz(0.92934004) q[0];
rz(-0.14840645) q[1];
sx q[1];
rz(-1.8971309) q[1];
sx q[1];
rz(2.94577) q[1];
rz(1.6918626) q[2];
sx q[2];
rz(-1.1101223) q[2];
sx q[2];
rz(2.3157673) q[2];
rz(1.8529057) q[3];
sx q[3];
rz(-0.28280453) q[3];
sx q[3];
rz(-0.40550532) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
