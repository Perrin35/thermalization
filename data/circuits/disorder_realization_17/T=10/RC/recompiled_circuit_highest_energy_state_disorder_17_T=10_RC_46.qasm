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
rz(1.849527) q[0];
sx q[0];
rz(5.4592291) q[0];
sx q[0];
rz(8.4973314) q[0];
rz(-1.0499586) q[1];
sx q[1];
rz(-2.1702622) q[1];
sx q[1];
rz(0.8144905) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096734418) q[0];
sx q[0];
rz(-2.0196901) q[0];
sx q[0];
rz(-0.4273703) q[0];
rz(-pi) q[1];
rz(-2.3454104) q[2];
sx q[2];
rz(-2.0416235) q[2];
sx q[2];
rz(-1.4126965) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8100639) q[1];
sx q[1];
rz(-1.074407) q[1];
sx q[1];
rz(-1.0277102) q[1];
rz(-pi) q[2];
rz(-0.67072224) q[3];
sx q[3];
rz(-0.94778396) q[3];
sx q[3];
rz(-2.5581806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3689975) q[2];
sx q[2];
rz(-1.7721704) q[2];
sx q[2];
rz(2.7102846) q[2];
rz(-1.6469693) q[3];
sx q[3];
rz(-1.5458115) q[3];
sx q[3];
rz(-2.3561884) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5539219) q[0];
sx q[0];
rz(-2.9449154) q[0];
sx q[0];
rz(1.3171296) q[0];
rz(-1.0493086) q[1];
sx q[1];
rz(-2.6056555) q[1];
sx q[1];
rz(2.5228693) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65370377) q[0];
sx q[0];
rz(-1.426322) q[0];
sx q[0];
rz(0.61019759) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.87845408) q[2];
sx q[2];
rz(-1.8900423) q[2];
sx q[2];
rz(-0.4934267) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2545027) q[1];
sx q[1];
rz(-1.9190333) q[1];
sx q[1];
rz(0.13910267) q[1];
rz(-0.57024184) q[3];
sx q[3];
rz(-2.5240353) q[3];
sx q[3];
rz(2.9887426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.72545663) q[2];
sx q[2];
rz(-1.8014896) q[2];
sx q[2];
rz(-0.49187342) q[2];
rz(-1.5623931) q[3];
sx q[3];
rz(-1.648936) q[3];
sx q[3];
rz(2.5621342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1270776) q[0];
sx q[0];
rz(-2.0836232) q[0];
sx q[0];
rz(-0.28011093) q[0];
rz(0.07490553) q[1];
sx q[1];
rz(-0.43498755) q[1];
sx q[1];
rz(-1.7710955) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6579191) q[0];
sx q[0];
rz(-2.3020805) q[0];
sx q[0];
rz(-0.28261225) q[0];
x q[1];
rz(-2.2417775) q[2];
sx q[2];
rz(-1.0428937) q[2];
sx q[2];
rz(1.4769276) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6447923) q[1];
sx q[1];
rz(-2.519517) q[1];
sx q[1];
rz(0.46228564) q[1];
x q[2];
rz(3.0093569) q[3];
sx q[3];
rz(-1.5894798) q[3];
sx q[3];
rz(1.736369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8568153) q[2];
sx q[2];
rz(-1.6723526) q[2];
sx q[2];
rz(2.3397297) q[2];
rz(2.2297468) q[3];
sx q[3];
rz(-2.3175479) q[3];
sx q[3];
rz(-2.1626933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3367679) q[0];
sx q[0];
rz(-2.473859) q[0];
sx q[0];
rz(-2.5208852) q[0];
rz(-0.2785109) q[1];
sx q[1];
rz(-1.6997489) q[1];
sx q[1];
rz(-1.0034358) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0517563) q[0];
sx q[0];
rz(-0.47430719) q[0];
sx q[0];
rz(0.82292212) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.070095991) q[2];
sx q[2];
rz(-1.150944) q[2];
sx q[2];
rz(-1.1481089) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.089503376) q[1];
sx q[1];
rz(-1.516253) q[1];
sx q[1];
rz(2.8724345) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5430319) q[3];
sx q[3];
rz(-1.7669456) q[3];
sx q[3];
rz(-2.83441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8371381) q[2];
sx q[2];
rz(-1.1123603) q[2];
sx q[2];
rz(0.066085286) q[2];
rz(-2.0153913) q[3];
sx q[3];
rz(-2.3812713) q[3];
sx q[3];
rz(-0.57536212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5758301) q[0];
sx q[0];
rz(-3.059721) q[0];
sx q[0];
rz(0.78731147) q[0];
rz(-0.85834223) q[1];
sx q[1];
rz(-0.97802496) q[1];
sx q[1];
rz(1.0816921) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4183732) q[0];
sx q[0];
rz(-1.0745418) q[0];
sx q[0];
rz(-2.4346274) q[0];
rz(-pi) q[1];
rz(0.85140793) q[2];
sx q[2];
rz(-0.59307304) q[2];
sx q[2];
rz(-2.3874758) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3216721) q[1];
sx q[1];
rz(-1.5770438) q[1];
sx q[1];
rz(0.28467463) q[1];
rz(-2.2899707) q[3];
sx q[3];
rz(-2.0083665) q[3];
sx q[3];
rz(-0.17862602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.16296884) q[2];
sx q[2];
rz(-1.1498412) q[2];
sx q[2];
rz(-2.2990189) q[2];
rz(0.27635559) q[3];
sx q[3];
rz(-2.1601951) q[3];
sx q[3];
rz(-1.3493376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2862947) q[0];
sx q[0];
rz(-1.8305625) q[0];
sx q[0];
rz(1.1899765) q[0];
rz(-1.2225993) q[1];
sx q[1];
rz(-1.2346376) q[1];
sx q[1];
rz(-2.55866) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5258497) q[0];
sx q[0];
rz(-1.2229563) q[0];
sx q[0];
rz(-0.49973947) q[0];
x q[1];
rz(0.013986258) q[2];
sx q[2];
rz(-0.99726935) q[2];
sx q[2];
rz(-1.397246) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.035810245) q[1];
sx q[1];
rz(-1.8086193) q[1];
sx q[1];
rz(1.6538804) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.510079) q[3];
sx q[3];
rz(-1.2119635) q[3];
sx q[3];
rz(-1.3541019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0588093) q[2];
sx q[2];
rz(-1.2502111) q[2];
sx q[2];
rz(-2.1018551) q[2];
rz(-2.5522363) q[3];
sx q[3];
rz(-1.4732692) q[3];
sx q[3];
rz(1.0065494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7650918) q[0];
sx q[0];
rz(-2.9917175) q[0];
sx q[0];
rz(-1.3366706) q[0];
rz(2.916015) q[1];
sx q[1];
rz(-1.0712737) q[1];
sx q[1];
rz(-2.371686) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8827646) q[0];
sx q[0];
rz(-2.6443091) q[0];
sx q[0];
rz(0.33035853) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3138397) q[2];
sx q[2];
rz(-1.079664) q[2];
sx q[2];
rz(0.14358768) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0437255) q[1];
sx q[1];
rz(-1.7683523) q[1];
sx q[1];
rz(1.6476589) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0287241) q[3];
sx q[3];
rz(-2.5333571) q[3];
sx q[3];
rz(-0.0054159482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5169107) q[2];
sx q[2];
rz(-1.3930438) q[2];
sx q[2];
rz(1.0153655) q[2];
rz(2.9376302) q[3];
sx q[3];
rz(-1.5488307) q[3];
sx q[3];
rz(1.3913733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90149752) q[0];
sx q[0];
rz(-2.5101341) q[0];
sx q[0];
rz(0.39328662) q[0];
rz(-2.6914864) q[1];
sx q[1];
rz(-1.7186807) q[1];
sx q[1];
rz(0.030287655) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5092376) q[0];
sx q[0];
rz(-1.3233174) q[0];
sx q[0];
rz(-1.6028274) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3322163) q[2];
sx q[2];
rz(-1.7580108) q[2];
sx q[2];
rz(0.16448122) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.62822484) q[1];
sx q[1];
rz(-1.0252044) q[1];
sx q[1];
rz(1.5978769) q[1];
x q[2];
rz(-2.6859517) q[3];
sx q[3];
rz(-0.29248387) q[3];
sx q[3];
rz(2.4258421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0344737) q[2];
sx q[2];
rz(-1.5849042) q[2];
sx q[2];
rz(2.0195473) q[2];
rz(1.0551039) q[3];
sx q[3];
rz(-0.3873581) q[3];
sx q[3];
rz(-1.8900227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17290641) q[0];
sx q[0];
rz(-2.1538669) q[0];
sx q[0];
rz(-0.77447844) q[0];
rz(-2.5075746) q[1];
sx q[1];
rz(-2.1114383) q[1];
sx q[1];
rz(1.0672306) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33444302) q[0];
sx q[0];
rz(-0.73842885) q[0];
sx q[0];
rz(0.025511857) q[0];
rz(3.0926782) q[2];
sx q[2];
rz(-1.8777784) q[2];
sx q[2];
rz(3.0119257) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6744819) q[1];
sx q[1];
rz(-1.8201105) q[1];
sx q[1];
rz(-1.8915908) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9036071) q[3];
sx q[3];
rz(-2.6389671) q[3];
sx q[3];
rz(-2.7183333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.60712236) q[2];
sx q[2];
rz(-0.12132135) q[2];
sx q[2];
rz(1.8758476) q[2];
rz(0.040146116) q[3];
sx q[3];
rz(-0.37467343) q[3];
sx q[3];
rz(-3.1025187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(0.48113111) q[0];
sx q[0];
rz(-2.2534695) q[0];
sx q[0];
rz(-1.2598502) q[0];
rz(-1.4866359) q[1];
sx q[1];
rz(-2.1075552) q[1];
sx q[1];
rz(-3.1390417) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77230763) q[0];
sx q[0];
rz(-0.17371674) q[0];
sx q[0];
rz(0.80887167) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4223406) q[2];
sx q[2];
rz(-0.064777834) q[2];
sx q[2];
rz(1.5671519) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8067183) q[1];
sx q[1];
rz(-2.2069227) q[1];
sx q[1];
rz(-1.0735372) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0052913) q[3];
sx q[3];
rz(-1.5948203) q[3];
sx q[3];
rz(2.0498118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.5551787) q[2];
sx q[2];
rz(-1.6224344) q[2];
sx q[2];
rz(2.4671538) q[2];
rz(-1.1174348) q[3];
sx q[3];
rz(-2.1148465) q[3];
sx q[3];
rz(2.8286772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.085070327) q[0];
sx q[0];
rz(-1.1310348) q[0];
sx q[0];
rz(-0.52666589) q[0];
rz(-2.7936735) q[1];
sx q[1];
rz(-1.9031453) q[1];
sx q[1];
rz(-1.3487945) q[1];
rz(-1.3580218) q[2];
sx q[2];
rz(-1.909303) q[2];
sx q[2];
rz(2.252966) q[2];
rz(-0.017853768) q[3];
sx q[3];
rz(-1.4244867) q[3];
sx q[3];
rz(-1.3467237) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
