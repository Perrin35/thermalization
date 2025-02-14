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
rz(-0.82395616) q[0];
sx q[0];
rz(-0.92744654) q[0];
rz(-4.1915512) q[1];
sx q[1];
rz(-0.97133049) q[1];
sx q[1];
rz(8.6102875) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.096734418) q[0];
sx q[0];
rz(-2.0196901) q[0];
sx q[0];
rz(-2.7142224) q[0];
x q[1];
rz(0.61887069) q[2];
sx q[2];
rz(-0.89779389) q[2];
sx q[2];
rz(0.25970632) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8100639) q[1];
sx q[1];
rz(-1.074407) q[1];
sx q[1];
rz(-2.1138825) q[1];
rz(-pi) q[2];
rz(-0.82858927) q[3];
sx q[3];
rz(-2.0999206) q[3];
sx q[3];
rz(0.55381004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.77259511) q[2];
sx q[2];
rz(-1.3694222) q[2];
sx q[2];
rz(0.43130809) q[2];
rz(-1.4946233) q[3];
sx q[3];
rz(-1.5458115) q[3];
sx q[3];
rz(2.3561884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5876708) q[0];
sx q[0];
rz(-0.19667721) q[0];
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
rz(2.4275682) q[0];
sx q[0];
rz(-2.51665) q[0];
sx q[0];
rz(-0.24863909) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0486062) q[2];
sx q[2];
rz(-2.3903762) q[2];
sx q[2];
rz(-0.71556811) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.26855384) q[1];
sx q[1];
rz(-1.4400926) q[1];
sx q[1];
rz(-1.9221582) q[1];
rz(-pi) q[2];
rz(-1.9369164) q[3];
sx q[3];
rz(-2.0799326) q[3];
sx q[3];
rz(0.81936554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.72545663) q[2];
sx q[2];
rz(-1.8014896) q[2];
sx q[2];
rz(-0.49187342) q[2];
rz(1.5623931) q[3];
sx q[3];
rz(-1.648936) q[3];
sx q[3];
rz(-2.5621342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1270776) q[0];
sx q[0];
rz(-1.0579695) q[0];
sx q[0];
rz(-0.28011093) q[0];
rz(0.07490553) q[1];
sx q[1];
rz(-2.7066051) q[1];
sx q[1];
rz(1.7710955) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48367354) q[0];
sx q[0];
rz(-0.83951211) q[0];
sx q[0];
rz(2.8589804) q[0];
rz(-pi) q[1];
rz(0.81746705) q[2];
sx q[2];
rz(-2.3139179) q[2];
sx q[2];
rz(0.65929123) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6447923) q[1];
sx q[1];
rz(-0.62207568) q[1];
sx q[1];
rz(0.46228564) q[1];
rz(-pi) q[2];
rz(0.13223572) q[3];
sx q[3];
rz(-1.5894798) q[3];
sx q[3];
rz(1.4052237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.28477731) q[2];
sx q[2];
rz(-1.4692401) q[2];
sx q[2];
rz(-0.80186296) q[2];
rz(0.91184584) q[3];
sx q[3];
rz(-0.82404476) q[3];
sx q[3];
rz(-2.1626933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80482471) q[0];
sx q[0];
rz(-0.66773361) q[0];
sx q[0];
rz(2.5208852) q[0];
rz(-0.2785109) q[1];
sx q[1];
rz(-1.6997489) q[1];
sx q[1];
rz(-1.0034358) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28351706) q[0];
sx q[0];
rz(-1.9122313) q[0];
sx q[0];
rz(0.33591875) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9915646) q[2];
sx q[2];
rz(-1.506797) q[2];
sx q[2];
rz(-2.6902932) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4651854) q[1];
sx q[1];
rz(-2.8670951) q[1];
sx q[1];
rz(0.20250116) q[1];
x q[2];
rz(-2.5430319) q[3];
sx q[3];
rz(-1.374647) q[3];
sx q[3];
rz(0.3071827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8371381) q[2];
sx q[2];
rz(-1.1123603) q[2];
sx q[2];
rz(0.066085286) q[2];
rz(1.1262013) q[3];
sx q[3];
rz(-0.76032138) q[3];
sx q[3];
rz(0.57536212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5758301) q[0];
sx q[0];
rz(-3.059721) q[0];
sx q[0];
rz(-2.3542812) q[0];
rz(-2.2832504) q[1];
sx q[1];
rz(-2.1635677) q[1];
sx q[1];
rz(-2.0599005) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72321945) q[0];
sx q[0];
rz(-1.0745418) q[0];
sx q[0];
rz(0.70696522) q[0];
rz(-0.85140793) q[2];
sx q[2];
rz(-2.5485196) q[2];
sx q[2];
rz(0.75411686) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3216721) q[1];
sx q[1];
rz(-1.5770438) q[1];
sx q[1];
rz(0.28467463) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.85162195) q[3];
sx q[3];
rz(-1.1332261) q[3];
sx q[3];
rz(2.9629666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.16296884) q[2];
sx q[2];
rz(-1.1498412) q[2];
sx q[2];
rz(2.2990189) q[2];
rz(0.27635559) q[3];
sx q[3];
rz(-0.98139757) q[3];
sx q[3];
rz(-1.792255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.855298) q[0];
sx q[0];
rz(-1.8305625) q[0];
sx q[0];
rz(1.1899765) q[0];
rz(-1.2225993) q[1];
sx q[1];
rz(-1.2346376) q[1];
sx q[1];
rz(-2.55866) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5258497) q[0];
sx q[0];
rz(-1.2229563) q[0];
sx q[0];
rz(2.6418532) q[0];
x q[1];
rz(-1.5924443) q[2];
sx q[2];
rz(-2.5679143) q[2];
sx q[2];
rz(1.4230185) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.035810245) q[1];
sx q[1];
rz(-1.3329734) q[1];
sx q[1];
rz(-1.4877122) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0058026) q[3];
sx q[3];
rz(-0.9851176) q[3];
sx q[3];
rz(-0.034736573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0827834) q[2];
sx q[2];
rz(-1.8913816) q[2];
sx q[2];
rz(-2.1018551) q[2];
rz(0.58935634) q[3];
sx q[3];
rz(-1.6683234) q[3];
sx q[3];
rz(-1.0065494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37650087) q[0];
sx q[0];
rz(-2.9917175) q[0];
sx q[0];
rz(1.804922) q[0];
rz(-2.916015) q[1];
sx q[1];
rz(-2.070319) q[1];
sx q[1];
rz(-2.371686) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25882803) q[0];
sx q[0];
rz(-2.6443091) q[0];
sx q[0];
rz(-2.8112341) q[0];
rz(2.3138397) q[2];
sx q[2];
rz(-2.0619287) q[2];
sx q[2];
rz(-2.998005) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0437255) q[1];
sx q[1];
rz(-1.7683523) q[1];
sx q[1];
rz(-1.6476589) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0287241) q[3];
sx q[3];
rz(-0.60823554) q[3];
sx q[3];
rz(0.0054159482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5169107) q[2];
sx q[2];
rz(-1.3930438) q[2];
sx q[2];
rz(-1.0153655) q[2];
rz(-2.9376302) q[3];
sx q[3];
rz(-1.5488307) q[3];
sx q[3];
rz(-1.3913733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90149752) q[0];
sx q[0];
rz(-0.63145852) q[0];
sx q[0];
rz(2.748306) q[0];
rz(0.45010629) q[1];
sx q[1];
rz(-1.7186807) q[1];
sx q[1];
rz(-3.111305) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0878828) q[0];
sx q[0];
rz(-1.5397414) q[0];
sx q[0];
rz(2.8939918) q[0];
x q[1];
rz(1.8387536) q[2];
sx q[2];
rz(-2.3620124) q[2];
sx q[2];
rz(-1.5991581) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.68037625) q[1];
sx q[1];
rz(-0.54619563) q[1];
sx q[1];
rz(-0.044574634) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4390596) q[3];
sx q[3];
rz(-1.8326958) q[3];
sx q[3];
rz(-1.1887664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1071189) q[2];
sx q[2];
rz(-1.5566885) q[2];
sx q[2];
rz(-2.0195473) q[2];
rz(-2.0864887) q[3];
sx q[3];
rz(-0.3873581) q[3];
sx q[3];
rz(-1.8900227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17290641) q[0];
sx q[0];
rz(-0.98772573) q[0];
sx q[0];
rz(2.3671142) q[0];
rz(-0.6340181) q[1];
sx q[1];
rz(-1.0301544) q[1];
sx q[1];
rz(-2.074362) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8071496) q[0];
sx q[0];
rz(-2.4031638) q[0];
sx q[0];
rz(-0.025511857) q[0];
x q[1];
rz(-1.2634694) q[2];
sx q[2];
rz(-1.6174223) q[2];
sx q[2];
rz(1.4263375) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.53505234) q[1];
sx q[1];
rz(-0.4036223) q[1];
sx q[1];
rz(0.89151793) q[1];
rz(-2.9638941) q[3];
sx q[3];
rz(-2.0434994) q[3];
sx q[3];
rz(-0.79897579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5344703) q[2];
sx q[2];
rz(-0.12132135) q[2];
sx q[2];
rz(-1.2657451) q[2];
rz(3.1014465) q[3];
sx q[3];
rz(-0.37467343) q[3];
sx q[3];
rz(-0.039073959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48113111) q[0];
sx q[0];
rz(-2.2534695) q[0];
sx q[0];
rz(-1.8817425) q[0];
rz(1.6549567) q[1];
sx q[1];
rz(-2.1075552) q[1];
sx q[1];
rz(0.0025509603) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1387864) q[0];
sx q[0];
rz(-1.696179) q[0];
sx q[0];
rz(0.12055293) q[0];
rz(-pi) q[1];
rz(-3.1319982) q[2];
sx q[2];
rz(-1.506732) q[2];
sx q[2];
rz(1.715915) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.59460179) q[1];
sx q[1];
rz(-2.3560215) q[1];
sx q[1];
rz(-2.5681096) q[1];
rz(-1.5137767) q[3];
sx q[3];
rz(-2.7064763) q[3];
sx q[3];
rz(-2.7142937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.5551787) q[2];
sx q[2];
rz(-1.6224344) q[2];
sx q[2];
rz(-2.4671538) q[2];
rz(-1.1174348) q[3];
sx q[3];
rz(-2.1148465) q[3];
sx q[3];
rz(2.8286772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.085070327) q[0];
sx q[0];
rz(-2.0105579) q[0];
sx q[0];
rz(2.6149268) q[0];
rz(2.7936735) q[1];
sx q[1];
rz(-1.2384474) q[1];
sx q[1];
rz(1.7927982) q[1];
rz(2.7958776) q[2];
sx q[2];
rz(-1.7713265) q[2];
sx q[2];
rz(-2.3878018) q[2];
rz(-1.4244637) q[3];
sx q[3];
rz(-1.5531333) q[3];
sx q[3];
rz(0.2266758) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
