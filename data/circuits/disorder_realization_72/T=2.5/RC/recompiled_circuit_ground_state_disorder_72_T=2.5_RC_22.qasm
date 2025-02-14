OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2312343) q[0];
sx q[0];
rz(5.4141747) q[0];
sx q[0];
rz(10.5095) q[0];
rz(1.9864858) q[1];
sx q[1];
rz(-2.3218563) q[1];
sx q[1];
rz(-2.2302332) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40075743) q[0];
sx q[0];
rz(-2.3199953) q[0];
sx q[0];
rz(1.2375968) q[0];
rz(-pi) q[1];
rz(-1.4248104) q[2];
sx q[2];
rz(-2.2772191) q[2];
sx q[2];
rz(1.8391158) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6877796) q[1];
sx q[1];
rz(-2.0168843) q[1];
sx q[1];
rz(0.58961726) q[1];
rz(0.064685589) q[3];
sx q[3];
rz(-2.0035335) q[3];
sx q[3];
rz(2.329934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.34065166) q[2];
sx q[2];
rz(-2.030535) q[2];
sx q[2];
rz(-0.079455376) q[2];
rz(2.5331412) q[3];
sx q[3];
rz(-1.2234917) q[3];
sx q[3];
rz(-0.89103812) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4873753) q[0];
sx q[0];
rz(-0.86367622) q[0];
sx q[0];
rz(-1.7864216) q[0];
rz(-1.0379418) q[1];
sx q[1];
rz(-0.67960056) q[1];
sx q[1];
rz(0.80744809) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45189598) q[0];
sx q[0];
rz(-2.1679255) q[0];
sx q[0];
rz(0.35291617) q[0];
rz(-0.37952559) q[2];
sx q[2];
rz(-2.1132601) q[2];
sx q[2];
rz(1.2778953) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3088919) q[1];
sx q[1];
rz(-1.8101793) q[1];
sx q[1];
rz(-1.6390087) q[1];
x q[2];
rz(3.03504) q[3];
sx q[3];
rz(-1.4378743) q[3];
sx q[3];
rz(-0.35655856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.21343931) q[2];
sx q[2];
rz(-1.5779053) q[2];
sx q[2];
rz(1.6820924) q[2];
rz(-0.45808074) q[3];
sx q[3];
rz(-0.57772485) q[3];
sx q[3];
rz(-0.95988449) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9816575) q[0];
sx q[0];
rz(-1.0502879) q[0];
sx q[0];
rz(2.4666069) q[0];
rz(-0.58810294) q[1];
sx q[1];
rz(-2.2876078) q[1];
sx q[1];
rz(1.3311707) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0494306) q[0];
sx q[0];
rz(-1.1750918) q[0];
sx q[0];
rz(2.6108598) q[0];
rz(0.66548062) q[2];
sx q[2];
rz(-1.1534216) q[2];
sx q[2];
rz(2.5145384) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9362671) q[1];
sx q[1];
rz(-0.85204879) q[1];
sx q[1];
rz(-0.34013744) q[1];
rz(-pi) q[2];
rz(-1.210545) q[3];
sx q[3];
rz(-0.52740806) q[3];
sx q[3];
rz(-2.1338192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5171234) q[2];
sx q[2];
rz(-0.34999592) q[2];
sx q[2];
rz(2.9208753) q[2];
rz(-2.3366426) q[3];
sx q[3];
rz(-1.5419818) q[3];
sx q[3];
rz(1.8055003) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6014366) q[0];
sx q[0];
rz(-2.8995081) q[0];
sx q[0];
rz(-2.5228187) q[0];
rz(-3.0586808) q[1];
sx q[1];
rz(-0.34268788) q[1];
sx q[1];
rz(1.6212911) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0205958) q[0];
sx q[0];
rz(-1.474664) q[0];
sx q[0];
rz(-1.0445879) q[0];
rz(1.3341212) q[2];
sx q[2];
rz(-2.1370721) q[2];
sx q[2];
rz(-2.2353954) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.027923) q[1];
sx q[1];
rz(-1.4884559) q[1];
sx q[1];
rz(-1.0587949) q[1];
rz(-pi) q[2];
rz(1.3906562) q[3];
sx q[3];
rz(-1.9400175) q[3];
sx q[3];
rz(-0.1843017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3380022) q[2];
sx q[2];
rz(-3.0323196) q[2];
sx q[2];
rz(0.14536157) q[2];
rz(0.27211443) q[3];
sx q[3];
rz(-1.4067255) q[3];
sx q[3];
rz(1.9947778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.014932545) q[0];
sx q[0];
rz(-1.8155875) q[0];
sx q[0];
rz(-0.10619157) q[0];
rz(-0.95942489) q[1];
sx q[1];
rz(-0.54490772) q[1];
sx q[1];
rz(0.19439654) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6694227) q[0];
sx q[0];
rz(-1.6568156) q[0];
sx q[0];
rz(0.73366139) q[0];
rz(-pi) q[1];
rz(-1.6856367) q[2];
sx q[2];
rz(-1.5756338) q[2];
sx q[2];
rz(2.9034535) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0433309) q[1];
sx q[1];
rz(-2.1594467) q[1];
sx q[1];
rz(-1.5215988) q[1];
x q[2];
rz(0.077094519) q[3];
sx q[3];
rz(-1.1655679) q[3];
sx q[3];
rz(-0.73540686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.83547366) q[2];
sx q[2];
rz(-1.4543616) q[2];
sx q[2];
rz(0.56841889) q[2];
rz(3.0889555) q[3];
sx q[3];
rz(-1.6391552) q[3];
sx q[3];
rz(-3.0047825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1133872) q[0];
sx q[0];
rz(-0.55332342) q[0];
sx q[0];
rz(0.18948874) q[0];
rz(2.1151309) q[1];
sx q[1];
rz(-2.2920513) q[1];
sx q[1];
rz(2.386327) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7758691) q[0];
sx q[0];
rz(-2.2138174) q[0];
sx q[0];
rz(0.68277208) q[0];
rz(-0.83024518) q[2];
sx q[2];
rz(-0.70607042) q[2];
sx q[2];
rz(2.490807) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2456296) q[1];
sx q[1];
rz(-1.4019483) q[1];
sx q[1];
rz(0.16014512) q[1];
rz(-1.2047661) q[3];
sx q[3];
rz(-1.2377973) q[3];
sx q[3];
rz(-2.7878417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.75817627) q[2];
sx q[2];
rz(-1.5551609) q[2];
sx q[2];
rz(-1.4413393) q[2];
rz(-1.7656743) q[3];
sx q[3];
rz(-2.223189) q[3];
sx q[3];
rz(-2.4413696) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43948424) q[0];
sx q[0];
rz(-2.0937884) q[0];
sx q[0];
rz(1.1442319) q[0];
rz(-0.2598091) q[1];
sx q[1];
rz(-0.83994284) q[1];
sx q[1];
rz(-2.1536486) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.086083273) q[0];
sx q[0];
rz(-1.6702154) q[0];
sx q[0];
rz(1.5539597) q[0];
rz(-1.9960711) q[2];
sx q[2];
rz(-2.2605611) q[2];
sx q[2];
rz(-0.20343101) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.14079796) q[1];
sx q[1];
rz(-0.70604815) q[1];
sx q[1];
rz(-2.7435859) q[1];
x q[2];
rz(-0.99616237) q[3];
sx q[3];
rz(-1.8271433) q[3];
sx q[3];
rz(-2.4202895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.40954956) q[2];
sx q[2];
rz(-1.6403551) q[2];
sx q[2];
rz(0.6130971) q[2];
rz(-2.4260855) q[3];
sx q[3];
rz(-0.77176538) q[3];
sx q[3];
rz(-0.36090052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20329423) q[0];
sx q[0];
rz(-0.84380117) q[0];
sx q[0];
rz(-0.26813689) q[0];
rz(0.2306436) q[1];
sx q[1];
rz(-1.5430887) q[1];
sx q[1];
rz(-0.014009744) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62416158) q[0];
sx q[0];
rz(-0.78963477) q[0];
sx q[0];
rz(2.8454418) q[0];
rz(-1.2410937) q[2];
sx q[2];
rz(-2.6791141) q[2];
sx q[2];
rz(-1.7623189) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.18302984) q[1];
sx q[1];
rz(-1.6154449) q[1];
sx q[1];
rz(2.6041998) q[1];
rz(-1.7746065) q[3];
sx q[3];
rz(-2.9958354) q[3];
sx q[3];
rz(-2.3043119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9739428) q[2];
sx q[2];
rz(-1.8029982) q[2];
sx q[2];
rz(0.30926427) q[2];
rz(0.15657982) q[3];
sx q[3];
rz(-1.7370109) q[3];
sx q[3];
rz(-2.1046765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2579047) q[0];
sx q[0];
rz(-0.64249277) q[0];
sx q[0];
rz(0.25074348) q[0];
rz(-0.58654395) q[1];
sx q[1];
rz(-1.5778912) q[1];
sx q[1];
rz(2.3768545) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4737807) q[0];
sx q[0];
rz(-2.7062253) q[0];
sx q[0];
rz(-2.3171168) q[0];
x q[1];
rz(0.7711507) q[2];
sx q[2];
rz(-1.7131107) q[2];
sx q[2];
rz(-0.76616312) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2231317) q[1];
sx q[1];
rz(-1.4479965) q[1];
sx q[1];
rz(1.5801013) q[1];
x q[2];
rz(-1.5161726) q[3];
sx q[3];
rz(-1.9491073) q[3];
sx q[3];
rz(-2.9448933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5137198) q[2];
sx q[2];
rz(-0.79078117) q[2];
sx q[2];
rz(-0.17390832) q[2];
rz(2.3993313) q[3];
sx q[3];
rz(-2.3895013) q[3];
sx q[3];
rz(0.043005634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1327508) q[0];
sx q[0];
rz(-0.77021563) q[0];
sx q[0];
rz(2.8736864) q[0];
rz(-1.335089) q[1];
sx q[1];
rz(-0.91064149) q[1];
sx q[1];
rz(1.7693899) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89248449) q[0];
sx q[0];
rz(-1.3297446) q[0];
sx q[0];
rz(-2.0349099) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6716953) q[2];
sx q[2];
rz(-1.865662) q[2];
sx q[2];
rz(0.036525846) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6698996) q[1];
sx q[1];
rz(-1.0581543) q[1];
sx q[1];
rz(-2.1218845) q[1];
x q[2];
rz(0.80268152) q[3];
sx q[3];
rz(-0.44375989) q[3];
sx q[3];
rz(2.8783523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.09769663) q[2];
sx q[2];
rz(-1.0363657) q[2];
sx q[2];
rz(-1.7608661) q[2];
rz(1.1287639) q[3];
sx q[3];
rz(-2.1916316) q[3];
sx q[3];
rz(1.5855764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2632521) q[0];
sx q[0];
rz(-1.5338407) q[0];
sx q[0];
rz(-1.1080678) q[0];
rz(-2.4551328) q[1];
sx q[1];
rz(-1.7484799) q[1];
sx q[1];
rz(1.9304986) q[1];
rz(2.6348719) q[2];
sx q[2];
rz(-2.5891149) q[2];
sx q[2];
rz(2.0074809) q[2];
rz(-1.0330647) q[3];
sx q[3];
rz(-1.2119537) q[3];
sx q[3];
rz(-2.4338943) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
