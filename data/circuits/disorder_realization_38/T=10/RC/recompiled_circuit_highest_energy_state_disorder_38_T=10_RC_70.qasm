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
rz(2.8528557) q[0];
sx q[0];
rz(-0.69536916) q[0];
sx q[0];
rz(-2.8737972) q[0];
rz(0.42203045) q[1];
sx q[1];
rz(4.0623436) q[1];
sx q[1];
rz(10.698591) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86912495) q[0];
sx q[0];
rz(-2.8960138) q[0];
sx q[0];
rz(1.7564943) q[0];
x q[1];
rz(-2.535475) q[2];
sx q[2];
rz(-0.96783468) q[2];
sx q[2];
rz(0.22113344) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.63034025) q[1];
sx q[1];
rz(-0.66424886) q[1];
sx q[1];
rz(-0.11099191) q[1];
rz(-1.4273604) q[3];
sx q[3];
rz(-1.535245) q[3];
sx q[3];
rz(-1.2811023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6196809) q[2];
sx q[2];
rz(-0.64138594) q[2];
sx q[2];
rz(3.0989975) q[2];
rz(0.28111449) q[3];
sx q[3];
rz(-1.576985) q[3];
sx q[3];
rz(-0.59578305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5113145) q[0];
sx q[0];
rz(-1.1742641) q[0];
sx q[0];
rz(2.0654772) q[0];
rz(1.0307182) q[1];
sx q[1];
rz(-1.1117671) q[1];
sx q[1];
rz(1.8310742) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56633184) q[0];
sx q[0];
rz(-1.5916414) q[0];
sx q[0];
rz(2.2290609) q[0];
x q[1];
rz(1.1300226) q[2];
sx q[2];
rz(-2.3557969) q[2];
sx q[2];
rz(-2.1919029) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4736997) q[1];
sx q[1];
rz(-1.7381428) q[1];
sx q[1];
rz(1.6909864) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1635246) q[3];
sx q[3];
rz(-2.0680475) q[3];
sx q[3];
rz(-0.45491342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.92404667) q[2];
sx q[2];
rz(-2.3895538) q[2];
sx q[2];
rz(-2.1853866) q[2];
rz(-1.3077959) q[3];
sx q[3];
rz(-2.3266413) q[3];
sx q[3];
rz(3.0202878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8420551) q[0];
sx q[0];
rz(-2.6326023) q[0];
sx q[0];
rz(3.1053542) q[0];
rz(-0.80129519) q[1];
sx q[1];
rz(-1.5943297) q[1];
sx q[1];
rz(1.3005728) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1195923) q[0];
sx q[0];
rz(-1.4109857) q[0];
sx q[0];
rz(-0.10800604) q[0];
rz(-pi) q[1];
rz(-1.5074129) q[2];
sx q[2];
rz(-2.5434539) q[2];
sx q[2];
rz(-1.6063521) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8138995) q[1];
sx q[1];
rz(-2.3997953) q[1];
sx q[1];
rz(-2.2375536) q[1];
x q[2];
rz(0.68978975) q[3];
sx q[3];
rz(-1.8588603) q[3];
sx q[3];
rz(-0.40460872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7654045) q[2];
sx q[2];
rz(-1.3724962) q[2];
sx q[2];
rz(-0.21793951) q[2];
rz(-0.91207063) q[3];
sx q[3];
rz(-1.6488766) q[3];
sx q[3];
rz(0.85165858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4172149) q[0];
sx q[0];
rz(-2.0700924) q[0];
sx q[0];
rz(2.3260314) q[0];
rz(-1.0527481) q[1];
sx q[1];
rz(-2.2595854) q[1];
sx q[1];
rz(1.3515333) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.049168436) q[0];
sx q[0];
rz(-1.6480849) q[0];
sx q[0];
rz(-2.671675) q[0];
x q[1];
rz(-0.10064023) q[2];
sx q[2];
rz(-2.422524) q[2];
sx q[2];
rz(-1.6883862) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.4554418) q[1];
sx q[1];
rz(-2.0459922) q[1];
sx q[1];
rz(0.032217044) q[1];
rz(-pi) q[2];
rz(2.5791956) q[3];
sx q[3];
rz(-0.80005433) q[3];
sx q[3];
rz(2.7871015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1286596) q[2];
sx q[2];
rz(-1.1933051) q[2];
sx q[2];
rz(2.7093757) q[2];
rz(2.2991119) q[3];
sx q[3];
rz(-2.1079) q[3];
sx q[3];
rz(-0.99561083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1603482) q[0];
sx q[0];
rz(-0.46537414) q[0];
sx q[0];
rz(2.5904742) q[0];
rz(-0.61800686) q[1];
sx q[1];
rz(-1.2851241) q[1];
sx q[1];
rz(2.0726223) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38902125) q[0];
sx q[0];
rz(-2.4271936) q[0];
sx q[0];
rz(-2.8117489) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0238012) q[2];
sx q[2];
rz(-2.0108622) q[2];
sx q[2];
rz(-2.5534782) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3772419) q[1];
sx q[1];
rz(-0.83989401) q[1];
sx q[1];
rz(2.9725171) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.51286814) q[3];
sx q[3];
rz(-1.697043) q[3];
sx q[3];
rz(3.0579289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3507639) q[2];
sx q[2];
rz(-2.5574234) q[2];
sx q[2];
rz(0.9551777) q[2];
rz(2.5480934) q[3];
sx q[3];
rz(-2.1953526) q[3];
sx q[3];
rz(1.745863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3747568) q[0];
sx q[0];
rz(-3.0992442) q[0];
sx q[0];
rz(-1.391885) q[0];
rz(-1.1426686) q[1];
sx q[1];
rz(-1.7997768) q[1];
sx q[1];
rz(-1.8291738) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29762938) q[0];
sx q[0];
rz(-2.2052551) q[0];
sx q[0];
rz(-2.5353801) q[0];
x q[1];
rz(-0.67709558) q[2];
sx q[2];
rz(-0.99092084) q[2];
sx q[2];
rz(0.092965417) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.23199546) q[1];
sx q[1];
rz(-1.0568406) q[1];
sx q[1];
rz(-1.6327052) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.083701) q[3];
sx q[3];
rz(-1.8857919) q[3];
sx q[3];
rz(-2.0712899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7890847) q[2];
sx q[2];
rz(-2.3859873) q[2];
sx q[2];
rz(1.4136723) q[2];
rz(-0.46755725) q[3];
sx q[3];
rz(-2.4831725) q[3];
sx q[3];
rz(2.5889034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6308052) q[0];
sx q[0];
rz(-3.0785705) q[0];
sx q[0];
rz(-0.68156534) q[0];
rz(1.956578) q[1];
sx q[1];
rz(-0.76528913) q[1];
sx q[1];
rz(-2.3550745) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2241572) q[0];
sx q[0];
rz(-2.6154714) q[0];
sx q[0];
rz(-0.71147646) q[0];
x q[1];
rz(2.6441387) q[2];
sx q[2];
rz(-1.0062394) q[2];
sx q[2];
rz(1.4184784) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6326846) q[1];
sx q[1];
rz(-0.68636319) q[1];
sx q[1];
rz(-1.4014161) q[1];
rz(1.5249355) q[3];
sx q[3];
rz(-2.4099382) q[3];
sx q[3];
rz(1.589244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.53175348) q[2];
sx q[2];
rz(-0.24007758) q[2];
sx q[2];
rz(1.5698203) q[2];
rz(1.2872559) q[3];
sx q[3];
rz(-1.5232892) q[3];
sx q[3];
rz(-0.94304812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.542273) q[0];
sx q[0];
rz(-0.71745187) q[0];
sx q[0];
rz(1.9532816) q[0];
rz(0.020847281) q[1];
sx q[1];
rz(-2.3884845) q[1];
sx q[1];
rz(-2.5426224) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5064829) q[0];
sx q[0];
rz(-2.0744893) q[0];
sx q[0];
rz(-3.0290108) q[0];
rz(-pi) q[1];
rz(-2.7487603) q[2];
sx q[2];
rz(-2.246703) q[2];
sx q[2];
rz(2.3623938) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.61781949) q[1];
sx q[1];
rz(-2.0460772) q[1];
sx q[1];
rz(0.75264024) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3009572) q[3];
sx q[3];
rz(-1.908769) q[3];
sx q[3];
rz(2.4455051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8591566) q[2];
sx q[2];
rz(-1.891529) q[2];
sx q[2];
rz(0.51521987) q[2];
rz(-0.53269261) q[3];
sx q[3];
rz(-1.0849413) q[3];
sx q[3];
rz(2.2044619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1135947) q[0];
sx q[0];
rz(-2.4585215) q[0];
sx q[0];
rz(0.37044507) q[0];
rz(-2.0619552) q[1];
sx q[1];
rz(-2.7307983) q[1];
sx q[1];
rz(3.0830141) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1616042) q[0];
sx q[0];
rz(-1.4128886) q[0];
sx q[0];
rz(-0.61451332) q[0];
rz(-pi) q[1];
x q[1];
rz(0.04444261) q[2];
sx q[2];
rz(-1.7979597) q[2];
sx q[2];
rz(-2.4054766) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1778922) q[1];
sx q[1];
rz(-0.78989115) q[1];
sx q[1];
rz(0.13479418) q[1];
rz(-pi) q[2];
rz(-2.305916) q[3];
sx q[3];
rz(-2.0234851) q[3];
sx q[3];
rz(1.9122461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2728682) q[2];
sx q[2];
rz(-2.3385907) q[2];
sx q[2];
rz(-0.33099428) q[2];
rz(3.1281779) q[3];
sx q[3];
rz(-0.9534854) q[3];
sx q[3];
rz(1.0564055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57824221) q[0];
sx q[0];
rz(-2.712482) q[0];
sx q[0];
rz(-0.44813928) q[0];
rz(-3.0478802) q[1];
sx q[1];
rz(-0.30148503) q[1];
sx q[1];
rz(-1.9261446) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96998065) q[0];
sx q[0];
rz(-1.6002065) q[0];
sx q[0];
rz(2.0038811) q[0];
rz(0.81268572) q[2];
sx q[2];
rz(-1.865109) q[2];
sx q[2];
rz(3.091623) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3642973) q[1];
sx q[1];
rz(-2.5357995) q[1];
sx q[1];
rz(2.0531027) q[1];
rz(-pi) q[2];
rz(-1.2566725) q[3];
sx q[3];
rz(-2.7249955) q[3];
sx q[3];
rz(-2.9521717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.965968) q[2];
sx q[2];
rz(-1.5886687) q[2];
sx q[2];
rz(2.442181) q[2];
rz(-1.8661963) q[3];
sx q[3];
rz(-1.1558775) q[3];
sx q[3];
rz(1.8705961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4324343) q[0];
sx q[0];
rz(-0.86031886) q[0];
sx q[0];
rz(1.1501089) q[0];
rz(1.746183) q[1];
sx q[1];
rz(-1.2184873) q[1];
sx q[1];
rz(1.7582735) q[1];
rz(-0.10192656) q[2];
sx q[2];
rz(-0.66833767) q[2];
sx q[2];
rz(2.3705033) q[2];
rz(0.46635177) q[3];
sx q[3];
rz(-1.4644571) q[3];
sx q[3];
rz(-1.1639948) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
