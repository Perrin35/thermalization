OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6991601) q[0];
sx q[0];
rz(4.5259024) q[0];
sx q[0];
rz(10.685267) q[0];
rz(-1.0386382) q[1];
sx q[1];
rz(4.4903978) q[1];
sx q[1];
rz(8.5010565) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7605654) q[0];
sx q[0];
rz(-1.1228704) q[0];
sx q[0];
rz(-2.7678124) q[0];
rz(2.997424) q[2];
sx q[2];
rz(-1.8509794) q[2];
sx q[2];
rz(-0.35137128) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0337861) q[1];
sx q[1];
rz(-0.68192712) q[1];
sx q[1];
rz(2.4297907) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5039623) q[3];
sx q[3];
rz(-0.59083592) q[3];
sx q[3];
rz(-1.6319815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2279921) q[2];
sx q[2];
rz(-1.2486518) q[2];
sx q[2];
rz(-2.9795734) q[2];
rz(2.2062733) q[3];
sx q[3];
rz(-0.98615065) q[3];
sx q[3];
rz(-0.71301618) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0682003) q[0];
sx q[0];
rz(-0.22664264) q[0];
sx q[0];
rz(1.1967999) q[0];
rz(0.67990047) q[1];
sx q[1];
rz(-2.6459243) q[1];
sx q[1];
rz(-1.4555567) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50549492) q[0];
sx q[0];
rz(-2.539145) q[0];
sx q[0];
rz(1.2732182) q[0];
x q[1];
rz(-2.935264) q[2];
sx q[2];
rz(-0.38197877) q[2];
sx q[2];
rz(-2.1260335) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.204139) q[1];
sx q[1];
rz(-2.0935645) q[1];
sx q[1];
rz(0.015603113) q[1];
rz(-pi) q[2];
rz(-3.0807642) q[3];
sx q[3];
rz(-2.7292477) q[3];
sx q[3];
rz(2.6044248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7130647) q[2];
sx q[2];
rz(-1.45168) q[2];
sx q[2];
rz(1.7896174) q[2];
rz(2.9591566) q[3];
sx q[3];
rz(-2.1648516) q[3];
sx q[3];
rz(-0.3119719) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75333726) q[0];
sx q[0];
rz(-0.68080807) q[0];
sx q[0];
rz(-0.80048168) q[0];
rz(0.02877409) q[1];
sx q[1];
rz(-2.0859699) q[1];
sx q[1];
rz(1.172539) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4370255) q[0];
sx q[0];
rz(-0.66172681) q[0];
sx q[0];
rz(-0.5163124) q[0];
x q[1];
rz(-2.795479) q[2];
sx q[2];
rz(-2.3556404) q[2];
sx q[2];
rz(1.9217938) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.064627083) q[1];
sx q[1];
rz(-0.58276999) q[1];
sx q[1];
rz(2.0240192) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3200687) q[3];
sx q[3];
rz(-1.6276976) q[3];
sx q[3];
rz(-1.0292605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0744434) q[2];
sx q[2];
rz(-1.6643486) q[2];
sx q[2];
rz(0.91119901) q[2];
rz(0.95101142) q[3];
sx q[3];
rz(-0.8042897) q[3];
sx q[3];
rz(2.2495911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7610385) q[0];
sx q[0];
rz(-3.0111713) q[0];
sx q[0];
rz(-0.12810853) q[0];
rz(3.065486) q[1];
sx q[1];
rz(-1.9271306) q[1];
sx q[1];
rz(0.52350837) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9317755) q[0];
sx q[0];
rz(-2.1486001) q[0];
sx q[0];
rz(-2.0156167) q[0];
rz(-pi) q[1];
rz(1.4857616) q[2];
sx q[2];
rz(-1.8997314) q[2];
sx q[2];
rz(2.0563682) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.47141155) q[1];
sx q[1];
rz(-1.3172611) q[1];
sx q[1];
rz(-2.2815435) q[1];
x q[2];
rz(-0.13682271) q[3];
sx q[3];
rz(-2.1577583) q[3];
sx q[3];
rz(2.6329991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.52544242) q[2];
sx q[2];
rz(-1.5972861) q[2];
sx q[2];
rz(0.564044) q[2];
rz(2.8530252) q[3];
sx q[3];
rz(-2.7189062) q[3];
sx q[3];
rz(-0.55571663) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6600835) q[0];
sx q[0];
rz(-0.68843377) q[0];
sx q[0];
rz(1.6500641) q[0];
rz(0.87961698) q[1];
sx q[1];
rz(-1.8938226) q[1];
sx q[1];
rz(-2.1496444) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6447727) q[0];
sx q[0];
rz(-0.20475514) q[0];
sx q[0];
rz(-0.55069189) q[0];
rz(-0.0089346272) q[2];
sx q[2];
rz(-1.2679456) q[2];
sx q[2];
rz(-0.47700275) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.79210287) q[1];
sx q[1];
rz(-1.8794685) q[1];
sx q[1];
rz(-1.5285138) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.81327849) q[3];
sx q[3];
rz(-1.6653898) q[3];
sx q[3];
rz(0.41155848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1296967) q[2];
sx q[2];
rz(-2.7719438) q[2];
sx q[2];
rz(-0.27080718) q[2];
rz(-2.9233542) q[3];
sx q[3];
rz(-1.3202347) q[3];
sx q[3];
rz(2.9158084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-0.72702423) q[0];
sx q[0];
rz(-0.71838656) q[0];
sx q[0];
rz(-1.3487934) q[0];
rz(2.7596966) q[1];
sx q[1];
rz(-2.8254639) q[1];
sx q[1];
rz(1.4250925) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8558559) q[0];
sx q[0];
rz(-0.55186134) q[0];
sx q[0];
rz(-0.43666552) q[0];
x q[1];
rz(1.233333) q[2];
sx q[2];
rz(-0.26974264) q[2];
sx q[2];
rz(2.595682) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.52275601) q[1];
sx q[1];
rz(-0.4546051) q[1];
sx q[1];
rz(1.418581) q[1];
rz(-0.53374966) q[3];
sx q[3];
rz(-1.0373877) q[3];
sx q[3];
rz(0.48983869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1340593) q[2];
sx q[2];
rz(-2.7441661) q[2];
sx q[2];
rz(0.56387222) q[2];
rz(2.9610736) q[3];
sx q[3];
rz(-1.518395) q[3];
sx q[3];
rz(0.40294161) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6329704) q[0];
sx q[0];
rz(-2.9794725) q[0];
sx q[0];
rz(0.41931835) q[0];
rz(-1.58889) q[1];
sx q[1];
rz(-1.2607375) q[1];
sx q[1];
rz(0.82180506) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.069698378) q[0];
sx q[0];
rz(-2.0928203) q[0];
sx q[0];
rz(-0.72899039) q[0];
x q[1];
rz(0.49634883) q[2];
sx q[2];
rz(-2.2349572) q[2];
sx q[2];
rz(-1.6830483) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2591178) q[1];
sx q[1];
rz(-1.7977409) q[1];
sx q[1];
rz(2.4005753) q[1];
rz(-pi) q[2];
rz(1.9415226) q[3];
sx q[3];
rz(-1.0154187) q[3];
sx q[3];
rz(-2.6851418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8043148) q[2];
sx q[2];
rz(-2.3874805) q[2];
sx q[2];
rz(-0.24469963) q[2];
rz(-0.129536) q[3];
sx q[3];
rz(-1.1641538) q[3];
sx q[3];
rz(-1.6285508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.725175) q[0];
sx q[0];
rz(-3.1224407) q[0];
sx q[0];
rz(-2.3186671) q[0];
rz(-0.30934632) q[1];
sx q[1];
rz(-1.7495218) q[1];
sx q[1];
rz(1.3051422) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7628521) q[0];
sx q[0];
rz(-1.2438602) q[0];
sx q[0];
rz(-2.5255894) q[0];
rz(-pi) q[1];
rz(1.9403946) q[2];
sx q[2];
rz(-0.9005138) q[2];
sx q[2];
rz(3.0073462) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9040363) q[1];
sx q[1];
rz(-2.190553) q[1];
sx q[1];
rz(2.0380286) q[1];
rz(-3.1180624) q[3];
sx q[3];
rz(-1.2315893) q[3];
sx q[3];
rz(-2.6587405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0017073) q[2];
sx q[2];
rz(-1.7557764) q[2];
sx q[2];
rz(1.4902327) q[2];
rz(2.0643318) q[3];
sx q[3];
rz(-2.1765985) q[3];
sx q[3];
rz(0.13154496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(1.7548783) q[0];
sx q[0];
rz(-1.2972378) q[0];
sx q[0];
rz(-2.819678) q[0];
rz(-1.5362668) q[1];
sx q[1];
rz(-1.9202817) q[1];
sx q[1];
rz(-2.4386491) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87899938) q[0];
sx q[0];
rz(-1.9290553) q[0];
sx q[0];
rz(-0.4155638) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4998352) q[2];
sx q[2];
rz(-2.6312201) q[2];
sx q[2];
rz(-0.77529782) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.798588) q[1];
sx q[1];
rz(-1.4061905) q[1];
sx q[1];
rz(-2.1144457) q[1];
rz(-2.1861107) q[3];
sx q[3];
rz(-1.2947047) q[3];
sx q[3];
rz(-2.626782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2408509) q[2];
sx q[2];
rz(-0.27975953) q[2];
sx q[2];
rz(-1.3396324) q[2];
rz(2.83589) q[3];
sx q[3];
rz(-1.327508) q[3];
sx q[3];
rz(-1.3302749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9777578) q[0];
sx q[0];
rz(-0.75755388) q[0];
sx q[0];
rz(-1.2257858) q[0];
rz(2.2380791) q[1];
sx q[1];
rz(-0.61360306) q[1];
sx q[1];
rz(2.6729565) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9194473) q[0];
sx q[0];
rz(-1.3627909) q[0];
sx q[0];
rz(-2.749445) q[0];
rz(-pi) q[1];
rz(0.32613812) q[2];
sx q[2];
rz(-1.9378127) q[2];
sx q[2];
rz(-2.0740167) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.84506449) q[1];
sx q[1];
rz(-1.6390641) q[1];
sx q[1];
rz(1.8370085) q[1];
x q[2];
rz(-2.0300794) q[3];
sx q[3];
rz(-1.1772403) q[3];
sx q[3];
rz(-0.72898385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.48352155) q[2];
sx q[2];
rz(-1.2926241) q[2];
sx q[2];
rz(1.998385) q[2];
rz(3.0269567) q[3];
sx q[3];
rz(-0.95364037) q[3];
sx q[3];
rz(1.6121929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4951915) q[0];
sx q[0];
rz(-1.9468745) q[0];
sx q[0];
rz(-0.68328802) q[0];
rz(-0.62190965) q[1];
sx q[1];
rz(-1.6786631) q[1];
sx q[1];
rz(2.8181029) q[1];
rz(2.3399578) q[2];
sx q[2];
rz(-2.2710706) q[2];
sx q[2];
rz(2.0588277) q[2];
rz(2.2223496) q[3];
sx q[3];
rz(-1.5083434) q[3];
sx q[3];
rz(-1.2283243) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
