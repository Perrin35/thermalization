OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.360541) q[0];
sx q[0];
rz(-0.6147576) q[0];
sx q[0];
rz(-1.0714666) q[0];
rz(2.7331424) q[1];
sx q[1];
rz(-0.93745679) q[1];
sx q[1];
rz(-1.6931005) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67249304) q[0];
sx q[0];
rz(-1.1156293) q[0];
sx q[0];
rz(-0.44123347) q[0];
rz(1.9117457) q[2];
sx q[2];
rz(-1.023479) q[2];
sx q[2];
rz(-0.42423074) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9706124) q[1];
sx q[1];
rz(-0.17380789) q[1];
sx q[1];
rz(-2.5998678) q[1];
rz(-2.0656282) q[3];
sx q[3];
rz(-0.77407661) q[3];
sx q[3];
rz(-2.2426393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.94413269) q[2];
sx q[2];
rz(-0.84386533) q[2];
sx q[2];
rz(-2.3485363) q[2];
rz(-0.26886764) q[3];
sx q[3];
rz(-1.8157418) q[3];
sx q[3];
rz(0.9602921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31084138) q[0];
sx q[0];
rz(-2.7590867) q[0];
sx q[0];
rz(0.88019669) q[0];
rz(-2.510732) q[1];
sx q[1];
rz(-2.3289101) q[1];
sx q[1];
rz(1.0989443) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1213721) q[0];
sx q[0];
rz(-1.3929875) q[0];
sx q[0];
rz(-1.5656548) q[0];
rz(-pi) q[1];
x q[1];
rz(0.95261344) q[2];
sx q[2];
rz(-1.7064377) q[2];
sx q[2];
rz(-2.408825) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4447282) q[1];
sx q[1];
rz(-2.676739) q[1];
sx q[1];
rz(1.2018728) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2078832) q[3];
sx q[3];
rz(-1.9572581) q[3];
sx q[3];
rz(-0.704788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7546996) q[2];
sx q[2];
rz(-1.4852445) q[2];
sx q[2];
rz(-2.5999542) q[2];
rz(-2.4382639) q[3];
sx q[3];
rz(-2.7024305) q[3];
sx q[3];
rz(0.37731236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2331053) q[0];
sx q[0];
rz(-2.2623514) q[0];
sx q[0];
rz(-3.044686) q[0];
rz(-2.5540409) q[1];
sx q[1];
rz(-2.2511626) q[1];
sx q[1];
rz(0.60120916) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5542463) q[0];
sx q[0];
rz(-2.4317435) q[0];
sx q[0];
rz(2.8621833) q[0];
rz(-pi) q[1];
x q[1];
rz(1.634609) q[2];
sx q[2];
rz(-1.7327961) q[2];
sx q[2];
rz(2.8500003) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.7060617) q[1];
sx q[1];
rz(-1.5630504) q[1];
sx q[1];
rz(0.46024291) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0151626) q[3];
sx q[3];
rz(-2.4345102) q[3];
sx q[3];
rz(-2.9418687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2931557) q[2];
sx q[2];
rz(-1.6907254) q[2];
sx q[2];
rz(-2.4270774) q[2];
rz(1.4826639) q[3];
sx q[3];
rz(-0.27479333) q[3];
sx q[3];
rz(2.7627435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42245427) q[0];
sx q[0];
rz(-1.2353354) q[0];
sx q[0];
rz(1.6374913) q[0];
rz(2.0488886) q[1];
sx q[1];
rz(-1.8605109) q[1];
sx q[1];
rz(0.5853931) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4883746) q[0];
sx q[0];
rz(-1.8420353) q[0];
sx q[0];
rz(-2.491889) q[0];
rz(0.010527654) q[2];
sx q[2];
rz(-1.4665571) q[2];
sx q[2];
rz(-2.5679905) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.55951872) q[1];
sx q[1];
rz(-2.2041956) q[1];
sx q[1];
rz(2.5410209) q[1];
rz(-2.5204757) q[3];
sx q[3];
rz(-0.74291544) q[3];
sx q[3];
rz(-2.7662504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.66712159) q[2];
sx q[2];
rz(-0.67255628) q[2];
sx q[2];
rz(2.7453864) q[2];
rz(-2.951156) q[3];
sx q[3];
rz(-0.93947828) q[3];
sx q[3];
rz(-2.6207391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.025295479) q[0];
sx q[0];
rz(-2.681356) q[0];
sx q[0];
rz(0.39475557) q[0];
rz(-2.1341628) q[1];
sx q[1];
rz(-1.9837244) q[1];
sx q[1];
rz(1.5347068) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7028977) q[0];
sx q[0];
rz(-0.83897018) q[0];
sx q[0];
rz(-1.1870866) q[0];
rz(-pi) q[1];
rz(1.2053648) q[2];
sx q[2];
rz(-2.2534568) q[2];
sx q[2];
rz(-2.3548369) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.462127) q[1];
sx q[1];
rz(-1.7662342) q[1];
sx q[1];
rz(1.2567169) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3716912) q[3];
sx q[3];
rz(-0.60958977) q[3];
sx q[3];
rz(-1.23884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.69636238) q[2];
sx q[2];
rz(-2.4424489) q[2];
sx q[2];
rz(-0.17808476) q[2];
rz(2.1961424) q[3];
sx q[3];
rz(-0.33792308) q[3];
sx q[3];
rz(2.7054355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6984542) q[0];
sx q[0];
rz(-2.4789424) q[0];
sx q[0];
rz(-0.15289256) q[0];
rz(-1.0180417) q[1];
sx q[1];
rz(-0.94296229) q[1];
sx q[1];
rz(0.61000383) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65839889) q[0];
sx q[0];
rz(-0.81419277) q[0];
sx q[0];
rz(-1.7479595) q[0];
rz(-2.1446682) q[2];
sx q[2];
rz(-1.5776433) q[2];
sx q[2];
rz(-2.8459542) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.57358068) q[1];
sx q[1];
rz(-1.3369155) q[1];
sx q[1];
rz(2.3877445) q[1];
rz(0.69645564) q[3];
sx q[3];
rz(-0.12731931) q[3];
sx q[3];
rz(-2.6501973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4623922) q[2];
sx q[2];
rz(-1.2827164) q[2];
sx q[2];
rz(-2.882615) q[2];
rz(0.14608832) q[3];
sx q[3];
rz(-0.77272213) q[3];
sx q[3];
rz(-0.63905382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6530957) q[0];
sx q[0];
rz(-0.86576068) q[0];
sx q[0];
rz(-2.3897032) q[0];
rz(2.409626) q[1];
sx q[1];
rz(-1.3804133) q[1];
sx q[1];
rz(-2.6123349) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2471591) q[0];
sx q[0];
rz(-0.34810796) q[0];
sx q[0];
rz(2.0624119) q[0];
rz(-pi) q[1];
rz(3.0050982) q[2];
sx q[2];
rz(-2.0220827) q[2];
sx q[2];
rz(-0.1642483) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0433474) q[1];
sx q[1];
rz(-1.4895413) q[1];
sx q[1];
rz(-3.1286865) q[1];
rz(-2.1629754) q[3];
sx q[3];
rz(-2.6088723) q[3];
sx q[3];
rz(-2.5902093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8818714) q[2];
sx q[2];
rz(-1.0292116) q[2];
sx q[2];
rz(0.79505801) q[2];
rz(-1.9184387) q[3];
sx q[3];
rz(-1.7984248) q[3];
sx q[3];
rz(-2.3651626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20599468) q[0];
sx q[0];
rz(-1.5482276) q[0];
sx q[0];
rz(-3.026631) q[0];
rz(-0.089275442) q[1];
sx q[1];
rz(-0.99901366) q[1];
sx q[1];
rz(0.1549391) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0099338) q[0];
sx q[0];
rz(-1.4408875) q[0];
sx q[0];
rz(-0.1013426) q[0];
x q[1];
rz(-0.9047382) q[2];
sx q[2];
rz(-0.78287941) q[2];
sx q[2];
rz(-0.43987405) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.66471264) q[1];
sx q[1];
rz(-0.47213337) q[1];
sx q[1];
rz(2.5499198) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0691419) q[3];
sx q[3];
rz(-0.84699291) q[3];
sx q[3];
rz(-0.75954306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0577008) q[2];
sx q[2];
rz(-1.2419147) q[2];
sx q[2];
rz(-0.6883626) q[2];
rz(2.5734731) q[3];
sx q[3];
rz(-1.0024242) q[3];
sx q[3];
rz(0.68305558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6883009) q[0];
sx q[0];
rz(-0.81128565) q[0];
sx q[0];
rz(-1.2055093) q[0];
rz(-2.0506809) q[1];
sx q[1];
rz(-0.58735192) q[1];
sx q[1];
rz(-1.6039414) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6610049) q[0];
sx q[0];
rz(-1.93297) q[0];
sx q[0];
rz(2.037338) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5715661) q[2];
sx q[2];
rz(-0.56790295) q[2];
sx q[2];
rz(0.22837328) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.22828211) q[1];
sx q[1];
rz(-0.78332114) q[1];
sx q[1];
rz(1.4128774) q[1];
rz(-pi) q[2];
rz(0.64301771) q[3];
sx q[3];
rz(-0.70482358) q[3];
sx q[3];
rz(0.091136668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9566112) q[2];
sx q[2];
rz(-1.1143755) q[2];
sx q[2];
rz(-1.1448917) q[2];
rz(-1.4780809) q[3];
sx q[3];
rz(-2.3749115) q[3];
sx q[3];
rz(2.0293106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0803273) q[0];
sx q[0];
rz(-0.63617951) q[0];
sx q[0];
rz(1.435745) q[0];
rz(2.0044633) q[1];
sx q[1];
rz(-0.69157332) q[1];
sx q[1];
rz(-2.2045076) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81737751) q[0];
sx q[0];
rz(-0.31459168) q[0];
sx q[0];
rz(0.28633519) q[0];
x q[1];
rz(-1.1446196) q[2];
sx q[2];
rz(-1.1010896) q[2];
sx q[2];
rz(1.6324279) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2236589) q[1];
sx q[1];
rz(-1.5964701) q[1];
sx q[1];
rz(-0.77083807) q[1];
rz(2.6804018) q[3];
sx q[3];
rz(-2.4726925) q[3];
sx q[3];
rz(-1.2707641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3252141) q[2];
sx q[2];
rz(-2.4983695) q[2];
sx q[2];
rz(-2.8312259) q[2];
rz(0.54272932) q[3];
sx q[3];
rz(-0.92971814) q[3];
sx q[3];
rz(-2.824596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19001374) q[0];
sx q[0];
rz(-1.0931451) q[0];
sx q[0];
rz(-2.6059294) q[0];
rz(0.71470064) q[1];
sx q[1];
rz(-1.2477881) q[1];
sx q[1];
rz(1.4664149) q[1];
rz(-0.71356365) q[2];
sx q[2];
rz(-0.55222558) q[2];
sx q[2];
rz(-1.4520558) q[2];
rz(3.0737446) q[3];
sx q[3];
rz(-2.4115457) q[3];
sx q[3];
rz(-1.5822165) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
