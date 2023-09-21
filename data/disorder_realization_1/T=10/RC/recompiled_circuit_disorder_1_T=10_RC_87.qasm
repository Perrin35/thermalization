OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3172265) q[0];
sx q[0];
rz(-2.0269725) q[0];
sx q[0];
rz(-3.1414519) q[0];
rz(-1.8074942) q[1];
sx q[1];
rz(-0.9642095) q[1];
sx q[1];
rz(1.948184) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6407335) q[0];
sx q[0];
rz(-1.8574323) q[0];
sx q[0];
rz(-0.1851693) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5932181) q[2];
sx q[2];
rz(-1.3142685) q[2];
sx q[2];
rz(2.246849) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.27543435) q[1];
sx q[1];
rz(-0.94238102) q[1];
sx q[1];
rz(-2.1584312) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0317806) q[3];
sx q[3];
rz(-1.3545274) q[3];
sx q[3];
rz(-3.0307378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6821735) q[2];
sx q[2];
rz(-3.1176304) q[2];
sx q[2];
rz(-1.2288644) q[2];
rz(1.7284283) q[3];
sx q[3];
rz(-2.0404405) q[3];
sx q[3];
rz(1.4878954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6035778) q[0];
sx q[0];
rz(-1.6390272) q[0];
sx q[0];
rz(-2.1287825) q[0];
rz(3.1139328) q[1];
sx q[1];
rz(-2.467997) q[1];
sx q[1];
rz(-1.123463) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7682122) q[0];
sx q[0];
rz(-2.0313615) q[0];
sx q[0];
rz(0.064231355) q[0];
rz(-pi) q[1];
x q[1];
rz(0.77983071) q[2];
sx q[2];
rz(-1.5260328) q[2];
sx q[2];
rz(-2.0335399) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0151129) q[1];
sx q[1];
rz(-1.0547332) q[1];
sx q[1];
rz(0.59241398) q[1];
rz(-pi) q[2];
rz(0.65621891) q[3];
sx q[3];
rz(-1.1074293) q[3];
sx q[3];
rz(0.046063395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3479487) q[2];
sx q[2];
rz(-1.0898033) q[2];
sx q[2];
rz(2.222555) q[2];
rz(2.4675026) q[3];
sx q[3];
rz(-2.489311) q[3];
sx q[3];
rz(-1.526171) q[3];
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
rz(0.27750257) q[0];
sx q[0];
rz(-2.9798177) q[0];
sx q[0];
rz(1.2751689) q[0];
rz(0.69349849) q[1];
sx q[1];
rz(-1.8854515) q[1];
sx q[1];
rz(2.0085874) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5223761) q[0];
sx q[0];
rz(-1.4285061) q[0];
sx q[0];
rz(3.1382568) q[0];
x q[1];
rz(1.0622382) q[2];
sx q[2];
rz(-0.8478176) q[2];
sx q[2];
rz(-1.6194956) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9085711) q[1];
sx q[1];
rz(-1.034522) q[1];
sx q[1];
rz(-0.33645333) q[1];
rz(-0.77550019) q[3];
sx q[3];
rz(-2.3470078) q[3];
sx q[3];
rz(-0.19898181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.2514078) q[2];
sx q[2];
rz(-2.3501985) q[2];
sx q[2];
rz(-1.2934925) q[2];
rz(0.039316468) q[3];
sx q[3];
rz(-1.2189564) q[3];
sx q[3];
rz(1.2600651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.2599729) q[0];
sx q[0];
rz(-3.0631174) q[0];
sx q[0];
rz(-1.9807293) q[0];
rz(-0.89598957) q[1];
sx q[1];
rz(-1.7005824) q[1];
sx q[1];
rz(0.13555759) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82942688) q[0];
sx q[0];
rz(-0.78071785) q[0];
sx q[0];
rz(-0.45129816) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9549471) q[2];
sx q[2];
rz(-1.640056) q[2];
sx q[2];
rz(0.70089507) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.902066) q[1];
sx q[1];
rz(-0.24589989) q[1];
sx q[1];
rz(-0.86218254) q[1];
rz(-2.0849864) q[3];
sx q[3];
rz(-0.31294542) q[3];
sx q[3];
rz(1.3514951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9049412) q[2];
sx q[2];
rz(-2.1950978) q[2];
sx q[2];
rz(2.2616852) q[2];
rz(3.0974292) q[3];
sx q[3];
rz(-1.6396089) q[3];
sx q[3];
rz(-2.8529609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1039466) q[0];
sx q[0];
rz(-0.3750616) q[0];
sx q[0];
rz(-2.1283545) q[0];
rz(3.0918616) q[1];
sx q[1];
rz(-2.2278992) q[1];
sx q[1];
rz(-2.0577046) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5729382) q[0];
sx q[0];
rz(-1.3875811) q[0];
sx q[0];
rz(1.818548) q[0];
rz(-pi) q[1];
rz(-2.3124144) q[2];
sx q[2];
rz(-0.36311705) q[2];
sx q[2];
rz(-1.6104289) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0005972) q[1];
sx q[1];
rz(-2.6037569) q[1];
sx q[1];
rz(1.3586033) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0524592) q[3];
sx q[3];
rz(-0.51557487) q[3];
sx q[3];
rz(0.46407947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9087387) q[2];
sx q[2];
rz(-0.32662699) q[2];
sx q[2];
rz(-2.8971635) q[2];
rz(-2.7092253) q[3];
sx q[3];
rz(-1.3997388) q[3];
sx q[3];
rz(-0.50306815) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2844834) q[0];
sx q[0];
rz(-1.720022) q[0];
sx q[0];
rz(0.094141468) q[0];
rz(-0.17177467) q[1];
sx q[1];
rz(-1.1356907) q[1];
sx q[1];
rz(2.24618) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58127922) q[0];
sx q[0];
rz(-0.34438294) q[0];
sx q[0];
rz(3.0292105) q[0];
rz(-pi) q[1];
rz(0.32161153) q[2];
sx q[2];
rz(-2.4827637) q[2];
sx q[2];
rz(-2.4668601) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6678644) q[1];
sx q[1];
rz(-1.6147991) q[1];
sx q[1];
rz(1.4021137) q[1];
rz(-pi) q[2];
x q[2];
rz(1.496109) q[3];
sx q[3];
rz(-1.5684621) q[3];
sx q[3];
rz(0.96935779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0078997) q[2];
sx q[2];
rz(-2.7329625) q[2];
sx q[2];
rz(2.3383979) q[2];
rz(-1.9512272) q[3];
sx q[3];
rz(-1.2322216) q[3];
sx q[3];
rz(2.7289594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.068709277) q[0];
sx q[0];
rz(-2.9769653) q[0];
sx q[0];
rz(0.51914006) q[0];
rz(-2.5601162) q[1];
sx q[1];
rz(-2.0362208) q[1];
sx q[1];
rz(1.2566459) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.885868) q[0];
sx q[0];
rz(-0.2547383) q[0];
sx q[0];
rz(-1.7314918) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9314552) q[2];
sx q[2];
rz(-1.1147095) q[2];
sx q[2];
rz(-1.2880585) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7278442) q[1];
sx q[1];
rz(-1.8520978) q[1];
sx q[1];
rz(1.5037699) q[1];
rz(0.36862936) q[3];
sx q[3];
rz(-0.7437403) q[3];
sx q[3];
rz(-1.7891235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.98823035) q[2];
sx q[2];
rz(-2.1116657) q[2];
sx q[2];
rz(-1.3640277) q[2];
rz(2.2310232) q[3];
sx q[3];
rz(-1.1547337) q[3];
sx q[3];
rz(1.5301269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0664739) q[0];
sx q[0];
rz(-0.56448889) q[0];
sx q[0];
rz(-2.8334154) q[0];
rz(-3.0691052) q[1];
sx q[1];
rz(-2.1283573) q[1];
sx q[1];
rz(-2.7546308) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68650866) q[0];
sx q[0];
rz(-0.63502705) q[0];
sx q[0];
rz(-1.1839266) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2690291) q[2];
sx q[2];
rz(-1.1485032) q[2];
sx q[2];
rz(-2.1625105) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4657198) q[1];
sx q[1];
rz(-2.4657001) q[1];
sx q[1];
rz(3.0216316) q[1];
rz(-pi) q[2];
rz(-0.49658637) q[3];
sx q[3];
rz(-1.7644617) q[3];
sx q[3];
rz(2.1878939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5179634) q[2];
sx q[2];
rz(-1.7771746) q[2];
sx q[2];
rz(0.44000885) q[2];
rz(0.7157588) q[3];
sx q[3];
rz(-1.7093168) q[3];
sx q[3];
rz(-1.0796775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0004262) q[0];
sx q[0];
rz(-0.74581242) q[0];
sx q[0];
rz(-2.0429042) q[0];
rz(-0.72775841) q[1];
sx q[1];
rz(-2.7658503) q[1];
sx q[1];
rz(3.0922906) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5987199) q[0];
sx q[0];
rz(-2.1763986) q[0];
sx q[0];
rz(2.4332895) q[0];
rz(-pi) q[1];
rz(2.6463037) q[2];
sx q[2];
rz(-0.9501179) q[2];
sx q[2];
rz(2.4923313) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2652492) q[1];
sx q[1];
rz(-2.0298404) q[1];
sx q[1];
rz(-0.52258073) q[1];
rz(-1.0827761) q[3];
sx q[3];
rz(-2.2501695) q[3];
sx q[3];
rz(2.8454012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.07842841) q[2];
sx q[2];
rz(-2.5421263) q[2];
sx q[2];
rz(2.4196529) q[2];
rz(-2.1980964) q[3];
sx q[3];
rz(-0.75073457) q[3];
sx q[3];
rz(0.25434428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(0.14770517) q[0];
sx q[0];
rz(-1.1575971) q[0];
sx q[0];
rz(-1.0797427) q[0];
rz(2.0823157) q[1];
sx q[1];
rz(-2.9187027) q[1];
sx q[1];
rz(-1.4019029) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3831543) q[0];
sx q[0];
rz(-0.24807319) q[0];
sx q[0];
rz(1.0323348) q[0];
rz(-2.9548168) q[2];
sx q[2];
rz(-1.597607) q[2];
sx q[2];
rz(-1.1550127) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8733752) q[1];
sx q[1];
rz(-2.2443319) q[1];
sx q[1];
rz(0.39336494) q[1];
x q[2];
rz(1.7781156) q[3];
sx q[3];
rz(-0.94982409) q[3];
sx q[3];
rz(0.45351115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6587276) q[2];
sx q[2];
rz(-1.3663224) q[2];
sx q[2];
rz(-1.520291) q[2];
rz(0.55082095) q[3];
sx q[3];
rz(-2.3362624) q[3];
sx q[3];
rz(0.6974535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14810066) q[0];
sx q[0];
rz(-1.8363331) q[0];
sx q[0];
rz(1.6114417) q[0];
rz(0.91611721) q[1];
sx q[1];
rz(-0.59090186) q[1];
sx q[1];
rz(-0.59060243) q[1];
rz(2.3475636) q[2];
sx q[2];
rz(-0.90894884) q[2];
sx q[2];
rz(-0.85469645) q[2];
rz(1.0536853) q[3];
sx q[3];
rz(-1.5552945) q[3];
sx q[3];
rz(-1.8154715) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
