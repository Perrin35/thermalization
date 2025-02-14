OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7810516) q[0];
sx q[0];
rz(3.7563503) q[0];
sx q[0];
rz(10.496245) q[0];
rz(-0.40845025) q[1];
sx q[1];
rz(4.0790494) q[1];
sx q[1];
rz(11.117878) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9928701) q[0];
sx q[0];
rz(-2.5186484) q[0];
sx q[0];
rz(0.8533661) q[0];
rz(-pi) q[1];
rz(1.229847) q[2];
sx q[2];
rz(-2.1181137) q[2];
sx q[2];
rz(2.7173619) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.422157) q[1];
sx q[1];
rz(-1.7195175) q[1];
sx q[1];
rz(1.4805111) q[1];
rz(-0.86032773) q[3];
sx q[3];
rz(-1.9091879) q[3];
sx q[3];
rz(2.8380054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.94413269) q[2];
sx q[2];
rz(-2.2977273) q[2];
sx q[2];
rz(2.3485363) q[2];
rz(-0.26886764) q[3];
sx q[3];
rz(-1.3258508) q[3];
sx q[3];
rz(2.1813006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31084138) q[0];
sx q[0];
rz(-2.7590867) q[0];
sx q[0];
rz(-0.88019669) q[0];
rz(2.510732) q[1];
sx q[1];
rz(-2.3289101) q[1];
sx q[1];
rz(-1.0989443) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99115935) q[0];
sx q[0];
rz(-2.9637103) q[0];
sx q[0];
rz(-3.1129897) q[0];
rz(-2.1889792) q[2];
sx q[2];
rz(-1.435155) q[2];
sx q[2];
rz(-0.73276765) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6827755) q[1];
sx q[1];
rz(-1.4084245) q[1];
sx q[1];
rz(-2.008325) q[1];
rz(-pi) q[2];
rz(0.97084011) q[3];
sx q[3];
rz(-0.73093855) q[3];
sx q[3];
rz(0.39492861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.38689303) q[2];
sx q[2];
rz(-1.6563481) q[2];
sx q[2];
rz(-0.54163843) q[2];
rz(0.70332876) q[3];
sx q[3];
rz(-0.43916217) q[3];
sx q[3];
rz(2.7642803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.90848732) q[0];
sx q[0];
rz(-2.2623514) q[0];
sx q[0];
rz(-0.096906699) q[0];
rz(-2.5540409) q[1];
sx q[1];
rz(-0.89043003) q[1];
sx q[1];
rz(2.5403835) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1977229) q[0];
sx q[0];
rz(-1.3900779) q[0];
sx q[0];
rz(0.69036071) q[0];
rz(-1.634609) q[2];
sx q[2];
rz(-1.4087965) q[2];
sx q[2];
rz(2.8500003) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.280698) q[1];
sx q[1];
rz(-1.1105683) q[1];
sx q[1];
rz(1.5621508) q[1];
rz(2.2279068) q[3];
sx q[3];
rz(-1.2877712) q[3];
sx q[3];
rz(-2.117827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2931557) q[2];
sx q[2];
rz(-1.6907254) q[2];
sx q[2];
rz(-2.4270774) q[2];
rz(1.6589288) q[3];
sx q[3];
rz(-2.8667993) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7191384) q[0];
sx q[0];
rz(-1.9062573) q[0];
sx q[0];
rz(-1.6374913) q[0];
rz(-2.0488886) q[1];
sx q[1];
rz(-1.8605109) q[1];
sx q[1];
rz(-0.5853931) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4883746) q[0];
sx q[0];
rz(-1.8420353) q[0];
sx q[0];
rz(-2.491889) q[0];
rz(-pi) q[1];
x q[1];
rz(0.010527654) q[2];
sx q[2];
rz(-1.4665571) q[2];
sx q[2];
rz(0.57360211) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8429381) q[1];
sx q[1];
rz(-2.2982592) q[1];
sx q[1];
rz(-0.91489376) q[1];
x q[2];
rz(1.0799418) q[3];
sx q[3];
rz(-2.1532791) q[3];
sx q[3];
rz(-2.7459308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.66712159) q[2];
sx q[2];
rz(-2.4690364) q[2];
sx q[2];
rz(-2.7453864) q[2];
rz(2.951156) q[3];
sx q[3];
rz(-2.2021144) q[3];
sx q[3];
rz(-2.6207391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-0.025295479) q[0];
sx q[0];
rz(-0.46023661) q[0];
sx q[0];
rz(-2.7468371) q[0];
rz(2.1341628) q[1];
sx q[1];
rz(-1.1578683) q[1];
sx q[1];
rz(-1.6068858) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10478131) q[0];
sx q[0];
rz(-0.80950538) q[0];
sx q[0];
rz(-0.39489371) q[0];
x q[1];
rz(-1.2053648) q[2];
sx q[2];
rz(-2.2534568) q[2];
sx q[2];
rz(-0.78675573) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6794656) q[1];
sx q[1];
rz(-1.7662342) q[1];
sx q[1];
rz(1.2567169) q[1];
rz(-pi) q[2];
rz(-2.0232361) q[3];
sx q[3];
rz(-1.9944189) q[3];
sx q[3];
rz(2.771559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4452303) q[2];
sx q[2];
rz(-2.4424489) q[2];
sx q[2];
rz(2.9635079) q[2];
rz(2.1961424) q[3];
sx q[3];
rz(-0.33792308) q[3];
sx q[3];
rz(-0.43615714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6984542) q[0];
sx q[0];
rz(-2.4789424) q[0];
sx q[0];
rz(-0.15289256) q[0];
rz(-2.1235509) q[1];
sx q[1];
rz(-2.1986304) q[1];
sx q[1];
rz(-2.5315888) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7383282) q[0];
sx q[0];
rz(-0.77304196) q[0];
sx q[0];
rz(-2.9570262) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.583408) q[2];
sx q[2];
rz(-0.57390814) q[2];
sx q[2];
rz(-1.877026) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.568012) q[1];
sx q[1];
rz(-1.3369155) q[1];
sx q[1];
rz(2.3877445) q[1];
x q[2];
rz(0.69645564) q[3];
sx q[3];
rz(-0.12731931) q[3];
sx q[3];
rz(-2.6501973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4623922) q[2];
sx q[2];
rz(-1.2827164) q[2];
sx q[2];
rz(0.25897762) q[2];
rz(0.14608832) q[3];
sx q[3];
rz(-0.77272213) q[3];
sx q[3];
rz(2.5025388) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4884969) q[0];
sx q[0];
rz(-2.275832) q[0];
sx q[0];
rz(-0.75188941) q[0];
rz(2.409626) q[1];
sx q[1];
rz(-1.3804133) q[1];
sx q[1];
rz(-2.6123349) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.142685) q[0];
sx q[0];
rz(-1.7325263) q[0];
sx q[0];
rz(-1.2611752) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2970778) q[2];
sx q[2];
rz(-0.47011312) q[2];
sx q[2];
rz(2.672247) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.25594357) q[1];
sx q[1];
rz(-0.082271345) q[1];
sx q[1];
rz(1.72797) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1629754) q[3];
sx q[3];
rz(-0.53272034) q[3];
sx q[3];
rz(0.5513834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.25972128) q[2];
sx q[2];
rz(-2.112381) q[2];
sx q[2];
rz(0.79505801) q[2];
rz(-1.9184387) q[3];
sx q[3];
rz(-1.3431679) q[3];
sx q[3];
rz(-0.77643001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.935598) q[0];
sx q[0];
rz(-1.5482276) q[0];
sx q[0];
rz(-3.026631) q[0];
rz(3.0523172) q[1];
sx q[1];
rz(-2.142579) q[1];
sx q[1];
rz(2.9866536) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45231014) q[0];
sx q[0];
rz(-1.4703106) q[0];
sx q[0];
rz(-1.7013676) q[0];
x q[1];
rz(-0.55122113) q[2];
sx q[2];
rz(-2.1586426) q[2];
sx q[2];
rz(-1.8648704) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.66471264) q[1];
sx q[1];
rz(-2.6694593) q[1];
sx q[1];
rz(0.5916729) q[1];
x q[2];
rz(-0.84568923) q[3];
sx q[3];
rz(-1.5165303) q[3];
sx q[3];
rz(2.3783663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0838919) q[2];
sx q[2];
rz(-1.8996779) q[2];
sx q[2];
rz(0.6883626) q[2];
rz(-2.5734731) q[3];
sx q[3];
rz(-1.0024242) q[3];
sx q[3];
rz(-0.68305558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6883009) q[0];
sx q[0];
rz(-2.330307) q[0];
sx q[0];
rz(1.2055093) q[0];
rz(-2.0506809) q[1];
sx q[1];
rz(-0.58735192) q[1];
sx q[1];
rz(-1.6039414) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4805877) q[0];
sx q[0];
rz(-1.93297) q[0];
sx q[0];
rz(1.1042547) q[0];
rz(-pi) q[1];
rz(-3.1411016) q[2];
sx q[2];
rz(-2.1386991) q[2];
sx q[2];
rz(0.22746023) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4548757) q[1];
sx q[1];
rz(-1.6819961) q[1];
sx q[1];
rz(0.79373549) q[1];
x q[2];
rz(2.5439369) q[3];
sx q[3];
rz(-1.9697876) q[3];
sx q[3];
rz(-1.9982709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9566112) q[2];
sx q[2];
rz(-1.1143755) q[2];
sx q[2];
rz(1.1448917) q[2];
rz(-1.4780809) q[3];
sx q[3];
rz(-2.3749115) q[3];
sx q[3];
rz(-1.112282) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0612653) q[0];
sx q[0];
rz(-2.5054131) q[0];
sx q[0];
rz(1.7058477) q[0];
rz(-2.0044633) q[1];
sx q[1];
rz(-0.69157332) q[1];
sx q[1];
rz(-0.93708509) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0264056) q[0];
sx q[0];
rz(-1.6583024) q[0];
sx q[0];
rz(2.8390272) q[0];
rz(-pi) q[1];
rz(0.50855277) q[2];
sx q[2];
rz(-1.1932696) q[2];
sx q[2];
rz(0.14103954) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.62793362) q[1];
sx q[1];
rz(-0.80027831) q[1];
sx q[1];
rz(-1.5350128) q[1];
x q[2];
rz(2.6804018) q[3];
sx q[3];
rz(-0.66890016) q[3];
sx q[3];
rz(-1.8708285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.81637853) q[2];
sx q[2];
rz(-0.64322317) q[2];
sx q[2];
rz(-2.8312259) q[2];
rz(2.5988633) q[3];
sx q[3];
rz(-0.92971814) q[3];
sx q[3];
rz(2.824596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9515789) q[0];
sx q[0];
rz(-2.0484476) q[0];
sx q[0];
rz(0.53566326) q[0];
rz(-0.71470064) q[1];
sx q[1];
rz(-1.8938046) q[1];
sx q[1];
rz(-1.6751777) q[1];
rz(1.9541478) q[2];
sx q[2];
rz(-1.162983) q[2];
sx q[2];
rz(2.4833124) q[2];
rz(0.72890239) q[3];
sx q[3];
rz(-1.5255675) q[3];
sx q[3];
rz(-0.062011157) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
