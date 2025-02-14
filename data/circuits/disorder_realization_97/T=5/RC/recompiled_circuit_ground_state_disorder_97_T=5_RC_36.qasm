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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67249304) q[0];
sx q[0];
rz(-2.0259633) q[0];
sx q[0];
rz(0.44123347) q[0];
rz(-2.5675814) q[2];
sx q[2];
rz(-1.8603627) q[2];
sx q[2];
rz(-2.1776108) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.422157) q[1];
sx q[1];
rz(-1.4220752) q[1];
sx q[1];
rz(1.6610816) q[1];
rz(-0.86032773) q[3];
sx q[3];
rz(-1.9091879) q[3];
sx q[3];
rz(2.8380054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.94413269) q[2];
sx q[2];
rz(-2.2977273) q[2];
sx q[2];
rz(-0.79305631) q[2];
rz(-0.26886764) q[3];
sx q[3];
rz(-1.3258508) q[3];
sx q[3];
rz(-0.9602921) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31084138) q[0];
sx q[0];
rz(-0.38250592) q[0];
sx q[0];
rz(-2.261396) q[0];
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
rz(2.1213721) q[0];
sx q[0];
rz(-1.7486052) q[0];
sx q[0];
rz(-1.5656548) q[0];
x q[1];
rz(-2.1889792) q[2];
sx q[2];
rz(-1.7064377) q[2];
sx q[2];
rz(-2.408825) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6968644) q[1];
sx q[1];
rz(-2.676739) q[1];
sx q[1];
rz(-1.2018728) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2078832) q[3];
sx q[3];
rz(-1.1843345) q[3];
sx q[3];
rz(-0.704788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7546996) q[2];
sx q[2];
rz(-1.4852445) q[2];
sx q[2];
rz(-2.5999542) q[2];
rz(-2.4382639) q[3];
sx q[3];
rz(-0.43916217) q[3];
sx q[3];
rz(2.7642803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2331053) q[0];
sx q[0];
rz(-2.2623514) q[0];
sx q[0];
rz(-3.044686) q[0];
rz(2.5540409) q[1];
sx q[1];
rz(-2.2511626) q[1];
sx q[1];
rz(-0.60120916) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9158949) q[0];
sx q[0];
rz(-2.2477799) q[0];
sx q[0];
rz(-1.8034794) q[0];
x q[1];
rz(1.5069836) q[2];
sx q[2];
rz(-1.7327961) q[2];
sx q[2];
rz(-2.8500003) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.86089462) q[1];
sx q[1];
rz(-2.0310244) q[1];
sx q[1];
rz(1.5794419) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2279068) q[3];
sx q[3];
rz(-1.2877712) q[3];
sx q[3];
rz(1.0237657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.84843695) q[2];
sx q[2];
rz(-1.6907254) q[2];
sx q[2];
rz(0.71451521) q[2];
rz(1.6589288) q[3];
sx q[3];
rz(-2.8667993) q[3];
sx q[3];
rz(2.7627435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42245427) q[0];
sx q[0];
rz(-1.9062573) q[0];
sx q[0];
rz(1.6374913) q[0];
rz(-2.0488886) q[1];
sx q[1];
rz(-1.8605109) q[1];
sx q[1];
rz(-0.5853931) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2566168) q[0];
sx q[0];
rz(-0.69640771) q[0];
sx q[0];
rz(0.43088669) q[0];
rz(3.131065) q[2];
sx q[2];
rz(-1.4665571) q[2];
sx q[2];
rz(2.5679905) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.745137) q[1];
sx q[1];
rz(-1.0978699) q[1];
sx q[1];
rz(0.84348444) q[1];
rz(1.0799418) q[3];
sx q[3];
rz(-0.98831359) q[3];
sx q[3];
rz(-0.39566186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4744711) q[2];
sx q[2];
rz(-2.4690364) q[2];
sx q[2];
rz(-0.39620623) q[2];
rz(0.19043663) q[3];
sx q[3];
rz(-0.93947828) q[3];
sx q[3];
rz(0.52085352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1162972) q[0];
sx q[0];
rz(-2.681356) q[0];
sx q[0];
rz(2.7468371) q[0];
rz(1.0074298) q[1];
sx q[1];
rz(-1.1578683) q[1];
sx q[1];
rz(1.6068858) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7028977) q[0];
sx q[0];
rz(-0.83897018) q[0];
sx q[0];
rz(1.954506) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41409576) q[2];
sx q[2];
rz(-0.76024917) q[2];
sx q[2];
rz(2.9000226) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.462127) q[1];
sx q[1];
rz(-1.7662342) q[1];
sx q[1];
rz(-1.2567169) q[1];
x q[2];
rz(-2.6768489) q[3];
sx q[3];
rz(-1.9806974) q[3];
sx q[3];
rz(1.003554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.69636238) q[2];
sx q[2];
rz(-2.4424489) q[2];
sx q[2];
rz(0.17808476) q[2];
rz(2.1961424) q[3];
sx q[3];
rz(-0.33792308) q[3];
sx q[3];
rz(2.7054355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4431385) q[0];
sx q[0];
rz(-0.66265023) q[0];
sx q[0];
rz(-2.9887001) q[0];
rz(-2.1235509) q[1];
sx q[1];
rz(-0.94296229) q[1];
sx q[1];
rz(-0.61000383) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40326443) q[0];
sx q[0];
rz(-0.77304196) q[0];
sx q[0];
rz(0.18456642) q[0];
rz(-2.1446682) q[2];
sx q[2];
rz(-1.5776433) q[2];
sx q[2];
rz(0.29563841) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9023561) q[1];
sx q[1];
rz(-0.78241759) q[1];
sx q[1];
rz(-2.8066325) q[1];
rz(0.69645564) q[3];
sx q[3];
rz(-0.12731931) q[3];
sx q[3];
rz(-2.6501973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6792004) q[2];
sx q[2];
rz(-1.2827164) q[2];
sx q[2];
rz(0.25897762) q[2];
rz(0.14608832) q[3];
sx q[3];
rz(-2.3688705) q[3];
sx q[3];
rz(0.63905382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6530957) q[0];
sx q[0];
rz(-2.275832) q[0];
sx q[0];
rz(0.75188941) q[0];
rz(-0.73196661) q[1];
sx q[1];
rz(-1.3804133) q[1];
sx q[1];
rz(0.52925777) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2471591) q[0];
sx q[0];
rz(-0.34810796) q[0];
sx q[0];
rz(2.0624119) q[0];
rz(-pi) q[1];
rz(-3.0050982) q[2];
sx q[2];
rz(-1.1195099) q[2];
sx q[2];
rz(2.9773444) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4715035) q[1];
sx q[1];
rz(-1.5836599) q[1];
sx q[1];
rz(-1.4895346) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1629754) q[3];
sx q[3];
rz(-2.6088723) q[3];
sx q[3];
rz(2.5902093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.25972128) q[2];
sx q[2];
rz(-2.112381) q[2];
sx q[2];
rz(-0.79505801) q[2];
rz(1.9184387) q[3];
sx q[3];
rz(-1.7984248) q[3];
sx q[3];
rz(2.3651626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.935598) q[0];
sx q[0];
rz(-1.5933651) q[0];
sx q[0];
rz(0.11496168) q[0];
rz(0.089275442) q[1];
sx q[1];
rz(-2.142579) q[1];
sx q[1];
rz(0.1549391) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6754571) q[0];
sx q[0];
rz(-2.977006) q[0];
sx q[0];
rz(0.91187061) q[0];
rz(0.55122113) q[2];
sx q[2];
rz(-2.1586426) q[2];
sx q[2];
rz(-1.2767222) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7748002) q[1];
sx q[1];
rz(-1.3143365) q[1];
sx q[1];
rz(0.40089295) q[1];
x q[2];
rz(-0.84568923) q[3];
sx q[3];
rz(-1.6250623) q[3];
sx q[3];
rz(-2.3783663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0838919) q[2];
sx q[2];
rz(-1.8996779) q[2];
sx q[2];
rz(-2.4532301) q[2];
rz(-0.56811959) q[3];
sx q[3];
rz(-1.0024242) q[3];
sx q[3];
rz(-2.4585371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.6883009) q[0];
sx q[0];
rz(-0.81128565) q[0];
sx q[0];
rz(1.9360833) q[0];
rz(2.0506809) q[1];
sx q[1];
rz(-2.5542407) q[1];
sx q[1];
rz(1.5376512) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6610049) q[0];
sx q[0];
rz(-1.2086226) q[0];
sx q[0];
rz(-2.037338) q[0];
rz(1.5715661) q[2];
sx q[2];
rz(-2.5736897) q[2];
sx q[2];
rz(-0.22837328) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4548757) q[1];
sx q[1];
rz(-1.4595965) q[1];
sx q[1];
rz(-2.3478572) q[1];
rz(1.0991715) q[3];
sx q[3];
rz(-2.1159008) q[3];
sx q[3];
rz(0.68600149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.18498147) q[2];
sx q[2];
rz(-1.1143755) q[2];
sx q[2];
rz(1.9967009) q[2];
rz(-1.6635118) q[3];
sx q[3];
rz(-2.3749115) q[3];
sx q[3];
rz(-2.0293106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0612653) q[0];
sx q[0];
rz(-0.63617951) q[0];
sx q[0];
rz(1.7058477) q[0];
rz(-1.1371293) q[1];
sx q[1];
rz(-0.69157332) q[1];
sx q[1];
rz(0.93708509) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81737751) q[0];
sx q[0];
rz(-2.827001) q[0];
sx q[0];
rz(-2.8552575) q[0];
rz(-pi) q[1];
rz(1.1446196) q[2];
sx q[2];
rz(-2.0405031) q[2];
sx q[2];
rz(-1.5091648) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.67928606) q[1];
sx q[1];
rz(-2.3704154) q[1];
sx q[1];
rz(-3.1047526) q[1];
rz(-pi) q[2];
rz(-2.5255975) q[3];
sx q[3];
rz(-1.8503891) q[3];
sx q[3];
rz(0.071690138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.81637853) q[2];
sx q[2];
rz(-0.64322317) q[2];
sx q[2];
rz(0.31036672) q[2];
rz(-0.54272932) q[3];
sx q[3];
rz(-0.92971814) q[3];
sx q[3];
rz(-0.31699666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19001374) q[0];
sx q[0];
rz(-2.0484476) q[0];
sx q[0];
rz(0.53566326) q[0];
rz(0.71470064) q[1];
sx q[1];
rz(-1.2477881) q[1];
sx q[1];
rz(1.4664149) q[1];
rz(-0.43595366) q[2];
sx q[2];
rz(-1.2203078) q[2];
sx q[2];
rz(0.75390799) q[2];
rz(1.6313995) q[3];
sx q[3];
rz(-2.2987859) q[3];
sx q[3];
rz(1.4684341) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
