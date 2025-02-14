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
rz(0.49993604) q[0];
sx q[0];
rz(-1.191782) q[0];
sx q[0];
rz(1.7697822) q[0];
rz(1.6504047) q[1];
sx q[1];
rz(-0.86725441) q[1];
sx q[1];
rz(0.39630085) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4258408) q[0];
sx q[0];
rz(-1.4914728) q[0];
sx q[0];
rz(2.9315445) q[0];
x q[1];
rz(0.17163817) q[2];
sx q[2];
rz(-1.2002103) q[2];
sx q[2];
rz(2.5570392) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.635879) q[1];
sx q[1];
rz(-1.6649655) q[1];
sx q[1];
rz(2.0999184) q[1];
rz(0.39537786) q[3];
sx q[3];
rz(-1.8413183) q[3];
sx q[3];
rz(0.45365712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0013915) q[2];
sx q[2];
rz(-2.3346257) q[2];
sx q[2];
rz(1.3414471) q[2];
rz(-2.0724824) q[3];
sx q[3];
rz(-0.1461229) q[3];
sx q[3];
rz(1.3816381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57217252) q[0];
sx q[0];
rz(-0.91745806) q[0];
sx q[0];
rz(2.4826352) q[0];
rz(0.072703687) q[1];
sx q[1];
rz(-2.0854988) q[1];
sx q[1];
rz(-1.8812077) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9958268) q[0];
sx q[0];
rz(-0.46016177) q[0];
sx q[0];
rz(1.3868679) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.057823) q[2];
sx q[2];
rz(-1.1758846) q[2];
sx q[2];
rz(2.9213443) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1932133) q[1];
sx q[1];
rz(-2.8682289) q[1];
sx q[1];
rz(2.8990919) q[1];
rz(-pi) q[2];
rz(2.1815592) q[3];
sx q[3];
rz(-1.6690147) q[3];
sx q[3];
rz(1.7383505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.8947123) q[2];
sx q[2];
rz(-2.3983045) q[2];
sx q[2];
rz(-1.3512208) q[2];
rz(-0.61659914) q[3];
sx q[3];
rz(-2.1117881) q[3];
sx q[3];
rz(0.29115796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.63224822) q[0];
sx q[0];
rz(-0.46724874) q[0];
sx q[0];
rz(0.21009357) q[0];
rz(0.083077438) q[1];
sx q[1];
rz(-1.9030842) q[1];
sx q[1];
rz(-3.074379) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9191362) q[0];
sx q[0];
rz(-1.558715) q[0];
sx q[0];
rz(3.1280845) q[0];
rz(-2.1234197) q[2];
sx q[2];
rz(-0.45028307) q[2];
sx q[2];
rz(-1.0339111) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9720502) q[1];
sx q[1];
rz(-1.9055371) q[1];
sx q[1];
rz(0.85595815) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4806858) q[3];
sx q[3];
rz(-2.2743974) q[3];
sx q[3];
rz(1.4058324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3136966) q[2];
sx q[2];
rz(-1.4915497) q[2];
sx q[2];
rz(1.3649887) q[2];
rz(0.40417534) q[3];
sx q[3];
rz(-1.3173236) q[3];
sx q[3];
rz(1.1990168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-0.86554027) q[0];
sx q[0];
rz(-1.7409538) q[0];
sx q[0];
rz(-0.47819594) q[0];
rz(-2.1887691) q[1];
sx q[1];
rz(-0.15004221) q[1];
sx q[1];
rz(2.731954) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7573625) q[0];
sx q[0];
rz(-0.91529796) q[0];
sx q[0];
rz(-1.537498) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.19721376) q[2];
sx q[2];
rz(-0.46425113) q[2];
sx q[2];
rz(-2.0418389) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2634025) q[1];
sx q[1];
rz(-0.95118427) q[1];
sx q[1];
rz(0.81677572) q[1];
rz(-1.8919935) q[3];
sx q[3];
rz(-1.7068591) q[3];
sx q[3];
rz(-2.4709159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5459368) q[2];
sx q[2];
rz(-1.7031534) q[2];
sx q[2];
rz(-2.0036009) q[2];
rz(0.99500895) q[3];
sx q[3];
rz(-2.5840839) q[3];
sx q[3];
rz(0.086624302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.397641) q[0];
sx q[0];
rz(-1.8331563) q[0];
sx q[0];
rz(0.27807903) q[0];
rz(1.2644348) q[1];
sx q[1];
rz(-1.8587298) q[1];
sx q[1];
rz(-1.5778731) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20533097) q[0];
sx q[0];
rz(-0.38547984) q[0];
sx q[0];
rz(1.2275213) q[0];
rz(0.78863849) q[2];
sx q[2];
rz(-1.83162) q[2];
sx q[2];
rz(-0.80718416) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1335781) q[1];
sx q[1];
rz(-2.1431794) q[1];
sx q[1];
rz(-1.7126669) q[1];
x q[2];
rz(-0.092838959) q[3];
sx q[3];
rz(-1.412409) q[3];
sx q[3];
rz(3.1166589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9331253) q[2];
sx q[2];
rz(-1.5310023) q[2];
sx q[2];
rz(2.7739286) q[2];
rz(1.095088) q[3];
sx q[3];
rz(-2.118066) q[3];
sx q[3];
rz(0.17525214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4472189) q[0];
sx q[0];
rz(-0.29964888) q[0];
sx q[0];
rz(-1.8023941) q[0];
rz(0.12204349) q[1];
sx q[1];
rz(-1.782878) q[1];
sx q[1];
rz(0.36062127) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0568389) q[0];
sx q[0];
rz(-0.90843102) q[0];
sx q[0];
rz(-1.0785036) q[0];
x q[1];
rz(1.1507658) q[2];
sx q[2];
rz(-1.5345062) q[2];
sx q[2];
rz(1.0880214) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.079626849) q[1];
sx q[1];
rz(-2.6089416) q[1];
sx q[1];
rz(2.1497186) q[1];
rz(-pi) q[2];
rz(0.25288323) q[3];
sx q[3];
rz(-1.5759625) q[3];
sx q[3];
rz(0.21765003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.7755166) q[2];
sx q[2];
rz(-1.6785494) q[2];
sx q[2];
rz(-1.1386846) q[2];
rz(-0.25389296) q[3];
sx q[3];
rz(-0.64520276) q[3];
sx q[3];
rz(0.74786389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2426386) q[0];
sx q[0];
rz(-0.88976088) q[0];
sx q[0];
rz(2.1658072) q[0];
rz(1.1294533) q[1];
sx q[1];
rz(-0.44810805) q[1];
sx q[1];
rz(1.7879558) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0669374) q[0];
sx q[0];
rz(-2.2569509) q[0];
sx q[0];
rz(-1.9524491) q[0];
rz(-pi) q[1];
rz(1.7464306) q[2];
sx q[2];
rz(-0.43000107) q[2];
sx q[2];
rz(1.2295837) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5013988) q[1];
sx q[1];
rz(-0.36194776) q[1];
sx q[1];
rz(-2.7622403) q[1];
rz(-pi) q[2];
rz(-2.4134992) q[3];
sx q[3];
rz(-0.59641664) q[3];
sx q[3];
rz(-2.8567258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5740616) q[2];
sx q[2];
rz(-0.83063829) q[2];
sx q[2];
rz(2.2243824) q[2];
rz(-1.5926825) q[3];
sx q[3];
rz(-0.33532381) q[3];
sx q[3];
rz(2.3622021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2061283) q[0];
sx q[0];
rz(-1.520282) q[0];
sx q[0];
rz(-0.024854831) q[0];
rz(1.8553597) q[1];
sx q[1];
rz(-1.356946) q[1];
sx q[1];
rz(1.6110427) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0760374) q[0];
sx q[0];
rz(-2.5352056) q[0];
sx q[0];
rz(0.72640149) q[0];
rz(-1.3127682) q[2];
sx q[2];
rz(-2.240643) q[2];
sx q[2];
rz(0.19425288) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.81653412) q[1];
sx q[1];
rz(-1.5360502) q[1];
sx q[1];
rz(-0.28157061) q[1];
rz(-pi) q[2];
rz(-0.99175158) q[3];
sx q[3];
rz(-2.8914521) q[3];
sx q[3];
rz(2.5219805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1207235) q[2];
sx q[2];
rz(-2.1376231) q[2];
sx q[2];
rz(-1.8570159) q[2];
rz(0.32127109) q[3];
sx q[3];
rz(-0.64906859) q[3];
sx q[3];
rz(0.2215213) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.016841737) q[0];
sx q[0];
rz(-1.9945972) q[0];
sx q[0];
rz(0.91957134) q[0];
rz(-2.549767) q[1];
sx q[1];
rz(-1.0706736) q[1];
sx q[1];
rz(2.8048973) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29558173) q[0];
sx q[0];
rz(-2.8407556) q[0];
sx q[0];
rz(-1.1016125) q[0];
rz(-pi) q[1];
rz(2.8317802) q[2];
sx q[2];
rz(-2.7618119) q[2];
sx q[2];
rz(-1.7526232) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.10178653) q[1];
sx q[1];
rz(-2.617604) q[1];
sx q[1];
rz(-0.26079674) q[1];
rz(-pi) q[2];
rz(-1.3161427) q[3];
sx q[3];
rz(-1.4166204) q[3];
sx q[3];
rz(-1.7642989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0261953) q[2];
sx q[2];
rz(-0.43033174) q[2];
sx q[2];
rz(2.3587312) q[2];
rz(-3.03249) q[3];
sx q[3];
rz(-1.0166054) q[3];
sx q[3];
rz(1.8946064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1054909) q[0];
sx q[0];
rz(-1.6721268) q[0];
sx q[0];
rz(0.92392695) q[0];
rz(2.6149514) q[1];
sx q[1];
rz(-1.275332) q[1];
sx q[1];
rz(-2.142876) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.625542) q[0];
sx q[0];
rz(-1.7241417) q[0];
sx q[0];
rz(1.6244066) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2998215) q[2];
sx q[2];
rz(-2.8405115) q[2];
sx q[2];
rz(-1.4512514) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0622039) q[1];
sx q[1];
rz(-0.98323373) q[1];
sx q[1];
rz(-0.04157898) q[1];
rz(-pi) q[2];
rz(-2.4630354) q[3];
sx q[3];
rz(-1.2781837) q[3];
sx q[3];
rz(1.34174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0052789) q[2];
sx q[2];
rz(-2.2515991) q[2];
sx q[2];
rz(-2.8078553) q[2];
rz(-0.8606832) q[3];
sx q[3];
rz(-1.4849097) q[3];
sx q[3];
rz(3.1349643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
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
rz(-2.6776047) q[0];
sx q[0];
rz(-1.9515568) q[0];
sx q[0];
rz(-2.4757181) q[0];
rz(-0.55466501) q[1];
sx q[1];
rz(-1.278109) q[1];
sx q[1];
rz(0.14229933) q[1];
rz(0.43989506) q[2];
sx q[2];
rz(-1.7976201) q[2];
sx q[2];
rz(1.3646361) q[2];
rz(-2.4246115) q[3];
sx q[3];
rz(-2.3752799) q[3];
sx q[3];
rz(0.3175288) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
