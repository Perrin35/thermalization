OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5300712) q[0];
sx q[0];
rz(-2.0687215) q[0];
sx q[0];
rz(-1.7142417) q[0];
rz(-2.144835) q[1];
sx q[1];
rz(-0.61350322) q[1];
sx q[1];
rz(3.0761392) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50953853) q[0];
sx q[0];
rz(-0.51842123) q[0];
sx q[0];
rz(-2.753756) q[0];
rz(-pi) q[1];
rz(0.36875002) q[2];
sx q[2];
rz(-1.6095709) q[2];
sx q[2];
rz(-1.4132959) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1004384) q[1];
sx q[1];
rz(-1.8722495) q[1];
sx q[1];
rz(0.90561066) q[1];
x q[2];
rz(2.0569589) q[3];
sx q[3];
rz(-1.0118985) q[3];
sx q[3];
rz(-1.3317747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.87302152) q[2];
sx q[2];
rz(-0.34653386) q[2];
sx q[2];
rz(1.8447426) q[2];
rz(2.0348564) q[3];
sx q[3];
rz(-0.96702558) q[3];
sx q[3];
rz(-1.4510252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1065555) q[0];
sx q[0];
rz(-0.7283926) q[0];
sx q[0];
rz(-1.6075217) q[0];
rz(-0.84397498) q[1];
sx q[1];
rz(-1.6232099) q[1];
sx q[1];
rz(2.8704876) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8753238) q[0];
sx q[0];
rz(-1.8401339) q[0];
sx q[0];
rz(1.6104417) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9366802) q[2];
sx q[2];
rz(-0.26981631) q[2];
sx q[2];
rz(1.3617977) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9535189) q[1];
sx q[1];
rz(-1.9862576) q[1];
sx q[1];
rz(1.6759592) q[1];
x q[2];
rz(-0.45174349) q[3];
sx q[3];
rz(-2.9752033) q[3];
sx q[3];
rz(1.1903669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8035182) q[2];
sx q[2];
rz(-0.98793554) q[2];
sx q[2];
rz(-0.79338497) q[2];
rz(0.042898305) q[3];
sx q[3];
rz(-1.2267313) q[3];
sx q[3];
rz(-2.4734917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0637829) q[0];
sx q[0];
rz(-0.48352799) q[0];
sx q[0];
rz(0.10312816) q[0];
rz(0.25281301) q[1];
sx q[1];
rz(-1.4272855) q[1];
sx q[1];
rz(2.8275183) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.26337) q[0];
sx q[0];
rz(-0.65257699) q[0];
sx q[0];
rz(-0.67499365) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0843749) q[2];
sx q[2];
rz(-1.0891181) q[2];
sx q[2];
rz(1.526598) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1669877) q[1];
sx q[1];
rz(-1.0613975) q[1];
sx q[1];
rz(2.2231119) q[1];
rz(0.57468412) q[3];
sx q[3];
rz(-0.30467027) q[3];
sx q[3];
rz(1.6449354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0146497) q[2];
sx q[2];
rz(-2.5003771) q[2];
sx q[2];
rz(-2.3250697) q[2];
rz(2.4294295) q[3];
sx q[3];
rz(-1.4543337) q[3];
sx q[3];
rz(-0.85675353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4954869) q[0];
sx q[0];
rz(-3.0480338) q[0];
sx q[0];
rz(2.5529472) q[0];
rz(-0.37216392) q[1];
sx q[1];
rz(-1.2969505) q[1];
sx q[1];
rz(2.1968496) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78658453) q[0];
sx q[0];
rz(-1.8898801) q[0];
sx q[0];
rz(-3.0727076) q[0];
rz(-pi) q[1];
rz(1.107457) q[2];
sx q[2];
rz(-2.6012408) q[2];
sx q[2];
rz(-2.5519443) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2887047) q[1];
sx q[1];
rz(-1.0094704) q[1];
sx q[1];
rz(1.0594548) q[1];
x q[2];
rz(0.1907986) q[3];
sx q[3];
rz(-1.2114297) q[3];
sx q[3];
rz(-0.99890733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9113691) q[2];
sx q[2];
rz(-1.5350124) q[2];
sx q[2];
rz(1.4851419) q[2];
rz(2.7459775) q[3];
sx q[3];
rz(-1.6499358) q[3];
sx q[3];
rz(-3.0806471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7570801) q[0];
sx q[0];
rz(-0.39329305) q[0];
sx q[0];
rz(-1.8994037) q[0];
rz(3.0643265) q[1];
sx q[1];
rz(-1.9237513) q[1];
sx q[1];
rz(-0.088931106) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1031694) q[0];
sx q[0];
rz(-0.22617243) q[0];
sx q[0];
rz(-2.3875176) q[0];
rz(1.6709264) q[2];
sx q[2];
rz(-1.0271974) q[2];
sx q[2];
rz(1.3058654) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6731074) q[1];
sx q[1];
rz(-1.846828) q[1];
sx q[1];
rz(-1.0795781) q[1];
rz(-pi) q[2];
rz(-1.704028) q[3];
sx q[3];
rz(-1.4733088) q[3];
sx q[3];
rz(-2.9436802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.58101216) q[2];
sx q[2];
rz(-1.3835013) q[2];
sx q[2];
rz(-1.9182473) q[2];
rz(1.2286202) q[3];
sx q[3];
rz(-0.88594985) q[3];
sx q[3];
rz(-0.74738735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.05518797) q[0];
sx q[0];
rz(-2.3265657) q[0];
sx q[0];
rz(0.95274693) q[0];
rz(-1.791026) q[1];
sx q[1];
rz(-1.2341276) q[1];
sx q[1];
rz(-1.9780673) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6442808) q[0];
sx q[0];
rz(-1.4587879) q[0];
sx q[0];
rz(-0.73581672) q[0];
rz(-pi) q[1];
rz(1.8356531) q[2];
sx q[2];
rz(-0.9914757) q[2];
sx q[2];
rz(0.91897041) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.64505534) q[1];
sx q[1];
rz(-0.24894196) q[1];
sx q[1];
rz(-1.7824936) q[1];
rz(-0.58168488) q[3];
sx q[3];
rz(-0.69184985) q[3];
sx q[3];
rz(-2.8210616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.20840883) q[2];
sx q[2];
rz(-1.9219425) q[2];
sx q[2];
rz(-1.5737994) q[2];
rz(2.9257704) q[3];
sx q[3];
rz(-0.62041557) q[3];
sx q[3];
rz(-0.1568493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1167574) q[0];
sx q[0];
rz(-0.76501608) q[0];
sx q[0];
rz(-1.5417954) q[0];
rz(0.4190017) q[1];
sx q[1];
rz(-2.0332789) q[1];
sx q[1];
rz(-1.7760743) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1772108) q[0];
sx q[0];
rz(-2.1809077) q[0];
sx q[0];
rz(0.23680738) q[0];
rz(-pi) q[1];
rz(0.78132479) q[2];
sx q[2];
rz(-1.1622687) q[2];
sx q[2];
rz(0.21121001) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8518062) q[1];
sx q[1];
rz(-0.86871457) q[1];
sx q[1];
rz(0.21721022) q[1];
x q[2];
rz(-1.0700052) q[3];
sx q[3];
rz(-0.65088256) q[3];
sx q[3];
rz(0.18903519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2447723) q[2];
sx q[2];
rz(-0.12910566) q[2];
sx q[2];
rz(-1.4531892) q[2];
rz(-1.7425273) q[3];
sx q[3];
rz(-1.941967) q[3];
sx q[3];
rz(0.80756584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6145265) q[0];
sx q[0];
rz(-0.87166059) q[0];
sx q[0];
rz(-2.6283619) q[0];
rz(-2.1649583) q[1];
sx q[1];
rz(-2.4726424) q[1];
sx q[1];
rz(1.4248779) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0289405) q[0];
sx q[0];
rz(-1.4412349) q[0];
sx q[0];
rz(0.86973637) q[0];
x q[1];
rz(2.5930357) q[2];
sx q[2];
rz(-1.6335677) q[2];
sx q[2];
rz(1.5323764) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.91814525) q[1];
sx q[1];
rz(-0.53358101) q[1];
sx q[1];
rz(2.1672638) q[1];
rz(-pi) q[2];
rz(-3.023671) q[3];
sx q[3];
rz(-1.4050845) q[3];
sx q[3];
rz(-0.57512257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6500515) q[2];
sx q[2];
rz(-2.1341133) q[2];
sx q[2];
rz(0.36637351) q[2];
rz(-2.5036687) q[3];
sx q[3];
rz(-1.7584691) q[3];
sx q[3];
rz(1.5095507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3731308) q[0];
sx q[0];
rz(-0.32485425) q[0];
sx q[0];
rz(1.1573867) q[0];
rz(0.53818446) q[1];
sx q[1];
rz(-1.8073852) q[1];
sx q[1];
rz(-3.11788) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0483646) q[0];
sx q[0];
rz(-2.1914838) q[0];
sx q[0];
rz(-1.4976504) q[0];
rz(1.1309409) q[2];
sx q[2];
rz(-0.80200101) q[2];
sx q[2];
rz(-1.8126876) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7754375) q[1];
sx q[1];
rz(-2.0733655) q[1];
sx q[1];
rz(-1.9910046) q[1];
rz(-pi) q[2];
rz(1.7718763) q[3];
sx q[3];
rz(-1.6938801) q[3];
sx q[3];
rz(1.1986002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9862765) q[2];
sx q[2];
rz(-1.3331058) q[2];
sx q[2];
rz(-2.7434366) q[2];
rz(2.6768173) q[3];
sx q[3];
rz(-2.0538752) q[3];
sx q[3];
rz(1.8006511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5584797) q[0];
sx q[0];
rz(-0.3322424) q[0];
sx q[0];
rz(-0.94859052) q[0];
rz(0.96725431) q[1];
sx q[1];
rz(-1.9890669) q[1];
sx q[1];
rz(-2.1327877) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43320424) q[0];
sx q[0];
rz(-0.77022591) q[0];
sx q[0];
rz(-1.5608833) q[0];
rz(-pi) q[1];
rz(2.45935) q[2];
sx q[2];
rz(-2.3385404) q[2];
sx q[2];
rz(-2.19953) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.29459144) q[1];
sx q[1];
rz(-1.090511) q[1];
sx q[1];
rz(-1.7748407) q[1];
rz(-2.4835577) q[3];
sx q[3];
rz(-1.3962126) q[3];
sx q[3];
rz(0.10887077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.41445109) q[2];
sx q[2];
rz(-0.87697566) q[2];
sx q[2];
rz(0.021075185) q[2];
rz(0.4840788) q[3];
sx q[3];
rz(-2.0120072) q[3];
sx q[3];
rz(2.4626203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8503583) q[0];
sx q[0];
rz(-1.5342916) q[0];
sx q[0];
rz(1.7271484) q[0];
rz(2.5454632) q[1];
sx q[1];
rz(-0.78501751) q[1];
sx q[1];
rz(2.6529978) q[1];
rz(0.1487371) q[2];
sx q[2];
rz(-1.7784468) q[2];
sx q[2];
rz(-0.39685984) q[2];
rz(-3.0335043) q[3];
sx q[3];
rz(-1.3664403) q[3];
sx q[3];
rz(-1.4667778) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
