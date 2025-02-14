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
rz(2.0789335) q[0];
sx q[0];
rz(-2.3261429) q[0];
sx q[0];
rz(2.4315779) q[0];
rz(-0.31399909) q[1];
sx q[1];
rz(-0.93500885) q[1];
sx q[1];
rz(1.8097872) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3478617) q[0];
sx q[0];
rz(-1.8457185) q[0];
sx q[0];
rz(1.8040276) q[0];
x q[1];
rz(2.5153036) q[2];
sx q[2];
rz(-0.41538996) q[2];
sx q[2];
rz(-2.9363869) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0441452) q[1];
sx q[1];
rz(-2.4821979) q[1];
sx q[1];
rz(-1.161518) q[1];
rz(-2.136257) q[3];
sx q[3];
rz(-2.4255803) q[3];
sx q[3];
rz(-2.211192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4189202) q[2];
sx q[2];
rz(-1.924943) q[2];
sx q[2];
rz(2.326272) q[2];
rz(-0.80960649) q[3];
sx q[3];
rz(-1.4987192) q[3];
sx q[3];
rz(-2.0577551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0631436) q[0];
sx q[0];
rz(-2.0922631) q[0];
sx q[0];
rz(3.0310042) q[0];
rz(-0.94509205) q[1];
sx q[1];
rz(-2.7022305) q[1];
sx q[1];
rz(1.5623215) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54540578) q[0];
sx q[0];
rz(-0.89841398) q[0];
sx q[0];
rz(-0.59242146) q[0];
x q[1];
rz(2.7209372) q[2];
sx q[2];
rz(-1.1569945) q[2];
sx q[2];
rz(-2.5977573) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.72772273) q[1];
sx q[1];
rz(-1.4663823) q[1];
sx q[1];
rz(2.6563679) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0575297) q[3];
sx q[3];
rz(-0.81765122) q[3];
sx q[3];
rz(-0.14150894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9940146) q[2];
sx q[2];
rz(-0.14573228) q[2];
sx q[2];
rz(0.4501403) q[2];
rz(1.133793) q[3];
sx q[3];
rz(-1.6885875) q[3];
sx q[3];
rz(0.97761893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56871539) q[0];
sx q[0];
rz(-2.1935538) q[0];
sx q[0];
rz(0.28133389) q[0];
rz(-1.3866792) q[1];
sx q[1];
rz(-1.4298871) q[1];
sx q[1];
rz(2.5414355) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2413797) q[0];
sx q[0];
rz(-2.8474387) q[0];
sx q[0];
rz(-2.9808729) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.96397206) q[2];
sx q[2];
rz(-1.5247048) q[2];
sx q[2];
rz(-0.54296903) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.197613) q[1];
sx q[1];
rz(-0.67786067) q[1];
sx q[1];
rz(-2.6070836) q[1];
rz(-pi) q[2];
rz(-0.26140499) q[3];
sx q[3];
rz(-2.5525744) q[3];
sx q[3];
rz(-0.95688577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0967789) q[2];
sx q[2];
rz(-1.1324977) q[2];
sx q[2];
rz(1.0463932) q[2];
rz(0.3977631) q[3];
sx q[3];
rz(-2.4989276) q[3];
sx q[3];
rz(0.1325632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9699049) q[0];
sx q[0];
rz(-3.1145018) q[0];
sx q[0];
rz(1.0525674) q[0];
rz(1.196208) q[1];
sx q[1];
rz(-1.9320678) q[1];
sx q[1];
rz(2.722091) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2118111) q[0];
sx q[0];
rz(-2.503184) q[0];
sx q[0];
rz(-0.45244658) q[0];
rz(2.9151462) q[2];
sx q[2];
rz(-2.1927877) q[2];
sx q[2];
rz(-2.7055307) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4070523) q[1];
sx q[1];
rz(-2.190935) q[1];
sx q[1];
rz(0.46874115) q[1];
rz(-pi) q[2];
rz(-1.2798843) q[3];
sx q[3];
rz(-2.0613163) q[3];
sx q[3];
rz(1.377493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9516051) q[2];
sx q[2];
rz(-1.1185458) q[2];
sx q[2];
rz(-2.0513963) q[2];
rz(2.6103141) q[3];
sx q[3];
rz(-2.4254906) q[3];
sx q[3];
rz(-1.5268911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(1.4215609) q[0];
sx q[0];
rz(-1.5694542) q[0];
sx q[0];
rz(-1.3979727) q[0];
rz(-3.0086503) q[1];
sx q[1];
rz(-1.3016737) q[1];
sx q[1];
rz(-3.1403819) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5842297) q[0];
sx q[0];
rz(-1.5869291) q[0];
sx q[0];
rz(2.7610012) q[0];
rz(1.2852077) q[2];
sx q[2];
rz(-1.4175709) q[2];
sx q[2];
rz(1.1040579) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.397314) q[1];
sx q[1];
rz(-0.71956735) q[1];
sx q[1];
rz(2.6707763) q[1];
rz(1.7204956) q[3];
sx q[3];
rz(-1.5970917) q[3];
sx q[3];
rz(0.61248518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2289537) q[2];
sx q[2];
rz(-1.3372083) q[2];
sx q[2];
rz(2.9310628) q[2];
rz(-2.1052965) q[3];
sx q[3];
rz(-2.3573037) q[3];
sx q[3];
rz(-1.9265296) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9482816) q[0];
sx q[0];
rz(-1.9183777) q[0];
sx q[0];
rz(0.40570983) q[0];
rz(-1.7871208) q[1];
sx q[1];
rz(-2.3190119) q[1];
sx q[1];
rz(0.92232651) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2012537) q[0];
sx q[0];
rz(-0.30863276) q[0];
sx q[0];
rz(-1.2135394) q[0];
rz(-pi) q[1];
rz(-1.7013266) q[2];
sx q[2];
rz(-0.67085999) q[2];
sx q[2];
rz(2.6595569) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.74444425) q[1];
sx q[1];
rz(-1.6849298) q[1];
sx q[1];
rz(2.1836917) q[1];
rz(-pi) q[2];
rz(0.57862307) q[3];
sx q[3];
rz(-2.4994183) q[3];
sx q[3];
rz(-2.7901621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4313844) q[2];
sx q[2];
rz(-1.9269383) q[2];
sx q[2];
rz(-2.7396438) q[2];
rz(1.8534144) q[3];
sx q[3];
rz(-0.86789075) q[3];
sx q[3];
rz(1.8624381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5227018) q[0];
sx q[0];
rz(-1.1033449) q[0];
sx q[0];
rz(0.064362137) q[0];
rz(2.3035658) q[1];
sx q[1];
rz(-2.3152654) q[1];
sx q[1];
rz(1.3305957) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15312402) q[0];
sx q[0];
rz(-0.58725427) q[0];
sx q[0];
rz(-0.22756094) q[0];
x q[1];
rz(-1.072791) q[2];
sx q[2];
rz(-2.118022) q[2];
sx q[2];
rz(-2.4233873) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8992675) q[1];
sx q[1];
rz(-1.0501672) q[1];
sx q[1];
rz(-1.316687) q[1];
x q[2];
rz(-2.182992) q[3];
sx q[3];
rz(-1.7122088) q[3];
sx q[3];
rz(-2.8116016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7360721) q[2];
sx q[2];
rz(-1.5084927) q[2];
sx q[2];
rz(-2.4580477) q[2];
rz(-2.4077967) q[3];
sx q[3];
rz(-1.7283231) q[3];
sx q[3];
rz(2.2939513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7965294) q[0];
sx q[0];
rz(-2.0579484) q[0];
sx q[0];
rz(-1.1731359) q[0];
rz(2.7660811) q[1];
sx q[1];
rz(-2.651732) q[1];
sx q[1];
rz(-1.0367905) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79085892) q[0];
sx q[0];
rz(-0.030412523) q[0];
sx q[0];
rz(2.4615672) q[0];
rz(-2.9762245) q[2];
sx q[2];
rz(-1.2980808) q[2];
sx q[2];
rz(-0.80406666) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4040096) q[1];
sx q[1];
rz(-0.77464691) q[1];
sx q[1];
rz(-2.9402016) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6508207) q[3];
sx q[3];
rz(-1.0939071) q[3];
sx q[3];
rz(0.27974883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.57572395) q[2];
sx q[2];
rz(-2.2278991) q[2];
sx q[2];
rz(-0.73053989) q[2];
rz(-2.9324487) q[3];
sx q[3];
rz(-2.6181965) q[3];
sx q[3];
rz(1.8714347) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4850979) q[0];
sx q[0];
rz(-2.7182343) q[0];
sx q[0];
rz(-3.094161) q[0];
rz(0.221953) q[1];
sx q[1];
rz(-2.0675979) q[1];
sx q[1];
rz(-2.2304631) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.297356) q[0];
sx q[0];
rz(-0.043023303) q[0];
sx q[0];
rz(2.4639315) q[0];
x q[1];
rz(1.5510727) q[2];
sx q[2];
rz(-0.64177401) q[2];
sx q[2];
rz(0.80427792) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.641727) q[1];
sx q[1];
rz(-1.2033101) q[1];
sx q[1];
rz(1.5843452) q[1];
x q[2];
rz(0.20424517) q[3];
sx q[3];
rz(-1.0139272) q[3];
sx q[3];
rz(-2.000252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.86965108) q[2];
sx q[2];
rz(-0.38001529) q[2];
sx q[2];
rz(2.9730049) q[2];
rz(2.9224959) q[3];
sx q[3];
rz(-0.85182652) q[3];
sx q[3];
rz(1.3271837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.99437) q[0];
sx q[0];
rz(-0.34039012) q[0];
sx q[0];
rz(-0.45743531) q[0];
rz(-1.2262729) q[1];
sx q[1];
rz(-2.8044082) q[1];
sx q[1];
rz(1.3630684) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4922614) q[0];
sx q[0];
rz(-0.96423414) q[0];
sx q[0];
rz(-1.6442293) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.326637) q[2];
sx q[2];
rz(-0.95521046) q[2];
sx q[2];
rz(-2.2509839) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5536834) q[1];
sx q[1];
rz(-1.739555) q[1];
sx q[1];
rz(0.38270271) q[1];
rz(-3.0567597) q[3];
sx q[3];
rz(-2.2515577) q[3];
sx q[3];
rz(0.57455243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3739796) q[2];
sx q[2];
rz(-0.55459905) q[2];
sx q[2];
rz(0.45200959) q[2];
rz(2.7376145) q[3];
sx q[3];
rz(-1.890506) q[3];
sx q[3];
rz(-1.1261136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44611888) q[0];
sx q[0];
rz(-1.5344545) q[0];
sx q[0];
rz(-2.9547966) q[0];
rz(0.52234621) q[1];
sx q[1];
rz(-2.6112687) q[1];
sx q[1];
rz(-1.912259) q[1];
rz(-2.859237) q[2];
sx q[2];
rz(-1.3236125) q[2];
sx q[2];
rz(-0.57409928) q[2];
rz(2.6330144) q[3];
sx q[3];
rz(-2.8420035) q[3];
sx q[3];
rz(-3.0723078) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
