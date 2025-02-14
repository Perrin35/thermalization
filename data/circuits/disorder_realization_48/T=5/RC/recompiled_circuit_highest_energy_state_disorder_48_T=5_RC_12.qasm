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
rz(-2.2263865) q[0];
sx q[0];
rz(-2.6829166) q[0];
sx q[0];
rz(-0.31804481) q[0];
rz(1.3660499) q[1];
sx q[1];
rz(-2.4429758) q[1];
sx q[1];
rz(-3.0787992) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1284769) q[0];
sx q[0];
rz(-2.4139691) q[0];
sx q[0];
rz(-0.81650556) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.43831011) q[2];
sx q[2];
rz(-0.93899124) q[2];
sx q[2];
rz(-2.8614728) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.56927189) q[1];
sx q[1];
rz(-0.57479984) q[1];
sx q[1];
rz(1.040609) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1413689) q[3];
sx q[3];
rz(-2.4834342) q[3];
sx q[3];
rz(-2.6702945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.723145) q[2];
sx q[2];
rz(-1.3157088) q[2];
sx q[2];
rz(-1.0203699) q[2];
rz(-2.4127035) q[3];
sx q[3];
rz(-0.43292361) q[3];
sx q[3];
rz(-1.5322022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51434022) q[0];
sx q[0];
rz(-1.1662551) q[0];
sx q[0];
rz(-0.038272055) q[0];
rz(0.12380869) q[1];
sx q[1];
rz(-2.6688711) q[1];
sx q[1];
rz(-1.570787) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4306612) q[0];
sx q[0];
rz(-1.5607906) q[0];
sx q[0];
rz(-1.5512933) q[0];
rz(-pi) q[1];
rz(-0.10199396) q[2];
sx q[2];
rz(-0.77384678) q[2];
sx q[2];
rz(0.71600658) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1727708) q[1];
sx q[1];
rz(-1.5712196) q[1];
sx q[1];
rz(-1.569773) q[1];
rz(-pi) q[2];
rz(-2.6255076) q[3];
sx q[3];
rz(-1.6760964) q[3];
sx q[3];
rz(2.1886231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6191285) q[2];
sx q[2];
rz(-1.3294486) q[2];
sx q[2];
rz(-0.5788571) q[2];
rz(-1.045643) q[3];
sx q[3];
rz(-0.78543109) q[3];
sx q[3];
rz(-0.75700179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66913644) q[0];
sx q[0];
rz(-2.0643015) q[0];
sx q[0];
rz(-0.45397595) q[0];
rz(2.9767735) q[1];
sx q[1];
rz(-2.3655128) q[1];
sx q[1];
rz(-1.4705315) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5002086) q[0];
sx q[0];
rz(-2.7140896) q[0];
sx q[0];
rz(-1.1469202) q[0];
x q[1];
rz(-1.0452095) q[2];
sx q[2];
rz(-1.282935) q[2];
sx q[2];
rz(-0.20341104) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1118033) q[1];
sx q[1];
rz(-1.8743427) q[1];
sx q[1];
rz(0.70036594) q[1];
rz(-pi) q[2];
rz(2.4824247) q[3];
sx q[3];
rz(-1.9999095) q[3];
sx q[3];
rz(1.3339213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.37956023) q[2];
sx q[2];
rz(-1.6635868) q[2];
sx q[2];
rz(-1.3564159) q[2];
rz(-0.49946579) q[3];
sx q[3];
rz(-1.5288589) q[3];
sx q[3];
rz(-1.0511901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22719638) q[0];
sx q[0];
rz(-0.90407404) q[0];
sx q[0];
rz(2.0830578) q[0];
rz(-1.9805485) q[1];
sx q[1];
rz(-0.5608905) q[1];
sx q[1];
rz(3.0243691) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6446021) q[0];
sx q[0];
rz(-2.1128901) q[0];
sx q[0];
rz(-1.6165471) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2160596) q[2];
sx q[2];
rz(-1.7286073) q[2];
sx q[2];
rz(0.35343808) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8414054) q[1];
sx q[1];
rz(-2.3267728) q[1];
sx q[1];
rz(-2.7542265) q[1];
rz(-pi) q[2];
rz(-2.6113308) q[3];
sx q[3];
rz(-0.70928364) q[3];
sx q[3];
rz(1.5269296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.78493541) q[2];
sx q[2];
rz(-1.4484118) q[2];
sx q[2];
rz(-1.3680722) q[2];
rz(-1.28537) q[3];
sx q[3];
rz(-2.0338438) q[3];
sx q[3];
rz(0.72803298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0990937) q[0];
sx q[0];
rz(-0.92785257) q[0];
sx q[0];
rz(1.4121144) q[0];
rz(-2.1389351) q[1];
sx q[1];
rz(-0.24191562) q[1];
sx q[1];
rz(1.5826506) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6928685) q[0];
sx q[0];
rz(-2.4131615) q[0];
sx q[0];
rz(0.65653519) q[0];
x q[1];
rz(-3.0465165) q[2];
sx q[2];
rz(-2.4209967) q[2];
sx q[2];
rz(0.61340767) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.11410275) q[1];
sx q[1];
rz(-1.5472425) q[1];
sx q[1];
rz(1.0204169) q[1];
rz(-pi) q[2];
rz(1.847397) q[3];
sx q[3];
rz(-2.8081354) q[3];
sx q[3];
rz(1.5302304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.58854181) q[2];
sx q[2];
rz(-1.43247) q[2];
sx q[2];
rz(-0.39349619) q[2];
rz(0.95997512) q[3];
sx q[3];
rz(-0.67592755) q[3];
sx q[3];
rz(-2.1737607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.2734964) q[0];
sx q[0];
rz(-3.0143026) q[0];
sx q[0];
rz(-0.18336503) q[0];
rz(2.6496437) q[1];
sx q[1];
rz(-0.79194561) q[1];
sx q[1];
rz(-1.2194182) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26579612) q[0];
sx q[0];
rz(-1.6599433) q[0];
sx q[0];
rz(-1.4899979) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6492789) q[2];
sx q[2];
rz(-1.6011136) q[2];
sx q[2];
rz(2.226269) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9907836) q[1];
sx q[1];
rz(-1.9766625) q[1];
sx q[1];
rz(-1.957914) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0571396) q[3];
sx q[3];
rz(-1.1005529) q[3];
sx q[3];
rz(-3.0601644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0101953) q[2];
sx q[2];
rz(-1.4364028) q[2];
sx q[2];
rz(-0.08237002) q[2];
rz(2.2149337) q[3];
sx q[3];
rz(-0.72946531) q[3];
sx q[3];
rz(-1.6389219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(2.6019186) q[0];
sx q[0];
rz(-2.813756) q[0];
sx q[0];
rz(-2.4251921) q[0];
rz(-2.7377103) q[1];
sx q[1];
rz(-1.5682181) q[1];
sx q[1];
rz(-1.4366038) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2540909) q[0];
sx q[0];
rz(-1.9483101) q[0];
sx q[0];
rz(0.15477263) q[0];
x q[1];
rz(0.5282497) q[2];
sx q[2];
rz(-1.6852323) q[2];
sx q[2];
rz(-2.6334762) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.75175367) q[1];
sx q[1];
rz(-1.1734278) q[1];
sx q[1];
rz(1.6786871) q[1];
x q[2];
rz(2.3327504) q[3];
sx q[3];
rz(-1.0029965) q[3];
sx q[3];
rz(1.1274991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7776103) q[2];
sx q[2];
rz(-0.97475514) q[2];
sx q[2];
rz(-3.0762365) q[2];
rz(-2.2743547) q[3];
sx q[3];
rz(-1.6403653) q[3];
sx q[3];
rz(-1.3245827) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48677483) q[0];
sx q[0];
rz(-1.7205394) q[0];
sx q[0];
rz(2.4105893) q[0];
rz(-0.77516088) q[1];
sx q[1];
rz(-2.8072) q[1];
sx q[1];
rz(2.7269272) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.054798445) q[0];
sx q[0];
rz(-2.2761184) q[0];
sx q[0];
rz(0.84673015) q[0];
rz(-pi) q[1];
rz(-2.6111905) q[2];
sx q[2];
rz(-1.3866405) q[2];
sx q[2];
rz(0.14094409) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9379256) q[1];
sx q[1];
rz(-1.3199184) q[1];
sx q[1];
rz(2.832042) q[1];
rz(-pi) q[2];
rz(2.25423) q[3];
sx q[3];
rz(-0.21036869) q[3];
sx q[3];
rz(0.31926814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9097462) q[2];
sx q[2];
rz(-0.48603386) q[2];
sx q[2];
rz(3.0492142) q[2];
rz(2.3653024) q[3];
sx q[3];
rz(-1.679136) q[3];
sx q[3];
rz(-2.9379454) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48731503) q[0];
sx q[0];
rz(-1.646811) q[0];
sx q[0];
rz(0.017729433) q[0];
rz(-1.0461294) q[1];
sx q[1];
rz(-1.7283864) q[1];
sx q[1];
rz(2.0808992) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26955206) q[0];
sx q[0];
rz(-1.5020796) q[0];
sx q[0];
rz(-1.5567354) q[0];
rz(-0.25814105) q[2];
sx q[2];
rz(-2.1434868) q[2];
sx q[2];
rz(1.0024459) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.880278) q[1];
sx q[1];
rz(-2.3572817) q[1];
sx q[1];
rz(1.8159291) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0396427) q[3];
sx q[3];
rz(-2.6155439) q[3];
sx q[3];
rz(-0.40204429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9866508) q[2];
sx q[2];
rz(-0.81601802) q[2];
sx q[2];
rz(-2.4972534) q[2];
rz(-2.2996357) q[3];
sx q[3];
rz(-1.5717477) q[3];
sx q[3];
rz(2.2346066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-1.7613206) q[0];
sx q[0];
rz(-0.43571061) q[0];
sx q[0];
rz(-2.6897588) q[0];
rz(1.0093581) q[1];
sx q[1];
rz(-1.511907) q[1];
sx q[1];
rz(-1.570328) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6752351) q[0];
sx q[0];
rz(-1.5646828) q[0];
sx q[0];
rz(0.17019043) q[0];
rz(-pi) q[1];
rz(-3.0250835) q[2];
sx q[2];
rz(-2.0429869) q[2];
sx q[2];
rz(-2.0906868) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.81262368) q[1];
sx q[1];
rz(-1.4419893) q[1];
sx q[1];
rz(1.179943) q[1];
rz(-pi) q[2];
rz(-1.5785923) q[3];
sx q[3];
rz(-2.2931794) q[3];
sx q[3];
rz(-2.9786199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5181804) q[2];
sx q[2];
rz(-2.0698915) q[2];
sx q[2];
rz(-1.5706459) q[2];
rz(3.0359641) q[3];
sx q[3];
rz(-0.94625866) q[3];
sx q[3];
rz(-2.6251729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8138206) q[0];
sx q[0];
rz(-2.1229424) q[0];
sx q[0];
rz(0.13737296) q[0];
rz(0.2188006) q[1];
sx q[1];
rz(-1.4004424) q[1];
sx q[1];
rz(-0.91011824) q[1];
rz(-1.995718) q[2];
sx q[2];
rz(-0.47987249) q[2];
sx q[2];
rz(2.1940827) q[2];
rz(-2.4841251) q[3];
sx q[3];
rz(-1.8651265) q[3];
sx q[3];
rz(-0.65262564) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
