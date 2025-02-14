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
rz(2.7163765) q[0];
sx q[0];
rz(-1.334231) q[0];
sx q[0];
rz(0.38036007) q[0];
rz(1.7827787) q[1];
sx q[1];
rz(-0.18667297) q[1];
sx q[1];
rz(0.92460728) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0565489) q[0];
sx q[0];
rz(-2.2593772) q[0];
sx q[0];
rz(-2.0184645) q[0];
rz(1.0546646) q[2];
sx q[2];
rz(-1.4251815) q[2];
sx q[2];
rz(-2.9235385) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2362931) q[1];
sx q[1];
rz(-0.84262203) q[1];
sx q[1];
rz(0.052680101) q[1];
rz(-3.0136631) q[3];
sx q[3];
rz(-1.5618213) q[3];
sx q[3];
rz(-2.8913446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1251462) q[2];
sx q[2];
rz(-2.1852198) q[2];
sx q[2];
rz(-2.1535786) q[2];
rz(-0.2068578) q[3];
sx q[3];
rz(-1.7346953) q[3];
sx q[3];
rz(0.44001165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67966953) q[0];
sx q[0];
rz(-1.672687) q[0];
sx q[0];
rz(0.2521421) q[0];
rz(-3.120046) q[1];
sx q[1];
rz(-2.4485059) q[1];
sx q[1];
rz(0.85539877) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48777521) q[0];
sx q[0];
rz(-3.1145373) q[0];
sx q[0];
rz(-1.6271126) q[0];
x q[1];
rz(-0.55412519) q[2];
sx q[2];
rz(-1.9337855) q[2];
sx q[2];
rz(-2.3100694) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.508042) q[1];
sx q[1];
rz(-1.3505858) q[1];
sx q[1];
rz(-2.1160265) q[1];
rz(0.40327252) q[3];
sx q[3];
rz(-0.72411116) q[3];
sx q[3];
rz(1.6722752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1713193) q[2];
sx q[2];
rz(-0.0042985175) q[2];
sx q[2];
rz(-2.8805736) q[2];
rz(0.039693443) q[3];
sx q[3];
rz(-1.3986162) q[3];
sx q[3];
rz(-0.54350129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4259341) q[0];
sx q[0];
rz(-1.4733227) q[0];
sx q[0];
rz(0.69450992) q[0];
rz(-2.014324) q[1];
sx q[1];
rz(-0.85113168) q[1];
sx q[1];
rz(0.99004254) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8111525) q[0];
sx q[0];
rz(-1.2227204) q[0];
sx q[0];
rz(-1.6578107) q[0];
x q[1];
rz(2.9431413) q[2];
sx q[2];
rz(-2.1439512) q[2];
sx q[2];
rz(-2.1015374) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7496365) q[1];
sx q[1];
rz(-0.5881433) q[1];
sx q[1];
rz(-0.68444499) q[1];
rz(0.95616266) q[3];
sx q[3];
rz(-0.80397881) q[3];
sx q[3];
rz(-1.965919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8433044) q[2];
sx q[2];
rz(-0.97323209) q[2];
sx q[2];
rz(2.6514371) q[2];
rz(-2.7847024) q[3];
sx q[3];
rz(-2.7481952) q[3];
sx q[3];
rz(1.2662158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.056034293) q[0];
sx q[0];
rz(-1.1857251) q[0];
sx q[0];
rz(0.84247843) q[0];
rz(0.8017686) q[1];
sx q[1];
rz(-2.8623878) q[1];
sx q[1];
rz(-2.3559949) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7467037) q[0];
sx q[0];
rz(-0.75468844) q[0];
sx q[0];
rz(-2.7695023) q[0];
rz(-pi) q[1];
rz(-2.9688492) q[2];
sx q[2];
rz(-0.72451353) q[2];
sx q[2];
rz(-0.24562626) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.468049) q[1];
sx q[1];
rz(-1.0338514) q[1];
sx q[1];
rz(-1.931744) q[1];
x q[2];
rz(0.20468851) q[3];
sx q[3];
rz(-0.83016268) q[3];
sx q[3];
rz(1.9616753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.12513146) q[2];
sx q[2];
rz(-0.74840122) q[2];
sx q[2];
rz(-1.008519) q[2];
rz(-2.290001) q[3];
sx q[3];
rz(-2.3840756) q[3];
sx q[3];
rz(1.0803224) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1763024) q[0];
sx q[0];
rz(-2.1602614) q[0];
sx q[0];
rz(1.0705795) q[0];
rz(2.3640682) q[1];
sx q[1];
rz(-0.91064015) q[1];
sx q[1];
rz(-2.1308965) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1436737) q[0];
sx q[0];
rz(-2.5287313) q[0];
sx q[0];
rz(-2.4852703) q[0];
rz(-1.7456876) q[2];
sx q[2];
rz(-1.5064459) q[2];
sx q[2];
rz(-1.5140057) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.7091211) q[1];
sx q[1];
rz(-0.84408954) q[1];
sx q[1];
rz(2.5120887) q[1];
rz(-0.56212496) q[3];
sx q[3];
rz(-1.6549806) q[3];
sx q[3];
rz(-0.69417324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8396478) q[2];
sx q[2];
rz(-3.0366615) q[2];
sx q[2];
rz(-1.5860175) q[2];
rz(1.0157061) q[3];
sx q[3];
rz(-1.1414707) q[3];
sx q[3];
rz(1.4122081) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79065901) q[0];
sx q[0];
rz(-0.18336329) q[0];
sx q[0];
rz(0.80107981) q[0];
rz(0.13310295) q[1];
sx q[1];
rz(-1.3214654) q[1];
sx q[1];
rz(0.4020234) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6989845) q[0];
sx q[0];
rz(-0.8359209) q[0];
sx q[0];
rz(1.1135654) q[0];
rz(-1.5636251) q[2];
sx q[2];
rz(-2.8507887) q[2];
sx q[2];
rz(-0.70722843) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.61643663) q[1];
sx q[1];
rz(-2.0503938) q[1];
sx q[1];
rz(-2.8308771) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4834255) q[3];
sx q[3];
rz(-1.6699381) q[3];
sx q[3];
rz(1.0189354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.63783995) q[2];
sx q[2];
rz(-1.736234) q[2];
sx q[2];
rz(-2.3657738) q[2];
rz(1.4555629) q[3];
sx q[3];
rz(-1.9676696) q[3];
sx q[3];
rz(2.1337401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.044947226) q[0];
sx q[0];
rz(-2.2438887) q[0];
sx q[0];
rz(-2.4626379) q[0];
rz(-1.2848162) q[1];
sx q[1];
rz(-0.7130475) q[1];
sx q[1];
rz(-1.3791893) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5011936) q[0];
sx q[0];
rz(-1.7441347) q[0];
sx q[0];
rz(0.79357432) q[0];
x q[1];
rz(3.0581362) q[2];
sx q[2];
rz(-1.8493358) q[2];
sx q[2];
rz(-0.73501529) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0028982) q[1];
sx q[1];
rz(-1.0779253) q[1];
sx q[1];
rz(-1.0652131) q[1];
rz(-pi) q[2];
rz(-1.3016765) q[3];
sx q[3];
rz(-2.9342071) q[3];
sx q[3];
rz(0.086750448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3369559) q[2];
sx q[2];
rz(-0.75504428) q[2];
sx q[2];
rz(1.9026559) q[2];
rz(1.8094481) q[3];
sx q[3];
rz(-2.381031) q[3];
sx q[3];
rz(-2.8968887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7234583) q[0];
sx q[0];
rz(-1.9754388) q[0];
sx q[0];
rz(0.13846692) q[0];
rz(-2.9578517) q[1];
sx q[1];
rz(-0.67444363) q[1];
sx q[1];
rz(-0.26652452) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4999102) q[0];
sx q[0];
rz(-1.2453834) q[0];
sx q[0];
rz(-0.17123674) q[0];
rz(-pi) q[1];
rz(-1.0037854) q[2];
sx q[2];
rz(-1.8006174) q[2];
sx q[2];
rz(2.9858495) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9723818) q[1];
sx q[1];
rz(-1.5145497) q[1];
sx q[1];
rz(-1.1509658) q[1];
rz(0.72690225) q[3];
sx q[3];
rz(-1.57734) q[3];
sx q[3];
rz(0.18101276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0515685) q[2];
sx q[2];
rz(-1.2568306) q[2];
sx q[2];
rz(0.23452342) q[2];
rz(1.6656434) q[3];
sx q[3];
rz(-1.539295) q[3];
sx q[3];
rz(0.75638151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.804857) q[0];
sx q[0];
rz(-2.2754301) q[0];
sx q[0];
rz(0.79972237) q[0];
rz(2.5526478) q[1];
sx q[1];
rz(-0.84732333) q[1];
sx q[1];
rz(-2.934093) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.070415592) q[0];
sx q[0];
rz(-1.9635597) q[0];
sx q[0];
rz(-1.7123187) q[0];
rz(2.2830711) q[2];
sx q[2];
rz(-1.3337161) q[2];
sx q[2];
rz(1.3905987) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7991528) q[1];
sx q[1];
rz(-0.9121597) q[1];
sx q[1];
rz(2.6299146) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6274243) q[3];
sx q[3];
rz(-0.47041962) q[3];
sx q[3];
rz(-0.75477598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3346682) q[2];
sx q[2];
rz(-2.25756) q[2];
sx q[2];
rz(-0.36726382) q[2];
rz(1.4372829) q[3];
sx q[3];
rz(-1.9673012) q[3];
sx q[3];
rz(-1.9678763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.2631898) q[0];
sx q[0];
rz(-1.0028361) q[0];
sx q[0];
rz(-2.3314085) q[0];
rz(2.0267678) q[1];
sx q[1];
rz(-0.38433847) q[1];
sx q[1];
rz(-0.55327639) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1076521) q[0];
sx q[0];
rz(-2.2080511) q[0];
sx q[0];
rz(-3.0725543) q[0];
rz(-pi) q[1];
rz(2.3478339) q[2];
sx q[2];
rz(-1.1053893) q[2];
sx q[2];
rz(-0.88875729) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3590433) q[1];
sx q[1];
rz(-1.3842674) q[1];
sx q[1];
rz(-1.2839509) q[1];
rz(2.7891879) q[3];
sx q[3];
rz(-2.3239072) q[3];
sx q[3];
rz(1.529983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1222003) q[2];
sx q[2];
rz(-0.70851749) q[2];
sx q[2];
rz(2.9019287) q[2];
rz(1.3278809) q[3];
sx q[3];
rz(-2.0769104) q[3];
sx q[3];
rz(0.46752587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0620621) q[0];
sx q[0];
rz(-1.7876328) q[0];
sx q[0];
rz(-2.7476516) q[0];
rz(2.486034) q[1];
sx q[1];
rz(-1.9656904) q[1];
sx q[1];
rz(0.94737731) q[1];
rz(2.4782933) q[2];
sx q[2];
rz(-2.6313783) q[2];
sx q[2];
rz(0.22102697) q[2];
rz(-1.3169133) q[3];
sx q[3];
rz(-2.1364273) q[3];
sx q[3];
rz(1.9404503) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
