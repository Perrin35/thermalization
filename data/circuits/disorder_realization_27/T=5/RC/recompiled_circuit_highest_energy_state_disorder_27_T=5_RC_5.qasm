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
rz(1.2394387) q[0];
sx q[0];
rz(-1.8128938) q[0];
sx q[0];
rz(-0.21188307) q[0];
rz(-1.4172685) q[1];
sx q[1];
rz(-2.6091726) q[1];
sx q[1];
rz(0.37791696) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53055313) q[0];
sx q[0];
rz(-1.5612649) q[0];
sx q[0];
rz(3.1411489) q[0];
rz(-1.9340408) q[2];
sx q[2];
rz(-2.2170728) q[2];
sx q[2];
rz(-2.5763047) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7211044) q[1];
sx q[1];
rz(-1.2961565) q[1];
sx q[1];
rz(0.65914776) q[1];
rz(-pi) q[2];
rz(-0.7573646) q[3];
sx q[3];
rz(-0.93776449) q[3];
sx q[3];
rz(-2.4773134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8523031) q[2];
sx q[2];
rz(-1.0785582) q[2];
sx q[2];
rz(0.079785384) q[2];
rz(0.96528178) q[3];
sx q[3];
rz(-1.691317) q[3];
sx q[3];
rz(-2.4108346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-0.47736436) q[0];
sx q[0];
rz(-2.1381162) q[0];
sx q[0];
rz(3.0526414) q[0];
rz(-1.7695919) q[1];
sx q[1];
rz(-1.7141432) q[1];
sx q[1];
rz(-2.037183) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0029582214) q[0];
sx q[0];
rz(-1.3413652) q[0];
sx q[0];
rz(0.63284875) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.78271336) q[2];
sx q[2];
rz(-1.209895) q[2];
sx q[2];
rz(-0.95907839) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.6767533) q[1];
sx q[1];
rz(-1.1285121) q[1];
sx q[1];
rz(-0.16525903) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6703963) q[3];
sx q[3];
rz(-2.1079113) q[3];
sx q[3];
rz(0.84505862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.264512) q[2];
sx q[2];
rz(-1.1254213) q[2];
sx q[2];
rz(-1.6748927) q[2];
rz(3.1125715) q[3];
sx q[3];
rz(-2.1252188) q[3];
sx q[3];
rz(-0.69028729) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90917176) q[0];
sx q[0];
rz(-3.0747774) q[0];
sx q[0];
rz(-0.27798852) q[0];
rz(-1.7104644) q[1];
sx q[1];
rz(-2.2093096) q[1];
sx q[1];
rz(-2.7412282) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58996449) q[0];
sx q[0];
rz(-2.0144406) q[0];
sx q[0];
rz(-0.69728627) q[0];
x q[1];
rz(0.30817356) q[2];
sx q[2];
rz(-0.71012596) q[2];
sx q[2];
rz(-1.9911204) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.4246042) q[1];
sx q[1];
rz(-1.4115872) q[1];
sx q[1];
rz(-0.47674322) q[1];
rz(1.2648029) q[3];
sx q[3];
rz(-1.8974432) q[3];
sx q[3];
rz(-1.5375002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4572738) q[2];
sx q[2];
rz(-1.6088586) q[2];
sx q[2];
rz(2.686783) q[2];
rz(1.2416035) q[3];
sx q[3];
rz(-0.95507115) q[3];
sx q[3];
rz(0.72030592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95551816) q[0];
sx q[0];
rz(-1.3294514) q[0];
sx q[0];
rz(-3.1410826) q[0];
rz(2.5406802) q[1];
sx q[1];
rz(-2.3020703) q[1];
sx q[1];
rz(2.9972163) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7847608) q[0];
sx q[0];
rz(-0.5934754) q[0];
sx q[0];
rz(-1.9463825) q[0];
x q[1];
rz(-2.8881489) q[2];
sx q[2];
rz(-1.7161233) q[2];
sx q[2];
rz(3.1325454) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6295602) q[1];
sx q[1];
rz(-1.7981537) q[1];
sx q[1];
rz(0.65452202) q[1];
rz(1.6779283) q[3];
sx q[3];
rz(-2.1153573) q[3];
sx q[3];
rz(0.0060826172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.942975) q[2];
sx q[2];
rz(-1.696442) q[2];
sx q[2];
rz(0.96251881) q[2];
rz(1.5396384) q[3];
sx q[3];
rz(-1.736015) q[3];
sx q[3];
rz(0.34077728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9050423) q[0];
sx q[0];
rz(-1.4555229) q[0];
sx q[0];
rz(-1.1849674) q[0];
rz(-2.9153337) q[1];
sx q[1];
rz(-0.87892756) q[1];
sx q[1];
rz(1.0386946) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5796367) q[0];
sx q[0];
rz(-2.547894) q[0];
sx q[0];
rz(0.14621347) q[0];
rz(-pi) q[1];
rz(0.67579999) q[2];
sx q[2];
rz(-0.86188176) q[2];
sx q[2];
rz(-1.0647286) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1490399) q[1];
sx q[1];
rz(-2.2254308) q[1];
sx q[1];
rz(-0.66019411) q[1];
rz(-pi) q[2];
rz(1.709278) q[3];
sx q[3];
rz(-0.44381986) q[3];
sx q[3];
rz(-2.6276772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7096536) q[2];
sx q[2];
rz(-0.69810549) q[2];
sx q[2];
rz(2.8254438) q[2];
rz(-1.6759253) q[3];
sx q[3];
rz(-2.0114055) q[3];
sx q[3];
rz(1.7620311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6370711) q[0];
sx q[0];
rz(-0.9032473) q[0];
sx q[0];
rz(-1.1035408) q[0];
rz(0.48031131) q[1];
sx q[1];
rz(-0.66245285) q[1];
sx q[1];
rz(0.59741098) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7561853) q[0];
sx q[0];
rz(-2.8592337) q[0];
sx q[0];
rz(0.83448164) q[0];
x q[1];
rz(-1.811932) q[2];
sx q[2];
rz(-0.84000194) q[2];
sx q[2];
rz(-2.3868449) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.589041) q[1];
sx q[1];
rz(-1.8218166) q[1];
sx q[1];
rz(1.7415206) q[1];
rz(0.091986309) q[3];
sx q[3];
rz(-1.2013271) q[3];
sx q[3];
rz(1.8487768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.20208134) q[2];
sx q[2];
rz(-1.6263522) q[2];
sx q[2];
rz(2.0651979) q[2];
rz(1.9988029) q[3];
sx q[3];
rz(-0.79311526) q[3];
sx q[3];
rz(-2.6514261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48653212) q[0];
sx q[0];
rz(-1.9831816) q[0];
sx q[0];
rz(3.1031188) q[0];
rz(-0.064727457) q[1];
sx q[1];
rz(-1.3836626) q[1];
sx q[1];
rz(0.23385349) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.192894) q[0];
sx q[0];
rz(-1.7561551) q[0];
sx q[0];
rz(-2.9167487) q[0];
rz(0.87217561) q[2];
sx q[2];
rz(-2.6769014) q[2];
sx q[2];
rz(-1.8605491) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.55091399) q[1];
sx q[1];
rz(-1.4812902) q[1];
sx q[1];
rz(2.4924335) q[1];
rz(-pi) q[2];
rz(3.0117118) q[3];
sx q[3];
rz(-1.6280773) q[3];
sx q[3];
rz(-2.2735325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.48876277) q[2];
sx q[2];
rz(-1.5831999) q[2];
sx q[2];
rz(0.59824198) q[2];
rz(-0.11387842) q[3];
sx q[3];
rz(-1.3860393) q[3];
sx q[3];
rz(-2.2940476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4050196) q[0];
sx q[0];
rz(-0.84252715) q[0];
sx q[0];
rz(-0.2463499) q[0];
rz(1.8440638) q[1];
sx q[1];
rz(-2.0228701) q[1];
sx q[1];
rz(-0.49682239) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56820541) q[0];
sx q[0];
rz(-2.3620689) q[0];
sx q[0];
rz(0.60978344) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.70513983) q[2];
sx q[2];
rz(-1.6821096) q[2];
sx q[2];
rz(1.706858) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0441974) q[1];
sx q[1];
rz(-0.5438416) q[1];
sx q[1];
rz(2.0378276) q[1];
x q[2];
rz(-0.15487352) q[3];
sx q[3];
rz(-1.0836156) q[3];
sx q[3];
rz(-2.7497512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7184427) q[2];
sx q[2];
rz(-0.53360525) q[2];
sx q[2];
rz(-1.144484) q[2];
rz(-1.1540958) q[3];
sx q[3];
rz(-1.7172979) q[3];
sx q[3];
rz(0.19788876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4474739) q[0];
sx q[0];
rz(-2.9070774) q[0];
sx q[0];
rz(3.0754454) q[0];
rz(-1.1478395) q[1];
sx q[1];
rz(-1.3775974) q[1];
sx q[1];
rz(-2.537421) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87245524) q[0];
sx q[0];
rz(-2.7949484) q[0];
sx q[0];
rz(-0.64328648) q[0];
rz(-pi) q[1];
rz(-1.7275024) q[2];
sx q[2];
rz(-1.0494266) q[2];
sx q[2];
rz(-0.40796134) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7614038) q[1];
sx q[1];
rz(-2.1288925) q[1];
sx q[1];
rz(1.7978884) q[1];
rz(-pi) q[2];
x q[2];
rz(0.6571397) q[3];
sx q[3];
rz(-0.98987386) q[3];
sx q[3];
rz(-2.7681153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.26312795) q[2];
sx q[2];
rz(-2.2915514) q[2];
sx q[2];
rz(-0.43295941) q[2];
rz(1.912502) q[3];
sx q[3];
rz(-1.9357598) q[3];
sx q[3];
rz(-1.3273201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94806725) q[0];
sx q[0];
rz(-1.075241) q[0];
sx q[0];
rz(-1.8684335) q[0];
rz(-0.46514568) q[1];
sx q[1];
rz(-1.3659313) q[1];
sx q[1];
rz(-0.21496162) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5453542) q[0];
sx q[0];
rz(-2.4430877) q[0];
sx q[0];
rz(-0.39909382) q[0];
rz(-pi) q[1];
x q[1];
rz(0.5318435) q[2];
sx q[2];
rz(-1.1924679) q[2];
sx q[2];
rz(1.7500306) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7219639) q[1];
sx q[1];
rz(-0.89744324) q[1];
sx q[1];
rz(0.90760214) q[1];
rz(1.3137903) q[3];
sx q[3];
rz(-2.1177835) q[3];
sx q[3];
rz(-1.6099324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.8421858) q[2];
sx q[2];
rz(-2.3278548) q[2];
sx q[2];
rz(2.1827533) q[2];
rz(-1.9793319) q[3];
sx q[3];
rz(-1.4872888) q[3];
sx q[3];
rz(-1.6963262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6017629) q[0];
sx q[0];
rz(-1.5400664) q[0];
sx q[0];
rz(-1.6590317) q[0];
rz(0.70855793) q[1];
sx q[1];
rz(-2.9500912) q[1];
sx q[1];
rz(-0.74831829) q[1];
rz(2.5790527) q[2];
sx q[2];
rz(-1.2934791) q[2];
sx q[2];
rz(0.51359609) q[2];
rz(2.419653) q[3];
sx q[3];
rz(-1.7122713) q[3];
sx q[3];
rz(0.28817978) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
