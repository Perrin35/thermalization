OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6150317) q[0];
sx q[0];
rz(2.5685413) q[0];
sx q[0];
rz(11.723784) q[0];
rz(2.1057582) q[1];
sx q[1];
rz(-1.0993212) q[1];
sx q[1];
rz(1.458118) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6095088) q[0];
sx q[0];
rz(-2.3713787) q[0];
sx q[0];
rz(0.30755933) q[0];
x q[1];
rz(1.5904434) q[2];
sx q[2];
rz(-2.1845449) q[2];
sx q[2];
rz(-2.8710499) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2394489) q[1];
sx q[1];
rz(-1.4414756) q[1];
sx q[1];
rz(-0.37954482) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.17957844) q[3];
sx q[3];
rz(-0.60185963) q[3];
sx q[3];
rz(1.4048525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6136916) q[2];
sx q[2];
rz(-2.1353022) q[2];
sx q[2];
rz(-0.17949417) q[2];
rz(-1.2256631) q[3];
sx q[3];
rz(-1.7951199) q[3];
sx q[3];
rz(2.3195482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3935788) q[0];
sx q[0];
rz(-2.2606235) q[0];
sx q[0];
rz(-0.32546145) q[0];
rz(-1.7851967) q[1];
sx q[1];
rz(-1.0486832) q[1];
sx q[1];
rz(1.9869841) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1240631) q[0];
sx q[0];
rz(-1.5895827) q[0];
sx q[0];
rz(-1.5536874) q[0];
x q[1];
rz(1.0500533) q[2];
sx q[2];
rz(-0.69486952) q[2];
sx q[2];
rz(0.75887242) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6007538) q[1];
sx q[1];
rz(-2.3327017) q[1];
sx q[1];
rz(-1.4536588) q[1];
rz(-2.9566544) q[3];
sx q[3];
rz(-1.3122845) q[3];
sx q[3];
rz(1.0227433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4521728) q[2];
sx q[2];
rz(-1.2499115) q[2];
sx q[2];
rz(-0.88341218) q[2];
rz(2.6702821) q[3];
sx q[3];
rz(-1.4383957) q[3];
sx q[3];
rz(0.78770351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31323355) q[0];
sx q[0];
rz(-1.6468843) q[0];
sx q[0];
rz(-1.5154243) q[0];
rz(0.60107636) q[1];
sx q[1];
rz(-2.5939012) q[1];
sx q[1];
rz(1.0916969) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1900345) q[0];
sx q[0];
rz(-2.1265731) q[0];
sx q[0];
rz(0.65727289) q[0];
rz(-pi) q[1];
rz(1.0513564) q[2];
sx q[2];
rz(-0.72548496) q[2];
sx q[2];
rz(-1.5067593) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6669238) q[1];
sx q[1];
rz(-0.83819929) q[1];
sx q[1];
rz(1.0679507) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5063498) q[3];
sx q[3];
rz(-1.4430874) q[3];
sx q[3];
rz(2.8116022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8213356) q[2];
sx q[2];
rz(-0.50575033) q[2];
sx q[2];
rz(-2.2606405) q[2];
rz(-1.3736003) q[3];
sx q[3];
rz(-1.6146086) q[3];
sx q[3];
rz(1.0176456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83051935) q[0];
sx q[0];
rz(-1.7493462) q[0];
sx q[0];
rz(0.4367035) q[0];
rz(-0.23315915) q[1];
sx q[1];
rz(-1.2522839) q[1];
sx q[1];
rz(-2.8312347) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46582857) q[0];
sx q[0];
rz(-0.40256631) q[0];
sx q[0];
rz(2.7990544) q[0];
x q[1];
rz(-0.15375806) q[2];
sx q[2];
rz(-2.4506844) q[2];
sx q[2];
rz(1.5916057) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1102317) q[1];
sx q[1];
rz(-2.7831315) q[1];
sx q[1];
rz(1.6237753) q[1];
rz(1.2918606) q[3];
sx q[3];
rz(-3.012987) q[3];
sx q[3];
rz(2.0898553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.13005304) q[2];
sx q[2];
rz(-0.71802846) q[2];
sx q[2];
rz(2.0641573) q[2];
rz(0.056190101) q[3];
sx q[3];
rz(-0.63779938) q[3];
sx q[3];
rz(1.594054) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.189165) q[0];
sx q[0];
rz(-1.0452894) q[0];
sx q[0];
rz(-0.24965832) q[0];
rz(1.5646308) q[1];
sx q[1];
rz(-0.77762929) q[1];
sx q[1];
rz(-0.87019428) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72083005) q[0];
sx q[0];
rz(-1.2569068) q[0];
sx q[0];
rz(1.785196) q[0];
rz(-pi) q[1];
rz(-2.444961) q[2];
sx q[2];
rz(-1.4259035) q[2];
sx q[2];
rz(-0.87755132) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.40904564) q[1];
sx q[1];
rz(-2.2605719) q[1];
sx q[1];
rz(1.9802666) q[1];
rz(-pi) q[2];
rz(-1.3051885) q[3];
sx q[3];
rz(-2.5634273) q[3];
sx q[3];
rz(2.2418914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8683118) q[2];
sx q[2];
rz(-1.8153278) q[2];
sx q[2];
rz(2.4678521) q[2];
rz(-2.8379748) q[3];
sx q[3];
rz(-1.2250591) q[3];
sx q[3];
rz(1.3195066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7917787) q[0];
sx q[0];
rz(-2.204201) q[0];
sx q[0];
rz(-0.2579903) q[0];
rz(-2.7164283) q[1];
sx q[1];
rz(-0.95562569) q[1];
sx q[1];
rz(-1.649883) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30116044) q[0];
sx q[0];
rz(-1.668881) q[0];
sx q[0];
rz(-0.43099404) q[0];
x q[1];
rz(2.4620373) q[2];
sx q[2];
rz(-1.9050042) q[2];
sx q[2];
rz(-2.213775) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3867492) q[1];
sx q[1];
rz(-1.0068839) q[1];
sx q[1];
rz(-2.2350603) q[1];
x q[2];
rz(1.9353988) q[3];
sx q[3];
rz(-1.4758665) q[3];
sx q[3];
rz(1.7942384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1288746) q[2];
sx q[2];
rz(-2.1779163) q[2];
sx q[2];
rz(-3.0498665) q[2];
rz(-2.2979459) q[3];
sx q[3];
rz(-0.97976145) q[3];
sx q[3];
rz(0.89404026) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82350746) q[0];
sx q[0];
rz(-1.2473236) q[0];
sx q[0];
rz(-0.41123018) q[0];
rz(2.2757018) q[1];
sx q[1];
rz(-2.829268) q[1];
sx q[1];
rz(3.1076028) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5916409) q[0];
sx q[0];
rz(-2.3410428) q[0];
sx q[0];
rz(0.22649015) q[0];
x q[1];
rz(0.88862822) q[2];
sx q[2];
rz(-0.76105984) q[2];
sx q[2];
rz(1.6737446) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.781144) q[1];
sx q[1];
rz(-0.89372674) q[1];
sx q[1];
rz(-0.064868049) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9339318) q[3];
sx q[3];
rz(-2.5163979) q[3];
sx q[3];
rz(0.085426424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6035446) q[2];
sx q[2];
rz(-2.5431583) q[2];
sx q[2];
rz(-0.87654385) q[2];
rz(2.792568) q[3];
sx q[3];
rz(-1.2004431) q[3];
sx q[3];
rz(-0.14311895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8975163) q[0];
sx q[0];
rz(-1.4325457) q[0];
sx q[0];
rz(-2.7476655) q[0];
rz(0.36755964) q[1];
sx q[1];
rz(-1.7575248) q[1];
sx q[1];
rz(-1.4454909) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5515585) q[0];
sx q[0];
rz(-2.708205) q[0];
sx q[0];
rz(-1.5831468) q[0];
rz(-pi) q[1];
rz(2.3169575) q[2];
sx q[2];
rz(-2.3687009) q[2];
sx q[2];
rz(0.075721272) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6519424) q[1];
sx q[1];
rz(-1.7282657) q[1];
sx q[1];
rz(0.55649346) q[1];
rz(-pi) q[2];
rz(-0.74650713) q[3];
sx q[3];
rz(-0.81614796) q[3];
sx q[3];
rz(1.5619123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2945071) q[2];
sx q[2];
rz(-2.2512348) q[2];
sx q[2];
rz(2.7344446) q[2];
rz(1.6242705) q[3];
sx q[3];
rz(-1.1573236) q[3];
sx q[3];
rz(2.8919162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80609926) q[0];
sx q[0];
rz(-0.51500106) q[0];
sx q[0];
rz(1.8898213) q[0];
rz(-2.4720526) q[1];
sx q[1];
rz(-1.1839097) q[1];
sx q[1];
rz(0.30977419) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71000242) q[0];
sx q[0];
rz(-2.6012523) q[0];
sx q[0];
rz(-2.8938328) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3938815) q[2];
sx q[2];
rz(-1.5614911) q[2];
sx q[2];
rz(-2.6975346) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.138315) q[1];
sx q[1];
rz(-1.365108) q[1];
sx q[1];
rz(0.059271952) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5449949) q[3];
sx q[3];
rz(-1.784424) q[3];
sx q[3];
rz(-1.9006157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.70242515) q[2];
sx q[2];
rz(-2.4283786) q[2];
sx q[2];
rz(-1.9343728) q[2];
rz(1.0369982) q[3];
sx q[3];
rz(-1.2456649) q[3];
sx q[3];
rz(0.65565482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.050215125) q[0];
sx q[0];
rz(-1.8176879) q[0];
sx q[0];
rz(1.9357095) q[0];
rz(0.58569113) q[1];
sx q[1];
rz(-2.060545) q[1];
sx q[1];
rz(1.4996128) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9046281) q[0];
sx q[0];
rz(-1.4693854) q[0];
sx q[0];
rz(-1.8629835) q[0];
x q[1];
rz(-1.502938) q[2];
sx q[2];
rz(-2.3336764) q[2];
sx q[2];
rz(0.18005904) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3824532) q[1];
sx q[1];
rz(-1.1083974) q[1];
sx q[1];
rz(-2.9049302) q[1];
rz(-1.5564735) q[3];
sx q[3];
rz(-2.9687772) q[3];
sx q[3];
rz(-0.9848136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2575834) q[2];
sx q[2];
rz(-1.7929701) q[2];
sx q[2];
rz(2.1949027) q[2];
rz(0.36869129) q[3];
sx q[3];
rz(-1.576141) q[3];
sx q[3];
rz(-0.45599109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35836999) q[0];
sx q[0];
rz(-1.9932278) q[0];
sx q[0];
rz(2.7182462) q[0];
rz(3.070667) q[1];
sx q[1];
rz(-1.4535041) q[1];
sx q[1];
rz(2.8765875) q[1];
rz(0.46438607) q[2];
sx q[2];
rz(-1.1190363) q[2];
sx q[2];
rz(-2.6418532) q[2];
rz(0.72017097) q[3];
sx q[3];
rz(-1.9577033) q[3];
sx q[3];
rz(-2.4693558) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];