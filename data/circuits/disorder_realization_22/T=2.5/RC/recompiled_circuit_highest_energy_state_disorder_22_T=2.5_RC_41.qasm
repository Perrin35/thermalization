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
rz(0.15777388) q[0];
sx q[0];
rz(-1.1717492) q[0];
sx q[0];
rz(-1.8921312) q[0];
rz(-0.14021048) q[1];
sx q[1];
rz(-1.6970716) q[1];
sx q[1];
rz(0.061847774) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0842347) q[0];
sx q[0];
rz(-1.0751197) q[0];
sx q[0];
rz(2.2593014) q[0];
rz(-1.5711478) q[2];
sx q[2];
rz(-1.4803737) q[2];
sx q[2];
rz(3.0302559) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1143519) q[1];
sx q[1];
rz(-1.6682965) q[1];
sx q[1];
rz(2.2932294) q[1];
rz(-pi) q[2];
rz(-0.88229813) q[3];
sx q[3];
rz(-1.8940482) q[3];
sx q[3];
rz(2.0417449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.42741117) q[2];
sx q[2];
rz(-2.2990871) q[2];
sx q[2];
rz(-2.8851435) q[2];
rz(-2.9668258) q[3];
sx q[3];
rz(-1.1656961) q[3];
sx q[3];
rz(-0.83785653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3278219) q[0];
sx q[0];
rz(-0.34532663) q[0];
sx q[0];
rz(3.0178965) q[0];
rz(-0.23873121) q[1];
sx q[1];
rz(-0.78014603) q[1];
sx q[1];
rz(-0.10496584) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86228131) q[0];
sx q[0];
rz(-2.263592) q[0];
sx q[0];
rz(-2.9458118) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3388073) q[2];
sx q[2];
rz(-1.9769443) q[2];
sx q[2];
rz(0.37877235) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.87890307) q[1];
sx q[1];
rz(-0.57003747) q[1];
sx q[1];
rz(-1.0703342) q[1];
rz(-0.70521783) q[3];
sx q[3];
rz(-1.2949484) q[3];
sx q[3];
rz(2.5666756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3671942) q[2];
sx q[2];
rz(-1.6509193) q[2];
sx q[2];
rz(-1.4614089) q[2];
rz(1.2857619) q[3];
sx q[3];
rz(-1.4444084) q[3];
sx q[3];
rz(2.8694966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1478145) q[0];
sx q[0];
rz(-2.9312134) q[0];
sx q[0];
rz(-1.0149957) q[0];
rz(2.2442832) q[1];
sx q[1];
rz(-1.4941447) q[1];
sx q[1];
rz(-2.4651333) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7583026) q[0];
sx q[0];
rz(-1.5448017) q[0];
sx q[0];
rz(-0.95375188) q[0];
x q[1];
rz(-0.49336596) q[2];
sx q[2];
rz(-2.1378433) q[2];
sx q[2];
rz(2.9285448) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.23692098) q[1];
sx q[1];
rz(-2.7075012) q[1];
sx q[1];
rz(2.0660115) q[1];
rz(-pi) q[2];
rz(-3.1234497) q[3];
sx q[3];
rz(-2.3251495) q[3];
sx q[3];
rz(-1.7545561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.93495381) q[2];
sx q[2];
rz(-2.166344) q[2];
sx q[2];
rz(0.76033956) q[2];
rz(1.9350516) q[3];
sx q[3];
rz(-2.3307266) q[3];
sx q[3];
rz(0.84347239) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7059785) q[0];
sx q[0];
rz(-2.5137081) q[0];
sx q[0];
rz(-1.5950369) q[0];
rz(-0.71890038) q[1];
sx q[1];
rz(-1.3201821) q[1];
sx q[1];
rz(-2.9100606) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.084615413) q[0];
sx q[0];
rz(-1.1418793) q[0];
sx q[0];
rz(2.5713628) q[0];
rz(-1.5668554) q[2];
sx q[2];
rz(-2.3634499) q[2];
sx q[2];
rz(-0.63150327) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0445679) q[1];
sx q[1];
rz(-2.6480354) q[1];
sx q[1];
rz(0.19792168) q[1];
rz(-pi) q[2];
rz(2.5367686) q[3];
sx q[3];
rz(-1.2422274) q[3];
sx q[3];
rz(2.891345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4270758) q[2];
sx q[2];
rz(-0.38893739) q[2];
sx q[2];
rz(-2.5330949) q[2];
rz(-2.9454339) q[3];
sx q[3];
rz(-1.720263) q[3];
sx q[3];
rz(0.99854809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(1.8683559) q[0];
sx q[0];
rz(-2.4899857) q[0];
sx q[0];
rz(0.77907816) q[0];
rz(0.92998663) q[1];
sx q[1];
rz(-1.9735186) q[1];
sx q[1];
rz(1.9680061) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5835607) q[0];
sx q[0];
rz(-0.3380709) q[0];
sx q[0];
rz(-1.5826691) q[0];
rz(1.3814244) q[2];
sx q[2];
rz(-2.5177285) q[2];
sx q[2];
rz(0.80846918) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8820937) q[1];
sx q[1];
rz(-0.32392392) q[1];
sx q[1];
rz(1.7739576) q[1];
x q[2];
rz(0.39802528) q[3];
sx q[3];
rz(-0.46062352) q[3];
sx q[3];
rz(1.6470136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8194627) q[2];
sx q[2];
rz(-2.9746015) q[2];
sx q[2];
rz(0.67506153) q[2];
rz(2.2760462) q[3];
sx q[3];
rz(-1.5138488) q[3];
sx q[3];
rz(1.9338098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.753767) q[0];
sx q[0];
rz(-0.017711552) q[0];
sx q[0];
rz(-0.87131635) q[0];
rz(-2.9246092) q[1];
sx q[1];
rz(-1.5996409) q[1];
sx q[1];
rz(-1.3380231) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17413134) q[0];
sx q[0];
rz(-0.92087692) q[0];
sx q[0];
rz(1.4667257) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7291405) q[2];
sx q[2];
rz(-2.4096903) q[2];
sx q[2];
rz(-2.5661693) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2830401) q[1];
sx q[1];
rz(-1.6190908) q[1];
sx q[1];
rz(1.4273391) q[1];
x q[2];
rz(-0.48465259) q[3];
sx q[3];
rz(-1.7395057) q[3];
sx q[3];
rz(1.3107306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.013082144) q[2];
sx q[2];
rz(-1.4980114) q[2];
sx q[2];
rz(0.80005542) q[2];
rz(2.1498146) q[3];
sx q[3];
rz(-2.1938775) q[3];
sx q[3];
rz(-0.94314027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8868788) q[0];
sx q[0];
rz(-0.44097057) q[0];
sx q[0];
rz(2.9898341) q[0];
rz(-1.4166547) q[1];
sx q[1];
rz(-2.2817426) q[1];
sx q[1];
rz(-2.3053665) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72513103) q[0];
sx q[0];
rz(-0.90354474) q[0];
sx q[0];
rz(-0.50531549) q[0];
x q[1];
rz(0.012846666) q[2];
sx q[2];
rz(-0.57672665) q[2];
sx q[2];
rz(1.6614514) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6265833) q[1];
sx q[1];
rz(-1.3224416) q[1];
sx q[1];
rz(3.0529313) q[1];
x q[2];
rz(-1.7015905) q[3];
sx q[3];
rz(-2.0328641) q[3];
sx q[3];
rz(-0.82926428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4862711) q[2];
sx q[2];
rz(-3.0781015) q[2];
sx q[2];
rz(-0.07902321) q[2];
rz(2.4231353) q[3];
sx q[3];
rz(-1.6301165) q[3];
sx q[3];
rz(2.2169936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(2.786161) q[0];
sx q[0];
rz(-0.59816718) q[0];
sx q[0];
rz(-2.2736736) q[0];
rz(-0.68880853) q[1];
sx q[1];
rz(-0.30615607) q[1];
sx q[1];
rz(2.205663) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3132738) q[0];
sx q[0];
rz(-2.5848645) q[0];
sx q[0];
rz(2.7931045) q[0];
x q[1];
rz(1.9804968) q[2];
sx q[2];
rz(-1.2563224) q[2];
sx q[2];
rz(-2.1358228) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3080146) q[1];
sx q[1];
rz(-1.2460099) q[1];
sx q[1];
rz(-3.0851425) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3981188) q[3];
sx q[3];
rz(-1.5718565) q[3];
sx q[3];
rz(-0.27976945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9800637) q[2];
sx q[2];
rz(-0.65902013) q[2];
sx q[2];
rz(2.8847983) q[2];
rz(2.1234546) q[3];
sx q[3];
rz(-1.1192106) q[3];
sx q[3];
rz(-0.97091466) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37860206) q[0];
sx q[0];
rz(-2.584223) q[0];
sx q[0];
rz(0.85365224) q[0];
rz(-1.254982) q[1];
sx q[1];
rz(-1.1458784) q[1];
sx q[1];
rz(-0.34057158) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4494394) q[0];
sx q[0];
rz(-1.7652349) q[0];
sx q[0];
rz(2.9783335) q[0];
rz(-pi) q[1];
rz(2.8355424) q[2];
sx q[2];
rz(-0.69717625) q[2];
sx q[2];
rz(-2.9070284) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.53696991) q[1];
sx q[1];
rz(-2.3019362) q[1];
sx q[1];
rz(2.5174058) q[1];
x q[2];
rz(-0.27517516) q[3];
sx q[3];
rz(-2.0180297) q[3];
sx q[3];
rz(2.3768611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.42369947) q[2];
sx q[2];
rz(-1.7607949) q[2];
sx q[2];
rz(2.763486) q[2];
rz(2.9041491) q[3];
sx q[3];
rz(-1.0708555) q[3];
sx q[3];
rz(2.8606991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.045192748) q[0];
sx q[0];
rz(-2.0572331) q[0];
sx q[0];
rz(0.64724809) q[0];
rz(-1.9735533) q[1];
sx q[1];
rz(-2.2868575) q[1];
sx q[1];
rz(-2.7580269) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76825324) q[0];
sx q[0];
rz(-3.1041234) q[0];
sx q[0];
rz(2.5277407) q[0];
rz(-pi) q[1];
rz(1.2403383) q[2];
sx q[2];
rz(-1.0631732) q[2];
sx q[2];
rz(-0.79703813) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3473222) q[1];
sx q[1];
rz(-1.8350661) q[1];
sx q[1];
rz(1.7842954) q[1];
rz(-1.3829154) q[3];
sx q[3];
rz(-0.5231103) q[3];
sx q[3];
rz(-0.95978242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8770404) q[2];
sx q[2];
rz(-0.88901797) q[2];
sx q[2];
rz(-1.5846579) q[2];
rz(1.5133739) q[3];
sx q[3];
rz(-1.7729365) q[3];
sx q[3];
rz(-2.7148066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.3601111) q[0];
sx q[0];
rz(-1.4327015) q[0];
sx q[0];
rz(-2.844368) q[0];
rz(1.9602641) q[1];
sx q[1];
rz(-1.7971296) q[1];
sx q[1];
rz(2.8167579) q[1];
rz(2.0484106) q[2];
sx q[2];
rz(-1.7872767) q[2];
sx q[2];
rz(1.7959303) q[2];
rz(-1.6024797) q[3];
sx q[3];
rz(-1.0864232) q[3];
sx q[3];
rz(-0.3530799) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
