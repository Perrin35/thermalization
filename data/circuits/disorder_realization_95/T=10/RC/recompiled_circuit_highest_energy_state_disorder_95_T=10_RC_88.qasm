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
rz(0.79701841) q[0];
sx q[0];
rz(0.92172829) q[0];
sx q[0];
rz(14.447639) q[0];
rz(-2.4951275) q[1];
sx q[1];
rz(-0.87352455) q[1];
sx q[1];
rz(-2.8096107) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2811738) q[0];
sx q[0];
rz(-1.5462942) q[0];
sx q[0];
rz(-1.1479196) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.69692287) q[2];
sx q[2];
rz(-1.8644608) q[2];
sx q[2];
rz(0.98048335) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0003854) q[1];
sx q[1];
rz(-1.6024593) q[1];
sx q[1];
rz(2.382676) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0044972) q[3];
sx q[3];
rz(-0.44586102) q[3];
sx q[3];
rz(-0.23030989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.476764) q[2];
sx q[2];
rz(-2.2085184) q[2];
sx q[2];
rz(-0.19922166) q[2];
rz(-2.9557989) q[3];
sx q[3];
rz(-1.4150861) q[3];
sx q[3];
rz(-1.7004405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.324447) q[0];
sx q[0];
rz(-0.48878336) q[0];
sx q[0];
rz(-2.4031438) q[0];
rz(-1.6959408) q[1];
sx q[1];
rz(-1.7414469) q[1];
sx q[1];
rz(-2.6587528) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2525185) q[0];
sx q[0];
rz(-1.7216428) q[0];
sx q[0];
rz(0.14223994) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1170399) q[2];
sx q[2];
rz(-0.38402176) q[2];
sx q[2];
rz(-2.6817317) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.079025913) q[1];
sx q[1];
rz(-1.1084659) q[1];
sx q[1];
rz(2.5231925) q[1];
rz(-pi) q[2];
rz(-2.9779997) q[3];
sx q[3];
rz(-1.0181659) q[3];
sx q[3];
rz(-1.9528509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8740251) q[2];
sx q[2];
rz(-1.678062) q[2];
sx q[2];
rz(1.2919424) q[2];
rz(-1.1235631) q[3];
sx q[3];
rz(-0.70190391) q[3];
sx q[3];
rz(-1.7263713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32347754) q[0];
sx q[0];
rz(-0.98387843) q[0];
sx q[0];
rz(1.4418607) q[0];
rz(1.9135176) q[1];
sx q[1];
rz(-1.2308729) q[1];
sx q[1];
rz(2.0679811) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96127733) q[0];
sx q[0];
rz(-0.89160365) q[0];
sx q[0];
rz(2.017867) q[0];
rz(-pi) q[1];
rz(-0.19811689) q[2];
sx q[2];
rz(-0.92495944) q[2];
sx q[2];
rz(0.4465296) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.57697188) q[1];
sx q[1];
rz(-0.70319093) q[1];
sx q[1];
rz(-1.6220051) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.4573044) q[3];
sx q[3];
rz(-1.4239582) q[3];
sx q[3];
rz(-2.4350805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.29232612) q[2];
sx q[2];
rz(-2.8385415) q[2];
sx q[2];
rz(1.0230052) q[2];
rz(-0.060221378) q[3];
sx q[3];
rz(-1.9950208) q[3];
sx q[3];
rz(-2.4165418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4802454) q[0];
sx q[0];
rz(-2.9954973) q[0];
sx q[0];
rz(0.51937854) q[0];
rz(2.5410779) q[1];
sx q[1];
rz(-2.0527716) q[1];
sx q[1];
rz(-2.3868938) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31904991) q[0];
sx q[0];
rz(-2.5124335) q[0];
sx q[0];
rz(2.1436917) q[0];
x q[1];
rz(-2.2738931) q[2];
sx q[2];
rz(-1.461173) q[2];
sx q[2];
rz(-0.90849691) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.939173) q[1];
sx q[1];
rz(-1.8853097) q[1];
sx q[1];
rz(-2.0862723) q[1];
rz(-1.174332) q[3];
sx q[3];
rz(-1.045634) q[3];
sx q[3];
rz(-1.7523867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.34919136) q[2];
sx q[2];
rz(-2.3442522) q[2];
sx q[2];
rz(2.965773) q[2];
rz(-2.0154121) q[3];
sx q[3];
rz(-1.7837985) q[3];
sx q[3];
rz(0.064182909) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.675451) q[0];
sx q[0];
rz(-1.6459246) q[0];
sx q[0];
rz(0.58188907) q[0];
rz(-1.9582845) q[1];
sx q[1];
rz(-0.71154037) q[1];
sx q[1];
rz(-1.8875095) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81800705) q[0];
sx q[0];
rz(-0.76420751) q[0];
sx q[0];
rz(-2.0463405) q[0];
rz(-2.3603201) q[2];
sx q[2];
rz(-2.5102977) q[2];
sx q[2];
rz(2.6897813) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.075335927) q[1];
sx q[1];
rz(-1.7793223) q[1];
sx q[1];
rz(-1.2700524) q[1];
x q[2];
rz(-1.759911) q[3];
sx q[3];
rz(-2.0678036) q[3];
sx q[3];
rz(-0.86732098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1648272) q[2];
sx q[2];
rz(-2.3259614) q[2];
sx q[2];
rz(2.7867479) q[2];
rz(1.9918848) q[3];
sx q[3];
rz(-1.0747654) q[3];
sx q[3];
rz(-0.84158516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0717764) q[0];
sx q[0];
rz(-2.8700097) q[0];
sx q[0];
rz(-0.040104453) q[0];
rz(1.4647723) q[1];
sx q[1];
rz(-0.96950871) q[1];
sx q[1];
rz(0.47223314) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4667763) q[0];
sx q[0];
rz(-0.32186723) q[0];
sx q[0];
rz(-2.1842514) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5750999) q[2];
sx q[2];
rz(-2.488778) q[2];
sx q[2];
rz(-0.091856591) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.80010092) q[1];
sx q[1];
rz(-1.8566378) q[1];
sx q[1];
rz(1.2263917) q[1];
rz(-pi) q[2];
rz(-0.7615044) q[3];
sx q[3];
rz(-1.780394) q[3];
sx q[3];
rz(2.6246894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8670696) q[2];
sx q[2];
rz(-0.7663061) q[2];
sx q[2];
rz(-1.6609894) q[2];
rz(0.04145043) q[3];
sx q[3];
rz(-0.71235123) q[3];
sx q[3];
rz(0.97431549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90254766) q[0];
sx q[0];
rz(-2.3426265) q[0];
sx q[0];
rz(2.3833185) q[0];
rz(2.7039418) q[1];
sx q[1];
rz(-1.2145019) q[1];
sx q[1];
rz(-0.95380107) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24702913) q[0];
sx q[0];
rz(-1.5252542) q[0];
sx q[0];
rz(-1.9532502) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.72118451) q[2];
sx q[2];
rz(-1.0892592) q[2];
sx q[2];
rz(2.6480716) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.071324997) q[1];
sx q[1];
rz(-0.70601459) q[1];
sx q[1];
rz(1.474375) q[1];
x q[2];
rz(-0.28080583) q[3];
sx q[3];
rz(-1.971741) q[3];
sx q[3];
rz(1.3459447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4147676) q[2];
sx q[2];
rz(-2.813952) q[2];
sx q[2];
rz(0.3978351) q[2];
rz(-1.8757403) q[3];
sx q[3];
rz(-1.7222907) q[3];
sx q[3];
rz(-2.6771767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9809113) q[0];
sx q[0];
rz(-0.92781624) q[0];
sx q[0];
rz(-2.8371147) q[0];
rz(-2.0092087) q[1];
sx q[1];
rz(-2.4696923) q[1];
sx q[1];
rz(1.0736046) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7614339) q[0];
sx q[0];
rz(-3.112987) q[0];
sx q[0];
rz(2.2053115) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4377549) q[2];
sx q[2];
rz(-2.6148301) q[2];
sx q[2];
rz(-2.0671699) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.5636812) q[1];
sx q[1];
rz(-1.8133687) q[1];
sx q[1];
rz(0.96176186) q[1];
rz(-3.0746634) q[3];
sx q[3];
rz(-2.4968377) q[3];
sx q[3];
rz(0.88551846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4094746) q[2];
sx q[2];
rz(-1.2827337) q[2];
sx q[2];
rz(0.049962433) q[2];
rz(0.91160715) q[3];
sx q[3];
rz(-2.215569) q[3];
sx q[3];
rz(0.91134206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6066345) q[0];
sx q[0];
rz(-1.3953403) q[0];
sx q[0];
rz(0.362679) q[0];
rz(-0.72049385) q[1];
sx q[1];
rz(-2.3299005) q[1];
sx q[1];
rz(1.9627176) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6828109) q[0];
sx q[0];
rz(-2.1401554) q[0];
sx q[0];
rz(0.95627934) q[0];
rz(-pi) q[1];
rz(2.8209053) q[2];
sx q[2];
rz(-1.4960519) q[2];
sx q[2];
rz(2.8976669) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8803823) q[1];
sx q[1];
rz(-0.19621135) q[1];
sx q[1];
rz(1.0537149) q[1];
x q[2];
rz(-0.54525996) q[3];
sx q[3];
rz(-1.2125748) q[3];
sx q[3];
rz(-2.6521366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.0044340213) q[2];
sx q[2];
rz(-0.68778554) q[2];
sx q[2];
rz(1.2031817) q[2];
rz(3.1249937) q[3];
sx q[3];
rz(-1.4998452) q[3];
sx q[3];
rz(1.871292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
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
rz(0.59852973) q[0];
sx q[0];
rz(-0.5235343) q[0];
sx q[0];
rz(-1.4434927) q[0];
rz(0.47830018) q[1];
sx q[1];
rz(-2.4595478) q[1];
sx q[1];
rz(1.9313448) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0296625) q[0];
sx q[0];
rz(-1.0138982) q[0];
sx q[0];
rz(0.98390376) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4912676) q[2];
sx q[2];
rz(-1.4286094) q[2];
sx q[2];
rz(1.1917758) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7164439) q[1];
sx q[1];
rz(-1.8952491) q[1];
sx q[1];
rz(1.5157881) q[1];
rz(-pi) q[2];
rz(0.5860255) q[3];
sx q[3];
rz(-0.3908433) q[3];
sx q[3];
rz(1.217921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.18834194) q[2];
sx q[2];
rz(-2.5276999) q[2];
sx q[2];
rz(-2.4284412) q[2];
rz(-2.5085416) q[3];
sx q[3];
rz(-0.96708599) q[3];
sx q[3];
rz(2.2658394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0093832) q[0];
sx q[0];
rz(-1.8108981) q[0];
sx q[0];
rz(-0.43999408) q[0];
rz(-0.72883365) q[1];
sx q[1];
rz(-2.3757917) q[1];
sx q[1];
rz(1.052312) q[1];
rz(3.0265774) q[2];
sx q[2];
rz(-1.9585522) q[2];
sx q[2];
rz(1.1481651) q[2];
rz(0.36076693) q[3];
sx q[3];
rz(-2.3008892) q[3];
sx q[3];
rz(0.6194722) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
