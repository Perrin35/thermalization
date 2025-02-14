OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3044337) q[0];
sx q[0];
rz(-1.8456012) q[0];
sx q[0];
rz(1.0977828) q[0];
rz(-0.9552362) q[1];
sx q[1];
rz(-2.3999441) q[1];
sx q[1];
rz(2.9902966) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1681874) q[0];
sx q[0];
rz(-0.12790132) q[0];
sx q[0];
rz(-1.3692877) q[0];
rz(-pi) q[1];
rz(-1.382231) q[2];
sx q[2];
rz(-0.55743581) q[2];
sx q[2];
rz(-0.093890015) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.94231168) q[1];
sx q[1];
rz(-2.1122871) q[1];
sx q[1];
rz(-2.460206) q[1];
rz(-pi) q[2];
rz(0.80161867) q[3];
sx q[3];
rz(-1.773917) q[3];
sx q[3];
rz(2.5020848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0155045) q[2];
sx q[2];
rz(-1.8512923) q[2];
sx q[2];
rz(0.31910953) q[2];
rz(2.0305521) q[3];
sx q[3];
rz(-0.55912656) q[3];
sx q[3];
rz(-1.8266953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1183209) q[0];
sx q[0];
rz(-0.51902223) q[0];
sx q[0];
rz(0.58854377) q[0];
rz(0.59666807) q[1];
sx q[1];
rz(-1.8110954) q[1];
sx q[1];
rz(-2.9002424) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7241192) q[0];
sx q[0];
rz(-0.80097526) q[0];
sx q[0];
rz(-2.1885314) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0701635) q[2];
sx q[2];
rz(-1.0712581) q[2];
sx q[2];
rz(-2.3880434) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6728404) q[1];
sx q[1];
rz(-2.6939055) q[1];
sx q[1];
rz(-1.8022369) q[1];
rz(1.7123132) q[3];
sx q[3];
rz(-2.5114473) q[3];
sx q[3];
rz(0.41434789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0236464) q[2];
sx q[2];
rz(-1.4586552) q[2];
sx q[2];
rz(-0.092197593) q[2];
rz(-2.2495031) q[3];
sx q[3];
rz(-2.2690319) q[3];
sx q[3];
rz(1.0649072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43111619) q[0];
sx q[0];
rz(-1.5065864) q[0];
sx q[0];
rz(-3.0431252) q[0];
rz(-0.59421986) q[1];
sx q[1];
rz(-2.0829945) q[1];
sx q[1];
rz(-0.99064151) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4020549) q[0];
sx q[0];
rz(-1.4841813) q[0];
sx q[0];
rz(1.9461856) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1192693) q[2];
sx q[2];
rz(-2.2532941) q[2];
sx q[2];
rz(-0.85993953) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.65730598) q[1];
sx q[1];
rz(-2.4542744) q[1];
sx q[1];
rz(-1.9621852) q[1];
rz(-pi) q[2];
rz(0.22728592) q[3];
sx q[3];
rz(-1.1133988) q[3];
sx q[3];
rz(1.5279087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6777665) q[2];
sx q[2];
rz(-2.3550484) q[2];
sx q[2];
rz(-2.3393935) q[2];
rz(-1.5471316) q[3];
sx q[3];
rz(-1.0651257) q[3];
sx q[3];
rz(-2.2772363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7441854) q[0];
sx q[0];
rz(-0.35905251) q[0];
sx q[0];
rz(-1.8205951) q[0];
rz(-1.6162704) q[1];
sx q[1];
rz(-2.1883712) q[1];
sx q[1];
rz(-2.9811409) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65112075) q[0];
sx q[0];
rz(-1.3297538) q[0];
sx q[0];
rz(0.81220497) q[0];
rz(0.58893369) q[2];
sx q[2];
rz(-2.9351165) q[2];
sx q[2];
rz(1.4331417) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8453809) q[1];
sx q[1];
rz(-1.4176765) q[1];
sx q[1];
rz(-2.2232242) q[1];
rz(-0.53628199) q[3];
sx q[3];
rz(-1.2068401) q[3];
sx q[3];
rz(2.3126147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.82115951) q[2];
sx q[2];
rz(-1.6511788) q[2];
sx q[2];
rz(2.5725345) q[2];
rz(1.2218366) q[3];
sx q[3];
rz(-2.8025083) q[3];
sx q[3];
rz(3.0088185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64567599) q[0];
sx q[0];
rz(-2.3585632) q[0];
sx q[0];
rz(2.3107279) q[0];
rz(-1.2105385) q[1];
sx q[1];
rz(-1.4930864) q[1];
sx q[1];
rz(0.99162203) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0474393) q[0];
sx q[0];
rz(-1.2625361) q[0];
sx q[0];
rz(0.83931132) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2328202) q[2];
sx q[2];
rz(-2.4261279) q[2];
sx q[2];
rz(1.5275265) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6779532) q[1];
sx q[1];
rz(-1.4153743) q[1];
sx q[1];
rz(-2.3033106) q[1];
rz(-pi) q[2];
rz(-0.95007054) q[3];
sx q[3];
rz(-0.80214989) q[3];
sx q[3];
rz(-0.42250326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3811938) q[2];
sx q[2];
rz(-0.92279592) q[2];
sx q[2];
rz(-2.7823616) q[2];
rz(0.71581101) q[3];
sx q[3];
rz(-0.78407136) q[3];
sx q[3];
rz(-1.951096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0773709) q[0];
sx q[0];
rz(-1.8697898) q[0];
sx q[0];
rz(-2.6203058) q[0];
rz(2.298666) q[1];
sx q[1];
rz(-1.9631674) q[1];
sx q[1];
rz(3.0311323) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2409637) q[0];
sx q[0];
rz(-0.74375737) q[0];
sx q[0];
rz(2.7636011) q[0];
rz(-pi) q[1];
rz(-1.6104389) q[2];
sx q[2];
rz(-0.64293282) q[2];
sx q[2];
rz(-0.55243353) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0065436) q[1];
sx q[1];
rz(-2.035917) q[1];
sx q[1];
rz(0.22174447) q[1];
rz(-0.25909781) q[3];
sx q[3];
rz(-1.0631592) q[3];
sx q[3];
rz(2.6258084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9342186) q[2];
sx q[2];
rz(-1.5418345) q[2];
sx q[2];
rz(2.1477487) q[2];
rz(-2.5908616) q[3];
sx q[3];
rz(-0.5144853) q[3];
sx q[3];
rz(-1.5229092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9024502) q[0];
sx q[0];
rz(-0.37343326) q[0];
sx q[0];
rz(2.9504839) q[0];
rz(0.36901078) q[1];
sx q[1];
rz(-1.399682) q[1];
sx q[1];
rz(-0.63327995) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8275237) q[0];
sx q[0];
rz(-1.6681328) q[0];
sx q[0];
rz(-0.23989664) q[0];
rz(0.72004135) q[2];
sx q[2];
rz(-1.9332464) q[2];
sx q[2];
rz(1.8099648) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.25200462) q[1];
sx q[1];
rz(-1.377863) q[1];
sx q[1];
rz(1.7798406) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.19509372) q[3];
sx q[3];
rz(-0.70835241) q[3];
sx q[3];
rz(-0.1699902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9621027) q[2];
sx q[2];
rz(-1.3733764) q[2];
sx q[2];
rz(0.38086677) q[2];
rz(1.7535836) q[3];
sx q[3];
rz(-1.0264779) q[3];
sx q[3];
rz(-1.7109722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66985828) q[0];
sx q[0];
rz(-0.40100455) q[0];
sx q[0];
rz(1.5420472) q[0];
rz(-1.3767287) q[1];
sx q[1];
rz(-1.7574666) q[1];
sx q[1];
rz(-2.210604) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41808266) q[0];
sx q[0];
rz(-2.2049954) q[0];
sx q[0];
rz(-1.4618327) q[0];
rz(0.45371666) q[2];
sx q[2];
rz(-1.3435817) q[2];
sx q[2];
rz(-1.5308282) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.763995) q[1];
sx q[1];
rz(-0.74255172) q[1];
sx q[1];
rz(1.5686036) q[1];
rz(-pi) q[2];
rz(0.91522907) q[3];
sx q[3];
rz(-1.9977565) q[3];
sx q[3];
rz(-0.87414908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5652183) q[2];
sx q[2];
rz(-0.56052506) q[2];
sx q[2];
rz(0.068923846) q[2];
rz(-1.8779514) q[3];
sx q[3];
rz(-0.37216035) q[3];
sx q[3];
rz(0.69839683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-2.3987592) q[0];
sx q[0];
rz(-0.95785207) q[0];
sx q[0];
rz(3.0897019) q[0];
rz(2.9715111) q[1];
sx q[1];
rz(-0.64337987) q[1];
sx q[1];
rz(-0.35194078) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1387685) q[0];
sx q[0];
rz(-1.5252946) q[0];
sx q[0];
rz(-0.36624927) q[0];
rz(-0.96723084) q[2];
sx q[2];
rz(-1.033342) q[2];
sx q[2];
rz(3.1307901) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5173147) q[1];
sx q[1];
rz(-2.0950965) q[1];
sx q[1];
rz(-0.45714) q[1];
rz(-pi) q[2];
x q[2];
rz(0.073411302) q[3];
sx q[3];
rz(-1.2651099) q[3];
sx q[3];
rz(-2.2524407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6841782) q[2];
sx q[2];
rz(-0.64932051) q[2];
sx q[2];
rz(-0.43668288) q[2];
rz(-2.7501578) q[3];
sx q[3];
rz(-1.6792363) q[3];
sx q[3];
rz(0.88633886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6201685) q[0];
sx q[0];
rz(-2.4346209) q[0];
sx q[0];
rz(3.022505) q[0];
rz(-1.2991615) q[1];
sx q[1];
rz(-1.7487339) q[1];
sx q[1];
rz(1.7838759) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.085310809) q[0];
sx q[0];
rz(-1.1650024) q[0];
sx q[0];
rz(0.88078518) q[0];
rz(-pi) q[1];
rz(0.36116551) q[2];
sx q[2];
rz(-2.2884011) q[2];
sx q[2];
rz(-1.4360365) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4132781) q[1];
sx q[1];
rz(-2.4367146) q[1];
sx q[1];
rz(-2.610763) q[1];
rz(2.6068654) q[3];
sx q[3];
rz(-0.89294723) q[3];
sx q[3];
rz(-2.0224366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.208821) q[2];
sx q[2];
rz(-1.5522771) q[2];
sx q[2];
rz(-2.5271752) q[2];
rz(1.5116073) q[3];
sx q[3];
rz(-2.4698518) q[3];
sx q[3];
rz(3.0577799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83196249) q[0];
sx q[0];
rz(-2.2086668) q[0];
sx q[0];
rz(-2.8902239) q[0];
rz(2.108719) q[1];
sx q[1];
rz(-1.947247) q[1];
sx q[1];
rz(-1.4364545) q[1];
rz(2.15945) q[2];
sx q[2];
rz(-0.12599421) q[2];
sx q[2];
rz(-0.85859184) q[2];
rz(2.3336505) q[3];
sx q[3];
rz(-1.0026889) q[3];
sx q[3];
rz(-1.9849594) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
