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
rz(1.4690118) q[0];
sx q[0];
rz(5.3634085) q[0];
sx q[0];
rz(10.07977) q[0];
rz(1.6545777) q[1];
sx q[1];
rz(-1.7188641) q[1];
sx q[1];
rz(-1.3230327) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1974063) q[0];
sx q[0];
rz(-1.2639158) q[0];
sx q[0];
rz(-2.5635864) q[0];
x q[1];
rz(2.8236515) q[2];
sx q[2];
rz(-1.3481209) q[2];
sx q[2];
rz(-0.26938619) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.536024) q[1];
sx q[1];
rz(-1.4711416) q[1];
sx q[1];
rz(-1.2936564) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4134656) q[3];
sx q[3];
rz(-1.9890729) q[3];
sx q[3];
rz(2.3969216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.46372867) q[2];
sx q[2];
rz(-0.68631309) q[2];
sx q[2];
rz(0.65832552) q[2];
rz(-2.5074734) q[3];
sx q[3];
rz(-1.6987957) q[3];
sx q[3];
rz(1.770795) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089791678) q[0];
sx q[0];
rz(-2.8103516) q[0];
sx q[0];
rz(-0.077089699) q[0];
rz(0.58473051) q[1];
sx q[1];
rz(-1.8010151) q[1];
sx q[1];
rz(-1.9416521) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28040417) q[0];
sx q[0];
rz(-1.1139835) q[0];
sx q[0];
rz(2.9218319) q[0];
rz(-pi) q[1];
rz(1.3951833) q[2];
sx q[2];
rz(-2.9746911) q[2];
sx q[2];
rz(1.4457955) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.89948326) q[1];
sx q[1];
rz(-0.98190183) q[1];
sx q[1];
rz(-2.6857826) q[1];
rz(0.55910965) q[3];
sx q[3];
rz(-1.5190795) q[3];
sx q[3];
rz(0.24538183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.2526907) q[2];
sx q[2];
rz(-2.1609047) q[2];
sx q[2];
rz(-0.67548951) q[2];
rz(3.1318956) q[3];
sx q[3];
rz(-0.24295013) q[3];
sx q[3];
rz(3.1233136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6460687) q[0];
sx q[0];
rz(-1.4718453) q[0];
sx q[0];
rz(-0.75463265) q[0];
rz(0.010628788) q[1];
sx q[1];
rz(-1.3898712) q[1];
sx q[1];
rz(2.1307814) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6957798) q[0];
sx q[0];
rz(-0.58923972) q[0];
sx q[0];
rz(-2.1408098) q[0];
x q[1];
rz(-2.25242) q[2];
sx q[2];
rz(-1.592836) q[2];
sx q[2];
rz(0.88543788) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.39937035) q[1];
sx q[1];
rz(-2.3934747) q[1];
sx q[1];
rz(-3.0518603) q[1];
rz(-pi) q[2];
rz(0.63186462) q[3];
sx q[3];
rz(-0.87040983) q[3];
sx q[3];
rz(2.4720342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2406771) q[2];
sx q[2];
rz(-2.7132576) q[2];
sx q[2];
rz(-2.6348616) q[2];
rz(1.2725376) q[3];
sx q[3];
rz(-1.7839909) q[3];
sx q[3];
rz(0.27819628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5946567) q[0];
sx q[0];
rz(-1.9173859) q[0];
sx q[0];
rz(1.1154255) q[0];
rz(1.2660654) q[1];
sx q[1];
rz(-1.5770117) q[1];
sx q[1];
rz(-0.87475264) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47622555) q[0];
sx q[0];
rz(-1.3194808) q[0];
sx q[0];
rz(0.091500207) q[0];
rz(-pi) q[1];
rz(-2.2244637) q[2];
sx q[2];
rz(-1.5558331) q[2];
sx q[2];
rz(1.5666759) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.529244) q[1];
sx q[1];
rz(-1.9772286) q[1];
sx q[1];
rz(1.8608577) q[1];
x q[2];
rz(-2.9625499) q[3];
sx q[3];
rz(-0.87233018) q[3];
sx q[3];
rz(1.8836138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0227585) q[2];
sx q[2];
rz(-0.47670445) q[2];
sx q[2];
rz(1.215722) q[2];
rz(1.7440965) q[3];
sx q[3];
rz(-1.6179061) q[3];
sx q[3];
rz(-1.7106748) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27595156) q[0];
sx q[0];
rz(-0.10198274) q[0];
sx q[0];
rz(1.4707461) q[0];
rz(1.3784846) q[1];
sx q[1];
rz(-0.78796092) q[1];
sx q[1];
rz(-1.4102304) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2988457) q[0];
sx q[0];
rz(-1.0797459) q[0];
sx q[0];
rz(-0.23570717) q[0];
rz(0.67362154) q[2];
sx q[2];
rz(-1.7604019) q[2];
sx q[2];
rz(0.11826917) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3258563) q[1];
sx q[1];
rz(-1.5691901) q[1];
sx q[1];
rz(2.8291507) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.35378176) q[3];
sx q[3];
rz(-1.980034) q[3];
sx q[3];
rz(1.080846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0146279) q[2];
sx q[2];
rz(-0.85863272) q[2];
sx q[2];
rz(-1.3237759) q[2];
rz(-1.7823559) q[3];
sx q[3];
rz(-1.2991354) q[3];
sx q[3];
rz(-2.7544379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0303665) q[0];
sx q[0];
rz(-1.4596326) q[0];
sx q[0];
rz(-3.1193745) q[0];
rz(2.4413595) q[1];
sx q[1];
rz(-1.8902238) q[1];
sx q[1];
rz(0.95692316) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9807463) q[0];
sx q[0];
rz(-2.2637746) q[0];
sx q[0];
rz(1.7599544) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9238052) q[2];
sx q[2];
rz(-1.7484669) q[2];
sx q[2];
rz(-1.9363994) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.8808525) q[1];
sx q[1];
rz(-1.2965974) q[1];
sx q[1];
rz(-1.9874279) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.92659183) q[3];
sx q[3];
rz(-0.53172382) q[3];
sx q[3];
rz(-0.34429541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2172829) q[2];
sx q[2];
rz(-2.0192912) q[2];
sx q[2];
rz(1.699532) q[2];
rz(1.1066655) q[3];
sx q[3];
rz(-2.536074) q[3];
sx q[3];
rz(1.9895915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7149413) q[0];
sx q[0];
rz(-1.277667) q[0];
sx q[0];
rz(0.46698025) q[0];
rz(1.8309719) q[1];
sx q[1];
rz(-0.95342031) q[1];
sx q[1];
rz(2.7511645) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.098179558) q[0];
sx q[0];
rz(-1.5382086) q[0];
sx q[0];
rz(-3.0010953) q[0];
rz(-1.2762345) q[2];
sx q[2];
rz(-1.6737564) q[2];
sx q[2];
rz(-2.0077133) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5891287) q[1];
sx q[1];
rz(-2.2099582) q[1];
sx q[1];
rz(0.57560779) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.95818635) q[3];
sx q[3];
rz(-1.4496007) q[3];
sx q[3];
rz(2.6874128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.19162576) q[2];
sx q[2];
rz(-1.9248631) q[2];
sx q[2];
rz(-2.1245655) q[2];
rz(-0.93938604) q[3];
sx q[3];
rz(-2.2685969) q[3];
sx q[3];
rz(-2.2123607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4262714) q[0];
sx q[0];
rz(-0.66522288) q[0];
sx q[0];
rz(-1.2128879) q[0];
rz(-0.60316482) q[1];
sx q[1];
rz(-1.1319356) q[1];
sx q[1];
rz(0.42339465) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.903666) q[0];
sx q[0];
rz(-2.4366424) q[0];
sx q[0];
rz(-1.6909107) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7216484) q[2];
sx q[2];
rz(-1.4948339) q[2];
sx q[2];
rz(-1.1830038) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.054476995) q[1];
sx q[1];
rz(-1.7810139) q[1];
sx q[1];
rz(1.3940548) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0984108) q[3];
sx q[3];
rz(-0.16999741) q[3];
sx q[3];
rz(1.7640597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6508871) q[2];
sx q[2];
rz(-2.8287973) q[2];
sx q[2];
rz(0.97839626) q[2];
rz(-0.69495106) q[3];
sx q[3];
rz(-1.2000822) q[3];
sx q[3];
rz(1.2945226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.097718358) q[0];
sx q[0];
rz(-0.23625034) q[0];
sx q[0];
rz(0.043721113) q[0];
rz(-1.9678736) q[1];
sx q[1];
rz(-2.3269188) q[1];
sx q[1];
rz(2.3497605) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39003644) q[0];
sx q[0];
rz(-1.9851713) q[0];
sx q[0];
rz(0.66585559) q[0];
x q[1];
rz(-1.2709684) q[2];
sx q[2];
rz(-1.3770896) q[2];
sx q[2];
rz(1.1739588) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3085295) q[1];
sx q[1];
rz(-1.9030142) q[1];
sx q[1];
rz(2.019472) q[1];
rz(-pi) q[2];
x q[2];
rz(0.65483396) q[3];
sx q[3];
rz(-0.5780226) q[3];
sx q[3];
rz(-0.27930799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.23058471) q[2];
sx q[2];
rz(-2.8318996) q[2];
sx q[2];
rz(1.9045551) q[2];
rz(2.6593995) q[3];
sx q[3];
rz(-0.92415205) q[3];
sx q[3];
rz(-2.6528416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.024254) q[0];
sx q[0];
rz(-2.2145705) q[0];
sx q[0];
rz(1.8035969) q[0];
rz(-0.85211873) q[1];
sx q[1];
rz(-1.9133277) q[1];
sx q[1];
rz(-1.1504014) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83417621) q[0];
sx q[0];
rz(-2.0518853) q[0];
sx q[0];
rz(-0.9721945) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8878292) q[2];
sx q[2];
rz(-1.46508) q[2];
sx q[2];
rz(-0.52717613) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4070421) q[1];
sx q[1];
rz(-1.7104539) q[1];
sx q[1];
rz(-1.9662844) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1086333) q[3];
sx q[3];
rz(-0.67194429) q[3];
sx q[3];
rz(2.5047461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.368025) q[2];
sx q[2];
rz(-2.1358392) q[2];
sx q[2];
rz(0.31039882) q[2];
rz(2.5144905) q[3];
sx q[3];
rz(-0.98395625) q[3];
sx q[3];
rz(1.1712801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-1.8695759) q[0];
sx q[0];
rz(-1.5740812) q[0];
sx q[0];
rz(1.5480702) q[0];
rz(1.9238453) q[1];
sx q[1];
rz(-1.4532614) q[1];
sx q[1];
rz(-0.90167602) q[1];
rz(-0.2097585) q[2];
sx q[2];
rz(-1.0804313) q[2];
sx q[2];
rz(-2.6624138) q[2];
rz(-2.8195856) q[3];
sx q[3];
rz(-1.5125572) q[3];
sx q[3];
rz(-1.5199979) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
