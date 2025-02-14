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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0812936) q[0];
sx q[0];
rz(-0.64612389) q[0];
sx q[0];
rz(-2.6160014) q[0];
rz(-pi) q[1];
rz(-2.5147314) q[2];
sx q[2];
rz(-0.38598362) q[2];
sx q[2];
rz(-0.71039334) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.6055687) q[1];
sx q[1];
rz(-1.4711416) q[1];
sx q[1];
rz(-1.2936564) q[1];
rz(0.72812702) q[3];
sx q[3];
rz(-1.9890729) q[3];
sx q[3];
rz(2.3969216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.677864) q[2];
sx q[2];
rz(-2.4552796) q[2];
sx q[2];
rz(-0.65832552) q[2];
rz(-2.5074734) q[3];
sx q[3];
rz(-1.6987957) q[3];
sx q[3];
rz(-1.3707976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.051801) q[0];
sx q[0];
rz(-0.3312411) q[0];
sx q[0];
rz(-3.064503) q[0];
rz(-0.58473051) q[1];
sx q[1];
rz(-1.8010151) q[1];
sx q[1];
rz(1.9416521) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.752992) q[0];
sx q[0];
rz(-1.3738828) q[0];
sx q[0];
rz(2.0373175) q[0];
rz(0.029424981) q[2];
sx q[2];
rz(-1.7351073) q[2];
sx q[2];
rz(-1.8738333) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6226095) q[1];
sx q[1];
rz(-2.4138192) q[1];
sx q[1];
rz(0.98811291) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.55910965) q[3];
sx q[3];
rz(-1.6225132) q[3];
sx q[3];
rz(-2.8962108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8889019) q[2];
sx q[2];
rz(-2.1609047) q[2];
sx q[2];
rz(0.67548951) q[2];
rz(-0.0096970079) q[3];
sx q[3];
rz(-2.8986425) q[3];
sx q[3];
rz(0.01827904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6460687) q[0];
sx q[0];
rz(-1.4718453) q[0];
sx q[0];
rz(-0.75463265) q[0];
rz(-3.1309639) q[1];
sx q[1];
rz(-1.7517215) q[1];
sx q[1];
rz(-2.1307814) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6957798) q[0];
sx q[0];
rz(-2.5523529) q[0];
sx q[0];
rz(-2.1408098) q[0];
rz(1.5358244) q[2];
sx q[2];
rz(-0.681923) q[2];
sx q[2];
rz(0.65820314) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.39937035) q[1];
sx q[1];
rz(-0.74811799) q[1];
sx q[1];
rz(-3.0518603) q[1];
x q[2];
rz(-2.378024) q[3];
sx q[3];
rz(-2.0393622) q[3];
sx q[3];
rz(2.6811622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.90091554) q[2];
sx q[2];
rz(-2.7132576) q[2];
sx q[2];
rz(2.6348616) q[2];
rz(-1.2725376) q[3];
sx q[3];
rz(-1.7839909) q[3];
sx q[3];
rz(2.8633964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5469359) q[0];
sx q[0];
rz(-1.9173859) q[0];
sx q[0];
rz(1.1154255) q[0];
rz(1.8755272) q[1];
sx q[1];
rz(-1.5770117) q[1];
sx q[1];
rz(-2.26684) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6653671) q[0];
sx q[0];
rz(-1.3194808) q[0];
sx q[0];
rz(-0.091500207) q[0];
rz(-1.5953996) q[2];
sx q[2];
rz(-2.4877791) q[2];
sx q[2];
rz(0.023651274) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.61234867) q[1];
sx q[1];
rz(-1.9772286) q[1];
sx q[1];
rz(1.8608577) q[1];
rz(-pi) q[2];
rz(-0.17904277) q[3];
sx q[3];
rz(-0.87233018) q[3];
sx q[3];
rz(1.2579789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0227585) q[2];
sx q[2];
rz(-2.6648882) q[2];
sx q[2];
rz(1.215722) q[2];
rz(1.3974961) q[3];
sx q[3];
rz(-1.6179061) q[3];
sx q[3];
rz(-1.4309179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27595156) q[0];
sx q[0];
rz(-3.0396099) q[0];
sx q[0];
rz(1.6708466) q[0];
rz(1.3784846) q[1];
sx q[1];
rz(-2.3536317) q[1];
sx q[1];
rz(-1.7313622) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3138235) q[0];
sx q[0];
rz(-2.6010989) q[0];
sx q[0];
rz(-1.9825516) q[0];
rz(-pi) q[1];
x q[1];
rz(1.330014) q[2];
sx q[2];
rz(-2.2301939) q[2];
sx q[2];
rz(-1.3032152) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.2454589) q[1];
sx q[1];
rz(-1.8832379) q[1];
sx q[1];
rz(-1.5724843) q[1];
x q[2];
rz(1.1377119) q[3];
sx q[3];
rz(-1.2473462) q[3];
sx q[3];
rz(2.5057305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1269647) q[2];
sx q[2];
rz(-0.85863272) q[2];
sx q[2];
rz(-1.3237759) q[2];
rz(-1.7823559) q[3];
sx q[3];
rz(-1.8424572) q[3];
sx q[3];
rz(-0.38715473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1112261) q[0];
sx q[0];
rz(-1.4596326) q[0];
sx q[0];
rz(-3.1193745) q[0];
rz(2.4413595) q[1];
sx q[1];
rz(-1.8902238) q[1];
sx q[1];
rz(0.95692316) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45201225) q[0];
sx q[0];
rz(-2.4274024) q[0];
sx q[0];
rz(-0.22269188) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9238052) q[2];
sx q[2];
rz(-1.7484669) q[2];
sx q[2];
rz(-1.2051932) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2607402) q[1];
sx q[1];
rz(-1.2965974) q[1];
sx q[1];
rz(1.9874279) q[1];
x q[2];
rz(-2.2150008) q[3];
sx q[3];
rz(-2.6098688) q[3];
sx q[3];
rz(-0.34429541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2172829) q[2];
sx q[2];
rz(-2.0192912) q[2];
sx q[2];
rz(1.4420606) q[2];
rz(-2.0349272) q[3];
sx q[3];
rz(-0.60551867) q[3];
sx q[3];
rz(1.1520011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7149413) q[0];
sx q[0];
rz(-1.8639257) q[0];
sx q[0];
rz(-2.6746124) q[0];
rz(1.3106208) q[1];
sx q[1];
rz(-2.1881723) q[1];
sx q[1];
rz(2.7511645) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4680088) q[0];
sx q[0];
rz(-1.7112186) q[0];
sx q[0];
rz(1.6037081) q[0];
x q[1];
rz(-1.8653581) q[2];
sx q[2];
rz(-1.6737564) q[2];
sx q[2];
rz(2.0077133) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.55246395) q[1];
sx q[1];
rz(-2.2099582) q[1];
sx q[1];
rz(-2.5659849) q[1];
rz(-pi) q[2];
x q[2];
rz(0.14777811) q[3];
sx q[3];
rz(-2.1782603) q[3];
sx q[3];
rz(-2.1097418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9499669) q[2];
sx q[2];
rz(-1.9248631) q[2];
sx q[2];
rz(2.1245655) q[2];
rz(2.2022066) q[3];
sx q[3];
rz(-0.87299577) q[3];
sx q[3];
rz(2.2123607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4262714) q[0];
sx q[0];
rz(-0.66522288) q[0];
sx q[0];
rz(-1.9287047) q[0];
rz(0.60316482) q[1];
sx q[1];
rz(-1.1319356) q[1];
sx q[1];
rz(2.718198) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2379266) q[0];
sx q[0];
rz(-0.70495021) q[0];
sx q[0];
rz(-1.4506819) q[0];
rz(1.7216484) q[2];
sx q[2];
rz(-1.4948339) q[2];
sx q[2];
rz(1.1830038) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5535721) q[1];
sx q[1];
rz(-1.7436073) q[1];
sx q[1];
rz(2.9281479) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5633864) q[3];
sx q[3];
rz(-1.4009589) q[3];
sx q[3];
rz(1.8078723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6508871) q[2];
sx q[2];
rz(-0.31279534) q[2];
sx q[2];
rz(0.97839626) q[2];
rz(-0.69495106) q[3];
sx q[3];
rz(-1.9415104) q[3];
sx q[3];
rz(1.8470701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0438743) q[0];
sx q[0];
rz(-2.9053423) q[0];
sx q[0];
rz(3.0978715) q[0];
rz(1.9678736) q[1];
sx q[1];
rz(-0.81467384) q[1];
sx q[1];
rz(2.3497605) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4872015) q[0];
sx q[0];
rz(-0.7672317) q[0];
sx q[0];
rz(-2.5228398) q[0];
rz(-pi) q[1];
rz(-2.9390807) q[2];
sx q[2];
rz(-1.8648476) q[2];
sx q[2];
rz(2.685315) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9992912) q[1];
sx q[1];
rz(-0.55146927) q[1];
sx q[1];
rz(-0.89889185) q[1];
rz(-2.4867587) q[3];
sx q[3];
rz(-2.5635701) q[3];
sx q[3];
rz(0.27930799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.23058471) q[2];
sx q[2];
rz(-2.8318996) q[2];
sx q[2];
rz(1.9045551) q[2];
rz(0.48219314) q[3];
sx q[3];
rz(-0.92415205) q[3];
sx q[3];
rz(2.6528416) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.024254) q[0];
sx q[0];
rz(-2.2145705) q[0];
sx q[0];
rz(-1.3379958) q[0];
rz(0.85211873) q[1];
sx q[1];
rz(-1.2282649) q[1];
sx q[1];
rz(1.9911912) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8088732) q[0];
sx q[0];
rz(-2.3925663) q[0];
sx q[0];
rz(-0.82360424) q[0];
rz(-pi) q[1];
rz(-3.030376) q[2];
sx q[2];
rz(-1.885998) q[2];
sx q[2];
rz(-2.1325796) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.627189) q[1];
sx q[1];
rz(-0.41819388) q[1];
sx q[1];
rz(1.2209284) q[1];
x q[2];
rz(2.1701067) q[3];
sx q[3];
rz(-1.2462292) q[3];
sx q[3];
rz(0.49719405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.368025) q[2];
sx q[2];
rz(-1.0057534) q[2];
sx q[2];
rz(-0.31039882) q[2];
rz(2.5144905) q[3];
sx q[3];
rz(-2.1576364) q[3];
sx q[3];
rz(1.9703126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
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
rz(2.9318342) q[2];
sx q[2];
rz(-1.0804313) q[2];
sx q[2];
rz(-2.6624138) q[2];
rz(1.6321833) q[3];
sx q[3];
rz(-1.8922378) q[3];
sx q[3];
rz(0.031382244) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
