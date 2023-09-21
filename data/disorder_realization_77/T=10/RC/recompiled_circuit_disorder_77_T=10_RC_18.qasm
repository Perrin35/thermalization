OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.66184008) q[0];
sx q[0];
rz(-0.84364426) q[0];
sx q[0];
rz(-2.9736829) q[0];
rz(-1.9703938) q[1];
sx q[1];
rz(-0.29532239) q[1];
sx q[1];
rz(3.0854316) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1355609) q[0];
sx q[0];
rz(-0.52838415) q[0];
sx q[0];
rz(1.1392659) q[0];
rz(-0.21284717) q[2];
sx q[2];
rz(-2.2058862) q[2];
sx q[2];
rz(-1.1320621) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.010477) q[1];
sx q[1];
rz(-1.8382204) q[1];
sx q[1];
rz(-2.1285776) q[1];
x q[2];
rz(2.9763016) q[3];
sx q[3];
rz(-0.2632907) q[3];
sx q[3];
rz(-2.1804682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7636259) q[2];
sx q[2];
rz(-0.28181919) q[2];
sx q[2];
rz(-2.7089233) q[2];
rz(-1.9487322) q[3];
sx q[3];
rz(-1.9038707) q[3];
sx q[3];
rz(2.7584934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52779657) q[0];
sx q[0];
rz(-0.48848099) q[0];
sx q[0];
rz(1.8288076) q[0];
rz(0.20547543) q[1];
sx q[1];
rz(-0.97646362) q[1];
sx q[1];
rz(-1.1516494) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5702471) q[0];
sx q[0];
rz(-2.3803108) q[0];
sx q[0];
rz(2.0061357) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4685681) q[2];
sx q[2];
rz(-0.7012127) q[2];
sx q[2];
rz(-0.53403026) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0484867) q[1];
sx q[1];
rz(-1.9890607) q[1];
sx q[1];
rz(-0.47936819) q[1];
rz(-pi) q[2];
rz(1.1851951) q[3];
sx q[3];
rz(-1.10023) q[3];
sx q[3];
rz(0.040599559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1318704) q[2];
sx q[2];
rz(-1.6513609) q[2];
sx q[2];
rz(0.22182626) q[2];
rz(2.7644073) q[3];
sx q[3];
rz(-0.42755869) q[3];
sx q[3];
rz(-0.77243531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.8310228) q[0];
sx q[0];
rz(-3.0492058) q[0];
sx q[0];
rz(-3.1047399) q[0];
rz(-2.316078) q[1];
sx q[1];
rz(-1.8258391) q[1];
sx q[1];
rz(-3.085014) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53783137) q[0];
sx q[0];
rz(-2.1268401) q[0];
sx q[0];
rz(-2.6105196) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7043731) q[2];
sx q[2];
rz(-1.3165511) q[2];
sx q[2];
rz(2.996252) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1287071) q[1];
sx q[1];
rz(-1.371908) q[1];
sx q[1];
rz(-0.30602869) q[1];
rz(-0.64485456) q[3];
sx q[3];
rz(-1.8861024) q[3];
sx q[3];
rz(-0.23526084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8686707) q[2];
sx q[2];
rz(-1.6438831) q[2];
sx q[2];
rz(2.2154714) q[2];
rz(-0.55666322) q[3];
sx q[3];
rz(-0.29354468) q[3];
sx q[3];
rz(1.042897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91519231) q[0];
sx q[0];
rz(-2.4214348) q[0];
sx q[0];
rz(0.74209374) q[0];
rz(1.1391976) q[1];
sx q[1];
rz(-2.6622055) q[1];
sx q[1];
rz(2.6779968) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6620561) q[0];
sx q[0];
rz(-2.2502796) q[0];
sx q[0];
rz(0.38328538) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6323339) q[2];
sx q[2];
rz(-1.9029641) q[2];
sx q[2];
rz(-0.57519826) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0193034) q[1];
sx q[1];
rz(-0.79320723) q[1];
sx q[1];
rz(-1.3293468) q[1];
rz(-pi) q[2];
x q[2];
rz(0.87426825) q[3];
sx q[3];
rz(-0.6797176) q[3];
sx q[3];
rz(0.10248871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7745557) q[2];
sx q[2];
rz(-2.98428) q[2];
sx q[2];
rz(-0.049499361) q[2];
rz(0.1285304) q[3];
sx q[3];
rz(-1.5934207) q[3];
sx q[3];
rz(-3.0226504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0054935) q[0];
sx q[0];
rz(-2.7042784) q[0];
sx q[0];
rz(2.8438925) q[0];
rz(-2.659335) q[1];
sx q[1];
rz(-0.75459701) q[1];
sx q[1];
rz(0.94435) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0194861) q[0];
sx q[0];
rz(-1.5257611) q[0];
sx q[0];
rz(1.3615863) q[0];
x q[1];
rz(-1.1733426) q[2];
sx q[2];
rz(-2.3059418) q[2];
sx q[2];
rz(0.022692516) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3790834) q[1];
sx q[1];
rz(-3.0804539) q[1];
sx q[1];
rz(1.6054543) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.827042) q[3];
sx q[3];
rz(-1.991193) q[3];
sx q[3];
rz(-2.8375569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.258761) q[2];
sx q[2];
rz(-2.0488887) q[2];
sx q[2];
rz(-0.10822254) q[2];
rz(-3.1392858) q[3];
sx q[3];
rz(-1.6121515) q[3];
sx q[3];
rz(0.32430696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7047983) q[0];
sx q[0];
rz(-0.44610766) q[0];
sx q[0];
rz(-2.5571402) q[0];
rz(-0.8862409) q[1];
sx q[1];
rz(-0.61683547) q[1];
sx q[1];
rz(0.054919682) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14558218) q[0];
sx q[0];
rz(-1.3025563) q[0];
sx q[0];
rz(-0.085573816) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1228742) q[2];
sx q[2];
rz(-1.2476377) q[2];
sx q[2];
rz(1.9402372) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0262895) q[1];
sx q[1];
rz(-0.68740986) q[1];
sx q[1];
rz(0.61647146) q[1];
rz(0.46338007) q[3];
sx q[3];
rz(-2.1043092) q[3];
sx q[3];
rz(-2.2802071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3871258) q[2];
sx q[2];
rz(-3.0209164) q[2];
sx q[2];
rz(-1.0167271) q[2];
rz(2.5975442) q[3];
sx q[3];
rz(-2.7816911) q[3];
sx q[3];
rz(1.8090766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.465437) q[0];
sx q[0];
rz(-0.98452079) q[0];
sx q[0];
rz(2.8570535) q[0];
rz(-2.1971205) q[1];
sx q[1];
rz(-1.1962793) q[1];
sx q[1];
rz(-0.91032666) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7401687) q[0];
sx q[0];
rz(-1.5059885) q[0];
sx q[0];
rz(-0.054697371) q[0];
rz(1.1548642) q[2];
sx q[2];
rz(-1.2565194) q[2];
sx q[2];
rz(1.2736125) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0525166) q[1];
sx q[1];
rz(-2.0564582) q[1];
sx q[1];
rz(-1.3658701) q[1];
x q[2];
rz(2.4207553) q[3];
sx q[3];
rz(-1.1158873) q[3];
sx q[3];
rz(2.7618559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3900782) q[2];
sx q[2];
rz(-0.063515924) q[2];
sx q[2];
rz(-0.92203036) q[2];
rz(2.5743124) q[3];
sx q[3];
rz(-1.7276948) q[3];
sx q[3];
rz(-1.0197619) q[3];
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
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89408016) q[0];
sx q[0];
rz(-2.4649354) q[0];
sx q[0];
rz(0.12938736) q[0];
rz(2.5091876) q[1];
sx q[1];
rz(-2.1148732) q[1];
sx q[1];
rz(2.8410889) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.520641) q[0];
sx q[0];
rz(-1.9410987) q[0];
sx q[0];
rz(-0.83604367) q[0];
rz(-pi) q[1];
rz(2.1883165) q[2];
sx q[2];
rz(-2.2308908) q[2];
sx q[2];
rz(-2.2613139) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0771675) q[1];
sx q[1];
rz(-0.38591138) q[1];
sx q[1];
rz(-2.6618631) q[1];
rz(-pi) q[2];
rz(1.9037876) q[3];
sx q[3];
rz(-0.78562842) q[3];
sx q[3];
rz(-1.9612519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.58632103) q[2];
sx q[2];
rz(-0.99493146) q[2];
sx q[2];
rz(2.3596181) q[2];
rz(0.55109763) q[3];
sx q[3];
rz(-1.3827773) q[3];
sx q[3];
rz(-2.9836695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5683811) q[0];
sx q[0];
rz(-1.1567572) q[0];
sx q[0];
rz(-3.0138299) q[0];
rz(-2.5993775) q[1];
sx q[1];
rz(-2.1844889) q[1];
sx q[1];
rz(-2.382747) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42620537) q[0];
sx q[0];
rz(-1.1466768) q[0];
sx q[0];
rz(-3.12294) q[0];
rz(-2.8685832) q[2];
sx q[2];
rz(-0.62676478) q[2];
sx q[2];
rz(-0.11944709) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.54770494) q[1];
sx q[1];
rz(-0.87915671) q[1];
sx q[1];
rz(-1.3681075) q[1];
rz(-pi) q[2];
rz(-1.5522478) q[3];
sx q[3];
rz(-2.4759001) q[3];
sx q[3];
rz(-1.9513643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1252497) q[2];
sx q[2];
rz(-1.7813851) q[2];
sx q[2];
rz(2.6515567) q[2];
rz(-1.7193433) q[3];
sx q[3];
rz(-1.9336721) q[3];
sx q[3];
rz(1.0629883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7816417) q[0];
sx q[0];
rz(-2.521823) q[0];
sx q[0];
rz(-3.066257) q[0];
rz(-0.8967337) q[1];
sx q[1];
rz(-1.1743841) q[1];
sx q[1];
rz(0.60992253) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2595554) q[0];
sx q[0];
rz(-1.2872739) q[0];
sx q[0];
rz(-2.3638704) q[0];
x q[1];
rz(-2.5209849) q[2];
sx q[2];
rz(-0.45244103) q[2];
sx q[2];
rz(-1.1021745) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8049106) q[1];
sx q[1];
rz(-1.8098127) q[1];
sx q[1];
rz(-2.3906624) q[1];
x q[2];
rz(-0.53819733) q[3];
sx q[3];
rz(-0.81848577) q[3];
sx q[3];
rz(-1.324211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.23218368) q[2];
sx q[2];
rz(-0.82513088) q[2];
sx q[2];
rz(-0.71371901) q[2];
rz(2.7632726) q[3];
sx q[3];
rz(-0.49351966) q[3];
sx q[3];
rz(0.87987125) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.338035) q[0];
sx q[0];
rz(-1.150158) q[0];
sx q[0];
rz(-1.5858571) q[0];
rz(0.65080416) q[1];
sx q[1];
rz(-1.6497859) q[1];
sx q[1];
rz(-0.12129687) q[1];
rz(-1.3748319) q[2];
sx q[2];
rz(-1.5407731) q[2];
sx q[2];
rz(-2.4827448) q[2];
rz(-2.0402504) q[3];
sx q[3];
rz(-0.51870844) q[3];
sx q[3];
rz(-1.4389256) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
