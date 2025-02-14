OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4189897) q[0];
sx q[0];
rz(-0.64282066) q[0];
sx q[0];
rz(1.2195725) q[0];
rz(0.83528432) q[1];
sx q[1];
rz(-1.3568027) q[1];
sx q[1];
rz(0.9019444) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.565958) q[0];
sx q[0];
rz(-2.1778641) q[0];
sx q[0];
rz(-1.3835039) q[0];
rz(-pi) q[1];
x q[1];
rz(0.38036687) q[2];
sx q[2];
rz(-1.3970301) q[2];
sx q[2];
rz(1.1443537) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7985417) q[1];
sx q[1];
rz(-1.3968811) q[1];
sx q[1];
rz(-1.2760722) q[1];
rz(-pi) q[2];
rz(2.2252623) q[3];
sx q[3];
rz(-2.188465) q[3];
sx q[3];
rz(-0.82181069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.96282643) q[2];
sx q[2];
rz(-1.0054532) q[2];
sx q[2];
rz(-0.82628769) q[2];
rz(1.4055584) q[3];
sx q[3];
rz(-2.8345351) q[3];
sx q[3];
rz(-0.74339408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5590782) q[0];
sx q[0];
rz(-0.40575108) q[0];
sx q[0];
rz(1.7412809) q[0];
rz(-1.1236313) q[1];
sx q[1];
rz(-0.39389899) q[1];
sx q[1];
rz(1.0981015) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.708824) q[0];
sx q[0];
rz(-1.4147804) q[0];
sx q[0];
rz(1.4601329) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33771659) q[2];
sx q[2];
rz(-1.8840839) q[2];
sx q[2];
rz(-0.17853436) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.84867618) q[1];
sx q[1];
rz(-1.6458938) q[1];
sx q[1];
rz(-3.0203222) q[1];
rz(2.9750175) q[3];
sx q[3];
rz(-1.9431356) q[3];
sx q[3];
rz(0.53322116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0565679) q[2];
sx q[2];
rz(-0.3173863) q[2];
sx q[2];
rz(1.8918096) q[2];
rz(1.7791087) q[3];
sx q[3];
rz(-2.1407849) q[3];
sx q[3];
rz(-1.484163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63610858) q[0];
sx q[0];
rz(-0.69378575) q[0];
sx q[0];
rz(-2.8662477) q[0];
rz(-2.6823726) q[1];
sx q[1];
rz(-2.1870435) q[1];
sx q[1];
rz(3.088248) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33063525) q[0];
sx q[0];
rz(-2.914408) q[0];
sx q[0];
rz(2.7869528) q[0];
rz(-pi) q[1];
rz(-1.1098451) q[2];
sx q[2];
rz(-1.58733) q[2];
sx q[2];
rz(1.2993882) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5596603) q[1];
sx q[1];
rz(-1.5301063) q[1];
sx q[1];
rz(0.55198021) q[1];
rz(-1.2417913) q[3];
sx q[3];
rz(-1.3734372) q[3];
sx q[3];
rz(-2.1579735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.82103819) q[2];
sx q[2];
rz(-0.36131636) q[2];
sx q[2];
rz(-1.6374755) q[2];
rz(-1.5004246) q[3];
sx q[3];
rz(-2.0913561) q[3];
sx q[3];
rz(-2.8659081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(2.7212873) q[0];
sx q[0];
rz(-0.5834226) q[0];
sx q[0];
rz(-2.6570008) q[0];
rz(1.2424319) q[1];
sx q[1];
rz(-2.395605) q[1];
sx q[1];
rz(0.34908435) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1734245) q[0];
sx q[0];
rz(-1.5281046) q[0];
sx q[0];
rz(2.1195565) q[0];
x q[1];
rz(0.075320638) q[2];
sx q[2];
rz(-1.9131129) q[2];
sx q[2];
rz(1.7750193) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.610747) q[1];
sx q[1];
rz(-1.3813003) q[1];
sx q[1];
rz(-2.9086472) q[1];
rz(-pi) q[2];
rz(2.7233299) q[3];
sx q[3];
rz(-1.0111448) q[3];
sx q[3];
rz(2.5734316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.35859534) q[2];
sx q[2];
rz(-1.9390257) q[2];
sx q[2];
rz(-0.45004582) q[2];
rz(2.3027244) q[3];
sx q[3];
rz(-1.5939555) q[3];
sx q[3];
rz(-2.94256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49804509) q[0];
sx q[0];
rz(-1.6723375) q[0];
sx q[0];
rz(-0.099844649) q[0];
rz(-1.8210583) q[1];
sx q[1];
rz(-2.1456199) q[1];
sx q[1];
rz(-2.5340396) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6927101) q[0];
sx q[0];
rz(-1.4878055) q[0];
sx q[0];
rz(1.6458428) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.54536087) q[2];
sx q[2];
rz(-1.4358476) q[2];
sx q[2];
rz(3.1317186) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.93170184) q[1];
sx q[1];
rz(-1.7875515) q[1];
sx q[1];
rz(-2.5561129) q[1];
rz(-pi) q[2];
rz(-1.3123061) q[3];
sx q[3];
rz(-1.6655169) q[3];
sx q[3];
rz(-0.80918377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1943835) q[2];
sx q[2];
rz(-2.3781229) q[2];
sx q[2];
rz(-2.6109931) q[2];
rz(2.7001906) q[3];
sx q[3];
rz(-0.97359052) q[3];
sx q[3];
rz(-0.94943625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43668231) q[0];
sx q[0];
rz(-1.7969776) q[0];
sx q[0];
rz(0.5740903) q[0];
rz(-0.87999815) q[1];
sx q[1];
rz(-1.1209542) q[1];
sx q[1];
rz(-1.7990187) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3898252) q[0];
sx q[0];
rz(-0.51503599) q[0];
sx q[0];
rz(1.1549321) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3903177) q[2];
sx q[2];
rz(-2.5715368) q[2];
sx q[2];
rz(1.7715724) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9989661) q[1];
sx q[1];
rz(-1.9311096) q[1];
sx q[1];
rz(0.25288881) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9632666) q[3];
sx q[3];
rz(-2.470782) q[3];
sx q[3];
rz(1.3516155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4673956) q[2];
sx q[2];
rz(-2.8574222) q[2];
sx q[2];
rz(0.44710844) q[2];
rz(2.1111264) q[3];
sx q[3];
rz(-1.6298529) q[3];
sx q[3];
rz(-0.38465056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7717188) q[0];
sx q[0];
rz(-1.8444703) q[0];
sx q[0];
rz(2.7287927) q[0];
rz(1.075047) q[1];
sx q[1];
rz(-2.0037035) q[1];
sx q[1];
rz(2.3379751) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8014288) q[0];
sx q[0];
rz(-1.0699843) q[0];
sx q[0];
rz(2.1636971) q[0];
rz(-2.2140131) q[2];
sx q[2];
rz(-2.0394435) q[2];
sx q[2];
rz(1.2609294) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0688096) q[1];
sx q[1];
rz(-0.78662614) q[1];
sx q[1];
rz(-1.5907611) q[1];
x q[2];
rz(2.9622224) q[3];
sx q[3];
rz(-1.6458047) q[3];
sx q[3];
rz(1.4574432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2794118) q[2];
sx q[2];
rz(-1.2181686) q[2];
sx q[2];
rz(-1.5009521) q[2];
rz(2.5189279) q[3];
sx q[3];
rz(-2.3707135) q[3];
sx q[3];
rz(-1.4890495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72614661) q[0];
sx q[0];
rz(-1.5203238) q[0];
sx q[0];
rz(0.55794445) q[0];
rz(-2.5803512) q[1];
sx q[1];
rz(-1.3713501) q[1];
sx q[1];
rz(1.7254613) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6218044) q[0];
sx q[0];
rz(-1.1237696) q[0];
sx q[0];
rz(-0.31479737) q[0];
x q[1];
rz(-2.2165856) q[2];
sx q[2];
rz(-0.46621284) q[2];
sx q[2];
rz(0.031161873) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0812644) q[1];
sx q[1];
rz(-0.34354478) q[1];
sx q[1];
rz(1.2765803) q[1];
x q[2];
rz(0.21702311) q[3];
sx q[3];
rz(-2.3742834) q[3];
sx q[3];
rz(-1.8930774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.8884362) q[2];
sx q[2];
rz(-2.3641868) q[2];
sx q[2];
rz(-2.1412795) q[2];
rz(0.12217626) q[3];
sx q[3];
rz(-1.5001985) q[3];
sx q[3];
rz(2.9848671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4633453) q[0];
sx q[0];
rz(-1.1868917) q[0];
sx q[0];
rz(-0.004322411) q[0];
rz(0.16383544) q[1];
sx q[1];
rz(-2.7163353) q[1];
sx q[1];
rz(1.77553) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28986606) q[0];
sx q[0];
rz(-2.4938739) q[0];
sx q[0];
rz(1.6169548) q[0];
rz(-pi) q[1];
rz(-1.2605569) q[2];
sx q[2];
rz(-1.9266591) q[2];
sx q[2];
rz(0.17525338) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7228702) q[1];
sx q[1];
rz(-1.3197749) q[1];
sx q[1];
rz(-1.6081564) q[1];
rz(-pi) q[2];
x q[2];
rz(0.86413278) q[3];
sx q[3];
rz(-0.88418761) q[3];
sx q[3];
rz(-0.1089801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0020478) q[2];
sx q[2];
rz(-0.9674955) q[2];
sx q[2];
rz(-2.5223993) q[2];
rz(-0.55245095) q[3];
sx q[3];
rz(-1.8360454) q[3];
sx q[3];
rz(-0.14510554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7904952) q[0];
sx q[0];
rz(-2.9572697) q[0];
sx q[0];
rz(2.6628394) q[0];
rz(1.152773) q[1];
sx q[1];
rz(-2.8515127) q[1];
sx q[1];
rz(0.94720381) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21175948) q[0];
sx q[0];
rz(-0.18778983) q[0];
sx q[0];
rz(1.7540356) q[0];
rz(-pi) q[1];
rz(2.8962063) q[2];
sx q[2];
rz(-2.1446501) q[2];
sx q[2];
rz(-1.4978486) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6852764) q[1];
sx q[1];
rz(-2.3446977) q[1];
sx q[1];
rz(2.0048098) q[1];
rz(-pi) q[2];
rz(-2.5371505) q[3];
sx q[3];
rz(-0.42057188) q[3];
sx q[3];
rz(0.67422359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7572215) q[2];
sx q[2];
rz(-0.20558509) q[2];
sx q[2];
rz(0.26665404) q[2];
rz(2.0742778) q[3];
sx q[3];
rz(-1.2131194) q[3];
sx q[3];
rz(2.1783569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.08854475) q[0];
sx q[0];
rz(-1.6041258) q[0];
sx q[0];
rz(-1.5969101) q[0];
rz(-2.0883941) q[1];
sx q[1];
rz(-1.9269301) q[1];
sx q[1];
rz(-2.3763837) q[1];
rz(-0.011372671) q[2];
sx q[2];
rz(-2.8995902) q[2];
sx q[2];
rz(-2.216969) q[2];
rz(-0.57907818) q[3];
sx q[3];
rz(-1.3828207) q[3];
sx q[3];
rz(-2.4170392) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
