OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.49178034) q[0];
sx q[0];
rz(3.427504) q[0];
sx q[0];
rz(8.9094845) q[0];
rz(-1.7973068) q[1];
sx q[1];
rz(-0.15434115) q[1];
sx q[1];
rz(-0.57758346) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89864697) q[0];
sx q[0];
rz(-2.4624914) q[0];
sx q[0];
rz(-0.26429096) q[0];
x q[1];
rz(-2.0787813) q[2];
sx q[2];
rz(-2.0487818) q[2];
sx q[2];
rz(0.21131549) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5756302) q[1];
sx q[1];
rz(-1.168025) q[1];
sx q[1];
rz(0.29213841) q[1];
x q[2];
rz(-2.3732784) q[3];
sx q[3];
rz(-0.81211219) q[3];
sx q[3];
rz(-2.1384359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.704533) q[2];
sx q[2];
rz(-1.5960863) q[2];
sx q[2];
rz(0.68721592) q[2];
rz(2.1263188) q[3];
sx q[3];
rz(-1.3736558) q[3];
sx q[3];
rz(-0.12250531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9706443) q[0];
sx q[0];
rz(-2.0630554) q[0];
sx q[0];
rz(-1.2600391) q[0];
rz(2.1353703) q[1];
sx q[1];
rz(-2.1496014) q[1];
sx q[1];
rz(0.84567436) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0592247) q[0];
sx q[0];
rz(-2.0031643) q[0];
sx q[0];
rz(2.5655377) q[0];
x q[1];
rz(0.77180736) q[2];
sx q[2];
rz(-2.9559921) q[2];
sx q[2];
rz(2.0303357) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7686292) q[1];
sx q[1];
rz(-1.0781204) q[1];
sx q[1];
rz(1.782934) q[1];
x q[2];
rz(2.1917079) q[3];
sx q[3];
rz(-1.1511027) q[3];
sx q[3];
rz(-1.3148395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3530897) q[2];
sx q[2];
rz(-0.22558364) q[2];
sx q[2];
rz(0.4804002) q[2];
rz(-1.3530312) q[3];
sx q[3];
rz(-1.055911) q[3];
sx q[3];
rz(-1.1876748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1903494) q[0];
sx q[0];
rz(-0.23878637) q[0];
sx q[0];
rz(0.7730661) q[0];
rz(-3.0103325) q[1];
sx q[1];
rz(-1.8570329) q[1];
sx q[1];
rz(-2.0551596) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4341136) q[0];
sx q[0];
rz(-1.1093603) q[0];
sx q[0];
rz(-3.1397318) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2181182) q[2];
sx q[2];
rz(-1.4787276) q[2];
sx q[2];
rz(-2.0331241) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4194581) q[1];
sx q[1];
rz(-0.28885435) q[1];
sx q[1];
rz(-0.018317776) q[1];
x q[2];
rz(1.1658737) q[3];
sx q[3];
rz(-1.4110663) q[3];
sx q[3];
rz(2.3111642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1290258) q[2];
sx q[2];
rz(-1.4386703) q[2];
sx q[2];
rz(1.770299) q[2];
rz(-0.38315547) q[3];
sx q[3];
rz(-1.8846735) q[3];
sx q[3];
rz(-2.3390521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4521769) q[0];
sx q[0];
rz(-1.8912264) q[0];
sx q[0];
rz(0.048359811) q[0];
rz(2.9776749) q[1];
sx q[1];
rz(-0.36968958) q[1];
sx q[1];
rz(-1.4455459) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1235385) q[0];
sx q[0];
rz(-0.92045438) q[0];
sx q[0];
rz(1.7975848) q[0];
rz(0.21638685) q[2];
sx q[2];
rz(-1.4828223) q[2];
sx q[2];
rz(2.3073334) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2913937) q[1];
sx q[1];
rz(-1.7899917) q[1];
sx q[1];
rz(-1.6367957) q[1];
rz(1.0936071) q[3];
sx q[3];
rz(-0.573199) q[3];
sx q[3];
rz(-1.5447865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.066102862) q[2];
sx q[2];
rz(-1.7843856) q[2];
sx q[2];
rz(-1.0243105) q[2];
rz(1.6131489) q[3];
sx q[3];
rz(-1.6201092) q[3];
sx q[3];
rz(2.9083692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7016474) q[0];
sx q[0];
rz(-0.82413903) q[0];
sx q[0];
rz(-1.2874999) q[0];
rz(-2.8225186) q[1];
sx q[1];
rz(-1.5417475) q[1];
sx q[1];
rz(0.85420001) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6720649) q[0];
sx q[0];
rz(-0.85692353) q[0];
sx q[0];
rz(3.1307427) q[0];
x q[1];
rz(0.48798497) q[2];
sx q[2];
rz(-2.2099566) q[2];
sx q[2];
rz(-1.1346863) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2687208) q[1];
sx q[1];
rz(-2.1495021) q[1];
sx q[1];
rz(-0.89923664) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0075931) q[3];
sx q[3];
rz(-1.0922722) q[3];
sx q[3];
rz(-2.0625045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.95191082) q[2];
sx q[2];
rz(-0.59331912) q[2];
sx q[2];
rz(-2.5642776) q[2];
rz(-2.632085) q[3];
sx q[3];
rz(-0.43764344) q[3];
sx q[3];
rz(0.90665162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1381056) q[0];
sx q[0];
rz(-1.071799) q[0];
sx q[0];
rz(-3.0694718) q[0];
rz(-2.0347319) q[1];
sx q[1];
rz(-0.51263428) q[1];
sx q[1];
rz(3.0153826) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0059347) q[0];
sx q[0];
rz(-1.5392443) q[0];
sx q[0];
rz(-0.21261442) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3665479) q[2];
sx q[2];
rz(-1.2942874) q[2];
sx q[2];
rz(1.0724049) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.26423745) q[1];
sx q[1];
rz(-0.74062956) q[1];
sx q[1];
rz(-1.2901558) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0522271) q[3];
sx q[3];
rz(-0.37399451) q[3];
sx q[3];
rz(0.75043375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.90298992) q[2];
sx q[2];
rz(-2.949252) q[2];
sx q[2];
rz(2.3664756) q[2];
rz(-0.827968) q[3];
sx q[3];
rz(-2.8505846) q[3];
sx q[3];
rz(-2.0194139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5489952) q[0];
sx q[0];
rz(-2.7153375) q[0];
sx q[0];
rz(-3.0431842) q[0];
rz(-1.1920284) q[1];
sx q[1];
rz(-1.8076618) q[1];
sx q[1];
rz(2.5820406) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0347621) q[0];
sx q[0];
rz(-2.1666399) q[0];
sx q[0];
rz(1.725127) q[0];
rz(-pi) q[1];
rz(2.9497629) q[2];
sx q[2];
rz(-2.8755113) q[2];
sx q[2];
rz(-1.1508458) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0911078) q[1];
sx q[1];
rz(-0.15639601) q[1];
sx q[1];
rz(2.1662103) q[1];
rz(-pi) q[2];
rz(0.26569326) q[3];
sx q[3];
rz(-0.64722792) q[3];
sx q[3];
rz(2.2006187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6825535) q[2];
sx q[2];
rz(-1.295853) q[2];
sx q[2];
rz(2.7977978) q[2];
rz(0.5665468) q[3];
sx q[3];
rz(-0.44851258) q[3];
sx q[3];
rz(-0.47376537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.375181) q[0];
sx q[0];
rz(-1.3324998) q[0];
sx q[0];
rz(2.4108316) q[0];
rz(2.9991951) q[1];
sx q[1];
rz(-1.2700894) q[1];
sx q[1];
rz(-2.2699845) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9974737) q[0];
sx q[0];
rz(-1.561957) q[0];
sx q[0];
rz(-3.0661422) q[0];
rz(0.40558221) q[2];
sx q[2];
rz(-1.431258) q[2];
sx q[2];
rz(1.9987193) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3155047) q[1];
sx q[1];
rz(-1.9313889) q[1];
sx q[1];
rz(-2.4396067) q[1];
rz(3.001508) q[3];
sx q[3];
rz(-0.9730556) q[3];
sx q[3];
rz(0.41543451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.72620755) q[2];
sx q[2];
rz(-1.0691079) q[2];
sx q[2];
rz(2.0020206) q[2];
rz(1.4987882) q[3];
sx q[3];
rz(-0.39396861) q[3];
sx q[3];
rz(-2.22877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9386439) q[0];
sx q[0];
rz(-1.4511755) q[0];
sx q[0];
rz(1.9198445) q[0];
rz(2.9755759) q[1];
sx q[1];
rz(-1.8211726) q[1];
sx q[1];
rz(1.5244012) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.411392) q[0];
sx q[0];
rz(-0.087326614) q[0];
sx q[0];
rz(1.9090396) q[0];
rz(-1.379307) q[2];
sx q[2];
rz(-2.0404437) q[2];
sx q[2];
rz(0.42275235) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9164239) q[1];
sx q[1];
rz(-1.7290219) q[1];
sx q[1];
rz(-1.2374864) q[1];
rz(-1.8839621) q[3];
sx q[3];
rz(-2.1861665) q[3];
sx q[3];
rz(-0.19389158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.021585492) q[2];
sx q[2];
rz(-1.465613) q[2];
sx q[2];
rz(0.35153708) q[2];
rz(-2.0848138) q[3];
sx q[3];
rz(-2.6119699) q[3];
sx q[3];
rz(0.74469152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64365023) q[0];
sx q[0];
rz(-2.239776) q[0];
sx q[0];
rz(1.3051916) q[0];
rz(-2.7611043) q[1];
sx q[1];
rz(-1.0419798) q[1];
sx q[1];
rz(2.8881853) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2449269) q[0];
sx q[0];
rz(-1.9650808) q[0];
sx q[0];
rz(0.11163296) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.017574) q[2];
sx q[2];
rz(-0.3728711) q[2];
sx q[2];
rz(-2.1602221) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.39292654) q[1];
sx q[1];
rz(-1.8112507) q[1];
sx q[1];
rz(-2.7886831) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9990342) q[3];
sx q[3];
rz(-1.0034475) q[3];
sx q[3];
rz(2.464307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5499251) q[2];
sx q[2];
rz(-2.2558236) q[2];
sx q[2];
rz(2.6386476) q[2];
rz(0.89899603) q[3];
sx q[3];
rz(-1.8476202) q[3];
sx q[3];
rz(-1.1635273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1702561) q[0];
sx q[0];
rz(-1.6032871) q[0];
sx q[0];
rz(0.26300318) q[0];
rz(0.7111711) q[1];
sx q[1];
rz(-2.053459) q[1];
sx q[1];
rz(-1.4278535) q[1];
rz(-1.3118369) q[2];
sx q[2];
rz(-1.2987518) q[2];
sx q[2];
rz(-2.1546298) q[2];
rz(-2.3861804) q[3];
sx q[3];
rz(-1.0676386) q[3];
sx q[3];
rz(-0.044015351) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];