OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6498123) q[0];
sx q[0];
rz(-0.28591135) q[0];
sx q[0];
rz(-2.6262992) q[0];
rz(-1.7973068) q[1];
sx q[1];
rz(-0.15434115) q[1];
sx q[1];
rz(-0.57758346) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87969765) q[0];
sx q[0];
rz(-1.4059773) q[0];
sx q[0];
rz(2.4796955) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0628113) q[2];
sx q[2];
rz(-1.0928109) q[2];
sx q[2];
rz(0.21131549) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9211728) q[1];
sx q[1];
rz(-0.49282679) q[1];
sx q[1];
rz(-2.1652031) q[1];
x q[2];
rz(-0.9382117) q[3];
sx q[3];
rz(-2.1198366) q[3];
sx q[3];
rz(0.051018056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.704533) q[2];
sx q[2];
rz(-1.5960863) q[2];
sx q[2];
rz(-2.4543767) q[2];
rz(-1.0152738) q[3];
sx q[3];
rz(-1.3736558) q[3];
sx q[3];
rz(-0.12250531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17094831) q[0];
sx q[0];
rz(-2.0630554) q[0];
sx q[0];
rz(-1.8815536) q[0];
rz(-2.1353703) q[1];
sx q[1];
rz(-0.99199122) q[1];
sx q[1];
rz(0.84567436) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2226919) q[0];
sx q[0];
rz(-1.0534304) q[0];
sx q[0];
rz(2.0738686) q[0];
rz(1.7010062) q[2];
sx q[2];
rz(-1.4381593) q[2];
sx q[2];
rz(2.8108033) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0963124) q[1];
sx q[1];
rz(-1.3842061) q[1];
sx q[1];
rz(-2.6394096) q[1];
rz(2.2250697) q[3];
sx q[3];
rz(-2.407981) q[3];
sx q[3];
rz(-0.77378002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.78850293) q[2];
sx q[2];
rz(-0.22558364) q[2];
sx q[2];
rz(-2.6611924) q[2];
rz(1.3530312) q[3];
sx q[3];
rz(-1.055911) q[3];
sx q[3];
rz(-1.9539179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95124328) q[0];
sx q[0];
rz(-2.9028063) q[0];
sx q[0];
rz(0.7730661) q[0];
rz(-3.0103325) q[1];
sx q[1];
rz(-1.2845598) q[1];
sx q[1];
rz(-1.0864331) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13585424) q[0];
sx q[0];
rz(-1.5724626) q[0];
sx q[0];
rz(1.1093596) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.418872) q[2];
sx q[2];
rz(-0.65290367) q[2];
sx q[2];
rz(-0.34130794) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.8311027) q[1];
sx q[1];
rz(-1.576014) q[1];
sx q[1];
rz(2.8527841) q[1];
rz(-pi) q[2];
rz(1.1826035) q[3];
sx q[3];
rz(-0.43366323) q[3];
sx q[3];
rz(-2.7565623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1290258) q[2];
sx q[2];
rz(-1.4386703) q[2];
sx q[2];
rz(-1.770299) q[2];
rz(2.7584372) q[3];
sx q[3];
rz(-1.2569191) q[3];
sx q[3];
rz(-0.80254054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6894158) q[0];
sx q[0];
rz(-1.8912264) q[0];
sx q[0];
rz(0.048359811) q[0];
rz(2.9776749) q[1];
sx q[1];
rz(-2.7719031) q[1];
sx q[1];
rz(1.4455459) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7955129) q[0];
sx q[0];
rz(-0.68329408) q[0];
sx q[0];
rz(0.28738316) q[0];
rz(-pi) q[1];
rz(-2.7518026) q[2];
sx q[2];
rz(-2.9082657) q[2];
sx q[2];
rz(2.0248272) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9963035) q[1];
sx q[1];
rz(-2.9128296) q[1];
sx q[1];
rz(-2.8537675) q[1];
rz(-pi) q[2];
rz(-2.8533832) q[3];
sx q[3];
rz(-2.0734348) q[3];
sx q[3];
rz(0.99311815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0754898) q[2];
sx q[2];
rz(-1.3572071) q[2];
sx q[2];
rz(1.0243105) q[2];
rz(1.5284437) q[3];
sx q[3];
rz(-1.5214835) q[3];
sx q[3];
rz(2.9083692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7016474) q[0];
sx q[0];
rz(-0.82413903) q[0];
sx q[0];
rz(-1.8540927) q[0];
rz(-2.8225186) q[1];
sx q[1];
rz(-1.5998452) q[1];
sx q[1];
rz(-0.85420001) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.655495) q[0];
sx q[0];
rz(-2.4276519) q[0];
sx q[0];
rz(1.5582725) q[0];
rz(0.48798497) q[2];
sx q[2];
rz(-2.2099566) q[2];
sx q[2];
rz(2.0069063) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2367868) q[1];
sx q[1];
rz(-0.85610897) q[1];
sx q[1];
rz(-2.3805815) q[1];
rz(-pi) q[2];
rz(2.1339995) q[3];
sx q[3];
rz(-1.0922722) q[3];
sx q[3];
rz(2.0625045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1896818) q[2];
sx q[2];
rz(-0.59331912) q[2];
sx q[2];
rz(0.577315) q[2];
rz(-0.50950766) q[3];
sx q[3];
rz(-0.43764344) q[3];
sx q[3];
rz(-0.90665162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
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
rz(0.072120897) q[0];
rz(1.1068608) q[1];
sx q[1];
rz(-0.51263428) q[1];
sx q[1];
rz(3.0153826) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55805154) q[0];
sx q[0];
rz(-1.3582894) q[0];
sx q[0];
rz(1.6030747) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3665479) q[2];
sx q[2];
rz(-1.2942874) q[2];
sx q[2];
rz(2.0691878) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5050161) q[1];
sx q[1];
rz(-0.86537433) q[1];
sx q[1];
rz(0.24800639) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8993127) q[3];
sx q[3];
rz(-1.3887172) q[3];
sx q[3];
rz(1.3086705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.90298992) q[2];
sx q[2];
rz(-0.1923407) q[2];
sx q[2];
rz(0.77511707) q[2];
rz(-2.3136247) q[3];
sx q[3];
rz(-2.8505846) q[3];
sx q[3];
rz(-1.1221788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5489952) q[0];
sx q[0];
rz(-2.7153375) q[0];
sx q[0];
rz(3.0431842) q[0];
rz(1.1920284) q[1];
sx q[1];
rz(-1.8076618) q[1];
sx q[1];
rz(-2.5820406) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10683051) q[0];
sx q[0];
rz(-2.1666399) q[0];
sx q[0];
rz(-1.4164657) q[0];
x q[1];
rz(0.19182972) q[2];
sx q[2];
rz(-2.8755113) q[2];
sx q[2];
rz(1.1508458) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4899788) q[1];
sx q[1];
rz(-1.7001121) q[1];
sx q[1];
rz(3.053385) q[1];
rz(2.8758994) q[3];
sx q[3];
rz(-0.64722792) q[3];
sx q[3];
rz(0.94097394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.45903912) q[2];
sx q[2];
rz(-1.295853) q[2];
sx q[2];
rz(-0.34379488) q[2];
rz(-0.5665468) q[3];
sx q[3];
rz(-0.44851258) q[3];
sx q[3];
rz(-2.6678273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.375181) q[0];
sx q[0];
rz(-1.3324998) q[0];
sx q[0];
rz(-2.4108316) q[0];
rz(2.9991951) q[1];
sx q[1];
rz(-1.2700894) q[1];
sx q[1];
rz(0.87160814) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14411892) q[0];
sx q[0];
rz(-1.5796356) q[0];
sx q[0];
rz(-3.0661422) q[0];
rz(1.7224738) q[2];
sx q[2];
rz(-1.9722087) q[2];
sx q[2];
rz(0.36827189) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3155047) q[1];
sx q[1];
rz(-1.9313889) q[1];
sx q[1];
rz(-2.4396067) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7730764) q[3];
sx q[3];
rz(-2.5296122) q[3];
sx q[3];
rz(0.66093854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.72620755) q[2];
sx q[2];
rz(-2.0724847) q[2];
sx q[2];
rz(2.0020206) q[2];
rz(1.4987882) q[3];
sx q[3];
rz(-0.39396861) q[3];
sx q[3];
rz(0.9128226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9386439) q[0];
sx q[0];
rz(-1.6904172) q[0];
sx q[0];
rz(-1.2217481) q[0];
rz(-2.9755759) q[1];
sx q[1];
rz(-1.32042) q[1];
sx q[1];
rz(-1.6171914) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.411392) q[0];
sx q[0];
rz(-0.087326614) q[0];
sx q[0];
rz(1.9090396) q[0];
rz(-2.6644601) q[2];
sx q[2];
rz(-1.7413365) q[2];
sx q[2];
rz(1.0605304) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0601378) q[1];
sx q[1];
rz(-2.7739035) q[1];
sx q[1];
rz(-2.0245488) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.41096656) q[3];
sx q[3];
rz(-2.4604359) q[3];
sx q[3];
rz(-2.4364803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.021585492) q[2];
sx q[2];
rz(-1.465613) q[2];
sx q[2];
rz(2.7900556) q[2];
rz(-1.0567788) q[3];
sx q[3];
rz(-2.6119699) q[3];
sx q[3];
rz(-0.74469152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4979424) q[0];
sx q[0];
rz(-2.239776) q[0];
sx q[0];
rz(-1.836401) q[0];
rz(2.7611043) q[1];
sx q[1];
rz(-1.0419798) q[1];
sx q[1];
rz(-2.8881853) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2828335) q[0];
sx q[0];
rz(-1.4677605) q[0];
sx q[0];
rz(1.1742924) q[0];
x q[1];
rz(-1.1240187) q[2];
sx q[2];
rz(-2.7687216) q[2];
sx q[2];
rz(0.98137059) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0903783) q[1];
sx q[1];
rz(-1.2284632) q[1];
sx q[1];
rz(1.315209) q[1];
rz(-pi) q[2];
rz(1.1425584) q[3];
sx q[3];
rz(-1.0034475) q[3];
sx q[3];
rz(-2.464307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5499251) q[2];
sx q[2];
rz(-0.88576907) q[2];
sx q[2];
rz(-2.6386476) q[2];
rz(-0.89899603) q[3];
sx q[3];
rz(-1.2939724) q[3];
sx q[3];
rz(-1.1635273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9713365) q[0];
sx q[0];
rz(-1.5383056) q[0];
sx q[0];
rz(-2.8785895) q[0];
rz(-2.4304216) q[1];
sx q[1];
rz(-2.053459) q[1];
sx q[1];
rz(-1.4278535) q[1];
rz(1.3118369) q[2];
sx q[2];
rz(-1.8428409) q[2];
sx q[2];
rz(0.98696282) q[2];
rz(0.67646277) q[3];
sx q[3];
rz(-0.87920311) q[3];
sx q[3];
rz(1.9999947) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
