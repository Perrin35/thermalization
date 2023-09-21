OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5443213) q[0];
sx q[0];
rz(-2.7379524) q[0];
sx q[0];
rz(0.37024745) q[0];
rz(-2.9397842) q[1];
sx q[1];
rz(-1.8887853) q[1];
sx q[1];
rz(-1.7226146) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50492935) q[0];
sx q[0];
rz(-1.3267656) q[0];
sx q[0];
rz(1.7457477) q[0];
x q[1];
rz(0.51214829) q[2];
sx q[2];
rz(-1.8036267) q[2];
sx q[2];
rz(-1.6793959) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.943489) q[1];
sx q[1];
rz(-1.5183105) q[1];
sx q[1];
rz(1.0204888) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8098104) q[3];
sx q[3];
rz(-0.79091573) q[3];
sx q[3];
rz(-0.055981759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3966763) q[2];
sx q[2];
rz(-1.6884721) q[2];
sx q[2];
rz(-2.7837226) q[2];
rz(0.19168028) q[3];
sx q[3];
rz(-2.7087757) q[3];
sx q[3];
rz(-0.55309692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1382004) q[0];
sx q[0];
rz(-1.9906185) q[0];
sx q[0];
rz(-3.0733118) q[0];
rz(1.0066907) q[1];
sx q[1];
rz(-3.0196562) q[1];
sx q[1];
rz(0.65111792) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40819528) q[0];
sx q[0];
rz(-1.7719643) q[0];
sx q[0];
rz(0.64543076) q[0];
x q[1];
rz(-1.5468555) q[2];
sx q[2];
rz(-1.0451473) q[2];
sx q[2];
rz(2.0999694) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.843833) q[1];
sx q[1];
rz(-0.69121541) q[1];
sx q[1];
rz(-0.83696951) q[1];
rz(2.6504374) q[3];
sx q[3];
rz(-1.7633811) q[3];
sx q[3];
rz(2.0302949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6887001) q[2];
sx q[2];
rz(-2.9850027) q[2];
sx q[2];
rz(2.2382656) q[2];
rz(-2.3870758) q[3];
sx q[3];
rz(-1.4947596) q[3];
sx q[3];
rz(0.27404684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84905255) q[0];
sx q[0];
rz(-0.64340574) q[0];
sx q[0];
rz(-2.6480411) q[0];
rz(-1.6429398) q[1];
sx q[1];
rz(-2.7354) q[1];
sx q[1];
rz(-2.1123871) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.878359) q[0];
sx q[0];
rz(-1.9089111) q[0];
sx q[0];
rz(-2.6792206) q[0];
x q[1];
rz(2.7846787) q[2];
sx q[2];
rz(-1.9682353) q[2];
sx q[2];
rz(2.7473161) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5268847) q[1];
sx q[1];
rz(-1.7803368) q[1];
sx q[1];
rz(-2.8309612) q[1];
rz(-pi) q[2];
rz(-2.1152705) q[3];
sx q[3];
rz(-1.4635651) q[3];
sx q[3];
rz(-2.2162007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9323953) q[2];
sx q[2];
rz(-0.63964996) q[2];
sx q[2];
rz(1.2970682) q[2];
rz(0.57859892) q[3];
sx q[3];
rz(-1.9208627) q[3];
sx q[3];
rz(0.67563081) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2066752) q[0];
sx q[0];
rz(-2.4949555) q[0];
sx q[0];
rz(2.7329965) q[0];
rz(-1.714255) q[1];
sx q[1];
rz(-1.6427549) q[1];
sx q[1];
rz(-2.2600007) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13215412) q[0];
sx q[0];
rz(-2.3883925) q[0];
sx q[0];
rz(2.8892345) q[0];
rz(-2.2927631) q[2];
sx q[2];
rz(-1.0366882) q[2];
sx q[2];
rz(1.3918849) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2243005) q[1];
sx q[1];
rz(-2.2227193) q[1];
sx q[1];
rz(-0.93445458) q[1];
x q[2];
rz(-0.14933417) q[3];
sx q[3];
rz(-2.7200326) q[3];
sx q[3];
rz(-1.4967407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0980229) q[2];
sx q[2];
rz(-1.4901525) q[2];
sx q[2];
rz(2.0289452) q[2];
rz(2.5860795) q[3];
sx q[3];
rz(-1.7881309) q[3];
sx q[3];
rz(1.5295193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9920138) q[0];
sx q[0];
rz(-1.1554138) q[0];
sx q[0];
rz(-1.3647112) q[0];
rz(0.31750202) q[1];
sx q[1];
rz(-0.96173871) q[1];
sx q[1];
rz(-3.024335) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7754585) q[0];
sx q[0];
rz(-1.701208) q[0];
sx q[0];
rz(0.58701697) q[0];
rz(-pi) q[1];
rz(-1.7848894) q[2];
sx q[2];
rz(-0.97852409) q[2];
sx q[2];
rz(-1.585373) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.77814233) q[1];
sx q[1];
rz(-0.25097978) q[1];
sx q[1];
rz(-0.66448934) q[1];
rz(-pi) q[2];
rz(2.1579671) q[3];
sx q[3];
rz(-1.6348602) q[3];
sx q[3];
rz(-0.6928882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7964898) q[2];
sx q[2];
rz(-1.779665) q[2];
sx q[2];
rz(1.7618746) q[2];
rz(1.951925) q[3];
sx q[3];
rz(-0.15933557) q[3];
sx q[3];
rz(-0.074507944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-0.18774408) q[0];
sx q[0];
rz(-1.501361) q[0];
sx q[0];
rz(0.70621079) q[0];
rz(2.0300991) q[1];
sx q[1];
rz(-2.5128384) q[1];
sx q[1];
rz(3.0336753) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6861434) q[0];
sx q[0];
rz(-0.77308547) q[0];
sx q[0];
rz(-2.7129052) q[0];
rz(-1.5594257) q[2];
sx q[2];
rz(-1.5307384) q[2];
sx q[2];
rz(-2.1795321) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.94090998) q[1];
sx q[1];
rz(-2.4001277) q[1];
sx q[1];
rz(-2.7156668) q[1];
rz(-pi) q[2];
rz(2.2599254) q[3];
sx q[3];
rz(-1.1881184) q[3];
sx q[3];
rz(-1.6430287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.32249054) q[2];
sx q[2];
rz(-2.1671447) q[2];
sx q[2];
rz(-2.0325913) q[2];
rz(1.8479944) q[3];
sx q[3];
rz(-1.785991) q[3];
sx q[3];
rz(0.10722815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7727707) q[0];
sx q[0];
rz(-1.6596376) q[0];
sx q[0];
rz(3.1266881) q[0];
rz(-0.42124721) q[1];
sx q[1];
rz(-2.0887471) q[1];
sx q[1];
rz(-0.79963911) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22916238) q[0];
sx q[0];
rz(-1.8254571) q[0];
sx q[0];
rz(2.4486662) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8268379) q[2];
sx q[2];
rz(-1.2999279) q[2];
sx q[2];
rz(2.1115007) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.15496847) q[1];
sx q[1];
rz(-1.6264377) q[1];
sx q[1];
rz(-0.14543047) q[1];
rz(-1.3560489) q[3];
sx q[3];
rz(-0.575892) q[3];
sx q[3];
rz(-1.4240595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.556276) q[2];
sx q[2];
rz(-2.5053146) q[2];
sx q[2];
rz(-0.9220534) q[2];
rz(1.8317892) q[3];
sx q[3];
rz(-1.2049048) q[3];
sx q[3];
rz(-0.95782763) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9668982) q[0];
sx q[0];
rz(-1.0555203) q[0];
sx q[0];
rz(-2.877537) q[0];
rz(-1.7407725) q[1];
sx q[1];
rz(-1.4512647) q[1];
sx q[1];
rz(2.82428) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6415629) q[0];
sx q[0];
rz(-0.39425685) q[0];
sx q[0];
rz(-1.7612329) q[0];
rz(-pi) q[1];
rz(-0.1605026) q[2];
sx q[2];
rz(-1.2264894) q[2];
sx q[2];
rz(-0.80489327) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.322284) q[1];
sx q[1];
rz(-1.6849736) q[1];
sx q[1];
rz(1.6627922) q[1];
x q[2];
rz(-0.21934261) q[3];
sx q[3];
rz(-0.6570802) q[3];
sx q[3];
rz(-2.5430162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.43549609) q[2];
sx q[2];
rz(-0.93402445) q[2];
sx q[2];
rz(-1.6395817) q[2];
rz(2.8912985) q[3];
sx q[3];
rz(-1.7264265) q[3];
sx q[3];
rz(-1.9992453) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89649993) q[0];
sx q[0];
rz(-0.30830202) q[0];
sx q[0];
rz(-1.4244351) q[0];
rz(-2.6622488) q[1];
sx q[1];
rz(-1.4762676) q[1];
sx q[1];
rz(0.11553484) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4452267) q[0];
sx q[0];
rz(-1.0123024) q[0];
sx q[0];
rz(-0.88748705) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.3165452) q[2];
sx q[2];
rz(-2.4015421) q[2];
sx q[2];
rz(-2.9209602) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.1852473) q[1];
sx q[1];
rz(-2.6223409) q[1];
sx q[1];
rz(-2.1812056) q[1];
x q[2];
rz(-2.4453836) q[3];
sx q[3];
rz(-1.4783825) q[3];
sx q[3];
rz(2.948451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.66403786) q[2];
sx q[2];
rz(-2.0319735) q[2];
sx q[2];
rz(1.4257365) q[2];
rz(1.7410295) q[3];
sx q[3];
rz(-1.0481342) q[3];
sx q[3];
rz(0.28276309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7857159) q[0];
sx q[0];
rz(-0.73369217) q[0];
sx q[0];
rz(2.8724331) q[0];
rz(1.0956988) q[1];
sx q[1];
rz(-2.2311189) q[1];
sx q[1];
rz(-1.7620618) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7786498) q[0];
sx q[0];
rz(-1.4654219) q[0];
sx q[0];
rz(-1.4044579) q[0];
x q[1];
rz(1.3921521) q[2];
sx q[2];
rz(-2.7182455) q[2];
sx q[2];
rz(-1.4942126) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5386242) q[1];
sx q[1];
rz(-1.6011964) q[1];
sx q[1];
rz(-1.3974742) q[1];
rz(-pi) q[2];
x q[2];
rz(0.33196253) q[3];
sx q[3];
rz(-1.276607) q[3];
sx q[3];
rz(-2.4320137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.049008869) q[2];
sx q[2];
rz(-0.2139341) q[2];
sx q[2];
rz(-1.6258378) q[2];
rz(-1.1670636) q[3];
sx q[3];
rz(-1.556282) q[3];
sx q[3];
rz(1.0313755) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9243069) q[0];
sx q[0];
rz(-1.3047682) q[0];
sx q[0];
rz(-2.7346942) q[0];
rz(-0.45463195) q[1];
sx q[1];
rz(-1.1063207) q[1];
sx q[1];
rz(2.8938821) q[1];
rz(1.2417912) q[2];
sx q[2];
rz(-0.49578373) q[2];
sx q[2];
rz(1.0052581) q[2];
rz(1.8299673) q[3];
sx q[3];
rz(-1.6298686) q[3];
sx q[3];
rz(2.009404) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];