OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7862608) q[0];
sx q[0];
rz(-0.064602764) q[0];
sx q[0];
rz(0.021615418) q[0];
rz(2.1463483) q[1];
sx q[1];
rz(-1.8145476) q[1];
sx q[1];
rz(-1.8099161) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3992213) q[0];
sx q[0];
rz(-2.135015) q[0];
sx q[0];
rz(-2.7447634) q[0];
rz(-pi) q[1];
rz(0.10279074) q[2];
sx q[2];
rz(-1.3445026) q[2];
sx q[2];
rz(3.0068827) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5035489) q[1];
sx q[1];
rz(-0.81144864) q[1];
sx q[1];
rz(-1.391295) q[1];
rz(-pi) q[2];
rz(2.0256151) q[3];
sx q[3];
rz(-2.4664719) q[3];
sx q[3];
rz(1.7636553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.33064476) q[2];
sx q[2];
rz(-1.6487048) q[2];
sx q[2];
rz(-0.60418207) q[2];
rz(-2.1172681) q[3];
sx q[3];
rz(-1.9842792) q[3];
sx q[3];
rz(1.1272875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8169096) q[0];
sx q[0];
rz(-0.01318251) q[0];
sx q[0];
rz(-2.0781793) q[0];
rz(0.88513199) q[1];
sx q[1];
rz(-1.5848031) q[1];
sx q[1];
rz(-0.0016454776) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5981818) q[0];
sx q[0];
rz(-3.0953005) q[0];
sx q[0];
rz(1.8020736) q[0];
rz(-pi) q[1];
rz(2.8029289) q[2];
sx q[2];
rz(-0.67064697) q[2];
sx q[2];
rz(-1.8530958) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9284294) q[1];
sx q[1];
rz(-2.3191116) q[1];
sx q[1];
rz(2.806042) q[1];
x q[2];
rz(0.9886338) q[3];
sx q[3];
rz(-0.30246624) q[3];
sx q[3];
rz(-2.8305588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1559747) q[2];
sx q[2];
rz(-1.5094455) q[2];
sx q[2];
rz(-2.4334811) q[2];
rz(-1.0937141) q[3];
sx q[3];
rz(-1.8023068) q[3];
sx q[3];
rz(-2.1480907) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.920632) q[0];
sx q[0];
rz(-1.4039803) q[0];
sx q[0];
rz(-0.8272585) q[0];
rz(-0.0050841252) q[1];
sx q[1];
rz(-1.2132443) q[1];
sx q[1];
rz(-2.0522096) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1339061) q[0];
sx q[0];
rz(-2.0984762) q[0];
sx q[0];
rz(-0.49467996) q[0];
x q[1];
rz(-1.1459648) q[2];
sx q[2];
rz(-2.1263188) q[2];
sx q[2];
rz(-0.59031634) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8100064) q[1];
sx q[1];
rz(-2.7467105) q[1];
sx q[1];
rz(-1.3473131) q[1];
x q[2];
rz(0.79077625) q[3];
sx q[3];
rz(-1.1547935) q[3];
sx q[3];
rz(-0.63000597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0638782) q[2];
sx q[2];
rz(-2.2641247) q[2];
sx q[2];
rz(0.88469488) q[2];
rz(1.1832773) q[3];
sx q[3];
rz(-1.3180472) q[3];
sx q[3];
rz(-1.4900835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5692212) q[0];
sx q[0];
rz(-2.5771993) q[0];
sx q[0];
rz(0.72682056) q[0];
rz(0.80365333) q[1];
sx q[1];
rz(-2.0539961) q[1];
sx q[1];
rz(2.7817536) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7774178) q[0];
sx q[0];
rz(-0.67876498) q[0];
sx q[0];
rz(1.8666408) q[0];
rz(-2.4904576) q[2];
sx q[2];
rz(-1.6614117) q[2];
sx q[2];
rz(-0.091094253) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6872014) q[1];
sx q[1];
rz(-1.193421) q[1];
sx q[1];
rz(2.6881933) q[1];
rz(-pi) q[2];
rz(-2.8344645) q[3];
sx q[3];
rz(-1.9046475) q[3];
sx q[3];
rz(-3.0577554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.56746733) q[2];
sx q[2];
rz(-1.987477) q[2];
sx q[2];
rz(-0.7652258) q[2];
rz(0.75677538) q[3];
sx q[3];
rz(-2.5496343) q[3];
sx q[3];
rz(0.24100196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1446447) q[0];
sx q[0];
rz(-0.47676555) q[0];
sx q[0];
rz(1.7255406) q[0];
rz(-2.7744746) q[1];
sx q[1];
rz(-1.3857625) q[1];
sx q[1];
rz(-2.1062772) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1015548) q[0];
sx q[0];
rz(-2.3674175) q[0];
sx q[0];
rz(2.220962) q[0];
x q[1];
rz(-2.6501422) q[2];
sx q[2];
rz(-1.5361668) q[2];
sx q[2];
rz(-2.2143242) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.94381911) q[1];
sx q[1];
rz(-3.0316331) q[1];
sx q[1];
rz(0.52093671) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3607929) q[3];
sx q[3];
rz(-1.8404507) q[3];
sx q[3];
rz(-2.8402013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3570024) q[2];
sx q[2];
rz(-1.4214397) q[2];
sx q[2];
rz(-2.3727097) q[2];
rz(-0.33603493) q[3];
sx q[3];
rz(-2.3612645) q[3];
sx q[3];
rz(1.4679573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96200213) q[0];
sx q[0];
rz(-0.63350326) q[0];
sx q[0];
rz(-1.996421) q[0];
rz(-2.0369453) q[1];
sx q[1];
rz(-1.2303753) q[1];
sx q[1];
rz(-2.9343658) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0475537) q[0];
sx q[0];
rz(-2.5045966) q[0];
sx q[0];
rz(0.23547049) q[0];
x q[1];
rz(-2.5099498) q[2];
sx q[2];
rz(-2.4372059) q[2];
sx q[2];
rz(-2.7001675) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6991899) q[1];
sx q[1];
rz(-1.3283722) q[1];
sx q[1];
rz(1.8276023) q[1];
x q[2];
rz(2.3239922) q[3];
sx q[3];
rz(-1.451965) q[3];
sx q[3];
rz(-2.6314051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3383639) q[2];
sx q[2];
rz(-0.88023606) q[2];
sx q[2];
rz(0.72171372) q[2];
rz(-1.8668113) q[3];
sx q[3];
rz(-1.5694247) q[3];
sx q[3];
rz(0.26003626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24467829) q[0];
sx q[0];
rz(-1.5667863) q[0];
sx q[0];
rz(-0.73202837) q[0];
rz(-3.1320944) q[1];
sx q[1];
rz(-2.5932725) q[1];
sx q[1];
rz(-0.2917372) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1222629) q[0];
sx q[0];
rz(-2.7195192) q[0];
sx q[0];
rz(0.12515573) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2675915) q[2];
sx q[2];
rz(-1.6723987) q[2];
sx q[2];
rz(0.50819699) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.74624324) q[1];
sx q[1];
rz(-2.8031073) q[1];
sx q[1];
rz(-0.45792087) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8802059) q[3];
sx q[3];
rz(-2.072812) q[3];
sx q[3];
rz(1.3537784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.26414028) q[2];
sx q[2];
rz(-1.143127) q[2];
sx q[2];
rz(-2.3042802) q[2];
rz(1.9705747) q[3];
sx q[3];
rz(-1.5191017) q[3];
sx q[3];
rz(3.0795735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.0983122) q[0];
sx q[0];
rz(-2.4536112) q[0];
sx q[0];
rz(-0.6638546) q[0];
rz(3.035416) q[1];
sx q[1];
rz(-2.5352434) q[1];
sx q[1];
rz(0.95867872) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3631358) q[0];
sx q[0];
rz(-1.0582557) q[0];
sx q[0];
rz(0.97495671) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8616222) q[2];
sx q[2];
rz(-2.0915871) q[2];
sx q[2];
rz(-2.1998646) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9194591) q[1];
sx q[1];
rz(-2.6047915) q[1];
sx q[1];
rz(1.3608576) q[1];
rz(-1.4484552) q[3];
sx q[3];
rz(-2.1515176) q[3];
sx q[3];
rz(1.3337222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0083996) q[2];
sx q[2];
rz(-1.9063176) q[2];
sx q[2];
rz(-2.2224902) q[2];
rz(1.7637926) q[3];
sx q[3];
rz(-1.931124) q[3];
sx q[3];
rz(1.9581883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.33525) q[0];
sx q[0];
rz(-0.85334539) q[0];
sx q[0];
rz(2.7046955) q[0];
rz(0.70029744) q[1];
sx q[1];
rz(-1.7179787) q[1];
sx q[1];
rz(-2.6760496) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3532928) q[0];
sx q[0];
rz(-0.71209891) q[0];
sx q[0];
rz(-3.0112991) q[0];
rz(1.9829468) q[2];
sx q[2];
rz(-1.4530384) q[2];
sx q[2];
rz(-2.462537) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1411966) q[1];
sx q[1];
rz(-1.9138412) q[1];
sx q[1];
rz(-2.6617906) q[1];
rz(1.9302619) q[3];
sx q[3];
rz(-1.8065479) q[3];
sx q[3];
rz(-0.93709968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7312701) q[2];
sx q[2];
rz(-0.53933829) q[2];
sx q[2];
rz(1.643606) q[2];
rz(0.20478976) q[3];
sx q[3];
rz(-1.0054761) q[3];
sx q[3];
rz(2.8695316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4964504) q[0];
sx q[0];
rz(-2.5800939) q[0];
sx q[0];
rz(1.1219332) q[0];
rz(2.3902068) q[1];
sx q[1];
rz(-1.6151927) q[1];
sx q[1];
rz(2.5591992) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84401417) q[0];
sx q[0];
rz(-1.7739002) q[0];
sx q[0];
rz(-0.88589478) q[0];
rz(-1.8495314) q[2];
sx q[2];
rz(-2.6346452) q[2];
sx q[2];
rz(0.35001937) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0092587) q[1];
sx q[1];
rz(-1.7528105) q[1];
sx q[1];
rz(3.0218714) q[1];
x q[2];
rz(0.21721812) q[3];
sx q[3];
rz(-1.5434885) q[3];
sx q[3];
rz(0.56009968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0104684) q[2];
sx q[2];
rz(-0.98103395) q[2];
sx q[2];
rz(0.76254145) q[2];
rz(1.7307581) q[3];
sx q[3];
rz(-0.91791955) q[3];
sx q[3];
rz(-1.5677174) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0762155) q[0];
sx q[0];
rz(-2.1552754) q[0];
sx q[0];
rz(-1.7393204) q[0];
rz(-1.8021884) q[1];
sx q[1];
rz(-1.4504455) q[1];
sx q[1];
rz(1.6557678) q[1];
rz(0.81417685) q[2];
sx q[2];
rz(-1.3635175) q[2];
sx q[2];
rz(-1.2843532) q[2];
rz(-2.3253757) q[3];
sx q[3];
rz(-0.64707884) q[3];
sx q[3];
rz(-0.94221471) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
