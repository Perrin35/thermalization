OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3396575) q[0];
sx q[0];
rz(-0.8987838) q[0];
sx q[0];
rz(1.3026097) q[0];
rz(-3.0175735) q[1];
sx q[1];
rz(-1.7350585) q[1];
sx q[1];
rz(1.6279434) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2500656) q[0];
sx q[0];
rz(-0.45338777) q[0];
sx q[0];
rz(2.2946623) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3816125) q[2];
sx q[2];
rz(-0.70356762) q[2];
sx q[2];
rz(2.4576996) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.044643) q[1];
sx q[1];
rz(-1.5162807) q[1];
sx q[1];
rz(-1.2300371) q[1];
x q[2];
rz(2.7645993) q[3];
sx q[3];
rz(-0.97929472) q[3];
sx q[3];
rz(-1.5821004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0579494) q[2];
sx q[2];
rz(-1.3336072) q[2];
sx q[2];
rz(2.7570214) q[2];
rz(-2.7390076) q[3];
sx q[3];
rz(-1.7076924) q[3];
sx q[3];
rz(-2.3334077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3835417) q[0];
sx q[0];
rz(-1.533968) q[0];
sx q[0];
rz(2.5480399) q[0];
rz(-1.3213762) q[1];
sx q[1];
rz(-2.0621736) q[1];
sx q[1];
rz(0.039332565) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6932019) q[0];
sx q[0];
rz(-1.9422008) q[0];
sx q[0];
rz(-1.0813963) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.987542) q[2];
sx q[2];
rz(-2.4994435) q[2];
sx q[2];
rz(2.1352445) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6249121) q[1];
sx q[1];
rz(-2.6221226) q[1];
sx q[1];
rz(3.0416529) q[1];
x q[2];
rz(-1.5402628) q[3];
sx q[3];
rz(-2.1200006) q[3];
sx q[3];
rz(0.98992482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.84052691) q[2];
sx q[2];
rz(-1.4378005) q[2];
sx q[2];
rz(0.62704101) q[2];
rz(-0.4380694) q[3];
sx q[3];
rz(-1.4984683) q[3];
sx q[3];
rz(2.0048678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7196734) q[0];
sx q[0];
rz(-2.0543583) q[0];
sx q[0];
rz(1.7701953) q[0];
rz(-0.60945359) q[1];
sx q[1];
rz(-0.89027673) q[1];
sx q[1];
rz(-1.9166463) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.030239074) q[0];
sx q[0];
rz(-2.9269106) q[0];
sx q[0];
rz(2.6743321) q[0];
x q[1];
rz(2.248522) q[2];
sx q[2];
rz(-1.8883103) q[2];
sx q[2];
rz(1.7567321) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4546734) q[1];
sx q[1];
rz(-0.49263601) q[1];
sx q[1];
rz(-0.83322816) q[1];
rz(-1.8277934) q[3];
sx q[3];
rz(-2.42958) q[3];
sx q[3];
rz(-0.76699996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8372535) q[2];
sx q[2];
rz(-0.86696583) q[2];
sx q[2];
rz(0.85401094) q[2];
rz(0.52418661) q[3];
sx q[3];
rz(-1.5093404) q[3];
sx q[3];
rz(2.1715651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9147341) q[0];
sx q[0];
rz(-0.15758841) q[0];
sx q[0];
rz(0.7937113) q[0];
rz(0.0018250068) q[1];
sx q[1];
rz(-0.55627126) q[1];
sx q[1];
rz(0.8108286) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13769606) q[0];
sx q[0];
rz(-2.3783149) q[0];
sx q[0];
rz(1.3062817) q[0];
rz(0.89703441) q[2];
sx q[2];
rz(-1.22222) q[2];
sx q[2];
rz(-0.86212117) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6135679) q[1];
sx q[1];
rz(-0.9903637) q[1];
sx q[1];
rz(2.7356262) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.14489095) q[3];
sx q[3];
rz(-1.1456744) q[3];
sx q[3];
rz(-1.8505877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4578555) q[2];
sx q[2];
rz(-2.808414) q[2];
sx q[2];
rz(1.6443171) q[2];
rz(-1.3582683) q[3];
sx q[3];
rz(-2.0446916) q[3];
sx q[3];
rz(-1.3076521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4182279) q[0];
sx q[0];
rz(-1.4606553) q[0];
sx q[0];
rz(2.8826868) q[0];
rz(-0.12340165) q[1];
sx q[1];
rz(-2.5105208) q[1];
sx q[1];
rz(-0.48615989) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21307316) q[0];
sx q[0];
rz(-2.3384636) q[0];
sx q[0];
rz(-2.2691239) q[0];
rz(1.4347378) q[2];
sx q[2];
rz(-2.0199088) q[2];
sx q[2];
rz(-0.27829188) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7546995) q[1];
sx q[1];
rz(-0.4122977) q[1];
sx q[1];
rz(-2.0569494) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2083602) q[3];
sx q[3];
rz(-2.3321816) q[3];
sx q[3];
rz(0.903086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2316124) q[2];
sx q[2];
rz(-0.51311246) q[2];
sx q[2];
rz(1.7463589) q[2];
rz(-0.95476556) q[3];
sx q[3];
rz(-1.3937817) q[3];
sx q[3];
rz(-2.098293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-2.2024277) q[0];
sx q[0];
rz(-1.2347777) q[0];
sx q[0];
rz(2.3811316) q[0];
rz(2.3072534) q[1];
sx q[1];
rz(-1.1015247) q[1];
sx q[1];
rz(1.1044501) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6532947) q[0];
sx q[0];
rz(-2.1135931) q[0];
sx q[0];
rz(-2.8858868) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1708753) q[2];
sx q[2];
rz(-1.2515704) q[2];
sx q[2];
rz(0.17344061) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.901746) q[1];
sx q[1];
rz(-0.86565986) q[1];
sx q[1];
rz(-2.2589931) q[1];
rz(0.31046162) q[3];
sx q[3];
rz(-2.1230773) q[3];
sx q[3];
rz(2.3776835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1793648) q[2];
sx q[2];
rz(-0.70415512) q[2];
sx q[2];
rz(3.1000225) q[2];
rz(2.9465594) q[3];
sx q[3];
rz(-0.2427559) q[3];
sx q[3];
rz(-2.354505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4389909) q[0];
sx q[0];
rz(-0.85346237) q[0];
sx q[0];
rz(-0.75337291) q[0];
rz(1.0607177) q[1];
sx q[1];
rz(-0.80442387) q[1];
sx q[1];
rz(-0.20739584) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8327861) q[0];
sx q[0];
rz(-2.1784003) q[0];
sx q[0];
rz(-2.4101546) q[0];
rz(-pi) q[1];
x q[1];
rz(2.676187) q[2];
sx q[2];
rz(-1.3108746) q[2];
sx q[2];
rz(-2.9170819) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7809697) q[1];
sx q[1];
rz(-2.1163673) q[1];
sx q[1];
rz(-1.8459794) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5630521) q[3];
sx q[3];
rz(-1.5317347) q[3];
sx q[3];
rz(2.9970053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.14219001) q[2];
sx q[2];
rz(-1.157607) q[2];
sx q[2];
rz(3.0863777) q[2];
rz(-2.2986872) q[3];
sx q[3];
rz(-0.75777811) q[3];
sx q[3];
rz(2.260476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.371599) q[0];
sx q[0];
rz(-2.433625) q[0];
sx q[0];
rz(-1.1771033) q[0];
rz(-2.2616995) q[1];
sx q[1];
rz(-2.8113007) q[1];
sx q[1];
rz(-3.0807307) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.092120425) q[0];
sx q[0];
rz(-0.15639399) q[0];
sx q[0];
rz(2.2230287) q[0];
rz(1.251077) q[2];
sx q[2];
rz(-1.1100612) q[2];
sx q[2];
rz(-1.899903) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1424928) q[1];
sx q[1];
rz(-1.4542034) q[1];
sx q[1];
rz(0.57013843) q[1];
rz(-1.8350321) q[3];
sx q[3];
rz(-1.4995534) q[3];
sx q[3];
rz(0.83282214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0941144) q[2];
sx q[2];
rz(-0.55314174) q[2];
sx q[2];
rz(1.6806357) q[2];
rz(0.85584062) q[3];
sx q[3];
rz(-2.1688921) q[3];
sx q[3];
rz(-0.87876764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1020553) q[0];
sx q[0];
rz(-0.85298959) q[0];
sx q[0];
rz(3.1374186) q[0];
rz(1.2511085) q[1];
sx q[1];
rz(-3.0079542) q[1];
sx q[1];
rz(-2.4600696) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7345372) q[0];
sx q[0];
rz(-0.87293599) q[0];
sx q[0];
rz(2.7584793) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.90687554) q[2];
sx q[2];
rz(-1.0685759) q[2];
sx q[2];
rz(-2.2447667) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.2118069) q[1];
sx q[1];
rz(-0.34084596) q[1];
sx q[1];
rz(-2.5548359) q[1];
x q[2];
rz(2.8130346) q[3];
sx q[3];
rz(-0.78379455) q[3];
sx q[3];
rz(0.72600049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5847797) q[2];
sx q[2];
rz(-0.4759554) q[2];
sx q[2];
rz(1.7468096) q[2];
rz(-2.7252588) q[3];
sx q[3];
rz(-1.8463912) q[3];
sx q[3];
rz(-0.39286119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3861179) q[0];
sx q[0];
rz(-0.32805726) q[0];
sx q[0];
rz(-1.0020483) q[0];
rz(0.26456061) q[1];
sx q[1];
rz(-1.325565) q[1];
sx q[1];
rz(1.6511668) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0984983) q[0];
sx q[0];
rz(-1.7014456) q[0];
sx q[0];
rz(-2.5976973) q[0];
x q[1];
rz(1.7298431) q[2];
sx q[2];
rz(-0.76888409) q[2];
sx q[2];
rz(-1.8184219) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.24764168) q[1];
sx q[1];
rz(-0.14130768) q[1];
sx q[1];
rz(0.9402449) q[1];
rz(-pi) q[2];
rz(0.57820676) q[3];
sx q[3];
rz(-3.0319408) q[3];
sx q[3];
rz(-2.6771817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1313021) q[2];
sx q[2];
rz(-1.7346069) q[2];
sx q[2];
rz(-1.7245801) q[2];
rz(2.8085282) q[3];
sx q[3];
rz(-2.371623) q[3];
sx q[3];
rz(-1.7634348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3604816) q[0];
sx q[0];
rz(-1.8753373) q[0];
sx q[0];
rz(-1.2486096) q[0];
rz(-3.1151415) q[1];
sx q[1];
rz(-1.6117663) q[1];
sx q[1];
rz(1.5996006) q[1];
rz(-1.1556861) q[2];
sx q[2];
rz(-0.90183707) q[2];
sx q[2];
rz(3.1169008) q[2];
rz(-0.61458434) q[3];
sx q[3];
rz(-1.3086223) q[3];
sx q[3];
rz(2.3400459) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
