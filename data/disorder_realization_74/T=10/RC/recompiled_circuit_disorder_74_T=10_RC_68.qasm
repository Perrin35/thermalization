OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.4607853) q[0];
sx q[0];
rz(-2.1587125) q[0];
sx q[0];
rz(2.13184) q[0];
rz(2.9653964) q[1];
sx q[1];
rz(5.3806452) q[1];
sx q[1];
rz(7.5723958) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9292944) q[0];
sx q[0];
rz(-2.0295143) q[0];
sx q[0];
rz(-0.32231583) q[0];
rz(-1.5640366) q[2];
sx q[2];
rz(-1.5464371) q[2];
sx q[2];
rz(1.7560132) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2002441) q[1];
sx q[1];
rz(-1.779088) q[1];
sx q[1];
rz(1.3079446) q[1];
rz(-pi) q[2];
rz(3.0715277) q[3];
sx q[3];
rz(-1.3032039) q[3];
sx q[3];
rz(-2.23578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5499128) q[2];
sx q[2];
rz(-2.9788571) q[2];
sx q[2];
rz(-3.0060449) q[2];
rz(0.1401976) q[3];
sx q[3];
rz(-2.086816) q[3];
sx q[3];
rz(-0.31301096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41483375) q[0];
sx q[0];
rz(-2.2048075) q[0];
sx q[0];
rz(-0.78080368) q[0];
rz(2.8813598) q[1];
sx q[1];
rz(-0.5865016) q[1];
sx q[1];
rz(-1.3134726) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.27539) q[0];
sx q[0];
rz(-1.7821527) q[0];
sx q[0];
rz(-0.22558055) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0275231) q[2];
sx q[2];
rz(-0.55210219) q[2];
sx q[2];
rz(0.13763025) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4310415) q[1];
sx q[1];
rz(-1.4841054) q[1];
sx q[1];
rz(2.980742) q[1];
rz(-pi) q[2];
rz(-3.0659862) q[3];
sx q[3];
rz(-1.4652271) q[3];
sx q[3];
rz(2.3188863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.69497067) q[2];
sx q[2];
rz(-1.5335252) q[2];
sx q[2];
rz(2.1874645) q[2];
rz(-1.702884) q[3];
sx q[3];
rz(-1.3043159) q[3];
sx q[3];
rz(-0.71189705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0838098) q[0];
sx q[0];
rz(-0.6018146) q[0];
sx q[0];
rz(-2.4457248) q[0];
rz(-0.35481915) q[1];
sx q[1];
rz(-1.4777007) q[1];
sx q[1];
rz(2.9755039) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8649193) q[0];
sx q[0];
rz(-1.858798) q[0];
sx q[0];
rz(-0.15171451) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0679929) q[2];
sx q[2];
rz(-2.1192126) q[2];
sx q[2];
rz(1.3082248) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7754164) q[1];
sx q[1];
rz(-2.011538) q[1];
sx q[1];
rz(-0.52904769) q[1];
x q[2];
rz(1.8649678) q[3];
sx q[3];
rz(-2.4066381) q[3];
sx q[3];
rz(2.8923349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7819536) q[2];
sx q[2];
rz(-0.91527462) q[2];
sx q[2];
rz(-2.1598143) q[2];
rz(-2.0180457) q[3];
sx q[3];
rz(-2.8765364) q[3];
sx q[3];
rz(1.8410929) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.05900255) q[0];
sx q[0];
rz(-0.95472097) q[0];
sx q[0];
rz(2.847805) q[0];
rz(-2.7000973) q[1];
sx q[1];
rz(-1.0499294) q[1];
sx q[1];
rz(-2.3667483) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7168717) q[0];
sx q[0];
rz(-1.4972591) q[0];
sx q[0];
rz(1.5543235) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5261126) q[2];
sx q[2];
rz(-1.2031021) q[2];
sx q[2];
rz(1.1553264) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.47961071) q[1];
sx q[1];
rz(-1.6488711) q[1];
sx q[1];
rz(-1.9740723) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1500548) q[3];
sx q[3];
rz(-2.8354128) q[3];
sx q[3];
rz(-1.7508208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.0094770771) q[2];
sx q[2];
rz(-1.0813794) q[2];
sx q[2];
rz(2.5715128) q[2];
rz(-1.8360957) q[3];
sx q[3];
rz(-0.92731849) q[3];
sx q[3];
rz(-1.501804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.775979) q[0];
sx q[0];
rz(-1.403084) q[0];
sx q[0];
rz(-0.18138012) q[0];
rz(2.5326305) q[1];
sx q[1];
rz(-2.6700171) q[1];
sx q[1];
rz(-0.10770527) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62892249) q[0];
sx q[0];
rz(-1.9229917) q[0];
sx q[0];
rz(0.098350071) q[0];
rz(-pi) q[1];
rz(0.88476752) q[2];
sx q[2];
rz(-0.77026412) q[2];
sx q[2];
rz(-1.9212854) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9039893) q[1];
sx q[1];
rz(-0.9346107) q[1];
sx q[1];
rz(0.49828766) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0531353) q[3];
sx q[3];
rz(-0.62873757) q[3];
sx q[3];
rz(0.79013463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.057377664) q[2];
sx q[2];
rz(-1.6485018) q[2];
sx q[2];
rz(0.28953141) q[2];
rz(0.036751898) q[3];
sx q[3];
rz(-2.3993902) q[3];
sx q[3];
rz(-3.1307722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0832131) q[0];
sx q[0];
rz(-0.20861861) q[0];
sx q[0];
rz(-1.5104729) q[0];
rz(-2.0026813) q[1];
sx q[1];
rz(-1.4803671) q[1];
sx q[1];
rz(-2.1246134) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3973629) q[0];
sx q[0];
rz(-1.1083535) q[0];
sx q[0];
rz(0.37325333) q[0];
x q[1];
rz(-2.9404042) q[2];
sx q[2];
rz(-1.5941465) q[2];
sx q[2];
rz(-0.77229283) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0000227) q[1];
sx q[1];
rz(-2.0380028) q[1];
sx q[1];
rz(1.0228553) q[1];
rz(-pi) q[2];
rz(-2.6050657) q[3];
sx q[3];
rz(-1.121472) q[3];
sx q[3];
rz(-0.41538737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5974474) q[2];
sx q[2];
rz(-2.1400673) q[2];
sx q[2];
rz(2.5069359) q[2];
rz(-2.0641816) q[3];
sx q[3];
rz(-2.1118739) q[3];
sx q[3];
rz(2.6950148) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25062659) q[0];
sx q[0];
rz(-2.4996596) q[0];
sx q[0];
rz(2.2977258) q[0];
rz(1.5232874) q[1];
sx q[1];
rz(-2.4089676) q[1];
sx q[1];
rz(-1.2197781) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49918136) q[0];
sx q[0];
rz(-1.6089604) q[0];
sx q[0];
rz(3.0495318) q[0];
rz(-1.9856521) q[2];
sx q[2];
rz(-1.2659237) q[2];
sx q[2];
rz(2.4883303) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9417604) q[1];
sx q[1];
rz(-0.70711771) q[1];
sx q[1];
rz(1.9553528) q[1];
rz(-pi) q[2];
rz(1.7759833) q[3];
sx q[3];
rz(-0.41654166) q[3];
sx q[3];
rz(-0.52683559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5119778) q[2];
sx q[2];
rz(-2.2045366) q[2];
sx q[2];
rz(0.31663319) q[2];
rz(0.037242446) q[3];
sx q[3];
rz(-1.3469632) q[3];
sx q[3];
rz(2.3855456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.0379631) q[0];
sx q[0];
rz(-2.0541971) q[0];
sx q[0];
rz(-0.33139247) q[0];
rz(1.9695075) q[1];
sx q[1];
rz(-0.66292271) q[1];
sx q[1];
rz(-0.40922871) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0497826) q[0];
sx q[0];
rz(-0.70735332) q[0];
sx q[0];
rz(1.1298256) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8843083) q[2];
sx q[2];
rz(-2.0743309) q[2];
sx q[2];
rz(-1.7538278) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.85247773) q[1];
sx q[1];
rz(-0.58809885) q[1];
sx q[1];
rz(0.27681338) q[1];
rz(-pi) q[2];
rz(-1.3274566) q[3];
sx q[3];
rz(-2.0753324) q[3];
sx q[3];
rz(-0.53989053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1428712) q[2];
sx q[2];
rz(-1.2850392) q[2];
sx q[2];
rz(-0.55465737) q[2];
rz(2.5000642) q[3];
sx q[3];
rz(-2.9252164) q[3];
sx q[3];
rz(0.085111246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(0.91759578) q[0];
sx q[0];
rz(-1.6308835) q[0];
sx q[0];
rz(3.1316485) q[0];
rz(1.0558646) q[1];
sx q[1];
rz(-1.5751585) q[1];
sx q[1];
rz(-1.508629) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5976534) q[0];
sx q[0];
rz(-2.2542852) q[0];
sx q[0];
rz(-0.20792122) q[0];
rz(1.6998859) q[2];
sx q[2];
rz(-0.61868762) q[2];
sx q[2];
rz(2.8512851) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.846261) q[1];
sx q[1];
rz(-1.2519072) q[1];
sx q[1];
rz(1.0237414) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6375332) q[3];
sx q[3];
rz(-2.4373694) q[3];
sx q[3];
rz(-0.080554068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9987954) q[2];
sx q[2];
rz(-1.5974416) q[2];
sx q[2];
rz(-0.43668401) q[2];
rz(1.3302594) q[3];
sx q[3];
rz(-0.67970777) q[3];
sx q[3];
rz(1.0288303) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0325539) q[0];
sx q[0];
rz(-2.3333896) q[0];
sx q[0];
rz(0.56525266) q[0];
rz(2.8887707) q[1];
sx q[1];
rz(-1.0193453) q[1];
sx q[1];
rz(-1.3814829) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5436343) q[0];
sx q[0];
rz(-1.6126452) q[0];
sx q[0];
rz(2.6393487) q[0];
rz(-pi) q[1];
rz(-1.6818468) q[2];
sx q[2];
rz(-1.0990708) q[2];
sx q[2];
rz(1.9327088) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6227464) q[1];
sx q[1];
rz(-1.985605) q[1];
sx q[1];
rz(-0.82621375) q[1];
x q[2];
rz(2.284427) q[3];
sx q[3];
rz(-1.790159) q[3];
sx q[3];
rz(0.21881783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6910203) q[2];
sx q[2];
rz(-0.82442966) q[2];
sx q[2];
rz(2.5449469) q[2];
rz(0.48554844) q[3];
sx q[3];
rz(-0.89533007) q[3];
sx q[3];
rz(0.027464494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(0.52994603) q[0];
sx q[0];
rz(-1.9367138) q[0];
sx q[0];
rz(1.1116897) q[0];
rz(-1.8687517) q[1];
sx q[1];
rz(-1.0796937) q[1];
sx q[1];
rz(-0.83273522) q[1];
rz(3.0957072) q[2];
sx q[2];
rz(-1.8250209) q[2];
sx q[2];
rz(1.4179306) q[2];
rz(-1.4894555) q[3];
sx q[3];
rz(-2.079439) q[3];
sx q[3];
rz(0.54872201) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
