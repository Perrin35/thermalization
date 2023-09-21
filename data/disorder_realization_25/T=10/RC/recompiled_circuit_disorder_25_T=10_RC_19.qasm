OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.0619573) q[0];
sx q[0];
rz(-0.24793967) q[0];
sx q[0];
rz(2.3214582) q[0];
rz(0.040854383) q[1];
sx q[1];
rz(3.9305384) q[1];
sx q[1];
rz(9.5126704) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1782921) q[0];
sx q[0];
rz(-2.0110197) q[0];
sx q[0];
rz(-0.51142366) q[0];
rz(-pi) q[1];
rz(0.42638875) q[2];
sx q[2];
rz(-2.0585367) q[2];
sx q[2];
rz(-2.3418155) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3275798) q[1];
sx q[1];
rz(-2.6373133) q[1];
sx q[1];
rz(-0.0045279702) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5486693) q[3];
sx q[3];
rz(-2.8828354) q[3];
sx q[3];
rz(-2.4165137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1408046) q[2];
sx q[2];
rz(-1.452383) q[2];
sx q[2];
rz(-2.5722356) q[2];
rz(-1.5287483) q[3];
sx q[3];
rz(-2.5053535) q[3];
sx q[3];
rz(1.3585842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.59812087) q[0];
sx q[0];
rz(-2.9765029) q[0];
sx q[0];
rz(-2.5879481) q[0];
rz(1.2373295) q[1];
sx q[1];
rz(-1.7747223) q[1];
sx q[1];
rz(-1.9083317) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9995995) q[0];
sx q[0];
rz(-2.2337231) q[0];
sx q[0];
rz(-0.87447383) q[0];
rz(-pi) q[1];
rz(2.4670062) q[2];
sx q[2];
rz(-1.6709176) q[2];
sx q[2];
rz(-1.9763725) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.711449) q[1];
sx q[1];
rz(-1.5866382) q[1];
sx q[1];
rz(1.7769017) q[1];
x q[2];
rz(1.8096829) q[3];
sx q[3];
rz(-2.2242862) q[3];
sx q[3];
rz(-1.4373923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1374986) q[2];
sx q[2];
rz(-1.4592905) q[2];
sx q[2];
rz(-3.1393576) q[2];
rz(-0.8301174) q[3];
sx q[3];
rz(-2.3486962) q[3];
sx q[3];
rz(-1.2494276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24901351) q[0];
sx q[0];
rz(-0.61616388) q[0];
sx q[0];
rz(-0.85154831) q[0];
rz(-2.3705204) q[1];
sx q[1];
rz(-1.675019) q[1];
sx q[1];
rz(-0.071203701) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93807968) q[0];
sx q[0];
rz(-2.8305614) q[0];
sx q[0];
rz(2.2631893) q[0];
rz(-2.0729012) q[2];
sx q[2];
rz(-0.50902589) q[2];
sx q[2];
rz(-2.6976762) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5350266) q[1];
sx q[1];
rz(-0.92381322) q[1];
sx q[1];
rz(-1.4024629) q[1];
rz(2.2749388) q[3];
sx q[3];
rz(-0.94368499) q[3];
sx q[3];
rz(1.8303527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8292134) q[2];
sx q[2];
rz(-2.2177314) q[2];
sx q[2];
rz(2.6181347) q[2];
rz(2.8097025) q[3];
sx q[3];
rz(-3.0811946) q[3];
sx q[3];
rz(-0.50326842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4098542) q[0];
sx q[0];
rz(-0.24895746) q[0];
sx q[0];
rz(2.3160034) q[0];
rz(1.1666974) q[1];
sx q[1];
rz(-0.86667934) q[1];
sx q[1];
rz(1.4473787) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71890812) q[0];
sx q[0];
rz(-1.3410765) q[0];
sx q[0];
rz(2.9740888) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0091323098) q[2];
sx q[2];
rz(-1.3924686) q[2];
sx q[2];
rz(-2.669131) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4059658) q[1];
sx q[1];
rz(-1.7298797) q[1];
sx q[1];
rz(-0.87079485) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7778266) q[3];
sx q[3];
rz(-1.7978661) q[3];
sx q[3];
rz(2.098339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6772785) q[2];
sx q[2];
rz(-1.7769287) q[2];
sx q[2];
rz(-2.9966667) q[2];
rz(1.011301) q[3];
sx q[3];
rz(-0.62711182) q[3];
sx q[3];
rz(0.01468006) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1458364) q[0];
sx q[0];
rz(-0.053481426) q[0];
sx q[0];
rz(-0.68403912) q[0];
rz(-1.9513291) q[1];
sx q[1];
rz(-2.4763156) q[1];
sx q[1];
rz(1.9285944) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1367462) q[0];
sx q[0];
rz(-0.46127013) q[0];
sx q[0];
rz(0.61704163) q[0];
rz(-0.62972516) q[2];
sx q[2];
rz(-2.1246398) q[2];
sx q[2];
rz(0.50309203) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2781093) q[1];
sx q[1];
rz(-0.28243318) q[1];
sx q[1];
rz(0.43538283) q[1];
rz(-pi) q[2];
rz(-2.9641987) q[3];
sx q[3];
rz(-1.8284441) q[3];
sx q[3];
rz(-2.859476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6325536) q[2];
sx q[2];
rz(-1.8703439) q[2];
sx q[2];
rz(-0.63009134) q[2];
rz(0.57224327) q[3];
sx q[3];
rz(-2.4961491) q[3];
sx q[3];
rz(-2.2563289) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1155788) q[0];
sx q[0];
rz(-0.84340874) q[0];
sx q[0];
rz(2.8552326) q[0];
rz(-0.59965602) q[1];
sx q[1];
rz(-2.0136132) q[1];
sx q[1];
rz(2.5568331) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.023914) q[0];
sx q[0];
rz(-2.4524134) q[0];
sx q[0];
rz(0.50424772) q[0];
x q[1];
rz(-0.61442394) q[2];
sx q[2];
rz(-1.5493869) q[2];
sx q[2];
rz(1.1504088) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6962291) q[1];
sx q[1];
rz(-0.74659691) q[1];
sx q[1];
rz(-2.0018875) q[1];
rz(-pi) q[2];
rz(0.29295178) q[3];
sx q[3];
rz(-2.9679899) q[3];
sx q[3];
rz(1.1374744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9770603) q[2];
sx q[2];
rz(-0.88447276) q[2];
sx q[2];
rz(-1.6284846) q[2];
rz(-2.1785054) q[3];
sx q[3];
rz(-1.3685127) q[3];
sx q[3];
rz(-0.62455463) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0528089) q[0];
sx q[0];
rz(-1.4316906) q[0];
sx q[0];
rz(0.54876304) q[0];
rz(-3.0896297) q[1];
sx q[1];
rz(-2.6497662) q[1];
sx q[1];
rz(-0.60639492) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.711395) q[0];
sx q[0];
rz(-2.093962) q[0];
sx q[0];
rz(0.27336143) q[0];
rz(-pi) q[1];
rz(2.7658471) q[2];
sx q[2];
rz(-2.1990015) q[2];
sx q[2];
rz(1.9919765) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7402621) q[1];
sx q[1];
rz(-1.4583424) q[1];
sx q[1];
rz(-0.80660352) q[1];
x q[2];
rz(-2.2783979) q[3];
sx q[3];
rz(-2.6226225) q[3];
sx q[3];
rz(0.88245813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9397883) q[2];
sx q[2];
rz(-2.1494892) q[2];
sx q[2];
rz(-0.50841224) q[2];
rz(1.1172179) q[3];
sx q[3];
rz(-2.9682187) q[3];
sx q[3];
rz(-1.6569051) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8023119) q[0];
sx q[0];
rz(-2.7957714) q[0];
sx q[0];
rz(2.3876277) q[0];
rz(2.1915961) q[1];
sx q[1];
rz(-1.4818622) q[1];
sx q[1];
rz(-0.22769134) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84275093) q[0];
sx q[0];
rz(-1.9831053) q[0];
sx q[0];
rz(0.077667872) q[0];
rz(0.67543244) q[2];
sx q[2];
rz(-2.5917412) q[2];
sx q[2];
rz(-0.39381248) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.258192) q[1];
sx q[1];
rz(-2.0996143) q[1];
sx q[1];
rz(1.3099111) q[1];
rz(0.55898198) q[3];
sx q[3];
rz(-1.535927) q[3];
sx q[3];
rz(2.1502286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.078538744) q[2];
sx q[2];
rz(-2.9324014) q[2];
sx q[2];
rz(-2.3030346) q[2];
rz(-2.7770384) q[3];
sx q[3];
rz(-1.3804599) q[3];
sx q[3];
rz(0.84028876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39695981) q[0];
sx q[0];
rz(-1.7663706) q[0];
sx q[0];
rz(2.3378085) q[0];
rz(-2.0869758) q[1];
sx q[1];
rz(-2.5024253) q[1];
sx q[1];
rz(1.9931591) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52590695) q[0];
sx q[0];
rz(-1.4988168) q[0];
sx q[0];
rz(-0.023948897) q[0];
rz(-pi) q[1];
rz(2.7788413) q[2];
sx q[2];
rz(-1.0748378) q[2];
sx q[2];
rz(-0.19778684) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.941274) q[1];
sx q[1];
rz(-0.66220821) q[1];
sx q[1];
rz(-2.3108285) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5356482) q[3];
sx q[3];
rz(-1.2190378) q[3];
sx q[3];
rz(-2.7621321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.049872963) q[2];
sx q[2];
rz(-2.1240978) q[2];
sx q[2];
rz(0.77825528) q[2];
rz(2.9459279) q[3];
sx q[3];
rz(-1.0396495) q[3];
sx q[3];
rz(-1.0218609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8675999) q[0];
sx q[0];
rz(-0.99853981) q[0];
sx q[0];
rz(-2.8883949) q[0];
rz(-0.7397488) q[1];
sx q[1];
rz(-2.4052129) q[1];
sx q[1];
rz(-2.443312) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4822599) q[0];
sx q[0];
rz(-1.8730622) q[0];
sx q[0];
rz(1.4955826) q[0];
rz(-pi) q[1];
rz(1.4966428) q[2];
sx q[2];
rz(-1.5955131) q[2];
sx q[2];
rz(0.93862426) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.56863927) q[1];
sx q[1];
rz(-2.0268974) q[1];
sx q[1];
rz(-0.24104636) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2020338) q[3];
sx q[3];
rz(-0.59951111) q[3];
sx q[3];
rz(-0.60839073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0239821) q[2];
sx q[2];
rz(-0.88279804) q[2];
sx q[2];
rz(2.3013766) q[2];
rz(0.64030567) q[3];
sx q[3];
rz(-2.2655723) q[3];
sx q[3];
rz(-1.9742112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2587851) q[0];
sx q[0];
rz(-0.85492815) q[0];
sx q[0];
rz(1.7229765) q[0];
rz(0.029126833) q[1];
sx q[1];
rz(-0.11321414) q[1];
sx q[1];
rz(-1.3197457) q[1];
rz(2.6559033) q[2];
sx q[2];
rz(-2.4808241) q[2];
sx q[2];
rz(-0.11447699) q[2];
rz(0.15016951) q[3];
sx q[3];
rz(-2.2554845) q[3];
sx q[3];
rz(3.1156202) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];