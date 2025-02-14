OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.46848133) q[0];
sx q[0];
rz(-0.10804478) q[0];
sx q[0];
rz(0.98969069) q[0];
rz(-1.5762848) q[1];
sx q[1];
rz(-1.58374) q[1];
sx q[1];
rz(-1.6852112) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1030758) q[0];
sx q[0];
rz(-1.4355735) q[0];
sx q[0];
rz(1.200202) q[0];
rz(1.6077432) q[2];
sx q[2];
rz(-1.4375028) q[2];
sx q[2];
rz(2.5306866) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4850581) q[1];
sx q[1];
rz(-0.63737255) q[1];
sx q[1];
rz(2.0591048) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1478808) q[3];
sx q[3];
rz(-1.607393) q[3];
sx q[3];
rz(0.59714506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9911389) q[2];
sx q[2];
rz(-0.011971124) q[2];
sx q[2];
rz(1.0452622) q[2];
rz(-0.996905) q[3];
sx q[3];
rz(-3.1361339) q[3];
sx q[3];
rz(1.8256942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.588722) q[0];
sx q[0];
rz(-1.9010239) q[0];
sx q[0];
rz(1.3528104) q[0];
rz(3.1006587) q[1];
sx q[1];
rz(-1.9238238) q[1];
sx q[1];
rz(-1.5997684) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38986015) q[0];
sx q[0];
rz(-2.9417188) q[0];
sx q[0];
rz(0.78011192) q[0];
rz(-pi) q[1];
rz(1.5495318) q[2];
sx q[2];
rz(-1.5545115) q[2];
sx q[2];
rz(0.88693888) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.62470755) q[1];
sx q[1];
rz(-1.5398281) q[1];
sx q[1];
rz(2.1672441) q[1];
x q[2];
rz(0.38741855) q[3];
sx q[3];
rz(-1.0325047) q[3];
sx q[3];
rz(-1.1673934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6359977) q[2];
sx q[2];
rz(-0.032568585) q[2];
sx q[2];
rz(2.6231498) q[2];
rz(1.1238267) q[3];
sx q[3];
rz(-2.3321407) q[3];
sx q[3];
rz(2.7037485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6464624) q[0];
sx q[0];
rz(-3.0914682) q[0];
sx q[0];
rz(1.5396402) q[0];
rz(2.4165972) q[1];
sx q[1];
rz(-3.1097737) q[1];
sx q[1];
rz(2.4620788) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2452263) q[0];
sx q[0];
rz(-0.99649444) q[0];
sx q[0];
rz(2.0942874) q[0];
rz(-pi) q[1];
rz(0.21444397) q[2];
sx q[2];
rz(-1.5753928) q[2];
sx q[2];
rz(-1.7018122) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.589349) q[1];
sx q[1];
rz(-1.8830944) q[1];
sx q[1];
rz(2.8683788) q[1];
rz(-0.26247611) q[3];
sx q[3];
rz(-2.1427689) q[3];
sx q[3];
rz(-2.9140896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9069549) q[2];
sx q[2];
rz(-2.8997771) q[2];
sx q[2];
rz(-2.6406636) q[2];
rz(0.45589724) q[3];
sx q[3];
rz(-3.1152476) q[3];
sx q[3];
rz(-0.19550368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2723715) q[0];
sx q[0];
rz(-0.048373241) q[0];
sx q[0];
rz(-2.3424171) q[0];
rz(2.8436106) q[1];
sx q[1];
rz(-2.8831392) q[1];
sx q[1];
rz(-2.2395649) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74583861) q[0];
sx q[0];
rz(-0.82358783) q[0];
sx q[0];
rz(1.7584778) q[0];
rz(-1.41961) q[2];
sx q[2];
rz(-2.107321) q[2];
sx q[2];
rz(0.46162185) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.30107547) q[1];
sx q[1];
rz(-1.4231735) q[1];
sx q[1];
rz(-3.1319437) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4933735) q[3];
sx q[3];
rz(-1.5158049) q[3];
sx q[3];
rz(-2.1907326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.67348376) q[2];
sx q[2];
rz(-0.031688422) q[2];
sx q[2];
rz(1.6710949) q[2];
rz(0.18943131) q[3];
sx q[3];
rz(-3.093284) q[3];
sx q[3];
rz(-0.76907492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9162132) q[0];
sx q[0];
rz(-3.0140641) q[0];
sx q[0];
rz(-0.39176971) q[0];
rz(1.0852934) q[1];
sx q[1];
rz(-3.1317874) q[1];
sx q[1];
rz(2.772803) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.987326) q[0];
sx q[0];
rz(-1.0079103) q[0];
sx q[0];
rz(1.994654) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1283705) q[2];
sx q[2];
rz(-0.48071024) q[2];
sx q[2];
rz(1.7465357) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6435191) q[1];
sx q[1];
rz(-3.1349224) q[1];
sx q[1];
rz(-1.4590864) q[1];
rz(-pi) q[2];
x q[2];
rz(0.7097575) q[3];
sx q[3];
rz(-1.7492948) q[3];
sx q[3];
rz(2.3861726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.58086777) q[2];
sx q[2];
rz(-3.0256425) q[2];
sx q[2];
rz(-1.3690534) q[2];
rz(3.0154058) q[3];
sx q[3];
rz(-0.43712619) q[3];
sx q[3];
rz(1.060846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5419902) q[0];
sx q[0];
rz(-2.3805711) q[0];
sx q[0];
rz(1.5498932) q[0];
rz(-2.1688993) q[1];
sx q[1];
rz(-0.21316554) q[1];
sx q[1];
rz(1.0290283) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1465107) q[0];
sx q[0];
rz(-1.8712988) q[0];
sx q[0];
rz(-0.62733688) q[0];
rz(-0.23795655) q[2];
sx q[2];
rz(-1.4334442) q[2];
sx q[2];
rz(0.99983012) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2756535) q[1];
sx q[1];
rz(-1.5703619) q[1];
sx q[1];
rz(-1.666953) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.036110445) q[3];
sx q[3];
rz(-1.7381769) q[3];
sx q[3];
rz(0.13425628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.35919967) q[2];
sx q[2];
rz(-1.7155557) q[2];
sx q[2];
rz(0.4314118) q[2];
rz(-0.99724489) q[3];
sx q[3];
rz(-3.1044208) q[3];
sx q[3];
rz(0.55716151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7998841) q[0];
sx q[0];
rz(-2.6567617) q[0];
sx q[0];
rz(1.0579911) q[0];
rz(-2.3175088) q[1];
sx q[1];
rz(-1.7015142e-05) q[1];
sx q[1];
rz(0.81738671) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0213172) q[0];
sx q[0];
rz(-0.45545721) q[0];
sx q[0];
rz(-1.0767471) q[0];
x q[1];
rz(-1.5777052) q[2];
sx q[2];
rz(-1.5890317) q[2];
sx q[2];
rz(-1.4438786) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.89725607) q[1];
sx q[1];
rz(-2.896999) q[1];
sx q[1];
rz(-1.7331428) q[1];
rz(-pi) q[2];
rz(-2.3895719) q[3];
sx q[3];
rz(-1.9746426) q[3];
sx q[3];
rz(0.93896093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.99007964) q[2];
sx q[2];
rz(-1.238287) q[2];
sx q[2];
rz(-1.6286758) q[2];
rz(0.40142909) q[3];
sx q[3];
rz(-3.1135058) q[3];
sx q[3];
rz(-0.95429558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3681188) q[0];
sx q[0];
rz(-3.0768657) q[0];
sx q[0];
rz(1.7777959) q[0];
rz(3.0996481) q[1];
sx q[1];
rz(-0.13101235) q[1];
sx q[1];
rz(2.1379474) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.100228) q[0];
sx q[0];
rz(-1.4778839) q[0];
sx q[0];
rz(-2.5071457) q[0];
rz(0.0081459002) q[2];
sx q[2];
rz(-1.2334494) q[2];
sx q[2];
rz(-1.0376736) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.81627203) q[1];
sx q[1];
rz(-2.8311816) q[1];
sx q[1];
rz(-0.85122935) q[1];
rz(-pi) q[2];
rz(0.3866638) q[3];
sx q[3];
rz(-0.69299504) q[3];
sx q[3];
rz(-0.7469783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4369138) q[2];
sx q[2];
rz(-0.059160058) q[2];
sx q[2];
rz(-1.3979647) q[2];
rz(2.6513903) q[3];
sx q[3];
rz(-0.043488113) q[3];
sx q[3];
rz(1.1787666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.101508) q[0];
sx q[0];
rz(-0.12797102) q[0];
sx q[0];
rz(-2.9301933) q[0];
rz(-1.118411) q[1];
sx q[1];
rz(-3.1299997) q[1];
sx q[1];
rz(1.6308019) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.495962) q[0];
sx q[0];
rz(-1.8798774) q[0];
sx q[0];
rz(1.8739971) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1338992) q[2];
sx q[2];
rz(-2.1497576) q[2];
sx q[2];
rz(-1.4329662) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.37247) q[1];
sx q[1];
rz(-1.5536397) q[1];
sx q[1];
rz(-2.9869798) q[1];
x q[2];
rz(-0.49064891) q[3];
sx q[3];
rz(-2.9534441) q[3];
sx q[3];
rz(0.33056459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7819034) q[2];
sx q[2];
rz(-3.1236881) q[2];
sx q[2];
rz(-2.8622799) q[2];
rz(2.277788) q[3];
sx q[3];
rz(-3.1371208) q[3];
sx q[3];
rz(2.4590676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.913468) q[0];
sx q[0];
rz(-1.1044015) q[0];
sx q[0];
rz(-1.3155235) q[0];
rz(0.77519351) q[1];
sx q[1];
rz(-0.23861353) q[1];
sx q[1];
rz(-1.7528037) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9520541) q[0];
sx q[0];
rz(-1.4217958) q[0];
sx q[0];
rz(2.8556264) q[0];
rz(-pi) q[1];
x q[1];
rz(0.49225537) q[2];
sx q[2];
rz(-0.45212072) q[2];
sx q[2];
rz(-1.4398354) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.57143104) q[1];
sx q[1];
rz(-1.5699016) q[1];
sx q[1];
rz(-1.5714297) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7241012) q[3];
sx q[3];
rz(-1.1840828) q[3];
sx q[3];
rz(1.2717122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3751601) q[2];
sx q[2];
rz(-3.1258686) q[2];
sx q[2];
rz(1.6895705) q[2];
rz(-2.8826513) q[3];
sx q[3];
rz(-0.18766923) q[3];
sx q[3];
rz(-1.1567206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6257085) q[0];
sx q[0];
rz(-2.4242171) q[0];
sx q[0];
rz(-1.7051359) q[0];
rz(1.4868078) q[1];
sx q[1];
rz(-0.27324067) q[1];
sx q[1];
rz(-2.9437093) q[1];
rz(0.033587348) q[2];
sx q[2];
rz(-0.20034364) q[2];
sx q[2];
rz(0.26058254) q[2];
rz(-1.6428357) q[3];
sx q[3];
rz(-2.7938953) q[3];
sx q[3];
rz(0.10216879) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
