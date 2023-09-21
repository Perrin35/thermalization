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
rz(0.51529348) q[0];
rz(-1.7973068) q[1];
sx q[1];
rz(-0.15434115) q[1];
sx q[1];
rz(-0.57758346) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.261895) q[0];
sx q[0];
rz(-1.4059773) q[0];
sx q[0];
rz(-0.66189712) q[0];
rz(2.0787813) q[2];
sx q[2];
rz(-1.0928109) q[2];
sx q[2];
rz(0.21131549) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.88749332) q[1];
sx q[1];
rz(-1.3026397) q[1];
sx q[1];
rz(-1.1521794) q[1];
rz(-pi) q[2];
rz(-2.203381) q[3];
sx q[3];
rz(-2.1198366) q[3];
sx q[3];
rz(3.0905746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43705964) q[2];
sx q[2];
rz(-1.5960863) q[2];
sx q[2];
rz(-2.4543767) q[2];
rz(-2.1263188) q[3];
sx q[3];
rz(-1.7679368) q[3];
sx q[3];
rz(3.0190873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17094831) q[0];
sx q[0];
rz(-2.0630554) q[0];
sx q[0];
rz(-1.8815536) q[0];
rz(-1.0062224) q[1];
sx q[1];
rz(-2.1496014) q[1];
sx q[1];
rz(0.84567436) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0613522) q[0];
sx q[0];
rz(-2.4363359) q[0];
sx q[0];
rz(0.7028701) q[0];
x q[1];
rz(-2.3697853) q[2];
sx q[2];
rz(-0.1856005) q[2];
sx q[2];
rz(1.1112569) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.80026514) q[1];
sx q[1];
rz(-2.6086573) q[1];
sx q[1];
rz(-0.37377263) q[1];
rz(-pi) q[2];
x q[2];
rz(0.50176974) q[3];
sx q[3];
rz(-2.1309149) q[3];
sx q[3];
rz(3.113941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3530897) q[2];
sx q[2];
rz(-2.916009) q[2];
sx q[2];
rz(-2.6611924) q[2];
rz(1.7885615) q[3];
sx q[3];
rz(-1.055911) q[3];
sx q[3];
rz(-1.1876748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1903494) q[0];
sx q[0];
rz(-2.9028063) q[0];
sx q[0];
rz(-2.3685266) q[0];
rz(3.0103325) q[1];
sx q[1];
rz(-1.8570329) q[1];
sx q[1];
rz(-1.0864331) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4341136) q[0];
sx q[0];
rz(-2.0322324) q[0];
sx q[0];
rz(0.001860851) q[0];
x q[1];
rz(-2.2181182) q[2];
sx q[2];
rz(-1.4787276) q[2];
sx q[2];
rz(2.0331241) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4194581) q[1];
sx q[1];
rz(-2.8527383) q[1];
sx q[1];
rz(-3.1232749) q[1];
rz(-1.1658737) q[3];
sx q[3];
rz(-1.4110663) q[3];
sx q[3];
rz(-2.3111642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1290258) q[2];
sx q[2];
rz(-1.7029224) q[2];
sx q[2];
rz(1.770299) q[2];
rz(0.38315547) q[3];
sx q[3];
rz(-1.8846735) q[3];
sx q[3];
rz(2.3390521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.6894158) q[0];
sx q[0];
rz(-1.2503662) q[0];
sx q[0];
rz(3.0932328) q[0];
rz(0.16391779) q[1];
sx q[1];
rz(-0.36968958) q[1];
sx q[1];
rz(1.4455459) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.018054124) q[0];
sx q[0];
rz(-0.92045438) q[0];
sx q[0];
rz(-1.3440078) q[0];
x q[1];
rz(-1.4807329) q[2];
sx q[2];
rz(-1.3552595) q[2];
sx q[2];
rz(2.4243674) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.85019894) q[1];
sx q[1];
rz(-1.7899917) q[1];
sx q[1];
rz(-1.5047969) q[1];
rz(-2.8533832) q[3];
sx q[3];
rz(-2.0734348) q[3];
sx q[3];
rz(-2.1484745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.066102862) q[2];
sx q[2];
rz(-1.7843856) q[2];
sx q[2];
rz(2.1172822) q[2];
rz(-1.6131489) q[3];
sx q[3];
rz(-1.5214835) q[3];
sx q[3];
rz(-0.23322341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7016474) q[0];
sx q[0];
rz(-0.82413903) q[0];
sx q[0];
rz(1.2874999) q[0];
rz(-0.31907407) q[1];
sx q[1];
rz(-1.5998452) q[1];
sx q[1];
rz(0.85420001) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0474284) q[0];
sx q[0];
rz(-1.578997) q[0];
sx q[0];
rz(-0.8568944) q[0];
rz(-0.87128432) q[2];
sx q[2];
rz(-1.185002) q[2];
sx q[2];
rz(0.7427578) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.10798771) q[1];
sx q[1];
rz(-2.1186947) q[1];
sx q[1];
rz(-2.4461436) q[1];
x q[2];
rz(-0.55027996) q[3];
sx q[3];
rz(-2.0645421) q[3];
sx q[3];
rz(-0.20875904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1896818) q[2];
sx q[2];
rz(-2.5482735) q[2];
sx q[2];
rz(2.5642776) q[2];
rz(0.50950766) q[3];
sx q[3];
rz(-0.43764344) q[3];
sx q[3];
rz(-2.234941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0034870738) q[0];
sx q[0];
rz(-2.0697937) q[0];
sx q[0];
rz(0.072120897) q[0];
rz(1.1068608) q[1];
sx q[1];
rz(-2.6289584) q[1];
sx q[1];
rz(0.12621005) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.135658) q[0];
sx q[0];
rz(-1.5392443) q[0];
sx q[0];
rz(0.21261442) q[0];
rz(-pi) q[1];
rz(-1.1926786) q[2];
sx q[2];
rz(-2.3092804) q[2];
sx q[2];
rz(-2.9044915) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8773552) q[1];
sx q[1];
rz(-2.4009631) q[1];
sx q[1];
rz(1.8514368) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.19212171) q[3];
sx q[3];
rz(-1.2479094) q[3];
sx q[3];
rz(0.20048143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.90298992) q[2];
sx q[2];
rz(-2.949252) q[2];
sx q[2];
rz(-2.3664756) q[2];
rz(-0.827968) q[3];
sx q[3];
rz(-2.8505846) q[3];
sx q[3];
rz(1.1221788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5489952) q[0];
sx q[0];
rz(-0.42625517) q[0];
sx q[0];
rz(-0.098408498) q[0];
rz(-1.9495643) q[1];
sx q[1];
rz(-1.8076618) q[1];
sx q[1];
rz(0.55955204) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37724272) q[0];
sx q[0];
rz(-0.61315216) q[0];
sx q[0];
rz(-2.9186547) q[0];
x q[1];
rz(-0.19182972) q[2];
sx q[2];
rz(-0.26608135) q[2];
sx q[2];
rz(-1.9907469) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2110062) q[1];
sx q[1];
rz(-1.483327) q[1];
sx q[1];
rz(1.4409815) q[1];
x q[2];
rz(-0.26569326) q[3];
sx q[3];
rz(-2.4943647) q[3];
sx q[3];
rz(2.2006187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.45903912) q[2];
sx q[2];
rz(-1.8457396) q[2];
sx q[2];
rz(-2.7977978) q[2];
rz(-2.5750459) q[3];
sx q[3];
rz(-0.44851258) q[3];
sx q[3];
rz(2.6678273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7664117) q[0];
sx q[0];
rz(-1.8090929) q[0];
sx q[0];
rz(-2.4108316) q[0];
rz(0.14239755) q[1];
sx q[1];
rz(-1.8715033) q[1];
sx q[1];
rz(-2.2699845) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5985142) q[0];
sx q[0];
rz(-3.0656272) q[0];
sx q[0];
rz(0.11673467) q[0];
rz(-pi) q[1];
rz(0.34198728) q[2];
sx q[2];
rz(-0.42765289) q[2];
sx q[2];
rz(0.74117408) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3155047) q[1];
sx q[1];
rz(-1.9313889) q[1];
sx q[1];
rz(0.70198595) q[1];
rz(-0.96846795) q[3];
sx q[3];
rz(-1.6864711) q[3];
sx q[3];
rz(1.0761716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4153851) q[2];
sx q[2];
rz(-2.0724847) q[2];
sx q[2];
rz(2.0020206) q[2];
rz(-1.4987882) q[3];
sx q[3];
rz(-2.747624) q[3];
sx q[3];
rz(-2.22877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(0.20294872) q[0];
sx q[0];
rz(-1.6904172) q[0];
sx q[0];
rz(1.9198445) q[0];
rz(-0.16601673) q[1];
sx q[1];
rz(-1.32042) q[1];
sx q[1];
rz(1.6171914) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6451384) q[0];
sx q[0];
rz(-1.5997412) q[0];
sx q[0];
rz(-1.4883947) q[0];
x q[1];
rz(-2.7828214) q[2];
sx q[2];
rz(-0.504474) q[2];
sx q[2];
rz(-2.3141253) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.4001273) q[1];
sx q[1];
rz(-1.2418081) q[1];
sx q[1];
rz(0.16727438) q[1];
x q[2];
rz(-0.63906007) q[3];
sx q[3];
rz(-1.3165054) q[3];
sx q[3];
rz(1.1921079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1200072) q[2];
sx q[2];
rz(-1.6759796) q[2];
sx q[2];
rz(0.35153708) q[2];
rz(1.0567788) q[3];
sx q[3];
rz(-2.6119699) q[3];
sx q[3];
rz(-2.3969011) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4979424) q[0];
sx q[0];
rz(-0.90181667) q[0];
sx q[0];
rz(1.836401) q[0];
rz(0.38048831) q[1];
sx q[1];
rz(-1.0419798) q[1];
sx q[1];
rz(2.8881853) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8587592) q[0];
sx q[0];
rz(-1.6738322) q[0];
sx q[0];
rz(-1.1742924) q[0];
rz(-pi) q[1];
rz(-0.16742736) q[2];
sx q[2];
rz(-1.2360459) q[2];
sx q[2];
rz(1.6850922) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3897755) q[1];
sx q[1];
rz(-2.7174065) q[1];
sx q[1];
rz(-2.5245689) q[1];
rz(-pi) q[2];
x q[2];
rz(0.57755034) q[3];
sx q[3];
rz(-2.4452219) q[3];
sx q[3];
rz(-1.3814572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5499251) q[2];
sx q[2];
rz(-2.2558236) q[2];
sx q[2];
rz(-2.6386476) q[2];
rz(0.89899603) q[3];
sx q[3];
rz(-1.8476202) q[3];
sx q[3];
rz(1.9780654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9713365) q[0];
sx q[0];
rz(-1.5383056) q[0];
sx q[0];
rz(-2.8785895) q[0];
rz(-0.7111711) q[1];
sx q[1];
rz(-1.0881337) q[1];
sx q[1];
rz(1.7137391) q[1];
rz(-2.3989427) q[2];
sx q[2];
rz(-2.7682318) q[2];
sx q[2];
rz(-2.9329185) q[2];
rz(-0.75541227) q[3];
sx q[3];
rz(-2.073954) q[3];
sx q[3];
rz(3.0975773) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
