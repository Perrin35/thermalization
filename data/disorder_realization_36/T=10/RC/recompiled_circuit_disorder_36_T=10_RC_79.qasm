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
rz(-2.8556813) q[0];
sx q[0];
rz(-0.51529348) q[0];
rz(-1.7973068) q[1];
sx q[1];
rz(-0.15434115) q[1];
sx q[1];
rz(-0.57758346) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5776423) q[0];
sx q[0];
rz(-0.9194153) q[0];
sx q[0];
rz(-1.7786068) q[0];
x q[1];
rz(2.3876786) q[2];
sx q[2];
rz(-0.68281096) q[2];
sx q[2];
rz(2.4726601) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5756302) q[1];
sx q[1];
rz(-1.9735676) q[1];
sx q[1];
rz(2.8494542) q[1];
rz(-pi) q[2];
rz(-0.9382117) q[3];
sx q[3];
rz(-2.1198366) q[3];
sx q[3];
rz(-3.0905746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.43705964) q[2];
sx q[2];
rz(-1.5455064) q[2];
sx q[2];
rz(0.68721592) q[2];
rz(-1.0152738) q[3];
sx q[3];
rz(-1.3736558) q[3];
sx q[3];
rz(-0.12250531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9706443) q[0];
sx q[0];
rz(-2.0630554) q[0];
sx q[0];
rz(-1.8815536) q[0];
rz(2.1353703) q[1];
sx q[1];
rz(-2.1496014) q[1];
sx q[1];
rz(0.84567436) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0613522) q[0];
sx q[0];
rz(-2.4363359) q[0];
sx q[0];
rz(-2.4387226) q[0];
rz(-pi) q[1];
rz(-0.77180736) q[2];
sx q[2];
rz(-0.1856005) q[2];
sx q[2];
rz(2.0303357) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0452803) q[1];
sx q[1];
rz(-1.7573866) q[1];
sx q[1];
rz(-0.50218302) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2250697) q[3];
sx q[3];
rz(-0.73361165) q[3];
sx q[3];
rz(-0.77378002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3530897) q[2];
sx q[2];
rz(-0.22558364) q[2];
sx q[2];
rz(2.6611924) q[2];
rz(-1.7885615) q[3];
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
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1903494) q[0];
sx q[0];
rz(-0.23878637) q[0];
sx q[0];
rz(-0.7730661) q[0];
rz(0.13126016) q[1];
sx q[1];
rz(-1.8570329) q[1];
sx q[1];
rz(1.0864331) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7032996) q[0];
sx q[0];
rz(-0.46143954) q[0];
sx q[0];
rz(-1.567054) q[0];
rz(-pi) q[1];
rz(-2.2181182) q[2];
sx q[2];
rz(-1.662865) q[2];
sx q[2];
rz(-2.0331241) q[2];
rz(-pi/2) q[3];
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
rz(-0.28880854) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9757189) q[3];
sx q[3];
rz(-1.4110663) q[3];
sx q[3];
rz(-0.83042849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1290258) q[2];
sx q[2];
rz(-1.4386703) q[2];
sx q[2];
rz(-1.3712937) q[2];
rz(2.7584372) q[3];
sx q[3];
rz(-1.8846735) q[3];
sx q[3];
rz(-2.3390521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4521769) q[0];
sx q[0];
rz(-1.8912264) q[0];
sx q[0];
rz(-0.048359811) q[0];
rz(-2.9776749) q[1];
sx q[1];
rz(-0.36968958) q[1];
sx q[1];
rz(-1.6960467) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.018054124) q[0];
sx q[0];
rz(-0.92045438) q[0];
sx q[0];
rz(-1.3440078) q[0];
x q[1];
rz(-2.9252058) q[2];
sx q[2];
rz(-1.4828223) q[2];
sx q[2];
rz(-0.83425922) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4353664) q[1];
sx q[1];
rz(-1.6352143) q[1];
sx q[1];
rz(0.21965841) q[1];
rz(-0.28820949) q[3];
sx q[3];
rz(-1.0681579) q[3];
sx q[3];
rz(0.99311815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.066102862) q[2];
sx q[2];
rz(-1.3572071) q[2];
sx q[2];
rz(2.1172822) q[2];
rz(1.5284437) q[3];
sx q[3];
rz(-1.6201092) q[3];
sx q[3];
rz(0.23322341) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7016474) q[0];
sx q[0];
rz(-0.82413903) q[0];
sx q[0];
rz(1.8540927) q[0];
rz(-0.31907407) q[1];
sx q[1];
rz(-1.5998452) q[1];
sx q[1];
rz(0.85420001) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.094164205) q[0];
sx q[0];
rz(-1.578997) q[0];
sx q[0];
rz(-2.2846983) q[0];
rz(-pi) q[1];
rz(-1.008026) q[2];
sx q[2];
rz(-0.78283435) q[2];
sx q[2];
rz(-0.40751878) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2367868) q[1];
sx q[1];
rz(-0.85610897) q[1];
sx q[1];
rz(-0.76101117) q[1];
rz(-pi) q[2];
rz(2.5913127) q[3];
sx q[3];
rz(-2.0645421) q[3];
sx q[3];
rz(-0.20875904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.95191082) q[2];
sx q[2];
rz(-0.59331912) q[2];
sx q[2];
rz(2.5642776) q[2];
rz(-2.632085) q[3];
sx q[3];
rz(-0.43764344) q[3];
sx q[3];
rz(0.90665162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.0034870738) q[0];
sx q[0];
rz(-1.071799) q[0];
sx q[0];
rz(0.072120897) q[0];
rz(1.1068608) q[1];
sx q[1];
rz(-2.6289584) q[1];
sx q[1];
rz(-3.0153826) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5835411) q[0];
sx q[0];
rz(-1.3582894) q[0];
sx q[0];
rz(1.538518) q[0];
rz(-2.3665479) q[2];
sx q[2];
rz(-1.8473052) q[2];
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
rz(1.0969442) q[1];
sx q[1];
rz(-1.7587887) q[1];
sx q[1];
rz(0.85) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8993127) q[3];
sx q[3];
rz(-1.7528755) q[3];
sx q[3];
rz(-1.8329221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.90298992) q[2];
sx q[2];
rz(-2.949252) q[2];
sx q[2];
rz(-2.3664756) q[2];
rz(-2.3136247) q[3];
sx q[3];
rz(-2.8505846) q[3];
sx q[3];
rz(-1.1221788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59259748) q[0];
sx q[0];
rz(-2.7153375) q[0];
sx q[0];
rz(-0.098408498) q[0];
rz(1.9495643) q[1];
sx q[1];
rz(-1.8076618) q[1];
sx q[1];
rz(-0.55955204) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0347621) q[0];
sx q[0];
rz(-0.97495279) q[0];
sx q[0];
rz(1.725127) q[0];
rz(-pi) q[1];
rz(2.8801708) q[2];
sx q[2];
rz(-1.520642) q[2];
sx q[2];
rz(2.9068771) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.050484867) q[1];
sx q[1];
rz(-2.9851966) q[1];
sx q[1];
rz(-0.97538235) q[1];
rz(1.3748752) q[3];
sx q[3];
rz(-2.19176) q[3];
sx q[3];
rz(1.2697112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6825535) q[2];
sx q[2];
rz(-1.295853) q[2];
sx q[2];
rz(2.7977978) q[2];
rz(-0.5665468) q[3];
sx q[3];
rz(-2.6930801) q[3];
sx q[3];
rz(-0.47376537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.375181) q[0];
sx q[0];
rz(-1.8090929) q[0];
sx q[0];
rz(-0.73076105) q[0];
rz(0.14239755) q[1];
sx q[1];
rz(-1.8715033) q[1];
sx q[1];
rz(0.87160814) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14411892) q[0];
sx q[0];
rz(-1.5796356) q[0];
sx q[0];
rz(3.0661422) q[0];
rz(-pi) q[1];
rz(-1.4191188) q[2];
sx q[2];
rz(-1.1693839) q[2];
sx q[2];
rz(-0.36827189) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.65054446) q[1];
sx q[1];
rz(-2.3666413) q[1];
sx q[1];
rz(2.6130555) q[1];
x q[2];
rz(1.7730764) q[3];
sx q[3];
rz(-2.5296122) q[3];
sx q[3];
rz(2.4806541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.72620755) q[2];
sx q[2];
rz(-2.0724847) q[2];
sx q[2];
rz(-1.139572) q[2];
rz(-1.4987882) q[3];
sx q[3];
rz(-2.747624) q[3];
sx q[3];
rz(0.9128226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.9386439) q[0];
sx q[0];
rz(-1.6904172) q[0];
sx q[0];
rz(1.2217481) q[0];
rz(2.9755759) q[1];
sx q[1];
rz(-1.32042) q[1];
sx q[1];
rz(1.6171914) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0696408) q[0];
sx q[0];
rz(-1.6531634) q[0];
sx q[0];
rz(0.029043341) q[0];
x q[1];
rz(1.379307) q[2];
sx q[2];
rz(-1.101149) q[2];
sx q[2];
rz(0.42275235) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2251687) q[1];
sx q[1];
rz(-1.4125707) q[1];
sx q[1];
rz(-1.9041063) q[1];
x q[2];
rz(0.63906007) q[3];
sx q[3];
rz(-1.8250873) q[3];
sx q[3];
rz(1.1921079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.021585492) q[2];
sx q[2];
rz(-1.465613) q[2];
sx q[2];
rz(-0.35153708) q[2];
rz(1.0567788) q[3];
sx q[3];
rz(-0.52962279) q[3];
sx q[3];
rz(-0.74469152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4979424) q[0];
sx q[0];
rz(-0.90181667) q[0];
sx q[0];
rz(-1.3051916) q[0];
rz(-0.38048831) q[1];
sx q[1];
rz(-1.0419798) q[1];
sx q[1];
rz(-2.8881853) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5288552) q[0];
sx q[0];
rz(-0.40898541) q[0];
sx q[0];
rz(-1.3091875) q[0];
x q[1];
rz(-2.017574) q[2];
sx q[2];
rz(-2.7687216) q[2];
sx q[2];
rz(-0.98137059) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7518172) q[1];
sx q[1];
rz(-2.7174065) q[1];
sx q[1];
rz(0.61702375) q[1];
x q[2];
rz(-2.5640423) q[3];
sx q[3];
rz(-0.69637075) q[3];
sx q[3];
rz(1.3814572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5499251) q[2];
sx q[2];
rz(-2.2558236) q[2];
sx q[2];
rz(-0.5029451) q[2];
rz(2.2425966) q[3];
sx q[3];
rz(-1.2939724) q[3];
sx q[3];
rz(1.9780654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.74264991) q[2];
sx q[2];
rz(-0.37336083) q[2];
sx q[2];
rz(0.20867418) q[2];
rz(-2.4651299) q[3];
sx q[3];
rz(-0.87920311) q[3];
sx q[3];
rz(1.9999947) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
