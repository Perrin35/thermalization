OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.71896267) q[0];
sx q[0];
rz(-0.29932061) q[0];
sx q[0];
rz(-2.646995) q[0];
rz(-1.9994796) q[1];
sx q[1];
rz(-2.1358868) q[1];
sx q[1];
rz(-1.1297273) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8324234) q[0];
sx q[0];
rz(-2.4465843) q[0];
sx q[0];
rz(2.8929936) q[0];
rz(-pi) q[1];
rz(2.0689993) q[2];
sx q[2];
rz(-1.9747726) q[2];
sx q[2];
rz(-2.7369268) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2430302) q[1];
sx q[1];
rz(-0.77706017) q[1];
sx q[1];
rz(1.9352566) q[1];
x q[2];
rz(2.7143728) q[3];
sx q[3];
rz(-2.1764206) q[3];
sx q[3];
rz(2.1744414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8226681) q[2];
sx q[2];
rz(-1.319004) q[2];
sx q[2];
rz(0.85282105) q[2];
rz(-1.3301814) q[3];
sx q[3];
rz(-0.69555247) q[3];
sx q[3];
rz(0.046796355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80832076) q[0];
sx q[0];
rz(-3.0905753) q[0];
sx q[0];
rz(1.464123) q[0];
rz(1.6630215) q[1];
sx q[1];
rz(-1.9824948) q[1];
sx q[1];
rz(1.9333855) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0639275) q[0];
sx q[0];
rz(-2.0443235) q[0];
sx q[0];
rz(0.60654503) q[0];
x q[1];
rz(0.21734997) q[2];
sx q[2];
rz(-2.0361379) q[2];
sx q[2];
rz(0.67721043) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7866655) q[1];
sx q[1];
rz(-1.8372722) q[1];
sx q[1];
rz(3.0158426) q[1];
rz(-pi) q[2];
rz(-2.6868846) q[3];
sx q[3];
rz(-1.4147926) q[3];
sx q[3];
rz(1.7059513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6055484) q[2];
sx q[2];
rz(-0.41993419) q[2];
sx q[2];
rz(1.4844683) q[2];
rz(-3.1260955) q[3];
sx q[3];
rz(-1.9280547) q[3];
sx q[3];
rz(-1.9626455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.991796) q[0];
sx q[0];
rz(-1.8114256) q[0];
sx q[0];
rz(2.8699744) q[0];
rz(-2.244921) q[1];
sx q[1];
rz(-0.50183693) q[1];
sx q[1];
rz(1.2976049) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54416849) q[0];
sx q[0];
rz(-0.46314683) q[0];
sx q[0];
rz(1.5482272) q[0];
x q[1];
rz(-2.6372725) q[2];
sx q[2];
rz(-1.9078622) q[2];
sx q[2];
rz(2.6884109) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5601968) q[1];
sx q[1];
rz(-0.27276688) q[1];
sx q[1];
rz(-2.0062937) q[1];
rz(-0.0053169189) q[3];
sx q[3];
rz(-2.379619) q[3];
sx q[3];
rz(0.055082037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4012332) q[2];
sx q[2];
rz(-0.77784246) q[2];
sx q[2];
rz(3.0878301) q[2];
rz(-1.2939804) q[3];
sx q[3];
rz(-0.79837489) q[3];
sx q[3];
rz(0.23538858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(2.9409222) q[0];
sx q[0];
rz(-0.5680474) q[0];
sx q[0];
rz(-1.4642375) q[0];
rz(-1.1306521) q[1];
sx q[1];
rz(-0.73431763) q[1];
sx q[1];
rz(-3.0272223) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8371724) q[0];
sx q[0];
rz(-2.3995598) q[0];
sx q[0];
rz(-2.4993308) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61538265) q[2];
sx q[2];
rz(-0.75591959) q[2];
sx q[2];
rz(2.4494954) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3614157) q[1];
sx q[1];
rz(-0.89539278) q[1];
sx q[1];
rz(3.1076868) q[1];
rz(0.19802494) q[3];
sx q[3];
rz(-1.5643483) q[3];
sx q[3];
rz(-0.40606582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4734681) q[2];
sx q[2];
rz(-1.5359842) q[2];
sx q[2];
rz(0.075695666) q[2];
rz(-2.7068052) q[3];
sx q[3];
rz(-1.1554759) q[3];
sx q[3];
rz(-0.079631478) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39744034) q[0];
sx q[0];
rz(-1.8295153) q[0];
sx q[0];
rz(0.42006668) q[0];
rz(-2.9282667) q[1];
sx q[1];
rz(-1.7330287) q[1];
sx q[1];
rz(-1.2423645) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38753375) q[0];
sx q[0];
rz(-2.2926169) q[0];
sx q[0];
rz(-2.5639064) q[0];
rz(-pi) q[1];
rz(2.7105646) q[2];
sx q[2];
rz(-1.8371474) q[2];
sx q[2];
rz(-0.14542994) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2957569) q[1];
sx q[1];
rz(-1.6713665) q[1];
sx q[1];
rz(-3.1261338) q[1];
rz(-pi) q[2];
x q[2];
rz(2.931704) q[3];
sx q[3];
rz(-0.22612962) q[3];
sx q[3];
rz(-1.1864288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8215948) q[2];
sx q[2];
rz(-2.2264806) q[2];
sx q[2];
rz(0.37459174) q[2];
rz(-1.0125259) q[3];
sx q[3];
rz(-2.7495224) q[3];
sx q[3];
rz(0.64468002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0300765) q[0];
sx q[0];
rz(-0.498963) q[0];
sx q[0];
rz(2.7110355) q[0];
rz(0.14398362) q[1];
sx q[1];
rz(-1.0824243) q[1];
sx q[1];
rz(0.63124257) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86595097) q[0];
sx q[0];
rz(-1.5168204) q[0];
sx q[0];
rz(-0.60906054) q[0];
x q[1];
rz(-0.14603931) q[2];
sx q[2];
rz(-1.1542873) q[2];
sx q[2];
rz(0.064695875) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.65578016) q[1];
sx q[1];
rz(-0.98193278) q[1];
sx q[1];
rz(0.9649802) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.33511929) q[3];
sx q[3];
rz(-2.0010173) q[3];
sx q[3];
rz(2.8443401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1708019) q[2];
sx q[2];
rz(-1.0589212) q[2];
sx q[2];
rz(1.4405174) q[2];
rz(-1.3614281) q[3];
sx q[3];
rz(-2.0378588) q[3];
sx q[3];
rz(0.48719278) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.454527) q[0];
sx q[0];
rz(-2.8519958) q[0];
sx q[0];
rz(2.2173296) q[0];
rz(-0.29280064) q[1];
sx q[1];
rz(-1.9127138) q[1];
sx q[1];
rz(2.3407095) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7445114) q[0];
sx q[0];
rz(-1.9430706) q[0];
sx q[0];
rz(0.27639322) q[0];
rz(-2.2325781) q[2];
sx q[2];
rz(-1.2940727) q[2];
sx q[2];
rz(0.02034517) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.427721) q[1];
sx q[1];
rz(-1.5414951) q[1];
sx q[1];
rz(-3.0364365) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.8893251) q[3];
sx q[3];
rz(-3.1066537) q[3];
sx q[3];
rz(-0.068471758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.822829) q[2];
sx q[2];
rz(-1.2206565) q[2];
sx q[2];
rz(1.2804383) q[2];
rz(0.66796962) q[3];
sx q[3];
rz(-0.48798713) q[3];
sx q[3];
rz(2.7282696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25705826) q[0];
sx q[0];
rz(-1.9030544) q[0];
sx q[0];
rz(1.1827693) q[0];
rz(0.23712748) q[1];
sx q[1];
rz(-3.0393937) q[1];
sx q[1];
rz(-3.0822486) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5049669) q[0];
sx q[0];
rz(-1.6904181) q[0];
sx q[0];
rz(-1.8455532) q[0];
x q[1];
rz(0.58299033) q[2];
sx q[2];
rz(-0.81524476) q[2];
sx q[2];
rz(-1.5423403) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0036111) q[1];
sx q[1];
rz(-1.6754158) q[1];
sx q[1];
rz(-1.8358747) q[1];
rz(-1.7982676) q[3];
sx q[3];
rz(-2.8099647) q[3];
sx q[3];
rz(-2.22081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.55988971) q[2];
sx q[2];
rz(-0.54055944) q[2];
sx q[2];
rz(-1.3524559) q[2];
rz(-1.7743568) q[3];
sx q[3];
rz(-1.7188027) q[3];
sx q[3];
rz(0.8555612) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2818114) q[0];
sx q[0];
rz(-0.78524041) q[0];
sx q[0];
rz(-0.45968858) q[0];
rz(3.1135318) q[1];
sx q[1];
rz(-1.9919845) q[1];
sx q[1];
rz(1.185816) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7677886) q[0];
sx q[0];
rz(-1.8361183) q[0];
sx q[0];
rz(-1.6173043) q[0];
rz(-pi) q[1];
rz(1.6721647) q[2];
sx q[2];
rz(-1.4918054) q[2];
sx q[2];
rz(-1.422342) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.016077135) q[1];
sx q[1];
rz(-1.0927769) q[1];
sx q[1];
rz(2.1090871) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5236058) q[3];
sx q[3];
rz(-1.3387965) q[3];
sx q[3];
rz(2.8232676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6431553) q[2];
sx q[2];
rz(-2.2648621) q[2];
sx q[2];
rz(0.15677491) q[2];
rz(-2.5458941) q[3];
sx q[3];
rz(-0.34606338) q[3];
sx q[3];
rz(-3.116385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51782411) q[0];
sx q[0];
rz(-2.1110004) q[0];
sx q[0];
rz(-0.18381707) q[0];
rz(0.078016438) q[1];
sx q[1];
rz(-1.0573496) q[1];
sx q[1];
rz(-0.26430166) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83885709) q[0];
sx q[0];
rz(-2.502901) q[0];
sx q[0];
rz(-0.46052082) q[0];
rz(2.5223612) q[2];
sx q[2];
rz(-0.78699099) q[2];
sx q[2];
rz(0.84112043) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2956315) q[1];
sx q[1];
rz(-0.23989883) q[1];
sx q[1];
rz(-2.4227581) q[1];
rz(2.6464858) q[3];
sx q[3];
rz(-1.4864085) q[3];
sx q[3];
rz(1.5051248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2901624) q[2];
sx q[2];
rz(-2.1992079) q[2];
sx q[2];
rz(0.22872049) q[2];
rz(-1.4043572) q[3];
sx q[3];
rz(-0.93969932) q[3];
sx q[3];
rz(3.0586045) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9901154) q[0];
sx q[0];
rz(-0.84828068) q[0];
sx q[0];
rz(-1.620851) q[0];
rz(2.4304541) q[1];
sx q[1];
rz(-1.3093206) q[1];
sx q[1];
rz(0.73285229) q[1];
rz(-2.0696832) q[2];
sx q[2];
rz(-1.6098534) q[2];
sx q[2];
rz(0.44865566) q[2];
rz(-0.55599273) q[3];
sx q[3];
rz(-1.0992194) q[3];
sx q[3];
rz(1.9627375) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
