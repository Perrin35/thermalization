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
rz(3.427504) q[0];
sx q[0];
rz(8.9094845) q[0];
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
rz(0.89864697) q[0];
sx q[0];
rz(-2.4624914) q[0];
sx q[0];
rz(-0.26429096) q[0];
rz(-pi) q[1];
rz(-1.0628113) q[2];
sx q[2];
rz(-2.0487818) q[2];
sx q[2];
rz(-0.21131549) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.88749332) q[1];
sx q[1];
rz(-1.8389529) q[1];
sx q[1];
rz(1.9894132) q[1];
rz(0.9382117) q[3];
sx q[3];
rz(-2.1198366) q[3];
sx q[3];
rz(-0.051018056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.43705964) q[2];
sx q[2];
rz(-1.5960863) q[2];
sx q[2];
rz(-0.68721592) q[2];
rz(1.0152738) q[3];
sx q[3];
rz(-1.3736558) q[3];
sx q[3];
rz(-3.0190873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17094831) q[0];
sx q[0];
rz(-1.0785372) q[0];
sx q[0];
rz(-1.8815536) q[0];
rz(1.0062224) q[1];
sx q[1];
rz(-2.1496014) q[1];
sx q[1];
rz(-0.84567436) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0592247) q[0];
sx q[0];
rz(-1.1384283) q[0];
sx q[0];
rz(-2.5655377) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3697853) q[2];
sx q[2];
rz(-2.9559921) q[2];
sx q[2];
rz(-1.1112569) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.80026514) q[1];
sx q[1];
rz(-0.53293537) q[1];
sx q[1];
rz(0.37377263) q[1];
rz(-pi) q[2];
rz(2.6398229) q[3];
sx q[3];
rz(-2.1309149) q[3];
sx q[3];
rz(0.027651699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3530897) q[2];
sx q[2];
rz(-2.916009) q[2];
sx q[2];
rz(2.6611924) q[2];
rz(1.3530312) q[3];
sx q[3];
rz(-2.0856817) q[3];
sx q[3];
rz(1.9539179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.95124328) q[0];
sx q[0];
rz(-2.9028063) q[0];
sx q[0];
rz(0.7730661) q[0];
rz(0.13126016) q[1];
sx q[1];
rz(-1.2845598) q[1];
sx q[1];
rz(2.0551596) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4341136) q[0];
sx q[0];
rz(-1.1093603) q[0];
sx q[0];
rz(3.1397318) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.418872) q[2];
sx q[2];
rz(-2.488689) q[2];
sx q[2];
rz(-2.8002847) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.8311027) q[1];
sx q[1];
rz(-1.576014) q[1];
sx q[1];
rz(-2.8527841) q[1];
x q[2];
rz(-0.17351405) q[3];
sx q[3];
rz(-1.9702692) q[3];
sx q[3];
rz(-2.3331593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1290258) q[2];
sx q[2];
rz(-1.7029224) q[2];
sx q[2];
rz(1.770299) q[2];
rz(-0.38315547) q[3];
sx q[3];
rz(-1.2569191) q[3];
sx q[3];
rz(2.3390521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.4521769) q[0];
sx q[0];
rz(-1.2503662) q[0];
sx q[0];
rz(0.048359811) q[0];
rz(-2.9776749) q[1];
sx q[1];
rz(-0.36968958) q[1];
sx q[1];
rz(-1.6960467) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6915582) q[0];
sx q[0];
rz(-1.7507179) q[0];
sx q[0];
rz(2.4787089) q[0];
rz(-0.21638685) q[2];
sx q[2];
rz(-1.4828223) q[2];
sx q[2];
rz(0.83425922) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.85019894) q[1];
sx q[1];
rz(-1.7899917) q[1];
sx q[1];
rz(-1.6367957) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0936071) q[3];
sx q[3];
rz(-2.5683937) q[3];
sx q[3];
rz(-1.5968061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.066102862) q[2];
sx q[2];
rz(-1.3572071) q[2];
sx q[2];
rz(2.1172822) q[2];
rz(-1.6131489) q[3];
sx q[3];
rz(-1.6201092) q[3];
sx q[3];
rz(0.23322341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.7016474) q[0];
sx q[0];
rz(-0.82413903) q[0];
sx q[0];
rz(1.2874999) q[0];
rz(0.31907407) q[1];
sx q[1];
rz(-1.5998452) q[1];
sx q[1];
rz(-0.85420001) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4860977) q[0];
sx q[0];
rz(-2.4276519) q[0];
sx q[0];
rz(1.5582725) q[0];
rz(-pi) q[1];
x q[1];
rz(0.87128432) q[2];
sx q[2];
rz(-1.185002) q[2];
sx q[2];
rz(-0.7427578) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0336049) q[1];
sx q[1];
rz(-2.1186947) q[1];
sx q[1];
rz(-0.69544905) q[1];
x q[2];
rz(2.1339995) q[3];
sx q[3];
rz(-2.0493205) q[3];
sx q[3];
rz(1.0790881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.95191082) q[2];
sx q[2];
rz(-0.59331912) q[2];
sx q[2];
rz(-0.577315) q[2];
rz(0.50950766) q[3];
sx q[3];
rz(-0.43764344) q[3];
sx q[3];
rz(-2.234941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0034870738) q[0];
sx q[0];
rz(-1.071799) q[0];
sx q[0];
rz(3.0694718) q[0];
rz(1.1068608) q[1];
sx q[1];
rz(-2.6289584) q[1];
sx q[1];
rz(-3.0153826) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5835411) q[0];
sx q[0];
rz(-1.7833033) q[0];
sx q[0];
rz(-1.6030747) q[0];
rz(-pi) q[1];
x q[1];
rz(0.38527617) q[2];
sx q[2];
rz(-2.3284973) q[2];
sx q[2];
rz(2.3713881) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0969442) q[1];
sx q[1];
rz(-1.7587887) q[1];
sx q[1];
rz(-0.85) q[1];
x q[2];
rz(1.24228) q[3];
sx q[3];
rz(-1.3887172) q[3];
sx q[3];
rz(-1.8329221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2386027) q[2];
sx q[2];
rz(-2.949252) q[2];
sx q[2];
rz(0.77511707) q[2];
rz(0.827968) q[3];
sx q[3];
rz(-0.29100806) q[3];
sx q[3];
rz(-2.0194139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5489952) q[0];
sx q[0];
rz(-2.7153375) q[0];
sx q[0];
rz(-0.098408498) q[0];
rz(-1.9495643) q[1];
sx q[1];
rz(-1.8076618) q[1];
sx q[1];
rz(-2.5820406) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10683051) q[0];
sx q[0];
rz(-0.97495279) q[0];
sx q[0];
rz(1.4164657) q[0];
rz(-pi) q[1];
rz(-0.19182972) q[2];
sx q[2];
rz(-0.26608135) q[2];
sx q[2];
rz(-1.9907469) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.050484867) q[1];
sx q[1];
rz(-0.15639601) q[1];
sx q[1];
rz(-0.97538235) q[1];
rz(-2.8758994) q[3];
sx q[3];
rz(-0.64722792) q[3];
sx q[3];
rz(-0.94097394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.45903912) q[2];
sx q[2];
rz(-1.8457396) q[2];
sx q[2];
rz(2.7977978) q[2];
rz(-2.5750459) q[3];
sx q[3];
rz(-2.6930801) q[3];
sx q[3];
rz(-2.6678273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.375181) q[0];
sx q[0];
rz(-1.8090929) q[0];
sx q[0];
rz(2.4108316) q[0];
rz(-0.14239755) q[1];
sx q[1];
rz(-1.2700894) q[1];
sx q[1];
rz(0.87160814) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14411892) q[0];
sx q[0];
rz(-1.561957) q[0];
sx q[0];
rz(0.075450443) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4191188) q[2];
sx q[2];
rz(-1.9722087) q[2];
sx q[2];
rz(-0.36827189) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.034681319) q[1];
sx q[1];
rz(-2.2195663) q[1];
sx q[1];
rz(1.1120863) q[1];
rz(-1.3685162) q[3];
sx q[3];
rz(-0.61198046) q[3];
sx q[3];
rz(-2.4806541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4153851) q[2];
sx q[2];
rz(-2.0724847) q[2];
sx q[2];
rz(2.0020206) q[2];
rz(-1.6428044) q[3];
sx q[3];
rz(-0.39396861) q[3];
sx q[3];
rz(-2.22877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.9386439) q[0];
sx q[0];
rz(-1.6904172) q[0];
sx q[0];
rz(1.9198445) q[0];
rz(-0.16601673) q[1];
sx q[1];
rz(-1.8211726) q[1];
sx q[1];
rz(1.5244012) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.071951889) q[0];
sx q[0];
rz(-1.4884293) q[0];
sx q[0];
rz(-0.029043341) q[0];
rz(-pi) q[1];
rz(0.35877123) q[2];
sx q[2];
rz(-0.504474) q[2];
sx q[2];
rz(-2.3141253) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7414654) q[1];
sx q[1];
rz(-1.8997846) q[1];
sx q[1];
rz(2.9743183) q[1];
x q[2];
rz(0.41096656) q[3];
sx q[3];
rz(-0.68115679) q[3];
sx q[3];
rz(0.70511234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.021585492) q[2];
sx q[2];
rz(-1.6759796) q[2];
sx q[2];
rz(-0.35153708) q[2];
rz(-1.0567788) q[3];
sx q[3];
rz(-2.6119699) q[3];
sx q[3];
rz(-0.74469152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4979424) q[0];
sx q[0];
rz(-2.239776) q[0];
sx q[0];
rz(-1.836401) q[0];
rz(-2.7611043) q[1];
sx q[1];
rz(-2.0996129) q[1];
sx q[1];
rz(0.25340733) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2449269) q[0];
sx q[0];
rz(-1.9650808) q[0];
sx q[0];
rz(3.0299597) q[0];
rz(-pi) q[1];
rz(2.017574) q[2];
sx q[2];
rz(-0.3728711) q[2];
sx q[2];
rz(2.1602221) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3897755) q[1];
sx q[1];
rz(-2.7174065) q[1];
sx q[1];
rz(2.5245689) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9990342) q[3];
sx q[3];
rz(-2.1381452) q[3];
sx q[3];
rz(-0.67728562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.59166756) q[2];
sx q[2];
rz(-0.88576907) q[2];
sx q[2];
rz(-0.5029451) q[2];
rz(-2.2425966) q[3];
sx q[3];
rz(-1.8476202) q[3];
sx q[3];
rz(-1.1635273) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1702561) q[0];
sx q[0];
rz(-1.5383056) q[0];
sx q[0];
rz(-2.8785895) q[0];
rz(-2.4304216) q[1];
sx q[1];
rz(-2.053459) q[1];
sx q[1];
rz(-1.4278535) q[1];
rz(0.28094963) q[2];
sx q[2];
rz(-1.8200257) q[2];
sx q[2];
rz(-0.65489468) q[2];
rz(-2.2181702) q[3];
sx q[3];
rz(-0.9265201) q[3];
sx q[3];
rz(-2.0410782) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
