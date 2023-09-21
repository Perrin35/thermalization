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
rz(-0.87969765) q[0];
sx q[0];
rz(-1.4059773) q[0];
sx q[0];
rz(0.66189712) q[0];
rz(-0.75391407) q[2];
sx q[2];
rz(-2.4587817) q[2];
sx q[2];
rz(-2.4726601) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.56596245) q[1];
sx q[1];
rz(-1.168025) q[1];
sx q[1];
rz(0.29213841) q[1];
rz(-0.9382117) q[3];
sx q[3];
rz(-2.1198366) q[3];
sx q[3];
rz(-3.0905746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.43705964) q[2];
sx q[2];
rz(-1.5455064) q[2];
sx q[2];
rz(2.4543767) q[2];
rz(-1.0152738) q[3];
sx q[3];
rz(-1.7679368) q[3];
sx q[3];
rz(0.12250531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17094831) q[0];
sx q[0];
rz(-2.0630554) q[0];
sx q[0];
rz(-1.2600391) q[0];
rz(-2.1353703) q[1];
sx q[1];
rz(-2.1496014) q[1];
sx q[1];
rz(2.2959183) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0613522) q[0];
sx q[0];
rz(-2.4363359) q[0];
sx q[0];
rz(0.7028701) q[0];
rz(-pi) q[1];
rz(3.0078366) q[2];
sx q[2];
rz(-1.4417366) q[2];
sx q[2];
rz(-1.918902) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.37296346) q[1];
sx q[1];
rz(-1.0781204) q[1];
sx q[1];
rz(-1.3586587) q[1];
x q[2];
rz(0.50176974) q[3];
sx q[3];
rz(-2.1309149) q[3];
sx q[3];
rz(-0.027651699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.78850293) q[2];
sx q[2];
rz(-2.916009) q[2];
sx q[2];
rz(2.6611924) q[2];
rz(-1.3530312) q[3];
sx q[3];
rz(-1.055911) q[3];
sx q[3];
rz(1.9539179) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1903494) q[0];
sx q[0];
rz(-0.23878637) q[0];
sx q[0];
rz(0.7730661) q[0];
rz(-0.13126016) q[1];
sx q[1];
rz(-1.2845598) q[1];
sx q[1];
rz(1.0864331) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4341136) q[0];
sx q[0];
rz(-1.1093603) q[0];
sx q[0];
rz(0.001860851) q[0];
x q[1];
rz(-0.11523192) q[2];
sx q[2];
rz(-0.92667246) q[2];
sx q[2];
rz(-0.53172058) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4003488) q[1];
sx q[1];
rz(-1.2819918) q[1];
sx q[1];
rz(1.5762394) q[1];
x q[2];
rz(-2.9680786) q[3];
sx q[3];
rz(-1.1713235) q[3];
sx q[3];
rz(-2.3331593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1290258) q[2];
sx q[2];
rz(-1.4386703) q[2];
sx q[2];
rz(-1.3712937) q[2];
rz(-0.38315547) q[3];
sx q[3];
rz(-1.8846735) q[3];
sx q[3];
rz(0.80254054) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6894158) q[0];
sx q[0];
rz(-1.8912264) q[0];
sx q[0];
rz(-0.048359811) q[0];
rz(2.9776749) q[1];
sx q[1];
rz(-2.7719031) q[1];
sx q[1];
rz(1.4455459) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4500344) q[0];
sx q[0];
rz(-1.7507179) q[0];
sx q[0];
rz(0.66288373) q[0];
rz(-0.38979001) q[2];
sx q[2];
rz(-0.23332694) q[2];
sx q[2];
rz(2.0248272) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9963035) q[1];
sx q[1];
rz(-2.9128296) q[1];
sx q[1];
rz(0.28782515) q[1];
x q[2];
rz(2.0479855) q[3];
sx q[3];
rz(-2.5683937) q[3];
sx q[3];
rz(1.5968061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.066102862) q[2];
sx q[2];
rz(-1.7843856) q[2];
sx q[2];
rz(2.1172822) q[2];
rz(1.6131489) q[3];
sx q[3];
rz(-1.6201092) q[3];
sx q[3];
rz(-0.23322341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
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
rz(-1.5417475) q[1];
sx q[1];
rz(0.85420001) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6720649) q[0];
sx q[0];
rz(-0.85692353) q[0];
sx q[0];
rz(-3.1307427) q[0];
rz(-0.87128432) q[2];
sx q[2];
rz(-1.185002) q[2];
sx q[2];
rz(-2.3988349) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.90480587) q[1];
sx q[1];
rz(-2.2854837) q[1];
sx q[1];
rz(0.76101117) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1339995) q[3];
sx q[3];
rz(-1.0922722) q[3];
sx q[3];
rz(1.0790881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1896818) q[2];
sx q[2];
rz(-0.59331912) q[2];
sx q[2];
rz(-2.5642776) q[2];
rz(-2.632085) q[3];
sx q[3];
rz(-0.43764344) q[3];
sx q[3];
rz(-2.234941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0034870738) q[0];
sx q[0];
rz(-2.0697937) q[0];
sx q[0];
rz(-3.0694718) q[0];
rz(-1.1068608) q[1];
sx q[1];
rz(-2.6289584) q[1];
sx q[1];
rz(3.0153826) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70996767) q[0];
sx q[0];
rz(-2.9266848) q[0];
sx q[0];
rz(0.14847319) q[0];
x q[1];
rz(0.77504471) q[2];
sx q[2];
rz(-1.8473052) q[2];
sx q[2];
rz(1.0724049) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8773552) q[1];
sx q[1];
rz(-0.74062956) q[1];
sx q[1];
rz(1.8514368) q[1];
rz(0.19212171) q[3];
sx q[3];
rz(-1.8936833) q[3];
sx q[3];
rz(0.20048143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.90298992) q[2];
sx q[2];
rz(-2.949252) q[2];
sx q[2];
rz(-2.3664756) q[2];
rz(-0.827968) q[3];
sx q[3];
rz(-0.29100806) q[3];
sx q[3];
rz(2.0194139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(0.59259748) q[0];
sx q[0];
rz(-0.42625517) q[0];
sx q[0];
rz(0.098408498) q[0];
rz(-1.9495643) q[1];
sx q[1];
rz(-1.3339309) q[1];
sx q[1];
rz(-0.55955204) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10683051) q[0];
sx q[0];
rz(-0.97495279) q[0];
sx q[0];
rz(-1.725127) q[0];
x q[1];
rz(2.9497629) q[2];
sx q[2];
rz(-0.26608135) q[2];
sx q[2];
rz(-1.9907469) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4899788) q[1];
sx q[1];
rz(-1.4414806) q[1];
sx q[1];
rz(3.053385) q[1];
rz(2.5114602) q[3];
sx q[3];
rz(-1.4117985) q[3];
sx q[3];
rz(-0.4160479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.45903912) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7664117) q[0];
sx q[0];
rz(-1.8090929) q[0];
sx q[0];
rz(0.73076105) q[0];
rz(2.9991951) q[1];
sx q[1];
rz(-1.2700894) q[1];
sx q[1];
rz(-2.2699845) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9974737) q[0];
sx q[0];
rz(-1.561957) q[0];
sx q[0];
rz(-0.075450443) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7360104) q[2];
sx q[2];
rz(-1.7103346) q[2];
sx q[2];
rz(1.1428733) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4910482) q[1];
sx q[1];
rz(-0.77495134) q[1];
sx q[1];
rz(2.6130555) q[1];
rz(-3.001508) q[3];
sx q[3];
rz(-0.9730556) q[3];
sx q[3];
rz(-0.41543451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4153851) q[2];
sx q[2];
rz(-1.0691079) q[2];
sx q[2];
rz(2.0020206) q[2];
rz(1.6428044) q[3];
sx q[3];
rz(-0.39396861) q[3];
sx q[3];
rz(-0.9128226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
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
rz(-1.8211726) q[1];
sx q[1];
rz(-1.6171914) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.411392) q[0];
sx q[0];
rz(-0.087326614) q[0];
sx q[0];
rz(-1.9090396) q[0];
x q[1];
rz(-2.6644601) q[2];
sx q[2];
rz(-1.4002561) q[2];
sx q[2];
rz(2.0810623) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0601378) q[1];
sx q[1];
rz(-2.7739035) q[1];
sx q[1];
rz(1.1170438) q[1];
x q[2];
rz(-1.2576305) q[3];
sx q[3];
rz(-2.1861665) q[3];
sx q[3];
rz(0.19389158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.021585492) q[2];
sx q[2];
rz(-1.6759796) q[2];
sx q[2];
rz(-0.35153708) q[2];
rz(-1.0567788) q[3];
sx q[3];
rz(-2.6119699) q[3];
sx q[3];
rz(2.3969011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64365023) q[0];
sx q[0];
rz(-2.239776) q[0];
sx q[0];
rz(-1.836401) q[0];
rz(-2.7611043) q[1];
sx q[1];
rz(-2.0996129) q[1];
sx q[1];
rz(-2.8881853) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8587592) q[0];
sx q[0];
rz(-1.6738322) q[0];
sx q[0];
rz(-1.1742924) q[0];
rz(0.16742736) q[2];
sx q[2];
rz(-1.9055467) q[2];
sx q[2];
rz(-1.4565005) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3897755) q[1];
sx q[1];
rz(-2.7174065) q[1];
sx q[1];
rz(2.5245689) q[1];
rz(-1.1425584) q[3];
sx q[3];
rz(-2.1381452) q[3];
sx q[3];
rz(-2.464307) q[3];
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
rz(-0.89899603) q[3];
sx q[3];
rz(-1.8476202) q[3];
sx q[3];
rz(-1.9780654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9713365) q[0];
sx q[0];
rz(-1.6032871) q[0];
sx q[0];
rz(0.26300318) q[0];
rz(-0.7111711) q[1];
sx q[1];
rz(-1.0881337) q[1];
sx q[1];
rz(1.7137391) q[1];
rz(1.3118369) q[2];
sx q[2];
rz(-1.8428409) q[2];
sx q[2];
rz(0.98696282) q[2];
rz(-2.3861804) q[3];
sx q[3];
rz(-1.0676386) q[3];
sx q[3];
rz(-0.044015351) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
