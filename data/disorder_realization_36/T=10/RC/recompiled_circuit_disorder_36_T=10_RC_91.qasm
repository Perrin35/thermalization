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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56395036) q[0];
sx q[0];
rz(-0.9194153) q[0];
sx q[0];
rz(-1.3629859) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6063927) q[2];
sx q[2];
rz(-1.1241962) q[2];
sx q[2];
rz(-1.5314147) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9211728) q[1];
sx q[1];
rz(-0.49282679) q[1];
sx q[1];
rz(2.1652031) q[1];
rz(-pi) q[2];
x q[2];
rz(2.203381) q[3];
sx q[3];
rz(-1.0217561) q[3];
sx q[3];
rz(3.0905746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.43705964) q[2];
sx q[2];
rz(-1.5455064) q[2];
sx q[2];
rz(2.4543767) q[2];
rz(1.0152738) q[3];
sx q[3];
rz(-1.3736558) q[3];
sx q[3];
rz(0.12250531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9706443) q[0];
sx q[0];
rz(-2.0630554) q[0];
sx q[0];
rz(1.8815536) q[0];
rz(-2.1353703) q[1];
sx q[1];
rz(-2.1496014) q[1];
sx q[1];
rz(2.2959183) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0802404) q[0];
sx q[0];
rz(-2.4363359) q[0];
sx q[0];
rz(0.7028701) q[0];
rz(0.13375608) q[2];
sx q[2];
rz(-1.6998561) q[2];
sx q[2];
rz(1.2226906) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3413275) q[1];
sx q[1];
rz(-0.53293537) q[1];
sx q[1];
rz(2.76782) q[1];
x q[2];
rz(2.1917079) q[3];
sx q[3];
rz(-1.9904899) q[3];
sx q[3];
rz(-1.8267531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3530897) q[2];
sx q[2];
rz(-2.916009) q[2];
sx q[2];
rz(-2.6611924) q[2];
rz(-1.7885615) q[3];
sx q[3];
rz(-2.0856817) q[3];
sx q[3];
rz(-1.1876748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95124328) q[0];
sx q[0];
rz(-0.23878637) q[0];
sx q[0];
rz(-0.7730661) q[0];
rz(-3.0103325) q[1];
sx q[1];
rz(-1.8570329) q[1];
sx q[1];
rz(-2.0551596) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13585424) q[0];
sx q[0];
rz(-1.5724626) q[0];
sx q[0];
rz(-2.0322331) q[0];
rz(-pi) q[1];
rz(-1.418872) q[2];
sx q[2];
rz(-2.488689) q[2];
sx q[2];
rz(-2.8002847) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.31049) q[1];
sx q[1];
rz(-1.5655787) q[1];
sx q[1];
rz(-0.28880854) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1826035) q[3];
sx q[3];
rz(-2.7079294) q[3];
sx q[3];
rz(-2.7565623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0125668) q[2];
sx q[2];
rz(-1.7029224) q[2];
sx q[2];
rz(-1.3712937) q[2];
rz(-0.38315547) q[3];
sx q[3];
rz(-1.2569191) q[3];
sx q[3];
rz(2.3390521) q[3];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4521769) q[0];
sx q[0];
rz(-1.8912264) q[0];
sx q[0];
rz(-0.048359811) q[0];
rz(2.9776749) q[1];
sx q[1];
rz(-2.7719031) q[1];
sx q[1];
rz(1.4455459) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4500344) q[0];
sx q[0];
rz(-1.7507179) q[0];
sx q[0];
rz(-2.4787089) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6608597) q[2];
sx q[2];
rz(-1.3552595) q[2];
sx q[2];
rz(2.4243674) q[2];
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
rz(-pi) q[2];
rz(-1.0501782) q[3];
sx q[3];
rz(-1.319066) q[3];
sx q[3];
rz(-0.43581918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.066102862) q[2];
sx q[2];
rz(-1.3572071) q[2];
sx q[2];
rz(1.0243105) q[2];
rz(1.5284437) q[3];
sx q[3];
rz(-1.6201092) q[3];
sx q[3];
rz(-2.9083692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0474284) q[0];
sx q[0];
rz(-1.5625956) q[0];
sx q[0];
rz(-0.8568944) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2703083) q[2];
sx q[2];
rz(-1.9565906) q[2];
sx q[2];
rz(2.3988349) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.3805815) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5913127) q[3];
sx q[3];
rz(-1.0770505) q[3];
sx q[3];
rz(0.20875904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1896818) q[2];
sx q[2];
rz(-0.59331912) q[2];
sx q[2];
rz(0.577315) q[2];
rz(2.632085) q[3];
sx q[3];
rz(-2.7039492) q[3];
sx q[3];
rz(0.90665162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0034870738) q[0];
sx q[0];
rz(-1.071799) q[0];
sx q[0];
rz(-0.072120897) q[0];
rz(-2.0347319) q[1];
sx q[1];
rz(-2.6289584) q[1];
sx q[1];
rz(0.12621005) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.135658) q[0];
sx q[0];
rz(-1.6023484) q[0];
sx q[0];
rz(-2.9289782) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1926786) q[2];
sx q[2];
rz(-2.3092804) q[2];
sx q[2];
rz(-0.2371012) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.63657657) q[1];
sx q[1];
rz(-2.2762183) q[1];
sx q[1];
rz(-0.24800639) q[1];
rz(2.9494709) q[3];
sx q[3];
rz(-1.8936833) q[3];
sx q[3];
rz(2.9411112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2386027) q[2];
sx q[2];
rz(-0.1923407) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
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
rz(1.9495643) q[1];
sx q[1];
rz(-1.8076618) q[1];
sx q[1];
rz(-0.55955204) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37724272) q[0];
sx q[0];
rz(-2.5284405) q[0];
sx q[0];
rz(2.9186547) q[0];
rz(-pi) q[1];
rz(-2.9497629) q[2];
sx q[2];
rz(-2.8755113) q[2];
sx q[2];
rz(-1.9907469) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.65161381) q[1];
sx q[1];
rz(-1.4414806) q[1];
sx q[1];
rz(0.08820769) q[1];
x q[2];
rz(-0.63013245) q[3];
sx q[3];
rz(-1.7297941) q[3];
sx q[3];
rz(-2.7255448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6825535) q[2];
sx q[2];
rz(-1.295853) q[2];
sx q[2];
rz(-0.34379488) q[2];
rz(0.5665468) q[3];
sx q[3];
rz(-0.44851258) q[3];
sx q[3];
rz(2.6678273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
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
rz(-2.9991951) q[1];
sx q[1];
rz(-1.8715033) q[1];
sx q[1];
rz(0.87160814) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5985142) q[0];
sx q[0];
rz(-0.075965479) q[0];
sx q[0];
rz(3.024858) q[0];
rz(-pi) q[1];
rz(2.7996054) q[2];
sx q[2];
rz(-0.42765289) q[2];
sx q[2];
rz(2.4004186) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4910482) q[1];
sx q[1];
rz(-2.3666413) q[1];
sx q[1];
rz(2.6130555) q[1];
x q[2];
rz(-2.1731247) q[3];
sx q[3];
rz(-1.4551216) q[3];
sx q[3];
rz(1.0761716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.72620755) q[2];
sx q[2];
rz(-1.0691079) q[2];
sx q[2];
rz(1.139572) q[2];
rz(-1.4987882) q[3];
sx q[3];
rz(-0.39396861) q[3];
sx q[3];
rz(-0.9128226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20294872) q[0];
sx q[0];
rz(-1.4511755) q[0];
sx q[0];
rz(1.2217481) q[0];
rz(2.9755759) q[1];
sx q[1];
rz(-1.32042) q[1];
sx q[1];
rz(1.6171914) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6451384) q[0];
sx q[0];
rz(-1.5997412) q[0];
sx q[0];
rz(-1.653198) q[0];
x q[1];
rz(1.7622856) q[2];
sx q[2];
rz(-1.101149) q[2];
sx q[2];
rz(2.7188403) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.081454885) q[1];
sx q[1];
rz(-2.7739035) q[1];
sx q[1];
rz(1.1170438) q[1];
x q[2];
rz(2.7306261) q[3];
sx q[3];
rz(-2.4604359) q[3];
sx q[3];
rz(0.70511234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1200072) q[2];
sx q[2];
rz(-1.465613) q[2];
sx q[2];
rz(-0.35153708) q[2];
rz(2.0848138) q[3];
sx q[3];
rz(-0.52962279) q[3];
sx q[3];
rz(-2.3969011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64365023) q[0];
sx q[0];
rz(-2.239776) q[0];
sx q[0];
rz(1.836401) q[0];
rz(-2.7611043) q[1];
sx q[1];
rz(-1.0419798) q[1];
sx q[1];
rz(2.8881853) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8966658) q[0];
sx q[0];
rz(-1.9650808) q[0];
sx q[0];
rz(-0.11163296) q[0];
x q[1];
rz(0.16742736) q[2];
sx q[2];
rz(-1.2360459) q[2];
sx q[2];
rz(1.4565005) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7518172) q[1];
sx q[1];
rz(-0.42418617) q[1];
sx q[1];
rz(-0.61702375) q[1];
rz(-0.6110544) q[3];
sx q[3];
rz(-1.9285893) q[3];
sx q[3];
rz(-0.65294453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5499251) q[2];
sx q[2];
rz(-2.2558236) q[2];
sx q[2];
rz(2.6386476) q[2];
rz(2.2425966) q[3];
sx q[3];
rz(-1.8476202) q[3];
sx q[3];
rz(-1.9780654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1702561) q[0];
sx q[0];
rz(-1.6032871) q[0];
sx q[0];
rz(0.26300318) q[0];
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
