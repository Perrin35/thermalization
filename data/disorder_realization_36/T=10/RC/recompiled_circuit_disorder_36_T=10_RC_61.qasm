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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5776423) q[0];
sx q[0];
rz(-2.2221774) q[0];
sx q[0];
rz(1.3629859) q[0];
rz(-pi) q[1];
rz(2.3876786) q[2];
sx q[2];
rz(-2.4587817) q[2];
sx q[2];
rz(-2.4726601) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2540993) q[1];
sx q[1];
rz(-1.8389529) q[1];
sx q[1];
rz(1.9894132) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64896119) q[3];
sx q[3];
rz(-2.0994086) q[3];
sx q[3];
rz(-1.987207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.704533) q[2];
sx q[2];
rz(-1.5960863) q[2];
sx q[2];
rz(0.68721592) q[2];
rz(1.0152738) q[3];
sx q[3];
rz(-1.3736558) q[3];
sx q[3];
rz(-3.0190873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17094831) q[0];
sx q[0];
rz(-2.0630554) q[0];
sx q[0];
rz(-1.8815536) q[0];
rz(2.1353703) q[1];
sx q[1];
rz(-2.1496014) q[1];
sx q[1];
rz(0.84567436) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9189008) q[0];
sx q[0];
rz(-2.0881623) q[0];
sx q[0];
rz(1.067724) q[0];
rz(-pi) q[1];
x q[1];
rz(0.13375608) q[2];
sx q[2];
rz(-1.4417366) q[2];
sx q[2];
rz(1.918902) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.37296346) q[1];
sx q[1];
rz(-2.0634723) q[1];
sx q[1];
rz(1.782934) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2250697) q[3];
sx q[3];
rz(-2.407981) q[3];
sx q[3];
rz(-0.77378002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.78850293) q[2];
sx q[2];
rz(-2.916009) q[2];
sx q[2];
rz(-0.4804002) q[2];
rz(1.3530312) q[3];
sx q[3];
rz(-2.0856817) q[3];
sx q[3];
rz(-1.1876748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95124328) q[0];
sx q[0];
rz(-0.23878637) q[0];
sx q[0];
rz(0.7730661) q[0];
rz(-0.13126016) q[1];
sx q[1];
rz(-1.2845598) q[1];
sx q[1];
rz(-2.0551596) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.438293) q[0];
sx q[0];
rz(-2.6801531) q[0];
sx q[0];
rz(1.5745387) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.92347446) q[2];
sx q[2];
rz(-1.662865) q[2];
sx q[2];
rz(-1.1084686) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.8311027) q[1];
sx q[1];
rz(-1.5655787) q[1];
sx q[1];
rz(-2.8527841) q[1];
rz(-pi) q[2];
rz(1.9757189) q[3];
sx q[3];
rz(-1.7305264) q[3];
sx q[3];
rz(-0.83042849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1290258) q[2];
sx q[2];
rz(-1.7029224) q[2];
sx q[2];
rz(-1.770299) q[2];
rz(-0.38315547) q[3];
sx q[3];
rz(-1.2569191) q[3];
sx q[3];
rz(-0.80254054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6894158) q[0];
sx q[0];
rz(-1.2503662) q[0];
sx q[0];
rz(-0.048359811) q[0];
rz(-2.9776749) q[1];
sx q[1];
rz(-2.7719031) q[1];
sx q[1];
rz(-1.4455459) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1235385) q[0];
sx q[0];
rz(-2.2211383) q[0];
sx q[0];
rz(1.7975848) q[0];
rz(-pi) q[1];
rz(-0.21638685) q[2];
sx q[2];
rz(-1.6587703) q[2];
sx q[2];
rz(2.3073334) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2913937) q[1];
sx q[1];
rz(-1.7899917) q[1];
sx q[1];
rz(1.5047969) q[1];
rz(-pi) q[2];
rz(-1.0936071) q[3];
sx q[3];
rz(-0.573199) q[3];
sx q[3];
rz(1.5447865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0754898) q[2];
sx q[2];
rz(-1.7843856) q[2];
sx q[2];
rz(1.0243105) q[2];
rz(1.5284437) q[3];
sx q[3];
rz(-1.5214835) q[3];
sx q[3];
rz(-0.23322341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4399453) q[0];
sx q[0];
rz(-2.3174536) q[0];
sx q[0];
rz(-1.8540927) q[0];
rz(0.31907407) q[1];
sx q[1];
rz(-1.5417475) q[1];
sx q[1];
rz(0.85420001) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4860977) q[0];
sx q[0];
rz(-0.71394074) q[0];
sx q[0];
rz(-1.5833202) q[0];
rz(-1.008026) q[2];
sx q[2];
rz(-0.78283435) q[2];
sx q[2];
rz(-0.40751878) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2687208) q[1];
sx q[1];
rz(-2.1495021) q[1];
sx q[1];
rz(0.89923664) q[1];
rz(-2.1339995) q[3];
sx q[3];
rz(-2.0493205) q[3];
sx q[3];
rz(2.0625045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1896818) q[2];
sx q[2];
rz(-2.5482735) q[2];
sx q[2];
rz(-0.577315) q[2];
rz(2.632085) q[3];
sx q[3];
rz(-2.7039492) q[3];
sx q[3];
rz(-2.234941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.0034870738) q[0];
sx q[0];
rz(-1.071799) q[0];
sx q[0];
rz(-3.0694718) q[0];
rz(-2.0347319) q[1];
sx q[1];
rz(-0.51263428) q[1];
sx q[1];
rz(-0.12621005) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0059347) q[0];
sx q[0];
rz(-1.5392443) q[0];
sx q[0];
rz(2.9289782) q[0];
rz(-pi) q[1];
rz(-1.948914) q[2];
sx q[2];
rz(-0.83231229) q[2];
sx q[2];
rz(0.2371012) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0969442) q[1];
sx q[1];
rz(-1.3828039) q[1];
sx q[1];
rz(-0.85) q[1];
rz(-pi) q[2];
rz(2.0893655) q[3];
sx q[3];
rz(-0.37399451) q[3];
sx q[3];
rz(0.75043375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2386027) q[2];
sx q[2];
rz(-2.949252) q[2];
sx q[2];
rz(2.3664756) q[2];
rz(-0.827968) q[3];
sx q[3];
rz(-0.29100806) q[3];
sx q[3];
rz(-1.1221788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59259748) q[0];
sx q[0];
rz(-0.42625517) q[0];
sx q[0];
rz(3.0431842) q[0];
rz(1.9495643) q[1];
sx q[1];
rz(-1.3339309) q[1];
sx q[1];
rz(0.55955204) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3768809) q[0];
sx q[0];
rz(-1.6983713) q[0];
sx q[0];
rz(-0.60140951) q[0];
x q[1];
rz(-2.9497629) q[2];
sx q[2];
rz(-0.26608135) q[2];
sx q[2];
rz(-1.1508458) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2110062) q[1];
sx q[1];
rz(-1.6582656) q[1];
sx q[1];
rz(1.4409815) q[1];
x q[2];
rz(1.3748752) q[3];
sx q[3];
rz(-0.94983263) q[3];
sx q[3];
rz(1.8718815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.45903912) q[2];
sx q[2];
rz(-1.295853) q[2];
sx q[2];
rz(0.34379488) q[2];
rz(0.5665468) q[3];
sx q[3];
rz(-2.6930801) q[3];
sx q[3];
rz(0.47376537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.7664117) q[0];
sx q[0];
rz(-1.8090929) q[0];
sx q[0];
rz(-0.73076105) q[0];
rz(-0.14239755) q[1];
sx q[1];
rz(-1.8715033) q[1];
sx q[1];
rz(-0.87160814) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14411892) q[0];
sx q[0];
rz(-1.561957) q[0];
sx q[0];
rz(-0.075450443) q[0];
rz(-pi) q[1];
x q[1];
rz(0.34198728) q[2];
sx q[2];
rz(-0.42765289) q[2];
sx q[2];
rz(-2.4004186) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.65054446) q[1];
sx q[1];
rz(-2.3666413) q[1];
sx q[1];
rz(0.52853711) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96846795) q[3];
sx q[3];
rz(-1.6864711) q[3];
sx q[3];
rz(-2.065421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4153851) q[2];
sx q[2];
rz(-1.0691079) q[2];
sx q[2];
rz(-1.139572) q[2];
rz(-1.4987882) q[3];
sx q[3];
rz(-0.39396861) q[3];
sx q[3];
rz(2.22877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9386439) q[0];
sx q[0];
rz(-1.4511755) q[0];
sx q[0];
rz(-1.9198445) q[0];
rz(-0.16601673) q[1];
sx q[1];
rz(-1.8211726) q[1];
sx q[1];
rz(-1.6171914) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.411392) q[0];
sx q[0];
rz(-0.087326614) q[0];
sx q[0];
rz(1.2325531) q[0];
x q[1];
rz(0.47713251) q[2];
sx q[2];
rz(-1.7413365) q[2];
sx q[2];
rz(1.0605304) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9164239) q[1];
sx q[1];
rz(-1.7290219) q[1];
sx q[1];
rz(-1.2374864) q[1];
x q[2];
rz(0.41096656) q[3];
sx q[3];
rz(-0.68115679) q[3];
sx q[3];
rz(-2.4364803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.021585492) q[2];
sx q[2];
rz(-1.6759796) q[2];
sx q[2];
rz(0.35153708) q[2];
rz(-1.0567788) q[3];
sx q[3];
rz(-2.6119699) q[3];
sx q[3];
rz(-0.74469152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64365023) q[0];
sx q[0];
rz(-0.90181667) q[0];
sx q[0];
rz(1.3051916) q[0];
rz(2.7611043) q[1];
sx q[1];
rz(-1.0419798) q[1];
sx q[1];
rz(0.25340733) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6127374) q[0];
sx q[0];
rz(-2.7326072) q[0];
sx q[0];
rz(-1.3091875) q[0];
rz(-pi) q[1];
rz(-1.2316522) q[2];
sx q[2];
rz(-1.412743) q[2];
sx q[2];
rz(-2.9718285) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7518172) q[1];
sx q[1];
rz(-2.7174065) q[1];
sx q[1];
rz(-0.61702375) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1425584) q[3];
sx q[3];
rz(-1.0034475) q[3];
sx q[3];
rz(2.464307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5499251) q[2];
sx q[2];
rz(-2.2558236) q[2];
sx q[2];
rz(0.5029451) q[2];
rz(0.89899603) q[3];
sx q[3];
rz(-1.2939724) q[3];
sx q[3];
rz(1.1635273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1702561) q[0];
sx q[0];
rz(-1.5383056) q[0];
sx q[0];
rz(-2.8785895) q[0];
rz(-0.7111711) q[1];
sx q[1];
rz(-1.0881337) q[1];
sx q[1];
rz(1.7137391) q[1];
rz(2.860643) q[2];
sx q[2];
rz(-1.3215669) q[2];
sx q[2];
rz(2.486698) q[2];
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