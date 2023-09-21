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
rz(2.5640092) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89864697) q[0];
sx q[0];
rz(-2.4624914) q[0];
sx q[0];
rz(-2.8773017) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0787813) q[2];
sx q[2];
rz(-2.0487818) q[2];
sx q[2];
rz(-2.9302772) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2204199) q[1];
sx q[1];
rz(-2.6487659) q[1];
sx q[1];
rz(-0.9763896) q[1];
x q[2];
rz(0.76831423) q[3];
sx q[3];
rz(-0.81211219) q[3];
sx q[3];
rz(1.0031568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.43705964) q[2];
sx q[2];
rz(-1.5455064) q[2];
sx q[2];
rz(2.4543767) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9706443) q[0];
sx q[0];
rz(-2.0630554) q[0];
sx q[0];
rz(-1.2600391) q[0];
rz(2.1353703) q[1];
sx q[1];
rz(-2.1496014) q[1];
sx q[1];
rz(-2.2959183) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0802404) q[0];
sx q[0];
rz(-2.4363359) q[0];
sx q[0];
rz(2.4387226) q[0];
rz(-0.13375608) q[2];
sx q[2];
rz(-1.4417366) q[2];
sx q[2];
rz(-1.918902) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0963124) q[1];
sx q[1];
rz(-1.3842061) q[1];
sx q[1];
rz(0.50218302) q[1];
x q[2];
rz(-2.1917079) q[3];
sx q[3];
rz(-1.1511027) q[3];
sx q[3];
rz(-1.8267531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.78850293) q[2];
sx q[2];
rz(-0.22558364) q[2];
sx q[2];
rz(2.6611924) q[2];
rz(1.3530312) q[3];
sx q[3];
rz(-2.0856817) q[3];
sx q[3];
rz(-1.1876748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1903494) q[0];
sx q[0];
rz(-2.9028063) q[0];
sx q[0];
rz(0.7730661) q[0];
rz(-0.13126016) q[1];
sx q[1];
rz(-1.8570329) q[1];
sx q[1];
rz(2.0551596) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.438293) q[0];
sx q[0];
rz(-2.6801531) q[0];
sx q[0];
rz(-1.5745387) q[0];
rz(1.7227206) q[2];
sx q[2];
rz(-0.65290367) q[2];
sx q[2];
rz(-0.34130794) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.72213458) q[1];
sx q[1];
rz(-0.28885435) q[1];
sx q[1];
rz(3.1232749) q[1];
rz(-pi) q[2];
rz(2.9680786) q[3];
sx q[3];
rz(-1.9702692) q[3];
sx q[3];
rz(-2.3331593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0125668) q[2];
sx q[2];
rz(-1.4386703) q[2];
sx q[2];
rz(-1.770299) q[2];
rz(2.7584372) q[3];
sx q[3];
rz(-1.8846735) q[3];
sx q[3];
rz(0.80254054) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4521769) q[0];
sx q[0];
rz(-1.2503662) q[0];
sx q[0];
rz(0.048359811) q[0];
rz(-0.16391779) q[1];
sx q[1];
rz(-0.36968958) q[1];
sx q[1];
rz(1.6960467) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1235385) q[0];
sx q[0];
rz(-2.2211383) q[0];
sx q[0];
rz(-1.7975848) q[0];
x q[1];
rz(2.9252058) q[2];
sx q[2];
rz(-1.6587703) q[2];
sx q[2];
rz(-0.83425922) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2913937) q[1];
sx q[1];
rz(-1.351601) q[1];
sx q[1];
rz(-1.6367957) q[1];
rz(-2.0479855) q[3];
sx q[3];
rz(-0.573199) q[3];
sx q[3];
rz(1.5968061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0754898) q[2];
sx q[2];
rz(-1.3572071) q[2];
sx q[2];
rz(-2.1172822) q[2];
rz(1.6131489) q[3];
sx q[3];
rz(-1.5214835) q[3];
sx q[3];
rz(0.23322341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7016474) q[0];
sx q[0];
rz(-2.3174536) q[0];
sx q[0];
rz(-1.8540927) q[0];
rz(-0.31907407) q[1];
sx q[1];
rz(-1.5417475) q[1];
sx q[1];
rz(-0.85420001) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4860977) q[0];
sx q[0];
rz(-2.4276519) q[0];
sx q[0];
rz(-1.5833202) q[0];
rz(-0.48798497) q[2];
sx q[2];
rz(-2.2099566) q[2];
sx q[2];
rz(-2.0069063) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0336049) q[1];
sx q[1];
rz(-1.0228979) q[1];
sx q[1];
rz(2.4461436) q[1];
rz(-pi) q[2];
rz(-2.5913127) q[3];
sx q[3];
rz(-2.0645421) q[3];
sx q[3];
rz(0.20875904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1896818) q[2];
sx q[2];
rz(-0.59331912) q[2];
sx q[2];
rz(2.5642776) q[2];
rz(2.632085) q[3];
sx q[3];
rz(-2.7039492) q[3];
sx q[3];
rz(0.90665162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0034870738) q[0];
sx q[0];
rz(-2.0697937) q[0];
sx q[0];
rz(-0.072120897) q[0];
rz(2.0347319) q[1];
sx q[1];
rz(-0.51263428) q[1];
sx q[1];
rz(0.12621005) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70996767) q[0];
sx q[0];
rz(-2.9266848) q[0];
sx q[0];
rz(0.14847319) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3665479) q[2];
sx q[2];
rz(-1.8473052) q[2];
sx q[2];
rz(1.0724049) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5050161) q[1];
sx q[1];
rz(-0.86537433) q[1];
sx q[1];
rz(0.24800639) q[1];
rz(-pi) q[2];
rz(-1.0522271) q[3];
sx q[3];
rz(-0.37399451) q[3];
sx q[3];
rz(-2.3911589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2386027) q[2];
sx q[2];
rz(-2.949252) q[2];
sx q[2];
rz(0.77511707) q[2];
rz(-0.827968) q[3];
sx q[3];
rz(-0.29100806) q[3];
sx q[3];
rz(-1.1221788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5489952) q[0];
sx q[0];
rz(-2.7153375) q[0];
sx q[0];
rz(-3.0431842) q[0];
rz(-1.9495643) q[1];
sx q[1];
rz(-1.8076618) q[1];
sx q[1];
rz(0.55955204) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10683051) q[0];
sx q[0];
rz(-2.1666399) q[0];
sx q[0];
rz(-1.725127) q[0];
x q[1];
rz(1.6227116) q[2];
sx q[2];
rz(-1.3097109) q[2];
sx q[2];
rz(1.3494929) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0911078) q[1];
sx q[1];
rz(-0.15639601) q[1];
sx q[1];
rz(0.97538235) q[1];
x q[2];
rz(0.26569326) q[3];
sx q[3];
rz(-2.4943647) q[3];
sx q[3];
rz(-2.2006187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.45903912) q[2];
sx q[2];
rz(-1.8457396) q[2];
sx q[2];
rz(-0.34379488) q[2];
rz(-0.5665468) q[3];
sx q[3];
rz(-0.44851258) q[3];
sx q[3];
rz(0.47376537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.7664117) q[0];
sx q[0];
rz(-1.8090929) q[0];
sx q[0];
rz(0.73076105) q[0];
rz(2.9991951) q[1];
sx q[1];
rz(-1.8715033) q[1];
sx q[1];
rz(2.2699845) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7155834) q[0];
sx q[0];
rz(-1.4953488) q[0];
sx q[0];
rz(1.5619318) q[0];
rz(0.40558221) q[2];
sx q[2];
rz(-1.7103346) q[2];
sx q[2];
rz(1.1428733) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.65054446) q[1];
sx q[1];
rz(-0.77495134) q[1];
sx q[1];
rz(2.6130555) q[1];
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
sx q[1];
rz(-pi/2) q[1];
rz(-0.72620755) q[2];
sx q[2];
rz(-1.0691079) q[2];
sx q[2];
rz(1.139572) q[2];
rz(1.6428044) q[3];
sx q[3];
rz(-2.747624) q[3];
sx q[3];
rz(0.9128226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20294872) q[0];
sx q[0];
rz(-1.4511755) q[0];
sx q[0];
rz(-1.2217481) q[0];
rz(-0.16601673) q[1];
sx q[1];
rz(-1.32042) q[1];
sx q[1];
rz(1.6171914) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6451384) q[0];
sx q[0];
rz(-1.5997412) q[0];
sx q[0];
rz(1.653198) q[0];
x q[1];
rz(-1.379307) q[2];
sx q[2];
rz(-1.101149) q[2];
sx q[2];
rz(2.7188403) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0601378) q[1];
sx q[1];
rz(-2.7739035) q[1];
sx q[1];
rz(1.1170438) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2576305) q[3];
sx q[3];
rz(-2.1861665) q[3];
sx q[3];
rz(-0.19389158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1200072) q[2];
sx q[2];
rz(-1.6759796) q[2];
sx q[2];
rz(2.7900556) q[2];
rz(1.0567788) q[3];
sx q[3];
rz(-0.52962279) q[3];
sx q[3];
rz(2.3969011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64365023) q[0];
sx q[0];
rz(-0.90181667) q[0];
sx q[0];
rz(-1.3051916) q[0];
rz(-2.7611043) q[1];
sx q[1];
rz(-1.0419798) q[1];
sx q[1];
rz(-0.25340733) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8966658) q[0];
sx q[0];
rz(-1.9650808) q[0];
sx q[0];
rz(3.0299597) q[0];
rz(-pi) q[1];
rz(2.017574) q[2];
sx q[2];
rz(-2.7687216) q[2];
sx q[2];
rz(0.98137059) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.39292654) q[1];
sx q[1];
rz(-1.8112507) q[1];
sx q[1];
rz(2.7886831) q[1];
rz(1.1425584) q[3];
sx q[3];
rz(-2.1381452) q[3];
sx q[3];
rz(-0.67728562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5499251) q[2];
sx q[2];
rz(-0.88576907) q[2];
sx q[2];
rz(-0.5029451) q[2];
rz(-2.2425966) q[3];
sx q[3];
rz(-1.2939724) q[3];
sx q[3];
rz(1.1635273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(1.8297557) q[2];
sx q[2];
rz(-1.2987518) q[2];
sx q[2];
rz(-2.1546298) q[2];
rz(0.75541227) q[3];
sx q[3];
rz(-1.0676386) q[3];
sx q[3];
rz(-0.044015351) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];