OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7918474) q[0];
sx q[0];
rz(-2.1748769) q[0];
sx q[0];
rz(-0.36122394) q[0];
rz(-0.53520441) q[1];
sx q[1];
rz(-1.1517795) q[1];
sx q[1];
rz(1.8395543) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0326234) q[0];
sx q[0];
rz(-1.291467) q[0];
sx q[0];
rz(0.33786122) q[0];
x q[1];
rz(0.7187219) q[2];
sx q[2];
rz(-1.5437922) q[2];
sx q[2];
rz(0.57098963) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.227461) q[1];
sx q[1];
rz(-1.1644851) q[1];
sx q[1];
rz(-1.3106457) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5005712) q[3];
sx q[3];
rz(-1.3596137) q[3];
sx q[3];
rz(0.9259895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5754622) q[2];
sx q[2];
rz(-2.8464816) q[2];
sx q[2];
rz(-2.5556514) q[2];
rz(-1.5015886) q[3];
sx q[3];
rz(-1.274704) q[3];
sx q[3];
rz(1.2007825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-2.5457299) q[0];
sx q[0];
rz(-1.0965309) q[0];
sx q[0];
rz(-2.0134266) q[0];
rz(0.7496756) q[1];
sx q[1];
rz(-1.3355037) q[1];
sx q[1];
rz(-1.658176) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7455341) q[0];
sx q[0];
rz(-2.0307956) q[0];
sx q[0];
rz(0.7454304) q[0];
rz(-2.5407578) q[2];
sx q[2];
rz(-1.877583) q[2];
sx q[2];
rz(2.6400849) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7871889) q[1];
sx q[1];
rz(-1.8285079) q[1];
sx q[1];
rz(-0.54097213) q[1];
rz(-pi) q[2];
rz(-2.800249) q[3];
sx q[3];
rz(-1.1142379) q[3];
sx q[3];
rz(-1.8648636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.38773203) q[2];
sx q[2];
rz(-0.78247491) q[2];
sx q[2];
rz(2.1369047) q[2];
rz(-3.1333771) q[3];
sx q[3];
rz(-1.7693844) q[3];
sx q[3];
rz(0.43732244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32366556) q[0];
sx q[0];
rz(-1.265047) q[0];
sx q[0];
rz(-1.0311968) q[0];
rz(-2.9464856) q[1];
sx q[1];
rz(-0.61385265) q[1];
sx q[1];
rz(-2.2574147) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3613286) q[0];
sx q[0];
rz(-1.5705646) q[0];
sx q[0];
rz(-0.0032629691) q[0];
x q[1];
rz(0.43135204) q[2];
sx q[2];
rz(-2.6425053) q[2];
sx q[2];
rz(-1.8519397) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9337195) q[1];
sx q[1];
rz(-1.6657511) q[1];
sx q[1];
rz(-1.1594561) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9724794) q[3];
sx q[3];
rz(-2.5163262) q[3];
sx q[3];
rz(2.7234944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.6241793) q[2];
sx q[2];
rz(-1.8884582) q[2];
sx q[2];
rz(-2.6873028) q[2];
rz(-1.0034466) q[3];
sx q[3];
rz(-0.61383057) q[3];
sx q[3];
rz(1.7295037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2237332) q[0];
sx q[0];
rz(-0.35686019) q[0];
sx q[0];
rz(0.80765635) q[0];
rz(1.7478906) q[1];
sx q[1];
rz(-1.423577) q[1];
sx q[1];
rz(1.4180988) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7634228) q[0];
sx q[0];
rz(-2.6900107) q[0];
sx q[0];
rz(1.9096791) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.70576422) q[2];
sx q[2];
rz(-1.4250722) q[2];
sx q[2];
rz(0.55225295) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1313626) q[1];
sx q[1];
rz(-1.5829965) q[1];
sx q[1];
rz(2.7598445) q[1];
x q[2];
rz(2.0924166) q[3];
sx q[3];
rz(-1.9289013) q[3];
sx q[3];
rz(1.1529837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.35004804) q[2];
sx q[2];
rz(-1.2309265) q[2];
sx q[2];
rz(-0.13354224) q[2];
rz(-0.40889016) q[3];
sx q[3];
rz(-0.10433993) q[3];
sx q[3];
rz(-2.1321808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23311663) q[0];
sx q[0];
rz(-0.82670832) q[0];
sx q[0];
rz(-0.46434656) q[0];
rz(2.6404479) q[1];
sx q[1];
rz(-2.8906288) q[1];
sx q[1];
rz(-2.4026925) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2621612) q[0];
sx q[0];
rz(-1.3152243) q[0];
sx q[0];
rz(-1.68925) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7064507) q[2];
sx q[2];
rz(-2.1096482) q[2];
sx q[2];
rz(1.5340027) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1170144) q[1];
sx q[1];
rz(-1.1305594) q[1];
sx q[1];
rz(-1.3280895) q[1];
x q[2];
rz(-2.49315) q[3];
sx q[3];
rz(-2.01727) q[3];
sx q[3];
rz(0.97383271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.240694) q[2];
sx q[2];
rz(-1.7186761) q[2];
sx q[2];
rz(-2.9071729) q[2];
rz(-0.33469409) q[3];
sx q[3];
rz(-2.5321999) q[3];
sx q[3];
rz(-1.4332244) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1471106) q[0];
sx q[0];
rz(-0.88279498) q[0];
sx q[0];
rz(-0.91142803) q[0];
rz(1.8270739) q[1];
sx q[1];
rz(-0.85845566) q[1];
sx q[1];
rz(-1.7427157) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0056038) q[0];
sx q[0];
rz(-2.4244747) q[0];
sx q[0];
rz(0.75053458) q[0];
rz(-pi) q[1];
rz(-1.6009459) q[2];
sx q[2];
rz(-1.7099524) q[2];
sx q[2];
rz(2.2517532) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6487471) q[1];
sx q[1];
rz(-1.6688818) q[1];
sx q[1];
rz(2.0010038) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1460378) q[3];
sx q[3];
rz(-2.5638608) q[3];
sx q[3];
rz(2.7932624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2167902) q[2];
sx q[2];
rz(-2.3239467) q[2];
sx q[2];
rz(-1.7693046) q[2];
rz(-1.6670082) q[3];
sx q[3];
rz(-2.0788914) q[3];
sx q[3];
rz(0.32349989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40397662) q[0];
sx q[0];
rz(-2.1196899) q[0];
sx q[0];
rz(1.2603941) q[0];
rz(-2.8003287) q[1];
sx q[1];
rz(-0.9287467) q[1];
sx q[1];
rz(2.0492679) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72427801) q[0];
sx q[0];
rz(-1.5515741) q[0];
sx q[0];
rz(0.5165944) q[0];
rz(-0.50470074) q[2];
sx q[2];
rz(-2.7554858) q[2];
sx q[2];
rz(-2.8442755) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.87783972) q[1];
sx q[1];
rz(-1.3307376) q[1];
sx q[1];
rz(-0.93024436) q[1];
x q[2];
rz(-0.63082781) q[3];
sx q[3];
rz(-2.6371427) q[3];
sx q[3];
rz(0.40182477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1616538) q[2];
sx q[2];
rz(-1.7540163) q[2];
sx q[2];
rz(0.8650583) q[2];
rz(3.0480399) q[3];
sx q[3];
rz(-1.4253987) q[3];
sx q[3];
rz(3.0430999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7155257) q[0];
sx q[0];
rz(-0.97398296) q[0];
sx q[0];
rz(-1.4564212) q[0];
rz(0.23297019) q[1];
sx q[1];
rz(-1.3825682) q[1];
sx q[1];
rz(1.9219386) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9838966) q[0];
sx q[0];
rz(-2.2737163) q[0];
sx q[0];
rz(-1.9546175) q[0];
rz(0.43841534) q[2];
sx q[2];
rz(-1.9389956) q[2];
sx q[2];
rz(-2.6121998) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2623655) q[1];
sx q[1];
rz(-2.7272537) q[1];
sx q[1];
rz(-0.44466876) q[1];
x q[2];
rz(2.574) q[3];
sx q[3];
rz(-0.78062468) q[3];
sx q[3];
rz(-0.34263698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0435656) q[2];
sx q[2];
rz(-1.4412014) q[2];
sx q[2];
rz(-3.0230057) q[2];
rz(-0.81973997) q[3];
sx q[3];
rz(-2.6974862) q[3];
sx q[3];
rz(-2.8290101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0176625) q[0];
sx q[0];
rz(-0.95159641) q[0];
sx q[0];
rz(2.1753525) q[0];
rz(-2.7746334) q[1];
sx q[1];
rz(-0.96261135) q[1];
sx q[1];
rz(-0.56698322) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9241739) q[0];
sx q[0];
rz(-1.3710306) q[0];
sx q[0];
rz(-2.6471247) q[0];
rz(-pi) q[1];
rz(-2.1927715) q[2];
sx q[2];
rz(-0.89911956) q[2];
sx q[2];
rz(-2.8664884) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.39886576) q[1];
sx q[1];
rz(-2.4165771) q[1];
sx q[1];
rz(-0.30959495) q[1];
x q[2];
rz(2.0040599) q[3];
sx q[3];
rz(-0.71914395) q[3];
sx q[3];
rz(0.55717232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6734267) q[2];
sx q[2];
rz(-1.4665073) q[2];
sx q[2];
rz(-1.0183938) q[2];
rz(-0.23032019) q[3];
sx q[3];
rz(-0.46263328) q[3];
sx q[3];
rz(3.0751626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51909834) q[0];
sx q[0];
rz(-2.1343756) q[0];
sx q[0];
rz(-1.3315211) q[0];
rz(1.9271756) q[1];
sx q[1];
rz(-1.4446832) q[1];
sx q[1];
rz(2.725098) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40503026) q[0];
sx q[0];
rz(-2.2433537) q[0];
sx q[0];
rz(2.879438) q[0];
rz(-pi) q[1];
rz(1.824786) q[2];
sx q[2];
rz(-1.2653627) q[2];
sx q[2];
rz(-1.0930201) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.107376) q[1];
sx q[1];
rz(-1.1469081) q[1];
sx q[1];
rz(0.54383212) q[1];
x q[2];
rz(3.0399226) q[3];
sx q[3];
rz(-2.0196126) q[3];
sx q[3];
rz(3.085595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.74849621) q[2];
sx q[2];
rz(-1.8668207) q[2];
sx q[2];
rz(-1.3416802) q[2];
rz(2.0210576) q[3];
sx q[3];
rz(-1.7399961) q[3];
sx q[3];
rz(3.024658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31711598) q[0];
sx q[0];
rz(-1.8335637) q[0];
sx q[0];
rz(1.2807922) q[0];
rz(0.9723797) q[1];
sx q[1];
rz(-1.096611) q[1];
sx q[1];
rz(-0.70659804) q[1];
rz(2.7823051) q[2];
sx q[2];
rz(-2.4854598) q[2];
sx q[2];
rz(1.796915) q[2];
rz(1.6617358) q[3];
sx q[3];
rz(-1.5385689) q[3];
sx q[3];
rz(-3.0868019) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
