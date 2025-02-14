OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72260296) q[0];
sx q[0];
rz(-2.498772) q[0];
sx q[0];
rz(-1.2195725) q[0];
rz(0.83528432) q[1];
sx q[1];
rz(-1.3568027) q[1];
sx q[1];
rz(-2.2396483) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57563462) q[0];
sx q[0];
rz(-2.1778641) q[0];
sx q[0];
rz(1.3835039) q[0];
rz(1.7576394) q[2];
sx q[2];
rz(-1.9451491) q[2];
sx q[2];
rz(-0.49546212) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.29034258) q[1];
sx q[1];
rz(-2.800673) q[1];
sx q[1];
rz(-1.0268289) q[1];
rz(-pi) q[2];
rz(2.2252623) q[3];
sx q[3];
rz(-2.188465) q[3];
sx q[3];
rz(2.319782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1787662) q[2];
sx q[2];
rz(-1.0054532) q[2];
sx q[2];
rz(2.315305) q[2];
rz(1.4055584) q[3];
sx q[3];
rz(-2.8345351) q[3];
sx q[3];
rz(2.3981986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58251441) q[0];
sx q[0];
rz(-0.40575108) q[0];
sx q[0];
rz(1.7412809) q[0];
rz(-2.0179613) q[1];
sx q[1];
rz(-2.7476937) q[1];
sx q[1];
rz(1.0981015) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0535677) q[0];
sx q[0];
rz(-0.19101772) q[0];
sx q[0];
rz(-2.5294526) q[0];
rz(2.3674521) q[2];
sx q[2];
rz(-2.6850651) q[2];
sx q[2];
rz(2.4696333) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4103293) q[1];
sx q[1];
rz(-1.4498693) q[1];
sx q[1];
rz(-1.4951453) q[1];
rz(-2.9750175) q[3];
sx q[3];
rz(-1.9431356) q[3];
sx q[3];
rz(-0.53322116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0850247) q[2];
sx q[2];
rz(-0.3173863) q[2];
sx q[2];
rz(-1.2497831) q[2];
rz(1.362484) q[3];
sx q[3];
rz(-2.1407849) q[3];
sx q[3];
rz(-1.6574297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63610858) q[0];
sx q[0];
rz(-0.69378575) q[0];
sx q[0];
rz(-0.27534494) q[0];
rz(-0.45922008) q[1];
sx q[1];
rz(-2.1870435) q[1];
sx q[1];
rz(-3.088248) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33063525) q[0];
sx q[0];
rz(-2.914408) q[0];
sx q[0];
rz(2.7869528) q[0];
rz(-pi) q[1];
rz(1.1098451) q[2];
sx q[2];
rz(-1.58733) q[2];
sx q[2];
rz(-1.2993882) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.054823067) q[1];
sx q[1];
rz(-0.55332282) q[1];
sx q[1];
rz(0.077484681) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0166078) q[3];
sx q[3];
rz(-0.38180581) q[3];
sx q[3];
rz(-3.0757381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.82103819) q[2];
sx q[2];
rz(-0.36131636) q[2];
sx q[2];
rz(-1.5041171) q[2];
rz(1.5004246) q[3];
sx q[3];
rz(-2.0913561) q[3];
sx q[3];
rz(-0.27568451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7212873) q[0];
sx q[0];
rz(-0.5834226) q[0];
sx q[0];
rz(-0.48459184) q[0];
rz(-1.2424319) q[1];
sx q[1];
rz(-2.395605) q[1];
sx q[1];
rz(2.7925083) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9681681) q[0];
sx q[0];
rz(-1.613488) q[0];
sx q[0];
rz(-2.1195565) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.227581) q[2];
sx q[2];
rz(-1.6417393) q[2];
sx q[2];
rz(-0.22954839) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.0047093948) q[1];
sx q[1];
rz(-1.7994969) q[1];
sx q[1];
rz(1.7654224) q[1];
rz(0.99557568) q[3];
sx q[3];
rz(-2.4565426) q[3];
sx q[3];
rz(-3.0127061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7829973) q[2];
sx q[2];
rz(-1.9390257) q[2];
sx q[2];
rz(-2.6915468) q[2];
rz(2.3027244) q[3];
sx q[3];
rz(-1.5939555) q[3];
sx q[3];
rz(0.1990327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6435476) q[0];
sx q[0];
rz(-1.6723375) q[0];
sx q[0];
rz(-0.099844649) q[0];
rz(1.3205344) q[1];
sx q[1];
rz(-0.99597275) q[1];
sx q[1];
rz(2.5340396) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6927101) q[0];
sx q[0];
rz(-1.6537871) q[0];
sx q[0];
rz(-1.4957499) q[0];
rz(-pi) q[1];
x q[1];
rz(1.728292) q[2];
sx q[2];
rz(-2.11065) q[2];
sx q[2];
rz(1.4992204) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1888426) q[1];
sx q[1];
rz(-0.61988587) q[1];
sx q[1];
rz(-0.37921885) q[1];
rz(0.097955161) q[3];
sx q[3];
rz(-1.8281015) q[3];
sx q[3];
rz(-0.78661455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1943835) q[2];
sx q[2];
rz(-0.76346976) q[2];
sx q[2];
rz(2.6109931) q[2];
rz(-2.7001906) q[3];
sx q[3];
rz(-2.1680021) q[3];
sx q[3];
rz(2.1921564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7049103) q[0];
sx q[0];
rz(-1.3446151) q[0];
sx q[0];
rz(0.5740903) q[0];
rz(-0.87999815) q[1];
sx q[1];
rz(-2.0206385) q[1];
sx q[1];
rz(-1.3425739) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3898252) q[0];
sx q[0];
rz(-0.51503599) q[0];
sx q[0];
rz(-1.1549321) q[0];
x q[1];
rz(-1.1583352) q[2];
sx q[2];
rz(-1.165373) q[2];
sx q[2];
rz(-0.5328005) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9989661) q[1];
sx q[1];
rz(-1.9311096) q[1];
sx q[1];
rz(2.8887038) q[1];
rz(-pi) q[2];
x q[2];
rz(0.17832605) q[3];
sx q[3];
rz(-2.470782) q[3];
sx q[3];
rz(1.7899771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4673956) q[2];
sx q[2];
rz(-0.28417045) q[2];
sx q[2];
rz(0.44710844) q[2];
rz(-2.1111264) q[3];
sx q[3];
rz(-1.5117398) q[3];
sx q[3];
rz(2.7569421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36987385) q[0];
sx q[0];
rz(-1.2971224) q[0];
sx q[0];
rz(2.7287927) q[0];
rz(2.0665456) q[1];
sx q[1];
rz(-2.0037035) q[1];
sx q[1];
rz(0.80361754) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34016383) q[0];
sx q[0];
rz(-2.0716084) q[0];
sx q[0];
rz(-0.97789557) q[0];
rz(0.86974025) q[2];
sx q[2];
rz(-0.77574965) q[2];
sx q[2];
rz(2.9091331) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4839125) q[1];
sx q[1];
rz(-1.5849304) q[1];
sx q[1];
rz(2.3573228) q[1];
x q[2];
rz(-2.7429306) q[3];
sx q[3];
rz(-2.9473262) q[3];
sx q[3];
rz(2.8630321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2794118) q[2];
sx q[2];
rz(-1.923424) q[2];
sx q[2];
rz(-1.6406406) q[2];
rz(0.62266478) q[3];
sx q[3];
rz(-2.3707135) q[3];
sx q[3];
rz(-1.6525432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.415446) q[0];
sx q[0];
rz(-1.6212689) q[0];
sx q[0];
rz(-2.5836482) q[0];
rz(2.5803512) q[1];
sx q[1];
rz(-1.3713501) q[1];
sx q[1];
rz(1.4161313) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2673906) q[0];
sx q[0];
rz(-0.54055291) q[0];
sx q[0];
rz(-2.1442599) q[0];
x q[1];
rz(2.2165856) q[2];
sx q[2];
rz(-2.6753798) q[2];
sx q[2];
rz(-3.1104308) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0812644) q[1];
sx q[1];
rz(-2.7980479) q[1];
sx q[1];
rz(-1.8650123) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.21702311) q[3];
sx q[3];
rz(-0.76730928) q[3];
sx q[3];
rz(-1.8930774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.8884362) q[2];
sx q[2];
rz(-0.77740589) q[2];
sx q[2];
rz(-1.0003132) q[2];
rz(0.12217626) q[3];
sx q[3];
rz(-1.5001985) q[3];
sx q[3];
rz(2.9848671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.4633453) q[0];
sx q[0];
rz(-1.954701) q[0];
sx q[0];
rz(3.1372702) q[0];
rz(2.9777572) q[1];
sx q[1];
rz(-2.7163353) q[1];
sx q[1];
rz(1.3660627) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9095847) q[0];
sx q[0];
rz(-0.9238832) q[0];
sx q[0];
rz(-0.034897403) q[0];
rz(-pi) q[1];
rz(-2.4539656) q[2];
sx q[2];
rz(-0.4677217) q[2];
sx q[2];
rz(-2.2224838) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1613579) q[1];
sx q[1];
rz(-1.5346077) q[1];
sx q[1];
rz(2.8904032) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2774599) q[3];
sx q[3];
rz(-2.257405) q[3];
sx q[3];
rz(0.1089801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0020478) q[2];
sx q[2];
rz(-0.9674955) q[2];
sx q[2];
rz(2.5223993) q[2];
rz(0.55245095) q[3];
sx q[3];
rz(-1.8360454) q[3];
sx q[3];
rz(-2.9964871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7904952) q[0];
sx q[0];
rz(-0.18432291) q[0];
sx q[0];
rz(2.6628394) q[0];
rz(-1.152773) q[1];
sx q[1];
rz(-0.29007998) q[1];
sx q[1];
rz(0.94720381) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21175948) q[0];
sx q[0];
rz(-2.9538028) q[0];
sx q[0];
rz(1.387557) q[0];
rz(-1.2113153) q[2];
sx q[2];
rz(-2.5229342) q[2];
sx q[2];
rz(-1.930069) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.80120197) q[1];
sx q[1];
rz(-1.2653192) q[1];
sx q[1];
rz(0.82251541) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8197038) q[3];
sx q[3];
rz(-1.9134023) q[3];
sx q[3];
rz(3.1150888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7572215) q[2];
sx q[2];
rz(-0.20558509) q[2];
sx q[2];
rz(-2.8749386) q[2];
rz(-2.0742778) q[3];
sx q[3];
rz(-1.9284733) q[3];
sx q[3];
rz(2.1783569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.08854475) q[0];
sx q[0];
rz(-1.6041258) q[0];
sx q[0];
rz(-1.5969101) q[0];
rz(2.0883941) q[1];
sx q[1];
rz(-1.2146626) q[1];
sx q[1];
rz(0.76520898) q[1];
rz(1.5736035) q[2];
sx q[2];
rz(-1.8127828) q[2];
sx q[2];
rz(0.91290963) q[2];
rz(-1.7942747) q[3];
sx q[3];
rz(-2.1383994) q[3];
sx q[3];
rz(2.4169328) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
