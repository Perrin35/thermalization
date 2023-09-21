OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.59453073) q[0];
sx q[0];
rz(-1.1214331) q[0];
sx q[0];
rz(-2.9601331) q[0];
rz(-1.0815066) q[1];
sx q[1];
rz(-2.4681611) q[1];
sx q[1];
rz(-2.0884617) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5232789) q[0];
sx q[0];
rz(-1.7821454) q[0];
sx q[0];
rz(-2.0413415) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7396169) q[2];
sx q[2];
rz(-0.66265124) q[2];
sx q[2];
rz(1.0647917) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3095113) q[1];
sx q[1];
rz(-0.81334269) q[1];
sx q[1];
rz(-1.3002212) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9795322) q[3];
sx q[3];
rz(-1.9063623) q[3];
sx q[3];
rz(-3.0367362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.16333214) q[2];
sx q[2];
rz(-2.0330727) q[2];
sx q[2];
rz(1.7738316) q[2];
rz(2.1286428) q[3];
sx q[3];
rz(-2.2949341) q[3];
sx q[3];
rz(3.0701385) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0435836) q[0];
sx q[0];
rz(-1.6921035) q[0];
sx q[0];
rz(2.5464771) q[0];
rz(-2.0960506) q[1];
sx q[1];
rz(-1.4141934) q[1];
sx q[1];
rz(1.5140623) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13578116) q[0];
sx q[0];
rz(-1.0732871) q[0];
sx q[0];
rz(-0.46732975) q[0];
x q[1];
rz(-2.3081231) q[2];
sx q[2];
rz(-1.5499299) q[2];
sx q[2];
rz(1.2534864) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4501614) q[1];
sx q[1];
rz(-1.3988004) q[1];
sx q[1];
rz(-1.7683692) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4684832) q[3];
sx q[3];
rz(-1.521109) q[3];
sx q[3];
rz(1.4440086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4864768) q[2];
sx q[2];
rz(-2.1003508) q[2];
sx q[2];
rz(-0.79616037) q[2];
rz(0.97186175) q[3];
sx q[3];
rz(-0.704851) q[3];
sx q[3];
rz(0.17175737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7650836) q[0];
sx q[0];
rz(-2.3628545) q[0];
sx q[0];
rz(-0.064095108) q[0];
rz(-0.31072101) q[1];
sx q[1];
rz(-1.4704082) q[1];
sx q[1];
rz(1.6832738) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9653939) q[0];
sx q[0];
rz(-1.5237234) q[0];
sx q[0];
rz(1.4233627) q[0];
x q[1];
rz(-2.4887423) q[2];
sx q[2];
rz(-2.0712426) q[2];
sx q[2];
rz(2.1476538) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.26865667) q[1];
sx q[1];
rz(-1.6246512) q[1];
sx q[1];
rz(2.7320646) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0890373) q[3];
sx q[3];
rz(-2.8263546) q[3];
sx q[3];
rz(-1.0759575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1594499) q[2];
sx q[2];
rz(-2.2950256) q[2];
sx q[2];
rz(0.67908755) q[2];
rz(1.7689765) q[3];
sx q[3];
rz(-1.2981828) q[3];
sx q[3];
rz(-0.94846559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2963592) q[0];
sx q[0];
rz(-2.5722752) q[0];
sx q[0];
rz(1.1244208) q[0];
rz(1.9212978) q[1];
sx q[1];
rz(-1.9492457) q[1];
sx q[1];
rz(1.4845928) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2161908) q[0];
sx q[0];
rz(-1.6214217) q[0];
sx q[0];
rz(0.030406818) q[0];
rz(-pi) q[1];
rz(-1.053327) q[2];
sx q[2];
rz(-1.7382858) q[2];
sx q[2];
rz(2.6094764) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.9323737) q[1];
sx q[1];
rz(-2.0062431) q[1];
sx q[1];
rz(0.93851628) q[1];
rz(2.4184166) q[3];
sx q[3];
rz(-1.3190862) q[3];
sx q[3];
rz(0.46079208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0174039) q[2];
sx q[2];
rz(-1.8576531) q[2];
sx q[2];
rz(-2.1253288) q[2];
rz(-1.7381564) q[3];
sx q[3];
rz(-1.6442464) q[3];
sx q[3];
rz(-1.0884292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50773412) q[0];
sx q[0];
rz(-1.8554747) q[0];
sx q[0];
rz(2.741709) q[0];
rz(-1.1625066) q[1];
sx q[1];
rz(-1.8116654) q[1];
sx q[1];
rz(-0.17366017) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4325718) q[0];
sx q[0];
rz(-0.077843277) q[0];
sx q[0];
rz(-2.7709333) q[0];
rz(-pi) q[1];
rz(-0.47541754) q[2];
sx q[2];
rz(-0.43736514) q[2];
sx q[2];
rz(-1.4471444) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.793321) q[1];
sx q[1];
rz(-1.4655359) q[1];
sx q[1];
rz(3.1131016) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3624304) q[3];
sx q[3];
rz(-1.2456129) q[3];
sx q[3];
rz(-1.9735826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6614723) q[2];
sx q[2];
rz(-1.9813333) q[2];
sx q[2];
rz(-2.373467) q[2];
rz(2.2875732) q[3];
sx q[3];
rz(-1.7212399) q[3];
sx q[3];
rz(-0.97222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0356692) q[0];
sx q[0];
rz(-0.24374715) q[0];
sx q[0];
rz(1.6712028) q[0];
rz(-2.629783) q[1];
sx q[1];
rz(-2.6302331) q[1];
sx q[1];
rz(-1.9981729) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2238732) q[0];
sx q[0];
rz(-1.6544764) q[0];
sx q[0];
rz(2.9907945) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60322275) q[2];
sx q[2];
rz(-1.2957186) q[2];
sx q[2];
rz(-3.0657257) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.032635078) q[1];
sx q[1];
rz(-0.40954486) q[1];
sx q[1];
rz(1.1213379) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0877803) q[3];
sx q[3];
rz(-0.59637585) q[3];
sx q[3];
rz(-3.0697825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4986971) q[2];
sx q[2];
rz(-2.2968447) q[2];
sx q[2];
rz(-1.6112304) q[2];
rz(1.6879843) q[3];
sx q[3];
rz(-1.0691103) q[3];
sx q[3];
rz(2.8924275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-3.0221508) q[0];
sx q[0];
rz(-2.3972153) q[0];
sx q[0];
rz(-1.7013593) q[0];
rz(-0.72921905) q[1];
sx q[1];
rz(-1.1518642) q[1];
sx q[1];
rz(2.008332) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2204809) q[0];
sx q[0];
rz(-0.32520121) q[0];
sx q[0];
rz(3.0082506) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.015958162) q[2];
sx q[2];
rz(-0.71231132) q[2];
sx q[2];
rz(-2.6238837) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.28753528) q[1];
sx q[1];
rz(-0.71811986) q[1];
sx q[1];
rz(-1.1487886) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4134737) q[3];
sx q[3];
rz(-1.3120578) q[3];
sx q[3];
rz(1.3820005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7363654) q[2];
sx q[2];
rz(-2.0183125) q[2];
sx q[2];
rz(-2.6531632) q[2];
rz(-1.8296261) q[3];
sx q[3];
rz(-0.043881504) q[3];
sx q[3];
rz(-0.64129889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0768123) q[0];
sx q[0];
rz(-1.8490054) q[0];
sx q[0];
rz(3.0704165) q[0];
rz(3.1094303) q[1];
sx q[1];
rz(-1.3379438) q[1];
sx q[1];
rz(-1.2088998) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6541518) q[0];
sx q[0];
rz(-1.7983266) q[0];
sx q[0];
rz(2.3339416) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2279943) q[2];
sx q[2];
rz(-0.29619869) q[2];
sx q[2];
rz(1.9343455) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.80450088) q[1];
sx q[1];
rz(-0.13131222) q[1];
sx q[1];
rz(1.4485703) q[1];
rz(-pi) q[2];
x q[2];
rz(0.62187059) q[3];
sx q[3];
rz(-1.8727881) q[3];
sx q[3];
rz(-2.5627476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.20748392) q[2];
sx q[2];
rz(-0.19342962) q[2];
sx q[2];
rz(1.5709546) q[2];
rz(-0.87336826) q[3];
sx q[3];
rz(-1.4040754) q[3];
sx q[3];
rz(2.7895555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(2.1677925) q[0];
sx q[0];
rz(-1.6163102) q[0];
sx q[0];
rz(2.8299676) q[0];
rz(-0.82178003) q[1];
sx q[1];
rz(-0.59097925) q[1];
sx q[1];
rz(1.5100381) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8845997) q[0];
sx q[0];
rz(-1.25367) q[0];
sx q[0];
rz(-1.7811799) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0636343) q[2];
sx q[2];
rz(-2.7535451) q[2];
sx q[2];
rz(2.0575112) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1875302) q[1];
sx q[1];
rz(-0.64779753) q[1];
sx q[1];
rz(-2.8026583) q[1];
x q[2];
rz(-1.7543091) q[3];
sx q[3];
rz(-1.3418875) q[3];
sx q[3];
rz(0.51975091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2733549) q[2];
sx q[2];
rz(-2.6028825) q[2];
sx q[2];
rz(2.3256425) q[2];
rz(0.50968918) q[3];
sx q[3];
rz(-2.3563801) q[3];
sx q[3];
rz(1.2379439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8907392) q[0];
sx q[0];
rz(-1.4287404) q[0];
sx q[0];
rz(2.9113286) q[0];
rz(2.5157805) q[1];
sx q[1];
rz(-2.1964549) q[1];
sx q[1];
rz(0.65840107) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9394768) q[0];
sx q[0];
rz(-2.3411223) q[0];
sx q[0];
rz(-0.70864422) q[0];
rz(-2.7294331) q[2];
sx q[2];
rz(-2.4240652) q[2];
sx q[2];
rz(1.621643) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2645996) q[1];
sx q[1];
rz(-1.0942642) q[1];
sx q[1];
rz(-2.9546253) q[1];
rz(-0.55024054) q[3];
sx q[3];
rz(-2.5329258) q[3];
sx q[3];
rz(-0.88702162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4663503) q[2];
sx q[2];
rz(-0.8224951) q[2];
sx q[2];
rz(-0.33774439) q[2];
rz(-1.0234458) q[3];
sx q[3];
rz(-1.7815536) q[3];
sx q[3];
rz(1.9256928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52453775) q[0];
sx q[0];
rz(-0.90159566) q[0];
sx q[0];
rz(3.0774075) q[0];
rz(-1.9032003) q[1];
sx q[1];
rz(-1.7108142) q[1];
sx q[1];
rz(-1.6731813) q[1];
rz(1.3255618) q[2];
sx q[2];
rz(-0.64074466) q[2];
sx q[2];
rz(-1.7111039) q[2];
rz(0.099048793) q[3];
sx q[3];
rz(-0.5746114) q[3];
sx q[3];
rz(-0.050669908) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];