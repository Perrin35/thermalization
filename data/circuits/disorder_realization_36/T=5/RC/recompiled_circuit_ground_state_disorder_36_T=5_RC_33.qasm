OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.9517684) q[0];
sx q[0];
rz(3.8444001) q[0];
sx q[0];
rz(8.202717) q[0];
rz(3.089978) q[1];
sx q[1];
rz(-1.6366704) q[1];
sx q[1];
rz(1.6610891) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9771271) q[0];
sx q[0];
rz(-1.7344622) q[0];
sx q[0];
rz(2.2360691) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0323276) q[2];
sx q[2];
rz(-1.4022298) q[2];
sx q[2];
rz(1.3079461) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8424219) q[1];
sx q[1];
rz(-1.7452733) q[1];
sx q[1];
rz(2.7122696) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0546255) q[3];
sx q[3];
rz(-1.9619521) q[3];
sx q[3];
rz(-3.1227171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1808971) q[2];
sx q[2];
rz(-2.3263859) q[2];
sx q[2];
rz(-1.0308456) q[2];
rz(-0.24813063) q[3];
sx q[3];
rz(-1.7957567) q[3];
sx q[3];
rz(0.26721755) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2613075) q[0];
sx q[0];
rz(-1.945865) q[0];
sx q[0];
rz(-1.1241359) q[0];
rz(1.0719489) q[1];
sx q[1];
rz(-1.5644466) q[1];
sx q[1];
rz(-0.42676485) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65157344) q[0];
sx q[0];
rz(-2.314004) q[0];
sx q[0];
rz(3.1378967) q[0];
rz(-pi) q[1];
rz(1.3091121) q[2];
sx q[2];
rz(-1.905447) q[2];
sx q[2];
rz(1.7098914) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2093231) q[1];
sx q[1];
rz(-0.70581064) q[1];
sx q[1];
rz(-2.9747444) q[1];
rz(-pi) q[2];
rz(-1.1817478) q[3];
sx q[3];
rz(-1.7579683) q[3];
sx q[3];
rz(-1.0295008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.65853226) q[2];
sx q[2];
rz(-2.6250562) q[2];
sx q[2];
rz(-3.0413682) q[2];
rz(2.0641067) q[3];
sx q[3];
rz(-1.5651549) q[3];
sx q[3];
rz(-1.235435) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.539262) q[0];
sx q[0];
rz(-2.1644008) q[0];
sx q[0];
rz(-1.4343028) q[0];
rz(2.3151248) q[1];
sx q[1];
rz(-1.6441328) q[1];
sx q[1];
rz(0.43063393) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42303681) q[0];
sx q[0];
rz(-1.406975) q[0];
sx q[0];
rz(2.8680766) q[0];
rz(-3.011322) q[2];
sx q[2];
rz(-1.300087) q[2];
sx q[2];
rz(-2.475428) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0890177) q[1];
sx q[1];
rz(-0.7886501) q[1];
sx q[1];
rz(2.1441205) q[1];
x q[2];
rz(0.39726302) q[3];
sx q[3];
rz(-2.0797585) q[3];
sx q[3];
rz(2.404117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2677801) q[2];
sx q[2];
rz(-1.9766821) q[2];
sx q[2];
rz(-0.57360345) q[2];
rz(2.7576647) q[3];
sx q[3];
rz(-1.3150747) q[3];
sx q[3];
rz(0.65281502) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4189932) q[0];
sx q[0];
rz(-0.79177952) q[0];
sx q[0];
rz(-2.0844039) q[0];
rz(0.81622299) q[1];
sx q[1];
rz(-2.1482601) q[1];
sx q[1];
rz(2.9073471) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1907983) q[0];
sx q[0];
rz(-1.6802579) q[0];
sx q[0];
rz(-1.1534766) q[0];
rz(-2.7170638) q[2];
sx q[2];
rz(-2.0574951) q[2];
sx q[2];
rz(2.8667712) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1409451) q[1];
sx q[1];
rz(-1.06303) q[1];
sx q[1];
rz(0.80054466) q[1];
rz(1.0660421) q[3];
sx q[3];
rz(-1.1467993) q[3];
sx q[3];
rz(-3.0847065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8546042) q[2];
sx q[2];
rz(-1.8082666) q[2];
sx q[2];
rz(-2.2387779) q[2];
rz(-1.758894) q[3];
sx q[3];
rz(-0.32302502) q[3];
sx q[3];
rz(-2.2949016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46120241) q[0];
sx q[0];
rz(-1.8620551) q[0];
sx q[0];
rz(-0.70988208) q[0];
rz(-0.4711802) q[1];
sx q[1];
rz(-1.9147562) q[1];
sx q[1];
rz(-3.0772298) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36633766) q[0];
sx q[0];
rz(-1.7058055) q[0];
sx q[0];
rz(-2.058821) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.2336524) q[2];
sx q[2];
rz(-2.0548247) q[2];
sx q[2];
rz(2.7357227) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.91436323) q[1];
sx q[1];
rz(-1.0681947) q[1];
sx q[1];
rz(-1.0881626) q[1];
rz(-pi) q[2];
x q[2];
rz(0.10897763) q[3];
sx q[3];
rz(-0.087317467) q[3];
sx q[3];
rz(2.2918037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.71089661) q[2];
sx q[2];
rz(-2.5888093) q[2];
sx q[2];
rz(-0.95332471) q[2];
rz(-2.583875) q[3];
sx q[3];
rz(-1.4671114) q[3];
sx q[3];
rz(-2.6127156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3144658) q[0];
sx q[0];
rz(-1.0449266) q[0];
sx q[0];
rz(-2.1449828) q[0];
rz(-3.0820471) q[1];
sx q[1];
rz(-0.98809067) q[1];
sx q[1];
rz(1.3522118) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9837275) q[0];
sx q[0];
rz(-1.0520555) q[0];
sx q[0];
rz(2.5602976) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.77162051) q[2];
sx q[2];
rz(-1.3305656) q[2];
sx q[2];
rz(0.26773237) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.88207605) q[1];
sx q[1];
rz(-1.8456371) q[1];
sx q[1];
rz(-0.3104574) q[1];
x q[2];
rz(2.2289946) q[3];
sx q[3];
rz(-1.8479947) q[3];
sx q[3];
rz(2.4656221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.754564) q[2];
sx q[2];
rz(-1.2754385) q[2];
sx q[2];
rz(-2.582029) q[2];
rz(-0.83545056) q[3];
sx q[3];
rz(-2.7159034) q[3];
sx q[3];
rz(-1.8387509) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4556722) q[0];
sx q[0];
rz(-2.4113825) q[0];
sx q[0];
rz(-0.36668229) q[0];
rz(-2.2600251) q[1];
sx q[1];
rz(-2.5036948) q[1];
sx q[1];
rz(-1.7104023) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0013051) q[0];
sx q[0];
rz(-1.9594345) q[0];
sx q[0];
rz(-2.4508618) q[0];
rz(-1.9164247) q[2];
sx q[2];
rz(-2.2273438) q[2];
sx q[2];
rz(-2.569927) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.28364946) q[1];
sx q[1];
rz(-1.6873296) q[1];
sx q[1];
rz(1.5999419) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.335698) q[3];
sx q[3];
rz(-2.2155016) q[3];
sx q[3];
rz(2.7489782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.35393474) q[2];
sx q[2];
rz(-2.6255609) q[2];
sx q[2];
rz(0.77989522) q[2];
rz(2.6025313) q[3];
sx q[3];
rz(-1.6293679) q[3];
sx q[3];
rz(-0.61162925) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1270545) q[0];
sx q[0];
rz(-2.4712565) q[0];
sx q[0];
rz(-2.3624453) q[0];
rz(-2.6226793) q[1];
sx q[1];
rz(-1.2823558) q[1];
sx q[1];
rz(-0.40294495) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.084672734) q[0];
sx q[0];
rz(-2.2209327) q[0];
sx q[0];
rz(2.5385602) q[0];
x q[1];
rz(2.1261931) q[2];
sx q[2];
rz(-0.55335535) q[2];
sx q[2];
rz(0.28552688) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1225452) q[1];
sx q[1];
rz(-1.9662939) q[1];
sx q[1];
rz(2.4707153) q[1];
x q[2];
rz(-1.1717848) q[3];
sx q[3];
rz(-0.61099377) q[3];
sx q[3];
rz(2.4965167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0275823) q[2];
sx q[2];
rz(-1.1677914) q[2];
sx q[2];
rz(-1.1172509) q[2];
rz(-0.69636017) q[3];
sx q[3];
rz(-1.1098692) q[3];
sx q[3];
rz(1.2874359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(0.47525147) q[0];
sx q[0];
rz(-0.23463686) q[0];
sx q[0];
rz(-0.41467211) q[0];
rz(2.0798202) q[1];
sx q[1];
rz(-2.2551408) q[1];
sx q[1];
rz(0.83962238) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79796004) q[0];
sx q[0];
rz(-2.0164444) q[0];
sx q[0];
rz(0.4526505) q[0];
rz(-pi) q[1];
rz(2.5171385) q[2];
sx q[2];
rz(-2.3169998) q[2];
sx q[2];
rz(-0.39442447) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.190358) q[1];
sx q[1];
rz(-0.91427416) q[1];
sx q[1];
rz(-1.7503529) q[1];
x q[2];
rz(-1.8168338) q[3];
sx q[3];
rz(-1.2241505) q[3];
sx q[3];
rz(0.7905761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.35905579) q[2];
sx q[2];
rz(-1.3691207) q[2];
sx q[2];
rz(3.0699406) q[2];
rz(0.045529384) q[3];
sx q[3];
rz(-1.1979016) q[3];
sx q[3];
rz(-0.45904747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.5422106) q[0];
sx q[0];
rz(-2.5686503) q[0];
sx q[0];
rz(2.086916) q[0];
rz(0.032546267) q[1];
sx q[1];
rz(-2.1701505) q[1];
sx q[1];
rz(2.6194825) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.78799) q[0];
sx q[0];
rz(-2.4756469) q[0];
sx q[0];
rz(-0.5086201) q[0];
rz(-2.4049875) q[2];
sx q[2];
rz(-1.1140559) q[2];
sx q[2];
rz(-1.119594) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5473816) q[1];
sx q[1];
rz(-1.949233) q[1];
sx q[1];
rz(-2.7345279) q[1];
rz(-1.2100092) q[3];
sx q[3];
rz(-0.57769247) q[3];
sx q[3];
rz(-3.011774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1471499) q[2];
sx q[2];
rz(-2.9394737) q[2];
sx q[2];
rz(1.408255) q[2];
rz(-1.5857006) q[3];
sx q[3];
rz(-2.9097911) q[3];
sx q[3];
rz(2.2304992) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0178575) q[0];
sx q[0];
rz(-1.3029079) q[0];
sx q[0];
rz(-2.2444176) q[0];
rz(2.1625715) q[1];
sx q[1];
rz(-1.9985825) q[1];
sx q[1];
rz(-0.40252007) q[1];
rz(2.0918905) q[2];
sx q[2];
rz(-2.7626531) q[2];
sx q[2];
rz(0.74041453) q[2];
rz(2.4784563) q[3];
sx q[3];
rz(-1.4226154) q[3];
sx q[3];
rz(-2.0285574) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
