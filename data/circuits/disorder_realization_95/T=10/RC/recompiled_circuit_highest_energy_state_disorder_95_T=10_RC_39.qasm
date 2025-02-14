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
rz(0.79701841) q[0];
sx q[0];
rz(-2.2198644) q[0];
sx q[0];
rz(1.2603238) q[0];
rz(-2.4951275) q[1];
sx q[1];
rz(-0.87352455) q[1];
sx q[1];
rz(0.33198196) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2352358) q[0];
sx q[0];
rz(-0.42354326) q[0];
sx q[0];
rz(1.5111501) q[0];
x q[1];
rz(-1.9464363) q[2];
sx q[2];
rz(-0.90919288) q[2];
sx q[2];
rz(2.313569) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1412072) q[1];
sx q[1];
rz(-1.5391333) q[1];
sx q[1];
rz(-2.382676) q[1];
rz(-3.0044972) q[3];
sx q[3];
rz(-2.6957316) q[3];
sx q[3];
rz(0.23030989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6648286) q[2];
sx q[2];
rz(-2.2085184) q[2];
sx q[2];
rz(-0.19922166) q[2];
rz(-0.18579379) q[3];
sx q[3];
rz(-1.7265065) q[3];
sx q[3];
rz(-1.7004405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81714565) q[0];
sx q[0];
rz(-2.6528093) q[0];
sx q[0];
rz(2.4031438) q[0];
rz(-1.4456519) q[1];
sx q[1];
rz(-1.7414469) q[1];
sx q[1];
rz(2.6587528) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49119887) q[0];
sx q[0];
rz(-0.20696124) q[0];
sx q[0];
rz(-0.82024337) q[0];
rz(-pi) q[1];
rz(-1.2224169) q[2];
sx q[2];
rz(-1.7357706) q[2];
sx q[2];
rz(2.4553187) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.184471) q[1];
sx q[1];
rz(-1.0252684) q[1];
sx q[1];
rz(2.1197182) q[1];
rz(-2.1294503) q[3];
sx q[3];
rz(-1.4317272) q[3];
sx q[3];
rz(2.6731051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.26756755) q[2];
sx q[2];
rz(-1.4635307) q[2];
sx q[2];
rz(-1.2919424) q[2];
rz(2.0180295) q[3];
sx q[3];
rz(-2.4396887) q[3];
sx q[3];
rz(-1.4152214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32347754) q[0];
sx q[0];
rz(-2.1577142) q[0];
sx q[0];
rz(1.6997319) q[0];
rz(-1.9135176) q[1];
sx q[1];
rz(-1.2308729) q[1];
sx q[1];
rz(1.0736116) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31698727) q[0];
sx q[0];
rz(-1.2277216) q[0];
sx q[0];
rz(2.4113684) q[0];
rz(0.19811689) q[2];
sx q[2];
rz(-0.92495944) q[2];
sx q[2];
rz(-0.4465296) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5646208) q[1];
sx q[1];
rz(-2.4384017) q[1];
sx q[1];
rz(1.6220051) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4074247) q[3];
sx q[3];
rz(-2.022812) q[3];
sx q[3];
rz(-2.3491835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8492665) q[2];
sx q[2];
rz(-2.8385415) q[2];
sx q[2];
rz(1.0230052) q[2];
rz(0.060221378) q[3];
sx q[3];
rz(-1.1465719) q[3];
sx q[3];
rz(-2.4165418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.4802454) q[0];
sx q[0];
rz(-0.1460954) q[0];
sx q[0];
rz(0.51937854) q[0];
rz(2.5410779) q[1];
sx q[1];
rz(-1.0888211) q[1];
sx q[1];
rz(-0.75469887) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1491283) q[0];
sx q[0];
rz(-2.0880648) q[0];
sx q[0];
rz(0.37578342) q[0];
rz(-pi) q[1];
x q[1];
rz(0.86769957) q[2];
sx q[2];
rz(-1.6804196) q[2];
sx q[2];
rz(0.90849691) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2024196) q[1];
sx q[1];
rz(-1.8853097) q[1];
sx q[1];
rz(-2.0862723) q[1];
rz(0.56086991) q[3];
sx q[3];
rz(-1.2301233) q[3];
sx q[3];
rz(-2.7531227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.34919136) q[2];
sx q[2];
rz(-0.79734048) q[2];
sx q[2];
rz(-2.965773) q[2];
rz(2.0154121) q[3];
sx q[3];
rz(-1.3577941) q[3];
sx q[3];
rz(-3.0774097) q[3];
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
rz(-pi/2) q[3];
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
rz(-0.46614161) q[0];
sx q[0];
rz(-1.6459246) q[0];
sx q[0];
rz(-0.58188907) q[0];
rz(-1.9582845) q[1];
sx q[1];
rz(-2.4300523) q[1];
sx q[1];
rz(-1.2540832) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39685707) q[0];
sx q[0];
rz(-1.2484457) q[0];
sx q[0];
rz(-2.2765323) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6627867) q[2];
sx q[2];
rz(-1.1421912) q[2];
sx q[2];
rz(1.3474825) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5820133) q[1];
sx q[1];
rz(-1.2767643) q[1];
sx q[1];
rz(-0.21802417) q[1];
rz(-pi) q[2];
rz(-2.6370089) q[3];
sx q[3];
rz(-1.736803) q[3];
sx q[3];
rz(0.61246757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9767655) q[2];
sx q[2];
rz(-0.81563121) q[2];
sx q[2];
rz(-2.7867479) q[2];
rz(-1.1497078) q[3];
sx q[3];
rz(-1.0747654) q[3];
sx q[3];
rz(2.3000075) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0717764) q[0];
sx q[0];
rz(-0.27158296) q[0];
sx q[0];
rz(-0.040104453) q[0];
rz(-1.6768203) q[1];
sx q[1];
rz(-2.1720839) q[1];
sx q[1];
rz(-0.47223314) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.036333648) q[0];
sx q[0];
rz(-1.3091631) q[0];
sx q[0];
rz(-0.18966578) q[0];
rz(-pi) q[1];
rz(-1.5664928) q[2];
sx q[2];
rz(-0.65281463) q[2];
sx q[2];
rz(0.091856591) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.66989723) q[1];
sx q[1];
rz(-1.2409087) q[1];
sx q[1];
rz(-2.8389588) q[1];
rz(-0.29904265) q[3];
sx q[3];
rz(-0.78419331) q[3];
sx q[3];
rz(-1.2687781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.27452305) q[2];
sx q[2];
rz(-2.3752866) q[2];
sx q[2];
rz(1.6609894) q[2];
rz(-0.04145043) q[3];
sx q[3];
rz(-0.71235123) q[3];
sx q[3];
rz(2.1672772) q[3];
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
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90254766) q[0];
sx q[0];
rz(-2.3426265) q[0];
sx q[0];
rz(-2.3833185) q[0];
rz(2.7039418) q[1];
sx q[1];
rz(-1.2145019) q[1];
sx q[1];
rz(2.1877916) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8945635) q[0];
sx q[0];
rz(-1.6163384) q[0];
sx q[0];
rz(1.1883424) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.96289159) q[2];
sx q[2];
rz(-0.94568077) q[2];
sx q[2];
rz(1.4639548) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5686534) q[1];
sx q[1];
rz(-1.6332989) q[1];
sx q[1];
rz(2.2745132) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1498044) q[3];
sx q[3];
rz(-2.6564993) q[3];
sx q[3];
rz(-2.4320784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7268251) q[2];
sx q[2];
rz(-0.32764062) q[2];
sx q[2];
rz(-0.3978351) q[2];
rz(-1.8757403) q[3];
sx q[3];
rz(-1.7222907) q[3];
sx q[3];
rz(0.46441594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1606814) q[0];
sx q[0];
rz(-2.2137764) q[0];
sx q[0];
rz(0.30447793) q[0];
rz(-1.132384) q[1];
sx q[1];
rz(-2.4696923) q[1];
sx q[1];
rz(-1.0736046) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7454483) q[0];
sx q[0];
rz(-1.5477596) q[0];
sx q[0];
rz(3.1246326) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9307641) q[2];
sx q[2];
rz(-1.9641293) q[2];
sx q[2];
rz(-2.8434812) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1731092) q[1];
sx q[1];
rz(-2.1595528) q[1];
sx q[1];
rz(-0.29300479) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.066929265) q[3];
sx q[3];
rz(-0.64475497) q[3];
sx q[3];
rz(0.88551846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.73211804) q[2];
sx q[2];
rz(-1.858859) q[2];
sx q[2];
rz(-3.0916302) q[2];
rz(-2.2299855) q[3];
sx q[3];
rz(-0.92602366) q[3];
sx q[3];
rz(-0.91134206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53495812) q[0];
sx q[0];
rz(-1.7462523) q[0];
sx q[0];
rz(0.362679) q[0];
rz(-2.4210988) q[1];
sx q[1];
rz(-0.81169218) q[1];
sx q[1];
rz(-1.1788751) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6828109) q[0];
sx q[0];
rz(-1.0014373) q[0];
sx q[0];
rz(-0.95627934) q[0];
rz(0.32068738) q[2];
sx q[2];
rz(-1.4960519) q[2];
sx q[2];
rz(0.2439258) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8803823) q[1];
sx q[1];
rz(-2.9453813) q[1];
sx q[1];
rz(1.0537149) q[1];
x q[2];
rz(-1.1580771) q[3];
sx q[3];
rz(-1.0636119) q[3];
sx q[3];
rz(-0.87178236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1371586) q[2];
sx q[2];
rz(-0.68778554) q[2];
sx q[2];
rz(1.9384109) q[2];
rz(3.1249937) q[3];
sx q[3];
rz(-1.4998452) q[3];
sx q[3];
rz(-1.2703007) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59852973) q[0];
sx q[0];
rz(-2.6180584) q[0];
sx q[0];
rz(1.4434927) q[0];
rz(-0.47830018) q[1];
sx q[1];
rz(-0.6820448) q[1];
sx q[1];
rz(-1.2102478) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.011199) q[0];
sx q[0];
rz(-2.3558295) q[0];
sx q[0];
rz(-2.4146621) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4912676) q[2];
sx q[2];
rz(-1.7129833) q[2];
sx q[2];
rz(1.1917758) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5961884) q[1];
sx q[1];
rz(-2.8126723) q[1];
sx q[1];
rz(0.1620345) q[1];
x q[2];
rz(-1.7948512) q[3];
sx q[3];
rz(-1.2478078) q[3];
sx q[3];
rz(-2.546348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.18834194) q[2];
sx q[2];
rz(-0.61389273) q[2];
sx q[2];
rz(-2.4284412) q[2];
rz(-2.5085416) q[3];
sx q[3];
rz(-0.96708599) q[3];
sx q[3];
rz(2.2658394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1322094) q[0];
sx q[0];
rz(-1.3306946) q[0];
sx q[0];
rz(2.7015986) q[0];
rz(0.72883365) q[1];
sx q[1];
rz(-0.76580096) q[1];
sx q[1];
rz(-2.0892807) q[1];
rz(-1.1807146) q[2];
sx q[2];
rz(-1.6772391) q[2];
sx q[2];
rz(-0.46628484) q[2];
rz(-0.80753978) q[3];
sx q[3];
rz(-1.3046466) q[3];
sx q[3];
rz(-0.70481963) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
