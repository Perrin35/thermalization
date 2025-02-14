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
rz(0.077724783) q[0];
sx q[0];
rz(-2.1169777) q[0];
sx q[0];
rz(-1.2999363) q[0];
rz(2.6990702) q[1];
sx q[1];
rz(4.2808851) q[1];
sx q[1];
rz(12.130375) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76749698) q[0];
sx q[0];
rz(-2.7347235) q[0];
sx q[0];
rz(-1.8396729) q[0];
rz(-2.9098815) q[2];
sx q[2];
rz(-0.75318906) q[2];
sx q[2];
rz(-1.1540268) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.30937815) q[1];
sx q[1];
rz(-1.16799) q[1];
sx q[1];
rz(-1.554053) q[1];
rz(-pi) q[2];
rz(0.12650872) q[3];
sx q[3];
rz(-0.96130575) q[3];
sx q[3];
rz(1.786527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0787597) q[2];
sx q[2];
rz(-1.6597513) q[2];
sx q[2];
rz(-0.01595846) q[2];
rz(-1.4912841) q[3];
sx q[3];
rz(-1.242638) q[3];
sx q[3];
rz(-2.7119467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66351873) q[0];
sx q[0];
rz(-1.0915382) q[0];
sx q[0];
rz(-1.7953405) q[0];
rz(1.1675872) q[1];
sx q[1];
rz(-1.0941894) q[1];
sx q[1];
rz(-1.2487372) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6561476) q[0];
sx q[0];
rz(-2.1753575) q[0];
sx q[0];
rz(0.47550918) q[0];
rz(-pi) q[1];
rz(2.4820868) q[2];
sx q[2];
rz(-2.9444866) q[2];
sx q[2];
rz(1.3172305) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.389321) q[1];
sx q[1];
rz(-1.5157425) q[1];
sx q[1];
rz(-1.9717713) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.98442673) q[3];
sx q[3];
rz(-1.274144) q[3];
sx q[3];
rz(-3.0595317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.0060129082) q[2];
sx q[2];
rz(-0.69301444) q[2];
sx q[2];
rz(-2.6683624) q[2];
rz(3.0114975) q[3];
sx q[3];
rz(-1.7546763) q[3];
sx q[3];
rz(0.010802833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8148282) q[0];
sx q[0];
rz(-0.15693754) q[0];
sx q[0];
rz(-2.9275295) q[0];
rz(-1.7474984) q[1];
sx q[1];
rz(-0.97624818) q[1];
sx q[1];
rz(3.0551547) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8562216) q[0];
sx q[0];
rz(-0.32111327) q[0];
sx q[0];
rz(1.1620528) q[0];
rz(-pi) q[1];
rz(-2.2485224) q[2];
sx q[2];
rz(-1.1249591) q[2];
sx q[2];
rz(2.0634212) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9443892) q[1];
sx q[1];
rz(-1.3968588) q[1];
sx q[1];
rz(-2.7157213) q[1];
rz(-pi) q[2];
rz(-2.6085601) q[3];
sx q[3];
rz(-2.1500476) q[3];
sx q[3];
rz(-1.2377626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4702686) q[2];
sx q[2];
rz(-0.9321804) q[2];
sx q[2];
rz(1.1364802) q[2];
rz(-1.7911576) q[3];
sx q[3];
rz(-1.8108436) q[3];
sx q[3];
rz(2.7854846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7478624) q[0];
sx q[0];
rz(-2.6342454) q[0];
sx q[0];
rz(2.2407929) q[0];
rz(-1.9081217) q[1];
sx q[1];
rz(-0.8546468) q[1];
sx q[1];
rz(-2.5305117) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40070286) q[0];
sx q[0];
rz(-0.4755334) q[0];
sx q[0];
rz(-1.5551546) q[0];
rz(-2.99154) q[2];
sx q[2];
rz(-0.23153472) q[2];
sx q[2];
rz(0.48047149) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7984182) q[1];
sx q[1];
rz(-1.1504984) q[1];
sx q[1];
rz(-2.0139) q[1];
rz(-1.4540205) q[3];
sx q[3];
rz(-1.2216785) q[3];
sx q[3];
rz(-1.735183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5078807) q[2];
sx q[2];
rz(-2.4780126) q[2];
sx q[2];
rz(-3.0207685) q[2];
rz(-3.1145596) q[3];
sx q[3];
rz(-2.9398672) q[3];
sx q[3];
rz(0.87593186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
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
rz(0.19288572) q[0];
sx q[0];
rz(-2.3159733) q[0];
sx q[0];
rz(1.8008308) q[0];
rz(-1.5277398) q[1];
sx q[1];
rz(-1.3132881) q[1];
sx q[1];
rz(0.54940474) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2717239) q[0];
sx q[0];
rz(-0.74942526) q[0];
sx q[0];
rz(-0.11336993) q[0];
x q[1];
rz(-2.1020736) q[2];
sx q[2];
rz(-1.5434859) q[2];
sx q[2];
rz(-2.4940048) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0548965) q[1];
sx q[1];
rz(-2.6532996) q[1];
sx q[1];
rz(-1.4642621) q[1];
x q[2];
rz(-0.033614393) q[3];
sx q[3];
rz(-2.0319124) q[3];
sx q[3];
rz(-0.5881084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1250829) q[2];
sx q[2];
rz(-1.5089401) q[2];
sx q[2];
rz(-0.58270085) q[2];
rz(-1.6216283) q[3];
sx q[3];
rz(-0.71526066) q[3];
sx q[3];
rz(-1.8250072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80578605) q[0];
sx q[0];
rz(-1.0842706) q[0];
sx q[0];
rz(1.2544607) q[0];
rz(1.1955903) q[1];
sx q[1];
rz(-1.4072199) q[1];
sx q[1];
rz(-1.6244627) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82255581) q[0];
sx q[0];
rz(-1.6375038) q[0];
sx q[0];
rz(0.011724756) q[0];
x q[1];
rz(-2.1468494) q[2];
sx q[2];
rz(-2.2855504) q[2];
sx q[2];
rz(1.7896259) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2020009) q[1];
sx q[1];
rz(-0.67145214) q[1];
sx q[1];
rz(2.6670649) q[1];
rz(-1.5593525) q[3];
sx q[3];
rz(-1.0557888) q[3];
sx q[3];
rz(2.8449165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.967041) q[2];
sx q[2];
rz(-1.5341362) q[2];
sx q[2];
rz(-3.0957481) q[2];
rz(-0.52715078) q[3];
sx q[3];
rz(-2.1093301) q[3];
sx q[3];
rz(0.84120685) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6531649) q[0];
sx q[0];
rz(-0.66075745) q[0];
sx q[0];
rz(0.38594693) q[0];
rz(-1.376232) q[1];
sx q[1];
rz(-2.0538581) q[1];
sx q[1];
rz(1.8650581) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4285884) q[0];
sx q[0];
rz(-2.3516293) q[0];
sx q[0];
rz(-0.53623627) q[0];
rz(1.1982199) q[2];
sx q[2];
rz(-1.5261006) q[2];
sx q[2];
rz(0.17105779) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1305728) q[1];
sx q[1];
rz(-1.5351194) q[1];
sx q[1];
rz(3.1336354) q[1];
x q[2];
rz(-0.035865772) q[3];
sx q[3];
rz(-1.6350785) q[3];
sx q[3];
rz(2.4202297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4253) q[2];
sx q[2];
rz(-1.7159117) q[2];
sx q[2];
rz(1.3670134) q[2];
rz(-0.084271757) q[3];
sx q[3];
rz(-2.027498) q[3];
sx q[3];
rz(-3.058694) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0609584) q[0];
sx q[0];
rz(-3.0301889) q[0];
sx q[0];
rz(1.2782619) q[0];
rz(-0.06079611) q[1];
sx q[1];
rz(-0.95525974) q[1];
sx q[1];
rz(-1.0999058) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17717028) q[0];
sx q[0];
rz(-1.9510498) q[0];
sx q[0];
rz(2.9833262) q[0];
x q[1];
rz(2.3448337) q[2];
sx q[2];
rz(-0.9769868) q[2];
sx q[2];
rz(-1.5409005) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2437224) q[1];
sx q[1];
rz(-1.9258782) q[1];
sx q[1];
rz(1.8753858) q[1];
x q[2];
rz(2.0792316) q[3];
sx q[3];
rz(-2.1223801) q[3];
sx q[3];
rz(-2.036117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3595769) q[2];
sx q[2];
rz(-2.4961553) q[2];
sx q[2];
rz(-2.6673356) q[2];
rz(0.54840243) q[3];
sx q[3];
rz(-2.4689398) q[3];
sx q[3];
rz(-1.3801581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43535522) q[0];
sx q[0];
rz(-0.42689231) q[0];
sx q[0];
rz(1.2109582) q[0];
rz(-1.8636761) q[1];
sx q[1];
rz(-1.5958818) q[1];
sx q[1];
rz(-0.011215297) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0408913) q[0];
sx q[0];
rz(-2.5731003) q[0];
sx q[0];
rz(-2.9799653) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8008046) q[2];
sx q[2];
rz(-1.9267779) q[2];
sx q[2];
rz(-1.180151) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1234731) q[1];
sx q[1];
rz(-1.5168861) q[1];
sx q[1];
rz(-2.3192295) q[1];
rz(-pi) q[2];
rz(-0.09402676) q[3];
sx q[3];
rz(-1.8593899) q[3];
sx q[3];
rz(1.2666324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8212905) q[2];
sx q[2];
rz(-1.9052637) q[2];
sx q[2];
rz(-0.75766364) q[2];
rz(0.6997987) q[3];
sx q[3];
rz(-0.57707969) q[3];
sx q[3];
rz(0.9052161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59117544) q[0];
sx q[0];
rz(-1.9015522) q[0];
sx q[0];
rz(2.8072939) q[0];
rz(2.7475157) q[1];
sx q[1];
rz(-0.95029345) q[1];
sx q[1];
rz(-1.2324415) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.619433) q[0];
sx q[0];
rz(-1.8364826) q[0];
sx q[0];
rz(2.9152318) q[0];
x q[1];
rz(1.7712461) q[2];
sx q[2];
rz(-0.88267372) q[2];
sx q[2];
rz(1.3039948) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5345911) q[1];
sx q[1];
rz(-1.2950194) q[1];
sx q[1];
rz(2.9124583) q[1];
rz(-pi) q[2];
rz(-2.0471935) q[3];
sx q[3];
rz(-0.55022424) q[3];
sx q[3];
rz(-0.27980313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.95840994) q[2];
sx q[2];
rz(-0.56768688) q[2];
sx q[2];
rz(-0.063610345) q[2];
rz(2.375864) q[3];
sx q[3];
rz(-1.1252517) q[3];
sx q[3];
rz(3.0797899) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3381989) q[0];
sx q[0];
rz(-2.3173208) q[0];
sx q[0];
rz(1.4650387) q[0];
rz(1.4631396) q[1];
sx q[1];
rz(-1.2668162) q[1];
sx q[1];
rz(-0.75513671) q[1];
rz(-1.4037618) q[2];
sx q[2];
rz(-0.93334953) q[2];
sx q[2];
rz(2.2095528) q[2];
rz(-2.3018671) q[3];
sx q[3];
rz(-1.4473549) q[3];
sx q[3];
rz(2.2458129) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
