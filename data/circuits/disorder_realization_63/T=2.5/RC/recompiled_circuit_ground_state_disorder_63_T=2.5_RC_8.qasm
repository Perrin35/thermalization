OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.11290057) q[0];
sx q[0];
rz(-0.24554645) q[0];
sx q[0];
rz(3.088933) q[0];
rz(1.117299) q[1];
sx q[1];
rz(-2.8007562) q[1];
sx q[1];
rz(-2.053082) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8262622) q[0];
sx q[0];
rz(-1.966094) q[0];
sx q[0];
rz(0.068282337) q[0];
rz(-pi) q[1];
rz(0.91204466) q[2];
sx q[2];
rz(-0.17870644) q[2];
sx q[2];
rz(-2.2026874) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0170235) q[1];
sx q[1];
rz(-0.93734159) q[1];
sx q[1];
rz(-1.3658466) q[1];
x q[2];
rz(1.151122) q[3];
sx q[3];
rz(-1.5745145) q[3];
sx q[3];
rz(-0.17824717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5071621) q[2];
sx q[2];
rz(-1.3619225) q[2];
sx q[2];
rz(-0.6081028) q[2];
rz(1.1242584) q[3];
sx q[3];
rz(-1.2048771) q[3];
sx q[3];
rz(2.2500136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7267777) q[0];
sx q[0];
rz(-1.5685273) q[0];
sx q[0];
rz(2.538105) q[0];
rz(0.51015774) q[1];
sx q[1];
rz(-2.3699103) q[1];
sx q[1];
rz(2.1147494) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0870846) q[0];
sx q[0];
rz(-1.3985976) q[0];
sx q[0];
rz(1.728375) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2910094) q[2];
sx q[2];
rz(-2.2424978) q[2];
sx q[2];
rz(-1.6784061) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.660381) q[1];
sx q[1];
rz(-0.91970316) q[1];
sx q[1];
rz(-0.91459413) q[1];
rz(-pi) q[2];
rz(0.90791865) q[3];
sx q[3];
rz(-1.1343805) q[3];
sx q[3];
rz(-1.7493356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5828731) q[2];
sx q[2];
rz(-1.1483973) q[2];
sx q[2];
rz(-2.4951475) q[2];
rz(1.8630113) q[3];
sx q[3];
rz(-2.1918112) q[3];
sx q[3];
rz(1.4518552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2090787) q[0];
sx q[0];
rz(-3.0113853) q[0];
sx q[0];
rz(1.5694438) q[0];
rz(0.34700829) q[1];
sx q[1];
rz(-1.6667112) q[1];
sx q[1];
rz(-2.2775547) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61053172) q[0];
sx q[0];
rz(-0.60734925) q[0];
sx q[0];
rz(0.37215287) q[0];
rz(0.43564312) q[2];
sx q[2];
rz(-1.3416939) q[2];
sx q[2];
rz(-1.733591) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4260508) q[1];
sx q[1];
rz(-2.155039) q[1];
sx q[1];
rz(2.3795695) q[1];
x q[2];
rz(0.10580833) q[3];
sx q[3];
rz(-1.0472626) q[3];
sx q[3];
rz(2.2744655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0916834) q[2];
sx q[2];
rz(-0.20541643) q[2];
sx q[2];
rz(-1.4570215) q[2];
rz(2.125804) q[3];
sx q[3];
rz(-2.6088645) q[3];
sx q[3];
rz(0.31639019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9006573) q[0];
sx q[0];
rz(-0.37739402) q[0];
sx q[0];
rz(1.5628016) q[0];
rz(-1.0904301) q[1];
sx q[1];
rz(-0.54160392) q[1];
sx q[1];
rz(-1.9900367) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3061515) q[0];
sx q[0];
rz(-0.62850289) q[0];
sx q[0];
rz(2.7449904) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0959553) q[2];
sx q[2];
rz(-0.86855382) q[2];
sx q[2];
rz(-2.4202731) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.71660605) q[1];
sx q[1];
rz(-2.0946436) q[1];
sx q[1];
rz(0.040158466) q[1];
x q[2];
rz(3.1138469) q[3];
sx q[3];
rz(-0.96409269) q[3];
sx q[3];
rz(1.9880027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2942723) q[2];
sx q[2];
rz(-0.5449833) q[2];
sx q[2];
rz(-1.6418183) q[2];
rz(-0.13449399) q[3];
sx q[3];
rz(-1.8650863) q[3];
sx q[3];
rz(-2.3597778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0624369) q[0];
sx q[0];
rz(-2.2298614) q[0];
sx q[0];
rz(2.675918) q[0];
rz(1.8648719) q[1];
sx q[1];
rz(-1.0036422) q[1];
sx q[1];
rz(-1.0328971) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5860234) q[0];
sx q[0];
rz(-0.15842552) q[0];
sx q[0];
rz(1.6371631) q[0];
x q[1];
rz(3.0014702) q[2];
sx q[2];
rz(-1.6701785) q[2];
sx q[2];
rz(1.4750208) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3001662) q[1];
sx q[1];
rz(-0.71117095) q[1];
sx q[1];
rz(-2.7513642) q[1];
x q[2];
rz(1.0276658) q[3];
sx q[3];
rz(-1.9082205) q[3];
sx q[3];
rz(-1.3285445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5501962) q[2];
sx q[2];
rz(-1.1893716) q[2];
sx q[2];
rz(-2.2097394) q[2];
rz(3.0469117) q[3];
sx q[3];
rz(-0.90016142) q[3];
sx q[3];
rz(-1.3100821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4513627) q[0];
sx q[0];
rz(-0.19915038) q[0];
sx q[0];
rz(-3.0354101) q[0];
rz(0.55446082) q[1];
sx q[1];
rz(-2.3531871) q[1];
sx q[1];
rz(-1.5415446) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7035737) q[0];
sx q[0];
rz(-0.39123639) q[0];
sx q[0];
rz(-0.06322579) q[0];
rz(-2.0279858) q[2];
sx q[2];
rz(-1.3234252) q[2];
sx q[2];
rz(-2.8988225) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2682802) q[1];
sx q[1];
rz(-2.2619216) q[1];
sx q[1];
rz(1.7296438) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8951449) q[3];
sx q[3];
rz(-1.5165899) q[3];
sx q[3];
rz(1.1185631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9083378) q[2];
sx q[2];
rz(-0.7408064) q[2];
sx q[2];
rz(1.6439269) q[2];
rz(0.46105591) q[3];
sx q[3];
rz(-2.4888829) q[3];
sx q[3];
rz(1.8259995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6719565) q[0];
sx q[0];
rz(-0.75103432) q[0];
sx q[0];
rz(2.0627956) q[0];
rz(-2.5476593) q[1];
sx q[1];
rz(-0.97336951) q[1];
sx q[1];
rz(-0.22720164) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1329472) q[0];
sx q[0];
rz(-1.7194178) q[0];
sx q[0];
rz(0.11283837) q[0];
rz(-pi) q[1];
rz(3.0005713) q[2];
sx q[2];
rz(-1.2407836) q[2];
sx q[2];
rz(-0.38244707) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1006443) q[1];
sx q[1];
rz(-2.4923263) q[1];
sx q[1];
rz(-2.098564) q[1];
rz(-pi) q[2];
rz(3.041275) q[3];
sx q[3];
rz(-1.2004108) q[3];
sx q[3];
rz(-0.60587347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6195153) q[2];
sx q[2];
rz(-2.7232813) q[2];
sx q[2];
rz(-2.412001) q[2];
rz(1.4531762) q[3];
sx q[3];
rz(-1.4540693) q[3];
sx q[3];
rz(-2.5280473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5486117) q[0];
sx q[0];
rz(-0.91658968) q[0];
sx q[0];
rz(0.35162893) q[0];
rz(-2.8278606) q[1];
sx q[1];
rz(-1.2064563) q[1];
sx q[1];
rz(-1.3765913) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54395318) q[0];
sx q[0];
rz(-0.40297976) q[0];
sx q[0];
rz(0.12059327) q[0];
rz(-0.96504546) q[2];
sx q[2];
rz(-2.9403983) q[2];
sx q[2];
rz(-0.66058841) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9326068) q[1];
sx q[1];
rz(-0.58949215) q[1];
sx q[1];
rz(-2.8578651) q[1];
x q[2];
rz(1.9405211) q[3];
sx q[3];
rz(-1.4784357) q[3];
sx q[3];
rz(2.8522688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.54479638) q[2];
sx q[2];
rz(-2.4185541) q[2];
sx q[2];
rz(-1.5210305) q[2];
rz(-0.14791402) q[3];
sx q[3];
rz(-1.3444129) q[3];
sx q[3];
rz(1.3676876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77383298) q[0];
sx q[0];
rz(-1.8748883) q[0];
sx q[0];
rz(-1.897478) q[0];
rz(1.4338088) q[1];
sx q[1];
rz(-1.560874) q[1];
sx q[1];
rz(0.22020766) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51802902) q[0];
sx q[0];
rz(-1.9401778) q[0];
sx q[0];
rz(-0.9250888) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.979109) q[2];
sx q[2];
rz(-2.0318609) q[2];
sx q[2];
rz(1.7562255) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9088194) q[1];
sx q[1];
rz(-2.4202883) q[1];
sx q[1];
rz(0.90296332) q[1];
rz(-pi) q[2];
rz(3.0433398) q[3];
sx q[3];
rz(-1.6643401) q[3];
sx q[3];
rz(0.52631179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9571017) q[2];
sx q[2];
rz(-1.1471006) q[2];
sx q[2];
rz(1.6415049) q[2];
rz(3.0575276) q[3];
sx q[3];
rz(-1.2596687) q[3];
sx q[3];
rz(-3.0715004) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67749196) q[0];
sx q[0];
rz(-0.75834948) q[0];
sx q[0];
rz(0.41859928) q[0];
rz(-2.1981926) q[1];
sx q[1];
rz(-1.489233) q[1];
sx q[1];
rz(-2.8387866) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43234149) q[0];
sx q[0];
rz(-2.925207) q[0];
sx q[0];
rz(2.2585013) q[0];
rz(-pi) q[1];
x q[1];
rz(0.19089107) q[2];
sx q[2];
rz(-1.7451236) q[2];
sx q[2];
rz(-3.0534923) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8119366) q[1];
sx q[1];
rz(-1.6534263) q[1];
sx q[1];
rz(1.5272115) q[1];
rz(-1.2556158) q[3];
sx q[3];
rz(-0.36096301) q[3];
sx q[3];
rz(0.40004594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4960949) q[2];
sx q[2];
rz(-0.87254137) q[2];
sx q[2];
rz(-2.3910451) q[2];
rz(-1.6802855) q[3];
sx q[3];
rz(-0.76992005) q[3];
sx q[3];
rz(2.6245978) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0208329) q[0];
sx q[0];
rz(-1.5789541) q[0];
sx q[0];
rz(-1.5776237) q[0];
rz(1.9137406) q[1];
sx q[1];
rz(-0.49929437) q[1];
sx q[1];
rz(1.1566537) q[1];
rz(0.56671178) q[2];
sx q[2];
rz(-1.7253582) q[2];
sx q[2];
rz(2.2127989) q[2];
rz(-2.078656) q[3];
sx q[3];
rz(-2.776317) q[3];
sx q[3];
rz(-1.1767514) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
