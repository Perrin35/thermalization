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
rz(2.8960462) q[0];
sx q[0];
rz(9.4774376) q[0];
rz(-2.0242937) q[1];
sx q[1];
rz(-0.3408365) q[1];
sx q[1];
rz(-1.0885106) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4910925) q[0];
sx q[0];
rz(-0.40084565) q[0];
sx q[0];
rz(1.7328788) q[0];
rz(-pi) q[1];
rz(2.229548) q[2];
sx q[2];
rz(-2.9628862) q[2];
sx q[2];
rz(-2.2026874) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3237985) q[1];
sx q[1];
rz(-1.7355738) q[1];
sx q[1];
rz(2.4980157) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1375211) q[3];
sx q[3];
rz(-1.1511251) q[3];
sx q[3];
rz(-1.3908902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6344305) q[2];
sx q[2];
rz(-1.3619225) q[2];
sx q[2];
rz(-2.5334899) q[2];
rz(-1.1242584) q[3];
sx q[3];
rz(-1.2048771) q[3];
sx q[3];
rz(-2.2500136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7267777) q[0];
sx q[0];
rz(-1.5685273) q[0];
sx q[0];
rz(-2.538105) q[0];
rz(-2.6314349) q[1];
sx q[1];
rz(-0.77168232) q[1];
sx q[1];
rz(1.0268432) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3393766) q[0];
sx q[0];
rz(-0.23288865) q[0];
sx q[0];
rz(-0.73407895) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.449069) q[2];
sx q[2];
rz(-2.1998458) q[2];
sx q[2];
rz(0.72450996) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7945343) q[1];
sx q[1];
rz(-2.0774986) q[1];
sx q[1];
rz(-2.3757412) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2193416) q[3];
sx q[3];
rz(-0.77510683) q[3];
sx q[3];
rz(2.4668281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5828731) q[2];
sx q[2];
rz(-1.1483973) q[2];
sx q[2];
rz(0.64644512) q[2];
rz(1.8630113) q[3];
sx q[3];
rz(-2.1918112) q[3];
sx q[3];
rz(-1.6897374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2090787) q[0];
sx q[0];
rz(-3.0113853) q[0];
sx q[0];
rz(-1.5694438) q[0];
rz(-2.7945844) q[1];
sx q[1];
rz(-1.6667112) q[1];
sx q[1];
rz(-2.2775547) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65007675) q[0];
sx q[0];
rz(-1.7798316) q[0];
sx q[0];
rz(-2.5670426) q[0];
x q[1];
rz(-0.50484294) q[2];
sx q[2];
rz(-2.6528203) q[2];
sx q[2];
rz(-2.8504643) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7622704) q[1];
sx q[1];
rz(-2.218609) q[1];
sx q[1];
rz(0.76384441) q[1];
rz(-1.7517463) q[3];
sx q[3];
rz(-0.53314185) q[3];
sx q[3];
rz(2.0651434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0916834) q[2];
sx q[2];
rz(-0.20541643) q[2];
sx q[2];
rz(-1.6845711) q[2];
rz(2.125804) q[3];
sx q[3];
rz(-2.6088645) q[3];
sx q[3];
rz(0.31639019) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24093534) q[0];
sx q[0];
rz(-0.37739402) q[0];
sx q[0];
rz(-1.5628016) q[0];
rz(2.0511625) q[1];
sx q[1];
rz(-2.5999887) q[1];
sx q[1];
rz(1.9900367) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82842302) q[0];
sx q[0];
rz(-0.99762625) q[0];
sx q[0];
rz(1.8445119) q[0];
rz(0.49539613) q[2];
sx q[2];
rz(-2.317111) q[2];
sx q[2];
rz(1.7480614) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.87428625) q[1];
sx q[1];
rz(-1.5360254) q[1];
sx q[1];
rz(-2.094993) q[1];
x q[2];
rz(1.6107481) q[3];
sx q[3];
rz(-2.5343347) q[3];
sx q[3];
rz(1.1049529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2942723) q[2];
sx q[2];
rz(-2.5966094) q[2];
sx q[2];
rz(1.6418183) q[2];
rz(-3.0070987) q[3];
sx q[3];
rz(-1.2765063) q[3];
sx q[3];
rz(0.78181481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0624369) q[0];
sx q[0];
rz(-0.9117313) q[0];
sx q[0];
rz(2.675918) q[0];
rz(1.8648719) q[1];
sx q[1];
rz(-2.1379505) q[1];
sx q[1];
rz(-2.1086955) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0807652) q[0];
sx q[0];
rz(-1.5812591) q[0];
sx q[0];
rz(-1.4127138) q[0];
rz(-pi) q[1];
rz(1.470437) q[2];
sx q[2];
rz(-1.7102229) q[2];
sx q[2];
rz(-3.0318236) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1101428) q[1];
sx q[1];
rz(-1.3198766) q[1];
sx q[1];
rz(-0.67278426) q[1];
rz(-pi) q[2];
rz(0.38893969) q[3];
sx q[3];
rz(-2.0802214) q[3];
sx q[3];
rz(0.044993613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5501962) q[2];
sx q[2];
rz(-1.1893716) q[2];
sx q[2];
rz(-2.2097394) q[2];
rz(3.0469117) q[3];
sx q[3];
rz(-2.2414312) q[3];
sx q[3];
rz(-1.8315106) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4513627) q[0];
sx q[0];
rz(-0.19915038) q[0];
sx q[0];
rz(0.1061826) q[0];
rz(0.55446082) q[1];
sx q[1];
rz(-0.78840557) q[1];
sx q[1];
rz(1.5415446) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5063972) q[0];
sx q[0];
rz(-1.180384) q[0];
sx q[0];
rz(1.5447389) q[0];
rz(-pi) q[1];
rz(-1.0511398) q[2];
sx q[2];
rz(-0.51562947) q[2];
sx q[2];
rz(1.7898498) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5145077) q[1];
sx q[1];
rz(-0.70620757) q[1];
sx q[1];
rz(-2.9526564) q[1];
rz(-2.8951449) q[3];
sx q[3];
rz(-1.6250027) q[3];
sx q[3];
rz(2.0230296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9083378) q[2];
sx q[2];
rz(-0.7408064) q[2];
sx q[2];
rz(-1.6439269) q[2];
rz(-0.46105591) q[3];
sx q[3];
rz(-0.65270972) q[3];
sx q[3];
rz(-1.3155931) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4696362) q[0];
sx q[0];
rz(-0.75103432) q[0];
sx q[0];
rz(-2.0627956) q[0];
rz(-0.59393334) q[1];
sx q[1];
rz(-2.1682231) q[1];
sx q[1];
rz(-0.22720164) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5626635) q[0];
sx q[0];
rz(-1.6823856) q[0];
sx q[0];
rz(-1.4212378) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0005713) q[2];
sx q[2];
rz(-1.900809) q[2];
sx q[2];
rz(-0.38244707) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1006443) q[1];
sx q[1];
rz(-0.64926636) q[1];
sx q[1];
rz(1.0430286) q[1];
rz(-pi) q[2];
x q[2];
rz(0.10031767) q[3];
sx q[3];
rz(-1.2004108) q[3];
sx q[3];
rz(0.60587347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5220773) q[2];
sx q[2];
rz(-2.7232813) q[2];
sx q[2];
rz(0.72959161) q[2];
rz(-1.4531762) q[3];
sx q[3];
rz(-1.6875234) q[3];
sx q[3];
rz(0.61354536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59298092) q[0];
sx q[0];
rz(-2.225003) q[0];
sx q[0];
rz(-2.7899637) q[0];
rz(2.8278606) q[1];
sx q[1];
rz(-1.2064563) q[1];
sx q[1];
rz(-1.7650013) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91582725) q[0];
sx q[0];
rz(-1.5236014) q[0];
sx q[0];
rz(2.7412358) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4046762) q[2];
sx q[2];
rz(-1.4567647) q[2];
sx q[2];
rz(-1.5064552) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9326068) q[1];
sx q[1];
rz(-0.58949215) q[1];
sx q[1];
rz(2.8578651) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8217161) q[3];
sx q[3];
rz(-2.7610215) q[3];
sx q[3];
rz(-1.5151092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5967963) q[2];
sx q[2];
rz(-0.72303855) q[2];
sx q[2];
rz(1.5210305) q[2];
rz(0.14791402) q[3];
sx q[3];
rz(-1.3444129) q[3];
sx q[3];
rz(-1.3676876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77383298) q[0];
sx q[0];
rz(-1.2667043) q[0];
sx q[0];
rz(-1.2441147) q[0];
rz(1.4338088) q[1];
sx q[1];
rz(-1.560874) q[1];
sx q[1];
rz(-2.921385) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51802902) q[0];
sx q[0];
rz(-1.9401778) q[0];
sx q[0];
rz(0.9250888) q[0];
rz(-1.8856064) q[2];
sx q[2];
rz(-0.48690912) q[2];
sx q[2];
rz(-2.1092507) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.23277321) q[1];
sx q[1];
rz(-2.4202883) q[1];
sx q[1];
rz(0.90296332) q[1];
x q[2];
rz(1.4768019) q[3];
sx q[3];
rz(-1.6686182) q[3];
sx q[3];
rz(-2.0879012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1844909) q[2];
sx q[2];
rz(-1.9944921) q[2];
sx q[2];
rz(-1.5000878) q[2];
rz(-3.0575276) q[3];
sx q[3];
rz(-1.2596687) q[3];
sx q[3];
rz(3.0715004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67749196) q[0];
sx q[0];
rz(-2.3832432) q[0];
sx q[0];
rz(0.41859928) q[0];
rz(0.94340008) q[1];
sx q[1];
rz(-1.6523596) q[1];
sx q[1];
rz(-0.30280608) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26696268) q[0];
sx q[0];
rz(-1.4041252) q[0];
sx q[0];
rz(-0.13864362) q[0];
x q[1];
rz(-0.19089107) q[2];
sx q[2];
rz(-1.7451236) q[2];
sx q[2];
rz(-0.088100351) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.23754072) q[1];
sx q[1];
rz(-1.5273603) q[1];
sx q[1];
rz(-0.082708184) q[1];
x q[2];
rz(0.11649152) q[3];
sx q[3];
rz(-1.2283652) q[3];
sx q[3];
rz(-0.064700944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4960949) q[2];
sx q[2];
rz(-0.87254137) q[2];
sx q[2];
rz(-0.75054753) q[2];
rz(1.6802855) q[3];
sx q[3];
rz(-0.76992005) q[3];
sx q[3];
rz(-2.6245978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0208329) q[0];
sx q[0];
rz(-1.5789541) q[0];
sx q[0];
rz(-1.5776237) q[0];
rz(-1.9137406) q[1];
sx q[1];
rz(-2.6422983) q[1];
sx q[1];
rz(-1.9849389) q[1];
rz(2.5748809) q[2];
sx q[2];
rz(-1.4162345) q[2];
sx q[2];
rz(-0.92879374) q[2];
rz(-1.0629366) q[3];
sx q[3];
rz(-0.36527562) q[3];
sx q[3];
rz(1.9648413) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
