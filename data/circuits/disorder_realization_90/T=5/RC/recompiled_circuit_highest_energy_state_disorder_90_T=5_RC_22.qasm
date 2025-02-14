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
rz(2.4861205) q[0];
sx q[0];
rz(-0.95215005) q[0];
sx q[0];
rz(3.0478391) q[0];
rz(-1.7733511) q[1];
sx q[1];
rz(-1.5612839) q[1];
sx q[1];
rz(-2.4743647) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.102948) q[0];
sx q[0];
rz(-1.1955452) q[0];
sx q[0];
rz(1.0374336) q[0];
rz(-pi) q[1];
rz(0.25847659) q[2];
sx q[2];
rz(-1.3745068) q[2];
sx q[2];
rz(2.1854913) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.92541313) q[1];
sx q[1];
rz(-2.254018) q[1];
sx q[1];
rz(0.36551468) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5205403) q[3];
sx q[3];
rz(-2.1571311) q[3];
sx q[3];
rz(2.1042657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5253456) q[2];
sx q[2];
rz(-1.5364001) q[2];
sx q[2];
rz(3.1180535) q[2];
rz(2.8871138) q[3];
sx q[3];
rz(-1.1602217) q[3];
sx q[3];
rz(-2.3833073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1940521) q[0];
sx q[0];
rz(-1.7451311) q[0];
sx q[0];
rz(-0.61554712) q[0];
rz(0.60000348) q[1];
sx q[1];
rz(-1.2702076) q[1];
sx q[1];
rz(-2.792865) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.00086) q[0];
sx q[0];
rz(-1.6475186) q[0];
sx q[0];
rz(2.4984635) q[0];
rz(-pi) q[1];
rz(-2.8519657) q[2];
sx q[2];
rz(-1.2420601) q[2];
sx q[2];
rz(1.6148953) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.80417577) q[1];
sx q[1];
rz(-1.2722208) q[1];
sx q[1];
rz(-1.0945303) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0660961) q[3];
sx q[3];
rz(-1.3196006) q[3];
sx q[3];
rz(-0.78758729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.49150026) q[2];
sx q[2];
rz(-1.3050175) q[2];
sx q[2];
rz(0.65518641) q[2];
rz(0.16264597) q[3];
sx q[3];
rz(-0.9442257) q[3];
sx q[3];
rz(0.20774016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96838897) q[0];
sx q[0];
rz(-2.0252616) q[0];
sx q[0];
rz(0.094245687) q[0];
rz(1.3000129) q[1];
sx q[1];
rz(-0.899122) q[1];
sx q[1];
rz(-1.947044) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.081095) q[0];
sx q[0];
rz(-2.0682242) q[0];
sx q[0];
rz(1.917385) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0897806) q[2];
sx q[2];
rz(-1.2070884) q[2];
sx q[2];
rz(1.6339982) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.49385168) q[1];
sx q[1];
rz(-1.6147227) q[1];
sx q[1];
rz(-1.8274587) q[1];
rz(-pi) q[2];
rz(2.7500509) q[3];
sx q[3];
rz(-2.3394008) q[3];
sx q[3];
rz(-1.6959617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3208348) q[2];
sx q[2];
rz(-1.4471549) q[2];
sx q[2];
rz(-2.7336332) q[2];
rz(1.3784493) q[3];
sx q[3];
rz(-2.2429376) q[3];
sx q[3];
rz(-0.11817008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(2.5648062) q[0];
sx q[0];
rz(-1.7498359) q[0];
sx q[0];
rz(-0.80192649) q[0];
rz(0.39464125) q[1];
sx q[1];
rz(-2.092974) q[1];
sx q[1];
rz(3.1065497) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5749767) q[0];
sx q[0];
rz(-2.1144951) q[0];
sx q[0];
rz(-0.4369785) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3801244) q[2];
sx q[2];
rz(-1.241127) q[2];
sx q[2];
rz(-2.6063188) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1890244) q[1];
sx q[1];
rz(-0.060563001) q[1];
sx q[1];
rz(0.80068717) q[1];
rz(-pi) q[2];
rz(2.7862042) q[3];
sx q[3];
rz(-1.3878763) q[3];
sx q[3];
rz(1.8754043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0617712) q[2];
sx q[2];
rz(-1.062919) q[2];
sx q[2];
rz(1.7064095) q[2];
rz(1.4052514) q[3];
sx q[3];
rz(-1.0900213) q[3];
sx q[3];
rz(-2.0315571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6768796) q[0];
sx q[0];
rz(-1.1486624) q[0];
sx q[0];
rz(1.1418463) q[0];
rz(-2.6308718) q[1];
sx q[1];
rz(-1.9074214) q[1];
sx q[1];
rz(1.7150735) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1175809) q[0];
sx q[0];
rz(-2.3387351) q[0];
sx q[0];
rz(0.76501655) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2586589) q[2];
sx q[2];
rz(-1.3545879) q[2];
sx q[2];
rz(0.59528159) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.99014501) q[1];
sx q[1];
rz(-1.4583318) q[1];
sx q[1];
rz(-2.165035) q[1];
rz(-1.9782039) q[3];
sx q[3];
rz(-1.7247044) q[3];
sx q[3];
rz(2.8606382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4003754) q[2];
sx q[2];
rz(-1.8343238) q[2];
sx q[2];
rz(-0.2571787) q[2];
rz(-0.87279618) q[3];
sx q[3];
rz(-1.8200487) q[3];
sx q[3];
rz(2.0040373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0208825) q[0];
sx q[0];
rz(-2.6873984) q[0];
sx q[0];
rz(-2.0632451) q[0];
rz(-2.2348166) q[1];
sx q[1];
rz(-0.6858784) q[1];
sx q[1];
rz(0.041898601) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3495765) q[0];
sx q[0];
rz(-1.0984165) q[0];
sx q[0];
rz(-0.040332009) q[0];
rz(-pi) q[1];
rz(-0.12280913) q[2];
sx q[2];
rz(-0.90722668) q[2];
sx q[2];
rz(-1.0994436) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7028693) q[1];
sx q[1];
rz(-1.8449191) q[1];
sx q[1];
rz(-1.1945748) q[1];
rz(-pi) q[2];
x q[2];
rz(0.83591977) q[3];
sx q[3];
rz(-1.984388) q[3];
sx q[3];
rz(1.5121764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8580253) q[2];
sx q[2];
rz(-0.70245409) q[2];
sx q[2];
rz(0.89135998) q[2];
rz(-2.4274965) q[3];
sx q[3];
rz(-0.73914206) q[3];
sx q[3];
rz(-2.0785296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63356584) q[0];
sx q[0];
rz(-1.897568) q[0];
sx q[0];
rz(-0.6749534) q[0];
rz(-1.5114816) q[1];
sx q[1];
rz(-2.5210896) q[1];
sx q[1];
rz(0.57704467) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9723608) q[0];
sx q[0];
rz(-1.0857538) q[0];
sx q[0];
rz(0.26974704) q[0];
rz(-pi) q[1];
x q[1];
rz(0.21930947) q[2];
sx q[2];
rz(-2.6990934) q[2];
sx q[2];
rz(-2.2377917) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7865078) q[1];
sx q[1];
rz(-1.6526645) q[1];
sx q[1];
rz(-0.90555377) q[1];
rz(0.8612154) q[3];
sx q[3];
rz(-0.83038721) q[3];
sx q[3];
rz(1.9423241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.453489) q[2];
sx q[2];
rz(-2.0192396) q[2];
sx q[2];
rz(-2.2171059) q[2];
rz(1.9214572) q[3];
sx q[3];
rz(-0.28885463) q[3];
sx q[3];
rz(-3.0815304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7500551) q[0];
sx q[0];
rz(-2.4704762) q[0];
sx q[0];
rz(-1.7975988) q[0];
rz(1.9186107) q[1];
sx q[1];
rz(-1.5593301) q[1];
sx q[1];
rz(1.5498243) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.236178) q[0];
sx q[0];
rz(-1.2981725) q[0];
sx q[0];
rz(0.18698606) q[0];
rz(-1.6835021) q[2];
sx q[2];
rz(-1.1980044) q[2];
sx q[2];
rz(-1.2959803) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2651974) q[1];
sx q[1];
rz(-1.205447) q[1];
sx q[1];
rz(-2.4202034) q[1];
rz(-2.377691) q[3];
sx q[3];
rz(-1.502486) q[3];
sx q[3];
rz(-2.1602283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0595155) q[2];
sx q[2];
rz(-1.0249219) q[2];
sx q[2];
rz(0.22641851) q[2];
rz(-3.1021127) q[3];
sx q[3];
rz(-1.4791146) q[3];
sx q[3];
rz(-1.9421019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74462849) q[0];
sx q[0];
rz(-1.2122943) q[0];
sx q[0];
rz(-2.3126171) q[0];
rz(-0.12116155) q[1];
sx q[1];
rz(-2.1794901) q[1];
sx q[1];
rz(1.6896348) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4273116) q[0];
sx q[0];
rz(-1.4330919) q[0];
sx q[0];
rz(2.4192823) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0214861) q[2];
sx q[2];
rz(-0.7996489) q[2];
sx q[2];
rz(0.19935184) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.02579298) q[1];
sx q[1];
rz(-2.3069972) q[1];
sx q[1];
rz(-2.2195312) q[1];
rz(-pi) q[2];
rz(-2.1690458) q[3];
sx q[3];
rz(-2.498811) q[3];
sx q[3];
rz(1.169426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.58763233) q[2];
sx q[2];
rz(-2.3162737) q[2];
sx q[2];
rz(2.2770605) q[2];
rz(0.65166059) q[3];
sx q[3];
rz(-1.0385907) q[3];
sx q[3];
rz(3.1237176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.637735) q[0];
sx q[0];
rz(-2.2561769) q[0];
sx q[0];
rz(-0.46022415) q[0];
rz(1.3453311) q[1];
sx q[1];
rz(-0.63667744) q[1];
sx q[1];
rz(1.771079) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0287728) q[0];
sx q[0];
rz(-1.1826853) q[0];
sx q[0];
rz(0.32916268) q[0];
rz(1.2865169) q[2];
sx q[2];
rz(-0.95259579) q[2];
sx q[2];
rz(-0.48910352) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7975841) q[1];
sx q[1];
rz(-1.3329778) q[1];
sx q[1];
rz(-1.8348376) q[1];
x q[2];
rz(-3.0881365) q[3];
sx q[3];
rz(-1.7776907) q[3];
sx q[3];
rz(-0.89689909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3687849) q[2];
sx q[2];
rz(-1.3489172) q[2];
sx q[2];
rz(-2.9985912) q[2];
rz(-0.16608206) q[3];
sx q[3];
rz(-2.4487285) q[3];
sx q[3];
rz(-0.79296976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76846692) q[0];
sx q[0];
rz(-1.4906727) q[0];
sx q[0];
rz(1.3937108) q[0];
rz(1.1557747) q[1];
sx q[1];
rz(-1.118569) q[1];
sx q[1];
rz(-2.7973693) q[1];
rz(1.0350107) q[2];
sx q[2];
rz(-2.2314318) q[2];
sx q[2];
rz(-0.4690276) q[2];
rz(-2.0438016) q[3];
sx q[3];
rz(-0.65752959) q[3];
sx q[3];
rz(0.65617954) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
