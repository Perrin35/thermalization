OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.3990134) q[0];
sx q[0];
rz(1.9785545) q[0];
sx q[0];
rz(7.4217441) q[0];
rz(2.5055655) q[1];
sx q[1];
rz(-2.6316402) q[1];
sx q[1];
rz(-0.62503254) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2306657) q[0];
sx q[0];
rz(-1.6146442) q[0];
sx q[0];
rz(-3.0005428) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1109842) q[2];
sx q[2];
rz(-1.5884879) q[2];
sx q[2];
rz(2.4938581) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.33746943) q[1];
sx q[1];
rz(-1.9826681) q[1];
sx q[1];
rz(0.98442529) q[1];
rz(-2.3287541) q[3];
sx q[3];
rz(-1.715797) q[3];
sx q[3];
rz(-0.53792153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.92224377) q[2];
sx q[2];
rz(-2.2897828) q[2];
sx q[2];
rz(0.40199486) q[2];
rz(-0.829202) q[3];
sx q[3];
rz(-1.821937) q[3];
sx q[3];
rz(1.7792262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7336693) q[0];
sx q[0];
rz(-2.5897554) q[0];
sx q[0];
rz(-2.0358987) q[0];
rz(2.3528174) q[1];
sx q[1];
rz(-0.62216798) q[1];
sx q[1];
rz(2.6355991) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9049022) q[0];
sx q[0];
rz(-0.2474898) q[0];
sx q[0];
rz(-0.31544073) q[0];
rz(-pi) q[1];
rz(2.3621735) q[2];
sx q[2];
rz(-0.96923087) q[2];
sx q[2];
rz(2.0815408) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0243903) q[1];
sx q[1];
rz(-1.7941231) q[1];
sx q[1];
rz(-1.2093556) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6873029) q[3];
sx q[3];
rz(-1.807782) q[3];
sx q[3];
rz(-1.2302356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.904423) q[2];
sx q[2];
rz(-2.0776896) q[2];
sx q[2];
rz(-0.87192956) q[2];
rz(1.0540086) q[3];
sx q[3];
rz(-2.0619312) q[3];
sx q[3];
rz(1.4409298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0525518) q[0];
sx q[0];
rz(-2.4854923) q[0];
sx q[0];
rz(1.22714) q[0];
rz(-2.6350806) q[1];
sx q[1];
rz(-0.7917234) q[1];
sx q[1];
rz(2.2167749) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4445405) q[0];
sx q[0];
rz(-2.9612975) q[0];
sx q[0];
rz(0.96526115) q[0];
rz(-pi) q[1];
rz(-0.13972762) q[2];
sx q[2];
rz(-0.46002166) q[2];
sx q[2];
rz(2.7163028) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1113092) q[1];
sx q[1];
rz(-0.96007512) q[1];
sx q[1];
rz(-2.769681) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9496578) q[3];
sx q[3];
rz(-2.1964245) q[3];
sx q[3];
rz(1.4257087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0792803) q[2];
sx q[2];
rz(-1.5996876) q[2];
sx q[2];
rz(2.6463032) q[2];
rz(1.371572) q[3];
sx q[3];
rz(-2.3973231) q[3];
sx q[3];
rz(-1.3220968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0168734) q[0];
sx q[0];
rz(-1.158241) q[0];
sx q[0];
rz(2.25368) q[0];
rz(1.7601695) q[1];
sx q[1];
rz(-2.1235178) q[1];
sx q[1];
rz(-0.51753128) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.237898) q[0];
sx q[0];
rz(-0.26142739) q[0];
sx q[0];
rz(-1.625007) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58184718) q[2];
sx q[2];
rz(-0.80815017) q[2];
sx q[2];
rz(-0.86607546) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1932839) q[1];
sx q[1];
rz(-1.4726606) q[1];
sx q[1];
rz(-0.0098642439) q[1];
x q[2];
rz(-1.1903619) q[3];
sx q[3];
rz(-2.0110134) q[3];
sx q[3];
rz(0.39716431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0722644) q[2];
sx q[2];
rz(-2.0158517) q[2];
sx q[2];
rz(2.1605055) q[2];
rz(-0.4782933) q[3];
sx q[3];
rz(-2.7893453) q[3];
sx q[3];
rz(2.6134764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1961871) q[0];
sx q[0];
rz(-1.4237175) q[0];
sx q[0];
rz(-3.1086573) q[0];
rz(-1.6163274) q[1];
sx q[1];
rz(-0.5739637) q[1];
sx q[1];
rz(2.2059435) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0230377) q[0];
sx q[0];
rz(-1.5574677) q[0];
sx q[0];
rz(-2.6715735) q[0];
rz(-pi) q[1];
rz(-1.0887371) q[2];
sx q[2];
rz(-2.3080359) q[2];
sx q[2];
rz(2.478046) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5577573) q[1];
sx q[1];
rz(-1.6331065) q[1];
sx q[1];
rz(0.77737096) q[1];
rz(0.35772459) q[3];
sx q[3];
rz(-2.5265363) q[3];
sx q[3];
rz(-2.0004724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.66964883) q[2];
sx q[2];
rz(-2.7574597) q[2];
sx q[2];
rz(1.498339) q[2];
rz(-2.5390427) q[3];
sx q[3];
rz(-0.82585255) q[3];
sx q[3];
rz(1.92441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54685408) q[0];
sx q[0];
rz(-2.8475519) q[0];
sx q[0];
rz(-0.038851693) q[0];
rz(2.7815869) q[1];
sx q[1];
rz(-1.729676) q[1];
sx q[1];
rz(-2.9927599) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63753269) q[0];
sx q[0];
rz(-1.6967745) q[0];
sx q[0];
rz(3.0138663) q[0];
x q[1];
rz(-1.4360481) q[2];
sx q[2];
rz(-2.1361669) q[2];
sx q[2];
rz(-3.0627405) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.31735709) q[1];
sx q[1];
rz(-2.3075793) q[1];
sx q[1];
rz(-2.4537817) q[1];
rz(0.58601309) q[3];
sx q[3];
rz(-1.3162321) q[3];
sx q[3];
rz(2.3814122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0570602) q[2];
sx q[2];
rz(-0.82814211) q[2];
sx q[2];
rz(0.35489902) q[2];
rz(2.0830294) q[3];
sx q[3];
rz(-0.94616977) q[3];
sx q[3];
rz(-2.6088349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8232089) q[0];
sx q[0];
rz(-1.7192778) q[0];
sx q[0];
rz(0.80867714) q[0];
rz(-3.0478364) q[1];
sx q[1];
rz(-2.3010727) q[1];
sx q[1];
rz(-2.7309928) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5275972) q[0];
sx q[0];
rz(-1.1777012) q[0];
sx q[0];
rz(-0.40531203) q[0];
x q[1];
rz(1.8316531) q[2];
sx q[2];
rz(-2.2078288) q[2];
sx q[2];
rz(-0.9185248) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4098874) q[1];
sx q[1];
rz(-0.52919878) q[1];
sx q[1];
rz(-0.7284109) q[1];
x q[2];
rz(-2.477272) q[3];
sx q[3];
rz(-1.4748286) q[3];
sx q[3];
rz(0.98379083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.14313301) q[2];
sx q[2];
rz(-2.748558) q[2];
sx q[2];
rz(0.03579363) q[2];
rz(-1.834747) q[3];
sx q[3];
rz(-1.9948317) q[3];
sx q[3];
rz(-2.3329363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3150113) q[0];
sx q[0];
rz(-0.97493521) q[0];
sx q[0];
rz(-0.63631979) q[0];
rz(-0.66337216) q[1];
sx q[1];
rz(-1.1095108) q[1];
sx q[1];
rz(2.5136307) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7610891) q[0];
sx q[0];
rz(-3.0596943) q[0];
sx q[0];
rz(2.7837672) q[0];
x q[1];
rz(2.7444674) q[2];
sx q[2];
rz(-1.6518804) q[2];
sx q[2];
rz(-0.18111595) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.20013576) q[1];
sx q[1];
rz(-2.1268599) q[1];
sx q[1];
rz(-1.0764313) q[1];
rz(-pi) q[2];
rz(2.9598755) q[3];
sx q[3];
rz(-0.19051208) q[3];
sx q[3];
rz(-0.2976892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.58225727) q[2];
sx q[2];
rz(-3.0310243) q[2];
sx q[2];
rz(1.7250693) q[2];
rz(-2.0420117) q[3];
sx q[3];
rz(-1.0397592) q[3];
sx q[3];
rz(-0.85272461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.002554) q[0];
sx q[0];
rz(-0.89109963) q[0];
sx q[0];
rz(0.6148327) q[0];
rz(-0.50827208) q[1];
sx q[1];
rz(-0.95242396) q[1];
sx q[1];
rz(0.27613786) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3373703) q[0];
sx q[0];
rz(-1.7003248) q[0];
sx q[0];
rz(1.1418377) q[0];
x q[1];
rz(-1.040784) q[2];
sx q[2];
rz(-1.4333087) q[2];
sx q[2];
rz(1.9078209) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.74127585) q[1];
sx q[1];
rz(-1.1172677) q[1];
sx q[1];
rz(-0.67429115) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2302141) q[3];
sx q[3];
rz(-1.8448534) q[3];
sx q[3];
rz(2.1942133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.189956) q[2];
sx q[2];
rz(-1.0389682) q[2];
sx q[2];
rz(-3.0432126) q[2];
rz(1.240587) q[3];
sx q[3];
rz(-2.0799347) q[3];
sx q[3];
rz(3.0978751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28546178) q[0];
sx q[0];
rz(-2.0933445) q[0];
sx q[0];
rz(0.35828006) q[0];
rz(1.0048535) q[1];
sx q[1];
rz(-2.256911) q[1];
sx q[1];
rz(2.7994432) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.005363883) q[0];
sx q[0];
rz(-2.507302) q[0];
sx q[0];
rz(-2.7691288) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1120295) q[2];
sx q[2];
rz(-2.4839253) q[2];
sx q[2];
rz(-0.56305158) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4163936) q[1];
sx q[1];
rz(-2.1402845) q[1];
sx q[1];
rz(1.6983022) q[1];
rz(-pi) q[2];
rz(3.0383184) q[3];
sx q[3];
rz(-1.905613) q[3];
sx q[3];
rz(2.0044199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3442012) q[2];
sx q[2];
rz(-2.6579393) q[2];
sx q[2];
rz(-0.3717711) q[2];
rz(2.9504635) q[3];
sx q[3];
rz(-1.1647977) q[3];
sx q[3];
rz(-2.718486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7525748) q[0];
sx q[0];
rz(-0.49674635) q[0];
sx q[0];
rz(2.4865271) q[0];
rz(0.3599421) q[1];
sx q[1];
rz(-1.5725726) q[1];
sx q[1];
rz(-1.6307065) q[1];
rz(-1.2684939) q[2];
sx q[2];
rz(-1.5110895) q[2];
sx q[2];
rz(-1.4541077) q[2];
rz(3.1265518) q[3];
sx q[3];
rz(-2.3287366) q[3];
sx q[3];
rz(-0.47470075) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
