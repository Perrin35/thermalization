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
rz(-0.63602716) q[1];
sx q[1];
rz(5.7732328) q[1];
sx q[1];
rz(10.049811) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63955414) q[0];
sx q[0];
rz(-0.14766492) q[0];
sx q[0];
rz(0.30252151) q[0];
rz(-1.6051859) q[2];
sx q[2];
rz(-2.6011438) q[2];
sx q[2];
rz(2.2480223) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.6909862) q[1];
sx q[1];
rz(-2.4391973) q[1];
sx q[1];
rz(2.2390963) q[1];
rz(-pi) q[2];
rz(0.19842438) q[3];
sx q[3];
rz(-0.82272595) q[3];
sx q[3];
rz(-0.89694512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.92224377) q[2];
sx q[2];
rz(-0.85180989) q[2];
sx q[2];
rz(-0.40199486) q[2];
rz(0.829202) q[3];
sx q[3];
rz(-1.3196557) q[3];
sx q[3];
rz(-1.3623665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40792337) q[0];
sx q[0];
rz(-2.5897554) q[0];
sx q[0];
rz(-2.0358987) q[0];
rz(0.78877527) q[1];
sx q[1];
rz(-2.5194247) q[1];
sx q[1];
rz(2.6355991) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.087990847) q[0];
sx q[0];
rz(-1.3357541) q[0];
sx q[0];
rz(-1.4925692) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.77941915) q[2];
sx q[2];
rz(-0.96923087) q[2];
sx q[2];
rz(2.0815408) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2179394) q[1];
sx q[1];
rz(-0.42227498) q[1];
sx q[1];
rz(2.1417066) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6384505) q[3];
sx q[3];
rz(-0.50853679) q[3];
sx q[3];
rz(-3.0340305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.904423) q[2];
sx q[2];
rz(-2.0776896) q[2];
sx q[2];
rz(-2.2696631) q[2];
rz(2.0875841) q[3];
sx q[3];
rz(-1.0796615) q[3];
sx q[3];
rz(-1.7006629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.089040861) q[0];
sx q[0];
rz(-2.4854923) q[0];
sx q[0];
rz(-1.22714) q[0];
rz(0.50651208) q[1];
sx q[1];
rz(-2.3498693) q[1];
sx q[1];
rz(0.92481771) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8657579) q[0];
sx q[0];
rz(-1.6730437) q[0];
sx q[0];
rz(-1.7195549) q[0];
rz(-2.6854555) q[2];
sx q[2];
rz(-1.6326687) q[2];
sx q[2];
rz(-1.8707239) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1113092) q[1];
sx q[1];
rz(-0.96007512) q[1];
sx q[1];
rz(-0.3719117) q[1];
x q[2];
rz(-1.1919349) q[3];
sx q[3];
rz(-2.1964245) q[3];
sx q[3];
rz(1.4257087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0623124) q[2];
sx q[2];
rz(-1.5996876) q[2];
sx q[2];
rz(-0.49528948) q[2];
rz(1.371572) q[3];
sx q[3];
rz(-0.74426952) q[3];
sx q[3];
rz(-1.8194958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0168734) q[0];
sx q[0];
rz(-1.9833516) q[0];
sx q[0];
rz(0.88791263) q[0];
rz(1.3814231) q[1];
sx q[1];
rz(-1.0180749) q[1];
sx q[1];
rz(-0.51753128) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84758112) q[0];
sx q[0];
rz(-1.309762) q[0];
sx q[0];
rz(3.1270967) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5597455) q[2];
sx q[2];
rz(-0.80815017) q[2];
sx q[2];
rz(0.86607546) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1932839) q[1];
sx q[1];
rz(-1.6689321) q[1];
sx q[1];
rz(-0.0098642439) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1903619) q[3];
sx q[3];
rz(-1.1305792) q[3];
sx q[3];
rz(-2.7444283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0722644) q[2];
sx q[2];
rz(-1.1257409) q[2];
sx q[2];
rz(-0.98108712) q[2];
rz(2.6632994) q[3];
sx q[3];
rz(-0.35224733) q[3];
sx q[3];
rz(0.52811629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1961871) q[0];
sx q[0];
rz(-1.4237175) q[0];
sx q[0];
rz(-0.032935306) q[0];
rz(1.6163274) q[1];
sx q[1];
rz(-2.567629) q[1];
sx q[1];
rz(-0.93564916) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52152741) q[0];
sx q[0];
rz(-0.47019401) q[0];
sx q[0];
rz(3.1121701) q[0];
rz(-pi) q[1];
rz(2.0528556) q[2];
sx q[2];
rz(-2.3080359) q[2];
sx q[2];
rz(-0.6635467) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.050154479) q[1];
sx q[1];
rz(-0.77934105) q[1];
sx q[1];
rz(3.0528751) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8132951) q[3];
sx q[3];
rz(-2.1418013) q[3];
sx q[3];
rz(1.5703438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4719438) q[2];
sx q[2];
rz(-2.7574597) q[2];
sx q[2];
rz(-1.498339) q[2];
rz(-2.5390427) q[3];
sx q[3];
rz(-0.82585255) q[3];
sx q[3];
rz(1.92441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54685408) q[0];
sx q[0];
rz(-0.29404077) q[0];
sx q[0];
rz(-0.038851693) q[0];
rz(-0.36000577) q[1];
sx q[1];
rz(-1.729676) q[1];
sx q[1];
rz(0.14883277) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.50406) q[0];
sx q[0];
rz(-1.6967745) q[0];
sx q[0];
rz(-0.12772638) q[0];
x q[1];
rz(-2.932933) q[2];
sx q[2];
rz(-2.5620915) q[2];
sx q[2];
rz(-0.3267056) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.57362) q[1];
sx q[1];
rz(-0.96155969) q[1];
sx q[1];
rz(0.9602169) q[1];
x q[2];
rz(-0.43977387) q[3];
sx q[3];
rz(-2.5086702) q[3];
sx q[3];
rz(-1.9683815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0570602) q[2];
sx q[2];
rz(-2.3134505) q[2];
sx q[2];
rz(-0.35489902) q[2];
rz(-1.0585632) q[3];
sx q[3];
rz(-2.1954229) q[3];
sx q[3];
rz(-0.53275776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8232089) q[0];
sx q[0];
rz(-1.4223149) q[0];
sx q[0];
rz(0.80867714) q[0];
rz(-3.0478364) q[1];
sx q[1];
rz(-0.84051991) q[1];
sx q[1];
rz(-0.41059986) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79390271) q[0];
sx q[0];
rz(-1.1979894) q[0];
sx q[0];
rz(1.1469141) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8062079) q[2];
sx q[2];
rz(-0.68143564) q[2];
sx q[2];
rz(-0.49668703) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6462999) q[1];
sx q[1];
rz(-1.22806) q[1];
sx q[1];
rz(0.41151026) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4491354) q[3];
sx q[3];
rz(-2.231519) q[3];
sx q[3];
rz(2.4796951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.14313301) q[2];
sx q[2];
rz(-2.748558) q[2];
sx q[2];
rz(0.03579363) q[2];
rz(-1.834747) q[3];
sx q[3];
rz(-1.9948317) q[3];
sx q[3];
rz(0.80865639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3150113) q[0];
sx q[0];
rz(-0.97493521) q[0];
sx q[0];
rz(0.63631979) q[0];
rz(2.4782205) q[1];
sx q[1];
rz(-2.0320818) q[1];
sx q[1];
rz(-2.5136307) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0215752) q[0];
sx q[0];
rz(-1.4940959) q[0];
sx q[0];
rz(-1.5995366) q[0];
rz(-pi) q[1];
x q[1];
rz(0.2070932) q[2];
sx q[2];
rz(-2.7367055) q[2];
sx q[2];
rz(-1.5804497) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0934712) q[1];
sx q[1];
rz(-1.1560165) q[1];
sx q[1];
rz(-2.5268447) q[1];
rz(-pi) q[2];
rz(2.9541438) q[3];
sx q[3];
rz(-1.6050242) q[3];
sx q[3];
rz(1.6899861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5593354) q[2];
sx q[2];
rz(-3.0310243) q[2];
sx q[2];
rz(-1.4165233) q[2];
rz(2.0420117) q[3];
sx q[3];
rz(-2.1018335) q[3];
sx q[3];
rz(-0.85272461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.002554) q[0];
sx q[0];
rz(-2.250493) q[0];
sx q[0];
rz(-0.6148327) q[0];
rz(-0.50827208) q[1];
sx q[1];
rz(-2.1891687) q[1];
sx q[1];
rz(-0.27613786) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50867453) q[0];
sx q[0];
rz(-2.6946697) q[0];
sx q[0];
rz(1.267295) q[0];
x q[1];
rz(2.1008087) q[2];
sx q[2];
rz(-1.4333087) q[2];
sx q[2];
rz(1.9078209) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.74127585) q[1];
sx q[1];
rz(-2.0243249) q[1];
sx q[1];
rz(-0.67429115) q[1];
x q[2];
rz(-2.7998447) q[3];
sx q[3];
rz(-0.93999388) q[3];
sx q[3];
rz(-2.3113826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.95163661) q[2];
sx q[2];
rz(-1.0389682) q[2];
sx q[2];
rz(-3.0432126) q[2];
rz(1.240587) q[3];
sx q[3];
rz(-2.0799347) q[3];
sx q[3];
rz(-0.043717535) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8561309) q[0];
sx q[0];
rz(-1.0482482) q[0];
sx q[0];
rz(-0.35828006) q[0];
rz(1.0048535) q[1];
sx q[1];
rz(-0.88468164) q[1];
sx q[1];
rz(-2.7994432) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.005363883) q[0];
sx q[0];
rz(-2.507302) q[0];
sx q[0];
rz(-2.7691288) q[0];
rz(-pi) q[1];
rz(1.0295632) q[2];
sx q[2];
rz(-2.4839253) q[2];
sx q[2];
rz(0.56305158) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.725199) q[1];
sx q[1];
rz(-2.1402845) q[1];
sx q[1];
rz(1.4432905) q[1];
rz(0.10327424) q[3];
sx q[3];
rz(-1.2359796) q[3];
sx q[3];
rz(-1.1371727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7973914) q[2];
sx q[2];
rz(-2.6579393) q[2];
sx q[2];
rz(0.3717711) q[2];
rz(-0.19112912) q[3];
sx q[3];
rz(-1.976795) q[3];
sx q[3];
rz(-0.42310664) q[3];
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
rz(0.38901781) q[0];
sx q[0];
rz(-0.49674635) q[0];
sx q[0];
rz(2.4865271) q[0];
rz(2.7816506) q[1];
sx q[1];
rz(-1.5690201) q[1];
sx q[1];
rz(1.5108861) q[1];
rz(0.062535738) q[2];
sx q[2];
rz(-1.8725431) q[2];
sx q[2];
rz(0.098081577) q[2];
rz(0.015040811) q[3];
sx q[3];
rz(-0.81285601) q[3];
sx q[3];
rz(2.6668919) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
