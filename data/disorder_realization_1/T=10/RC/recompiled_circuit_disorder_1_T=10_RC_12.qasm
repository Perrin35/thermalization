OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.82436615) q[0];
sx q[0];
rz(-1.1146201) q[0];
sx q[0];
rz(-0.00014076509) q[0];
rz(-1.8074942) q[1];
sx q[1];
rz(-0.9642095) q[1];
sx q[1];
rz(-1.1934086) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1228468) q[0];
sx q[0];
rz(-1.393264) q[0];
sx q[0];
rz(1.8621423) q[0];
rz(-pi) q[1];
rz(0.5483746) q[2];
sx q[2];
rz(-1.3142685) q[2];
sx q[2];
rz(-2.246849) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8661583) q[1];
sx q[1];
rz(-2.1992116) q[1];
sx q[1];
rz(-0.98316146) q[1];
x q[2];
rz(2.0334929) q[3];
sx q[3];
rz(-0.24216147) q[3];
sx q[3];
rz(0.36377454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.45941916) q[2];
sx q[2];
rz(-0.023962263) q[2];
sx q[2];
rz(1.2288644) q[2];
rz(1.4131644) q[3];
sx q[3];
rz(-1.1011522) q[3];
sx q[3];
rz(-1.6536973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5380149) q[0];
sx q[0];
rz(-1.6390272) q[0];
sx q[0];
rz(1.0128101) q[0];
rz(0.027659841) q[1];
sx q[1];
rz(-0.67359567) q[1];
sx q[1];
rz(-1.123463) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22966188) q[0];
sx q[0];
rz(-0.46470416) q[0];
sx q[0];
rz(-1.4421411) q[0];
x q[1];
rz(-0.77983071) q[2];
sx q[2];
rz(-1.5260328) q[2];
sx q[2];
rz(-1.1080527) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1883205) q[1];
sx q[1];
rz(-0.76474944) q[1];
sx q[1];
rz(-0.79337593) q[1];
rz(-pi) q[2];
rz(-2.4853737) q[3];
sx q[3];
rz(-1.1074293) q[3];
sx q[3];
rz(-3.0955293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3479487) q[2];
sx q[2];
rz(-1.0898033) q[2];
sx q[2];
rz(0.91903764) q[2];
rz(2.4675026) q[3];
sx q[3];
rz(-2.489311) q[3];
sx q[3];
rz(-1.526171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8640901) q[0];
sx q[0];
rz(-0.16177495) q[0];
sx q[0];
rz(1.8664237) q[0];
rz(2.4480942) q[1];
sx q[1];
rz(-1.2561412) q[1];
sx q[1];
rz(-1.1330053) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95110675) q[0];
sx q[0];
rz(-1.5674942) q[0];
sx q[0];
rz(-1.7130873) q[0];
rz(-pi) q[1];
rz(0.50425597) q[2];
sx q[2];
rz(-2.2849053) q[2];
sx q[2];
rz(-2.2222663) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.36724597) q[1];
sx q[1];
rz(-0.62421747) q[1];
sx q[1];
rz(-1.063785) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3660925) q[3];
sx q[3];
rz(-0.79458487) q[3];
sx q[3];
rz(-0.19898181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.2514078) q[2];
sx q[2];
rz(-0.79139411) q[2];
sx q[2];
rz(1.2934925) q[2];
rz(0.039316468) q[3];
sx q[3];
rz(-1.2189564) q[3];
sx q[3];
rz(-1.8815276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8816198) q[0];
sx q[0];
rz(-0.078475229) q[0];
sx q[0];
rz(1.1608634) q[0];
rz(-2.2456031) q[1];
sx q[1];
rz(-1.7005824) q[1];
sx q[1];
rz(-0.13555759) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7318152) q[0];
sx q[0];
rz(-1.8827794) q[0];
sx q[0];
rz(-0.72809763) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7538257) q[2];
sx q[2];
rz(-2.7515538) q[2];
sx q[2];
rz(2.4412465) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.902066) q[1];
sx q[1];
rz(-0.24589989) q[1];
sx q[1];
rz(-2.2794101) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9837708) q[3];
sx q[3];
rz(-1.2994248) q[3];
sx q[3];
rz(-1.8872758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.23665145) q[2];
sx q[2];
rz(-0.94649482) q[2];
sx q[2];
rz(0.87990749) q[2];
rz(0.044163477) q[3];
sx q[3];
rz(-1.6396089) q[3];
sx q[3];
rz(2.8529609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0376461) q[0];
sx q[0];
rz(-2.7665311) q[0];
sx q[0];
rz(-2.1283545) q[0];
rz(-0.049731072) q[1];
sx q[1];
rz(-0.91369349) q[1];
sx q[1];
rz(2.0577046) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3736553) q[0];
sx q[0];
rz(-2.8345788) q[0];
sx q[0];
rz(-0.92371773) q[0];
x q[1];
rz(-0.25123698) q[2];
sx q[2];
rz(-1.305797) q[2];
sx q[2];
rz(-2.3064409) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0005972) q[1];
sx q[1];
rz(-0.53783572) q[1];
sx q[1];
rz(-1.7829893) q[1];
x q[2];
rz(-0.089133457) q[3];
sx q[3];
rz(-2.6260178) q[3];
sx q[3];
rz(-2.6775132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9087387) q[2];
sx q[2];
rz(-2.8149657) q[2];
sx q[2];
rz(-2.8971635) q[2];
rz(0.43236732) q[3];
sx q[3];
rz(-1.7418539) q[3];
sx q[3];
rz(-2.6385245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2844834) q[0];
sx q[0];
rz(-1.4215707) q[0];
sx q[0];
rz(0.094141468) q[0];
rz(2.969818) q[1];
sx q[1];
rz(-1.1356907) q[1];
sx q[1];
rz(-0.89541268) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5603134) q[0];
sx q[0];
rz(-0.34438294) q[0];
sx q[0];
rz(3.0292105) q[0];
rz(-pi) q[1];
rz(-2.5080639) q[2];
sx q[2];
rz(-1.765536) q[2];
sx q[2];
rz(2.5031236) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0370334) q[1];
sx q[1];
rz(-1.7393141) q[1];
sx q[1];
rz(0.044635459) q[1];
rz(-pi) q[2];
rz(0.0023407188) q[3];
sx q[3];
rz(-1.4961092) q[3];
sx q[3];
rz(2.5403288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.133693) q[2];
sx q[2];
rz(-0.40863016) q[2];
sx q[2];
rz(-0.80319476) q[2];
rz(1.1903654) q[3];
sx q[3];
rz(-1.2322216) q[3];
sx q[3];
rz(2.7289594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.068709277) q[0];
sx q[0];
rz(-0.16462737) q[0];
sx q[0];
rz(0.51914006) q[0];
rz(0.58147645) q[1];
sx q[1];
rz(-1.1053718) q[1];
sx q[1];
rz(1.8849467) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089766895) q[0];
sx q[0];
rz(-1.3194114) q[0];
sx q[0];
rz(-0.041640394) q[0];
x q[1];
rz(1.16876) q[2];
sx q[2];
rz(-2.6425344) q[2];
sx q[2];
rz(-1.4025584) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9650967) q[1];
sx q[1];
rz(-2.8526222) q[1];
sx q[1];
rz(0.22775905) q[1];
x q[2];
rz(0.7092181) q[3];
sx q[3];
rz(-1.8172482) q[3];
sx q[3];
rz(-2.646288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.98823035) q[2];
sx q[2];
rz(-2.1116657) q[2];
sx q[2];
rz(-1.777565) q[2];
rz(0.91056943) q[3];
sx q[3];
rz(-1.986859) q[3];
sx q[3];
rz(-1.6114657) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0664739) q[0];
sx q[0];
rz(-0.56448889) q[0];
sx q[0];
rz(-0.30817729) q[0];
rz(3.0691052) q[1];
sx q[1];
rz(-2.1283573) q[1];
sx q[1];
rz(-0.38696188) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21803074) q[0];
sx q[0];
rz(-2.1523928) q[0];
sx q[0];
rz(0.27115718) q[0];
x q[1];
rz(2.1808124) q[2];
sx q[2];
rz(-0.79723251) q[2];
sx q[2];
rz(-2.0955992) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3124233) q[1];
sx q[1];
rz(-2.2409391) q[1];
sx q[1];
rz(-1.4751242) q[1];
rz(1.3513226) q[3];
sx q[3];
rz(-2.05728) q[3];
sx q[3];
rz(-0.51318491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.62362921) q[2];
sx q[2];
rz(-1.7771746) q[2];
sx q[2];
rz(-2.7015838) q[2];
rz(2.4258339) q[3];
sx q[3];
rz(-1.7093168) q[3];
sx q[3];
rz(-2.0619152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0004262) q[0];
sx q[0];
rz(-0.74581242) q[0];
sx q[0];
rz(1.0986885) q[0];
rz(2.4138342) q[1];
sx q[1];
rz(-0.37574238) q[1];
sx q[1];
rz(-3.0922906) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54287275) q[0];
sx q[0];
rz(-0.96519404) q[0];
sx q[0];
rz(-0.70830317) q[0];
rz(-pi) q[1];
rz(0.88843139) q[2];
sx q[2];
rz(-1.9677791) q[2];
sx q[2];
rz(1.2259442) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.038736131) q[1];
sx q[1];
rz(-2.4604081) q[1];
sx q[1];
rz(0.78050651) q[1];
rz(-0.52600577) q[3];
sx q[3];
rz(-2.3283539) q[3];
sx q[3];
rz(-2.7362652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.07842841) q[2];
sx q[2];
rz(-0.59946632) q[2];
sx q[2];
rz(0.72193974) q[2];
rz(-2.1980964) q[3];
sx q[3];
rz(-0.75073457) q[3];
sx q[3];
rz(0.25434428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14770517) q[0];
sx q[0];
rz(-1.9839956) q[0];
sx q[0];
rz(-1.0797427) q[0];
rz(1.059277) q[1];
sx q[1];
rz(-2.9187027) q[1];
sx q[1];
rz(-1.7396897) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.429075) q[0];
sx q[0];
rz(-1.6970465) q[0];
sx q[0];
rz(-1.7849126) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9548168) q[2];
sx q[2];
rz(-1.597607) q[2];
sx q[2];
rz(1.9865799) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.460443) q[1];
sx q[1];
rz(-0.76421684) q[1];
sx q[1];
rz(-1.1230254) q[1];
rz(0.28016443) q[3];
sx q[3];
rz(-2.4912842) q[3];
sx q[3];
rz(-3.0349817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.4828651) q[2];
sx q[2];
rz(-1.7752703) q[2];
sx q[2];
rz(1.6213017) q[2];
rz(0.55082095) q[3];
sx q[3];
rz(-2.3362624) q[3];
sx q[3];
rz(-2.4441392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14810066) q[0];
sx q[0];
rz(-1.3052595) q[0];
sx q[0];
rz(-1.530151) q[0];
rz(-2.2254754) q[1];
sx q[1];
rz(-0.59090186) q[1];
sx q[1];
rz(-0.59060243) q[1];
rz(2.312071) q[2];
sx q[2];
rz(-0.98486949) q[2];
sx q[2];
rz(-2.9688901) q[2];
rz(-2.0879073) q[3];
sx q[3];
rz(-1.5552945) q[3];
sx q[3];
rz(-1.8154715) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];