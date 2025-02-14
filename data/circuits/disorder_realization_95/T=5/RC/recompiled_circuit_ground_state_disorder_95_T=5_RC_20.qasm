OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.65741444) q[0];
sx q[0];
rz(1.2793469) q[0];
sx q[0];
rz(10.068324) q[0];
rz(0.14856385) q[1];
sx q[1];
rz(3.7879877) q[1];
sx q[1];
rz(11.990379) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9849562) q[0];
sx q[0];
rz(-2.505777) q[0];
sx q[0];
rz(1.030907) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9690016) q[2];
sx q[2];
rz(-1.2786037) q[2];
sx q[2];
rz(-1.9167231) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.483584) q[1];
sx q[1];
rz(-1.5414667) q[1];
sx q[1];
rz(1.6701677) q[1];
rz(-2.3131274) q[3];
sx q[3];
rz(-1.7196619) q[3];
sx q[3];
rz(0.016691301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5129471) q[2];
sx q[2];
rz(-0.60428667) q[2];
sx q[2];
rz(-2.4224572) q[2];
rz(-0.61270815) q[3];
sx q[3];
rz(-0.49379525) q[3];
sx q[3];
rz(-1.4601716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4145819) q[0];
sx q[0];
rz(-2.5717323) q[0];
sx q[0];
rz(1.7852831) q[0];
rz(2.5510229) q[1];
sx q[1];
rz(-1.5868203) q[1];
sx q[1];
rz(-0.86033598) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7588494) q[0];
sx q[0];
rz(-0.68419391) q[0];
sx q[0];
rz(1.280715) q[0];
rz(0.55540076) q[2];
sx q[2];
rz(-1.4216219) q[2];
sx q[2];
rz(-2.5859216) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.115827) q[1];
sx q[1];
rz(-1.888926) q[1];
sx q[1];
rz(0.33237388) q[1];
rz(-pi) q[2];
x q[2];
rz(0.72655118) q[3];
sx q[3];
rz(-1.6457215) q[3];
sx q[3];
rz(2.5963432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7343973) q[2];
sx q[2];
rz(-1.1043786) q[2];
sx q[2];
rz(-2.8779136) q[2];
rz(0.10144357) q[3];
sx q[3];
rz(-1.0976617) q[3];
sx q[3];
rz(1.6775848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9145255) q[0];
sx q[0];
rz(-0.83221808) q[0];
sx q[0];
rz(1.2365923) q[0];
rz(-0.49691686) q[1];
sx q[1];
rz(-2.2221815) q[1];
sx q[1];
rz(-0.20981728) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.041417) q[0];
sx q[0];
rz(-0.13209535) q[0];
sx q[0];
rz(0.14762525) q[0];
rz(-1.1205925) q[2];
sx q[2];
rz(-1.6944024) q[2];
sx q[2];
rz(-0.093401366) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1637437) q[1];
sx q[1];
rz(-1.6075378) q[1];
sx q[1];
rz(-2.0746465) q[1];
rz(2.3043892) q[3];
sx q[3];
rz(-1.4542237) q[3];
sx q[3];
rz(-1.0511412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1257552) q[2];
sx q[2];
rz(-1.1344942) q[2];
sx q[2];
rz(2.8718359) q[2];
rz(2.6848327) q[3];
sx q[3];
rz(-0.41958198) q[3];
sx q[3];
rz(-2.8381798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1064827) q[0];
sx q[0];
rz(-2.7484317) q[0];
sx q[0];
rz(0.092057236) q[0];
rz(2.5606142) q[1];
sx q[1];
rz(-0.61994225) q[1];
sx q[1];
rz(2.7041025) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96707805) q[0];
sx q[0];
rz(-1.0096435) q[0];
sx q[0];
rz(-2.9516417) q[0];
rz(-1.7873437) q[2];
sx q[2];
rz(-1.6833954) q[2];
sx q[2];
rz(-1.4527904) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6981004) q[1];
sx q[1];
rz(-2.6671706) q[1];
sx q[1];
rz(2.1441133) q[1];
rz(-2.8169223) q[3];
sx q[3];
rz(-2.2483147) q[3];
sx q[3];
rz(3.0538525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.94365326) q[2];
sx q[2];
rz(-1.7655756) q[2];
sx q[2];
rz(0.1680689) q[2];
rz(-0.69593143) q[3];
sx q[3];
rz(-1.9441425) q[3];
sx q[3];
rz(-1.7456938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54774061) q[0];
sx q[0];
rz(-0.20741367) q[0];
sx q[0];
rz(0.70163027) q[0];
rz(-2.592109) q[1];
sx q[1];
rz(-1.0797078) q[1];
sx q[1];
rz(0.93542498) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.815274) q[0];
sx q[0];
rz(-1.9978412) q[0];
sx q[0];
rz(0.27384211) q[0];
x q[1];
rz(-0.1232202) q[2];
sx q[2];
rz(-1.5180523) q[2];
sx q[2];
rz(-3.0739917) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0282057) q[1];
sx q[1];
rz(-1.0019687) q[1];
sx q[1];
rz(-2.6527106) q[1];
x q[2];
rz(-1.6356556) q[3];
sx q[3];
rz(-1.6422513) q[3];
sx q[3];
rz(0.91141975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6286185) q[2];
sx q[2];
rz(-1.3209891) q[2];
sx q[2];
rz(2.4690907) q[2];
rz(2.3682112) q[3];
sx q[3];
rz(-2.0956109) q[3];
sx q[3];
rz(-0.33236233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8124939) q[0];
sx q[0];
rz(-2.2342873) q[0];
sx q[0];
rz(-1.6727653) q[0];
rz(-0.86917296) q[1];
sx q[1];
rz(-2.4425127) q[1];
sx q[1];
rz(-0.3259784) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3968626) q[0];
sx q[0];
rz(-0.32856634) q[0];
sx q[0];
rz(-2.7129125) q[0];
x q[1];
rz(2.5766854) q[2];
sx q[2];
rz(-2.0015026) q[2];
sx q[2];
rz(0.6503833) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1010998) q[1];
sx q[1];
rz(-2.0475351) q[1];
sx q[1];
rz(-0.65495305) q[1];
rz(1.7160837) q[3];
sx q[3];
rz(-2.709148) q[3];
sx q[3];
rz(-2.7029981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.088923067) q[2];
sx q[2];
rz(-1.576705) q[2];
sx q[2];
rz(-0.79724533) q[2];
rz(3.0439324) q[3];
sx q[3];
rz(-0.74334136) q[3];
sx q[3];
rz(1.2231479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(0.866211) q[0];
sx q[0];
rz(-2.1480063) q[0];
sx q[0];
rz(-2.1589101) q[0];
rz(-1.4879701) q[1];
sx q[1];
rz(-2.5762312) q[1];
sx q[1];
rz(2.8178094) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57162962) q[0];
sx q[0];
rz(-1.6545719) q[0];
sx q[0];
rz(-0.28599583) q[0];
rz(-pi) q[1];
rz(-2.1828142) q[2];
sx q[2];
rz(-1.3347776) q[2];
sx q[2];
rz(2.4356206) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.946879) q[1];
sx q[1];
rz(-0.66821287) q[1];
sx q[1];
rz(-2.9577423) q[1];
rz(-1.9110095) q[3];
sx q[3];
rz(-1.0939301) q[3];
sx q[3];
rz(-0.25097706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4299778) q[2];
sx q[2];
rz(-0.51402503) q[2];
sx q[2];
rz(0.035710486) q[2];
rz(-2.4014373) q[3];
sx q[3];
rz(-1.6868846) q[3];
sx q[3];
rz(-1.9945701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3586054) q[0];
sx q[0];
rz(-2.7645223) q[0];
sx q[0];
rz(0.20088917) q[0];
rz(-0.56328493) q[1];
sx q[1];
rz(-1.3689684) q[1];
sx q[1];
rz(-0.39931452) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5168229) q[0];
sx q[0];
rz(-1.4218569) q[0];
sx q[0];
rz(-0.14601222) q[0];
x q[1];
rz(-1.6496193) q[2];
sx q[2];
rz(-0.48187253) q[2];
sx q[2];
rz(-2.3367867) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6714564) q[1];
sx q[1];
rz(-2.5920007) q[1];
sx q[1];
rz(-0.84322815) q[1];
rz(-pi) q[2];
rz(-2.3763555) q[3];
sx q[3];
rz(-1.8730193) q[3];
sx q[3];
rz(1.1208432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1839972) q[2];
sx q[2];
rz(-0.60773578) q[2];
sx q[2];
rz(-3.0228534) q[2];
rz(2.8710098) q[3];
sx q[3];
rz(-0.99004531) q[3];
sx q[3];
rz(2.9873007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1110693) q[0];
sx q[0];
rz(-1.1708165) q[0];
sx q[0];
rz(0.4554553) q[0];
rz(0.2332553) q[1];
sx q[1];
rz(-2.3979954) q[1];
sx q[1];
rz(3.1366248) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7457897) q[0];
sx q[0];
rz(-0.088287778) q[0];
sx q[0];
rz(2.7448065) q[0];
rz(-1.8894779) q[2];
sx q[2];
rz(-1.9405085) q[2];
sx q[2];
rz(1.2546509) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.78902728) q[1];
sx q[1];
rz(-2.0754756) q[1];
sx q[1];
rz(-0.53225174) q[1];
rz(-pi) q[2];
rz(1.6065323) q[3];
sx q[3];
rz(-2.7155345) q[3];
sx q[3];
rz(-0.96340513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7838955) q[2];
sx q[2];
rz(-1.8507439) q[2];
sx q[2];
rz(2.2515187) q[2];
rz(-0.87738532) q[3];
sx q[3];
rz(-0.93534094) q[3];
sx q[3];
rz(0.64524159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9599065) q[0];
sx q[0];
rz(-0.077597685) q[0];
sx q[0];
rz(2.087387) q[0];
rz(2.8554845) q[1];
sx q[1];
rz(-2.095463) q[1];
sx q[1];
rz(1.8792763) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0204476) q[0];
sx q[0];
rz(-0.95709267) q[0];
sx q[0];
rz(-1.2200939) q[0];
rz(1.3081864) q[2];
sx q[2];
rz(-1.76909) q[2];
sx q[2];
rz(-1.8526371) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9156277) q[1];
sx q[1];
rz(-2.0061992) q[1];
sx q[1];
rz(1.7724724) q[1];
rz(-pi) q[2];
rz(1.0159053) q[3];
sx q[3];
rz(-1.4994748) q[3];
sx q[3];
rz(-1.2268905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.012764843) q[2];
sx q[2];
rz(-0.78170693) q[2];
sx q[2];
rz(-0.26415602) q[2];
rz(-2.9014897) q[3];
sx q[3];
rz(-1.0485579) q[3];
sx q[3];
rz(-0.31606328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8964597) q[0];
sx q[0];
rz(-0.75440732) q[0];
sx q[0];
rz(-0.9040133) q[0];
rz(-0.28884197) q[1];
sx q[1];
rz(-1.387351) q[1];
sx q[1];
rz(3.0659061) q[1];
rz(2.8223128) q[2];
sx q[2];
rz(-1.3360595) q[2];
sx q[2];
rz(1.3315249) q[2];
rz(-3.0700923) q[3];
sx q[3];
rz(-1.6115474) q[3];
sx q[3];
rz(-1.1781884) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
