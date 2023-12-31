OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.26602715) q[0];
sx q[0];
rz(2.6063483) q[0];
sx q[0];
rz(8.6707414) q[0];
rz(-2.3513878) q[1];
sx q[1];
rz(-1.9145929) q[1];
sx q[1];
rz(1.1608646) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95522674) q[0];
sx q[0];
rz(-0.47954924) q[0];
sx q[0];
rz(0.064155302) q[0];
x q[1];
rz(-1.8765175) q[2];
sx q[2];
rz(-0.49760488) q[2];
sx q[2];
rz(1.9203609) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7272211) q[1];
sx q[1];
rz(-1.1611658) q[1];
sx q[1];
rz(-2.1280671) q[1];
rz(-pi) q[2];
rz(-1.7369032) q[3];
sx q[3];
rz(-1.7475834) q[3];
sx q[3];
rz(-3.0151031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6291818) q[2];
sx q[2];
rz(-1.2966195) q[2];
sx q[2];
rz(-0.66649246) q[2];
rz(-0.50981057) q[3];
sx q[3];
rz(-1.9911659) q[3];
sx q[3];
rz(-1.2734909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3614685) q[0];
sx q[0];
rz(-1.7391917) q[0];
sx q[0];
rz(-2.0289039) q[0];
rz(0.15377046) q[1];
sx q[1];
rz(-2.230481) q[1];
sx q[1];
rz(1.8033093) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1424944) q[0];
sx q[0];
rz(-0.72744838) q[0];
sx q[0];
rz(-2.8508458) q[0];
rz(-pi) q[1];
rz(1.0988702) q[2];
sx q[2];
rz(-0.79052351) q[2];
sx q[2];
rz(-2.3015442) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.2087304) q[1];
sx q[1];
rz(-1.6881697) q[1];
sx q[1];
rz(-0.37456234) q[1];
rz(0.58740939) q[3];
sx q[3];
rz(-1.4360768) q[3];
sx q[3];
rz(2.2692674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.082836941) q[2];
sx q[2];
rz(-0.78410167) q[2];
sx q[2];
rz(-1.178297) q[2];
rz(0.96238771) q[3];
sx q[3];
rz(-1.0932086) q[3];
sx q[3];
rz(-2.4760831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96175471) q[0];
sx q[0];
rz(-2.6340155) q[0];
sx q[0];
rz(1.5270365) q[0];
rz(2.4987192) q[1];
sx q[1];
rz(-1.0765272) q[1];
sx q[1];
rz(-0.33338526) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3231455) q[0];
sx q[0];
rz(-2.4172757) q[0];
sx q[0];
rz(-1.2562654) q[0];
x q[1];
rz(0.1336735) q[2];
sx q[2];
rz(-1.1867656) q[2];
sx q[2];
rz(3.0093699) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.15409878) q[1];
sx q[1];
rz(-1.0990267) q[1];
sx q[1];
rz(-3.0497453) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7942579) q[3];
sx q[3];
rz(-2.1922605) q[3];
sx q[3];
rz(0.37924757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.12144111) q[2];
sx q[2];
rz(-1.8276428) q[2];
sx q[2];
rz(0.80231673) q[2];
rz(0.24117593) q[3];
sx q[3];
rz(-2.4427588) q[3];
sx q[3];
rz(3.0407217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0872831) q[0];
sx q[0];
rz(-1.5232975) q[0];
sx q[0];
rz(0.56418443) q[0];
rz(0.57812771) q[1];
sx q[1];
rz(-1.4923948) q[1];
sx q[1];
rz(-2.633458) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.098658853) q[0];
sx q[0];
rz(-1.9693976) q[0];
sx q[0];
rz(-2.0304785) q[0];
rz(-0.63799413) q[2];
sx q[2];
rz(-1.5533414) q[2];
sx q[2];
rz(-0.88776112) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9307738) q[1];
sx q[1];
rz(-2.3465119) q[1];
sx q[1];
rz(0.37097431) q[1];
rz(-pi) q[2];
rz(-2.1448137) q[3];
sx q[3];
rz(-0.32696163) q[3];
sx q[3];
rz(-2.5329778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4529139) q[2];
sx q[2];
rz(-0.33334175) q[2];
sx q[2];
rz(-0.63344669) q[2];
rz(-0.59988919) q[3];
sx q[3];
rz(-1.9918631) q[3];
sx q[3];
rz(-1.5002804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3863581) q[0];
sx q[0];
rz(-0.30325493) q[0];
sx q[0];
rz(-0.19609837) q[0];
rz(-1.261699) q[1];
sx q[1];
rz(-2.321107) q[1];
sx q[1];
rz(2.0702147) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58760658) q[0];
sx q[0];
rz(-1.1452132) q[0];
sx q[0];
rz(-2.2527184) q[0];
rz(-2.0256261) q[2];
sx q[2];
rz(-2.358846) q[2];
sx q[2];
rz(-0.46403971) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2744659) q[1];
sx q[1];
rz(-1.9879838) q[1];
sx q[1];
rz(2.8309114) q[1];
rz(-2.3134872) q[3];
sx q[3];
rz(-0.58907408) q[3];
sx q[3];
rz(0.38618726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1065958) q[2];
sx q[2];
rz(-1.3369766) q[2];
sx q[2];
rz(0.98199797) q[2];
rz(2.9563831) q[3];
sx q[3];
rz(-0.84398142) q[3];
sx q[3];
rz(-1.265032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5380602) q[0];
sx q[0];
rz(-0.91941994) q[0];
sx q[0];
rz(-2.61125) q[0];
rz(-1.2999339) q[1];
sx q[1];
rz(-1.329774) q[1];
sx q[1];
rz(-2.9249654) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7055571) q[0];
sx q[0];
rz(-0.85058054) q[0];
sx q[0];
rz(0.35465045) q[0];
x q[1];
rz(-1.8703307) q[2];
sx q[2];
rz(-0.96200633) q[2];
sx q[2];
rz(1.8744206) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.12946373) q[1];
sx q[1];
rz(-0.33540091) q[1];
sx q[1];
rz(1.5957551) q[1];
rz(-pi) q[2];
rz(0.19200237) q[3];
sx q[3];
rz(-0.90930206) q[3];
sx q[3];
rz(2.0249174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4914322) q[2];
sx q[2];
rz(-1.6314793) q[2];
sx q[2];
rz(2.1177297) q[2];
rz(2.2653545) q[3];
sx q[3];
rz(-0.70228464) q[3];
sx q[3];
rz(1.8700301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.4890471) q[0];
sx q[0];
rz(-1.9788195) q[0];
sx q[0];
rz(0.52893692) q[0];
rz(1.6128929) q[1];
sx q[1];
rz(-1.1922319) q[1];
sx q[1];
rz(2.0524009) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76656094) q[0];
sx q[0];
rz(-1.2732817) q[0];
sx q[0];
rz(-3.1356698) q[0];
x q[1];
rz(-2.5838926) q[2];
sx q[2];
rz(-2.2039362) q[2];
sx q[2];
rz(-1.5932839) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.293922) q[1];
sx q[1];
rz(-1.5772181) q[1];
sx q[1];
rz(-0.033172219) q[1];
rz(3.0508556) q[3];
sx q[3];
rz(-1.4416579) q[3];
sx q[3];
rz(1.7307626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0163394) q[2];
sx q[2];
rz(-2.3374127) q[2];
sx q[2];
rz(-1.9160697) q[2];
rz(1.6493753) q[3];
sx q[3];
rz(-1.7047313) q[3];
sx q[3];
rz(-0.15587458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4847223) q[0];
sx q[0];
rz(-3.0796034) q[0];
sx q[0];
rz(-0.86762506) q[0];
rz(0.067226974) q[1];
sx q[1];
rz(-1.0311238) q[1];
sx q[1];
rz(-2.9464088) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7282043) q[0];
sx q[0];
rz(-1.2767316) q[0];
sx q[0];
rz(-1.3608027) q[0];
rz(0.033069177) q[2];
sx q[2];
rz(-0.40133921) q[2];
sx q[2];
rz(-0.44905845) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6703549) q[1];
sx q[1];
rz(-1.1413304) q[1];
sx q[1];
rz(0.045713748) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.94594749) q[3];
sx q[3];
rz(-0.56560707) q[3];
sx q[3];
rz(-2.2212976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8937257) q[2];
sx q[2];
rz(-0.21394955) q[2];
sx q[2];
rz(-2.2360738) q[2];
rz(1.9838105) q[3];
sx q[3];
rz(-1.1685305) q[3];
sx q[3];
rz(0.15914966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4147707) q[0];
sx q[0];
rz(-1.0960217) q[0];
sx q[0];
rz(-1.0409521) q[0];
rz(-0.078661593) q[1];
sx q[1];
rz(-0.18053308) q[1];
sx q[1];
rz(2.7862766) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8126292) q[0];
sx q[0];
rz(-1.2041429) q[0];
sx q[0];
rz(3.0366412) q[0];
rz(2.9012868) q[2];
sx q[2];
rz(-1.086364) q[2];
sx q[2];
rz(0.74953178) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4081817) q[1];
sx q[1];
rz(-1.2473885) q[1];
sx q[1];
rz(-0.059493382) q[1];
rz(-pi) q[2];
rz(-2.5411685) q[3];
sx q[3];
rz(-1.7250337) q[3];
sx q[3];
rz(2.141181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.67655247) q[2];
sx q[2];
rz(-0.6074473) q[2];
sx q[2];
rz(-1.4036277) q[2];
rz(-0.0062395652) q[3];
sx q[3];
rz(-1.3429567) q[3];
sx q[3];
rz(-1.5374373) q[3];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2733317) q[0];
sx q[0];
rz(-1.9737759) q[0];
sx q[0];
rz(-1.8433174) q[0];
rz(2.6955993) q[1];
sx q[1];
rz(-2.0263717) q[1];
sx q[1];
rz(2.840852) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8200127) q[0];
sx q[0];
rz(-1.0386779) q[0];
sx q[0];
rz(2.8768455) q[0];
x q[1];
rz(-0.54603521) q[2];
sx q[2];
rz(-1.9249501) q[2];
sx q[2];
rz(-1.491577) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3831331) q[1];
sx q[1];
rz(-2.8499095) q[1];
sx q[1];
rz(-0.056706927) q[1];
rz(2.7281076) q[3];
sx q[3];
rz(-0.95260145) q[3];
sx q[3];
rz(-1.9181044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5131502) q[2];
sx q[2];
rz(-2.0437888) q[2];
sx q[2];
rz(-2.002031) q[2];
rz(1.3509753) q[3];
sx q[3];
rz(-0.59777483) q[3];
sx q[3];
rz(-1.1736419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6931077) q[0];
sx q[0];
rz(-1.2379452) q[0];
sx q[0];
rz(-2.2647279) q[0];
rz(-1.4032455) q[1];
sx q[1];
rz(-1.2259903) q[1];
sx q[1];
rz(-1.2798053) q[1];
rz(0.030056603) q[2];
sx q[2];
rz(-0.056014225) q[2];
sx q[2];
rz(0.87607486) q[2];
rz(-0.26294796) q[3];
sx q[3];
rz(-2.1750952) q[3];
sx q[3];
rz(1.7261214) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
