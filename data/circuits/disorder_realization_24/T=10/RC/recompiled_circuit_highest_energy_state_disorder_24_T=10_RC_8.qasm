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
rz(-0.10676323) q[0];
sx q[0];
rz(-1.8101298) q[0];
sx q[0];
rz(-1.3504299) q[0];
rz(-2.8332233) q[1];
sx q[1];
rz(-2.8928533) q[1];
sx q[1];
rz(-2.1623936) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.072974394) q[0];
sx q[0];
rz(-2.3280297) q[0];
sx q[0];
rz(-1.4692151) q[0];
rz(-pi) q[1];
rz(-2.0151695) q[2];
sx q[2];
rz(-0.70894402) q[2];
sx q[2];
rz(-2.5042748) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0268271) q[1];
sx q[1];
rz(-1.9239167) q[1];
sx q[1];
rz(1.4985282) q[1];
x q[2];
rz(-0.37110801) q[3];
sx q[3];
rz(-1.4083997) q[3];
sx q[3];
rz(1.1134256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8637434) q[2];
sx q[2];
rz(-1.9305482) q[2];
sx q[2];
rz(-2.5085874) q[2];
rz(-0.64166075) q[3];
sx q[3];
rz(-2.6810665) q[3];
sx q[3];
rz(-2.3708926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4577456) q[0];
sx q[0];
rz(-2.1774543) q[0];
sx q[0];
rz(-0.82474166) q[0];
rz(1.9174346) q[1];
sx q[1];
rz(-0.85547248) q[1];
sx q[1];
rz(2.8214084) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.510957) q[0];
sx q[0];
rz(-0.63586006) q[0];
sx q[0];
rz(0.73131928) q[0];
rz(-1.0616706) q[2];
sx q[2];
rz(-1.7560282) q[2];
sx q[2];
rz(-0.19134609) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0517438) q[1];
sx q[1];
rz(-1.4898572) q[1];
sx q[1];
rz(-1.5481987) q[1];
rz(1.3489487) q[3];
sx q[3];
rz(-1.6288405) q[3];
sx q[3];
rz(2.0741418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0188633) q[2];
sx q[2];
rz(-1.2168695) q[2];
sx q[2];
rz(1.8625721) q[2];
rz(-3.0458798) q[3];
sx q[3];
rz(-1.0749823) q[3];
sx q[3];
rz(1.8224645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1529481) q[0];
sx q[0];
rz(-2.4929292) q[0];
sx q[0];
rz(1.0308107) q[0];
rz(2.5910494) q[1];
sx q[1];
rz(-2.2830453) q[1];
sx q[1];
rz(1.2446838) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6129788) q[0];
sx q[0];
rz(-1.6242149) q[0];
sx q[0];
rz(2.9821755) q[0];
rz(2.1878476) q[2];
sx q[2];
rz(-1.707952) q[2];
sx q[2];
rz(-2.0037665) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.66164368) q[1];
sx q[1];
rz(-0.91001653) q[1];
sx q[1];
rz(-2.8815072) q[1];
rz(-pi) q[2];
rz(-0.81826415) q[3];
sx q[3];
rz(-1.2845367) q[3];
sx q[3];
rz(0.32901007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2391669) q[2];
sx q[2];
rz(-2.5776165) q[2];
sx q[2];
rz(-1.172056) q[2];
rz(0.82550448) q[3];
sx q[3];
rz(-1.2647311) q[3];
sx q[3];
rz(-2.4765305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.091752) q[0];
sx q[0];
rz(-1.6850152) q[0];
sx q[0];
rz(-2.7276584) q[0];
rz(2.6992758) q[1];
sx q[1];
rz(-1.0041753) q[1];
sx q[1];
rz(0.54534674) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1236432) q[0];
sx q[0];
rz(-1.845903) q[0];
sx q[0];
rz(-1.6226044) q[0];
x q[1];
rz(-0.92939922) q[2];
sx q[2];
rz(-2.2365733) q[2];
sx q[2];
rz(-2.5432472) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1283577) q[1];
sx q[1];
rz(-1.2942794) q[1];
sx q[1];
rz(1.0324296) q[1];
x q[2];
rz(1.3096626) q[3];
sx q[3];
rz(-2.6178837) q[3];
sx q[3];
rz(-1.5098438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.44117323) q[2];
sx q[2];
rz(-1.2607231) q[2];
sx q[2];
rz(-2.9803661) q[2];
rz(1.553933) q[3];
sx q[3];
rz(-0.59993887) q[3];
sx q[3];
rz(-2.5509295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.081838354) q[0];
sx q[0];
rz(-1.2677544) q[0];
sx q[0];
rz(2.4181714) q[0];
rz(1.3149892) q[1];
sx q[1];
rz(-1.864121) q[1];
sx q[1];
rz(2.1037197) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68907952) q[0];
sx q[0];
rz(-0.85001365) q[0];
sx q[0];
rz(1.6759273) q[0];
x q[1];
rz(-2.5332301) q[2];
sx q[2];
rz(-2.3181097) q[2];
sx q[2];
rz(1.057511) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2813869) q[1];
sx q[1];
rz(-1.1347927) q[1];
sx q[1];
rz(1.866775) q[1];
rz(-pi) q[2];
x q[2];
rz(0.54609046) q[3];
sx q[3];
rz(-1.6973064) q[3];
sx q[3];
rz(1.5777335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3967241) q[2];
sx q[2];
rz(-2.0266271) q[2];
sx q[2];
rz(0.57207668) q[2];
rz(2.5165596) q[3];
sx q[3];
rz(-1.0038989) q[3];
sx q[3];
rz(-1.2264138) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52142757) q[0];
sx q[0];
rz(-3.1075952) q[0];
sx q[0];
rz(-0.52325621) q[0];
rz(-1.8190067) q[1];
sx q[1];
rz(-1.6736284) q[1];
sx q[1];
rz(-2.5214213) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5145743) q[0];
sx q[0];
rz(-0.73643273) q[0];
sx q[0];
rz(0.66177841) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.46862009) q[2];
sx q[2];
rz(-1.5053362) q[2];
sx q[2];
rz(0.58131389) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5413618) q[1];
sx q[1];
rz(-0.63284749) q[1];
sx q[1];
rz(-1.4221627) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.2175339) q[3];
sx q[3];
rz(-1.399588) q[3];
sx q[3];
rz(1.8678766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.72914499) q[2];
sx q[2];
rz(-0.77249384) q[2];
sx q[2];
rz(1.283851) q[2];
rz(-0.9681975) q[3];
sx q[3];
rz(-1.7627534) q[3];
sx q[3];
rz(1.5865954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-0.17484434) q[0];
sx q[0];
rz(-2.0059678) q[0];
sx q[0];
rz(1.9298166) q[0];
rz(0.9696331) q[1];
sx q[1];
rz(-1.8200834) q[1];
sx q[1];
rz(1.2744354) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2454555) q[0];
sx q[0];
rz(-1.7217759) q[0];
sx q[0];
rz(-2.2666988) q[0];
x q[1];
rz(1.4757122) q[2];
sx q[2];
rz(-0.78721744) q[2];
sx q[2];
rz(-1.1457535) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3434747) q[1];
sx q[1];
rz(-2.1785582) q[1];
sx q[1];
rz(2.7704984) q[1];
x q[2];
rz(-2.681143) q[3];
sx q[3];
rz(-1.2795953) q[3];
sx q[3];
rz(1.2429383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4659287) q[2];
sx q[2];
rz(-1.6589386) q[2];
sx q[2];
rz(-0.3234123) q[2];
rz(0.0023500738) q[3];
sx q[3];
rz(-1.6833143) q[3];
sx q[3];
rz(2.5200444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6482786) q[0];
sx q[0];
rz(-1.4819773) q[0];
sx q[0];
rz(-1.7907273) q[0];
rz(1.3510652) q[1];
sx q[1];
rz(-1.8565145) q[1];
sx q[1];
rz(1.2877119) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5918811) q[0];
sx q[0];
rz(-1.2882107) q[0];
sx q[0];
rz(3.0192399) q[0];
rz(1.9016198) q[2];
sx q[2];
rz(-0.74218732) q[2];
sx q[2];
rz(1.0617219) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8069597) q[1];
sx q[1];
rz(-1.8939085) q[1];
sx q[1];
rz(1.7503225) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.34890449) q[3];
sx q[3];
rz(-1.5530905) q[3];
sx q[3];
rz(2.4573093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.52435851) q[2];
sx q[2];
rz(-1.1523767) q[2];
sx q[2];
rz(2.1481245) q[2];
rz(1.0348381) q[3];
sx q[3];
rz(-0.60641685) q[3];
sx q[3];
rz(-2.9724227) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33191037) q[0];
sx q[0];
rz(-1.9669635) q[0];
sx q[0];
rz(-1.6140953) q[0];
rz(-0.88036674) q[1];
sx q[1];
rz(-2.7208734) q[1];
sx q[1];
rz(2.8299433) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0904734) q[0];
sx q[0];
rz(-1.7684121) q[0];
sx q[0];
rz(-1.2725426) q[0];
rz(-pi) q[1];
rz(-1.8021605) q[2];
sx q[2];
rz(-1.5867434) q[2];
sx q[2];
rz(-2.0248264) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5376725) q[1];
sx q[1];
rz(-0.89960704) q[1];
sx q[1];
rz(2.0449941) q[1];
rz(-pi) q[2];
rz(1.1470471) q[3];
sx q[3];
rz(-2.4422788) q[3];
sx q[3];
rz(2.5831403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6722022) q[2];
sx q[2];
rz(-2.2180874) q[2];
sx q[2];
rz(1.7507318) q[2];
rz(1.6145128) q[3];
sx q[3];
rz(-1.047784) q[3];
sx q[3];
rz(2.0670149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5020318) q[0];
sx q[0];
rz(-1.9968411) q[0];
sx q[0];
rz(-0.61258739) q[0];
rz(2.1514429) q[1];
sx q[1];
rz(-1.9691111) q[1];
sx q[1];
rz(-1.490907) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51524893) q[0];
sx q[0];
rz(-0.83413163) q[0];
sx q[0];
rz(0.5818999) q[0];
rz(-pi) q[1];
rz(0.25509768) q[2];
sx q[2];
rz(-2.2146985) q[2];
sx q[2];
rz(-1.7718499) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9673493) q[1];
sx q[1];
rz(-0.16160204) q[1];
sx q[1];
rz(-1.6095649) q[1];
rz(-pi) q[2];
rz(-2.4435639) q[3];
sx q[3];
rz(-1.7637611) q[3];
sx q[3];
rz(-2.2730227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7154197) q[2];
sx q[2];
rz(-2.3121068) q[2];
sx q[2];
rz(-3.1066011) q[2];
rz(1.70111) q[3];
sx q[3];
rz(-0.8927497) q[3];
sx q[3];
rz(-0.54097241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4424425) q[0];
sx q[0];
rz(-2.2461666) q[0];
sx q[0];
rz(2.9641892) q[0];
rz(-1.9433446) q[1];
sx q[1];
rz(-1.5180963) q[1];
sx q[1];
rz(0.95388283) q[1];
rz(-2.162355) q[2];
sx q[2];
rz(-2.1605347) q[2];
sx q[2];
rz(0.84183358) q[2];
rz(0.51552897) q[3];
sx q[3];
rz(-2.3524108) q[3];
sx q[3];
rz(-2.4084211) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
