OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.55968094) q[0];
sx q[0];
rz(2.0547325) q[0];
sx q[0];
rz(7.6261043) q[0];
rz(7.8328447) q[1];
sx q[1];
rz(3.0631493) q[1];
sx q[1];
rz(6.9258239) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6519878) q[0];
sx q[0];
rz(-1.542924) q[0];
sx q[0];
rz(-2.6023988) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6334555) q[2];
sx q[2];
rz(-1.7388441) q[2];
sx q[2];
rz(-0.44067581) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5832311) q[1];
sx q[1];
rz(-1.6915295) q[1];
sx q[1];
rz(-1.6912778) q[1];
rz(-0.1518941) q[3];
sx q[3];
rz(-0.23670247) q[3];
sx q[3];
rz(-3.0083002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6671483) q[2];
sx q[2];
rz(-2.4953304) q[2];
sx q[2];
rz(1.2791963) q[2];
rz(-2.4228418) q[3];
sx q[3];
rz(-1.5712534) q[3];
sx q[3];
rz(0.32354245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6779125) q[0];
sx q[0];
rz(-1.3543411) q[0];
sx q[0];
rz(2.1610778) q[0];
rz(-2.9837043) q[1];
sx q[1];
rz(-1.4973367) q[1];
sx q[1];
rz(-0.78871361) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6533587) q[0];
sx q[0];
rz(-1.6846859) q[0];
sx q[0];
rz(1.4466404) q[0];
rz(-pi) q[1];
rz(-2.8854495) q[2];
sx q[2];
rz(-1.6015341) q[2];
sx q[2];
rz(-2.5275633) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2930254) q[1];
sx q[1];
rz(-1.8752521) q[1];
sx q[1];
rz(0.57363631) q[1];
x q[2];
rz(0.73123587) q[3];
sx q[3];
rz(-0.9160708) q[3];
sx q[3];
rz(0.54192858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.366189) q[2];
sx q[2];
rz(-2.0663694) q[2];
sx q[2];
rz(2.9555087) q[2];
rz(-2.4880593) q[3];
sx q[3];
rz(-1.7553522) q[3];
sx q[3];
rz(-3.0236566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2114975) q[0];
sx q[0];
rz(-0.94267541) q[0];
sx q[0];
rz(-2.511456) q[0];
rz(3.0139626) q[1];
sx q[1];
rz(-2.4928513) q[1];
sx q[1];
rz(-0.72174597) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5436514) q[0];
sx q[0];
rz(-1.8231892) q[0];
sx q[0];
rz(0.6681722) q[0];
rz(-2.4475054) q[2];
sx q[2];
rz(-1.9352479) q[2];
sx q[2];
rz(1.2823766) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7194781) q[1];
sx q[1];
rz(-1.1269224) q[1];
sx q[1];
rz(-1.6133973) q[1];
rz(-pi) q[2];
rz(-0.33629041) q[3];
sx q[3];
rz(-0.97067562) q[3];
sx q[3];
rz(2.1093413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3815986) q[2];
sx q[2];
rz(-0.95278946) q[2];
sx q[2];
rz(2.2325113) q[2];
rz(-2.6233853) q[3];
sx q[3];
rz(-2.0776694) q[3];
sx q[3];
rz(2.853493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58406126) q[0];
sx q[0];
rz(-2.4052305) q[0];
sx q[0];
rz(0.042536143) q[0];
rz(0.78035367) q[1];
sx q[1];
rz(-2.6413554) q[1];
sx q[1];
rz(2.3775878) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7351748) q[0];
sx q[0];
rz(-1.58332) q[0];
sx q[0];
rz(-0.89212117) q[0];
x q[1];
rz(-0.1673844) q[2];
sx q[2];
rz(-2.4902654) q[2];
sx q[2];
rz(1.960388) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2149787) q[1];
sx q[1];
rz(-1.1264631) q[1];
sx q[1];
rz(0.39919969) q[1];
x q[2];
rz(0.56960168) q[3];
sx q[3];
rz(-2.3929425) q[3];
sx q[3];
rz(2.1514055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2970695) q[2];
sx q[2];
rz(-1.2244747) q[2];
sx q[2];
rz(-0.53304535) q[2];
rz(2.879203) q[3];
sx q[3];
rz(-2.5049987) q[3];
sx q[3];
rz(-0.46245241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2247291) q[0];
sx q[0];
rz(-0.30039766) q[0];
sx q[0];
rz(-1.2438783) q[0];
rz(-2.9149756) q[1];
sx q[1];
rz(-0.77426568) q[1];
sx q[1];
rz(2.7382543) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6809083) q[0];
sx q[0];
rz(-1.7366689) q[0];
sx q[0];
rz(-2.9537863) q[0];
rz(0.68571217) q[2];
sx q[2];
rz(-1.2843686) q[2];
sx q[2];
rz(3.0428257) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.22002815) q[1];
sx q[1];
rz(-2.0953396) q[1];
sx q[1];
rz(1.3925874) q[1];
rz(-pi) q[2];
rz(-1.0395398) q[3];
sx q[3];
rz(-0.96635339) q[3];
sx q[3];
rz(-1.8464551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8205745) q[2];
sx q[2];
rz(-1.0514739) q[2];
sx q[2];
rz(0.78249758) q[2];
rz(1.1123505) q[3];
sx q[3];
rz(-2.7045043) q[3];
sx q[3];
rz(-2.1330244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.430442) q[0];
sx q[0];
rz(-0.40200457) q[0];
sx q[0];
rz(0.45561403) q[0];
rz(-3.0420711) q[1];
sx q[1];
rz(-1.1385304) q[1];
sx q[1];
rz(-3.1351556) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59505263) q[0];
sx q[0];
rz(-1.8227302) q[0];
sx q[0];
rz(-2.3810054) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1941031) q[2];
sx q[2];
rz(-1.6695392) q[2];
sx q[2];
rz(-1.8284947) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2575063) q[1];
sx q[1];
rz(-1.3370561) q[1];
sx q[1];
rz(1.1463548) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3901859) q[3];
sx q[3];
rz(-0.91114984) q[3];
sx q[3];
rz(-0.71099647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7727938) q[2];
sx q[2];
rz(-1.5051179) q[2];
sx q[2];
rz(2.2951365) q[2];
rz(-0.99772292) q[3];
sx q[3];
rz(-2.3322451) q[3];
sx q[3];
rz(1.6830106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0434175) q[0];
sx q[0];
rz(-0.88413969) q[0];
sx q[0];
rz(0.055710677) q[0];
rz(-2.3588691) q[1];
sx q[1];
rz(-1.1306154) q[1];
sx q[1];
rz(-1.9810716) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6768778) q[0];
sx q[0];
rz(-0.58266312) q[0];
sx q[0];
rz(0.90375264) q[0];
rz(-pi) q[1];
rz(-2.2968282) q[2];
sx q[2];
rz(-0.65462199) q[2];
sx q[2];
rz(0.95915937) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6351663) q[1];
sx q[1];
rz(-2.4392358) q[1];
sx q[1];
rz(-2.173645) q[1];
rz(-pi) q[2];
rz(2.7619744) q[3];
sx q[3];
rz(-1.2109204) q[3];
sx q[3];
rz(-2.4740263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.065585) q[2];
sx q[2];
rz(-2.2202754) q[2];
sx q[2];
rz(-0.78545061) q[2];
rz(0.75585946) q[3];
sx q[3];
rz(-1.8463585) q[3];
sx q[3];
rz(0.41539645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8287559) q[0];
sx q[0];
rz(-0.090933196) q[0];
sx q[0];
rz(1.1428517) q[0];
rz(1.3061334) q[1];
sx q[1];
rz(-1.1885234) q[1];
sx q[1];
rz(0.41608861) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.04836719) q[0];
sx q[0];
rz(-2.0443633) q[0];
sx q[0];
rz(-2.918539) q[0];
rz(-2.9990254) q[2];
sx q[2];
rz(-2.0349742) q[2];
sx q[2];
rz(-0.67205059) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5059698) q[1];
sx q[1];
rz(-1.1715679) q[1];
sx q[1];
rz(-1.3969621) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11085005) q[3];
sx q[3];
rz(-2.0409611) q[3];
sx q[3];
rz(-1.6925616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0351506) q[2];
sx q[2];
rz(-1.4539377) q[2];
sx q[2];
rz(-0.32315928) q[2];
rz(-2.9390826) q[3];
sx q[3];
rz(-1.2812307) q[3];
sx q[3];
rz(2.5642853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37339661) q[0];
sx q[0];
rz(-2.51238) q[0];
sx q[0];
rz(1.1449822) q[0];
rz(1.1960944) q[1];
sx q[1];
rz(-2.9856666) q[1];
sx q[1];
rz(2.6224565) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7767169) q[0];
sx q[0];
rz(-1.4380699) q[0];
sx q[0];
rz(1.9507292) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9197308) q[2];
sx q[2];
rz(-1.8010745) q[2];
sx q[2];
rz(0.78267539) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5362894) q[1];
sx q[1];
rz(-1.481206) q[1];
sx q[1];
rz(0.78219608) q[1];
x q[2];
rz(2.9037644) q[3];
sx q[3];
rz(-1.3999709) q[3];
sx q[3];
rz(1.886614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8273932) q[2];
sx q[2];
rz(-1.9284356) q[2];
sx q[2];
rz(-3.1006151) q[2];
rz(-0.86769062) q[3];
sx q[3];
rz(-2.4882081) q[3];
sx q[3];
rz(-2.6303671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7460019) q[0];
sx q[0];
rz(-0.87930167) q[0];
sx q[0];
rz(1.6145153) q[0];
rz(-1.7136259) q[1];
sx q[1];
rz(-2.7797996) q[1];
sx q[1];
rz(1.6428927) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7259827) q[0];
sx q[0];
rz(-1.590953) q[0];
sx q[0];
rz(-1.6462506) q[0];
x q[1];
rz(-3.0936436) q[2];
sx q[2];
rz(-1.7257294) q[2];
sx q[2];
rz(-0.63873728) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8720819) q[1];
sx q[1];
rz(-1.2334358) q[1];
sx q[1];
rz(-2.3245287) q[1];
rz(-1.4478217) q[3];
sx q[3];
rz(-1.7719367) q[3];
sx q[3];
rz(-0.5826544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3830118) q[2];
sx q[2];
rz(-1.0772394) q[2];
sx q[2];
rz(-0.14653462) q[2];
rz(-0.81418973) q[3];
sx q[3];
rz(-2.345293) q[3];
sx q[3];
rz(-2.7248785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9183337) q[0];
sx q[0];
rz(-1.2467361) q[0];
sx q[0];
rz(0.99714739) q[0];
rz(-1.2659484) q[1];
sx q[1];
rz(-1.0026149) q[1];
sx q[1];
rz(1.2276585) q[1];
rz(-2.0064034) q[2];
sx q[2];
rz(-0.26293593) q[2];
sx q[2];
rz(1.5756366) q[2];
rz(-0.72600611) q[3];
sx q[3];
rz(-0.86961679) q[3];
sx q[3];
rz(1.1179954) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];