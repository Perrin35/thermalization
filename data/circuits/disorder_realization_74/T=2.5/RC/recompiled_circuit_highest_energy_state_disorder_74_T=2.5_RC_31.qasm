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
rz(-2.3251301) q[0];
sx q[0];
rz(-0.10179585) q[0];
sx q[0];
rz(2.6020004) q[0];
rz(3.6379023) q[1];
sx q[1];
rz(3.451347) q[1];
sx q[1];
rz(5.7440905) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4350454) q[0];
sx q[0];
rz(-0.4420949) q[0];
sx q[0];
rz(-3.0821783) q[0];
rz(-pi) q[1];
rz(-2.7560948) q[2];
sx q[2];
rz(-0.38542569) q[2];
sx q[2];
rz(0.62776923) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4161168) q[1];
sx q[1];
rz(-1.892964) q[1];
sx q[1];
rz(-0.14741082) q[1];
rz(-pi) q[2];
x q[2];
rz(0.92181262) q[3];
sx q[3];
rz(-1.6361437) q[3];
sx q[3];
rz(-1.1999038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.76217905) q[2];
sx q[2];
rz(-0.87612408) q[2];
sx q[2];
rz(-1.9525105) q[2];
rz(-1.1842229) q[3];
sx q[3];
rz(-2.2370179) q[3];
sx q[3];
rz(-1.5466461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9963843) q[0];
sx q[0];
rz(-1.5897911) q[0];
sx q[0];
rz(1.8183964) q[0];
rz(2.6595751) q[1];
sx q[1];
rz(-2.2299485) q[1];
sx q[1];
rz(2.1673896) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0486748) q[0];
sx q[0];
rz(-0.67103025) q[0];
sx q[0];
rz(1.8604408) q[0];
x q[1];
rz(-0.85136885) q[2];
sx q[2];
rz(-0.15002827) q[2];
sx q[2];
rz(-2.0025557) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1434553) q[1];
sx q[1];
rz(-2.1051198) q[1];
sx q[1];
rz(-1.8443395) q[1];
rz(-pi) q[2];
rz(-2.7747068) q[3];
sx q[3];
rz(-2.5436333) q[3];
sx q[3];
rz(2.6794527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1366068) q[2];
sx q[2];
rz(-2.5544781) q[2];
sx q[2];
rz(-0.74964398) q[2];
rz(-2.7096115) q[3];
sx q[3];
rz(-1.1155198) q[3];
sx q[3];
rz(1.8384793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5169446) q[0];
sx q[0];
rz(-2.8693146) q[0];
sx q[0];
rz(-2.5352449) q[0];
rz(-1.3336522) q[1];
sx q[1];
rz(-1.2140112) q[1];
sx q[1];
rz(-2.8820754) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0344275) q[0];
sx q[0];
rz(-1.3239064) q[0];
sx q[0];
rz(-2.9811923) q[0];
rz(2.9664842) q[2];
sx q[2];
rz(-1.7705743) q[2];
sx q[2];
rz(1.4005043) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1705679) q[1];
sx q[1];
rz(-1.9924506) q[1];
sx q[1];
rz(0.41280156) q[1];
rz(0.43478888) q[3];
sx q[3];
rz(-2.090095) q[3];
sx q[3];
rz(-3.0360707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8848662) q[2];
sx q[2];
rz(-1.3639516) q[2];
sx q[2];
rz(1.5941031) q[2];
rz(-2.3571842) q[3];
sx q[3];
rz(-1.514785) q[3];
sx q[3];
rz(-1.5659531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6467658) q[0];
sx q[0];
rz(-1.0972728) q[0];
sx q[0];
rz(0.25949091) q[0];
rz(1.9208113) q[1];
sx q[1];
rz(-0.44993284) q[1];
sx q[1];
rz(1.4422013) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9913015) q[0];
sx q[0];
rz(-1.1027005) q[0];
sx q[0];
rz(-2.1670114) q[0];
x q[1];
rz(-0.97564189) q[2];
sx q[2];
rz(-1.4624498) q[2];
sx q[2];
rz(-2.8841022) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4220548) q[1];
sx q[1];
rz(-1.9214848) q[1];
sx q[1];
rz(1.3077875) q[1];
rz(-pi) q[2];
rz(0.070017858) q[3];
sx q[3];
rz(-0.88980955) q[3];
sx q[3];
rz(0.79047608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0356902) q[2];
sx q[2];
rz(-1.3373969) q[2];
sx q[2];
rz(2.7316459) q[2];
rz(1.2116872) q[3];
sx q[3];
rz(-0.68710059) q[3];
sx q[3];
rz(1.3338026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7164417) q[0];
sx q[0];
rz(-2.4254159) q[0];
sx q[0];
rz(2.8675365) q[0];
rz(2.7033499) q[1];
sx q[1];
rz(-1.2934877) q[1];
sx q[1];
rz(-2.4058707) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0276707) q[0];
sx q[0];
rz(-0.7589853) q[0];
sx q[0];
rz(3.0149197) q[0];
rz(-pi) q[1];
rz(-1.7154791) q[2];
sx q[2];
rz(-1.4229454) q[2];
sx q[2];
rz(1.3593909) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2949038) q[1];
sx q[1];
rz(-1.3219993) q[1];
sx q[1];
rz(-1.873068) q[1];
x q[2];
rz(-2.6295119) q[3];
sx q[3];
rz(-2.1554073) q[3];
sx q[3];
rz(-1.1325815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.938544) q[2];
sx q[2];
rz(-1.4578578) q[2];
sx q[2];
rz(-0.61384821) q[2];
rz(-0.40890536) q[3];
sx q[3];
rz(-2.1730065) q[3];
sx q[3];
rz(-2.3992505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78855377) q[0];
sx q[0];
rz(-2.5034294) q[0];
sx q[0];
rz(1.9045389) q[0];
rz(-1.6963814) q[1];
sx q[1];
rz(-2.684869) q[1];
sx q[1];
rz(2.9663185) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0171368) q[0];
sx q[0];
rz(-1.5457898) q[0];
sx q[0];
rz(-1.3978005) q[0];
rz(-pi) q[1];
rz(0.66368033) q[2];
sx q[2];
rz(-0.77518565) q[2];
sx q[2];
rz(-0.41942393) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3289017) q[1];
sx q[1];
rz(-1.1393424) q[1];
sx q[1];
rz(-0.56568362) q[1];
rz(-pi) q[2];
rz(2.184559) q[3];
sx q[3];
rz(-2.0599457) q[3];
sx q[3];
rz(3.1135984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1767629) q[2];
sx q[2];
rz(-1.8076597) q[2];
sx q[2];
rz(-2.1691587) q[2];
rz(-1.4771627) q[3];
sx q[3];
rz(-1.9921314) q[3];
sx q[3];
rz(2.3798063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85457388) q[0];
sx q[0];
rz(-3.119097) q[0];
sx q[0];
rz(2.7884685) q[0];
rz(-1.0278206) q[1];
sx q[1];
rz(-1.2007583) q[1];
sx q[1];
rz(2.1305398) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42732692) q[0];
sx q[0];
rz(-1.9280701) q[0];
sx q[0];
rz(-1.5861804) q[0];
rz(-pi) q[1];
rz(-2.2765144) q[2];
sx q[2];
rz(-2.0524244) q[2];
sx q[2];
rz(0.31246802) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8098896) q[1];
sx q[1];
rz(-1.9192438) q[1];
sx q[1];
rz(0.16683677) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8500438) q[3];
sx q[3];
rz(-1.8729775) q[3];
sx q[3];
rz(-0.19552375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.31356835) q[2];
sx q[2];
rz(-2.821533) q[2];
sx q[2];
rz(1.3671406) q[2];
rz(-1.2683055) q[3];
sx q[3];
rz(-1.6627848) q[3];
sx q[3];
rz(0.84793276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1139514) q[0];
sx q[0];
rz(-0.60751644) q[0];
sx q[0];
rz(-2.0116346) q[0];
rz(-2.9226411) q[1];
sx q[1];
rz(-1.7349225) q[1];
sx q[1];
rz(0.91046441) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.509678) q[0];
sx q[0];
rz(-0.70075894) q[0];
sx q[0];
rz(-0.52978446) q[0];
rz(2.6633419) q[2];
sx q[2];
rz(-1.646981) q[2];
sx q[2];
rz(-2.8752799) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3748951) q[1];
sx q[1];
rz(-1.3037221) q[1];
sx q[1];
rz(-2.1007936) q[1];
rz(-0.72515709) q[3];
sx q[3];
rz(-1.8817543) q[3];
sx q[3];
rz(2.6127315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.8255446) q[2];
sx q[2];
rz(-0.80222183) q[2];
sx q[2];
rz(0.82553378) q[2];
rz(-0.46142203) q[3];
sx q[3];
rz(-1.8749219) q[3];
sx q[3];
rz(0.44574827) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6398741) q[0];
sx q[0];
rz(-1.3383144) q[0];
sx q[0];
rz(-2.3097532) q[0];
rz(0.72744751) q[1];
sx q[1];
rz(-1.4201545) q[1];
sx q[1];
rz(0.096253455) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1883586) q[0];
sx q[0];
rz(-1.112845) q[0];
sx q[0];
rz(-1.0012549) q[0];
rz(-1.6321858) q[2];
sx q[2];
rz(-1.824531) q[2];
sx q[2];
rz(2.369427) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.61273328) q[1];
sx q[1];
rz(-1.7816356) q[1];
sx q[1];
rz(3.0490521) q[1];
x q[2];
rz(-1.3745802) q[3];
sx q[3];
rz(-1.726578) q[3];
sx q[3];
rz(0.84583144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40621296) q[2];
sx q[2];
rz(-1.7071743) q[2];
sx q[2];
rz(-1.5489138) q[2];
rz(-2.5088572) q[3];
sx q[3];
rz(-1.9097208) q[3];
sx q[3];
rz(-2.8922141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7525472) q[0];
sx q[0];
rz(-2.1010375) q[0];
sx q[0];
rz(2.7885875) q[0];
rz(2.8189335) q[1];
sx q[1];
rz(-2.8162075) q[1];
sx q[1];
rz(2.7696612) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7854561) q[0];
sx q[0];
rz(-1.8710941) q[0];
sx q[0];
rz(-0.2072643) q[0];
rz(0.32970365) q[2];
sx q[2];
rz(-1.7404544) q[2];
sx q[2];
rz(-0.28126954) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0363716) q[1];
sx q[1];
rz(-0.37215713) q[1];
sx q[1];
rz(-2.6879265) q[1];
rz(-pi) q[2];
x q[2];
rz(0.84623611) q[3];
sx q[3];
rz(-1.5414943) q[3];
sx q[3];
rz(1.3361564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.17452621) q[2];
sx q[2];
rz(-1.1897949) q[2];
sx q[2];
rz(1.8444427) q[2];
rz(-0.8477115) q[3];
sx q[3];
rz(-1.1573557) q[3];
sx q[3];
rz(-2.1367836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90947718) q[0];
sx q[0];
rz(-1.3790601) q[0];
sx q[0];
rz(-2.7098304) q[0];
rz(1.3879981) q[1];
sx q[1];
rz(-1.9276062) q[1];
sx q[1];
rz(3.066317) q[1];
rz(-2.3222011) q[2];
sx q[2];
rz(-2.1234305) q[2];
sx q[2];
rz(2.7101868) q[2];
rz(-2.344178) q[3];
sx q[3];
rz(-1.8731464) q[3];
sx q[3];
rz(-2.1220589) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
