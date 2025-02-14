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
rz(-1.8104115) q[0];
sx q[0];
rz(-1.6348732) q[0];
sx q[0];
rz(2.9659086) q[0];
rz(1.0653347) q[1];
sx q[1];
rz(4.4720286) q[1];
sx q[1];
rz(8.7695697) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86533812) q[0];
sx q[0];
rz(-1.0840345) q[0];
sx q[0];
rz(-0.36051644) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6263481) q[2];
sx q[2];
rz(-1.4719988) q[2];
sx q[2];
rz(-0.21888079) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.53319327) q[1];
sx q[1];
rz(-3.0113868) q[1];
sx q[1];
rz(-0.59334468) q[1];
x q[2];
rz(1.1635861) q[3];
sx q[3];
rz(-2.4904479) q[3];
sx q[3];
rz(0.69738392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4483999) q[2];
sx q[2];
rz(-3.0312263) q[2];
sx q[2];
rz(2.8006862) q[2];
rz(2.36813) q[3];
sx q[3];
rz(-2.1229459) q[3];
sx q[3];
rz(3.036518) q[3];
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
rz(1.676749) q[0];
sx q[0];
rz(-2.8571547) q[0];
sx q[0];
rz(-2.8819717) q[0];
rz(0.78850857) q[1];
sx q[1];
rz(-0.85616833) q[1];
sx q[1];
rz(1.8964918) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7126078) q[0];
sx q[0];
rz(-2.1855506) q[0];
sx q[0];
rz(0.12354596) q[0];
rz(-pi) q[1];
rz(-0.2971895) q[2];
sx q[2];
rz(-2.8126394) q[2];
sx q[2];
rz(-2.7847814) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.81345612) q[1];
sx q[1];
rz(-0.81902796) q[1];
sx q[1];
rz(1.0830031) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8471242) q[3];
sx q[3];
rz(-2.614902) q[3];
sx q[3];
rz(-1.1931488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.239324) q[2];
sx q[2];
rz(-2.6332492) q[2];
sx q[2];
rz(2.2321205) q[2];
rz(-2.4746312) q[3];
sx q[3];
rz(-1.4374377) q[3];
sx q[3];
rz(3.0392569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7143836) q[0];
sx q[0];
rz(-0.41218555) q[0];
sx q[0];
rz(-1.3877731) q[0];
rz(2.1448403) q[1];
sx q[1];
rz(-0.69147384) q[1];
sx q[1];
rz(-0.48666418) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62138689) q[0];
sx q[0];
rz(-0.4356111) q[0];
sx q[0];
rz(1.9389047) q[0];
rz(-pi) q[1];
rz(2.4805807) q[2];
sx q[2];
rz(-2.1901988) q[2];
sx q[2];
rz(-2.7255779) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6909161) q[1];
sx q[1];
rz(-1.449739) q[1];
sx q[1];
rz(-2.6208505) q[1];
rz(-pi) q[2];
x q[2];
rz(0.9439982) q[3];
sx q[3];
rz(-1.4508012) q[3];
sx q[3];
rz(-1.7606408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.77728689) q[2];
sx q[2];
rz(-1.8824258) q[2];
sx q[2];
rz(-0.49059179) q[2];
rz(-0.81654882) q[3];
sx q[3];
rz(-2.3994763) q[3];
sx q[3];
rz(-1.0402927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33646026) q[0];
sx q[0];
rz(-1.8959624) q[0];
sx q[0];
rz(-2.1601987) q[0];
rz(2.7384752) q[1];
sx q[1];
rz(-0.50794452) q[1];
sx q[1];
rz(2.6037762) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6526529) q[0];
sx q[0];
rz(-1.8171628) q[0];
sx q[0];
rz(-0.15423488) q[0];
x q[1];
rz(-0.13366416) q[2];
sx q[2];
rz(-0.85459177) q[2];
sx q[2];
rz(1.4203526) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.06925) q[1];
sx q[1];
rz(-1.9543271) q[1];
sx q[1];
rz(1.0479397) q[1];
rz(-pi) q[2];
rz(-1.9547525) q[3];
sx q[3];
rz(-2.4297485) q[3];
sx q[3];
rz(1.8432338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2526523) q[2];
sx q[2];
rz(-2.6344968) q[2];
sx q[2];
rz(1.0520774) q[2];
rz(2.5893411) q[3];
sx q[3];
rz(-1.735732) q[3];
sx q[3];
rz(1.2205869) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.030516142) q[0];
sx q[0];
rz(-1.3330326) q[0];
sx q[0];
rz(3.1040177) q[0];
rz(-0.75282085) q[1];
sx q[1];
rz(-1.5767187) q[1];
sx q[1];
rz(3.0573696) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60606979) q[0];
sx q[0];
rz(-1.3227934) q[0];
sx q[0];
rz(2.1169784) q[0];
rz(1.2758314) q[2];
sx q[2];
rz(-1.9199445) q[2];
sx q[2];
rz(-0.63284208) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3854691) q[1];
sx q[1];
rz(-1.1471355) q[1];
sx q[1];
rz(-0.045457289) q[1];
rz(-pi) q[2];
rz(3.0952206) q[3];
sx q[3];
rz(-1.4102458) q[3];
sx q[3];
rz(0.35453803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.61421824) q[2];
sx q[2];
rz(-0.97443333) q[2];
sx q[2];
rz(0.43240377) q[2];
rz(1.6480986) q[3];
sx q[3];
rz(-1.7688388) q[3];
sx q[3];
rz(-2.4026332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37345165) q[0];
sx q[0];
rz(-1.5274763) q[0];
sx q[0];
rz(0.76171184) q[0];
rz(1.0234045) q[1];
sx q[1];
rz(-2.6506212) q[1];
sx q[1];
rz(-2.8010119) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0878076) q[0];
sx q[0];
rz(-1.6594961) q[0];
sx q[0];
rz(-0.78194639) q[0];
rz(-pi) q[1];
x q[1];
rz(0.24575524) q[2];
sx q[2];
rz(-0.57502103) q[2];
sx q[2];
rz(-1.0005282) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5300776) q[1];
sx q[1];
rz(-2.7002091) q[1];
sx q[1];
rz(0.89479052) q[1];
rz(-1.7238925) q[3];
sx q[3];
rz(-2.3153437) q[3];
sx q[3];
rz(1.561867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9684101) q[2];
sx q[2];
rz(-2.0081655) q[2];
sx q[2];
rz(2.3859712) q[2];
rz(-2.7544045) q[3];
sx q[3];
rz(-1.7500992) q[3];
sx q[3];
rz(-0.79425991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25283915) q[0];
sx q[0];
rz(-3.0536953) q[0];
sx q[0];
rz(1.9320236) q[0];
rz(-1.6190489) q[1];
sx q[1];
rz(-0.77170283) q[1];
sx q[1];
rz(-1.3873772) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98859859) q[0];
sx q[0];
rz(-1.8805046) q[0];
sx q[0];
rz(-1.214908) q[0];
rz(-pi) q[1];
rz(-2.8262994) q[2];
sx q[2];
rz(-1.0206501) q[2];
sx q[2];
rz(1.2619051) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8232812) q[1];
sx q[1];
rz(-2.2522587) q[1];
sx q[1];
rz(1.0266515) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.58297662) q[3];
sx q[3];
rz(-2.600353) q[3];
sx q[3];
rz(1.4547294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9428955) q[2];
sx q[2];
rz(-0.2936475) q[2];
sx q[2];
rz(-2.7315308) q[2];
rz(0.53747082) q[3];
sx q[3];
rz(-1.042807) q[3];
sx q[3];
rz(-1.0679831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.033567) q[0];
sx q[0];
rz(-1.9060059) q[0];
sx q[0];
rz(2.1349452) q[0];
rz(-2.0058696) q[1];
sx q[1];
rz(-2.6508811) q[1];
sx q[1];
rz(2.8699285) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1170809) q[0];
sx q[0];
rz(-2.4371464) q[0];
sx q[0];
rz(-2.4570877) q[0];
rz(0.92046787) q[2];
sx q[2];
rz(-1.0050541) q[2];
sx q[2];
rz(2.5793902) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2073839) q[1];
sx q[1];
rz(-1.0107688) q[1];
sx q[1];
rz(2.4337342) q[1];
x q[2];
rz(1.3471421) q[3];
sx q[3];
rz(-0.8394548) q[3];
sx q[3];
rz(2.770693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2037619) q[2];
sx q[2];
rz(-2.2088642) q[2];
sx q[2];
rz(-2.0358098) q[2];
rz(1.743099) q[3];
sx q[3];
rz(-2.1574557) q[3];
sx q[3];
rz(0.85247803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49649134) q[0];
sx q[0];
rz(-0.95106769) q[0];
sx q[0];
rz(2.7880461) q[0];
rz(-1.722466) q[1];
sx q[1];
rz(-1.6058233) q[1];
sx q[1];
rz(3.0790192) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1997027) q[0];
sx q[0];
rz(-1.0257922) q[0];
sx q[0];
rz(-0.99129321) q[0];
rz(0.75980564) q[2];
sx q[2];
rz(-1.9927597) q[2];
sx q[2];
rz(-2.9874731) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.313301) q[1];
sx q[1];
rz(-0.57692674) q[1];
sx q[1];
rz(-0.58842701) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.479019) q[3];
sx q[3];
rz(-1.4972613) q[3];
sx q[3];
rz(0.92360332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.320437) q[2];
sx q[2];
rz(-1.7675567) q[2];
sx q[2];
rz(0.42018166) q[2];
rz(0.34475103) q[3];
sx q[3];
rz(-0.77778608) q[3];
sx q[3];
rz(-0.93728089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.47278136) q[0];
sx q[0];
rz(-2.9720699) q[0];
sx q[0];
rz(-1.2813168) q[0];
rz(-1.5538813) q[1];
sx q[1];
rz(-2.3322767) q[1];
sx q[1];
rz(2.4344427) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67884655) q[0];
sx q[0];
rz(-0.94043676) q[0];
sx q[0];
rz(-0.014883777) q[0];
rz(-2.1936839) q[2];
sx q[2];
rz(-0.68205183) q[2];
sx q[2];
rz(2.1932604) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.36701074) q[1];
sx q[1];
rz(-0.38644192) q[1];
sx q[1];
rz(-2.5143753) q[1];
x q[2];
rz(-2.3570746) q[3];
sx q[3];
rz(-2.3459179) q[3];
sx q[3];
rz(-2.8592542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3215434) q[2];
sx q[2];
rz(-2.0115435) q[2];
sx q[2];
rz(-0.496544) q[2];
rz(-1.2921565) q[3];
sx q[3];
rz(-2.8437331) q[3];
sx q[3];
rz(2.7636102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15688607) q[0];
sx q[0];
rz(-2.8732185) q[0];
sx q[0];
rz(2.1068841) q[0];
rz(-2.5490419) q[1];
sx q[1];
rz(-2.4609346) q[1];
sx q[1];
rz(1.5998283) q[1];
rz(1.7405657) q[2];
sx q[2];
rz(-0.23525722) q[2];
sx q[2];
rz(2.7973882) q[2];
rz(2.5377688) q[3];
sx q[3];
rz(-1.3841938) q[3];
sx q[3];
rz(-2.9740372) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
