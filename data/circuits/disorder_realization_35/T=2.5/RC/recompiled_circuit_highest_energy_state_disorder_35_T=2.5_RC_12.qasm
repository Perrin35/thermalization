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
rz(2.2284722) q[0];
sx q[0];
rz(4.1780317) q[0];
sx q[0];
rz(10.960624) q[0];
rz(0.66134557) q[1];
sx q[1];
rz(-1.1868492) q[1];
sx q[1];
rz(-2.7050736) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3876037) q[0];
sx q[0];
rz(-1.4719098) q[0];
sx q[0];
rz(1.8900699) q[0];
rz(-pi) q[1];
rz(-0.65459697) q[2];
sx q[2];
rz(-0.67383781) q[2];
sx q[2];
rz(-1.0985653) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4492717) q[1];
sx q[1];
rz(-0.8682478) q[1];
sx q[1];
rz(0.0535123) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9645676) q[3];
sx q[3];
rz(-2.2469871) q[3];
sx q[3];
rz(2.2023622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.46301699) q[2];
sx q[2];
rz(-2.7028694) q[2];
sx q[2];
rz(-2.2185745) q[2];
rz(0.065741278) q[3];
sx q[3];
rz(-1.3275423) q[3];
sx q[3];
rz(-1.340723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1317257) q[0];
sx q[0];
rz(-2.0781524) q[0];
sx q[0];
rz(3.104082) q[0];
rz(-2.5610979) q[1];
sx q[1];
rz(-2.4721804) q[1];
sx q[1];
rz(0.69211778) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0789675) q[0];
sx q[0];
rz(-0.37165549) q[0];
sx q[0];
rz(1.070648) q[0];
x q[1];
rz(-0.28022556) q[2];
sx q[2];
rz(-1.0839274) q[2];
sx q[2];
rz(-3.1161199) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.925732) q[1];
sx q[1];
rz(-2.070148) q[1];
sx q[1];
rz(0.66470048) q[1];
rz(-0.18499891) q[3];
sx q[3];
rz(-0.76055148) q[3];
sx q[3];
rz(-2.4724677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6633501) q[2];
sx q[2];
rz(-1.9448091) q[2];
sx q[2];
rz(-1.3956068) q[2];
rz(2.9546402) q[3];
sx q[3];
rz(-1.170265) q[3];
sx q[3];
rz(2.601534) q[3];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4726987) q[0];
sx q[0];
rz(-2.5085594) q[0];
sx q[0];
rz(-0.5589267) q[0];
rz(2.3129297) q[1];
sx q[1];
rz(-1.5507973) q[1];
sx q[1];
rz(-0.25159803) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64564451) q[0];
sx q[0];
rz(-0.48316607) q[0];
sx q[0];
rz(0.42074411) q[0];
rz(1.0423123) q[2];
sx q[2];
rz(-2.5229075) q[2];
sx q[2];
rz(-1.2109735) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4797552) q[1];
sx q[1];
rz(-2.2244456) q[1];
sx q[1];
rz(0.33262555) q[1];
rz(-pi) q[2];
rz(2.2280327) q[3];
sx q[3];
rz(-2.1373478) q[3];
sx q[3];
rz(-0.65627894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.33637968) q[2];
sx q[2];
rz(-2.7153335) q[2];
sx q[2];
rz(-1.4317929) q[2];
rz(-0.21670565) q[3];
sx q[3];
rz(-1.5313989) q[3];
sx q[3];
rz(1.7823035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25877229) q[0];
sx q[0];
rz(-0.50676218) q[0];
sx q[0];
rz(1.850542) q[0];
rz(0.82921118) q[1];
sx q[1];
rz(-2.0334838) q[1];
sx q[1];
rz(2.1270027) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1433346) q[0];
sx q[0];
rz(-0.59276544) q[0];
sx q[0];
rz(0.44500764) q[0];
rz(-pi) q[1];
rz(2.5780771) q[2];
sx q[2];
rz(-1.2402099) q[2];
sx q[2];
rz(-2.3553768) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.18568745) q[1];
sx q[1];
rz(-0.57480156) q[1];
sx q[1];
rz(-1.5622524) q[1];
rz(0.15249522) q[3];
sx q[3];
rz(-2.3806981) q[3];
sx q[3];
rz(-0.12013474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0097051) q[2];
sx q[2];
rz(-1.9405126) q[2];
sx q[2];
rz(-0.76060549) q[2];
rz(0.46918121) q[3];
sx q[3];
rz(-1.4696591) q[3];
sx q[3];
rz(0.73452264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0483265) q[0];
sx q[0];
rz(-2.6272197) q[0];
sx q[0];
rz(2.560428) q[0];
rz(-1.5900853) q[1];
sx q[1];
rz(-2.4947512) q[1];
sx q[1];
rz(-0.69492984) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7030554) q[0];
sx q[0];
rz(-1.3378654) q[0];
sx q[0];
rz(-1.8700048) q[0];
rz(-1.476711) q[2];
sx q[2];
rz(-1.2796937) q[2];
sx q[2];
rz(1.7895607) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.123053) q[1];
sx q[1];
rz(-0.54529069) q[1];
sx q[1];
rz(-3.0569949) q[1];
x q[2];
rz(-0.65443734) q[3];
sx q[3];
rz(-0.29464606) q[3];
sx q[3];
rz(-1.0279931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.86357) q[2];
sx q[2];
rz(-1.5333971) q[2];
sx q[2];
rz(0.30280534) q[2];
rz(1.7152202) q[3];
sx q[3];
rz(-0.36077603) q[3];
sx q[3];
rz(2.5572131) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4645709) q[0];
sx q[0];
rz(-2.7320221) q[0];
sx q[0];
rz(2.4466178) q[0];
rz(2.8547844) q[1];
sx q[1];
rz(-2.5359055) q[1];
sx q[1];
rz(1.8668176) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5999278) q[0];
sx q[0];
rz(-2.554691) q[0];
sx q[0];
rz(0.98595263) q[0];
rz(1.1026682) q[2];
sx q[2];
rz(-2.063836) q[2];
sx q[2];
rz(1.4725034) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7949603) q[1];
sx q[1];
rz(-0.92466507) q[1];
sx q[1];
rz(2.6898795) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5596095) q[3];
sx q[3];
rz(-1.7714388) q[3];
sx q[3];
rz(-1.1021922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0810762) q[2];
sx q[2];
rz(-2.923107) q[2];
sx q[2];
rz(-2.1742353) q[2];
rz(-0.71990144) q[3];
sx q[3];
rz(-0.94341174) q[3];
sx q[3];
rz(1.2758183) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28441456) q[0];
sx q[0];
rz(-3.1107749) q[0];
sx q[0];
rz(1.6946633) q[0];
rz(-0.46724304) q[1];
sx q[1];
rz(-1.8485565) q[1];
sx q[1];
rz(-1.1955059) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23654437) q[0];
sx q[0];
rz(-1.1130739) q[0];
sx q[0];
rz(2.6603856) q[0];
rz(-2.823916) q[2];
sx q[2];
rz(-2.3001849) q[2];
sx q[2];
rz(2.9175772) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.171716) q[1];
sx q[1];
rz(-0.50285733) q[1];
sx q[1];
rz(1.4420532) q[1];
x q[2];
rz(1.5348263) q[3];
sx q[3];
rz(-1.5580873) q[3];
sx q[3];
rz(-2.980924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3168827) q[2];
sx q[2];
rz(-0.83493817) q[2];
sx q[2];
rz(0.52428025) q[2];
rz(2.8892062) q[3];
sx q[3];
rz(-0.10656825) q[3];
sx q[3];
rz(3.0805123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97279945) q[0];
sx q[0];
rz(-1.1058829) q[0];
sx q[0];
rz(-1.9148781) q[0];
rz(-2.3854158) q[1];
sx q[1];
rz(-2.1935479) q[1];
sx q[1];
rz(-2.7511168) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9928241) q[0];
sx q[0];
rz(-0.52339866) q[0];
sx q[0];
rz(2.5958909) q[0];
rz(-pi) q[1];
x q[1];
rz(0.080733732) q[2];
sx q[2];
rz(-0.44534031) q[2];
sx q[2];
rz(-0.75917086) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.272534) q[1];
sx q[1];
rz(-2.7400937) q[1];
sx q[1];
rz(-2.1436006) q[1];
rz(-pi) q[2];
x q[2];
rz(2.212575) q[3];
sx q[3];
rz(-1.9784156) q[3];
sx q[3];
rz(2.3315725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8172424) q[2];
sx q[2];
rz(-2.4411026) q[2];
sx q[2];
rz(1.3520757) q[2];
rz(-0.17524854) q[3];
sx q[3];
rz(-1.3388355) q[3];
sx q[3];
rz(-1.2043183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7290203) q[0];
sx q[0];
rz(-2.6798798) q[0];
sx q[0];
rz(-0.66226688) q[0];
rz(-0.63016713) q[1];
sx q[1];
rz(-1.8616734) q[1];
sx q[1];
rz(2.9060649) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61155897) q[0];
sx q[0];
rz(-1.3970102) q[0];
sx q[0];
rz(-0.21428971) q[0];
rz(-pi) q[1];
x q[1];
rz(1.477219) q[2];
sx q[2];
rz(-2.1380599) q[2];
sx q[2];
rz(-0.34514375) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.77068608) q[1];
sx q[1];
rz(-1.5944434) q[1];
sx q[1];
rz(1.6321983) q[1];
x q[2];
rz(-1.730758) q[3];
sx q[3];
rz(-0.70288881) q[3];
sx q[3];
rz(0.95332586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0887289) q[2];
sx q[2];
rz(-2.1910618) q[2];
sx q[2];
rz(-2.5843184) q[2];
rz(-2.3142464) q[3];
sx q[3];
rz(-1.5118303) q[3];
sx q[3];
rz(-1.6415589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92852229) q[0];
sx q[0];
rz(-2.3445573) q[0];
sx q[0];
rz(-2.573977) q[0];
rz(-0.40191832) q[1];
sx q[1];
rz(-2.4979976) q[1];
sx q[1];
rz(-1.7793122) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0716996) q[0];
sx q[0];
rz(-1.9687459) q[0];
sx q[0];
rz(-2.0073298) q[0];
rz(-pi) q[1];
rz(2.1902607) q[2];
sx q[2];
rz(-1.7835604) q[2];
sx q[2];
rz(1.484953) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6475923) q[1];
sx q[1];
rz(-0.63557445) q[1];
sx q[1];
rz(1.7032436) q[1];
rz(-pi) q[2];
rz(-2.3188842) q[3];
sx q[3];
rz(-2.1629982) q[3];
sx q[3];
rz(-1.4440735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7207429) q[2];
sx q[2];
rz(-1.2030615) q[2];
sx q[2];
rz(-2.8813349) q[2];
rz(-0.69089729) q[3];
sx q[3];
rz(-0.8553718) q[3];
sx q[3];
rz(-1.5232085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8061168) q[0];
sx q[0];
rz(-1.5609043) q[0];
sx q[0];
rz(1.5326395) q[0];
rz(0.43329049) q[1];
sx q[1];
rz(-1.3229803) q[1];
sx q[1];
rz(1.9569474) q[1];
rz(0.717502) q[2];
sx q[2];
rz(-2.2926471) q[2];
sx q[2];
rz(-2.9903464) q[2];
rz(2.3512202) q[3];
sx q[3];
rz(-1.6675259) q[3];
sx q[3];
rz(-0.70601757) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
