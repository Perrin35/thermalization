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
rz(-2.1051536) q[0];
sx q[0];
rz(1.5358465) q[0];
rz(-2.4802471) q[1];
sx q[1];
rz(-1.9547434) q[1];
sx q[1];
rz(-0.43651906) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.614994) q[0];
sx q[0];
rz(-2.8078571) q[0];
sx q[0];
rz(-1.2646535) q[0];
rz(-2.0233193) q[2];
sx q[2];
rz(-1.0529537) q[2];
sx q[2];
rz(-2.81942) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.84387842) q[1];
sx q[1];
rz(-1.6116287) q[1];
sx q[1];
rz(-0.86754129) q[1];
rz(-pi) q[2];
rz(-2.2546886) q[3];
sx q[3];
rz(-1.4330079) q[3];
sx q[3];
rz(0.52007127) q[3];
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
rz(-3.0758514) q[3];
sx q[3];
rz(-1.3275423) q[3];
sx q[3];
rz(-1.340723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0098669212) q[0];
sx q[0];
rz(-1.0634402) q[0];
sx q[0];
rz(-3.104082) q[0];
rz(-0.58049479) q[1];
sx q[1];
rz(-0.66941222) q[1];
sx q[1];
rz(-2.4494749) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6094309) q[0];
sx q[0];
rz(-1.2464644) q[0];
sx q[0];
rz(-0.18478365) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0742998) q[2];
sx q[2];
rz(-1.8177336) q[2];
sx q[2];
rz(1.4114789) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.925732) q[1];
sx q[1];
rz(-2.070148) q[1];
sx q[1];
rz(2.4768922) q[1];
rz(0.18499891) q[3];
sx q[3];
rz(-0.76055148) q[3];
sx q[3];
rz(2.4724677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.47824255) q[2];
sx q[2];
rz(-1.9448091) q[2];
sx q[2];
rz(1.7459858) q[2];
rz(-2.9546402) q[3];
sx q[3];
rz(-1.9713277) q[3];
sx q[3];
rz(-0.54005867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-0.66889399) q[0];
sx q[0];
rz(-0.63303328) q[0];
sx q[0];
rz(2.582666) q[0];
rz(-0.82866296) q[1];
sx q[1];
rz(-1.5507973) q[1];
sx q[1];
rz(-0.25159803) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54788816) q[0];
sx q[0];
rz(-1.761709) q[0];
sx q[0];
rz(2.6949791) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0195658) q[2];
sx q[2];
rz(-1.2740268) q[2];
sx q[2];
rz(0.084159764) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.66183749) q[1];
sx q[1];
rz(-2.2244456) q[1];
sx q[1];
rz(0.33262555) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.76520701) q[3];
sx q[3];
rz(-2.3022989) q[3];
sx q[3];
rz(-2.8347903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.805213) q[2];
sx q[2];
rz(-0.42625913) q[2];
sx q[2];
rz(1.4317929) q[2];
rz(-0.21670565) q[3];
sx q[3];
rz(-1.5313989) q[3];
sx q[3];
rz(-1.3592892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8828204) q[0];
sx q[0];
rz(-2.6348305) q[0];
sx q[0];
rz(1.850542) q[0];
rz(-0.82921118) q[1];
sx q[1];
rz(-1.1081089) q[1];
sx q[1];
rz(2.1270027) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9982581) q[0];
sx q[0];
rz(-0.59276544) q[0];
sx q[0];
rz(-2.696585) q[0];
x q[1];
rz(-1.185174) q[2];
sx q[2];
rz(-1.0411388) q[2];
sx q[2];
rz(-0.5822863) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9559052) q[1];
sx q[1];
rz(-2.5667911) q[1];
sx q[1];
rz(1.5793403) q[1];
x q[2];
rz(0.15249522) q[3];
sx q[3];
rz(-2.3806981) q[3];
sx q[3];
rz(-0.12013474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0097051) q[2];
sx q[2];
rz(-1.9405126) q[2];
sx q[2];
rz(0.76060549) q[2];
rz(0.46918121) q[3];
sx q[3];
rz(-1.6719336) q[3];
sx q[3];
rz(2.40707) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0932662) q[0];
sx q[0];
rz(-2.6272197) q[0];
sx q[0];
rz(-0.58116466) q[0];
rz(-1.5515074) q[1];
sx q[1];
rz(-2.4947512) q[1];
sx q[1];
rz(-2.4466628) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43853727) q[0];
sx q[0];
rz(-1.8037272) q[0];
sx q[0];
rz(1.2715879) q[0];
x q[1];
rz(0.29232358) q[2];
sx q[2];
rz(-1.6609123) q[2];
sx q[2];
rz(-2.9499049) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.1174005) q[1];
sx q[1];
rz(-2.1139189) q[1];
sx q[1];
rz(1.622011) q[1];
x q[2];
rz(-1.7534689) q[3];
sx q[3];
rz(-1.338306) q[3];
sx q[3];
rz(-2.7894265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.27802262) q[2];
sx q[2];
rz(-1.6081955) q[2];
sx q[2];
rz(0.30280534) q[2];
rz(1.4263724) q[3];
sx q[3];
rz(-0.36077603) q[3];
sx q[3];
rz(0.58437955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67702174) q[0];
sx q[0];
rz(-0.40957054) q[0];
sx q[0];
rz(0.69497481) q[0];
rz(0.28680828) q[1];
sx q[1];
rz(-2.5359055) q[1];
sx q[1];
rz(-1.8668176) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0114637) q[0];
sx q[0];
rz(-2.0507567) q[0];
sx q[0];
rz(-2.7897054) q[0];
rz(-pi) q[1];
rz(-1.1026682) q[2];
sx q[2];
rz(-1.0777567) q[2];
sx q[2];
rz(1.4725034) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6332153) q[1];
sx q[1];
rz(-1.214809) q[1];
sx q[1];
rz(0.87320019) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.58198317) q[3];
sx q[3];
rz(-1.3701539) q[3];
sx q[3];
rz(2.0394005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0605165) q[2];
sx q[2];
rz(-0.21848564) q[2];
sx q[2];
rz(2.1742353) q[2];
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
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28441456) q[0];
sx q[0];
rz(-3.1107749) q[0];
sx q[0];
rz(-1.4469294) q[0];
rz(0.46724304) q[1];
sx q[1];
rz(-1.8485565) q[1];
sx q[1];
rz(-1.9460868) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0366096) q[0];
sx q[0];
rz(-0.65147841) q[0];
sx q[0];
rz(-2.3250513) q[0];
rz(-pi) q[1];
rz(-0.81584064) q[2];
sx q[2];
rz(-1.3357478) q[2];
sx q[2];
rz(1.5624969) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.96987665) q[1];
sx q[1];
rz(-0.50285733) q[1];
sx q[1];
rz(-1.4420532) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5348263) q[3];
sx q[3];
rz(-1.5580873) q[3];
sx q[3];
rz(2.980924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3168827) q[2];
sx q[2];
rz(-0.83493817) q[2];
sx q[2];
rz(-0.52428025) q[2];
rz(-2.8892062) q[3];
sx q[3];
rz(-3.0350244) q[3];
sx q[3];
rz(-0.061080385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97279945) q[0];
sx q[0];
rz(-2.0357098) q[0];
sx q[0];
rz(1.9148781) q[0];
rz(-0.75617689) q[1];
sx q[1];
rz(-0.94804472) q[1];
sx q[1];
rz(-2.7511168) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9928241) q[0];
sx q[0];
rz(-2.618194) q[0];
sx q[0];
rz(-0.54570178) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5323213) q[2];
sx q[2];
rz(-2.0145825) q[2];
sx q[2];
rz(2.2930068) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8836964) q[1];
sx q[1];
rz(-1.9054277) q[1];
sx q[1];
rz(-2.9154214) q[1];
rz(-pi) q[2];
rz(-0.92901765) q[3];
sx q[3];
rz(-1.1631771) q[3];
sx q[3];
rz(-2.3315725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.32435027) q[2];
sx q[2];
rz(-0.70049006) q[2];
sx q[2];
rz(1.7895169) q[2];
rz(0.17524854) q[3];
sx q[3];
rz(-1.8027571) q[3];
sx q[3];
rz(1.9372743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4125724) q[0];
sx q[0];
rz(-2.6798798) q[0];
sx q[0];
rz(-2.4793258) q[0];
rz(0.63016713) q[1];
sx q[1];
rz(-1.8616734) q[1];
sx q[1];
rz(-2.9060649) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6306289) q[0];
sx q[0];
rz(-0.27505829) q[0];
sx q[0];
rz(-0.69013005) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9959683) q[2];
sx q[2];
rz(-0.57409414) q[2];
sx q[2];
rz(0.51806322) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9743226) q[1];
sx q[1];
rz(-3.0758) q[1];
sx q[1];
rz(1.9386824) q[1];
rz(-pi) q[2];
rz(0.8742378) q[3];
sx q[3];
rz(-1.6739419) q[3];
sx q[3];
rz(0.49498765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0528637) q[2];
sx q[2];
rz(-0.95053089) q[2];
sx q[2];
rz(0.55727422) q[2];
rz(2.3142464) q[3];
sx q[3];
rz(-1.6297623) q[3];
sx q[3];
rz(1.5000337) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2130704) q[0];
sx q[0];
rz(-0.7970354) q[0];
sx q[0];
rz(-2.573977) q[0];
rz(-0.40191832) q[1];
sx q[1];
rz(-0.64359507) q[1];
sx q[1];
rz(-1.3622805) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6797736) q[0];
sx q[0];
rz(-1.1704233) q[0];
sx q[0];
rz(-2.7072565) q[0];
rz(-1.2145813) q[2];
sx q[2];
rz(-0.65041861) q[2];
sx q[2];
rz(0.37373558) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6475923) q[1];
sx q[1];
rz(-2.5060182) q[1];
sx q[1];
rz(-1.7032436) q[1];
x q[2];
rz(2.3990223) q[3];
sx q[3];
rz(-0.97108632) q[3];
sx q[3];
rz(2.7903583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7207429) q[2];
sx q[2];
rz(-1.2030615) q[2];
sx q[2];
rz(0.26025772) q[2];
rz(-0.69089729) q[3];
sx q[3];
rz(-2.2862209) q[3];
sx q[3];
rz(1.5232085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33547587) q[0];
sx q[0];
rz(-1.5806883) q[0];
sx q[0];
rz(-1.6089532) q[0];
rz(0.43329049) q[1];
sx q[1];
rz(-1.3229803) q[1];
sx q[1];
rz(1.9569474) q[1];
rz(2.4337089) q[2];
sx q[2];
rz(-1.0546726) q[2];
sx q[2];
rz(2.2451014) q[2];
rz(1.4337486) q[3];
sx q[3];
rz(-0.78513405) q[3];
sx q[3];
rz(0.76754192) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
