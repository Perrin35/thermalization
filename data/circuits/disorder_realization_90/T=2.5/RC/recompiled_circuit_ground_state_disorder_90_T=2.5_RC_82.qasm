OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6381792) q[0];
sx q[0];
rz(-1.8621651) q[0];
sx q[0];
rz(2.3655565) q[0];
rz(-2.4322721) q[1];
sx q[1];
rz(-1.5172989) q[1];
sx q[1];
rz(0.59901839) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8255492) q[0];
sx q[0];
rz(-2.6514396) q[0];
sx q[0];
rz(1.8328299) q[0];
x q[1];
rz(-0.28223306) q[2];
sx q[2];
rz(-0.90105173) q[2];
sx q[2];
rz(-1.2362267) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6769961) q[1];
sx q[1];
rz(-2.0647767) q[1];
sx q[1];
rz(1.9736273) q[1];
x q[2];
rz(-2.7077449) q[3];
sx q[3];
rz(-2.3126855) q[3];
sx q[3];
rz(-0.69218194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1385931) q[2];
sx q[2];
rz(-2.004576) q[2];
sx q[2];
rz(0.19621672) q[2];
rz(1.044322) q[3];
sx q[3];
rz(-0.82572562) q[3];
sx q[3];
rz(-0.31061068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0394548) q[0];
sx q[0];
rz(-2.4882443) q[0];
sx q[0];
rz(2.2838604) q[0];
rz(2.6013382) q[1];
sx q[1];
rz(-2.0876355) q[1];
sx q[1];
rz(2.3589755) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1038541) q[0];
sx q[0];
rz(-1.4840992) q[0];
sx q[0];
rz(-2.187665) q[0];
x q[1];
rz(-1.247632) q[2];
sx q[2];
rz(-0.64715451) q[2];
sx q[2];
rz(-2.0314856) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9036213) q[1];
sx q[1];
rz(-0.97717932) q[1];
sx q[1];
rz(1.3376544) q[1];
x q[2];
rz(2.2563062) q[3];
sx q[3];
rz(-2.0205704) q[3];
sx q[3];
rz(1.7974082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.1648078) q[2];
sx q[2];
rz(-0.79443496) q[2];
sx q[2];
rz(2.5505193) q[2];
rz(-2.1137386) q[3];
sx q[3];
rz(-0.63575345) q[3];
sx q[3];
rz(3.0052321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22981055) q[0];
sx q[0];
rz(-2.6214143) q[0];
sx q[0];
rz(-2.0853364) q[0];
rz(-2.1504869) q[1];
sx q[1];
rz(-1.3737563) q[1];
sx q[1];
rz(2.3371005) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4933518) q[0];
sx q[0];
rz(-1.8364779) q[0];
sx q[0];
rz(0.65346395) q[0];
rz(-pi) q[1];
rz(-1.7011004) q[2];
sx q[2];
rz(-0.63281239) q[2];
sx q[2];
rz(-1.1820584) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.04136297) q[1];
sx q[1];
rz(-1.3303262) q[1];
sx q[1];
rz(-2.9036456) q[1];
rz(0.64739703) q[3];
sx q[3];
rz(-1.7248099) q[3];
sx q[3];
rz(-1.9141509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6911917) q[2];
sx q[2];
rz(-2.2914026) q[2];
sx q[2];
rz(-2.5737393) q[2];
rz(-1.19207) q[3];
sx q[3];
rz(-0.53693938) q[3];
sx q[3];
rz(0.15538628) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.49517) q[0];
sx q[0];
rz(-1.3986873) q[0];
sx q[0];
rz(1.1871185) q[0];
rz(-2.8861956) q[1];
sx q[1];
rz(-1.7385769) q[1];
sx q[1];
rz(0.38937169) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0856649) q[0];
sx q[0];
rz(-1.4367668) q[0];
sx q[0];
rz(0.22179752) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.537048) q[2];
sx q[2];
rz(-0.72740388) q[2];
sx q[2];
rz(-0.45798618) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7070605) q[1];
sx q[1];
rz(-1.7599306) q[1];
sx q[1];
rz(2.2347336) q[1];
rz(1.7674753) q[3];
sx q[3];
rz(-1.4405319) q[3];
sx q[3];
rz(-1.2127339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7369467) q[2];
sx q[2];
rz(-0.32361042) q[2];
sx q[2];
rz(1.7741989) q[2];
rz(-2.6020452) q[3];
sx q[3];
rz(-1.5522141) q[3];
sx q[3];
rz(-1.3246983) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0668199) q[0];
sx q[0];
rz(-2.9434151) q[0];
sx q[0];
rz(-2.5812126) q[0];
rz(0.21891521) q[1];
sx q[1];
rz(-2.0603265) q[1];
sx q[1];
rz(2.970649) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8290299) q[0];
sx q[0];
rz(-1.1033022) q[0];
sx q[0];
rz(-0.74689052) q[0];
rz(1.4291228) q[2];
sx q[2];
rz(-1.7729086) q[2];
sx q[2];
rz(1.0712445) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5569607) q[1];
sx q[1];
rz(-2.0804394) q[1];
sx q[1];
rz(0.72344785) q[1];
rz(-pi) q[2];
rz(1.6690977) q[3];
sx q[3];
rz(-0.62255961) q[3];
sx q[3];
rz(-0.05427256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7818452) q[2];
sx q[2];
rz(-2.013194) q[2];
sx q[2];
rz(2.183059) q[2];
rz(2.9790699) q[3];
sx q[3];
rz(-0.84238094) q[3];
sx q[3];
rz(-1.5925647) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.573134) q[0];
sx q[0];
rz(-3.077226) q[0];
sx q[0];
rz(2.3934613) q[0];
rz(2.023078) q[1];
sx q[1];
rz(-2.7225814) q[1];
sx q[1];
rz(-2.8245139) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90077315) q[0];
sx q[0];
rz(-2.8860807) q[0];
sx q[0];
rz(-0.27285646) q[0];
x q[1];
rz(-2.9002764) q[2];
sx q[2];
rz(-2.9394627) q[2];
sx q[2];
rz(1.8692819) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.78817716) q[1];
sx q[1];
rz(-1.4995575) q[1];
sx q[1];
rz(-0.80435462) q[1];
rz(-pi) q[2];
rz(3.1012332) q[3];
sx q[3];
rz(-0.77733126) q[3];
sx q[3];
rz(-0.07829994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.76724425) q[2];
sx q[2];
rz(-0.92864645) q[2];
sx q[2];
rz(1.413215) q[2];
rz(3.0610541) q[3];
sx q[3];
rz(-1.3486515) q[3];
sx q[3];
rz(2.9197689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0291075) q[0];
sx q[0];
rz(-1.8928098) q[0];
sx q[0];
rz(1.8674194) q[0];
rz(1.3451276) q[1];
sx q[1];
rz(-0.74256623) q[1];
sx q[1];
rz(-1.3412195) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2224436) q[0];
sx q[0];
rz(-2.0608927) q[0];
sx q[0];
rz(0.63700139) q[0];
rz(-pi) q[1];
x q[1];
rz(0.51711443) q[2];
sx q[2];
rz(-0.54143751) q[2];
sx q[2];
rz(2.0064029) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.53645027) q[1];
sx q[1];
rz(-1.4172232) q[1];
sx q[1];
rz(-1.3773514) q[1];
rz(-pi) q[2];
rz(0.11734924) q[3];
sx q[3];
rz(-0.92399358) q[3];
sx q[3];
rz(-0.58157677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.93366569) q[2];
sx q[2];
rz(-1.1873446) q[2];
sx q[2];
rz(-2.6742317) q[2];
rz(-0.76006877) q[3];
sx q[3];
rz(-1.4773388) q[3];
sx q[3];
rz(-1.2490341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46045983) q[0];
sx q[0];
rz(-1.0204027) q[0];
sx q[0];
rz(-1.4469752) q[0];
rz(0.32577062) q[1];
sx q[1];
rz(-0.21369801) q[1];
sx q[1];
rz(1.5209341) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50624079) q[0];
sx q[0];
rz(-1.304848) q[0];
sx q[0];
rz(-2.6145934) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.09163945) q[2];
sx q[2];
rz(-1.5704201) q[2];
sx q[2];
rz(-0.74972744) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7829166) q[1];
sx q[1];
rz(-0.39729983) q[1];
sx q[1];
rz(2.0372169) q[1];
rz(-pi) q[2];
rz(1.8160024) q[3];
sx q[3];
rz(-0.71906656) q[3];
sx q[3];
rz(-1.2326425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1763566) q[2];
sx q[2];
rz(-0.74702817) q[2];
sx q[2];
rz(-0.61526862) q[2];
rz(1.8687013) q[3];
sx q[3];
rz(-0.26510173) q[3];
sx q[3];
rz(-2.7191775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(-1.0203005) q[0];
sx q[0];
rz(-1.2739807) q[0];
sx q[0];
rz(-2.2913388) q[0];
rz(-2.0203159) q[1];
sx q[1];
rz(-1.3669776) q[1];
sx q[1];
rz(2.3698295) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61311713) q[0];
sx q[0];
rz(-0.69418797) q[0];
sx q[0];
rz(-1.5334237) q[0];
x q[1];
rz(-2.775101) q[2];
sx q[2];
rz(-2.9482609) q[2];
sx q[2];
rz(1.1761348) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5215251) q[1];
sx q[1];
rz(-1.8142533) q[1];
sx q[1];
rz(-3.0848178) q[1];
x q[2];
rz(-0.87987514) q[3];
sx q[3];
rz(-0.41112435) q[3];
sx q[3];
rz(-0.087669186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0261592) q[2];
sx q[2];
rz(-1.5966281) q[2];
sx q[2];
rz(2.3010632) q[2];
rz(-0.080951512) q[3];
sx q[3];
rz(-0.5778802) q[3];
sx q[3];
rz(-2.8988083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42846546) q[0];
sx q[0];
rz(-2.9436538) q[0];
sx q[0];
rz(1.9848829) q[0];
rz(2.1139862) q[1];
sx q[1];
rz(-1.7637858) q[1];
sx q[1];
rz(-0.81046945) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9543957) q[0];
sx q[0];
rz(-0.82747298) q[0];
sx q[0];
rz(-1.8660956) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9058305) q[2];
sx q[2];
rz(-1.050569) q[2];
sx q[2];
rz(-0.65641415) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1406882) q[1];
sx q[1];
rz(-1.8546805) q[1];
sx q[1];
rz(-1.8769916) q[1];
rz(-pi) q[2];
rz(1.719172) q[3];
sx q[3];
rz(-2.1273489) q[3];
sx q[3];
rz(-0.30239964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.79601866) q[2];
sx q[2];
rz(-1.2556262) q[2];
sx q[2];
rz(-1.5239117) q[2];
rz(2.0412622) q[3];
sx q[3];
rz(-1.1805781) q[3];
sx q[3];
rz(0.17980096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1179467) q[0];
sx q[0];
rz(-1.7272341) q[0];
sx q[0];
rz(-0.9247307) q[0];
rz(-3.0991411) q[1];
sx q[1];
rz(-1.1141384) q[1];
sx q[1];
rz(1.3420807) q[1];
rz(-1.4412075) q[2];
sx q[2];
rz(-2.7322506) q[2];
sx q[2];
rz(2.6553497) q[2];
rz(2.6841738) q[3];
sx q[3];
rz(-1.538496) q[3];
sx q[3];
rz(-2.9067007) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
