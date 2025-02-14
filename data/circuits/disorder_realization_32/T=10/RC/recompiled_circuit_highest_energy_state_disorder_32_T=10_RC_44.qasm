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
rz(-1.929317) q[0];
sx q[0];
rz(-1.1317929) q[0];
sx q[0];
rz(1.7687067) q[0];
rz(2.0139439) q[1];
sx q[1];
rz(4.5166587) q[1];
sx q[1];
rz(5.4969129) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36390135) q[0];
sx q[0];
rz(-0.4253952) q[0];
sx q[0];
rz(2.0869654) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2271871) q[2];
sx q[2];
rz(-2.0255651) q[2];
sx q[2];
rz(2.9637873) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8189803) q[1];
sx q[1];
rz(-2.8507345) q[1];
sx q[1];
rz(-1.6751272) q[1];
rz(-pi) q[2];
rz(-0.97162515) q[3];
sx q[3];
rz(-1.9135981) q[3];
sx q[3];
rz(2.9620217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0218411) q[2];
sx q[2];
rz(-1.4234875) q[2];
sx q[2];
rz(-1.2500259) q[2];
rz(2.405808) q[3];
sx q[3];
rz(-0.17446987) q[3];
sx q[3];
rz(-1.638394) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4982872) q[0];
sx q[0];
rz(-1.9685638) q[0];
sx q[0];
rz(-0.071320891) q[0];
rz(1.1530863) q[1];
sx q[1];
rz(-0.76140296) q[1];
sx q[1];
rz(2.6128795) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77069717) q[0];
sx q[0];
rz(-2.5636854) q[0];
sx q[0];
rz(0.96630567) q[0];
rz(-pi) q[1];
rz(-2.9527093) q[2];
sx q[2];
rz(-1.3391277) q[2];
sx q[2];
rz(2.846039) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.15859998) q[1];
sx q[1];
rz(-2.0865177) q[1];
sx q[1];
rz(0.92293585) q[1];
rz(-1.1558481) q[3];
sx q[3];
rz(-1.7768716) q[3];
sx q[3];
rz(-1.2690488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4802287) q[2];
sx q[2];
rz(-2.7429548) q[2];
sx q[2];
rz(3.1129692) q[2];
rz(-2.7875767) q[3];
sx q[3];
rz(-0.75494868) q[3];
sx q[3];
rz(2.5825175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70327586) q[0];
sx q[0];
rz(-2.1385312) q[0];
sx q[0];
rz(0.89135528) q[0];
rz(-0.50093961) q[1];
sx q[1];
rz(-1.5762065) q[1];
sx q[1];
rz(-0.058813728) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0100493) q[0];
sx q[0];
rz(-1.6185404) q[0];
sx q[0];
rz(0.086011925) q[0];
x q[1];
rz(0.28532736) q[2];
sx q[2];
rz(-1.0614191) q[2];
sx q[2];
rz(-0.90778186) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2477639) q[1];
sx q[1];
rz(-2.2732452) q[1];
sx q[1];
rz(0.6599627) q[1];
rz(-pi) q[2];
rz(-2.9887588) q[3];
sx q[3];
rz(-2.9590769) q[3];
sx q[3];
rz(-2.07978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7848876) q[2];
sx q[2];
rz(-0.55861837) q[2];
sx q[2];
rz(0.2365665) q[2];
rz(0.92075721) q[3];
sx q[3];
rz(-1.0509138) q[3];
sx q[3];
rz(1.7304272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22223602) q[0];
sx q[0];
rz(-1.8574497) q[0];
sx q[0];
rz(2.2663569) q[0];
rz(-1.5454166) q[1];
sx q[1];
rz(-0.86219209) q[1];
sx q[1];
rz(3.0962871) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52462477) q[0];
sx q[0];
rz(-3.0201206) q[0];
sx q[0];
rz(1.0018906) q[0];
rz(-pi) q[1];
rz(-2.2384134) q[2];
sx q[2];
rz(-1.982769) q[2];
sx q[2];
rz(2.2087086) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.873535) q[1];
sx q[1];
rz(-2.2378439) q[1];
sx q[1];
rz(-0.11063834) q[1];
rz(-0.1940345) q[3];
sx q[3];
rz(-2.1450419) q[3];
sx q[3];
rz(1.6358835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8594325) q[2];
sx q[2];
rz(-2.1789447) q[2];
sx q[2];
rz(-1.947594) q[2];
rz(0.89030877) q[3];
sx q[3];
rz(-1.2229536) q[3];
sx q[3];
rz(-2.1931026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0328261) q[0];
sx q[0];
rz(-0.81619167) q[0];
sx q[0];
rz(2.4591675) q[0];
rz(1.2884864) q[1];
sx q[1];
rz(-1.5623743) q[1];
sx q[1];
rz(1.8103745) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3159861) q[0];
sx q[0];
rz(-2.2974112) q[0];
sx q[0];
rz(3.1196703) q[0];
rz(-1.9648026) q[2];
sx q[2];
rz(-1.3298243) q[2];
sx q[2];
rz(0.35277982) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9115548) q[1];
sx q[1];
rz(-1.4304407) q[1];
sx q[1];
rz(-2.7799003) q[1];
x q[2];
rz(0.45068897) q[3];
sx q[3];
rz(-2.0540049) q[3];
sx q[3];
rz(-0.86408981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1149301) q[2];
sx q[2];
rz(-0.59938359) q[2];
sx q[2];
rz(0.66693532) q[2];
rz(-0.55495787) q[3];
sx q[3];
rz(-0.68735492) q[3];
sx q[3];
rz(2.8426898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8126467) q[0];
sx q[0];
rz(-1.8829367) q[0];
sx q[0];
rz(0.70372787) q[0];
rz(2.9723889) q[1];
sx q[1];
rz(-1.5237944) q[1];
sx q[1];
rz(0.15883787) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.279351) q[0];
sx q[0];
rz(-2.9547174) q[0];
sx q[0];
rz(-0.86743768) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6936982) q[2];
sx q[2];
rz(-1.4180866) q[2];
sx q[2];
rz(-2.0982775) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.006268) q[1];
sx q[1];
rz(-0.49990955) q[1];
sx q[1];
rz(-2.8630524) q[1];
rz(-pi) q[2];
x q[2];
rz(0.63780906) q[3];
sx q[3];
rz(-1.7553698) q[3];
sx q[3];
rz(2.0710433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.77330294) q[2];
sx q[2];
rz(-2.1998019) q[2];
sx q[2];
rz(0.85757315) q[2];
rz(-1.1693906) q[3];
sx q[3];
rz(-1.8604859) q[3];
sx q[3];
rz(-2.9331971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17539772) q[0];
sx q[0];
rz(-0.76681391) q[0];
sx q[0];
rz(-2.4229557) q[0];
rz(2.8596558) q[1];
sx q[1];
rz(-1.7666631) q[1];
sx q[1];
rz(-0.66863543) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0575585) q[0];
sx q[0];
rz(-2.4120122) q[0];
sx q[0];
rz(1.0490824) q[0];
rz(-pi) q[1];
x q[1];
rz(0.35923194) q[2];
sx q[2];
rz(-0.97373913) q[2];
sx q[2];
rz(0.12166858) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9397647) q[1];
sx q[1];
rz(-1.9376905) q[1];
sx q[1];
rz(0.35122996) q[1];
x q[2];
rz(2.5103522) q[3];
sx q[3];
rz(-1.5641912) q[3];
sx q[3];
rz(-0.94830482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4947027) q[2];
sx q[2];
rz(-1.1855482) q[2];
sx q[2];
rz(-2.1051712) q[2];
rz(-2.5070665) q[3];
sx q[3];
rz(-2.1536638) q[3];
sx q[3];
rz(0.79276597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6805639) q[0];
sx q[0];
rz(-0.11756086) q[0];
sx q[0];
rz(0.72599894) q[0];
rz(1.1622102) q[1];
sx q[1];
rz(-1.1770959) q[1];
sx q[1];
rz(1.9451709) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0045753) q[0];
sx q[0];
rz(-1.5423623) q[0];
sx q[0];
rz(1.7445376) q[0];
rz(-pi) q[1];
rz(3.093946) q[2];
sx q[2];
rz(-1.8777913) q[2];
sx q[2];
rz(-2.7757211) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6641621) q[1];
sx q[1];
rz(-0.71738418) q[1];
sx q[1];
rz(0.75023164) q[1];
rz(2.7128588) q[3];
sx q[3];
rz(-2.4247243) q[3];
sx q[3];
rz(0.69566876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.90765816) q[2];
sx q[2];
rz(-1.5479167) q[2];
sx q[2];
rz(-0.93351239) q[2];
rz(0.61698169) q[3];
sx q[3];
rz(-1.0022997) q[3];
sx q[3];
rz(-0.84701076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9762978) q[0];
sx q[0];
rz(-3.0369861) q[0];
sx q[0];
rz(-3.0883375) q[0];
rz(0.28757295) q[1];
sx q[1];
rz(-0.92929274) q[1];
sx q[1];
rz(0.25921777) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0236146) q[0];
sx q[0];
rz(-0.058246944) q[0];
sx q[0];
rz(1.9969166) q[0];
x q[1];
rz(-1.0561814) q[2];
sx q[2];
rz(-0.26517235) q[2];
sx q[2];
rz(2.3725703) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.99154918) q[1];
sx q[1];
rz(-2.3880312) q[1];
sx q[1];
rz(-1.6607453) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7227371) q[3];
sx q[3];
rz(-0.81063089) q[3];
sx q[3];
rz(-1.6068271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.58664924) q[2];
sx q[2];
rz(-1.8879075) q[2];
sx q[2];
rz(2.9602236) q[2];
rz(-0.3950611) q[3];
sx q[3];
rz(-2.919988) q[3];
sx q[3];
rz(-2.1612397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7419389) q[0];
sx q[0];
rz(-1.548865) q[0];
sx q[0];
rz(-2.7594866) q[0];
rz(-2.8451879) q[1];
sx q[1];
rz(-1.0045241) q[1];
sx q[1];
rz(1.872725) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73055501) q[0];
sx q[0];
rz(-2.0417622) q[0];
sx q[0];
rz(2.6593326) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2440654) q[2];
sx q[2];
rz(-1.466905) q[2];
sx q[2];
rz(-1.2938703) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8386744) q[1];
sx q[1];
rz(-1.5494635) q[1];
sx q[1];
rz(-1.231074) q[1];
rz(-2.6878854) q[3];
sx q[3];
rz(-1.4611011) q[3];
sx q[3];
rz(-0.060893313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2962013) q[2];
sx q[2];
rz(-0.56075823) q[2];
sx q[2];
rz(-0.38983795) q[2];
rz(-1.2975533) q[3];
sx q[3];
rz(-1.9997948) q[3];
sx q[3];
rz(-2.2977184) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1525477) q[0];
sx q[0];
rz(-1.4023517) q[0];
sx q[0];
rz(-1.5040816) q[0];
rz(-1.7851495) q[1];
sx q[1];
rz(-2.103613) q[1];
sx q[1];
rz(1.6115859) q[1];
rz(0.98567166) q[2];
sx q[2];
rz(-2.4092842) q[2];
sx q[2];
rz(0.53436919) q[2];
rz(0.33792321) q[3];
sx q[3];
rz(-2.3281964) q[3];
sx q[3];
rz(-2.7184814) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
