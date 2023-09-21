OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3192531) q[0];
sx q[0];
rz(-2.1847794) q[0];
sx q[0];
rz(1.7083038) q[0];
rz(-0.23694555) q[1];
sx q[1];
rz(-1.911093) q[1];
sx q[1];
rz(1.9967611) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5231397) q[0];
sx q[0];
rz(-1.3292399) q[0];
sx q[0];
rz(-0.068058204) q[0];
rz(-1.1510017) q[2];
sx q[2];
rz(-1.3090054) q[2];
sx q[2];
rz(-0.37444886) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1022542) q[1];
sx q[1];
rz(-1.5896817) q[1];
sx q[1];
rz(-1.1182055) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8688698) q[3];
sx q[3];
rz(-1.0469336) q[3];
sx q[3];
rz(-3.1071172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.11674374) q[2];
sx q[2];
rz(-1.8086834) q[2];
sx q[2];
rz(-1.2343181) q[2];
rz(-1.0553137) q[3];
sx q[3];
rz(-2.0827259) q[3];
sx q[3];
rz(1.1528667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18594436) q[0];
sx q[0];
rz(-1.9444822) q[0];
sx q[0];
rz(2.3117075) q[0];
rz(0.25289598) q[1];
sx q[1];
rz(-2.3650832) q[1];
sx q[1];
rz(0.84709644) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69764187) q[0];
sx q[0];
rz(-1.9662494) q[0];
sx q[0];
rz(2.6614499) q[0];
rz(-pi) q[1];
rz(-1.2752091) q[2];
sx q[2];
rz(-0.81381961) q[2];
sx q[2];
rz(0.57397599) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6716799) q[1];
sx q[1];
rz(-2.5922853) q[1];
sx q[1];
rz(-3.0104907) q[1];
x q[2];
rz(-1.133425) q[3];
sx q[3];
rz(-1.2542033) q[3];
sx q[3];
rz(-2.7377627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0236686) q[2];
sx q[2];
rz(-0.37332049) q[2];
sx q[2];
rz(0.71933293) q[2];
rz(1.2891278) q[3];
sx q[3];
rz(-0.23580655) q[3];
sx q[3];
rz(-1.4891362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
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
rz(1.8711202) q[0];
sx q[0];
rz(-1.2501165) q[0];
sx q[0];
rz(0.57139325) q[0];
rz(-2.1748623) q[1];
sx q[1];
rz(-1.8833908) q[1];
sx q[1];
rz(0.5273231) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29887154) q[0];
sx q[0];
rz(-1.5776331) q[0];
sx q[0];
rz(2.8323035) q[0];
rz(-pi) q[1];
rz(-0.40076077) q[2];
sx q[2];
rz(-2.4279865) q[2];
sx q[2];
rz(2.1378627) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.69818245) q[1];
sx q[1];
rz(-1.3122307) q[1];
sx q[1];
rz(-2.3727388) q[1];
rz(-pi) q[2];
rz(-1.2445883) q[3];
sx q[3];
rz(-0.86498125) q[3];
sx q[3];
rz(-1.0969485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8335235) q[2];
sx q[2];
rz(-1.4493194) q[2];
sx q[2];
rz(-2.4148338) q[2];
rz(-2.1440992) q[3];
sx q[3];
rz(-0.62524978) q[3];
sx q[3];
rz(1.5184901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24546394) q[0];
sx q[0];
rz(-1.6587057) q[0];
sx q[0];
rz(-2.1380651) q[0];
rz(3.1009122) q[1];
sx q[1];
rz(-2.0122806) q[1];
sx q[1];
rz(-2.3153268) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47047939) q[0];
sx q[0];
rz(-0.90958909) q[0];
sx q[0];
rz(-2.3040422) q[0];
x q[1];
rz(-0.063886558) q[2];
sx q[2];
rz(-0.70110828) q[2];
sx q[2];
rz(2.137616) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.367825) q[1];
sx q[1];
rz(-2.7959679) q[1];
sx q[1];
rz(-1.7961545) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1816979) q[3];
sx q[3];
rz(-0.97550052) q[3];
sx q[3];
rz(2.1967595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5740009) q[2];
sx q[2];
rz(-1.8969798) q[2];
sx q[2];
rz(2.2272002) q[2];
rz(-0.10522035) q[3];
sx q[3];
rz(-1.5172232) q[3];
sx q[3];
rz(0.58925327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88702622) q[0];
sx q[0];
rz(-1.6396602) q[0];
sx q[0];
rz(-1.1337093) q[0];
rz(-0.0056313593) q[1];
sx q[1];
rz(-1.7040323) q[1];
sx q[1];
rz(-0.61757913) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78695801) q[0];
sx q[0];
rz(-1.600453) q[0];
sx q[0];
rz(-1.6164854) q[0];
rz(-pi) q[1];
x q[1];
rz(0.64230201) q[2];
sx q[2];
rz(-1.9335434) q[2];
sx q[2];
rz(1.158266) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.51296556) q[1];
sx q[1];
rz(-1.8497397) q[1];
sx q[1];
rz(-2.1668424) q[1];
x q[2];
rz(-1.5836309) q[3];
sx q[3];
rz(-1.7463297) q[3];
sx q[3];
rz(0.83735355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3474943) q[2];
sx q[2];
rz(-1.8811052) q[2];
sx q[2];
rz(-2.9079672) q[2];
rz(0.99308333) q[3];
sx q[3];
rz(-2.1916094) q[3];
sx q[3];
rz(2.5036507) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1722906) q[0];
sx q[0];
rz(-3.0811221) q[0];
sx q[0];
rz(-3.0517975) q[0];
rz(0.88109294) q[1];
sx q[1];
rz(-2.6125364) q[1];
sx q[1];
rz(2.9615013) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0952547) q[0];
sx q[0];
rz(-2.1349847) q[0];
sx q[0];
rz(2.6020223) q[0];
rz(-pi) q[1];
rz(0.96713709) q[2];
sx q[2];
rz(-1.4146311) q[2];
sx q[2];
rz(0.16317633) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5345726) q[1];
sx q[1];
rz(-1.2828476) q[1];
sx q[1];
rz(-1.401591) q[1];
rz(-pi) q[2];
rz(-1.0780219) q[3];
sx q[3];
rz(-1.587095) q[3];
sx q[3];
rz(-1.8340045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.2580516) q[2];
sx q[2];
rz(-1.5774612) q[2];
sx q[2];
rz(1.1759261) q[2];
rz(-3.0684493) q[3];
sx q[3];
rz(-2.4172343) q[3];
sx q[3];
rz(-2.6426962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.183855) q[0];
sx q[0];
rz(-0.65559214) q[0];
sx q[0];
rz(1.7234329) q[0];
rz(1.9765967) q[1];
sx q[1];
rz(-2.5823451) q[1];
sx q[1];
rz(-3.1013536) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0979472) q[0];
sx q[0];
rz(-2.8890434) q[0];
sx q[0];
rz(1.9027684) q[0];
rz(-pi) q[1];
rz(-0.046182403) q[2];
sx q[2];
rz(-1.9765515) q[2];
sx q[2];
rz(-1.501776) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.38899598) q[1];
sx q[1];
rz(-2.1850359) q[1];
sx q[1];
rz(-2.9700301) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8242565) q[3];
sx q[3];
rz(-1.279497) q[3];
sx q[3];
rz(-2.3288162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0030901) q[2];
sx q[2];
rz(-2.3562727) q[2];
sx q[2];
rz(2.9677532) q[2];
rz(1.396817) q[3];
sx q[3];
rz(-1.4649748) q[3];
sx q[3];
rz(-0.85913908) q[3];
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
rz(-2.023507) q[0];
sx q[0];
rz(-2.9243587) q[0];
sx q[0];
rz(-1.404495) q[0];
rz(-1.2069758) q[1];
sx q[1];
rz(-1.6563621) q[1];
sx q[1];
rz(-1.6361902) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0907785) q[0];
sx q[0];
rz(-2.6065126) q[0];
sx q[0];
rz(0.65061609) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5152399) q[2];
sx q[2];
rz(-0.92712958) q[2];
sx q[2];
rz(-0.49382892) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.38547388) q[1];
sx q[1];
rz(-1.6033152) q[1];
sx q[1];
rz(0.69617747) q[1];
rz(-pi) q[2];
x q[2];
rz(0.32780111) q[3];
sx q[3];
rz(-0.9947239) q[3];
sx q[3];
rz(2.2097261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.76190844) q[2];
sx q[2];
rz(-0.48019335) q[2];
sx q[2];
rz(-2.9910679) q[2];
rz(-1.5395509) q[3];
sx q[3];
rz(-1.1952885) q[3];
sx q[3];
rz(-2.5185744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4432916) q[0];
sx q[0];
rz(-2.0069831) q[0];
sx q[0];
rz(-2.495893) q[0];
rz(-2.6422016) q[1];
sx q[1];
rz(-1.4285409) q[1];
sx q[1];
rz(1.2649149) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0255233) q[0];
sx q[0];
rz(-2.2795296) q[0];
sx q[0];
rz(0.52225964) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0090946) q[2];
sx q[2];
rz(-1.0602789) q[2];
sx q[2];
rz(-0.27062182) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7281108) q[1];
sx q[1];
rz(-2.9915447) q[1];
sx q[1];
rz(2.8996182) q[1];
rz(-pi) q[2];
rz(2.7154269) q[3];
sx q[3];
rz(-2.6821972) q[3];
sx q[3];
rz(-0.029749425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.50679961) q[2];
sx q[2];
rz(-2.1719077) q[2];
sx q[2];
rz(-0.52465049) q[2];
rz(-0.55650416) q[3];
sx q[3];
rz(-1.5126712) q[3];
sx q[3];
rz(2.7867253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-1.5234579) q[0];
sx q[0];
rz(-2.3045461) q[0];
sx q[0];
rz(0.24542228) q[0];
rz(-2.5852809) q[1];
sx q[1];
rz(-1.5309155) q[1];
sx q[1];
rz(-1.4642749) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20477644) q[0];
sx q[0];
rz(-1.6813155) q[0];
sx q[0];
rz(-2.9226581) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.98430888) q[2];
sx q[2];
rz(-1.7585635) q[2];
sx q[2];
rz(3.0629326) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7732685) q[1];
sx q[1];
rz(-1.3001633) q[1];
sx q[1];
rz(-2.2776105) q[1];
rz(1.3947992) q[3];
sx q[3];
rz(-1.6547852) q[3];
sx q[3];
rz(1.849911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2422553) q[2];
sx q[2];
rz(-1.9031886) q[2];
sx q[2];
rz(-2.396092) q[2];
rz(1.7177104) q[3];
sx q[3];
rz(-0.29885492) q[3];
sx q[3];
rz(-2.0132813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2765008) q[0];
sx q[0];
rz(-1.4958953) q[0];
sx q[0];
rz(-2.4687742) q[0];
rz(-1.2561692) q[1];
sx q[1];
rz(-2.3354463) q[1];
sx q[1];
rz(-1.068402) q[1];
rz(-2.0586661) q[2];
sx q[2];
rz(-2.2219873) q[2];
sx q[2];
rz(1.7359003) q[2];
rz(1.9477378) q[3];
sx q[3];
rz(-1.0300763) q[3];
sx q[3];
rz(-1.1618617) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
