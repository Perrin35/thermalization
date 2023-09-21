OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.82233959) q[0];
sx q[0];
rz(-0.95681325) q[0];
sx q[0];
rz(1.4332888) q[0];
rz(2.9046471) q[1];
sx q[1];
rz(-1.2304996) q[1];
sx q[1];
rz(-1.9967611) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.618453) q[0];
sx q[0];
rz(-1.3292399) q[0];
sx q[0];
rz(3.0735344) q[0];
rz(1.1510017) q[2];
sx q[2];
rz(-1.8325873) q[2];
sx q[2];
rz(2.7671438) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.634234) q[1];
sx q[1];
rz(-0.45295742) q[1];
sx q[1];
rz(1.6139612) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2727229) q[3];
sx q[3];
rz(-1.0469336) q[3];
sx q[3];
rz(-3.1071172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0248489) q[2];
sx q[2];
rz(-1.8086834) q[2];
sx q[2];
rz(1.9072745) q[2];
rz(2.0862789) q[3];
sx q[3];
rz(-2.0827259) q[3];
sx q[3];
rz(-1.9887259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18594436) q[0];
sx q[0];
rz(-1.1971104) q[0];
sx q[0];
rz(0.82988513) q[0];
rz(2.8886967) q[1];
sx q[1];
rz(-0.77650944) q[1];
sx q[1];
rz(0.84709644) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0704437) q[0];
sx q[0];
rz(-1.1304454) q[0];
sx q[0];
rz(-1.1308934) q[0];
rz(-0.77913021) q[2];
sx q[2];
rz(-1.7841633) q[2];
sx q[2];
rz(1.2029635) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6716799) q[1];
sx q[1];
rz(-2.5922853) q[1];
sx q[1];
rz(3.0104907) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.133425) q[3];
sx q[3];
rz(-1.2542033) q[3];
sx q[3];
rz(0.40382995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.117924) q[2];
sx q[2];
rz(-0.37332049) q[2];
sx q[2];
rz(2.4222597) q[2];
rz(-1.8524648) q[3];
sx q[3];
rz(-0.23580655) q[3];
sx q[3];
rz(1.6524564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8711202) q[0];
sx q[0];
rz(-1.8914762) q[0];
sx q[0];
rz(0.57139325) q[0];
rz(2.1748623) q[1];
sx q[1];
rz(-1.8833908) q[1];
sx q[1];
rz(2.6142696) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8910599) q[0];
sx q[0];
rz(-2.8322304) q[0];
sx q[0];
rz(-3.1191349) q[0];
rz(-pi) q[1];
rz(-0.6730404) q[2];
sx q[2];
rz(-1.3125784) q[2];
sx q[2];
rz(2.8845127) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1151162) q[1];
sx q[1];
rz(-2.3079702) q[1];
sx q[1];
rz(1.923418) q[1];
x q[2];
rz(2.7819488) q[3];
sx q[3];
rz(-2.375964) q[3];
sx q[3];
rz(-1.5639203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8335235) q[2];
sx q[2];
rz(-1.4493194) q[2];
sx q[2];
rz(-0.72675881) q[2];
rz(-0.99749342) q[3];
sx q[3];
rz(-2.5163429) q[3];
sx q[3];
rz(-1.6231026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-2.8961287) q[0];
sx q[0];
rz(-1.4828869) q[0];
sx q[0];
rz(-1.0035275) q[0];
rz(-3.1009122) q[1];
sx q[1];
rz(-2.0122806) q[1];
sx q[1];
rz(-0.8262659) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6711133) q[0];
sx q[0];
rz(-0.90958909) q[0];
sx q[0];
rz(0.83755042) q[0];
rz(-pi) q[1];
x q[1];
rz(0.063886558) q[2];
sx q[2];
rz(-0.70110828) q[2];
sx q[2];
rz(-2.137616) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0127718) q[1];
sx q[1];
rz(-1.9073309) q[1];
sx q[1];
rz(-3.0613042) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1816979) q[3];
sx q[3];
rz(-2.1660921) q[3];
sx q[3];
rz(2.1967595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5740009) q[2];
sx q[2];
rz(-1.8969798) q[2];
sx q[2];
rz(0.91439247) q[2];
rz(-0.10522035) q[3];
sx q[3];
rz(-1.5172232) q[3];
sx q[3];
rz(0.58925327) q[3];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88702622) q[0];
sx q[0];
rz(-1.5019324) q[0];
sx q[0];
rz(1.1337093) q[0];
rz(0.0056313593) q[1];
sx q[1];
rz(-1.4375604) q[1];
sx q[1];
rz(2.5240135) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3591101) q[0];
sx q[0];
rz(-1.5251274) q[0];
sx q[0];
rz(-3.111905) q[0];
rz(-pi) q[1];
x q[1];
rz(0.64230201) q[2];
sx q[2];
rz(-1.2080492) q[2];
sx q[2];
rz(-1.158266) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4434112) q[1];
sx q[1];
rz(-0.65084208) q[1];
sx q[1];
rz(1.0990259) q[1];
rz(-pi) q[2];
x q[2];
rz(0.17554749) q[3];
sx q[3];
rz(-1.5834337) q[3];
sx q[3];
rz(2.4059084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3474943) q[2];
sx q[2];
rz(-1.8811052) q[2];
sx q[2];
rz(-0.23362544) q[2];
rz(0.99308333) q[3];
sx q[3];
rz(-2.1916094) q[3];
sx q[3];
rz(-0.63794199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1722906) q[0];
sx q[0];
rz(-3.0811221) q[0];
sx q[0];
rz(-0.0897952) q[0];
rz(2.2604997) q[1];
sx q[1];
rz(-2.6125364) q[1];
sx q[1];
rz(-2.9615013) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0952547) q[0];
sx q[0];
rz(-2.1349847) q[0];
sx q[0];
rz(-0.5395704) q[0];
rz(-2.9526268) q[2];
sx q[2];
rz(-2.1660888) q[2];
sx q[2];
rz(1.300786) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.60702) q[1];
sx q[1];
rz(-1.8587451) q[1];
sx q[1];
rz(-1.7400017) q[1];
x q[2];
rz(0.018499231) q[3];
sx q[3];
rz(-2.0634994) q[3];
sx q[3];
rz(0.25445709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.2580516) q[2];
sx q[2];
rz(-1.5641314) q[2];
sx q[2];
rz(1.1759261) q[2];
rz(-0.073143395) q[3];
sx q[3];
rz(-0.72435838) q[3];
sx q[3];
rz(0.49889645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.183855) q[0];
sx q[0];
rz(-0.65559214) q[0];
sx q[0];
rz(1.4181597) q[0];
rz(-1.1649959) q[1];
sx q[1];
rz(-2.5823451) q[1];
sx q[1];
rz(-3.1013536) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7559164) q[0];
sx q[0];
rz(-1.809281) q[0];
sx q[0];
rz(0.083906108) q[0];
rz(-pi) q[1];
rz(-3.0954103) q[2];
sx q[2];
rz(-1.1650411) q[2];
sx q[2];
rz(1.6398167) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.38899598) q[1];
sx q[1];
rz(-2.1850359) q[1];
sx q[1];
rz(-2.9700301) q[1];
x q[2];
rz(-1.2651029) q[3];
sx q[3];
rz(-1.2672658) q[3];
sx q[3];
rz(0.85206735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0030901) q[2];
sx q[2];
rz(-2.3562727) q[2];
sx q[2];
rz(-2.9677532) q[2];
rz(-1.396817) q[3];
sx q[3];
rz(-1.4649748) q[3];
sx q[3];
rz(0.85913908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1180856) q[0];
sx q[0];
rz(-2.9243587) q[0];
sx q[0];
rz(1.404495) q[0];
rz(1.2069758) q[1];
sx q[1];
rz(-1.4852306) q[1];
sx q[1];
rz(1.5054024) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36639402) q[0];
sx q[0];
rz(-1.1530071) q[0];
sx q[0];
rz(1.226107) q[0];
rz(2.2340441) q[2];
sx q[2];
rz(-2.2758256) q[2];
sx q[2];
rz(2.7570587) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7561188) q[1];
sx q[1];
rz(-1.5382775) q[1];
sx q[1];
rz(-0.69617747) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1721341) q[3];
sx q[3];
rz(-1.2974032) q[3];
sx q[3];
rz(-0.45575842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.76190844) q[2];
sx q[2];
rz(-0.48019335) q[2];
sx q[2];
rz(-2.9910679) q[2];
rz(-1.6020417) q[3];
sx q[3];
rz(-1.9463041) q[3];
sx q[3];
rz(-2.5185744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
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
rz(-0.49939108) q[1];
sx q[1];
rz(-1.4285409) q[1];
sx q[1];
rz(-1.2649149) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1160693) q[0];
sx q[0];
rz(-2.2795296) q[0];
sx q[0];
rz(2.619333) q[0];
rz(2.5876797) q[2];
sx q[2];
rz(-1.1914807) q[2];
sx q[2];
rz(-1.5253138) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7281108) q[1];
sx q[1];
rz(-0.15004798) q[1];
sx q[1];
rz(2.8996182) q[1];
rz(-2.7154269) q[3];
sx q[3];
rz(-0.45939547) q[3];
sx q[3];
rz(-0.029749425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.50679961) q[2];
sx q[2];
rz(-0.96968499) q[2];
sx q[2];
rz(-0.52465049) q[2];
rz(2.5850885) q[3];
sx q[3];
rz(-1.5126712) q[3];
sx q[3];
rz(2.7867253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5234579) q[0];
sx q[0];
rz(-0.83704656) q[0];
sx q[0];
rz(2.8961704) q[0];
rz(0.55631176) q[1];
sx q[1];
rz(-1.5309155) q[1];
sx q[1];
rz(1.6773178) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9368162) q[0];
sx q[0];
rz(-1.4602772) q[0];
sx q[0];
rz(0.21893455) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9173031) q[2];
sx q[2];
rz(-0.99594342) q[2];
sx q[2];
rz(-1.526051) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6357395) q[1];
sx q[1];
rz(-0.74843279) q[1];
sx q[1];
rz(1.1670508) q[1];
rz(1.7467935) q[3];
sx q[3];
rz(-1.4868075) q[3];
sx q[3];
rz(1.849911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2422553) q[2];
sx q[2];
rz(-1.238404) q[2];
sx q[2];
rz(-0.74550068) q[2];
rz(-1.4238822) q[3];
sx q[3];
rz(-2.8427377) q[3];
sx q[3];
rz(2.0132813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-1.8650919) q[0];
sx q[0];
rz(-1.6456974) q[0];
sx q[0];
rz(0.67281848) q[0];
rz(-1.8854234) q[1];
sx q[1];
rz(-0.80614631) q[1];
sx q[1];
rz(2.0731906) q[1];
rz(0.55143572) q[2];
sx q[2];
rz(-2.3497992) q[2];
sx q[2];
rz(-2.1247911) q[2];
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