OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5175944) q[0];
sx q[0];
rz(-2.6156293) q[0];
sx q[0];
rz(2.9232803) q[0];
rz(-1.6649618) q[1];
sx q[1];
rz(-0.47774878) q[1];
sx q[1];
rz(0.55396095) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71611878) q[0];
sx q[0];
rz(-2.0332575) q[0];
sx q[0];
rz(1.0008971) q[0];
rz(-pi) q[1];
rz(3.0402115) q[2];
sx q[2];
rz(-2.6104365) q[2];
sx q[2];
rz(0.44975933) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.40383717) q[1];
sx q[1];
rz(-0.74450508) q[1];
sx q[1];
rz(1.8831403) q[1];
rz(-pi) q[2];
rz(-1.9945108) q[3];
sx q[3];
rz(-2.0689575) q[3];
sx q[3];
rz(-2.8776282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7654984) q[2];
sx q[2];
rz(-0.29355294) q[2];
sx q[2];
rz(-1.2676839) q[2];
rz(2.3119161) q[3];
sx q[3];
rz(-1.6206348) q[3];
sx q[3];
rz(0.42384306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.053452881) q[0];
sx q[0];
rz(-1.3351048) q[0];
sx q[0];
rz(-1.394519) q[0];
rz(-2.246619) q[1];
sx q[1];
rz(-1.0286237) q[1];
sx q[1];
rz(1.5590394) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6532063) q[0];
sx q[0];
rz(-2.1108642) q[0];
sx q[0];
rz(-2.3652943) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1505914) q[2];
sx q[2];
rz(-1.8981427) q[2];
sx q[2];
rz(0.98132747) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8445721) q[1];
sx q[1];
rz(-0.75648844) q[1];
sx q[1];
rz(-2.6761961) q[1];
x q[2];
rz(-0.30541916) q[3];
sx q[3];
rz(-1.1608184) q[3];
sx q[3];
rz(2.3596606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.47722694) q[2];
sx q[2];
rz(-1.619728) q[2];
sx q[2];
rz(-3.1108372) q[2];
rz(2.6886046) q[3];
sx q[3];
rz(-2.9121297) q[3];
sx q[3];
rz(-0.17791137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40193108) q[0];
sx q[0];
rz(-3.0406096) q[0];
sx q[0];
rz(-2.3133551) q[0];
rz(3.0939057) q[1];
sx q[1];
rz(-2.2779155) q[1];
sx q[1];
rz(1.9140859) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0704884) q[0];
sx q[0];
rz(-1.0336735) q[0];
sx q[0];
rz(2.0893196) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.21593185) q[2];
sx q[2];
rz(-1.9683483) q[2];
sx q[2];
rz(2.2550607) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8611698) q[1];
sx q[1];
rz(-2.6211328) q[1];
sx q[1];
rz(2.1153482) q[1];
x q[2];
rz(-0.84215409) q[3];
sx q[3];
rz(-1.531736) q[3];
sx q[3];
rz(-1.0674764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0487655) q[2];
sx q[2];
rz(-2.2266677) q[2];
sx q[2];
rz(-1.1009334) q[2];
rz(0.01384211) q[3];
sx q[3];
rz(-1.3661386) q[3];
sx q[3];
rz(-2.3069416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58532995) q[0];
sx q[0];
rz(-1.7381373) q[0];
sx q[0];
rz(-2.8072667) q[0];
rz(-0.73257929) q[1];
sx q[1];
rz(-2.1926447) q[1];
sx q[1];
rz(3.0308731) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3363655) q[0];
sx q[0];
rz(-2.3799161) q[0];
sx q[0];
rz(-0.021999981) q[0];
rz(-pi) q[1];
x q[1];
rz(0.51038536) q[2];
sx q[2];
rz(-2.9523627) q[2];
sx q[2];
rz(-0.58923474) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.91023129) q[1];
sx q[1];
rz(-1.0139324) q[1];
sx q[1];
rz(-0.75480144) q[1];
rz(-pi) q[2];
rz(-0.48583123) q[3];
sx q[3];
rz(-0.83668349) q[3];
sx q[3];
rz(-2.776718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.683814) q[2];
sx q[2];
rz(-0.28202287) q[2];
sx q[2];
rz(-1.4078183) q[2];
rz(1.1553361) q[3];
sx q[3];
rz(-1.1563533) q[3];
sx q[3];
rz(2.9467764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69283501) q[0];
sx q[0];
rz(-0.91788569) q[0];
sx q[0];
rz(-2.1441929) q[0];
rz(-1.6150486) q[1];
sx q[1];
rz(-2.5020182) q[1];
sx q[1];
rz(1.4195199) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90516312) q[0];
sx q[0];
rz(-1.7337203) q[0];
sx q[0];
rz(-1.7927367) q[0];
rz(-1.7260584) q[2];
sx q[2];
rz(-2.3387032) q[2];
sx q[2];
rz(0.93198317) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.75052634) q[1];
sx q[1];
rz(-1.8896034) q[1];
sx q[1];
rz(-0.90853779) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0180198) q[3];
sx q[3];
rz(-1.6032748) q[3];
sx q[3];
rz(0.44636727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.851696) q[2];
sx q[2];
rz(-0.09446129) q[2];
sx q[2];
rz(-0.32290253) q[2];
rz(-2.0276535) q[3];
sx q[3];
rz(-1.0909922) q[3];
sx q[3];
rz(2.7285301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7776529) q[0];
sx q[0];
rz(-2.8033065) q[0];
sx q[0];
rz(2.3349578) q[0];
rz(-0.58397645) q[1];
sx q[1];
rz(-2.0239425) q[1];
sx q[1];
rz(1.1899828) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34790137) q[0];
sx q[0];
rz(-0.37037235) q[0];
sx q[0];
rz(1.9166975) q[0];
rz(0.34752589) q[2];
sx q[2];
rz(-2.579877) q[2];
sx q[2];
rz(-0.39081854) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8998225) q[1];
sx q[1];
rz(-2.9651151) q[1];
sx q[1];
rz(1.8679669) q[1];
rz(-1.4011995) q[3];
sx q[3];
rz(-2.2940594) q[3];
sx q[3];
rz(-2.4270428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.40818647) q[2];
sx q[2];
rz(-1.6383645) q[2];
sx q[2];
rz(-2.9564986) q[2];
rz(1.5271651) q[3];
sx q[3];
rz(-1.4027169) q[3];
sx q[3];
rz(-2.4534498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9516893) q[0];
sx q[0];
rz(-0.95174319) q[0];
sx q[0];
rz(0.21251799) q[0];
rz(2.3566133) q[1];
sx q[1];
rz(-1.6467983) q[1];
sx q[1];
rz(2.3627538) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76701421) q[0];
sx q[0];
rz(-0.85142577) q[0];
sx q[0];
rz(-2.3423751) q[0];
rz(2.5539407) q[2];
sx q[2];
rz(-1.8853123) q[2];
sx q[2];
rz(1.8314198) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0743588) q[1];
sx q[1];
rz(-1.7869084) q[1];
sx q[1];
rz(-1.8317779) q[1];
rz(-pi) q[2];
rz(2.1206585) q[3];
sx q[3];
rz(-2.3519197) q[3];
sx q[3];
rz(-2.0146973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6276041) q[2];
sx q[2];
rz(-2.6476761) q[2];
sx q[2];
rz(-0.40840515) q[2];
rz(2.4397395) q[3];
sx q[3];
rz(-2.2097094) q[3];
sx q[3];
rz(1.8823889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7931165) q[0];
sx q[0];
rz(-1.2154673) q[0];
sx q[0];
rz(-3.075573) q[0];
rz(-1.5090212) q[1];
sx q[1];
rz(-1.0450109) q[1];
sx q[1];
rz(-0.95796934) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5292266) q[0];
sx q[0];
rz(-2.5924396) q[0];
sx q[0];
rz(-2.5757838) q[0];
x q[1];
rz(-2.1195378) q[2];
sx q[2];
rz(-2.3081452) q[2];
sx q[2];
rz(1.1481783) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9878475) q[1];
sx q[1];
rz(-1.2586795) q[1];
sx q[1];
rz(-0.69674833) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1470966) q[3];
sx q[3];
rz(-2.9589783) q[3];
sx q[3];
rz(-1.2507358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3079188) q[2];
sx q[2];
rz(-1.6180399) q[2];
sx q[2];
rz(-1.7306805) q[2];
rz(-0.69502568) q[3];
sx q[3];
rz(-1.6618988) q[3];
sx q[3];
rz(3.1237349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0787635) q[0];
sx q[0];
rz(-1.6595027) q[0];
sx q[0];
rz(-2.1516946) q[0];
rz(0.46317378) q[1];
sx q[1];
rz(-1.4521234) q[1];
sx q[1];
rz(1.2407726) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8340943) q[0];
sx q[0];
rz(-1.8655329) q[0];
sx q[0];
rz(0.060262279) q[0];
rz(2.0133205) q[2];
sx q[2];
rz(-1.754369) q[2];
sx q[2];
rz(1.5300446) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.040673469) q[1];
sx q[1];
rz(-1.2881491) q[1];
sx q[1];
rz(-1.2987178) q[1];
x q[2];
rz(-0.45205493) q[3];
sx q[3];
rz(-0.99957672) q[3];
sx q[3];
rz(-2.6719928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3339633) q[2];
sx q[2];
rz(-1.7619851) q[2];
sx q[2];
rz(-2.4228952) q[2];
rz(1.9888318) q[3];
sx q[3];
rz(-1.7199687) q[3];
sx q[3];
rz(-1.0144455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54866791) q[0];
sx q[0];
rz(-0.34807006) q[0];
sx q[0];
rz(-2.3396709) q[0];
rz(-1.0879263) q[1];
sx q[1];
rz(-1.3366924) q[1];
sx q[1];
rz(0.71802872) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.021026) q[0];
sx q[0];
rz(-2.2533198) q[0];
sx q[0];
rz(-1.1521856) q[0];
x q[1];
rz(0.52387107) q[2];
sx q[2];
rz(-1.1351368) q[2];
sx q[2];
rz(-0.73843971) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.57783571) q[1];
sx q[1];
rz(-2.8056762) q[1];
sx q[1];
rz(1.4116686) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4122756) q[3];
sx q[3];
rz(-0.98776796) q[3];
sx q[3];
rz(2.5860525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3512909) q[2];
sx q[2];
rz(-0.21778926) q[2];
sx q[2];
rz(2.6289319) q[2];
rz(-0.74448186) q[3];
sx q[3];
rz(-1.0267886) q[3];
sx q[3];
rz(0.7640394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9164593) q[0];
sx q[0];
rz(-1.8712578) q[0];
sx q[0];
rz(1.901392) q[0];
rz(-0.71612877) q[1];
sx q[1];
rz(-0.54812535) q[1];
sx q[1];
rz(-2.5352238) q[1];
rz(-3.0307583) q[2];
sx q[2];
rz(-1.7946984) q[2];
sx q[2];
rz(2.6738965) q[2];
rz(-1.8567139) q[3];
sx q[3];
rz(-0.34964041) q[3];
sx q[3];
rz(1.6007363) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
