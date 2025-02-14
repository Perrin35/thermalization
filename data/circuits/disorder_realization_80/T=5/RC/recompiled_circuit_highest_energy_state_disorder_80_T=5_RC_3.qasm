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
rz(0.39245519) q[0];
sx q[0];
rz(-0.28372228) q[0];
sx q[0];
rz(1.1878045) q[0];
rz(2.6999733) q[1];
sx q[1];
rz(-1.266357) q[1];
sx q[1];
rz(-1.1511572) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2806704) q[0];
sx q[0];
rz(-0.09083561) q[0];
sx q[0];
rz(-1.3263741) q[0];
rz(-1.2032937) q[2];
sx q[2];
rz(-1.1137059) q[2];
sx q[2];
rz(-0.45705308) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.77598876) q[1];
sx q[1];
rz(-1.9950486) q[1];
sx q[1];
rz(1.4308903) q[1];
x q[2];
rz(-2.0602442) q[3];
sx q[3];
rz(-0.57751211) q[3];
sx q[3];
rz(0.58693991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.89897951) q[2];
sx q[2];
rz(-0.8967163) q[2];
sx q[2];
rz(-2.8420281) q[2];
rz(-0.35563955) q[3];
sx q[3];
rz(-2.1483597) q[3];
sx q[3];
rz(-0.43660823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(-0.75293175) q[0];
sx q[0];
rz(-0.54838538) q[0];
sx q[0];
rz(-2.8652628) q[0];
rz(-0.07946864) q[1];
sx q[1];
rz(-0.53229585) q[1];
sx q[1];
rz(0.4963378) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2329335) q[0];
sx q[0];
rz(-2.56224) q[0];
sx q[0];
rz(2.0412372) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4935819) q[2];
sx q[2];
rz(-1.2334012) q[2];
sx q[2];
rz(1.5959306) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4071634) q[1];
sx q[1];
rz(-1.2447272) q[1];
sx q[1];
rz(2.3666275) q[1];
x q[2];
rz(-3.1383508) q[3];
sx q[3];
rz(-1.7326983) q[3];
sx q[3];
rz(1.1461794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.77208272) q[2];
sx q[2];
rz(-1.3630867) q[2];
sx q[2];
rz(2.962964) q[2];
rz(1.5313088) q[3];
sx q[3];
rz(-0.9261927) q[3];
sx q[3];
rz(-2.3791651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4302706) q[0];
sx q[0];
rz(-1.0866168) q[0];
sx q[0];
rz(1.7311199) q[0];
rz(-1.3871644) q[1];
sx q[1];
rz(-1.6491363) q[1];
sx q[1];
rz(-1.656146) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58830611) q[0];
sx q[0];
rz(-1.6035514) q[0];
sx q[0];
rz(2.4577228) q[0];
rz(-pi) q[1];
x q[1];
rz(0.40992592) q[2];
sx q[2];
rz(-2.2927444) q[2];
sx q[2];
rz(0.97463911) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3130114) q[1];
sx q[1];
rz(-1.269956) q[1];
sx q[1];
rz(-1.1238696) q[1];
rz(-pi) q[2];
rz(-2.0250506) q[3];
sx q[3];
rz(-2.4245533) q[3];
sx q[3];
rz(-1.3644791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9401271) q[2];
sx q[2];
rz(-2.1374233) q[2];
sx q[2];
rz(-2.6489769) q[2];
rz(2.7228739) q[3];
sx q[3];
rz(-1.0227572) q[3];
sx q[3];
rz(0.59188265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.017460499) q[0];
sx q[0];
rz(-5/(13*pi)) q[0];
sx q[0];
rz(0.19805743) q[0];
rz(-0.65912229) q[1];
sx q[1];
rz(-2.5810869) q[1];
sx q[1];
rz(-2.9445599) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7871088) q[0];
sx q[0];
rz(-1.0403883) q[0];
sx q[0];
rz(1.6097989) q[0];
rz(-pi) q[1];
rz(-0.23931673) q[2];
sx q[2];
rz(-1.2916358) q[2];
sx q[2];
rz(-2.5365732) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2391165) q[1];
sx q[1];
rz(-0.65985876) q[1];
sx q[1];
rz(0.59999864) q[1];
rz(-pi) q[2];
rz(-1.6313309) q[3];
sx q[3];
rz(-0.54577561) q[3];
sx q[3];
rz(2.7013403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.78843242) q[2];
sx q[2];
rz(-1.3015231) q[2];
sx q[2];
rz(1.8678467) q[2];
rz(-1.5984009) q[3];
sx q[3];
rz(-1.1929932) q[3];
sx q[3];
rz(0.25632349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0565599) q[0];
sx q[0];
rz(-1.7340478) q[0];
sx q[0];
rz(-0.10277596) q[0];
rz(0.13725266) q[1];
sx q[1];
rz(-0.44932258) q[1];
sx q[1];
rz(-1.4385673) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53046561) q[0];
sx q[0];
rz(-1.0004064) q[0];
sx q[0];
rz(0.80369759) q[0];
x q[1];
rz(-0.46311997) q[2];
sx q[2];
rz(-1.6119405) q[2];
sx q[2];
rz(2.9477313) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4754191) q[1];
sx q[1];
rz(-1.8295145) q[1];
sx q[1];
rz(-2.9889876) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6827312) q[3];
sx q[3];
rz(-1.4482968) q[3];
sx q[3];
rz(-1.4089597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.30456257) q[2];
sx q[2];
rz(-1.6670767) q[2];
sx q[2];
rz(2.4414818) q[2];
rz(0.89811283) q[3];
sx q[3];
rz(-0.63106314) q[3];
sx q[3];
rz(-1.496544) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20823088) q[0];
sx q[0];
rz(-0.58657402) q[0];
sx q[0];
rz(1.5775648) q[0];
rz(-0.6598407) q[1];
sx q[1];
rz(-2.3520062) q[1];
sx q[1];
rz(2.1671364) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3269269) q[0];
sx q[0];
rz(-1.5788933) q[0];
sx q[0];
rz(0.011122313) q[0];
x q[1];
rz(0.31807138) q[2];
sx q[2];
rz(-2.632395) q[2];
sx q[2];
rz(1.731002) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.240307) q[1];
sx q[1];
rz(-1.069233) q[1];
sx q[1];
rz(1.5585414) q[1];
rz(-1.5100689) q[3];
sx q[3];
rz(-2.1083045) q[3];
sx q[3];
rz(-2.2030061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4727891) q[2];
sx q[2];
rz(-0.69837022) q[2];
sx q[2];
rz(-0.93639708) q[2];
rz(-2.2514553) q[3];
sx q[3];
rz(-1.3606768) q[3];
sx q[3];
rz(2.9357125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-1.5278006) q[0];
sx q[0];
rz(-0.19659909) q[0];
sx q[0];
rz(3.0278681) q[0];
rz(-1.3971036) q[1];
sx q[1];
rz(-2.0753658) q[1];
sx q[1];
rz(0.024959175) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67201383) q[0];
sx q[0];
rz(-0.22525283) q[0];
sx q[0];
rz(1.8173056) q[0];
rz(-pi) q[1];
rz(1.3106842) q[2];
sx q[2];
rz(-0.94365135) q[2];
sx q[2];
rz(-0.24218923) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0164707) q[1];
sx q[1];
rz(-0.85713398) q[1];
sx q[1];
rz(-0.086507052) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9437441) q[3];
sx q[3];
rz(-1.6258844) q[3];
sx q[3];
rz(0.69869374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6156442) q[2];
sx q[2];
rz(-1.4333466) q[2];
sx q[2];
rz(-2.136266) q[2];
rz(2.2510236) q[3];
sx q[3];
rz(-1.4054479) q[3];
sx q[3];
rz(2.4864206) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6152076) q[0];
sx q[0];
rz(-0.40281519) q[0];
sx q[0];
rz(2.030754) q[0];
rz(3.0530744) q[1];
sx q[1];
rz(-2.7169777) q[1];
sx q[1];
rz(-1.7875338) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7294069) q[0];
sx q[0];
rz(-0.78347396) q[0];
sx q[0];
rz(-1.7307348) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8432003) q[2];
sx q[2];
rz(-1.7188311) q[2];
sx q[2];
rz(-1.4349951) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9971294) q[1];
sx q[1];
rz(-1.6465997) q[1];
sx q[1];
rz(2.8949225) q[1];
x q[2];
rz(2.5398387) q[3];
sx q[3];
rz(-2.2838998) q[3];
sx q[3];
rz(2.8181638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2985349) q[2];
sx q[2];
rz(-0.83710805) q[2];
sx q[2];
rz(2.7561772) q[2];
rz(-0.75000969) q[3];
sx q[3];
rz(-0.95732006) q[3];
sx q[3];
rz(0.065940417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0176004) q[0];
sx q[0];
rz(-2.0662859) q[0];
sx q[0];
rz(0.83220926) q[0];
rz(1.0819134) q[1];
sx q[1];
rz(-0.81580201) q[1];
sx q[1];
rz(-1.6197416) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8502626) q[0];
sx q[0];
rz(-1.1801774) q[0];
sx q[0];
rz(-1.0608835) q[0];
rz(-pi) q[1];
rz(2.3004901) q[2];
sx q[2];
rz(-0.37174598) q[2];
sx q[2];
rz(1.339762) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3491813) q[1];
sx q[1];
rz(-0.99073016) q[1];
sx q[1];
rz(-2.4522454) q[1];
rz(-pi) q[2];
x q[2];
rz(0.24116349) q[3];
sx q[3];
rz(-2.2134288) q[3];
sx q[3];
rz(1.3295028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1443783) q[2];
sx q[2];
rz(-1.8961366) q[2];
sx q[2];
rz(1.2388371) q[2];
rz(1.8501806) q[3];
sx q[3];
rz(-1.8700799) q[3];
sx q[3];
rz(3.0102357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.130403) q[0];
sx q[0];
rz(-0.953453) q[0];
sx q[0];
rz(0.74380547) q[0];
rz(3.0939843) q[1];
sx q[1];
rz(-2.036939) q[1];
sx q[1];
rz(1.1415175) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4726758) q[0];
sx q[0];
rz(-2.0341691) q[0];
sx q[0];
rz(-2.1457304) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9594934) q[2];
sx q[2];
rz(-2.1517589) q[2];
sx q[2];
rz(-2.3240391) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0298668) q[1];
sx q[1];
rz(-2.2891785) q[1];
sx q[1];
rz(0.54775441) q[1];
rz(-pi) q[2];
rz(-0.23903592) q[3];
sx q[3];
rz(-0.65863684) q[3];
sx q[3];
rz(-2.4995668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.49012524) q[2];
sx q[2];
rz(-2.9059124) q[2];
sx q[2];
rz(-2.8362595) q[2];
rz(-1.4281645) q[3];
sx q[3];
rz(-1.7686663) q[3];
sx q[3];
rz(1.292424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.74844985) q[0];
sx q[0];
rz(-0.23101692) q[0];
sx q[0];
rz(1.1363181) q[0];
rz(0.5592067) q[1];
sx q[1];
rz(-1.8103841) q[1];
sx q[1];
rz(0.24787535) q[1];
rz(1.2216907) q[2];
sx q[2];
rz(-0.93358718) q[2];
sx q[2];
rz(-2.6503218) q[2];
rz(-1.705966) q[3];
sx q[3];
rz(-1.4528989) q[3];
sx q[3];
rz(-1.0624878) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
