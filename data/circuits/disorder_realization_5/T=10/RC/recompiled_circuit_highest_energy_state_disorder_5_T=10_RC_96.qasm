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
rz(-1.9755826) q[0];
sx q[0];
rz(-2.2628885) q[0];
sx q[0];
rz(2.9455844) q[0];
rz(-1.4930383) q[1];
sx q[1];
rz(-0.9382481) q[1];
sx q[1];
rz(0.87747639) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24841079) q[0];
sx q[0];
rz(-2.0201004) q[0];
sx q[0];
rz(2.9025159) q[0];
x q[1];
rz(2.7174905) q[2];
sx q[2];
rz(-1.8908894) q[2];
sx q[2];
rz(0.84148511) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6266105) q[1];
sx q[1];
rz(-2.0143854) q[1];
sx q[1];
rz(-2.842998) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0095413) q[3];
sx q[3];
rz(-0.51902229) q[3];
sx q[3];
rz(-1.2281017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.33285546) q[2];
sx q[2];
rz(-0.70781195) q[2];
sx q[2];
rz(3.1044002) q[2];
rz(-0.91341364) q[3];
sx q[3];
rz(-2.1593058) q[3];
sx q[3];
rz(-1.6325525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-0.34441242) q[0];
sx q[0];
rz(-2.0672815) q[0];
sx q[0];
rz(2.9194226) q[0];
rz(0.19368681) q[1];
sx q[1];
rz(-1.8733459) q[1];
sx q[1];
rz(0.2598612) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7714789) q[0];
sx q[0];
rz(-0.51379097) q[0];
sx q[0];
rz(-1.0593848) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1104576) q[2];
sx q[2];
rz(-0.20763131) q[2];
sx q[2];
rz(-2.5467444) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4187964) q[1];
sx q[1];
rz(-1.7823311) q[1];
sx q[1];
rz(1.3822894) q[1];
rz(2.3823795) q[3];
sx q[3];
rz(-1.9369643) q[3];
sx q[3];
rz(-3.1389583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.76693177) q[2];
sx q[2];
rz(-2.5587475) q[2];
sx q[2];
rz(-0.15677162) q[2];
rz(0.80592704) q[3];
sx q[3];
rz(-2.0511878) q[3];
sx q[3];
rz(2.6167615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2369775) q[0];
sx q[0];
rz(-3.0293063) q[0];
sx q[0];
rz(-1.9554546) q[0];
rz(0.063749464) q[1];
sx q[1];
rz(-1.5260162) q[1];
sx q[1];
rz(2.729111) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.257791) q[0];
sx q[0];
rz(-2.3915419) q[0];
sx q[0];
rz(-0.81410637) q[0];
rz(1.3774728) q[2];
sx q[2];
rz(-0.68916048) q[2];
sx q[2];
rz(1.2528407) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8195678) q[1];
sx q[1];
rz(-1.7734999) q[1];
sx q[1];
rz(0.2804108) q[1];
rz(-0.90610151) q[3];
sx q[3];
rz(-2.0482488) q[3];
sx q[3];
rz(-0.64732823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5731262) q[2];
sx q[2];
rz(-0.44034475) q[2];
sx q[2];
rz(1.9754515) q[2];
rz(-0.79683534) q[3];
sx q[3];
rz(-1.3692057) q[3];
sx q[3];
rz(0.78127512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(1.2260988) q[0];
sx q[0];
rz(-1.6330566) q[0];
sx q[0];
rz(1.728212) q[0];
rz(0.49890292) q[1];
sx q[1];
rz(-1.3830802) q[1];
sx q[1];
rz(1.8477207) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5099611) q[0];
sx q[0];
rz(-0.089279739) q[0];
sx q[0];
rz(1.2651612) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2433238) q[2];
sx q[2];
rz(-1.4988384) q[2];
sx q[2];
rz(1.5639992) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.41999757) q[1];
sx q[1];
rz(-1.6328057) q[1];
sx q[1];
rz(-1.5739784) q[1];
rz(0.067469723) q[3];
sx q[3];
rz(-1.6327792) q[3];
sx q[3];
rz(0.27674473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.52900195) q[2];
sx q[2];
rz(-0.68098536) q[2];
sx q[2];
rz(2.9717818) q[2];
rz(-2.7134907) q[3];
sx q[3];
rz(-1.7120818) q[3];
sx q[3];
rz(-2.9812109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.661628) q[0];
sx q[0];
rz(-0.39230883) q[0];
sx q[0];
rz(2.3062134) q[0];
rz(-2.5021878) q[1];
sx q[1];
rz(-2.4765922) q[1];
sx q[1];
rz(1.0909874) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8919286) q[0];
sx q[0];
rz(-0.87452379) q[0];
sx q[0];
rz(-1.5868756) q[0];
rz(-pi) q[1];
rz(-1.7382938) q[2];
sx q[2];
rz(-1.8665627) q[2];
sx q[2];
rz(-2.5633321) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9613232) q[1];
sx q[1];
rz(-1.9641479) q[1];
sx q[1];
rz(2.250227) q[1];
rz(2.8150878) q[3];
sx q[3];
rz(-2.0707664) q[3];
sx q[3];
rz(0.96159354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.368448) q[2];
sx q[2];
rz(-1.9545363) q[2];
sx q[2];
rz(-0.0040357987) q[2];
rz(-1.076738) q[3];
sx q[3];
rz(-0.69759798) q[3];
sx q[3];
rz(-0.18690404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8232329) q[0];
sx q[0];
rz(-1.3613181) q[0];
sx q[0];
rz(-0.63364345) q[0];
rz(2.7404495) q[1];
sx q[1];
rz(-2.432343) q[1];
sx q[1];
rz(-0.69222442) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.322583) q[0];
sx q[0];
rz(-0.87912175) q[0];
sx q[0];
rz(2.5965967) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.879331) q[2];
sx q[2];
rz(-0.82595982) q[2];
sx q[2];
rz(-2.8070297) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8068829) q[1];
sx q[1];
rz(-0.63155424) q[1];
sx q[1];
rz(1.1401661) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4603087) q[3];
sx q[3];
rz(-0.23698209) q[3];
sx q[3];
rz(-2.9221688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.20574084) q[2];
sx q[2];
rz(-0.78448272) q[2];
sx q[2];
rz(-1.240823) q[2];
rz(-1.9567018) q[3];
sx q[3];
rz(-1.3013867) q[3];
sx q[3];
rz(-1.5968116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1320837) q[0];
sx q[0];
rz(-1.8385831) q[0];
sx q[0];
rz(1.7708154) q[0];
rz(2.7235183) q[1];
sx q[1];
rz(-1.6590174) q[1];
sx q[1];
rz(2.6549227) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3330227) q[0];
sx q[0];
rz(-1.4736703) q[0];
sx q[0];
rz(2.1311231) q[0];
rz(-pi) q[1];
rz(-1.0008294) q[2];
sx q[2];
rz(-2.5116133) q[2];
sx q[2];
rz(-2.2905245) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3102973) q[1];
sx q[1];
rz(-2.8716209) q[1];
sx q[1];
rz(-0.84862535) q[1];
rz(1.5237996) q[3];
sx q[3];
rz(-2.040427) q[3];
sx q[3];
rz(2.3835973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.439094) q[2];
sx q[2];
rz(-1.0489901) q[2];
sx q[2];
rz(-2.9443963) q[2];
rz(-0.50944734) q[3];
sx q[3];
rz(-0.26064894) q[3];
sx q[3];
rz(2.313224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2800804) q[0];
sx q[0];
rz(-1.0741638) q[0];
sx q[0];
rz(1.3664838) q[0];
rz(0.81661433) q[1];
sx q[1];
rz(-2.0520703) q[1];
sx q[1];
rz(0.28269592) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3168516) q[0];
sx q[0];
rz(-2.5159989) q[0];
sx q[0];
rz(-1.8815142) q[0];
rz(-pi) q[1];
x q[1];
rz(0.31158547) q[2];
sx q[2];
rz(-1.6610314) q[2];
sx q[2];
rz(0.923522) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7273409) q[1];
sx q[1];
rz(-1.7569572) q[1];
sx q[1];
rz(0.8826137) q[1];
x q[2];
rz(1.9603086) q[3];
sx q[3];
rz(-1.0822902) q[3];
sx q[3];
rz(-2.8581885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.035159811) q[2];
sx q[2];
rz(-1.0059493) q[2];
sx q[2];
rz(2.1584568) q[2];
rz(-1.2784917) q[3];
sx q[3];
rz(-0.92517868) q[3];
sx q[3];
rz(-2.2127051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4814257) q[0];
sx q[0];
rz(-1.3496512) q[0];
sx q[0];
rz(2.8283258) q[0];
rz(-0.52472862) q[1];
sx q[1];
rz(-2.8786761) q[1];
sx q[1];
rz(1.8892586) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1850644) q[0];
sx q[0];
rz(-1.4502429) q[0];
sx q[0];
rz(0.53114364) q[0];
x q[1];
rz(2.5400794) q[2];
sx q[2];
rz(-2.2819715) q[2];
sx q[2];
rz(-1.2486697) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5460427) q[1];
sx q[1];
rz(-0.79985207) q[1];
sx q[1];
rz(-2.355666) q[1];
x q[2];
rz(-0.081324025) q[3];
sx q[3];
rz(-2.3235851) q[3];
sx q[3];
rz(2.6214056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.19920443) q[2];
sx q[2];
rz(-0.62493268) q[2];
sx q[2];
rz(1.3814629) q[2];
rz(-2.840461) q[3];
sx q[3];
rz(-0.98265177) q[3];
sx q[3];
rz(-0.38100955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23406601) q[0];
sx q[0];
rz(-2.1387687) q[0];
sx q[0];
rz(1.3475077) q[0];
rz(-2.2566336) q[1];
sx q[1];
rz(-2.5159409) q[1];
sx q[1];
rz(1.7431097) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2725582) q[0];
sx q[0];
rz(-1.5814651) q[0];
sx q[0];
rz(-1.19633) q[0];
rz(-1.3560881) q[2];
sx q[2];
rz(-0.71678679) q[2];
sx q[2];
rz(1.4625664) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.32244476) q[1];
sx q[1];
rz(-1.6080753) q[1];
sx q[1];
rz(3.0611866) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2169525) q[3];
sx q[3];
rz(-0.5914591) q[3];
sx q[3];
rz(0.22218695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0884023) q[2];
sx q[2];
rz(-1.8209499) q[2];
sx q[2];
rz(-0.94584805) q[2];
rz(2.451402) q[3];
sx q[3];
rz(-1.918957) q[3];
sx q[3];
rz(1.3010196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.239554) q[0];
sx q[0];
rz(-1.6033462) q[0];
sx q[0];
rz(-0.70815804) q[0];
rz(-0.061307727) q[1];
sx q[1];
rz(-2.5262482) q[1];
sx q[1];
rz(1.6940438) q[1];
rz(-1.6534783) q[2];
sx q[2];
rz(-1.2808242) q[2];
sx q[2];
rz(-2.6416953) q[2];
rz(-1.9504299) q[3];
sx q[3];
rz(-2.9480231) q[3];
sx q[3];
rz(2.775016) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
