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
rz(-0.51802975) q[0];
sx q[0];
rz(-2.0002444) q[0];
sx q[0];
rz(3.1014693) q[0];
rz(0.24815458) q[1];
sx q[1];
rz(-1.0623955) q[1];
sx q[1];
rz(-0.084029347) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7617159) q[0];
sx q[0];
rz(-1.5119702) q[0];
sx q[0];
rz(-0.071278871) q[0];
rz(-pi) q[1];
rz(-1.3925926) q[2];
sx q[2];
rz(-1.1786818) q[2];
sx q[2];
rz(2.7313679) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2457468) q[1];
sx q[1];
rz(-1.7406643) q[1];
sx q[1];
rz(2.6636811) q[1];
rz(-pi) q[2];
rz(-2.5474882) q[3];
sx q[3];
rz(-2.4770774) q[3];
sx q[3];
rz(-0.49418172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6930801) q[2];
sx q[2];
rz(-1.7757519) q[2];
sx q[2];
rz(2.8903294) q[2];
rz(1.1723899) q[3];
sx q[3];
rz(-2.0386212) q[3];
sx q[3];
rz(-0.049886543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4548816) q[0];
sx q[0];
rz(-2.659681) q[0];
sx q[0];
rz(-1.7676109) q[0];
rz(-0.33085597) q[1];
sx q[1];
rz(-1.7287946) q[1];
sx q[1];
rz(-2.8673598) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2811063) q[0];
sx q[0];
rz(-1.2546859) q[0];
sx q[0];
rz(1.4962026) q[0];
x q[1];
rz(-3.129446) q[2];
sx q[2];
rz(-2.6082905) q[2];
sx q[2];
rz(0.27517327) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.17473681) q[1];
sx q[1];
rz(-1.8951594) q[1];
sx q[1];
rz(1.3423389) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0586309) q[3];
sx q[3];
rz(-1.1534235) q[3];
sx q[3];
rz(2.7132332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8809044) q[2];
sx q[2];
rz(-1.8131249) q[2];
sx q[2];
rz(-1.0880967) q[2];
rz(1.043383) q[3];
sx q[3];
rz(-0.16845307) q[3];
sx q[3];
rz(-2.8823631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-2.1328717) q[0];
sx q[0];
rz(-0.50478029) q[0];
sx q[0];
rz(1.1016499) q[0];
rz(2.3654826) q[1];
sx q[1];
rz(-1.2932237) q[1];
sx q[1];
rz(-0.64220846) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61313907) q[0];
sx q[0];
rz(-2.0106843) q[0];
sx q[0];
rz(-0.47275193) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8611254) q[2];
sx q[2];
rz(-0.73958042) q[2];
sx q[2];
rz(2.6698339) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.57802618) q[1];
sx q[1];
rz(-1.1456923) q[1];
sx q[1];
rz(3.0277287) q[1];
rz(-pi) q[2];
x q[2];
rz(0.65861838) q[3];
sx q[3];
rz(-2.082654) q[3];
sx q[3];
rz(0.030125253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3099826) q[2];
sx q[2];
rz(-1.9881366) q[2];
sx q[2];
rz(0.93008053) q[2];
rz(-2.3275404) q[3];
sx q[3];
rz(-1.1907153) q[3];
sx q[3];
rz(1.9109776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.880068) q[0];
sx q[0];
rz(-1.1898758) q[0];
sx q[0];
rz(-2.6370908) q[0];
rz(-2.2241459) q[1];
sx q[1];
rz(-0.73926273) q[1];
sx q[1];
rz(2.9333072) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.298131) q[0];
sx q[0];
rz(-1.7195452) q[0];
sx q[0];
rz(2.896033) q[0];
rz(-1.850205) q[2];
sx q[2];
rz(-2.2710751) q[2];
sx q[2];
rz(0.4513739) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.54305824) q[1];
sx q[1];
rz(-2.0428791) q[1];
sx q[1];
rz(-2.2170539) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8288021) q[3];
sx q[3];
rz(-2.1157221) q[3];
sx q[3];
rz(1.2598002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.31671277) q[2];
sx q[2];
rz(-2.0471639) q[2];
sx q[2];
rz(1.1701976) q[2];
rz(2.1757388) q[3];
sx q[3];
rz(-1.0713157) q[3];
sx q[3];
rz(-0.14537183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55249864) q[0];
sx q[0];
rz(-0.67411244) q[0];
sx q[0];
rz(-0.66459769) q[0];
rz(-0.26847863) q[1];
sx q[1];
rz(-1.6842027) q[1];
sx q[1];
rz(-2.1427515) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12832175) q[0];
sx q[0];
rz(-1.7230464) q[0];
sx q[0];
rz(0.98441846) q[0];
rz(-0.13370698) q[2];
sx q[2];
rz(-1.1503476) q[2];
sx q[2];
rz(0.22841098) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2951569) q[1];
sx q[1];
rz(-1.4017516) q[1];
sx q[1];
rz(-3.0205932) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.596101) q[3];
sx q[3];
rz(-0.69195834) q[3];
sx q[3];
rz(1.4922752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9773679) q[2];
sx q[2];
rz(-2.2013142) q[2];
sx q[2];
rz(0.11717907) q[2];
rz(-3.0818648) q[3];
sx q[3];
rz(-1.9220587) q[3];
sx q[3];
rz(-0.81502771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77434671) q[0];
sx q[0];
rz(-0.4137488) q[0];
sx q[0];
rz(1.2600979) q[0];
rz(0.11481181) q[1];
sx q[1];
rz(-1.3453307) q[1];
sx q[1];
rz(0.095349163) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1079102) q[0];
sx q[0];
rz(-0.20295197) q[0];
sx q[0];
rz(3.0686707) q[0];
rz(-1.0667858) q[2];
sx q[2];
rz(-1.4246449) q[2];
sx q[2];
rz(0.53737568) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.29987511) q[1];
sx q[1];
rz(-0.63427329) q[1];
sx q[1];
rz(0.81171616) q[1];
rz(-pi) q[2];
x q[2];
rz(0.010706832) q[3];
sx q[3];
rz(-2.0981776) q[3];
sx q[3];
rz(2.7675865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.054333869) q[2];
sx q[2];
rz(-0.91780535) q[2];
sx q[2];
rz(2.4799662) q[2];
rz(2.3719487) q[3];
sx q[3];
rz(-1.4504284) q[3];
sx q[3];
rz(-0.86696398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0041644) q[0];
sx q[0];
rz(-0.34808174) q[0];
sx q[0];
rz(0.75781703) q[0];
rz(-3.1163395) q[1];
sx q[1];
rz(-2.7460637) q[1];
sx q[1];
rz(-1.1753561) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34285082) q[0];
sx q[0];
rz(-1.2161868) q[0];
sx q[0];
rz(0.50574982) q[0];
rz(-2.1532159) q[2];
sx q[2];
rz(-2.3727131) q[2];
sx q[2];
rz(2.0748367) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3010294) q[1];
sx q[1];
rz(-1.9604654) q[1];
sx q[1];
rz(2.3430166) q[1];
rz(-pi) q[2];
rz(2.712599) q[3];
sx q[3];
rz(-1.5499664) q[3];
sx q[3];
rz(-0.57383895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.56388277) q[2];
sx q[2];
rz(-1.5241728) q[2];
sx q[2];
rz(1.2152524) q[2];
rz(0.75872129) q[3];
sx q[3];
rz(-1.4164378) q[3];
sx q[3];
rz(1.5846213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6587081) q[0];
sx q[0];
rz(-1.8526798) q[0];
sx q[0];
rz(-2.7523852) q[0];
rz(0.088118531) q[1];
sx q[1];
rz(-1.4382818) q[1];
sx q[1];
rz(1.5315936) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.599278) q[0];
sx q[0];
rz(-1.3345584) q[0];
sx q[0];
rz(1.3482984) q[0];
rz(-pi) q[1];
rz(2.0247718) q[2];
sx q[2];
rz(-2.1810227) q[2];
sx q[2];
rz(2.7895841) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.89985049) q[1];
sx q[1];
rz(-2.012737) q[1];
sx q[1];
rz(1.5133218) q[1];
rz(-pi) q[2];
rz(1.1346899) q[3];
sx q[3];
rz(-0.75911555) q[3];
sx q[3];
rz(-1.5576943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4072998) q[2];
sx q[2];
rz(-2.5278842) q[2];
sx q[2];
rz(-1.8180004) q[2];
rz(1.343441) q[3];
sx q[3];
rz(-0.7898134) q[3];
sx q[3];
rz(0.33671236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6650498) q[0];
sx q[0];
rz(-0.90710586) q[0];
sx q[0];
rz(-2.6764349) q[0];
rz(3.0910885) q[1];
sx q[1];
rz(-2.2513794) q[1];
sx q[1];
rz(2.6506298) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28117958) q[0];
sx q[0];
rz(-1.2523012) q[0];
sx q[0];
rz(2.7943576) q[0];
rz(-pi) q[1];
rz(-2.8341836) q[2];
sx q[2];
rz(-1.4239862) q[2];
sx q[2];
rz(1.4585782) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.16737398) q[1];
sx q[1];
rz(-2.189296) q[1];
sx q[1];
rz(2.2688686) q[1];
rz(2.773953) q[3];
sx q[3];
rz(-2.841137) q[3];
sx q[3];
rz(2.5641172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4831873) q[2];
sx q[2];
rz(-0.5961954) q[2];
sx q[2];
rz(-2.2779706) q[2];
rz(0.18381271) q[3];
sx q[3];
rz(-0.96537662) q[3];
sx q[3];
rz(2.1553154) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7154536) q[0];
sx q[0];
rz(-1.621839) q[0];
sx q[0];
rz(2.5157628) q[0];
rz(0.65678701) q[1];
sx q[1];
rz(-0.96581179) q[1];
sx q[1];
rz(1.9331369) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4759373) q[0];
sx q[0];
rz(-1.6724419) q[0];
sx q[0];
rz(-1.4370421) q[0];
rz(-pi) q[1];
rz(0.65932806) q[2];
sx q[2];
rz(-1.314333) q[2];
sx q[2];
rz(-1.9224482) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0872709) q[1];
sx q[1];
rz(-0.086471237) q[1];
sx q[1];
rz(0.030589624) q[1];
rz(-0.27140433) q[3];
sx q[3];
rz(-2.5529666) q[3];
sx q[3];
rz(-1.7388104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5829696) q[2];
sx q[2];
rz(-1.8331336) q[2];
sx q[2];
rz(-0.51897955) q[2];
rz(-0.44273043) q[3];
sx q[3];
rz(-1.3230007) q[3];
sx q[3];
rz(-1.7980661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2124355) q[0];
sx q[0];
rz(-2.1514308) q[0];
sx q[0];
rz(1.9707752) q[0];
rz(-0.15405542) q[1];
sx q[1];
rz(-0.93692056) q[1];
sx q[1];
rz(-0.9137203) q[1];
rz(0.072617037) q[2];
sx q[2];
rz(-2.6238863) q[2];
sx q[2];
rz(-3.1095502) q[2];
rz(-0.34482205) q[3];
sx q[3];
rz(-2.214684) q[3];
sx q[3];
rz(-0.82292258) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
