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
rz(2.6235629) q[0];
sx q[0];
rz(-1.1413483) q[0];
sx q[0];
rz(-3.1014693) q[0];
rz(-2.8934381) q[1];
sx q[1];
rz(-2.0791972) q[1];
sx q[1];
rz(0.084029347) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7617159) q[0];
sx q[0];
rz(-1.5119702) q[0];
sx q[0];
rz(3.0703138) q[0];
rz(-0.40496396) q[2];
sx q[2];
rz(-2.7128007) q[2];
sx q[2];
rz(0.85067174) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5539816) q[1];
sx q[1];
rz(-2.0412672) q[1];
sx q[1];
rz(1.7616097) q[1];
rz(2.5474882) q[3];
sx q[3];
rz(-2.4770774) q[3];
sx q[3];
rz(-2.6474109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4485126) q[2];
sx q[2];
rz(-1.3658407) q[2];
sx q[2];
rz(-2.8903294) q[2];
rz(1.9692028) q[3];
sx q[3];
rz(-2.0386212) q[3];
sx q[3];
rz(0.049886543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4548816) q[0];
sx q[0];
rz(-0.48191163) q[0];
sx q[0];
rz(1.3739817) q[0];
rz(0.33085597) q[1];
sx q[1];
rz(-1.7287946) q[1];
sx q[1];
rz(-0.27423283) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2811063) q[0];
sx q[0];
rz(-1.8869068) q[0];
sx q[0];
rz(-1.64539) q[0];
rz(-pi) q[1];
rz(0.012146668) q[2];
sx q[2];
rz(-2.6082905) q[2];
sx q[2];
rz(-2.8664194) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.17473681) q[1];
sx q[1];
rz(-1.2464332) q[1];
sx q[1];
rz(1.7992538) q[1];
rz(-1.1521454) q[3];
sx q[3];
rz(-1.6466221) q[3];
sx q[3];
rz(-2.0328498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8809044) q[2];
sx q[2];
rz(-1.8131249) q[2];
sx q[2];
rz(1.0880967) q[2];
rz(1.043383) q[3];
sx q[3];
rz(-2.9731396) q[3];
sx q[3];
rz(2.8823631) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.008721) q[0];
sx q[0];
rz(-0.50478029) q[0];
sx q[0];
rz(2.0399427) q[0];
rz(-2.3654826) q[1];
sx q[1];
rz(-1.2932237) q[1];
sx q[1];
rz(-2.4993842) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6519847) q[0];
sx q[0];
rz(-2.5074158) q[0];
sx q[0];
rz(2.339667) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4218161) q[2];
sx q[2];
rz(-1.7584561) q[2];
sx q[2];
rz(-1.3087147) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2930729) q[1];
sx q[1];
rz(-0.43918931) q[1];
sx q[1];
rz(1.8166914) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.95313426) q[3];
sx q[3];
rz(-1.0079621) q[3];
sx q[3];
rz(1.9632393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8316101) q[2];
sx q[2];
rz(-1.1534561) q[2];
sx q[2];
rz(0.93008053) q[2];
rz(0.81405226) q[3];
sx q[3];
rz(-1.9508773) q[3];
sx q[3];
rz(-1.9109776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.880068) q[0];
sx q[0];
rz(-1.1898758) q[0];
sx q[0];
rz(-0.50450182) q[0];
rz(-2.2241459) q[1];
sx q[1];
rz(-0.73926273) q[1];
sx q[1];
rz(2.9333072) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4513826) q[0];
sx q[0];
rz(-1.8135895) q[0];
sx q[0];
rz(-1.7240748) q[0];
rz(-2.8253352) q[2];
sx q[2];
rz(-0.74511792) q[2];
sx q[2];
rz(-0.032501246) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5985344) q[1];
sx q[1];
rz(-2.0428791) q[1];
sx q[1];
rz(2.2170539) q[1];
x q[2];
rz(-0.55996079) q[3];
sx q[3];
rz(-1.7907639) q[3];
sx q[3];
rz(0.44693963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8248799) q[2];
sx q[2];
rz(-1.0944288) q[2];
sx q[2];
rz(-1.1701976) q[2];
rz(-2.1757388) q[3];
sx q[3];
rz(-1.0713157) q[3];
sx q[3];
rz(-2.9962208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.589094) q[0];
sx q[0];
rz(-2.4674802) q[0];
sx q[0];
rz(0.66459769) q[0];
rz(-2.873114) q[1];
sx q[1];
rz(-1.45739) q[1];
sx q[1];
rz(0.99884117) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3420606) q[0];
sx q[0];
rz(-2.149509) q[0];
sx q[0];
rz(-0.18216746) q[0];
rz(-pi) q[1];
rz(1.2810318) q[2];
sx q[2];
rz(-2.701607) q[2];
sx q[2];
rz(0.54674613) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.70390648) q[1];
sx q[1];
rz(-1.6900628) q[1];
sx q[1];
rz(1.7410623) q[1];
rz(-pi) q[2];
rz(-3.1206297) q[3];
sx q[3];
rz(-0.87910324) q[3];
sx q[3];
rz(1.5251336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9773679) q[2];
sx q[2];
rz(-2.2013142) q[2];
sx q[2];
rz(-0.11717907) q[2];
rz(0.059727877) q[3];
sx q[3];
rz(-1.2195339) q[3];
sx q[3];
rz(0.81502771) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77434671) q[0];
sx q[0];
rz(-0.4137488) q[0];
sx q[0];
rz(1.2600979) q[0];
rz(-0.11481181) q[1];
sx q[1];
rz(-1.7962619) q[1];
sx q[1];
rz(0.095349163) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0336825) q[0];
sx q[0];
rz(-2.9386407) q[0];
sx q[0];
rz(-3.0686707) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8666519) q[2];
sx q[2];
rz(-2.618578) q[2];
sx q[2];
rz(2.366334) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9748677) q[1];
sx q[1];
rz(-2.0151867) q[1];
sx q[1];
rz(-2.6728898) q[1];
rz(-pi) q[2];
x q[2];
rz(0.010706832) q[3];
sx q[3];
rz(-2.0981776) q[3];
sx q[3];
rz(-0.37400613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.054333869) q[2];
sx q[2];
rz(-2.2237873) q[2];
sx q[2];
rz(0.66162649) q[2];
rz(2.3719487) q[3];
sx q[3];
rz(-1.4504284) q[3];
sx q[3];
rz(2.2746287) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0041644) q[0];
sx q[0];
rz(-2.7935109) q[0];
sx q[0];
rz(-2.3837756) q[0];
rz(-3.1163395) q[1];
sx q[1];
rz(-2.7460637) q[1];
sx q[1];
rz(1.9662366) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7879774) q[0];
sx q[0];
rz(-0.60871658) q[0];
sx q[0];
rz(0.65257163) q[0];
x q[1];
rz(-2.6525431) q[2];
sx q[2];
rz(-0.95121517) q[2];
sx q[2];
rz(-1.3330158) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.84056329) q[1];
sx q[1];
rz(-1.9604654) q[1];
sx q[1];
rz(-0.79857608) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0915499) q[3];
sx q[3];
rz(-2.712125) q[3];
sx q[3];
rz(2.1901401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.56388277) q[2];
sx q[2];
rz(-1.6174199) q[2];
sx q[2];
rz(-1.9263402) q[2];
rz(2.3828714) q[3];
sx q[3];
rz(-1.4164378) q[3];
sx q[3];
rz(-1.5846213) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6587081) q[0];
sx q[0];
rz(-1.2889129) q[0];
sx q[0];
rz(-2.7523852) q[0];
rz(0.088118531) q[1];
sx q[1];
rz(-1.7033109) q[1];
sx q[1];
rz(1.6099991) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9154598) q[0];
sx q[0];
rz(-2.8185039) q[0];
sx q[0];
rz(-2.3996572) q[0];
rz(-2.4803703) q[2];
sx q[2];
rz(-1.2031789) q[2];
sx q[2];
rz(1.4914545) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1080146) q[1];
sx q[1];
rz(-2.696175) q[1];
sx q[1];
rz(3.0207795) q[1];
rz(-pi) q[2];
x q[2];
rz(0.86054389) q[3];
sx q[3];
rz(-1.8657953) q[3];
sx q[3];
rz(-0.33911946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7342928) q[2];
sx q[2];
rz(-2.5278842) q[2];
sx q[2];
rz(1.3235922) q[2];
rz(-1.343441) q[3];
sx q[3];
rz(-0.7898134) q[3];
sx q[3];
rz(2.8048803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
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
rz(1.6650498) q[0];
sx q[0];
rz(-2.2344868) q[0];
sx q[0];
rz(2.6764349) q[0];
rz(3.0910885) q[1];
sx q[1];
rz(-2.2513794) q[1];
sx q[1];
rz(-0.49096289) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4024593) q[0];
sx q[0];
rz(-1.8998892) q[0];
sx q[0];
rz(-1.2335445) q[0];
rz(1.4168763) q[2];
sx q[2];
rz(-1.2668005) q[2];
sx q[2];
rz(2.9829669) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.16737398) q[1];
sx q[1];
rz(-2.189296) q[1];
sx q[1];
rz(2.2688686) q[1];
x q[2];
rz(1.6816986) q[3];
sx q[3];
rz(-1.2909799) q[3];
sx q[3];
rz(-2.1808437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4831873) q[2];
sx q[2];
rz(-2.5453973) q[2];
sx q[2];
rz(-0.8636221) q[2];
rz(-0.18381271) q[3];
sx q[3];
rz(-2.176216) q[3];
sx q[3];
rz(-0.98627728) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7154536) q[0];
sx q[0];
rz(-1.621839) q[0];
sx q[0];
rz(-2.5157628) q[0];
rz(-0.65678701) q[1];
sx q[1];
rz(-2.1757809) q[1];
sx q[1];
rz(-1.2084557) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10851181) q[0];
sx q[0];
rz(-1.7038561) q[0];
sx q[0];
rz(-3.0390374) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7371128) q[2];
sx q[2];
rz(-2.4411429) q[2];
sx q[2];
rz(3.1062448) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0565666) q[1];
sx q[1];
rz(-1.657227) q[1];
sx q[1];
rz(-1.5734476) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3937065) q[3];
sx q[3];
rz(-1.0064126) q[3];
sx q[3];
rz(1.0799112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5829696) q[2];
sx q[2];
rz(-1.308459) q[2];
sx q[2];
rz(-0.51897955) q[2];
rz(-0.44273043) q[3];
sx q[3];
rz(-1.818592) q[3];
sx q[3];
rz(-1.3435266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9291572) q[0];
sx q[0];
rz(-0.99016187) q[0];
sx q[0];
rz(-1.1708175) q[0];
rz(0.15405542) q[1];
sx q[1];
rz(-2.2046721) q[1];
sx q[1];
rz(2.2278723) q[1];
rz(3.0689756) q[2];
sx q[2];
rz(-0.5177064) q[2];
sx q[2];
rz(0.032042423) q[2];
rz(1.1476573) q[3];
sx q[3];
rz(-0.71865766) q[3];
sx q[3];
rz(-1.3620993) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
