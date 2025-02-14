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
rz(-0.80225575) q[0];
sx q[0];
rz(-1.7576317) q[0];
sx q[0];
rz(1.2686165) q[0];
rz(0.4624548) q[1];
sx q[1];
rz(-3.0825244) q[1];
sx q[1];
rz(-1.7360092) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2911513) q[0];
sx q[0];
rz(-1.1014525) q[0];
sx q[0];
rz(2.6525081) q[0];
x q[1];
rz(1.8008158) q[2];
sx q[2];
rz(-2.2284796) q[2];
sx q[2];
rz(2.7000526) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2785223) q[1];
sx q[1];
rz(-2.2691326) q[1];
sx q[1];
rz(2.5935943) q[1];
rz(2.2470993) q[3];
sx q[3];
rz(-1.749012) q[3];
sx q[3];
rz(1.3823929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9872226) q[2];
sx q[2];
rz(-1.9742249) q[2];
sx q[2];
rz(2.3517189) q[2];
rz(-3.1176944) q[3];
sx q[3];
rz(-1.3393211) q[3];
sx q[3];
rz(0.53643119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8928878) q[0];
sx q[0];
rz(-1.6510115) q[0];
sx q[0];
rz(1.8875246) q[0];
rz(1.6528543) q[1];
sx q[1];
rz(-2.2658927) q[1];
sx q[1];
rz(-2.658433) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89630187) q[0];
sx q[0];
rz(-1.5190796) q[0];
sx q[0];
rz(-1.381676) q[0];
x q[1];
rz(-2.0873726) q[2];
sx q[2];
rz(-0.90460515) q[2];
sx q[2];
rz(1.7214799) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.70139685) q[1];
sx q[1];
rz(-0.45397511) q[1];
sx q[1];
rz(2.5099436) q[1];
rz(-2.8511413) q[3];
sx q[3];
rz(-1.402352) q[3];
sx q[3];
rz(-0.12663933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2529605) q[2];
sx q[2];
rz(-1.4560207) q[2];
sx q[2];
rz(1.5144833) q[2];
rz(-3.0857981) q[3];
sx q[3];
rz(-1.5972219) q[3];
sx q[3];
rz(1.6519215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2273939) q[0];
sx q[0];
rz(-0.44454235) q[0];
sx q[0];
rz(-3.0134873) q[0];
rz(-2.5654492) q[1];
sx q[1];
rz(-2.4100401) q[1];
sx q[1];
rz(0.76549706) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5680117) q[0];
sx q[0];
rz(-0.67493248) q[0];
sx q[0];
rz(-0.57008596) q[0];
x q[1];
rz(0.68674318) q[2];
sx q[2];
rz(-1.8675065) q[2];
sx q[2];
rz(-2.3627757) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.24352267) q[1];
sx q[1];
rz(-1.0869893) q[1];
sx q[1];
rz(2.0953728) q[1];
rz(-0.14344826) q[3];
sx q[3];
rz(-1.5960403) q[3];
sx q[3];
rz(2.6526566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8219882) q[2];
sx q[2];
rz(-0.87279785) q[2];
sx q[2];
rz(-0.86307159) q[2];
rz(0.34267628) q[3];
sx q[3];
rz(-0.42049146) q[3];
sx q[3];
rz(-1.4884523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2494025) q[0];
sx q[0];
rz(-0.96977314) q[0];
sx q[0];
rz(0.48200193) q[0];
rz(2.8328698) q[1];
sx q[1];
rz(-2.8161507) q[1];
sx q[1];
rz(-0.6122922) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4541439) q[0];
sx q[0];
rz(-2.4061235) q[0];
sx q[0];
rz(-0.81540458) q[0];
x q[1];
rz(-0.014043645) q[2];
sx q[2];
rz(-1.6545452) q[2];
sx q[2];
rz(-1.75911) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4548851) q[1];
sx q[1];
rz(-2.5338182) q[1];
sx q[1];
rz(-1.6493504) q[1];
rz(1.5641698) q[3];
sx q[3];
rz(-1.7200791) q[3];
sx q[3];
rz(1.3914668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4838532) q[2];
sx q[2];
rz(-0.95119363) q[2];
sx q[2];
rz(-1.8335906) q[2];
rz(-2.7075503) q[3];
sx q[3];
rz(-1.3515892) q[3];
sx q[3];
rz(1.8724117) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16172116) q[0];
sx q[0];
rz(-1.5525818) q[0];
sx q[0];
rz(2.8671434) q[0];
rz(1.6203923) q[1];
sx q[1];
rz(-2.0976286) q[1];
sx q[1];
rz(1.8843947) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50097695) q[0];
sx q[0];
rz(-1.3269935) q[0];
sx q[0];
rz(-1.8200726) q[0];
rz(-pi) q[1];
rz(2.1047701) q[2];
sx q[2];
rz(-1.8920915) q[2];
sx q[2];
rz(0.23201135) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.91433954) q[1];
sx q[1];
rz(-1.7861331) q[1];
sx q[1];
rz(1.017175) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2755695) q[3];
sx q[3];
rz(-1.8187722) q[3];
sx q[3];
rz(0.06125227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7460798) q[2];
sx q[2];
rz(-1.505625) q[2];
sx q[2];
rz(2.5858322) q[2];
rz(2.3505576) q[3];
sx q[3];
rz(-1.0020703) q[3];
sx q[3];
rz(-1.5033495) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.898734) q[0];
sx q[0];
rz(-1.4396311) q[0];
sx q[0];
rz(0.71204251) q[0];
rz(-0.60246077) q[1];
sx q[1];
rz(-2.7280877) q[1];
sx q[1];
rz(-0.079040225) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3369737) q[0];
sx q[0];
rz(-0.078402407) q[0];
sx q[0];
rz(0.96167643) q[0];
rz(-0.075435813) q[2];
sx q[2];
rz(-1.4860538) q[2];
sx q[2];
rz(0.85979474) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.329603) q[1];
sx q[1];
rz(-1.7895849) q[1];
sx q[1];
rz(2.3892774) q[1];
x q[2];
rz(0.2980663) q[3];
sx q[3];
rz(-1.4722927) q[3];
sx q[3];
rz(0.57386604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1003803) q[2];
sx q[2];
rz(-2.2104287) q[2];
sx q[2];
rz(0.97990123) q[2];
rz(1.6117217) q[3];
sx q[3];
rz(-1.9001222) q[3];
sx q[3];
rz(-2.8821168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3945382) q[0];
sx q[0];
rz(-1.8600445) q[0];
sx q[0];
rz(-0.10093149) q[0];
rz(2.586567) q[1];
sx q[1];
rz(-0.9340159) q[1];
sx q[1];
rz(-1.2947882) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7332563) q[0];
sx q[0];
rz(-2.2332158) q[0];
sx q[0];
rz(2.8444879) q[0];
rz(-pi) q[1];
rz(-0.00097043911) q[2];
sx q[2];
rz(-1.8651267) q[2];
sx q[2];
rz(-2.7270339) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.76483549) q[1];
sx q[1];
rz(-1.7143814) q[1];
sx q[1];
rz(0.43068703) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0247308) q[3];
sx q[3];
rz(-2.6203558) q[3];
sx q[3];
rz(1.2574399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1780221) q[2];
sx q[2];
rz(-0.87248412) q[2];
sx q[2];
rz(-0.24173173) q[2];
rz(-2.669615) q[3];
sx q[3];
rz(-0.21502544) q[3];
sx q[3];
rz(1.5836466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3374775) q[0];
sx q[0];
rz(-0.98144704) q[0];
sx q[0];
rz(0.61019439) q[0];
rz(-0.21774165) q[1];
sx q[1];
rz(-2.1861031) q[1];
sx q[1];
rz(2.311923) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0141896) q[0];
sx q[0];
rz(-0.9547736) q[0];
sx q[0];
rz(2.4418529) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.83135624) q[2];
sx q[2];
rz(-2.6717253) q[2];
sx q[2];
rz(-2.510315) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1983932) q[1];
sx q[1];
rz(-1.7560609) q[1];
sx q[1];
rz(0.41364248) q[1];
rz(1.7046961) q[3];
sx q[3];
rz(-2.4256676) q[3];
sx q[3];
rz(-1.1435103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3024451) q[2];
sx q[2];
rz(-1.7586917) q[2];
sx q[2];
rz(-1.1350606) q[2];
rz(-0.48842397) q[3];
sx q[3];
rz(-1.4819744) q[3];
sx q[3];
rz(1.2357014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9957073) q[0];
sx q[0];
rz(-2.0304401) q[0];
sx q[0];
rz(-2.1671894) q[0];
rz(-2.3459332) q[1];
sx q[1];
rz(-2.7641422) q[1];
sx q[1];
rz(-2.5239351) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5461676) q[0];
sx q[0];
rz(-1.4873532) q[0];
sx q[0];
rz(-1.5682778) q[0];
rz(1.9699924) q[2];
sx q[2];
rz(-2.0621057) q[2];
sx q[2];
rz(-1.0688865) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.61231919) q[1];
sx q[1];
rz(-0.63748432) q[1];
sx q[1];
rz(2.6785664) q[1];
x q[2];
rz(2.8116954) q[3];
sx q[3];
rz(-2.0198235) q[3];
sx q[3];
rz(2.2638418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1254897) q[2];
sx q[2];
rz(-1.2697271) q[2];
sx q[2];
rz(1.7529091) q[2];
rz(-1.3056508) q[3];
sx q[3];
rz(-2.4180222) q[3];
sx q[3];
rz(1.8587662) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0131123) q[0];
sx q[0];
rz(-2.4025669) q[0];
sx q[0];
rz(-0.59984961) q[0];
rz(0.57303095) q[1];
sx q[1];
rz(-2.532798) q[1];
sx q[1];
rz(-1.0240239) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.065627873) q[0];
sx q[0];
rz(-0.5813404) q[0];
sx q[0];
rz(-0.79728787) q[0];
rz(3.0863831) q[2];
sx q[2];
rz(-0.6533567) q[2];
sx q[2];
rz(3.0467767) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8166067) q[1];
sx q[1];
rz(-1.4868951) q[1];
sx q[1];
rz(-0.3971667) q[1];
rz(-1.7217595) q[3];
sx q[3];
rz(-1.9486041) q[3];
sx q[3];
rz(1.0491766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2188501) q[2];
sx q[2];
rz(-1.3877733) q[2];
sx q[2];
rz(-0.77524033) q[2];
rz(-1.5357337) q[3];
sx q[3];
rz(-1.1865059) q[3];
sx q[3];
rz(1.2576013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6105462) q[0];
sx q[0];
rz(-2.1903867) q[0];
sx q[0];
rz(1.2280986) q[0];
rz(-1.3950521) q[1];
sx q[1];
rz(-1.4735305) q[1];
sx q[1];
rz(2.7458618) q[1];
rz(1.878123) q[2];
sx q[2];
rz(-2.0057445) q[2];
sx q[2];
rz(2.4705171) q[2];
rz(-1.6041605) q[3];
sx q[3];
rz(-0.76904528) q[3];
sx q[3];
rz(2.4861102) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
