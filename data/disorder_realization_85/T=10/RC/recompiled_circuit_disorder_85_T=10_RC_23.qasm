OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.934259) q[0];
sx q[0];
rz(-0.59036314) q[0];
sx q[0];
rz(-2.7705749) q[0];
rz(2.7603005) q[1];
sx q[1];
rz(-2.5420904) q[1];
sx q[1];
rz(-1.376027) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8126412) q[0];
sx q[0];
rz(-2.0094123) q[0];
sx q[0];
rz(-2.2493275) q[0];
rz(-pi) q[1];
rz(2.5901428) q[2];
sx q[2];
rz(-2.3460238) q[2];
sx q[2];
rz(1.431682) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8810597) q[1];
sx q[1];
rz(-2.7785289) q[1];
sx q[1];
rz(2.0102324) q[1];
rz(-pi) q[2];
rz(-0.052470603) q[3];
sx q[3];
rz(-2.3215508) q[3];
sx q[3];
rz(2.6657871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.084289) q[2];
sx q[2];
rz(-2.7412582) q[2];
sx q[2];
rz(-0.98891813) q[2];
rz(-0.75254285) q[3];
sx q[3];
rz(-1.9957333) q[3];
sx q[3];
rz(2.3108216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9343524) q[0];
sx q[0];
rz(-0.11226421) q[0];
sx q[0];
rz(-1.1799312) q[0];
rz(0.99769366) q[1];
sx q[1];
rz(-1.2832063) q[1];
sx q[1];
rz(-2.4172799) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2253101) q[0];
sx q[0];
rz(-2.3819469) q[0];
sx q[0];
rz(-2.0347974) q[0];
rz(-pi) q[1];
rz(2.7472277) q[2];
sx q[2];
rz(-2.9260203) q[2];
sx q[2];
rz(2.8087316) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.69246768) q[1];
sx q[1];
rz(-1.3652703) q[1];
sx q[1];
rz(-1.2114026) q[1];
rz(-2.904326) q[3];
sx q[3];
rz(-1.7627343) q[3];
sx q[3];
rz(-0.81258472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2362242) q[2];
sx q[2];
rz(-1.8111818) q[2];
sx q[2];
rz(-0.36188564) q[2];
rz(-0.13606717) q[3];
sx q[3];
rz(-0.55570221) q[3];
sx q[3];
rz(-3.0959685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1419462) q[0];
sx q[0];
rz(-1.0936341) q[0];
sx q[0];
rz(-1.746159) q[0];
rz(-2.6793001) q[1];
sx q[1];
rz(-2.7170083) q[1];
sx q[1];
rz(-1.9225072) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8273979) q[0];
sx q[0];
rz(-2.6959246) q[0];
sx q[0];
rz(-0.92339869) q[0];
rz(-pi) q[1];
rz(2.2723324) q[2];
sx q[2];
rz(-1.3420891) q[2];
sx q[2];
rz(1.2094091) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1278652) q[1];
sx q[1];
rz(-1.3294819) q[1];
sx q[1];
rz(-1.024854) q[1];
x q[2];
rz(1.5145281) q[3];
sx q[3];
rz(-2.8731822) q[3];
sx q[3];
rz(1.0750107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.42157713) q[2];
sx q[2];
rz(-2.4013459) q[2];
sx q[2];
rz(1.6960309) q[2];
rz(-2.5727663) q[3];
sx q[3];
rz(-2.2918662) q[3];
sx q[3];
rz(0.11051699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22359426) q[0];
sx q[0];
rz(-0.44819865) q[0];
sx q[0];
rz(-0.51825994) q[0];
rz(-0.7154243) q[1];
sx q[1];
rz(-2.0253069) q[1];
sx q[1];
rz(2.3148361) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2402156) q[0];
sx q[0];
rz(-0.6745406) q[0];
sx q[0];
rz(0.69112372) q[0];
x q[1];
rz(-0.13917285) q[2];
sx q[2];
rz(-0.77251245) q[2];
sx q[2];
rz(-3.006209) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0485059) q[1];
sx q[1];
rz(-2.6566681) q[1];
sx q[1];
rz(-2.4946458) q[1];
rz(-2.4171962) q[3];
sx q[3];
rz(-1.8076234) q[3];
sx q[3];
rz(-0.5326007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9892019) q[2];
sx q[2];
rz(-2.9426136) q[2];
sx q[2];
rz(-1.7626804) q[2];
rz(-0.072323024) q[3];
sx q[3];
rz(-0.81243378) q[3];
sx q[3];
rz(-1.4962083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3146661) q[0];
sx q[0];
rz(-2.5214654) q[0];
sx q[0];
rz(2.0157053) q[0];
rz(0.90244883) q[1];
sx q[1];
rz(-2.1676962) q[1];
sx q[1];
rz(0.28516969) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0641159) q[0];
sx q[0];
rz(-0.47668326) q[0];
sx q[0];
rz(-1.478273) q[0];
x q[1];
rz(-0.6638078) q[2];
sx q[2];
rz(-1.8823349) q[2];
sx q[2];
rz(0.89154348) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.73357108) q[1];
sx q[1];
rz(-1.5703778) q[1];
sx q[1];
rz(-1.8838521) q[1];
x q[2];
rz(-2.9167446) q[3];
sx q[3];
rz(-0.41089155) q[3];
sx q[3];
rz(-1.2332682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.29331648) q[2];
sx q[2];
rz(-0.50555503) q[2];
sx q[2];
rz(0.53945333) q[2];
rz(-2.8347677) q[3];
sx q[3];
rz(-0.8845194) q[3];
sx q[3];
rz(-2.6873798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8834615) q[0];
sx q[0];
rz(-2.3911609) q[0];
sx q[0];
rz(2.9845797) q[0];
rz(2.4482588) q[1];
sx q[1];
rz(-2.2608829) q[1];
sx q[1];
rz(-1.7745811) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.052664) q[0];
sx q[0];
rz(-0.0757218) q[0];
sx q[0];
rz(-2.0220387) q[0];
rz(-0.579367) q[2];
sx q[2];
rz(-0.89951347) q[2];
sx q[2];
rz(2.6210149) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.64703343) q[1];
sx q[1];
rz(-1.651598) q[1];
sx q[1];
rz(-1.2809491) q[1];
rz(1.6744162) q[3];
sx q[3];
rz(-1.2195671) q[3];
sx q[3];
rz(-2.2675089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8453025) q[2];
sx q[2];
rz(-0.85001105) q[2];
sx q[2];
rz(0.40346754) q[2];
rz(-2.6599595) q[3];
sx q[3];
rz(-1.0721595) q[3];
sx q[3];
rz(2.6223555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24213174) q[0];
sx q[0];
rz(-2.2583028) q[0];
sx q[0];
rz(-2.2677299) q[0];
rz(-0.44772398) q[1];
sx q[1];
rz(-2.402585) q[1];
sx q[1];
rz(1.9708995) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0872333) q[0];
sx q[0];
rz(-1.5938252) q[0];
sx q[0];
rz(1.5769632) q[0];
x q[1];
rz(3.0253719) q[2];
sx q[2];
rz(-1.6214317) q[2];
sx q[2];
rz(-1.2223787) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.891174) q[1];
sx q[1];
rz(-2.182057) q[1];
sx q[1];
rz(-0.82959081) q[1];
rz(-pi) q[2];
rz(-2.6074334) q[3];
sx q[3];
rz(-1.8338406) q[3];
sx q[3];
rz(2.0260889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0968904) q[2];
sx q[2];
rz(-2.5723852) q[2];
sx q[2];
rz(2.5308385) q[2];
rz(2.6664873) q[3];
sx q[3];
rz(-2.0510309) q[3];
sx q[3];
rz(-0.92774123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8996745) q[0];
sx q[0];
rz(-3.0245259) q[0];
sx q[0];
rz(2.8444667) q[0];
rz(-1.7469453) q[1];
sx q[1];
rz(-1.9906094) q[1];
sx q[1];
rz(-0.64613211) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0028249) q[0];
sx q[0];
rz(-1.3647623) q[0];
sx q[0];
rz(-0.10319184) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6693194) q[2];
sx q[2];
rz(-1.7211282) q[2];
sx q[2];
rz(-2.3538102) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.38478002) q[1];
sx q[1];
rz(-1.9128748) q[1];
sx q[1];
rz(-2.6095819) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5313247) q[3];
sx q[3];
rz(-2.2329997) q[3];
sx q[3];
rz(2.3074647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.064676553) q[2];
sx q[2];
rz(-2.1885394) q[2];
sx q[2];
rz(0.45483744) q[2];
rz(-2.440195) q[3];
sx q[3];
rz(-2.111179) q[3];
sx q[3];
rz(2.0075683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4421473) q[0];
sx q[0];
rz(-3*pi/16) q[0];
sx q[0];
rz(-2.3440857) q[0];
rz(0.51756716) q[1];
sx q[1];
rz(-2.3289754) q[1];
sx q[1];
rz(0.10841766) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90005504) q[0];
sx q[0];
rz(-2.2336322) q[0];
sx q[0];
rz(-0.073297757) q[0];
rz(-1.8234532) q[2];
sx q[2];
rz(-0.18441072) q[2];
sx q[2];
rz(-0.84469634) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.79418102) q[1];
sx q[1];
rz(-1.3550183) q[1];
sx q[1];
rz(2.9706035) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.18687825) q[3];
sx q[3];
rz(-2.2664321) q[3];
sx q[3];
rz(-0.38284341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0124399) q[2];
sx q[2];
rz(-1.7947349) q[2];
sx q[2];
rz(-0.33995315) q[2];
rz(-0.41839504) q[3];
sx q[3];
rz(-2.5451626) q[3];
sx q[3];
rz(2.4160014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6091992) q[0];
sx q[0];
rz(-2.7476855) q[0];
sx q[0];
rz(-2.4627731) q[0];
rz(0.36418307) q[1];
sx q[1];
rz(-1.4441676) q[1];
sx q[1];
rz(-3.0864339) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3653152) q[0];
sx q[0];
rz(-2.758983) q[0];
sx q[0];
rz(-0.17392735) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7807547) q[2];
sx q[2];
rz(-1.0210438) q[2];
sx q[2];
rz(-1.2325665) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.058051) q[1];
sx q[1];
rz(-1.0293048) q[1];
sx q[1];
rz(0.70152775) q[1];
rz(-1.6454562) q[3];
sx q[3];
rz(-0.45711043) q[3];
sx q[3];
rz(0.32170579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1577592) q[2];
sx q[2];
rz(-2.1201717) q[2];
sx q[2];
rz(2.6514163) q[2];
rz(0.13752078) q[3];
sx q[3];
rz(-1.0995882) q[3];
sx q[3];
rz(-0.93808758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4162083) q[0];
sx q[0];
rz(-1.2513456) q[0];
sx q[0];
rz(-0.67847897) q[0];
rz(-0.2086808) q[1];
sx q[1];
rz(-1.122767) q[1];
sx q[1];
rz(-1.541419) q[1];
rz(1.4738884) q[2];
sx q[2];
rz(-1.2599535) q[2];
sx q[2];
rz(-2.8355666) q[2];
rz(-1.6155852) q[3];
sx q[3];
rz(-1.3309892) q[3];
sx q[3];
rz(0.25467024) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
