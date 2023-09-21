OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.71198553) q[0];
sx q[0];
rz(3.5482121) q[0];
sx q[0];
rz(9.1756048) q[0];
rz(3.0781526) q[1];
sx q[1];
rz(-0.97172207) q[1];
sx q[1];
rz(-0.5501774) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8367856) q[0];
sx q[0];
rz(-1.3512304) q[0];
sx q[0];
rz(-3.1138793) q[0];
rz(-pi) q[1];
rz(0.43843856) q[2];
sx q[2];
rz(-2.0619832) q[2];
sx q[2];
rz(0.05471281) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5397415) q[1];
sx q[1];
rz(-1.9508233) q[1];
sx q[1];
rz(-0.7976346) q[1];
x q[2];
rz(0.95991858) q[3];
sx q[3];
rz(-1.6392518) q[3];
sx q[3];
rz(-1.5271036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7258519) q[2];
sx q[2];
rz(-0.44885138) q[2];
sx q[2];
rz(2.501781) q[2];
rz(2.2885627) q[3];
sx q[3];
rz(-0.60522389) q[3];
sx q[3];
rz(-0.38133347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8752276) q[0];
sx q[0];
rz(-1.1667644) q[0];
sx q[0];
rz(0.27045989) q[0];
rz(2.4282783) q[1];
sx q[1];
rz(-1.0353054) q[1];
sx q[1];
rz(-1.6289904) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4166116) q[0];
sx q[0];
rz(-2.1583301) q[0];
sx q[0];
rz(-2.7532817) q[0];
rz(0.28378758) q[2];
sx q[2];
rz(-1.5740526) q[2];
sx q[2];
rz(-2.6736205) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.70817972) q[1];
sx q[1];
rz(-2.4509894) q[1];
sx q[1];
rz(-1.3604926) q[1];
x q[2];
rz(1.8566425) q[3];
sx q[3];
rz(-1.3405521) q[3];
sx q[3];
rz(0.86499351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2018532) q[2];
sx q[2];
rz(-2.8994603) q[2];
sx q[2];
rz(-0.80336037) q[2];
rz(1.057829) q[3];
sx q[3];
rz(-1.6488896) q[3];
sx q[3];
rz(0.025618205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2597044) q[0];
sx q[0];
rz(-1.0972247) q[0];
sx q[0];
rz(0.29552466) q[0];
rz(-2.9064536) q[1];
sx q[1];
rz(-1.4328911) q[1];
sx q[1];
rz(-0.74584109) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64622067) q[0];
sx q[0];
rz(-1.4797858) q[0];
sx q[0];
rz(-1.792684) q[0];
rz(-pi) q[1];
rz(1.4065811) q[2];
sx q[2];
rz(-2.0914372) q[2];
sx q[2];
rz(0.28902136) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5513788) q[1];
sx q[1];
rz(-2.5775902) q[1];
sx q[1];
rz(0.53932921) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6882719) q[3];
sx q[3];
rz(-1.3461777) q[3];
sx q[3];
rz(0.9325222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.088034257) q[2];
sx q[2];
rz(-0.025667889) q[2];
sx q[2];
rz(-0.68874613) q[2];
rz(0.050343242) q[3];
sx q[3];
rz(-0.91209948) q[3];
sx q[3];
rz(1.5215727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.922309) q[0];
sx q[0];
rz(-1.0204717) q[0];
sx q[0];
rz(-3.0396089) q[0];
rz(-3.02137) q[1];
sx q[1];
rz(-0.52821237) q[1];
sx q[1];
rz(-0.27329683) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8768339) q[0];
sx q[0];
rz(-2.8059373) q[0];
sx q[0];
rz(-1.2980952) q[0];
x q[1];
rz(-2.1347413) q[2];
sx q[2];
rz(-1.0014921) q[2];
sx q[2];
rz(3.0727) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4218688) q[1];
sx q[1];
rz(-0.27742741) q[1];
sx q[1];
rz(-1.006702) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4920739) q[3];
sx q[3];
rz(-1.4276541) q[3];
sx q[3];
rz(-2.3424145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.79167241) q[2];
sx q[2];
rz(-2.4035954) q[2];
sx q[2];
rz(-2.6376574) q[2];
rz(0.079581633) q[3];
sx q[3];
rz(-1.9888398) q[3];
sx q[3];
rz(-2.7664405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85916096) q[0];
sx q[0];
rz(-1.0520042) q[0];
sx q[0];
rz(0.082745634) q[0];
rz(-2.4619608) q[1];
sx q[1];
rz(-1.4928763) q[1];
sx q[1];
rz(2.1544429) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25492451) q[0];
sx q[0];
rz(-2.3912171) q[0];
sx q[0];
rz(2.1819941) q[0];
x q[1];
rz(1.5613902) q[2];
sx q[2];
rz(-2.214589) q[2];
sx q[2];
rz(-2.8807092) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8122711) q[1];
sx q[1];
rz(-1.2580039) q[1];
sx q[1];
rz(-0.66004628) q[1];
rz(-0.042111245) q[3];
sx q[3];
rz(-1.9633506) q[3];
sx q[3];
rz(0.3596572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8905028) q[2];
sx q[2];
rz(-2.3513998) q[2];
sx q[2];
rz(2.9303072) q[2];
rz(-2.7206897) q[3];
sx q[3];
rz(-0.55287164) q[3];
sx q[3];
rz(3.0781854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6565276) q[0];
sx q[0];
rz(-1.8873029) q[0];
sx q[0];
rz(-2.3497537) q[0];
rz(0.99545288) q[1];
sx q[1];
rz(-0.96770006) q[1];
sx q[1];
rz(-1.8621559) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6364481) q[0];
sx q[0];
rz(-2.8594058) q[0];
sx q[0];
rz(1.5178773) q[0];
x q[1];
rz(1.2811786) q[2];
sx q[2];
rz(-0.98266232) q[2];
sx q[2];
rz(-2.9786125) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.93950677) q[1];
sx q[1];
rz(-2.0923951) q[1];
sx q[1];
rz(3.107198) q[1];
x q[2];
rz(2.2313314) q[3];
sx q[3];
rz(-0.22781867) q[3];
sx q[3];
rz(-1.6662625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8906158) q[2];
sx q[2];
rz(-2.1112517) q[2];
sx q[2];
rz(-2.3708169) q[2];
rz(1.4701014) q[3];
sx q[3];
rz(-2.7272868) q[3];
sx q[3];
rz(0.18338403) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32271785) q[0];
sx q[0];
rz(-0.6495496) q[0];
sx q[0];
rz(-0.055667002) q[0];
rz(-2.9260013) q[1];
sx q[1];
rz(-2.3781653) q[1];
sx q[1];
rz(0.0035704426) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63821793) q[0];
sx q[0];
rz(-2.5711381) q[0];
sx q[0];
rz(0.72871448) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7262906) q[2];
sx q[2];
rz(-0.79681444) q[2];
sx q[2];
rz(2.08564) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.53686245) q[1];
sx q[1];
rz(-2.4559048) q[1];
sx q[1];
rz(1.0889978) q[1];
rz(-pi) q[2];
rz(1.2083863) q[3];
sx q[3];
rz(-2.9299195) q[3];
sx q[3];
rz(0.025312245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.59721649) q[2];
sx q[2];
rz(-2.1901013) q[2];
sx q[2];
rz(-2.2214831) q[2];
rz(-0.19872935) q[3];
sx q[3];
rz(-1.9072429) q[3];
sx q[3];
rz(0.41771093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51628095) q[0];
sx q[0];
rz(-1.6315062) q[0];
sx q[0];
rz(-1.4177119) q[0];
rz(0.40813804) q[1];
sx q[1];
rz(-1.0393655) q[1];
sx q[1];
rz(0.67869854) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1269826) q[0];
sx q[0];
rz(-1.8982366) q[0];
sx q[0];
rz(-1.8493269) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0661725) q[2];
sx q[2];
rz(-2.3373211) q[2];
sx q[2];
rz(1.417516) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.24819599) q[1];
sx q[1];
rz(-2.0940603) q[1];
sx q[1];
rz(0.072695331) q[1];
rz(-pi) q[2];
rz(-0.42550605) q[3];
sx q[3];
rz(-1.7153499) q[3];
sx q[3];
rz(1.0021462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.17710182) q[2];
sx q[2];
rz(-1.0937546) q[2];
sx q[2];
rz(-2.2237681) q[2];
rz(-1.142189) q[3];
sx q[3];
rz(-0.85421383) q[3];
sx q[3];
rz(0.93723047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98638242) q[0];
sx q[0];
rz(-2.4676403) q[0];
sx q[0];
rz(-2.881799) q[0];
rz(-2.4329176) q[1];
sx q[1];
rz(-0.27888137) q[1];
sx q[1];
rz(-0.52694595) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5848815) q[0];
sx q[0];
rz(-1.3538133) q[0];
sx q[0];
rz(2.7816539) q[0];
rz(0.2692659) q[2];
sx q[2];
rz(-2.2109647) q[2];
sx q[2];
rz(-2.7562041) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.38339112) q[1];
sx q[1];
rz(-1.8806629) q[1];
sx q[1];
rz(1.8499225) q[1];
x q[2];
rz(-1.8606887) q[3];
sx q[3];
rz(-1.999951) q[3];
sx q[3];
rz(-3.1115301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8156585) q[2];
sx q[2];
rz(-2.2951173) q[2];
sx q[2];
rz(-0.79088598) q[2];
rz(2.7549426) q[3];
sx q[3];
rz(-2.1270042) q[3];
sx q[3];
rz(-0.33716831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7336422) q[0];
sx q[0];
rz(-2.9691417) q[0];
sx q[0];
rz(0.98544425) q[0];
rz(2.573029) q[1];
sx q[1];
rz(-1.111258) q[1];
sx q[1];
rz(0.3607761) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6919149) q[0];
sx q[0];
rz(-1.5957818) q[0];
sx q[0];
rz(1.7382394) q[0];
rz(-0.12014328) q[2];
sx q[2];
rz(-1.7842245) q[2];
sx q[2];
rz(-0.10211589) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4611778) q[1];
sx q[1];
rz(-2.3056681) q[1];
sx q[1];
rz(1.7635036) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1694069) q[3];
sx q[3];
rz(-1.7161887) q[3];
sx q[3];
rz(-1.8190847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0439904) q[2];
sx q[2];
rz(-2.4520935) q[2];
sx q[2];
rz(0.94341755) q[2];
rz(-0.57389456) q[3];
sx q[3];
rz(-0.4807764) q[3];
sx q[3];
rz(-1.7452314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5744793) q[0];
sx q[0];
rz(-1.6945101) q[0];
sx q[0];
rz(2.2816818) q[0];
rz(1.7815331) q[1];
sx q[1];
rz(-2.139745) q[1];
sx q[1];
rz(2.3812961) q[1];
rz(2.0104682) q[2];
sx q[2];
rz(-0.22618539) q[2];
sx q[2];
rz(-2.4629081) q[2];
rz(2.3902262) q[3];
sx q[3];
rz(-2.1204227) q[3];
sx q[3];
rz(-1.942996) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];