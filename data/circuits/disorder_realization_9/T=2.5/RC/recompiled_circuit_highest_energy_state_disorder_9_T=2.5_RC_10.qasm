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
rz(3.0577793) q[0];
sx q[0];
rz(-0.35141355) q[0];
sx q[0];
rz(1.0454398) q[0];
rz(-5.4391556) q[1];
sx q[1];
rz(0.56517711) q[1];
sx q[1];
rz(11.204389) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0365252) q[0];
sx q[0];
rz(-2.1224667) q[0];
sx q[0];
rz(1.9776634) q[0];
x q[1];
rz(-1.7704813) q[2];
sx q[2];
rz(-2.7058995) q[2];
sx q[2];
rz(-1.0267804) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.61698422) q[1];
sx q[1];
rz(-1.1676663) q[1];
sx q[1];
rz(1.8126081) q[1];
rz(-pi) q[2];
rz(2.1426422) q[3];
sx q[3];
rz(-0.57724726) q[3];
sx q[3];
rz(1.5791864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7094946) q[2];
sx q[2];
rz(-0.39958909) q[2];
sx q[2];
rz(0.70322767) q[2];
rz(0.75913366) q[3];
sx q[3];
rz(-2.1742564) q[3];
sx q[3];
rz(0.67559344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9264939) q[0];
sx q[0];
rz(-0.7248942) q[0];
sx q[0];
rz(2.3469927) q[0];
rz(2.3130747) q[1];
sx q[1];
rz(-1.0732032) q[1];
sx q[1];
rz(2.6603096) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3743418) q[0];
sx q[0];
rz(-1.4515948) q[0];
sx q[0];
rz(1.5990555) q[0];
x q[1];
rz(-2.6120466) q[2];
sx q[2];
rz(-0.82223071) q[2];
sx q[2];
rz(-1.032743) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1254639) q[1];
sx q[1];
rz(-1.4939346) q[1];
sx q[1];
rz(1.8015339) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0167502) q[3];
sx q[3];
rz(-0.9940799) q[3];
sx q[3];
rz(-2.5050636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1543697) q[2];
sx q[2];
rz(-2.1798446) q[2];
sx q[2];
rz(2.6340384) q[2];
rz(-2.4842723) q[3];
sx q[3];
rz(-2.1828914) q[3];
sx q[3];
rz(1.2348194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2051314) q[0];
sx q[0];
rz(-1.5345804) q[0];
sx q[0];
rz(-2.8954647) q[0];
rz(2.4283465) q[1];
sx q[1];
rz(-1.617022) q[1];
sx q[1];
rz(2.2770142) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1522035) q[0];
sx q[0];
rz(-1.8789904) q[0];
sx q[0];
rz(-0.21580036) q[0];
x q[1];
rz(2.6106493) q[2];
sx q[2];
rz(-2.1111672) q[2];
sx q[2];
rz(-2.2744389) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9731701) q[1];
sx q[1];
rz(-1.6515141) q[1];
sx q[1];
rz(1.1087085) q[1];
rz(-pi) q[2];
rz(-2.9900486) q[3];
sx q[3];
rz(-2.3484548) q[3];
sx q[3];
rz(1.053699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7264709) q[2];
sx q[2];
rz(-0.39602009) q[2];
sx q[2];
rz(-1.9351287) q[2];
rz(0.11780277) q[3];
sx q[3];
rz(-0.67540568) q[3];
sx q[3];
rz(0.12731586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29663157) q[0];
sx q[0];
rz(-1.4609818) q[0];
sx q[0];
rz(-2.6764181) q[0];
rz(-2.2370715) q[1];
sx q[1];
rz(-2.1280839) q[1];
sx q[1];
rz(1.4314338) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2003882) q[0];
sx q[0];
rz(-2.6058307) q[0];
sx q[0];
rz(-0.8922116) q[0];
rz(1.2063556) q[2];
sx q[2];
rz(-0.52204692) q[2];
sx q[2];
rz(-2.4625157) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5618973) q[1];
sx q[1];
rz(-2.2219147) q[1];
sx q[1];
rz(-1.1267594) q[1];
rz(-1.9172098) q[3];
sx q[3];
rz(-0.87408057) q[3];
sx q[3];
rz(-2.5908366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5005834) q[2];
sx q[2];
rz(-2.9312129) q[2];
sx q[2];
rz(-1.496199) q[2];
rz(0.03838852) q[3];
sx q[3];
rz(-1.8206785) q[3];
sx q[3];
rz(-0.7588318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0043623) q[0];
sx q[0];
rz(-0.69843233) q[0];
sx q[0];
rz(-1.5414365) q[0];
rz(-1.6372708) q[1];
sx q[1];
rz(-1.5802822) q[1];
sx q[1];
rz(-1.5024332) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7961032) q[0];
sx q[0];
rz(-1.1123344) q[0];
sx q[0];
rz(-1.0572516) q[0];
rz(-pi) q[1];
rz(0.2404332) q[2];
sx q[2];
rz(-0.56686646) q[2];
sx q[2];
rz(0.64449233) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0509661) q[1];
sx q[1];
rz(-0.6399261) q[1];
sx q[1];
rz(-0.76211318) q[1];
rz(-pi) q[2];
x q[2];
rz(0.10879119) q[3];
sx q[3];
rz(-1.5034893) q[3];
sx q[3];
rz(1.045026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1372823) q[2];
sx q[2];
rz(-2.8305125) q[2];
sx q[2];
rz(-2.7610371) q[2];
rz(-0.54008326) q[3];
sx q[3];
rz(-1.249908) q[3];
sx q[3];
rz(1.165747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49749097) q[0];
sx q[0];
rz(-3.0971425) q[0];
sx q[0];
rz(-0.38462001) q[0];
rz(-3.0907471) q[1];
sx q[1];
rz(-0.47639242) q[1];
sx q[1];
rz(-1.6772038) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58005262) q[0];
sx q[0];
rz(-2.0140352) q[0];
sx q[0];
rz(1.1053941) q[0];
rz(-pi) q[1];
rz(-2.2273686) q[2];
sx q[2];
rz(-2.3434635) q[2];
sx q[2];
rz(-1.4887336) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.205464) q[1];
sx q[1];
rz(-1.4805321) q[1];
sx q[1];
rz(1.4858264) q[1];
rz(-pi) q[2];
rz(-1.3648273) q[3];
sx q[3];
rz(-2.245083) q[3];
sx q[3];
rz(-2.8482361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9836318) q[2];
sx q[2];
rz(-0.76405683) q[2];
sx q[2];
rz(0.91510406) q[2];
rz(-0.30465952) q[3];
sx q[3];
rz(-0.72158486) q[3];
sx q[3];
rz(1.4753531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-1.7810998) q[0];
sx q[0];
rz(-2.2486794) q[0];
sx q[0];
rz(1.1156981) q[0];
rz(-0.64104331) q[1];
sx q[1];
rz(-1.2572925) q[1];
sx q[1];
rz(-1.8152274) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4531739) q[0];
sx q[0];
rz(-2.6368615) q[0];
sx q[0];
rz(-0.27413989) q[0];
rz(0.4037598) q[2];
sx q[2];
rz(-0.40382622) q[2];
sx q[2];
rz(-0.60690875) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.54394447) q[1];
sx q[1];
rz(-0.94606384) q[1];
sx q[1];
rz(0.90170699) q[1];
rz(-2.2828045) q[3];
sx q[3];
rz(-0.988171) q[3];
sx q[3];
rz(-1.4020385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5850087) q[2];
sx q[2];
rz(-0.98238397) q[2];
sx q[2];
rz(-1.8741685) q[2];
rz(-0.061138717) q[3];
sx q[3];
rz(-1.2631402) q[3];
sx q[3];
rz(0.88501969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1214445) q[0];
sx q[0];
rz(-0.76524884) q[0];
sx q[0];
rz(-2.6159317) q[0];
rz(-1.4875745) q[1];
sx q[1];
rz(-2.0272389) q[1];
sx q[1];
rz(0.18255998) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17815264) q[0];
sx q[0];
rz(-1.9885673) q[0];
sx q[0];
rz(-1.3304292) q[0];
rz(-pi) q[1];
rz(2.0261954) q[2];
sx q[2];
rz(-1.2898852) q[2];
sx q[2];
rz(-1.8473491) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9801431) q[1];
sx q[1];
rz(-1.7421466) q[1];
sx q[1];
rz(-0.93659393) q[1];
rz(-pi) q[2];
rz(-2.0917039) q[3];
sx q[3];
rz(-2.6080797) q[3];
sx q[3];
rz(-0.81797879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.21060264) q[2];
sx q[2];
rz(-1.3347722) q[2];
sx q[2];
rz(-0.1087428) q[2];
rz(-3.0336618) q[3];
sx q[3];
rz(-2.7863672) q[3];
sx q[3];
rz(1.3817374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0906319) q[0];
sx q[0];
rz(-1.9405631) q[0];
sx q[0];
rz(0.65318024) q[0];
rz(1.2921035) q[1];
sx q[1];
rz(-2.8958246) q[1];
sx q[1];
rz(-0.076315708) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6298593) q[0];
sx q[0];
rz(-2.7329067) q[0];
sx q[0];
rz(-2.2124394) q[0];
x q[1];
rz(0.033360783) q[2];
sx q[2];
rz(-1.4028869) q[2];
sx q[2];
rz(-2.1176934) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2321978) q[1];
sx q[1];
rz(-1.4320717) q[1];
sx q[1];
rz(-0.75896427) q[1];
rz(-pi) q[2];
rz(-1.1006484) q[3];
sx q[3];
rz(-2.4565434) q[3];
sx q[3];
rz(-2.4026826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.896758) q[2];
sx q[2];
rz(-1.1847757) q[2];
sx q[2];
rz(-2.1484788) q[2];
rz(1.7874329) q[3];
sx q[3];
rz(-0.5391776) q[3];
sx q[3];
rz(-1.7206934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83574522) q[0];
sx q[0];
rz(-1.5033686) q[0];
sx q[0];
rz(0.29356965) q[0];
rz(1.0319895) q[1];
sx q[1];
rz(-0.76621619) q[1];
sx q[1];
rz(-2.8242677) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.246884) q[0];
sx q[0];
rz(-0.63011677) q[0];
sx q[0];
rz(0.21960857) q[0];
rz(-pi) q[1];
rz(-0.34658771) q[2];
sx q[2];
rz(-2.5591597) q[2];
sx q[2];
rz(2.1726441) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.369243) q[1];
sx q[1];
rz(-1.1398106) q[1];
sx q[1];
rz(1.1915156) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6280319) q[3];
sx q[3];
rz(-1.2737107) q[3];
sx q[3];
rz(-0.12606584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.62852922) q[2];
sx q[2];
rz(-2.8184012) q[2];
sx q[2];
rz(0.31488669) q[2];
rz(3.1079187) q[3];
sx q[3];
rz(-1.0957054) q[3];
sx q[3];
rz(-1.7033738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2122129) q[0];
sx q[0];
rz(-2.0912981) q[0];
sx q[0];
rz(-0.97846497) q[0];
rz(2.9411511) q[1];
sx q[1];
rz(-1.0841752) q[1];
sx q[1];
rz(-0.60842327) q[1];
rz(0.048439518) q[2];
sx q[2];
rz(-1.6002527) q[2];
sx q[2];
rz(-3.0839828) q[2];
rz(-1.8157806) q[3];
sx q[3];
rz(-1.0272588) q[3];
sx q[3];
rz(2.5486265) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
