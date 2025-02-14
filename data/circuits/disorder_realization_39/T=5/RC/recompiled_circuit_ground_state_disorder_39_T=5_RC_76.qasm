OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.47867632) q[0];
sx q[0];
rz(7.8362099) q[0];
sx q[0];
rz(10.683164) q[0];
rz(-5.0212669) q[1];
sx q[1];
rz(6.8016383) q[1];
sx q[1];
rz(6.3432884) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10190554) q[0];
sx q[0];
rz(-1.6015669) q[0];
sx q[0];
rz(-3.0512848) q[0];
x q[1];
rz(2.8608039) q[2];
sx q[2];
rz(-2.6323689) q[2];
sx q[2];
rz(1.4851242) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2138252) q[1];
sx q[1];
rz(-2.7672184) q[1];
sx q[1];
rz(-1.8273749) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46334907) q[3];
sx q[3];
rz(-1.0245205) q[3];
sx q[3];
rz(-1.1322556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.95187676) q[2];
sx q[2];
rz(-0.78236255) q[2];
sx q[2];
rz(-2.961109) q[2];
rz(2.8372724) q[3];
sx q[3];
rz(-0.92420998) q[3];
sx q[3];
rz(2.9144104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75905281) q[0];
sx q[0];
rz(-2.5681684) q[0];
sx q[0];
rz(0.10391129) q[0];
rz(0.78481627) q[1];
sx q[1];
rz(-1.3094614) q[1];
sx q[1];
rz(-2.5596502) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73954813) q[0];
sx q[0];
rz(-0.53240896) q[0];
sx q[0];
rz(-0.043921434) q[0];
rz(0.50931661) q[2];
sx q[2];
rz(-1.881372) q[2];
sx q[2];
rz(0.79307014) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2375029) q[1];
sx q[1];
rz(-2.6640035) q[1];
sx q[1];
rz(1.6271923) q[1];
rz(-pi) q[2];
x q[2];
rz(2.373967) q[3];
sx q[3];
rz(-2.5968142) q[3];
sx q[3];
rz(-1.9035853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4265784) q[2];
sx q[2];
rz(-0.38807401) q[2];
sx q[2];
rz(-1.0283872) q[2];
rz(-0.55108023) q[3];
sx q[3];
rz(-1.9469399) q[3];
sx q[3];
rz(-2.6330131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5904163) q[0];
sx q[0];
rz(-1.6352147) q[0];
sx q[0];
rz(-1.2580385) q[0];
rz(2.538077) q[1];
sx q[1];
rz(-1.8305093) q[1];
sx q[1];
rz(-2.9579732) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70442536) q[0];
sx q[0];
rz(-1.9806304) q[0];
sx q[0];
rz(2.8528575) q[0];
x q[1];
rz(-0.48086353) q[2];
sx q[2];
rz(-1.4521952) q[2];
sx q[2];
rz(3.1080217) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0407942) q[1];
sx q[1];
rz(-1.9996627) q[1];
sx q[1];
rz(2.3593192) q[1];
x q[2];
rz(-0.19422005) q[3];
sx q[3];
rz(-2.9306539) q[3];
sx q[3];
rz(2.6489182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.52643481) q[2];
sx q[2];
rz(-1.3765455) q[2];
sx q[2];
rz(0.60205013) q[2];
rz(0.43271068) q[3];
sx q[3];
rz(-1.5475464) q[3];
sx q[3];
rz(-1.4759147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3848307) q[0];
sx q[0];
rz(-0.47465208) q[0];
sx q[0];
rz(-2.5153644) q[0];
rz(1.2681883) q[1];
sx q[1];
rz(-1.405193) q[1];
sx q[1];
rz(0.38571206) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2572266) q[0];
sx q[0];
rz(-0.8943087) q[0];
sx q[0];
rz(2.1348597) q[0];
rz(2.0247632) q[2];
sx q[2];
rz(-1.0359284) q[2];
sx q[2];
rz(-1.3605538) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0929132) q[1];
sx q[1];
rz(-1.5041222) q[1];
sx q[1];
rz(-2.9408022) q[1];
rz(-pi) q[2];
rz(-3.0195531) q[3];
sx q[3];
rz(-1.0309848) q[3];
sx q[3];
rz(0.5468747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6993774) q[2];
sx q[2];
rz(-1.135004) q[2];
sx q[2];
rz(-0.29083148) q[2];
rz(-1.95131) q[3];
sx q[3];
rz(-0.61605993) q[3];
sx q[3];
rz(-2.0975838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46665835) q[0];
sx q[0];
rz(-0.062677296) q[0];
sx q[0];
rz(1.4516996) q[0];
rz(-0.85707227) q[1];
sx q[1];
rz(-1.2985726) q[1];
sx q[1];
rz(0.42795408) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65954094) q[0];
sx q[0];
rz(-2.6845346) q[0];
sx q[0];
rz(1.2808475) q[0];
rz(-0.73889795) q[2];
sx q[2];
rz(-1.9618926) q[2];
sx q[2];
rz(1.0223688) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.89558307) q[1];
sx q[1];
rz(-1.9253139) q[1];
sx q[1];
rz(-0.73486967) q[1];
rz(2.1367504) q[3];
sx q[3];
rz(-2.1814846) q[3];
sx q[3];
rz(0.37328675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2032623) q[2];
sx q[2];
rz(-2.0945956) q[2];
sx q[2];
rz(-0.14673512) q[2];
rz(-0.19715582) q[3];
sx q[3];
rz(-1.4898172) q[3];
sx q[3];
rz(2.3728235) q[3];
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
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7305304) q[0];
sx q[0];
rz(-1.0163607) q[0];
sx q[0];
rz(-2.9732669) q[0];
rz(2.7492145) q[1];
sx q[1];
rz(-0.43586755) q[1];
sx q[1];
rz(-1.6900774) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.096286721) q[0];
sx q[0];
rz(-2.0074962) q[0];
sx q[0];
rz(3.1179881) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8658357) q[2];
sx q[2];
rz(-2.3988798) q[2];
sx q[2];
rz(-2.4969522) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.299356) q[1];
sx q[1];
rz(-0.82932392) q[1];
sx q[1];
rz(1.7620128) q[1];
x q[2];
rz(-2.3528162) q[3];
sx q[3];
rz(-1.2433854) q[3];
sx q[3];
rz(-0.57306266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3671941) q[2];
sx q[2];
rz(-0.44782475) q[2];
sx q[2];
rz(1.9642584) q[2];
rz(-2.2199953) q[3];
sx q[3];
rz(-0.84601837) q[3];
sx q[3];
rz(-2.6109076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0717936) q[0];
sx q[0];
rz(-1.8853747) q[0];
sx q[0];
rz(2.0846833) q[0];
rz(2.9131556) q[1];
sx q[1];
rz(-0.7799131) q[1];
sx q[1];
rz(-2.8905919) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91755501) q[0];
sx q[0];
rz(-2.7537808) q[0];
sx q[0];
rz(2.4075131) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2784581) q[2];
sx q[2];
rz(-2.302711) q[2];
sx q[2];
rz(-1.2108491) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9343379) q[1];
sx q[1];
rz(-1.7060857) q[1];
sx q[1];
rz(2.9771718) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5369473) q[3];
sx q[3];
rz(-0.6007768) q[3];
sx q[3];
rz(2.0913817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.935427) q[2];
sx q[2];
rz(-0.82038227) q[2];
sx q[2];
rz(0.41681918) q[2];
rz(-0.56944141) q[3];
sx q[3];
rz(-1.1007525) q[3];
sx q[3];
rz(2.4429564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5702629) q[0];
sx q[0];
rz(-1.2705734) q[0];
sx q[0];
rz(0.50462333) q[0];
rz(2.8847671) q[1];
sx q[1];
rz(-2.3378614) q[1];
sx q[1];
rz(-1.9907192) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4537612) q[0];
sx q[0];
rz(-2.6021299) q[0];
sx q[0];
rz(-1.1629172) q[0];
rz(-pi) q[1];
x q[1];
rz(0.35607349) q[2];
sx q[2];
rz(-0.96637539) q[2];
sx q[2];
rz(-0.52456896) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0432277) q[1];
sx q[1];
rz(-0.72611626) q[1];
sx q[1];
rz(-0.67732201) q[1];
x q[2];
rz(-0.41696291) q[3];
sx q[3];
rz(-1.7156258) q[3];
sx q[3];
rz(-2.8944998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5095832) q[2];
sx q[2];
rz(-2.7183967) q[2];
sx q[2];
rz(0.43935856) q[2];
rz(-1.1257233) q[3];
sx q[3];
rz(-1.537354) q[3];
sx q[3];
rz(1.841898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2706962) q[0];
sx q[0];
rz(-2.3211711) q[0];
sx q[0];
rz(2.6743555) q[0];
rz(-1.0386508) q[1];
sx q[1];
rz(-0.50913441) q[1];
sx q[1];
rz(-1.4899303) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.576173) q[0];
sx q[0];
rz(-1.7459946) q[0];
sx q[0];
rz(-2.6433667) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9949399) q[2];
sx q[2];
rz(-1.6645558) q[2];
sx q[2];
rz(2.3921426) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6257621) q[1];
sx q[1];
rz(-2.3861679) q[1];
sx q[1];
rz(-1.6095407) q[1];
x q[2];
rz(-2.5760004) q[3];
sx q[3];
rz(-1.1882236) q[3];
sx q[3];
rz(2.3747834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.44522875) q[2];
sx q[2];
rz(-2.7859521) q[2];
sx q[2];
rz(-1.5544372) q[2];
rz(-2.1571531) q[3];
sx q[3];
rz(-2.2319904) q[3];
sx q[3];
rz(1.4136774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.5928891) q[0];
sx q[0];
rz(-2.2991572) q[0];
sx q[0];
rz(-0.62826759) q[0];
rz(-2.9382622) q[1];
sx q[1];
rz(-2.0275828) q[1];
sx q[1];
rz(0.59439739) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0258163) q[0];
sx q[0];
rz(-2.5133142) q[0];
sx q[0];
rz(-1.9253057) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7901232) q[2];
sx q[2];
rz(-1.0426874) q[2];
sx q[2];
rz(2.2161432) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0815142) q[1];
sx q[1];
rz(-1.8412672) q[1];
sx q[1];
rz(-0.31281506) q[1];
rz(1.6420308) q[3];
sx q[3];
rz(-0.95893493) q[3];
sx q[3];
rz(0.16779403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2079042) q[2];
sx q[2];
rz(-2.4173357) q[2];
sx q[2];
rz(-2.9653463) q[2];
rz(1.3147973) q[3];
sx q[3];
rz(-2.3254471) q[3];
sx q[3];
rz(-0.76752457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15277302) q[0];
sx q[0];
rz(-1.7012699) q[0];
sx q[0];
rz(-1.758601) q[0];
rz(3.026961) q[1];
sx q[1];
rz(-0.75420598) q[1];
sx q[1];
rz(2.4892714) q[1];
rz(2.9822275) q[2];
sx q[2];
rz(-2.9014316) q[2];
sx q[2];
rz(0.17221774) q[2];
rz(2.7588941) q[3];
sx q[3];
rz(-2.0321587) q[3];
sx q[3];
rz(1.7176499) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
