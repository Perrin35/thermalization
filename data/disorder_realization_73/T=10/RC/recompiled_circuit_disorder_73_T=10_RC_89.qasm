OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1459382) q[0];
sx q[0];
rz(-2.6383658) q[0];
sx q[0];
rz(0.72416645) q[0];
rz(0.63996285) q[1];
sx q[1];
rz(-0.53007403) q[1];
sx q[1];
rz(-0.78483265) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9871702) q[0];
sx q[0];
rz(-1.3599456) q[0];
sx q[0];
rz(0.2984557) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7230524) q[2];
sx q[2];
rz(-1.4782018) q[2];
sx q[2];
rz(1.4946403) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9468294) q[1];
sx q[1];
rz(-2.0954872) q[1];
sx q[1];
rz(2.8835433) q[1];
x q[2];
rz(-0.033693245) q[3];
sx q[3];
rz(-1.1564621) q[3];
sx q[3];
rz(-1.8309483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5518387) q[2];
sx q[2];
rz(-1.4171615) q[2];
sx q[2];
rz(-3.0736249) q[2];
rz(-3.0170278) q[3];
sx q[3];
rz(-0.3228651) q[3];
sx q[3];
rz(1.386806) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92000604) q[0];
sx q[0];
rz(-0.13555549) q[0];
sx q[0];
rz(2.8979229) q[0];
rz(-2.5098353) q[1];
sx q[1];
rz(-1.4032204) q[1];
sx q[1];
rz(1.7858645) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5403554) q[0];
sx q[0];
rz(-1.3155126) q[0];
sx q[0];
rz(-3.113494) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.182105) q[2];
sx q[2];
rz(-1.9027862) q[2];
sx q[2];
rz(1.7322844) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4888549) q[1];
sx q[1];
rz(-1.4421717) q[1];
sx q[1];
rz(-2.0205523) q[1];
x q[2];
rz(2.5494266) q[3];
sx q[3];
rz(-1.5321931) q[3];
sx q[3];
rz(1.6174699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0791066) q[2];
sx q[2];
rz(-2.1495543) q[2];
sx q[2];
rz(-2.8919354) q[2];
rz(2.6349973) q[3];
sx q[3];
rz(-1.5157615) q[3];
sx q[3];
rz(0.33199582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24519414) q[0];
sx q[0];
rz(-1.2250552) q[0];
sx q[0];
rz(-2.2431592) q[0];
rz(-1.8067182) q[1];
sx q[1];
rz(-1.9060262) q[1];
sx q[1];
rz(1.2737087) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1184517) q[0];
sx q[0];
rz(-1.8225192) q[0];
sx q[0];
rz(1.9385507) q[0];
rz(1.4726228) q[2];
sx q[2];
rz(-1.86491) q[2];
sx q[2];
rz(1.836118) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9540005) q[1];
sx q[1];
rz(-2.2224269) q[1];
sx q[1];
rz(-1.1263532) q[1];
x q[2];
rz(2.5898315) q[3];
sx q[3];
rz(-2.3380087) q[3];
sx q[3];
rz(1.4078275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0818103) q[2];
sx q[2];
rz(-2.5368097) q[2];
sx q[2];
rz(-0.95345062) q[2];
rz(3.1070784) q[3];
sx q[3];
rz(-0.78648609) q[3];
sx q[3];
rz(-0.22687337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8614486) q[0];
sx q[0];
rz(-2.9252958) q[0];
sx q[0];
rz(0.24818534) q[0];
rz(2.10363) q[1];
sx q[1];
rz(-1.1231517) q[1];
sx q[1];
rz(3.0674556) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.603133) q[0];
sx q[0];
rz(-2.5939301) q[0];
sx q[0];
rz(-1.0663701) q[0];
x q[1];
rz(0.3481625) q[2];
sx q[2];
rz(-2.2194038) q[2];
sx q[2];
rz(1.6096889) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.76064199) q[1];
sx q[1];
rz(-1.4529072) q[1];
sx q[1];
rz(1.2537434) q[1];
x q[2];
rz(-1.3965963) q[3];
sx q[3];
rz(-1.8074236) q[3];
sx q[3];
rz(-0.06121204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8903824) q[2];
sx q[2];
rz(-0.40428287) q[2];
sx q[2];
rz(0.038643535) q[2];
rz(2.1679227) q[3];
sx q[3];
rz(-2.645851) q[3];
sx q[3];
rz(-2.8715449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53428179) q[0];
sx q[0];
rz(-1.6058291) q[0];
sx q[0];
rz(-1.779153) q[0];
rz(-0.81659395) q[1];
sx q[1];
rz(-1.2885619) q[1];
sx q[1];
rz(-1.1626676) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8602596) q[0];
sx q[0];
rz(-1.5409924) q[0];
sx q[0];
rz(-1.4593065) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2478998) q[2];
sx q[2];
rz(-1.731551) q[2];
sx q[2];
rz(1.9184743) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3711277) q[1];
sx q[1];
rz(-0.41509291) q[1];
sx q[1];
rz(2.0626555) q[1];
rz(-0.43443067) q[3];
sx q[3];
rz(-2.0475004) q[3];
sx q[3];
rz(1.3694976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6267307) q[2];
sx q[2];
rz(-2.0613487) q[2];
sx q[2];
rz(-0.14222063) q[2];
rz(2.2375315) q[3];
sx q[3];
rz(-1.3198493) q[3];
sx q[3];
rz(0.18946762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11480039) q[0];
sx q[0];
rz(-2.8650706) q[0];
sx q[0];
rz(-1.6739155) q[0];
rz(2.5698075) q[1];
sx q[1];
rz(-2.7829058) q[1];
sx q[1];
rz(-2.8335559) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.98826) q[0];
sx q[0];
rz(-1.4446265) q[0];
sx q[0];
rz(-1.6660965) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.61200895) q[2];
sx q[2];
rz(-1.0527805) q[2];
sx q[2];
rz(0.81105622) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8923556) q[1];
sx q[1];
rz(-2.4294469) q[1];
sx q[1];
rz(1.874079) q[1];
rz(-pi) q[2];
rz(1.4792535) q[3];
sx q[3];
rz(-1.2237826) q[3];
sx q[3];
rz(0.43859827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3623111) q[2];
sx q[2];
rz(-1.4004536) q[2];
sx q[2];
rz(-1.1479088) q[2];
rz(2.4273196) q[3];
sx q[3];
rz(-0.89759421) q[3];
sx q[3];
rz(2.4961297) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66184735) q[0];
sx q[0];
rz(-0.84091887) q[0];
sx q[0];
rz(3.0116144) q[0];
rz(-0.030844363) q[1];
sx q[1];
rz(-1.8519311) q[1];
sx q[1];
rz(-0.67108363) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2790047) q[0];
sx q[0];
rz(-0.064084856) q[0];
sx q[0];
rz(-2.5628753) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1110711) q[2];
sx q[2];
rz(-1.3866716) q[2];
sx q[2];
rz(-0.39779278) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1783501) q[1];
sx q[1];
rz(-1.9415783) q[1];
sx q[1];
rz(3.0313655) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.047366553) q[3];
sx q[3];
rz(-0.8960552) q[3];
sx q[3];
rz(-0.55344068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0188296) q[2];
sx q[2];
rz(-1.5315703) q[2];
sx q[2];
rz(0.57787952) q[2];
rz(3.1130062) q[3];
sx q[3];
rz(-1.280602) q[3];
sx q[3];
rz(1.8813429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-1.1639444) q[0];
sx q[0];
rz(-2.2387235) q[0];
sx q[0];
rz(0.41241616) q[0];
rz(-1.4498129) q[1];
sx q[1];
rz(-1.7990566) q[1];
sx q[1];
rz(-1.1669881) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1056571) q[0];
sx q[0];
rz(-1.3577537) q[0];
sx q[0];
rz(-0.94353326) q[0];
x q[1];
rz(-0.27600482) q[2];
sx q[2];
rz(-0.88062421) q[2];
sx q[2];
rz(-1.3921757) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3762714) q[1];
sx q[1];
rz(-1.5404535) q[1];
sx q[1];
rz(-2.9529851) q[1];
rz(-pi) q[2];
rz(0.1498296) q[3];
sx q[3];
rz(-2.095788) q[3];
sx q[3];
rz(2.6168952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2010487) q[2];
sx q[2];
rz(-0.87934914) q[2];
sx q[2];
rz(-2.9525625) q[2];
rz(-2.9947301) q[3];
sx q[3];
rz(-0.18460128) q[3];
sx q[3];
rz(-1.7485025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1259595) q[0];
sx q[0];
rz(-1.8305612) q[0];
sx q[0];
rz(-0.92700672) q[0];
rz(-1.3828297) q[1];
sx q[1];
rz(-0.60950509) q[1];
sx q[1];
rz(-1.6519201) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4091858) q[0];
sx q[0];
rz(-0.86940765) q[0];
sx q[0];
rz(-0.62774815) q[0];
rz(-1.4719047) q[2];
sx q[2];
rz(-2.6675468) q[2];
sx q[2];
rz(2.3941819) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.97396321) q[1];
sx q[1];
rz(-1.7427674) q[1];
sx q[1];
rz(2.9457438) q[1];
x q[2];
rz(-2.1065815) q[3];
sx q[3];
rz(-1.1318558) q[3];
sx q[3];
rz(2.5715695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5902517) q[2];
sx q[2];
rz(-0.44289032) q[2];
sx q[2];
rz(1.7112973) q[2];
rz(-2.5643505) q[3];
sx q[3];
rz(-2.2539299) q[3];
sx q[3];
rz(1.3841217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5883314) q[0];
sx q[0];
rz(-1.7734779) q[0];
sx q[0];
rz(-0.28840315) q[0];
rz(-2.6092031) q[1];
sx q[1];
rz(-0.45982292) q[1];
sx q[1];
rz(0.14702252) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0518236) q[0];
sx q[0];
rz(-0.71634403) q[0];
sx q[0];
rz(-0.96869529) q[0];
rz(-pi) q[1];
x q[1];
rz(1.64098) q[2];
sx q[2];
rz(-1.3314221) q[2];
sx q[2];
rz(-1.0603051) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8204931) q[1];
sx q[1];
rz(-1.5309146) q[1];
sx q[1];
rz(0.95221968) q[1];
x q[2];
rz(0.23707323) q[3];
sx q[3];
rz(-2.5519538) q[3];
sx q[3];
rz(-1.8001363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0599351) q[2];
sx q[2];
rz(-0.43490484) q[2];
sx q[2];
rz(-2.3975513) q[2];
rz(-2.384281) q[3];
sx q[3];
rz(-1.7777187) q[3];
sx q[3];
rz(-1.3967167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1159146) q[0];
sx q[0];
rz(-2.0712576) q[0];
sx q[0];
rz(2.0448137) q[0];
rz(2.3241282) q[1];
sx q[1];
rz(-1.2066963) q[1];
sx q[1];
rz(-0.6304601) q[1];
rz(0.4521162) q[2];
sx q[2];
rz(-1.631626) q[2];
sx q[2];
rz(0.22125868) q[2];
rz(2.1327303) q[3];
sx q[3];
rz(-1.4603793) q[3];
sx q[3];
rz(0.59752656) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
