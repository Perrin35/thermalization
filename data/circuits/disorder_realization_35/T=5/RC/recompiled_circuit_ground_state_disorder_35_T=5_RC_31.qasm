OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72260296) q[0];
sx q[0];
rz(-2.498772) q[0];
sx q[0];
rz(8.2052054) q[0];
rz(0.83528432) q[1];
sx q[1];
rz(-1.3568027) q[1];
sx q[1];
rz(-2.2396483) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.886698) q[0];
sx q[0];
rz(-0.6318081) q[0];
sx q[0];
rz(2.8796701) q[0];
rz(-pi) q[1];
rz(-0.44166126) q[2];
sx q[2];
rz(-0.41641372) q[2];
sx q[2];
rz(-3.1232782) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.28022596) q[1];
sx q[1];
rz(-1.8609443) q[1];
sx q[1];
rz(-2.9600083) q[1];
rz(-pi) q[2];
rz(-0.91633032) q[3];
sx q[3];
rz(-0.95312762) q[3];
sx q[3];
rz(0.82181069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.96282643) q[2];
sx q[2];
rz(-2.1361394) q[2];
sx q[2];
rz(-2.315305) q[2];
rz(-1.4055584) q[3];
sx q[3];
rz(-2.8345351) q[3];
sx q[3];
rz(-2.3981986) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58251441) q[0];
sx q[0];
rz(-2.7358416) q[0];
sx q[0];
rz(-1.4003117) q[0];
rz(-1.1236313) q[1];
sx q[1];
rz(-0.39389899) q[1];
sx q[1];
rz(-2.0434911) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.708824) q[0];
sx q[0];
rz(-1.4147804) q[0];
sx q[0];
rz(-1.6814598) q[0];
x q[1];
rz(0.33771659) q[2];
sx q[2];
rz(-1.8840839) q[2];
sx q[2];
rz(-0.17853436) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.73126331) q[1];
sx q[1];
rz(-1.6917233) q[1];
sx q[1];
rz(-1.6464473) q[1];
x q[2];
rz(1.9722749) q[3];
sx q[3];
rz(-0.40630925) q[3];
sx q[3];
rz(3.0412948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0565679) q[2];
sx q[2];
rz(-0.3173863) q[2];
sx q[2];
rz(-1.2497831) q[2];
rz(-1.7791087) q[3];
sx q[3];
rz(-1.0008078) q[3];
sx q[3];
rz(-1.484163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.5054841) q[0];
sx q[0];
rz(-2.4478069) q[0];
sx q[0];
rz(0.27534494) q[0];
rz(-0.45922008) q[1];
sx q[1];
rz(-0.95454916) q[1];
sx q[1];
rz(3.088248) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8109574) q[0];
sx q[0];
rz(-0.22718469) q[0];
sx q[0];
rz(-2.7869528) q[0];
rz(-pi) q[1];
rz(-0.01845999) q[2];
sx q[2];
rz(-2.0316796) q[2];
sx q[2];
rz(2.8783952) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.054823067) q[1];
sx q[1];
rz(-0.55332282) q[1];
sx q[1];
rz(-3.064108) q[1];
rz(-pi) q[2];
rz(1.0166078) q[3];
sx q[3];
rz(-2.7597868) q[3];
sx q[3];
rz(3.0757381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3205545) q[2];
sx q[2];
rz(-2.7802763) q[2];
sx q[2];
rz(1.6374755) q[2];
rz(-1.5004246) q[3];
sx q[3];
rz(-2.0913561) q[3];
sx q[3];
rz(-2.8659081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7212873) q[0];
sx q[0];
rz(-2.5581701) q[0];
sx q[0];
rz(-2.6570008) q[0];
rz(-1.8991607) q[1];
sx q[1];
rz(-2.395605) q[1];
sx q[1];
rz(-2.7925083) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9681681) q[0];
sx q[0];
rz(-1.5281046) q[0];
sx q[0];
rz(-1.0220362) q[0];
x q[1];
rz(3.066272) q[2];
sx q[2];
rz(-1.9131129) q[2];
sx q[2];
rz(-1.7750193) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5308456) q[1];
sx q[1];
rz(-1.3813003) q[1];
sx q[1];
rz(-0.23294542) q[1];
rz(2.1717668) q[3];
sx q[3];
rz(-1.9221913) q[3];
sx q[3];
rz(1.9072208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7829973) q[2];
sx q[2];
rz(-1.2025669) q[2];
sx q[2];
rz(0.45004582) q[2];
rz(-2.3027244) q[3];
sx q[3];
rz(-1.5476371) q[3];
sx q[3];
rz(-2.94256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49804509) q[0];
sx q[0];
rz(-1.4692551) q[0];
sx q[0];
rz(-3.041748) q[0];
rz(-1.3205344) q[1];
sx q[1];
rz(-2.1456199) q[1];
sx q[1];
rz(-0.60755306) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6927101) q[0];
sx q[0];
rz(-1.6537871) q[0];
sx q[0];
rz(-1.6458428) q[0];
rz(1.4133006) q[2];
sx q[2];
rz(-1.0309426) q[2];
sx q[2];
rz(-1.6423722) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.4974571) q[1];
sx q[1];
rz(-2.1408399) q[1];
sx q[1];
rz(1.8291147) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0436375) q[3];
sx q[3];
rz(-1.3134911) q[3];
sx q[3];
rz(2.3549781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9472092) q[2];
sx q[2];
rz(-2.3781229) q[2];
sx q[2];
rz(-0.53059951) q[2];
rz(0.44140205) q[3];
sx q[3];
rz(-2.1680021) q[3];
sx q[3];
rz(-0.94943625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7049103) q[0];
sx q[0];
rz(-1.7969776) q[0];
sx q[0];
rz(0.5740903) q[0];
rz(-2.2615945) q[1];
sx q[1];
rz(-2.0206385) q[1];
sx q[1];
rz(-1.7990187) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2821776) q[0];
sx q[0];
rz(-1.1033774) q[0];
sx q[0];
rz(0.22478454) q[0];
rz(-0.75127496) q[2];
sx q[2];
rz(-2.5715368) q[2];
sx q[2];
rz(-1.7715724) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8042829) q[1];
sx q[1];
rz(-1.334467) q[1];
sx q[1];
rz(-1.1996891) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.430949) q[3];
sx q[3];
rz(-2.2290843) q[3];
sx q[3];
rz(-2.0161395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4673956) q[2];
sx q[2];
rz(-2.8574222) q[2];
sx q[2];
rz(2.6944842) q[2];
rz(-1.0304662) q[3];
sx q[3];
rz(-1.6298529) q[3];
sx q[3];
rz(-0.38465056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36987385) q[0];
sx q[0];
rz(-1.8444703) q[0];
sx q[0];
rz(-0.41279992) q[0];
rz(2.0665456) q[1];
sx q[1];
rz(-2.0037035) q[1];
sx q[1];
rz(0.80361754) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5435073) q[0];
sx q[0];
rz(-2.0830562) q[0];
sx q[0];
rz(2.5582163) q[0];
rz(-0.86974025) q[2];
sx q[2];
rz(-2.365843) q[2];
sx q[2];
rz(-0.23245959) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0688096) q[1];
sx q[1];
rz(-2.3549665) q[1];
sx q[1];
rz(1.5907611) q[1];
rz(-pi) q[2];
rz(0.39866205) q[3];
sx q[3];
rz(-0.19426647) q[3];
sx q[3];
rz(-2.8630321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2794118) q[2];
sx q[2];
rz(-1.2181686) q[2];
sx q[2];
rz(-1.5009521) q[2];
rz(2.5189279) q[3];
sx q[3];
rz(-0.77087918) q[3];
sx q[3];
rz(-1.6525432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.72614661) q[0];
sx q[0];
rz(-1.6212689) q[0];
sx q[0];
rz(-0.55794445) q[0];
rz(0.56124148) q[1];
sx q[1];
rz(-1.3713501) q[1];
sx q[1];
rz(1.7254613) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87420207) q[0];
sx q[0];
rz(-2.6010397) q[0];
sx q[0];
rz(-0.99733277) q[0];
rz(-pi) q[1];
rz(2.8475262) q[2];
sx q[2];
rz(-1.9379788) q[2];
sx q[2];
rz(-2.4096556) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.76991612) q[1];
sx q[1];
rz(-1.8990108) q[1];
sx q[1];
rz(0.10336831) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3660356) q[3];
sx q[3];
rz(-2.3157218) q[3];
sx q[3];
rz(1.5457758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2531565) q[2];
sx q[2];
rz(-2.3641868) q[2];
sx q[2];
rz(1.0003132) q[2];
rz(0.12217626) q[3];
sx q[3];
rz(-1.6413942) q[3];
sx q[3];
rz(0.15672556) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67824739) q[0];
sx q[0];
rz(-1.954701) q[0];
sx q[0];
rz(0.004322411) q[0];
rz(-0.16383544) q[1];
sx q[1];
rz(-2.7163353) q[1];
sx q[1];
rz(1.3660627) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28986606) q[0];
sx q[0];
rz(-2.4938739) q[0];
sx q[0];
rz(-1.5246379) q[0];
rz(1.8810358) q[2];
sx q[2];
rz(-1.2149335) q[2];
sx q[2];
rz(2.9663393) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7228702) q[1];
sx q[1];
rz(-1.8218177) q[1];
sx q[1];
rz(-1.5334362) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.66994623) q[3];
sx q[3];
rz(-0.94196999) q[3];
sx q[3];
rz(-0.82312246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0020478) q[2];
sx q[2];
rz(-0.9674955) q[2];
sx q[2];
rz(0.61919332) q[2];
rz(-2.5891417) q[3];
sx q[3];
rz(-1.8360454) q[3];
sx q[3];
rz(-2.9964871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3510975) q[0];
sx q[0];
rz(-2.9572697) q[0];
sx q[0];
rz(-2.6628394) q[0];
rz(-1.152773) q[1];
sx q[1];
rz(-0.29007998) q[1];
sx q[1];
rz(0.94720381) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21175948) q[0];
sx q[0];
rz(-0.18778983) q[0];
sx q[0];
rz(-1.7540356) q[0];
x q[1];
rz(0.2453863) q[2];
sx q[2];
rz(-2.1446501) q[2];
sx q[2];
rz(1.4978486) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.041872) q[1];
sx q[1];
rz(-0.86468378) q[1];
sx q[1];
rz(0.40634917) q[1];
x q[2];
rz(-1.3218888) q[3];
sx q[3];
rz(-1.2281903) q[3];
sx q[3];
rz(-0.026503868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3843711) q[2];
sx q[2];
rz(-2.9360076) q[2];
sx q[2];
rz(0.26665404) q[2];
rz(1.0673149) q[3];
sx q[3];
rz(-1.9284733) q[3];
sx q[3];
rz(2.1783569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0530479) q[0];
sx q[0];
rz(-1.5374669) q[0];
sx q[0];
rz(1.5446825) q[0];
rz(1.0531986) q[1];
sx q[1];
rz(-1.9269301) q[1];
sx q[1];
rz(-2.3763837) q[1];
rz(1.5736035) q[2];
sx q[2];
rz(-1.8127828) q[2];
sx q[2];
rz(0.91290963) q[2];
rz(-0.33452928) q[3];
sx q[3];
rz(-0.60548895) q[3];
sx q[3];
rz(2.0169712) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
