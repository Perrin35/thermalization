OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6239983) q[0];
sx q[0];
rz(-0.52596337) q[0];
sx q[0];
rz(0.21831231) q[0];
rz(-1.6649618) q[1];
sx q[1];
rz(-0.47774878) q[1];
sx q[1];
rz(0.55396095) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57620063) q[0];
sx q[0];
rz(-1.0668653) q[0];
sx q[0];
rz(0.53458696) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6126819) q[2];
sx q[2];
rz(-1.6220835) q[2];
sx q[2];
rz(-1.9330658) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2077929) q[1];
sx q[1];
rz(-1.7805532) q[1];
sx q[1];
rz(-2.2906474) q[1];
rz(-2.4943236) q[3];
sx q[3];
rz(-2.4993901) q[3];
sx q[3];
rz(2.6489779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.37609425) q[2];
sx q[2];
rz(-0.29355294) q[2];
sx q[2];
rz(1.2676839) q[2];
rz(2.3119161) q[3];
sx q[3];
rz(-1.5209578) q[3];
sx q[3];
rz(2.7177496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.053452881) q[0];
sx q[0];
rz(-1.3351048) q[0];
sx q[0];
rz(-1.7470737) q[0];
rz(0.8949737) q[1];
sx q[1];
rz(-1.0286237) q[1];
sx q[1];
rz(-1.5825533) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6000309) q[0];
sx q[0];
rz(-0.91225925) q[0];
sx q[0];
rz(0.70777871) q[0];
rz(-pi) q[1];
rz(-1.1505914) q[2];
sx q[2];
rz(-1.24345) q[2];
sx q[2];
rz(-0.98132747) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.30732111) q[1];
sx q[1];
rz(-2.2311192) q[1];
sx q[1];
rz(-1.1701533) q[1];
x q[2];
rz(-0.96554324) q[3];
sx q[3];
rz(-2.6355834) q[3];
sx q[3];
rz(-3.0288896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.47722694) q[2];
sx q[2];
rz(-1.5218647) q[2];
sx q[2];
rz(3.1108372) q[2];
rz(-2.6886046) q[3];
sx q[3];
rz(-0.22946295) q[3];
sx q[3];
rz(2.9636813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40193108) q[0];
sx q[0];
rz(-0.10098305) q[0];
sx q[0];
rz(0.82823753) q[0];
rz(-0.047686934) q[1];
sx q[1];
rz(-0.86367718) q[1];
sx q[1];
rz(-1.9140859) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21624529) q[0];
sx q[0];
rz(-1.1309393) q[0];
sx q[0];
rz(0.60103215) q[0];
rz(-pi) q[1];
rz(0.21593185) q[2];
sx q[2];
rz(-1.1732444) q[2];
sx q[2];
rz(2.2550607) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8125898) q[1];
sx q[1];
rz(-1.131445) q[1];
sx q[1];
rz(-2.8529608) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.84215409) q[3];
sx q[3];
rz(-1.6098566) q[3];
sx q[3];
rz(-2.0741163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0487655) q[2];
sx q[2];
rz(-0.91492492) q[2];
sx q[2];
rz(2.0406593) q[2];
rz(3.1277505) q[3];
sx q[3];
rz(-1.775454) q[3];
sx q[3];
rz(-2.3069416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58532995) q[0];
sx q[0];
rz(-1.4034554) q[0];
sx q[0];
rz(0.334326) q[0];
rz(-0.73257929) q[1];
sx q[1];
rz(-0.94894797) q[1];
sx q[1];
rz(-3.0308731) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3601013) q[0];
sx q[0];
rz(-1.555614) q[0];
sx q[0];
rz(-0.76155565) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.16558318) q[2];
sx q[2];
rz(-1.6628169) q[2];
sx q[2];
rz(1.657287) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9423805) q[1];
sx q[1];
rz(-0.95006493) q[1];
sx q[1];
rz(-0.86356461) q[1];
rz(-2.3663051) q[3];
sx q[3];
rz(-1.2167769) q[3];
sx q[3];
rz(0.86590761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.683814) q[2];
sx q[2];
rz(-2.8595698) q[2];
sx q[2];
rz(1.4078183) q[2];
rz(-1.1553361) q[3];
sx q[3];
rz(-1.9852394) q[3];
sx q[3];
rz(2.9467764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69283501) q[0];
sx q[0];
rz(-0.91788569) q[0];
sx q[0];
rz(0.99739972) q[0];
rz(1.5265441) q[1];
sx q[1];
rz(-0.63957447) q[1];
sx q[1];
rz(1.7220727) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2888694) q[0];
sx q[0];
rz(-2.8670681) q[0];
sx q[0];
rz(-0.92936607) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.15879719) q[2];
sx q[2];
rz(-2.3613075) q[2];
sx q[2];
rz(2.4312388) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.75052634) q[1];
sx q[1];
rz(-1.8896034) q[1];
sx q[1];
rz(0.90853779) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0180198) q[3];
sx q[3];
rz(-1.6032748) q[3];
sx q[3];
rz(-0.44636727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2898966) q[2];
sx q[2];
rz(-3.0471314) q[2];
sx q[2];
rz(-2.8186901) q[2];
rz(-1.1139392) q[3];
sx q[3];
rz(-2.0506004) q[3];
sx q[3];
rz(-0.41306257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7776529) q[0];
sx q[0];
rz(-2.8033065) q[0];
sx q[0];
rz(0.80663484) q[0];
rz(-0.58397645) q[1];
sx q[1];
rz(-1.1176502) q[1];
sx q[1];
rz(-1.1899828) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71680741) q[0];
sx q[0];
rz(-1.9182701) q[0];
sx q[0];
rz(0.13089546) q[0];
rz(-pi) q[1];
rz(-1.3596542) q[2];
sx q[2];
rz(-2.0953155) q[2];
sx q[2];
rz(-0.013583029) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2014043) q[1];
sx q[1];
rz(-1.7394627) q[1];
sx q[1];
rz(-0.052171295) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9526918) q[3];
sx q[3];
rz(-2.402225) q[3];
sx q[3];
rz(-0.96773558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.40818647) q[2];
sx q[2];
rz(-1.5032282) q[2];
sx q[2];
rz(-2.9564986) q[2];
rz(1.6144276) q[3];
sx q[3];
rz(-1.7388758) q[3];
sx q[3];
rz(0.68814284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9516893) q[0];
sx q[0];
rz(-2.1898495) q[0];
sx q[0];
rz(2.9290747) q[0];
rz(0.7849794) q[1];
sx q[1];
rz(-1.6467983) q[1];
sx q[1];
rz(-2.3627538) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76701421) q[0];
sx q[0];
rz(-2.2901669) q[0];
sx q[0];
rz(-0.79921754) q[0];
x q[1];
rz(-0.58765192) q[2];
sx q[2];
rz(-1.2562804) q[2];
sx q[2];
rz(-1.8314198) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.067233868) q[1];
sx q[1];
rz(-1.7869084) q[1];
sx q[1];
rz(1.8317779) q[1];
x q[2];
rz(0.4850579) q[3];
sx q[3];
rz(-2.2211255) q[3];
sx q[3];
rz(-1.8431839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.51398858) q[2];
sx q[2];
rz(-0.49391654) q[2];
sx q[2];
rz(2.7331875) q[2];
rz(0.70185316) q[3];
sx q[3];
rz(-2.2097094) q[3];
sx q[3];
rz(1.2592038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34847611) q[0];
sx q[0];
rz(-1.9261253) q[0];
sx q[0];
rz(3.075573) q[0];
rz(-1.6325715) q[1];
sx q[1];
rz(-2.0965818) q[1];
sx q[1];
rz(-0.95796934) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45487472) q[0];
sx q[0];
rz(-1.8544079) q[0];
sx q[0];
rz(-0.4768178) q[0];
rz(-pi) q[1];
rz(2.3249469) q[2];
sx q[2];
rz(-1.1743675) q[2];
sx q[2];
rz(-0.81260896) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0762679) q[1];
sx q[1];
rz(-0.75265127) q[1];
sx q[1];
rz(-2.6757043) q[1];
rz(-1.7244206) q[3];
sx q[3];
rz(-1.6699176) q[3];
sx q[3];
rz(-0.88872582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8336739) q[2];
sx q[2];
rz(-1.5235528) q[2];
sx q[2];
rz(1.7306805) q[2];
rz(-2.446567) q[3];
sx q[3];
rz(-1.4796939) q[3];
sx q[3];
rz(-0.01785774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0787635) q[0];
sx q[0];
rz(-1.48209) q[0];
sx q[0];
rz(-2.1516946) q[0];
rz(0.46317378) q[1];
sx q[1];
rz(-1.4521234) q[1];
sx q[1];
rz(-1.90082) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1027066) q[0];
sx q[0];
rz(-0.30065824) q[0];
sx q[0];
rz(1.7666398) q[0];
x q[1];
rz(-1.979901) q[2];
sx q[2];
rz(-2.6648403) q[2];
sx q[2];
rz(0.40845218) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1009192) q[1];
sx q[1];
rz(-1.8534436) q[1];
sx q[1];
rz(-1.8428749) q[1];
rz(0.97384805) q[3];
sx q[3];
rz(-2.429109) q[3];
sx q[3];
rz(1.2013916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8076294) q[2];
sx q[2];
rz(-1.7619851) q[2];
sx q[2];
rz(0.71869746) q[2];
rz(-1.1527609) q[3];
sx q[3];
rz(-1.4216239) q[3];
sx q[3];
rz(-2.1271472) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5929247) q[0];
sx q[0];
rz(-2.7935226) q[0];
sx q[0];
rz(2.3396709) q[0];
rz(1.0879263) q[1];
sx q[1];
rz(-1.8049003) q[1];
sx q[1];
rz(0.71802872) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3182326) q[0];
sx q[0];
rz(-1.8917068) q[0];
sx q[0];
rz(-2.4145187) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3921591) q[2];
sx q[2];
rz(-0.66808703) q[2];
sx q[2];
rz(-1.4632478) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1433318) q[1];
sx q[1];
rz(-1.5185396) q[1];
sx q[1];
rz(1.2388171) q[1];
rz(-pi) q[2];
rz(2.3613648) q[3];
sx q[3];
rz(-0.89892584) q[3];
sx q[3];
rz(1.5742009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.7903018) q[2];
sx q[2];
rz(-2.9238034) q[2];
sx q[2];
rz(-0.51266074) q[2];
rz(0.74448186) q[3];
sx q[3];
rz(-2.1148041) q[3];
sx q[3];
rz(0.7640394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22513334) q[0];
sx q[0];
rz(-1.8712578) q[0];
sx q[0];
rz(1.901392) q[0];
rz(-0.71612877) q[1];
sx q[1];
rz(-0.54812535) q[1];
sx q[1];
rz(-2.5352238) q[1];
rz(3.0307583) q[2];
sx q[2];
rz(-1.3468942) q[2];
sx q[2];
rz(-0.46769618) q[2];
rz(-1.2848787) q[3];
sx q[3];
rz(-2.7919522) q[3];
sx q[3];
rz(-1.5408564) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
