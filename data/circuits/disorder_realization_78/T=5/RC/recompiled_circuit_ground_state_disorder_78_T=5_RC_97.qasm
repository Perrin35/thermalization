OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8371589) q[0];
sx q[0];
rz(4.9871939) q[0];
sx q[0];
rz(11.468588) q[0];
rz(2.1863565) q[1];
sx q[1];
rz(-0.74164852) q[1];
sx q[1];
rz(0.15129605) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79729743) q[0];
sx q[0];
rz(-1.5963285) q[0];
sx q[0];
rz(-1.6961369) q[0];
rz(-pi) q[1];
rz(1.382231) q[2];
sx q[2];
rz(-2.5841568) q[2];
sx q[2];
rz(-0.093890015) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9469493) q[1];
sx q[1];
rz(-0.8425172) q[1];
sx q[1];
rz(0.76232736) q[1];
rz(-pi) q[2];
rz(1.2829078) q[3];
sx q[3];
rz(-0.79018738) q[3];
sx q[3];
rz(2.0048646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0155045) q[2];
sx q[2];
rz(-1.2903004) q[2];
sx q[2];
rz(-2.8224831) q[2];
rz(-2.0305521) q[3];
sx q[3];
rz(-2.5824661) q[3];
sx q[3];
rz(1.3148974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.023271712) q[0];
sx q[0];
rz(-0.51902223) q[0];
sx q[0];
rz(2.5530489) q[0];
rz(-2.5449246) q[1];
sx q[1];
rz(-1.8110954) q[1];
sx q[1];
rz(-2.9002424) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7241192) q[0];
sx q[0];
rz(-2.3406174) q[0];
sx q[0];
rz(-2.1885314) q[0];
rz(-pi) q[1];
rz(1.4407519) q[2];
sx q[2];
rz(-2.6374014) q[2];
sx q[2];
rz(0.90182226) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7244563) q[1];
sx q[1];
rz(-1.1358741) q[1];
sx q[1];
rz(0.10970727) q[1];
rz(-pi) q[2];
rz(2.1961658) q[3];
sx q[3];
rz(-1.6540048) q[3];
sx q[3];
rz(2.0997467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0236464) q[2];
sx q[2];
rz(-1.4586552) q[2];
sx q[2];
rz(3.0493951) q[2];
rz(-2.2495031) q[3];
sx q[3];
rz(-2.2690319) q[3];
sx q[3];
rz(1.0649072) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43111619) q[0];
sx q[0];
rz(-1.5065864) q[0];
sx q[0];
rz(-0.098467501) q[0];
rz(0.59421986) q[1];
sx q[1];
rz(-1.0585982) q[1];
sx q[1];
rz(-0.99064151) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7567609) q[0];
sx q[0];
rz(-0.38479003) q[0];
sx q[0];
rz(-1.8033474) q[0];
rz(-2.5712396) q[2];
sx q[2];
rz(-2.2945171) q[2];
sx q[2];
rz(0.090426771) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2222683) q[1];
sx q[1];
rz(-1.8152555) q[1];
sx q[1];
rz(0.9217086) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1414755) q[3];
sx q[3];
rz(-0.50715441) q[3];
sx q[3];
rz(-1.0454659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.46382612) q[2];
sx q[2];
rz(-2.3550484) q[2];
sx q[2];
rz(2.3393935) q[2];
rz(1.5471316) q[3];
sx q[3];
rz(-2.0764669) q[3];
sx q[3];
rz(-2.2772363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7441854) q[0];
sx q[0];
rz(-0.35905251) q[0];
sx q[0];
rz(1.3209976) q[0];
rz(-1.6162704) q[1];
sx q[1];
rz(-2.1883712) q[1];
sx q[1];
rz(-2.9811409) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4904719) q[0];
sx q[0];
rz(-1.3297538) q[0];
sx q[0];
rz(0.81220497) q[0];
rz(-1.686626) q[2];
sx q[2];
rz(-1.7421075) q[2];
sx q[2];
rz(-1.1095604) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.39059475) q[1];
sx q[1];
rz(-0.92727755) q[1];
sx q[1];
rz(-2.949763) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9879278) q[3];
sx q[3];
rz(-1.0730181) q[3];
sx q[3];
rz(-2.1912632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.82115951) q[2];
sx q[2];
rz(-1.4904138) q[2];
sx q[2];
rz(-0.56905812) q[2];
rz(1.2218366) q[3];
sx q[3];
rz(-2.8025083) q[3];
sx q[3];
rz(-0.1327742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4959167) q[0];
sx q[0];
rz(-2.3585632) q[0];
sx q[0];
rz(2.3107279) q[0];
rz(-1.9310541) q[1];
sx q[1];
rz(-1.4930864) q[1];
sx q[1];
rz(2.1499706) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.338991) q[0];
sx q[0];
rz(-2.3590238) q[0];
sx q[0];
rz(-2.0156142) q[0];
rz(-pi) q[1];
rz(1.9087725) q[2];
sx q[2];
rz(-0.71546474) q[2];
sx q[2];
rz(1.6140661) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2455006) q[1];
sx q[1];
rz(-0.84907167) q[1];
sx q[1];
rz(-0.20770276) q[1];
rz(-pi) q[2];
rz(0.5414821) q[3];
sx q[3];
rz(-2.1953744) q[3];
sx q[3];
rz(1.9196212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7603989) q[2];
sx q[2];
rz(-2.2187967) q[2];
sx q[2];
rz(0.35923108) q[2];
rz(0.71581101) q[3];
sx q[3];
rz(-0.78407136) q[3];
sx q[3];
rz(1.1904967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0773709) q[0];
sx q[0];
rz(-1.8697898) q[0];
sx q[0];
rz(2.6203058) q[0];
rz(-0.84292665) q[1];
sx q[1];
rz(-1.9631674) q[1];
sx q[1];
rz(3.0311323) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4058286) q[0];
sx q[0];
rz(-0.89019201) q[0];
sx q[0];
rz(1.243478) q[0];
rz(-pi) q[1];
rz(-1.5311538) q[2];
sx q[2];
rz(-0.64293282) q[2];
sx q[2];
rz(0.55243353) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4723052) q[1];
sx q[1];
rz(-0.51176039) q[1];
sx q[1];
rz(1.157758) q[1];
rz(-pi) q[2];
rz(1.13917) q[3];
sx q[3];
rz(-0.56474287) q[3];
sx q[3];
rz(0.016591681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9342186) q[2];
sx q[2];
rz(-1.5997581) q[2];
sx q[2];
rz(-2.1477487) q[2];
rz(2.5908616) q[3];
sx q[3];
rz(-2.6271074) q[3];
sx q[3];
rz(1.6186835) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9024502) q[0];
sx q[0];
rz(-0.37343326) q[0];
sx q[0];
rz(2.9504839) q[0];
rz(-0.36901078) q[1];
sx q[1];
rz(-1.7419107) q[1];
sx q[1];
rz(-0.63327995) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12145081) q[0];
sx q[0];
rz(-2.883054) q[0];
sx q[0];
rz(-0.38991897) q[0];
rz(-pi) q[1];
rz(0.52187829) q[2];
sx q[2];
rz(-2.3503135) q[2];
sx q[2];
rz(-0.14497862) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.889588) q[1];
sx q[1];
rz(-1.7637296) q[1];
sx q[1];
rz(-1.7798406) q[1];
x q[2];
rz(-1.7353667) q[3];
sx q[3];
rz(-2.2630082) q[3];
sx q[3];
rz(0.084567955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9621027) q[2];
sx q[2];
rz(-1.7682163) q[2];
sx q[2];
rz(-0.38086677) q[2];
rz(-1.3880091) q[3];
sx q[3];
rz(-1.0264779) q[3];
sx q[3];
rz(1.4306205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4717344) q[0];
sx q[0];
rz(-2.7405881) q[0];
sx q[0];
rz(1.5420472) q[0];
rz(-1.764864) q[1];
sx q[1];
rz(-1.7574666) q[1];
sx q[1];
rz(-0.9309887) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5409398) q[0];
sx q[0];
rz(-2.4993745) q[0];
sx q[0];
rz(0.14677958) q[0];
rz(-pi) q[1];
rz(2.6561894) q[2];
sx q[2];
rz(-2.6377262) q[2];
sx q[2];
rz(0.47270838) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.763995) q[1];
sx q[1];
rz(-2.3990409) q[1];
sx q[1];
rz(1.5686036) q[1];
rz(2.2119207) q[3];
sx q[3];
rz(-2.3768209) q[3];
sx q[3];
rz(1.9509893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.57637438) q[2];
sx q[2];
rz(-0.56052506) q[2];
sx q[2];
rz(-0.068923846) q[2];
rz(-1.2636412) q[3];
sx q[3];
rz(-0.37216035) q[3];
sx q[3];
rz(2.4431958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-0.7428335) q[0];
sx q[0];
rz(-0.95785207) q[0];
sx q[0];
rz(0.051890705) q[0];
rz(-0.17008153) q[1];
sx q[1];
rz(-0.64337987) q[1];
sx q[1];
rz(2.7896519) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6860117) q[0];
sx q[0];
rz(-0.36893836) q[0];
sx q[0];
rz(0.12646778) q[0];
rz(-2.1743618) q[2];
sx q[2];
rz(-1.033342) q[2];
sx q[2];
rz(0.010802566) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.70507732) q[1];
sx q[1];
rz(-1.962858) q[1];
sx q[1];
rz(-2.1432671) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.073411302) q[3];
sx q[3];
rz(-1.2651099) q[3];
sx q[3];
rz(2.2524407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.45741442) q[2];
sx q[2];
rz(-0.64932051) q[2];
sx q[2];
rz(2.7049098) q[2];
rz(2.7501578) q[3];
sx q[3];
rz(-1.4623564) q[3];
sx q[3];
rz(-2.2552538) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52142414) q[0];
sx q[0];
rz(-2.4346209) q[0];
sx q[0];
rz(-3.022505) q[0];
rz(-1.8424312) q[1];
sx q[1];
rz(-1.7487339) q[1];
sx q[1];
rz(-1.7838759) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8004476) q[0];
sx q[0];
rz(-2.1954932) q[0];
sx q[0];
rz(-0.50826061) q[0];
rz(-pi) q[1];
rz(1.18612) q[2];
sx q[2];
rz(-0.78868491) q[2];
sx q[2];
rz(0.91463156) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3849302) q[1];
sx q[1];
rz(-2.1637056) q[1];
sx q[1];
rz(1.1641527) q[1];
rz(-pi) q[2];
rz(-2.1351027) q[3];
sx q[3];
rz(-2.3052633) q[3];
sx q[3];
rz(1.8759954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.208821) q[2];
sx q[2];
rz(-1.5893156) q[2];
sx q[2];
rz(0.61441747) q[2];
rz(1.5116073) q[3];
sx q[3];
rz(-0.67174086) q[3];
sx q[3];
rz(-3.0577799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-2.3096302) q[0];
sx q[0];
rz(-2.2086668) q[0];
sx q[0];
rz(-2.8902239) q[0];
rz(2.108719) q[1];
sx q[1];
rz(-1.947247) q[1];
sx q[1];
rz(-1.4364545) q[1];
rz(1.6757552) q[2];
sx q[2];
rz(-1.5009673) q[2];
sx q[2];
rz(-3.0143723) q[2];
rz(0.80794215) q[3];
sx q[3];
rz(-2.1389037) q[3];
sx q[3];
rz(1.1566333) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
