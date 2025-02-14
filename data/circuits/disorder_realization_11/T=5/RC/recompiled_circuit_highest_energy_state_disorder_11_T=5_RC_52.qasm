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
rz(-1.7186681) q[0];
sx q[0];
rz(-1.0942425) q[0];
sx q[0];
rz(-2.8835468) q[0];
rz(2.1482422) q[1];
sx q[1];
rz(1.8408096) q[1];
sx q[1];
rz(9.169133) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6171744) q[0];
sx q[0];
rz(-0.99992311) q[0];
sx q[0];
rz(0.85321315) q[0];
rz(2.9886888) q[2];
sx q[2];
rz(-2.7488378) q[2];
sx q[2];
rz(-0.037234779) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.62754831) q[1];
sx q[1];
rz(-1.6395634) q[1];
sx q[1];
rz(1.3876983) q[1];
rz(1.7840476) q[3];
sx q[3];
rz(-2.0590326) q[3];
sx q[3];
rz(0.31682107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6306182) q[2];
sx q[2];
rz(-1.3773842) q[2];
sx q[2];
rz(2.0261436) q[2];
rz(-3.1274146) q[3];
sx q[3];
rz(-1.3445798) q[3];
sx q[3];
rz(0.43629638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(1.2071335) q[0];
sx q[0];
rz(-0.62196982) q[0];
sx q[0];
rz(-1.2472664) q[0];
rz(2.6994052) q[1];
sx q[1];
rz(-1.3860044) q[1];
sx q[1];
rz(0.8173379) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7572875) q[0];
sx q[0];
rz(-1.7624859) q[0];
sx q[0];
rz(2.1038281) q[0];
rz(0.27702443) q[2];
sx q[2];
rz(-1.1454586) q[2];
sx q[2];
rz(2.9556731) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4358959) q[1];
sx q[1];
rz(-1.7361987) q[1];
sx q[1];
rz(0.64343217) q[1];
rz(-2.2226187) q[3];
sx q[3];
rz(-2.2187244) q[3];
sx q[3];
rz(1.0046665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3680129) q[2];
sx q[2];
rz(-2.7648338) q[2];
sx q[2];
rz(-2.9212941) q[2];
rz(-1.9593272) q[3];
sx q[3];
rz(-1.9672829) q[3];
sx q[3];
rz(-3.0784472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0910864) q[0];
sx q[0];
rz(-1.8244705) q[0];
sx q[0];
rz(-2.0029946) q[0];
rz(0.41172045) q[1];
sx q[1];
rz(-0.76985923) q[1];
sx q[1];
rz(2.128111) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.103055) q[0];
sx q[0];
rz(-2.8881209) q[0];
sx q[0];
rz(0.059132476) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0092952) q[2];
sx q[2];
rz(-2.0084642) q[2];
sx q[2];
rz(2.3486111) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.895993) q[1];
sx q[1];
rz(-2.2618544) q[1];
sx q[1];
rz(-2.5697903) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0115836) q[3];
sx q[3];
rz(-2.293236) q[3];
sx q[3];
rz(2.5353081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.98075214) q[2];
sx q[2];
rz(-1.9858805) q[2];
sx q[2];
rz(1.2551003) q[2];
rz(-2.8694425) q[3];
sx q[3];
rz(-2.3641219) q[3];
sx q[3];
rz(3.0009771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18144064) q[0];
sx q[0];
rz(-0.93638268) q[0];
sx q[0];
rz(-2.5312359) q[0];
rz(-0.70862526) q[1];
sx q[1];
rz(-2.1253864) q[1];
sx q[1];
rz(-1.5708539) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8023636) q[0];
sx q[0];
rz(-1.7238364) q[0];
sx q[0];
rz(1.7101076) q[0];
rz(-pi) q[1];
rz(-0.28684692) q[2];
sx q[2];
rz(-1.45668) q[2];
sx q[2];
rz(-0.6335887) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.849087) q[1];
sx q[1];
rz(-1.6189112) q[1];
sx q[1];
rz(2.8699257) q[1];
rz(-2.7486984) q[3];
sx q[3];
rz(-0.84019444) q[3];
sx q[3];
rz(3.068416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9307956) q[2];
sx q[2];
rz(-1.0127298) q[2];
sx q[2];
rz(-2.5861758) q[2];
rz(-0.72426116) q[3];
sx q[3];
rz(-1.1490425) q[3];
sx q[3];
rz(-2.485062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.637218) q[0];
sx q[0];
rz(-1.1906304) q[0];
sx q[0];
rz(-0.36886886) q[0];
rz(-0.24636191) q[1];
sx q[1];
rz(-1.3280222) q[1];
sx q[1];
rz(1.7074283) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95043889) q[0];
sx q[0];
rz(-0.75685793) q[0];
sx q[0];
rz(-0.0075154742) q[0];
x q[1];
rz(-2.6961961) q[2];
sx q[2];
rz(-1.4120308) q[2];
sx q[2];
rz(2.7743055) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.37232698) q[1];
sx q[1];
rz(-0.40180909) q[1];
sx q[1];
rz(0.010784464) q[1];
rz(1.1617817) q[3];
sx q[3];
rz(-2.8944765) q[3];
sx q[3];
rz(-2.1638526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4904334) q[2];
sx q[2];
rz(-1.8305402) q[2];
sx q[2];
rz(1.0820214) q[2];
rz(-0.79484445) q[3];
sx q[3];
rz(-1.669603) q[3];
sx q[3];
rz(-1.2580416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.865888) q[0];
sx q[0];
rz(-0.10144932) q[0];
sx q[0];
rz(-0.24359447) q[0];
rz(0.98681915) q[1];
sx q[1];
rz(-2.3934264) q[1];
sx q[1];
rz(-0.93596828) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83171885) q[0];
sx q[0];
rz(-1.9666934) q[0];
sx q[0];
rz(-2.7926366) q[0];
x q[1];
rz(1.1963821) q[2];
sx q[2];
rz(-0.78305675) q[2];
sx q[2];
rz(-0.20314344) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0985581) q[1];
sx q[1];
rz(-1.1698616) q[1];
sx q[1];
rz(0.75717302) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1706287) q[3];
sx q[3];
rz(-1.5884627) q[3];
sx q[3];
rz(-2.3223557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.662107) q[2];
sx q[2];
rz(-1.138849) q[2];
sx q[2];
rz(-0.34995079) q[2];
rz(-0.55772603) q[3];
sx q[3];
rz(-1.2099268) q[3];
sx q[3];
rz(0.85132712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-1.6996985) q[0];
sx q[0];
rz(-0.63755578) q[0];
sx q[0];
rz(0.30174524) q[0];
rz(2.8217577) q[1];
sx q[1];
rz(-1.539307) q[1];
sx q[1];
rz(-1.1357657) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2280884) q[0];
sx q[0];
rz(-0.64675179) q[0];
sx q[0];
rz(0.089368377) q[0];
rz(-2.2280424) q[2];
sx q[2];
rz(-0.74591178) q[2];
sx q[2];
rz(2.8532956) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3150683) q[1];
sx q[1];
rz(-0.60501912) q[1];
sx q[1];
rz(1.4304763) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0301827) q[3];
sx q[3];
rz(-0.79395959) q[3];
sx q[3];
rz(-2.4536228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4096421) q[2];
sx q[2];
rz(-2.4739517) q[2];
sx q[2];
rz(0.14275924) q[2];
rz(-1.3879294) q[3];
sx q[3];
rz(-1.9526491) q[3];
sx q[3];
rz(-0.68283844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1801572) q[0];
sx q[0];
rz(-0.016594369) q[0];
sx q[0];
rz(-0.55602443) q[0];
rz(-2.9122638) q[1];
sx q[1];
rz(-1.2622958) q[1];
sx q[1];
rz(-1.75846) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97998226) q[0];
sx q[0];
rz(-1.5357657) q[0];
sx q[0];
rz(-0.73765124) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6298619) q[2];
sx q[2];
rz(-2.0286273) q[2];
sx q[2];
rz(-2.3659467) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.759093) q[1];
sx q[1];
rz(-1.9803932) q[1];
sx q[1];
rz(0.31633693) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1874299) q[3];
sx q[3];
rz(-1.1463506) q[3];
sx q[3];
rz(0.92680537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6250299) q[2];
sx q[2];
rz(-0.27250686) q[2];
sx q[2];
rz(1.9005091) q[2];
rz(-0.062601335) q[3];
sx q[3];
rz(-1.4227941) q[3];
sx q[3];
rz(-0.82714287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1271707) q[0];
sx q[0];
rz(-1.4143455) q[0];
sx q[0];
rz(-2.993809) q[0];
rz(2.6328909) q[1];
sx q[1];
rz(-1.769442) q[1];
sx q[1];
rz(2.7630189) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5145549) q[0];
sx q[0];
rz(-0.48279027) q[0];
sx q[0];
rz(-0.52282368) q[0];
rz(2.0835593) q[2];
sx q[2];
rz(-2.4840925) q[2];
sx q[2];
rz(0.12557827) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.3264824) q[1];
sx q[1];
rz(-1.2604144) q[1];
sx q[1];
rz(1.2182477) q[1];
rz(-pi) q[2];
rz(0.13033615) q[3];
sx q[3];
rz(-2.3763083) q[3];
sx q[3];
rz(1.9627067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.96616894) q[2];
sx q[2];
rz(-0.95169008) q[2];
sx q[2];
rz(1.2790722) q[2];
rz(1.7133948) q[3];
sx q[3];
rz(-1.0019852) q[3];
sx q[3];
rz(-2.9960347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3928423) q[0];
sx q[0];
rz(-0.57806438) q[0];
sx q[0];
rz(-3.0774935) q[0];
rz(2.6129258) q[1];
sx q[1];
rz(-1.8730947) q[1];
sx q[1];
rz(0.97506964) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.753016) q[0];
sx q[0];
rz(-2.574769) q[0];
sx q[0];
rz(0.50811572) q[0];
rz(-pi) q[1];
x q[1];
rz(0.56985241) q[2];
sx q[2];
rz(-2.0790711) q[2];
sx q[2];
rz(-0.13392553) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.84555039) q[1];
sx q[1];
rz(-0.49440372) q[1];
sx q[1];
rz(2.666996) q[1];
rz(0.24419489) q[3];
sx q[3];
rz(-1.7932442) q[3];
sx q[3];
rz(-0.73058587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9958682) q[2];
sx q[2];
rz(-0.66103649) q[2];
sx q[2];
rz(-2.5346942) q[2];
rz(-2.8912344) q[3];
sx q[3];
rz(-1.7253877) q[3];
sx q[3];
rz(-0.50601602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.09457) q[0];
sx q[0];
rz(-1.9245514) q[0];
sx q[0];
rz(-1.2522329) q[0];
rz(0.49867123) q[1];
sx q[1];
rz(-1.55232) q[1];
sx q[1];
rz(1.3943863) q[1];
rz(-2.5468536) q[2];
sx q[2];
rz(-0.95437106) q[2];
sx q[2];
rz(0.30827733) q[2];
rz(1.3648894) q[3];
sx q[3];
rz(-1.4510703) q[3];
sx q[3];
rz(2.5616796) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
