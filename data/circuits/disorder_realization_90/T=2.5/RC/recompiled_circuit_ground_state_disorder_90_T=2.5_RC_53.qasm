OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5034135) q[0];
sx q[0];
rz(-1.2794275) q[0];
sx q[0];
rz(0.77603618) q[0];
rz(-2.4322721) q[1];
sx q[1];
rz(-1.5172989) q[1];
sx q[1];
rz(-2.5425743) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31604345) q[0];
sx q[0];
rz(-0.49015309) q[0];
sx q[0];
rz(1.8328299) q[0];
rz(-pi) q[1];
rz(-0.28223306) q[2];
sx q[2];
rz(-0.90105173) q[2];
sx q[2];
rz(1.9053659) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4645965) q[1];
sx q[1];
rz(-1.076816) q[1];
sx q[1];
rz(-1.1679653) q[1];
x q[2];
rz(2.0008068) q[3];
sx q[3];
rz(-0.83807349) q[3];
sx q[3];
rz(-0.091146745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1385931) q[2];
sx q[2];
rz(-2.004576) q[2];
sx q[2];
rz(-2.9453759) q[2];
rz(-2.0972706) q[3];
sx q[3];
rz(-2.315867) q[3];
sx q[3];
rz(0.31061068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1021378) q[0];
sx q[0];
rz(-2.4882443) q[0];
sx q[0];
rz(0.85773221) q[0];
rz(-0.54025447) q[1];
sx q[1];
rz(-2.0876355) q[1];
sx q[1];
rz(-0.78261715) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7359813) q[0];
sx q[0];
rz(-2.1850039) q[0];
sx q[0];
rz(0.10615291) q[0];
rz(-0.94901325) q[2];
sx q[2];
rz(-1.7634541) q[2];
sx q[2];
rz(0.72177835) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2007717) q[1];
sx q[1];
rz(-1.7634974) q[1];
sx q[1];
rz(-0.60639834) q[1];
x q[2];
rz(2.5839543) q[3];
sx q[3];
rz(-2.1774315) q[3];
sx q[3];
rz(0.56822694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9767849) q[2];
sx q[2];
rz(-0.79443496) q[2];
sx q[2];
rz(0.59107333) q[2];
rz(1.0278541) q[3];
sx q[3];
rz(-0.63575345) q[3];
sx q[3];
rz(3.0052321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9117821) q[0];
sx q[0];
rz(-2.6214143) q[0];
sx q[0];
rz(-1.0562563) q[0];
rz(0.99110574) q[1];
sx q[1];
rz(-1.7678363) q[1];
sx q[1];
rz(-2.3371005) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0206385) q[0];
sx q[0];
rz(-2.197663) q[0];
sx q[0];
rz(-1.9009717) q[0];
rz(-0.94203888) q[2];
sx q[2];
rz(-1.4938746) q[2];
sx q[2];
rz(-0.49401894) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5544598) q[1];
sx q[1];
rz(-1.3398223) q[1];
sx q[1];
rz(1.323631) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3785754) q[3];
sx q[3];
rz(-0.9323191) q[3];
sx q[3];
rz(2.9137127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.450401) q[2];
sx q[2];
rz(-2.2914026) q[2];
sx q[2];
rz(0.56785339) q[2];
rz(-1.19207) q[3];
sx q[3];
rz(-2.6046533) q[3];
sx q[3];
rz(2.9862064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.49517) q[0];
sx q[0];
rz(-1.7429054) q[0];
sx q[0];
rz(1.1871185) q[0];
rz(2.8861956) q[1];
sx q[1];
rz(-1.7385769) q[1];
sx q[1];
rz(-0.38937169) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5449937) q[0];
sx q[0];
rz(-1.7905718) q[0];
sx q[0];
rz(1.4334428) q[0];
rz(-1.6045447) q[2];
sx q[2];
rz(-2.4141888) q[2];
sx q[2];
rz(-0.45798618) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.90044777) q[1];
sx q[1];
rz(-0.68640781) q[1];
sx q[1];
rz(-1.871984) q[1];
rz(-pi) q[2];
rz(0.13279543) q[3];
sx q[3];
rz(-1.7657874) q[3];
sx q[3];
rz(-2.8094069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4046459) q[2];
sx q[2];
rz(-0.32361042) q[2];
sx q[2];
rz(1.7741989) q[2];
rz(-2.6020452) q[3];
sx q[3];
rz(-1.5522141) q[3];
sx q[3];
rz(-1.3246983) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0747727) q[0];
sx q[0];
rz(-2.9434151) q[0];
sx q[0];
rz(0.5603801) q[0];
rz(-2.9226774) q[1];
sx q[1];
rz(-2.0603265) q[1];
sx q[1];
rz(2.970649) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31256277) q[0];
sx q[0];
rz(-1.1033022) q[0];
sx q[0];
rz(-0.74689052) q[0];
rz(-pi) q[1];
rz(0.20410164) q[2];
sx q[2];
rz(-1.7095672) q[2];
sx q[2];
rz(2.6134174) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7485908) q[1];
sx q[1];
rz(-2.1868863) q[1];
sx q[1];
rz(-2.2115117) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6690977) q[3];
sx q[3];
rz(-0.62255961) q[3];
sx q[3];
rz(0.05427256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.35974744) q[2];
sx q[2];
rz(-2.013194) q[2];
sx q[2];
rz(-0.95853364) q[2];
rz(-2.9790699) q[3];
sx q[3];
rz(-2.2992117) q[3];
sx q[3];
rz(-1.5925647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5684587) q[0];
sx q[0];
rz(-0.064366654) q[0];
sx q[0];
rz(-0.74813133) q[0];
rz(1.1185147) q[1];
sx q[1];
rz(-0.41901127) q[1];
sx q[1];
rz(-2.8245139) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1823163) q[0];
sx q[0];
rz(-1.3249363) q[0];
sx q[0];
rz(1.5005174) q[0];
x q[1];
rz(0.24131624) q[2];
sx q[2];
rz(-0.20212999) q[2];
sx q[2];
rz(1.2723107) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.78817716) q[1];
sx q[1];
rz(-1.4995575) q[1];
sx q[1];
rz(-2.337238) q[1];
rz(-pi) q[2];
rz(-0.040359453) q[3];
sx q[3];
rz(-2.3642614) q[3];
sx q[3];
rz(-3.0632927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3743484) q[2];
sx q[2];
rz(-2.2129462) q[2];
sx q[2];
rz(1.7283776) q[2];
rz(0.080538571) q[3];
sx q[3];
rz(-1.7929411) q[3];
sx q[3];
rz(2.9197689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0291075) q[0];
sx q[0];
rz(-1.8928098) q[0];
sx q[0];
rz(-1.2741733) q[0];
rz(-1.796465) q[1];
sx q[1];
rz(-0.74256623) q[1];
sx q[1];
rz(1.8003731) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31650922) q[0];
sx q[0];
rz(-1.0183471) q[0];
sx q[0];
rz(-0.98487206) q[0];
rz(-2.6599081) q[2];
sx q[2];
rz(-1.8284214) q[2];
sx q[2];
rz(-3.1236908) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0772895) q[1];
sx q[1];
rz(-1.3796564) q[1];
sx q[1];
rz(2.9851476) q[1];
rz(-0.92067155) q[3];
sx q[3];
rz(-1.4772282) q[3];
sx q[3];
rz(2.2232995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.93366569) q[2];
sx q[2];
rz(-1.9542481) q[2];
sx q[2];
rz(2.6742317) q[2];
rz(0.76006877) q[3];
sx q[3];
rz(-1.6642539) q[3];
sx q[3];
rz(-1.2490341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6811328) q[0];
sx q[0];
rz(-1.0204027) q[0];
sx q[0];
rz(-1.6946174) q[0];
rz(2.815822) q[1];
sx q[1];
rz(-0.21369801) q[1];
sx q[1];
rz(1.6206585) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6353519) q[0];
sx q[0];
rz(-1.8367447) q[0];
sx q[0];
rz(-0.5269993) q[0];
rz(-pi) q[1];
rz(-0.09163945) q[2];
sx q[2];
rz(-1.5711725) q[2];
sx q[2];
rz(-2.3918652) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7829166) q[1];
sx q[1];
rz(-2.7442928) q[1];
sx q[1];
rz(-1.1043758) q[1];
rz(-pi) q[2];
x q[2];
rz(0.86674848) q[3];
sx q[3];
rz(-1.4102077) q[3];
sx q[3];
rz(2.989547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1763566) q[2];
sx q[2];
rz(-0.74702817) q[2];
sx q[2];
rz(0.61526862) q[2];
rz(1.2728914) q[3];
sx q[3];
rz(-0.26510173) q[3];
sx q[3];
rz(2.7191775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1212921) q[0];
sx q[0];
rz(-1.2739807) q[0];
sx q[0];
rz(2.2913388) q[0];
rz(-2.0203159) q[1];
sx q[1];
rz(-1.774615) q[1];
sx q[1];
rz(-2.3698295) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92895011) q[0];
sx q[0];
rz(-1.54689) q[0];
sx q[0];
rz(-0.8769518) q[0];
rz(-1.5007559) q[2];
sx q[2];
rz(-1.7511427) q[2];
sx q[2];
rz(-1.548962) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.38852019) q[1];
sx q[1];
rz(-2.8917312) q[1];
sx q[1];
rz(1.7953963) q[1];
rz(-pi) q[2];
rz(2.8706067) q[3];
sx q[3];
rz(-1.8838722) q[3];
sx q[3];
rz(-2.3199758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0261592) q[2];
sx q[2];
rz(-1.5449646) q[2];
sx q[2];
rz(2.3010632) q[2];
rz(-3.0606411) q[3];
sx q[3];
rz(-0.5778802) q[3];
sx q[3];
rz(-0.24278434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42846546) q[0];
sx q[0];
rz(-0.19793887) q[0];
sx q[0];
rz(-1.1567098) q[0];
rz(1.0276065) q[1];
sx q[1];
rz(-1.3778069) q[1];
sx q[1];
rz(2.3311232) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76475984) q[0];
sx q[0];
rz(-0.78928052) q[0];
sx q[0];
rz(-2.8349702) q[0];
rz(-1.2357622) q[2];
sx q[2];
rz(-1.050569) q[2];
sx q[2];
rz(-0.65641415) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0009044) q[1];
sx q[1];
rz(-1.2869121) q[1];
sx q[1];
rz(1.264601) q[1];
rz(-pi) q[2];
rz(1.4224206) q[3];
sx q[3];
rz(-2.1273489) q[3];
sx q[3];
rz(0.30239964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.345574) q[2];
sx q[2];
rz(-1.8859665) q[2];
sx q[2];
rz(1.6176809) q[2];
rz(-2.0412622) q[3];
sx q[3];
rz(-1.9610145) q[3];
sx q[3];
rz(-2.9617917) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1179467) q[0];
sx q[0];
rz(-1.4143586) q[0];
sx q[0];
rz(2.216862) q[0];
rz(-0.042451518) q[1];
sx q[1];
rz(-2.0274542) q[1];
sx q[1];
rz(-1.7995119) q[1];
rz(1.9770728) q[2];
sx q[2];
rz(-1.6222519) q[2];
sx q[2];
rz(-1.9380515) q[2];
rz(2.6841738) q[3];
sx q[3];
rz(-1.538496) q[3];
sx q[3];
rz(-2.9067007) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
