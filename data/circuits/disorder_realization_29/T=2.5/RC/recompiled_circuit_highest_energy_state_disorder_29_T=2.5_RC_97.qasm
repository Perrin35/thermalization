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
rz(1.2443378) q[0];
sx q[0];
rz(-0.79255784) q[0];
sx q[0];
rz(-0.25294024) q[0];
rz(-0.77415544) q[1];
sx q[1];
rz(-0.3781265) q[1];
sx q[1];
rz(0.3869431) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16432504) q[0];
sx q[0];
rz(-1.8935793) q[0];
sx q[0];
rz(0.27077814) q[0];
rz(-pi) q[1];
rz(2.5637881) q[2];
sx q[2];
rz(-1.4981604) q[2];
sx q[2];
rz(0.60608038) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7751037) q[1];
sx q[1];
rz(-1.4867884) q[1];
sx q[1];
rz(1.6506399) q[1];
rz(-2.4693842) q[3];
sx q[3];
rz(-2.432193) q[3];
sx q[3];
rz(-1.0649452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.046254961) q[2];
sx q[2];
rz(-0.75901186) q[2];
sx q[2];
rz(1.7389899) q[2];
rz(-1.4143573) q[3];
sx q[3];
rz(-2.4284913) q[3];
sx q[3];
rz(-0.49899092) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4614748) q[0];
sx q[0];
rz(-2.6533227) q[0];
sx q[0];
rz(-0.71037355) q[0];
rz(1.0164227) q[1];
sx q[1];
rz(-1.8417532) q[1];
sx q[1];
rz(-0.76510915) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.13694) q[0];
sx q[0];
rz(-2.1015133) q[0];
sx q[0];
rz(1.1999446) q[0];
x q[1];
rz(-2.3813407) q[2];
sx q[2];
rz(-2.2905802) q[2];
sx q[2];
rz(0.30530294) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4087832) q[1];
sx q[1];
rz(-1.8955232) q[1];
sx q[1];
rz(2.2546683) q[1];
rz(-pi) q[2];
x q[2];
rz(0.14987544) q[3];
sx q[3];
rz(-1.9016163) q[3];
sx q[3];
rz(-2.5811623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0244828) q[2];
sx q[2];
rz(-2.4958002) q[2];
sx q[2];
rz(2.4369241) q[2];
rz(-2.9674528) q[3];
sx q[3];
rz(-0.87248674) q[3];
sx q[3];
rz(-2.7599938) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0524549) q[0];
sx q[0];
rz(-1.8280886) q[0];
sx q[0];
rz(0.88743368) q[0];
rz(1.8916091) q[1];
sx q[1];
rz(-1.2108112) q[1];
sx q[1];
rz(-1.3608305) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0075390752) q[0];
sx q[0];
rz(-1.0666532) q[0];
sx q[0];
rz(-2.0043122) q[0];
rz(-pi) q[1];
x q[1];
rz(1.046692) q[2];
sx q[2];
rz(-2.9827318) q[2];
sx q[2];
rz(2.9422974) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2584784) q[1];
sx q[1];
rz(-0.7683903) q[1];
sx q[1];
rz(-1.7151296) q[1];
x q[2];
rz(-3.0988337) q[3];
sx q[3];
rz(-1.1819289) q[3];
sx q[3];
rz(2.1270909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.36797324) q[2];
sx q[2];
rz(-2.7757288) q[2];
sx q[2];
rz(1.6486637) q[2];
rz(0.44935539) q[3];
sx q[3];
rz(-1.2949233) q[3];
sx q[3];
rz(-1.9213283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0538977) q[0];
sx q[0];
rz(-2.5515285) q[0];
sx q[0];
rz(1.7328523) q[0];
rz(2.159481) q[1];
sx q[1];
rz(-1.5487919) q[1];
sx q[1];
rz(1.4062448) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81534751) q[0];
sx q[0];
rz(-0.24213386) q[0];
sx q[0];
rz(-1.0905488) q[0];
x q[1];
rz(-2.2135229) q[2];
sx q[2];
rz(-3.0756223) q[2];
sx q[2];
rz(0.82974354) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.92038735) q[1];
sx q[1];
rz(-1.8169329) q[1];
sx q[1];
rz(2.1046776) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2194524) q[3];
sx q[3];
rz(-1.3548082) q[3];
sx q[3];
rz(-2.7060425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1534319) q[2];
sx q[2];
rz(-1.016523) q[2];
sx q[2];
rz(-0.15596685) q[2];
rz(-2.4912452) q[3];
sx q[3];
rz(-1.9675083) q[3];
sx q[3];
rz(-0.11421886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0735737) q[0];
sx q[0];
rz(-1.2696215) q[0];
sx q[0];
rz(-2.9408348) q[0];
rz(-1.1772032) q[1];
sx q[1];
rz(-0.8526082) q[1];
sx q[1];
rz(-1.9815365) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8070712) q[0];
sx q[0];
rz(-1.40309) q[0];
sx q[0];
rz(1.7981862) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5963836) q[2];
sx q[2];
rz(-3.0070947) q[2];
sx q[2];
rz(2.4328977) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0811942) q[1];
sx q[1];
rz(-1.8738998) q[1];
sx q[1];
rz(-2.7737507) q[1];
x q[2];
rz(-1.2358642) q[3];
sx q[3];
rz(-0.55509243) q[3];
sx q[3];
rz(-2.5344283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7572299) q[2];
sx q[2];
rz(-2.195916) q[2];
sx q[2];
rz(-0.45822701) q[2];
rz(-0.630817) q[3];
sx q[3];
rz(-1.2354555) q[3];
sx q[3];
rz(-2.6423776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45796564) q[0];
sx q[0];
rz(-2.2891335) q[0];
sx q[0];
rz(-3.1347347) q[0];
rz(0.231617) q[1];
sx q[1];
rz(-1.3791142) q[1];
sx q[1];
rz(-2.0793656) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4214913) q[0];
sx q[0];
rz(-1.694223) q[0];
sx q[0];
rz(2.8272259) q[0];
rz(-0.3756282) q[2];
sx q[2];
rz(-2.6115719) q[2];
sx q[2];
rz(-2.9935868) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0343218) q[1];
sx q[1];
rz(-2.2244029) q[1];
sx q[1];
rz(-1.8922217) q[1];
x q[2];
rz(-1.7627349) q[3];
sx q[3];
rz(-2.4232691) q[3];
sx q[3];
rz(-3.0305221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.14083938) q[2];
sx q[2];
rz(-1.2522298) q[2];
sx q[2];
rz(-1.8737277) q[2];
rz(-0.86152348) q[3];
sx q[3];
rz(-2.0778766) q[3];
sx q[3];
rz(-1.4388194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
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
rz(0.93290257) q[0];
sx q[0];
rz(-2.8080495) q[0];
sx q[0];
rz(2.9260337) q[0];
rz(-2.6453099) q[1];
sx q[1];
rz(-1.9535306) q[1];
sx q[1];
rz(-2.9821679) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8529786) q[0];
sx q[0];
rz(-2.5067177) q[0];
sx q[0];
rz(-1.575241) q[0];
x q[1];
rz(1.764545) q[2];
sx q[2];
rz(-2.0341349) q[2];
sx q[2];
rz(-3.0510356) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1856975) q[1];
sx q[1];
rz(-1.5296827) q[1];
sx q[1];
rz(1.3354882) q[1];
rz(-pi) q[2];
rz(1.0219021) q[3];
sx q[3];
rz(-1.1682142) q[3];
sx q[3];
rz(0.91820133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7364007) q[2];
sx q[2];
rz(-2.0173732) q[2];
sx q[2];
rz(2.6250725) q[2];
rz(-1.9308331) q[3];
sx q[3];
rz(-1.6885933) q[3];
sx q[3];
rz(1.7621382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6571758) q[0];
sx q[0];
rz(-0.076198904) q[0];
sx q[0];
rz(-0.54022378) q[0];
rz(-1.6172488) q[1];
sx q[1];
rz(-2.1397619) q[1];
sx q[1];
rz(-1.2670955) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0277242) q[0];
sx q[0];
rz(-2.0696215) q[0];
sx q[0];
rz(-1.3214825) q[0];
x q[1];
rz(-0.77719633) q[2];
sx q[2];
rz(-2.1146449) q[2];
sx q[2];
rz(1.0842619) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9766751) q[1];
sx q[1];
rz(-1.4265713) q[1];
sx q[1];
rz(-0.37101908) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1096051) q[3];
sx q[3];
rz(-0.69503419) q[3];
sx q[3];
rz(-1.1170596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2544864) q[2];
sx q[2];
rz(-1.1595414) q[2];
sx q[2];
rz(-2.9642588) q[2];
rz(-1.0698498) q[3];
sx q[3];
rz(-2.0900574) q[3];
sx q[3];
rz(0.18226084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3607445) q[0];
sx q[0];
rz(-1.1540664) q[0];
sx q[0];
rz(0.10072197) q[0];
rz(1.1489457) q[1];
sx q[1];
rz(-1.9457685) q[1];
sx q[1];
rz(1.1740059) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74885313) q[0];
sx q[0];
rz(-2.1707188) q[0];
sx q[0];
rz(-1.3412745) q[0];
rz(-pi) q[1];
rz(0.10020013) q[2];
sx q[2];
rz(-0.76876193) q[2];
sx q[2];
rz(0.59610808) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2795978) q[1];
sx q[1];
rz(-1.5513541) q[1];
sx q[1];
rz(1.6085546) q[1];
rz(0.10281201) q[3];
sx q[3];
rz(-0.88231707) q[3];
sx q[3];
rz(-1.7243054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9775057) q[2];
sx q[2];
rz(-1.4848494) q[2];
sx q[2];
rz(0.40317765) q[2];
rz(2.4781503) q[3];
sx q[3];
rz(-0.57938975) q[3];
sx q[3];
rz(-3.1373533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4106301) q[0];
sx q[0];
rz(-1.0799438) q[0];
sx q[0];
rz(-0.078911111) q[0];
rz(-0.90947378) q[1];
sx q[1];
rz(-1.0458922) q[1];
sx q[1];
rz(-0.39631072) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3672597) q[0];
sx q[0];
rz(-3.1235187) q[0];
sx q[0];
rz(-2.0967183) q[0];
rz(3.134507) q[2];
sx q[2];
rz(-2.3478697) q[2];
sx q[2];
rz(0.34863472) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7381607) q[1];
sx q[1];
rz(-0.53817525) q[1];
sx q[1];
rz(-0.47026431) q[1];
x q[2];
rz(-0.5509866) q[3];
sx q[3];
rz(-1.0240304) q[3];
sx q[3];
rz(-1.3133698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4497455) q[2];
sx q[2];
rz(-2.2578466) q[2];
sx q[2];
rz(0.79908243) q[2];
rz(-0.07829047) q[3];
sx q[3];
rz(-1.2543863) q[3];
sx q[3];
rz(-2.4962795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4770724) q[0];
sx q[0];
rz(-1.399853) q[0];
sx q[0];
rz(0.54367263) q[0];
rz(-1.9480582) q[1];
sx q[1];
rz(-1.8652893) q[1];
sx q[1];
rz(0.51020772) q[1];
rz(1.1134182) q[2];
sx q[2];
rz(-1.7241679) q[2];
sx q[2];
rz(0.67571251) q[2];
rz(-0.54936784) q[3];
sx q[3];
rz(-2.3372169) q[3];
sx q[3];
rz(-1.611471) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
