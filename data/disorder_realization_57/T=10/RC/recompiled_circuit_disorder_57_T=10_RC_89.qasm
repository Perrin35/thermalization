OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.52743071) q[0];
sx q[0];
rz(-2.3318113) q[0];
sx q[0];
rz(0.53139395) q[0];
rz(0.2375138) q[1];
sx q[1];
rz(-1.7778492) q[1];
sx q[1];
rz(1.9385424) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8578313) q[0];
sx q[0];
rz(-2.8537822) q[0];
sx q[0];
rz(2.3057111) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4389624) q[2];
sx q[2];
rz(-2.812817) q[2];
sx q[2];
rz(1.6871014) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0359218) q[1];
sx q[1];
rz(-2.8997313) q[1];
sx q[1];
rz(-0.24005228) q[1];
rz(-pi) q[2];
x q[2];
rz(0.12227998) q[3];
sx q[3];
rz(-0.29502007) q[3];
sx q[3];
rz(2.2515841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8649341) q[2];
sx q[2];
rz(-1.1316789) q[2];
sx q[2];
rz(-0.19908389) q[2];
rz(-1.52786) q[3];
sx q[3];
rz(-2.644643) q[3];
sx q[3];
rz(-2.4075107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99877015) q[0];
sx q[0];
rz(-3.0144189) q[0];
sx q[0];
rz(2.3221827) q[0];
rz(0.2858513) q[1];
sx q[1];
rz(-2.2276623) q[1];
sx q[1];
rz(1.2664638) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0430543) q[0];
sx q[0];
rz(-0.96999723) q[0];
sx q[0];
rz(-1.1061125) q[0];
rz(-pi) q[1];
rz(-1.9752713) q[2];
sx q[2];
rz(-1.7911439) q[2];
sx q[2];
rz(-1.599556) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.61050615) q[1];
sx q[1];
rz(-0.55263457) q[1];
sx q[1];
rz(-2.3081739) q[1];
rz(-pi) q[2];
rz(-0.5760848) q[3];
sx q[3];
rz(-1.1844716) q[3];
sx q[3];
rz(-1.9379804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9849898) q[2];
sx q[2];
rz(-1.2133602) q[2];
sx q[2];
rz(1.9796237) q[2];
rz(-0.087163838) q[3];
sx q[3];
rz(-1.5036843) q[3];
sx q[3];
rz(3.0676837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6067628) q[0];
sx q[0];
rz(-2.02089) q[0];
sx q[0];
rz(-2.8379922) q[0];
rz(1.3820232) q[1];
sx q[1];
rz(-1.8654114) q[1];
sx q[1];
rz(0.64750013) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8436463) q[0];
sx q[0];
rz(-0.93695153) q[0];
sx q[0];
rz(1.514939) q[0];
rz(-pi) q[1];
rz(2.5325003) q[2];
sx q[2];
rz(-0.81431544) q[2];
sx q[2];
rz(2.7012205) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.69933575) q[1];
sx q[1];
rz(-1.934821) q[1];
sx q[1];
rz(1.0815716) q[1];
rz(-pi) q[2];
rz(-2.2494499) q[3];
sx q[3];
rz(-0.8222848) q[3];
sx q[3];
rz(-3.126614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2993762) q[2];
sx q[2];
rz(-0.78263485) q[2];
sx q[2];
rz(-1.5608609) q[2];
rz(2.1598699) q[3];
sx q[3];
rz(-1.862062) q[3];
sx q[3];
rz(-0.64341199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9075539) q[0];
sx q[0];
rz(-0.83775318) q[0];
sx q[0];
rz(2.2696944) q[0];
rz(-2.7543228) q[1];
sx q[1];
rz(-0.64278066) q[1];
sx q[1];
rz(0.29104582) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1211348) q[0];
sx q[0];
rz(-0.99636787) q[0];
sx q[0];
rz(2.8681884) q[0];
rz(0.29291885) q[2];
sx q[2];
rz(-1.8790763) q[2];
sx q[2];
rz(2.3777547) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2070771) q[1];
sx q[1];
rz(-2.0743437) q[1];
sx q[1];
rz(0.7657004) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8157186) q[3];
sx q[3];
rz(-0.59026736) q[3];
sx q[3];
rz(-0.11412379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3691833) q[2];
sx q[2];
rz(-0.79493752) q[2];
sx q[2];
rz(2.5975361) q[2];
rz(-2.6323075) q[3];
sx q[3];
rz(-1.0638758) q[3];
sx q[3];
rz(-2.2756186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13957025) q[0];
sx q[0];
rz(-1.0263356) q[0];
sx q[0];
rz(-2.496526) q[0];
rz(2.5158665) q[1];
sx q[1];
rz(-1.9179683) q[1];
sx q[1];
rz(-0.63017875) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.038534855) q[0];
sx q[0];
rz(-1.4487421) q[0];
sx q[0];
rz(-0.5794258) q[0];
x q[1];
rz(-2.3216256) q[2];
sx q[2];
rz(-1.2756639) q[2];
sx q[2];
rz(-0.68857869) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.6781792) q[1];
sx q[1];
rz(-0.58600512) q[1];
sx q[1];
rz(1.9834118) q[1];
rz(-pi) q[2];
rz(2.583902) q[3];
sx q[3];
rz(-1.7002749) q[3];
sx q[3];
rz(-1.492471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.75382918) q[2];
sx q[2];
rz(-1.3239653) q[2];
sx q[2];
rz(-1.3641664) q[2];
rz(-1.8917313) q[3];
sx q[3];
rz(-1.9259689) q[3];
sx q[3];
rz(-2.2201339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0867778) q[0];
sx q[0];
rz(-0.58371109) q[0];
sx q[0];
rz(0.71682799) q[0];
rz(-0.43771276) q[1];
sx q[1];
rz(-2.4611459) q[1];
sx q[1];
rz(-0.81370083) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5649367) q[0];
sx q[0];
rz(-1.7840958) q[0];
sx q[0];
rz(-3.0773786) q[0];
x q[1];
rz(0.72197638) q[2];
sx q[2];
rz(-2.0275896) q[2];
sx q[2];
rz(1.7320088) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5972283) q[1];
sx q[1];
rz(-1.9042943) q[1];
sx q[1];
rz(-0.79798214) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2171451) q[3];
sx q[3];
rz(-1.8600703) q[3];
sx q[3];
rz(-1.7464964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8852691) q[2];
sx q[2];
rz(-2.7041114) q[2];
sx q[2];
rz(-1.9539333) q[2];
rz(-2.2096283) q[3];
sx q[3];
rz(-2.0075802) q[3];
sx q[3];
rz(-2.7704172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1535783) q[0];
sx q[0];
rz(-1.0289285) q[0];
sx q[0];
rz(-0.6860835) q[0];
rz(1.6185121) q[1];
sx q[1];
rz(-2.1673817) q[1];
sx q[1];
rz(-2.4783321) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8406869) q[0];
sx q[0];
rz(-1.5468883) q[0];
sx q[0];
rz(2.2211214) q[0];
rz(1.7887605) q[2];
sx q[2];
rz(-2.1969165) q[2];
sx q[2];
rz(1.2499714) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2901193) q[1];
sx q[1];
rz(-2.0514601) q[1];
sx q[1];
rz(0.22766797) q[1];
x q[2];
rz(-2.456326) q[3];
sx q[3];
rz(-2.4589834) q[3];
sx q[3];
rz(0.21041378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.66701165) q[2];
sx q[2];
rz(-0.8374927) q[2];
sx q[2];
rz(3.1265124) q[2];
rz(1.0117426) q[3];
sx q[3];
rz(-2.0254617) q[3];
sx q[3];
rz(2.570178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0886154) q[0];
sx q[0];
rz(-2.4089854) q[0];
sx q[0];
rz(-0.10738871) q[0];
rz(2.836851) q[1];
sx q[1];
rz(-1.823103) q[1];
sx q[1];
rz(0.4577786) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4360355) q[0];
sx q[0];
rz(-1.2495263) q[0];
sx q[0];
rz(2.2309488) q[0];
rz(0.56577487) q[2];
sx q[2];
rz(-2.3104295) q[2];
sx q[2];
rz(0.6643675) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0929872) q[1];
sx q[1];
rz(-2.7106206) q[1];
sx q[1];
rz(0.98577568) q[1];
x q[2];
rz(2.5910208) q[3];
sx q[3];
rz(-0.62567657) q[3];
sx q[3];
rz(-1.5284556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.36499873) q[2];
sx q[2];
rz(-1.685131) q[2];
sx q[2];
rz(-1.7101074) q[2];
rz(-1.7163904) q[3];
sx q[3];
rz(-2.5148354) q[3];
sx q[3];
rz(1.2861929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9872221) q[0];
sx q[0];
rz(-2.6619338) q[0];
sx q[0];
rz(2.1881058) q[0];
rz(1.3257239) q[1];
sx q[1];
rz(-2.6486501) q[1];
sx q[1];
rz(0.58473933) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23823243) q[0];
sx q[0];
rz(-2.0582709) q[0];
sx q[0];
rz(1.6977298) q[0];
rz(-pi) q[1];
rz(-1.2051177) q[2];
sx q[2];
rz(-2.0846016) q[2];
sx q[2];
rz(-0.50525451) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.016420267) q[1];
sx q[1];
rz(-2.6965045) q[1];
sx q[1];
rz(2.9801286) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3034079) q[3];
sx q[3];
rz(-1.4328453) q[3];
sx q[3];
rz(1.4064058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8087625) q[2];
sx q[2];
rz(-2.3213883) q[2];
sx q[2];
rz(-0.27627036) q[2];
rz(-0.55109465) q[3];
sx q[3];
rz(-1.3994183) q[3];
sx q[3];
rz(0.38366693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0722512) q[0];
sx q[0];
rz(-2.6273917) q[0];
sx q[0];
rz(2.4861091) q[0];
rz(1.9845225) q[1];
sx q[1];
rz(-1.7600704) q[1];
sx q[1];
rz(1.6190593) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5009506) q[0];
sx q[0];
rz(-2.2146533) q[0];
sx q[0];
rz(1.180091) q[0];
rz(-pi) q[1];
rz(2.4535975) q[2];
sx q[2];
rz(-1.3032459) q[2];
sx q[2];
rz(-0.5984532) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9138227) q[1];
sx q[1];
rz(-1.7350983) q[1];
sx q[1];
rz(2.5096202) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7195815) q[3];
sx q[3];
rz(-1.9415628) q[3];
sx q[3];
rz(1.6251723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.15554252) q[2];
sx q[2];
rz(-0.4077929) q[2];
sx q[2];
rz(-2.299451) q[2];
rz(1.6048253) q[3];
sx q[3];
rz(-1.7611046) q[3];
sx q[3];
rz(3.1150637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90606541) q[0];
sx q[0];
rz(-1.1885831) q[0];
sx q[0];
rz(2.6901235) q[0];
rz(-0.72262598) q[1];
sx q[1];
rz(-2.2650748) q[1];
sx q[1];
rz(1.5320019) q[1];
rz(1.7799829) q[2];
sx q[2];
rz(-1.5968512) q[2];
sx q[2];
rz(-3.130393) q[2];
rz(-3.0620861) q[3];
sx q[3];
rz(-0.88149298) q[3];
sx q[3];
rz(0.029824921) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
