OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2774529) q[0];
sx q[0];
rz(4.6946445) q[0];
sx q[0];
rz(10.932218) q[0];
rz(1.5965257) q[1];
sx q[1];
rz(-0.59626055) q[1];
sx q[1];
rz(-2.526386) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9278487) q[0];
sx q[0];
rz(-2.1670177) q[0];
sx q[0];
rz(0.56791373) q[0];
x q[1];
rz(-2.7841714) q[2];
sx q[2];
rz(-1.7454141) q[2];
sx q[2];
rz(-1.071196) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2088036) q[1];
sx q[1];
rz(-2.1425769) q[1];
sx q[1];
rz(0.51704452) q[1];
rz(-pi) q[2];
rz(2.8350713) q[3];
sx q[3];
rz(-1.3984826) q[3];
sx q[3];
rz(-2.4337208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.4404099) q[2];
sx q[2];
rz(-1.5298693) q[2];
sx q[2];
rz(0.33828503) q[2];
rz(-1.7017378) q[3];
sx q[3];
rz(-2.2262636) q[3];
sx q[3];
rz(-2.2556944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97025362) q[0];
sx q[0];
rz(-2.4304424) q[0];
sx q[0];
rz(-3.1112444) q[0];
rz(-3.0753823) q[1];
sx q[1];
rz(-0.98774424) q[1];
sx q[1];
rz(-1.5240086) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0328007) q[0];
sx q[0];
rz(-1.5691225) q[0];
sx q[0];
rz(-1.770442) q[0];
rz(-pi) q[1];
rz(1.1163887) q[2];
sx q[2];
rz(-1.589236) q[2];
sx q[2];
rz(-1.6905284) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6807032) q[1];
sx q[1];
rz(-1.5017121) q[1];
sx q[1];
rz(-2.1417888) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1094692) q[3];
sx q[3];
rz(-0.20878775) q[3];
sx q[3];
rz(-1.3267645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7559738) q[2];
sx q[2];
rz(-0.9884584) q[2];
sx q[2];
rz(-1.9937817) q[2];
rz(1.8148445) q[3];
sx q[3];
rz(-1.3245405) q[3];
sx q[3];
rz(0.23708788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1266992) q[0];
sx q[0];
rz(-0.48148695) q[0];
sx q[0];
rz(2.8258064) q[0];
rz(-2.2029927) q[1];
sx q[1];
rz(-1.6789852) q[1];
sx q[1];
rz(-2.8895203) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7092428) q[0];
sx q[0];
rz(-0.81626695) q[0];
sx q[0];
rz(-0.66246756) q[0];
x q[1];
rz(-2.2601068) q[2];
sx q[2];
rz(-0.93886095) q[2];
sx q[2];
rz(1.880868) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7604916) q[1];
sx q[1];
rz(-0.64844202) q[1];
sx q[1];
rz(-1.2566503) q[1];
rz(-2.4113703) q[3];
sx q[3];
rz(-0.50958868) q[3];
sx q[3];
rz(-0.93917055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1217653) q[2];
sx q[2];
rz(-2.6941507) q[2];
sx q[2];
rz(-0.034051731) q[2];
rz(-0.017459067) q[3];
sx q[3];
rz(-1.358946) q[3];
sx q[3];
rz(-1.0954558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62717342) q[0];
sx q[0];
rz(-1.592941) q[0];
sx q[0];
rz(1.6148286) q[0];
rz(-1.0871672) q[1];
sx q[1];
rz(-0.68030578) q[1];
sx q[1];
rz(0.70708752) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3084761) q[0];
sx q[0];
rz(-1.1261228) q[0];
sx q[0];
rz(-1.1348669) q[0];
rz(0.25755067) q[2];
sx q[2];
rz(-2.6817245) q[2];
sx q[2];
rz(0.79007733) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6882119) q[1];
sx q[1];
rz(-2.4773295) q[1];
sx q[1];
rz(1.7654256) q[1];
rz(-pi) q[2];
rz(1.8791734) q[3];
sx q[3];
rz(-1.7763419) q[3];
sx q[3];
rz(1.8431078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.84918555) q[2];
sx q[2];
rz(-1.9533998) q[2];
sx q[2];
rz(-0.46009955) q[2];
rz(1.7442616) q[3];
sx q[3];
rz(-1.5887235) q[3];
sx q[3];
rz(-2.8989255) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8206772) q[0];
sx q[0];
rz(-2.9292332) q[0];
sx q[0];
rz(1.7472349) q[0];
rz(-2.0460515) q[1];
sx q[1];
rz(-1.6004326) q[1];
sx q[1];
rz(0.25462338) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6760315) q[0];
sx q[0];
rz(-1.6499632) q[0];
sx q[0];
rz(0.013750793) q[0];
rz(1.5484372) q[2];
sx q[2];
rz(-2.5639113) q[2];
sx q[2];
rz(-0.81800848) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4659781) q[1];
sx q[1];
rz(-1.8120159) q[1];
sx q[1];
rz(2.4720008) q[1];
rz(2.280064) q[3];
sx q[3];
rz(-1.2332752) q[3];
sx q[3];
rz(0.57818128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1191117) q[2];
sx q[2];
rz(-0.20038651) q[2];
sx q[2];
rz(-1.3767892) q[2];
rz(1.4962176) q[3];
sx q[3];
rz(-1.6330556) q[3];
sx q[3];
rz(-1.0866603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45143932) q[0];
sx q[0];
rz(-2.4017161) q[0];
sx q[0];
rz(2.8421463) q[0];
rz(2.1014138) q[1];
sx q[1];
rz(-1.6957915) q[1];
sx q[1];
rz(0.20656955) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9402007) q[0];
sx q[0];
rz(-2.3513146) q[0];
sx q[0];
rz(1.9076365) q[0];
x q[1];
rz(-1.0429522) q[2];
sx q[2];
rz(-1.9163418) q[2];
sx q[2];
rz(0.31574677) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5094604) q[1];
sx q[1];
rz(-1.7473979) q[1];
sx q[1];
rz(-0.86061865) q[1];
rz(-pi) q[2];
rz(1.9566831) q[3];
sx q[3];
rz(-2.4486802) q[3];
sx q[3];
rz(0.087156765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.841659) q[2];
sx q[2];
rz(-1.724023) q[2];
sx q[2];
rz(-2.858813) q[2];
rz(0.81280604) q[3];
sx q[3];
rz(-2.7225284) q[3];
sx q[3];
rz(0.16684428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92641002) q[0];
sx q[0];
rz(-1.0535425) q[0];
sx q[0];
rz(2.7600631) q[0];
rz(-2.5577257) q[1];
sx q[1];
rz(-0.54324141) q[1];
sx q[1];
rz(-1.8136224) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5751942) q[0];
sx q[0];
rz(-2.0090721) q[0];
sx q[0];
rz(-3.0153494) q[0];
rz(-pi) q[1];
rz(2.1930201) q[2];
sx q[2];
rz(-1.9554536) q[2];
sx q[2];
rz(-2.4228061) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2139637) q[1];
sx q[1];
rz(-2.2447526) q[1];
sx q[1];
rz(-0.2143292) q[1];
rz(-pi) q[2];
rz(1.1442723) q[3];
sx q[3];
rz(-2.0445619) q[3];
sx q[3];
rz(1.1748479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.61838377) q[2];
sx q[2];
rz(-2.0760459) q[2];
sx q[2];
rz(1.3605114) q[2];
rz(-1.7112188) q[3];
sx q[3];
rz(-1.0083219) q[3];
sx q[3];
rz(0.84806228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9946063) q[0];
sx q[0];
rz(-1.9817579) q[0];
sx q[0];
rz(-2.9597136) q[0];
rz(-2.6673642) q[1];
sx q[1];
rz(-1.0206181) q[1];
sx q[1];
rz(2.1906733) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61688214) q[0];
sx q[0];
rz(-1.2843772) q[0];
sx q[0];
rz(2.8911203) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0792564) q[2];
sx q[2];
rz(-1.8364292) q[2];
sx q[2];
rz(3.0692284) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2434395) q[1];
sx q[1];
rz(-2.9228518) q[1];
sx q[1];
rz(2.8380413) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3202053) q[3];
sx q[3];
rz(-2.6139724) q[3];
sx q[3];
rz(-1.9963095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9178847) q[2];
sx q[2];
rz(-2.6612838) q[2];
sx q[2];
rz(0.075909464) q[2];
rz(2.5935796) q[3];
sx q[3];
rz(-1.3242105) q[3];
sx q[3];
rz(0.63265911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(0.52371812) q[0];
sx q[0];
rz(-2.0848367) q[0];
sx q[0];
rz(-1.3762208) q[0];
rz(-0.41704047) q[1];
sx q[1];
rz(-1.4191671) q[1];
sx q[1];
rz(2.4818647) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2840246) q[0];
sx q[0];
rz(-2.095247) q[0];
sx q[0];
rz(-2.2542623) q[0];
x q[1];
rz(1.7890036) q[2];
sx q[2];
rz(-2.3404684) q[2];
sx q[2];
rz(-1.4584691) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8068741) q[1];
sx q[1];
rz(-1.8611307) q[1];
sx q[1];
rz(0.90805407) q[1];
rz(-0.057007313) q[3];
sx q[3];
rz(-2.1087286) q[3];
sx q[3];
rz(-1.2008592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.62347162) q[2];
sx q[2];
rz(-0.76449624) q[2];
sx q[2];
rz(-1.0260322) q[2];
rz(-2.9927411) q[3];
sx q[3];
rz(-2.1089349) q[3];
sx q[3];
rz(3.1159475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.898107) q[0];
sx q[0];
rz(-0.74111104) q[0];
sx q[0];
rz(-1.8359258) q[0];
rz(1.1765515) q[1];
sx q[1];
rz(-1.2780317) q[1];
sx q[1];
rz(2.1059039) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64338387) q[0];
sx q[0];
rz(-1.9783101) q[0];
sx q[0];
rz(1.9252752) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5752441) q[2];
sx q[2];
rz(-2.1841335) q[2];
sx q[2];
rz(1.7042421) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4020821) q[1];
sx q[1];
rz(-1.2716736) q[1];
sx q[1];
rz(-1.8226536) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3003179) q[3];
sx q[3];
rz(-0.83993739) q[3];
sx q[3];
rz(-0.65281103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5130561) q[2];
sx q[2];
rz(-0.81130242) q[2];
sx q[2];
rz(1.9899842) q[2];
rz(2.628905) q[3];
sx q[3];
rz(-1.0980462) q[3];
sx q[3];
rz(-2.776896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.37968996) q[0];
sx q[0];
rz(-1.7871478) q[0];
sx q[0];
rz(-2.3085069) q[0];
rz(-1.5079386) q[1];
sx q[1];
rz(-2.5588551) q[1];
sx q[1];
rz(-0.48164639) q[1];
rz(-0.41082906) q[2];
sx q[2];
rz(-1.5317691) q[2];
sx q[2];
rz(-0.017824235) q[2];
rz(2.5232265) q[3];
sx q[3];
rz(-1.5636087) q[3];
sx q[3];
rz(-2.9611361) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
