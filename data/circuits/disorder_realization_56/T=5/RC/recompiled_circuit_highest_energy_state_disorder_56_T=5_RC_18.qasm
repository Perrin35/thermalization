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
rz(1.4752969) q[0];
sx q[0];
rz(-1.2694321) q[0];
sx q[0];
rz(-2.4694634) q[0];
rz(0.44828662) q[1];
sx q[1];
rz(-1.5223794) q[1];
sx q[1];
rz(0.7575922) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8920994) q[0];
sx q[0];
rz(-1.4612911) q[0];
sx q[0];
rz(0.076105781) q[0];
rz(-pi) q[1];
rz(0.50067164) q[2];
sx q[2];
rz(-1.2462052) q[2];
sx q[2];
rz(0.90510923) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1266844) q[1];
sx q[1];
rz(-1.9796625) q[1];
sx q[1];
rz(1.8168628) q[1];
rz(2.4663062) q[3];
sx q[3];
rz(-0.8199586) q[3];
sx q[3];
rz(-0.91778558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0932833) q[2];
sx q[2];
rz(-1.8123241) q[2];
sx q[2];
rz(-1.7993571) q[2];
rz(-1.7736769) q[3];
sx q[3];
rz(-1.0932837) q[3];
sx q[3];
rz(-1.5526519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9369478) q[0];
sx q[0];
rz(-1.0138252) q[0];
sx q[0];
rz(-0.94888765) q[0];
rz(0.10143796) q[1];
sx q[1];
rz(-1.0600435) q[1];
sx q[1];
rz(-2.2116275) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6008015) q[0];
sx q[0];
rz(-1.7834429) q[0];
sx q[0];
rz(1.6746816) q[0];
x q[1];
rz(-1.5632767) q[2];
sx q[2];
rz(-2.0759656) q[2];
sx q[2];
rz(-1.9746321) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8023493) q[1];
sx q[1];
rz(-1.7319458) q[1];
sx q[1];
rz(-1.6768964) q[1];
x q[2];
rz(-2.5245776) q[3];
sx q[3];
rz(-0.91147826) q[3];
sx q[3];
rz(2.3544745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0019504) q[2];
sx q[2];
rz(-1.4652239) q[2];
sx q[2];
rz(0.742221) q[2];
rz(2.6774075) q[3];
sx q[3];
rz(-1.3451385) q[3];
sx q[3];
rz(1.5884885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6716229) q[0];
sx q[0];
rz(-0.18593423) q[0];
sx q[0];
rz(-1.6999014) q[0];
rz(0.16547671) q[1];
sx q[1];
rz(-2.4022357) q[1];
sx q[1];
rz(-2.2742719) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9242223) q[0];
sx q[0];
rz(-1.570756) q[0];
sx q[0];
rz(1.5702308) q[0];
rz(-pi) q[1];
rz(-1.0652415) q[2];
sx q[2];
rz(-1.3320413) q[2];
sx q[2];
rz(2.3468034) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8738385) q[1];
sx q[1];
rz(-1.2902564) q[1];
sx q[1];
rz(-0.33919097) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3300784) q[3];
sx q[3];
rz(-1.1574739) q[3];
sx q[3];
rz(2.1275525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9637588) q[2];
sx q[2];
rz(-1.9671665) q[2];
sx q[2];
rz(-2.3826694) q[2];
rz(2.3376076) q[3];
sx q[3];
rz(-2.1920429) q[3];
sx q[3];
rz(0.72833958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3114965) q[0];
sx q[0];
rz(-1.281597) q[0];
sx q[0];
rz(-0.48402825) q[0];
rz(-1.5361891) q[1];
sx q[1];
rz(-1.4151298) q[1];
sx q[1];
rz(1.7162292) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1041906) q[0];
sx q[0];
rz(-2.0811102) q[0];
sx q[0];
rz(-1.181382) q[0];
rz(-pi) q[1];
x q[1];
rz(0.31230782) q[2];
sx q[2];
rz(-1.3538401) q[2];
sx q[2];
rz(-1.7157451) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0565383) q[1];
sx q[1];
rz(-2.1966272) q[1];
sx q[1];
rz(-0.27295785) q[1];
rz(-pi) q[2];
rz(1.1953765) q[3];
sx q[3];
rz(-1.1470801) q[3];
sx q[3];
rz(2.8590607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.385685) q[2];
sx q[2];
rz(-1.0681095) q[2];
sx q[2];
rz(-0.29328406) q[2];
rz(-0.048132345) q[3];
sx q[3];
rz(-1.8354514) q[3];
sx q[3];
rz(2.8369246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3047979) q[0];
sx q[0];
rz(-1.7339107) q[0];
sx q[0];
rz(-1.1037214) q[0];
rz(-0.91148218) q[1];
sx q[1];
rz(-1.6171004) q[1];
sx q[1];
rz(0.10890659) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86065642) q[0];
sx q[0];
rz(-1.0901648) q[0];
sx q[0];
rz(2.3970766) q[0];
x q[1];
rz(0.53498603) q[2];
sx q[2];
rz(-1.4145383) q[2];
sx q[2];
rz(1.5076758) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0518347) q[1];
sx q[1];
rz(-1.4076774) q[1];
sx q[1];
rz(2.240313) q[1];
rz(1.2690684) q[3];
sx q[3];
rz(-1.0142039) q[3];
sx q[3];
rz(-2.3315786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.64985046) q[2];
sx q[2];
rz(-1.3567341) q[2];
sx q[2];
rz(0.25807992) q[2];
rz(-2.7598925) q[3];
sx q[3];
rz(-2.2038867) q[3];
sx q[3];
rz(0.55735731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7447164) q[0];
sx q[0];
rz(-0.98137403) q[0];
sx q[0];
rz(-1.1599524) q[0];
rz(-0.57811919) q[1];
sx q[1];
rz(-1.4719529) q[1];
sx q[1];
rz(0.32346183) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2787298) q[0];
sx q[0];
rz(-0.6375618) q[0];
sx q[0];
rz(0.80018534) q[0];
x q[1];
rz(-0.49655621) q[2];
sx q[2];
rz(-1.7944031) q[2];
sx q[2];
rz(-1.8211435) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.19666323) q[1];
sx q[1];
rz(-0.62360901) q[1];
sx q[1];
rz(2.7476279) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6810535) q[3];
sx q[3];
rz(-1.2804739) q[3];
sx q[3];
rz(3.0825305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9385927) q[2];
sx q[2];
rz(-2.3679569) q[2];
sx q[2];
rz(-2.2933551) q[2];
rz(-1.2674468) q[3];
sx q[3];
rz(-1.7402612) q[3];
sx q[3];
rz(2.447824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0271725) q[0];
sx q[0];
rz(-2.4878451) q[0];
sx q[0];
rz(2.2681336) q[0];
rz(1.9050441) q[1];
sx q[1];
rz(-2.1658587) q[1];
sx q[1];
rz(-0.083560856) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1655052) q[0];
sx q[0];
rz(-0.3437316) q[0];
sx q[0];
rz(-2.7718839) q[0];
rz(-pi) q[1];
rz(-1.0413076) q[2];
sx q[2];
rz(-1.7779967) q[2];
sx q[2];
rz(2.4251314) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4518731) q[1];
sx q[1];
rz(-2.3102133) q[1];
sx q[1];
rz(-1.0240133) q[1];
rz(-2.7807981) q[3];
sx q[3];
rz(-1.2865781) q[3];
sx q[3];
rz(0.36297114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7053335) q[2];
sx q[2];
rz(-1.6062364) q[2];
sx q[2];
rz(-1.3667038) q[2];
rz(-2.8691835) q[3];
sx q[3];
rz(-1.0206157) q[3];
sx q[3];
rz(-0.75850707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90683872) q[0];
sx q[0];
rz(-2.315157) q[0];
sx q[0];
rz(0.93836623) q[0];
rz(0.80288184) q[1];
sx q[1];
rz(-2.1736841) q[1];
sx q[1];
rz(-0.82069194) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7308189) q[0];
sx q[0];
rz(-1.2026498) q[0];
sx q[0];
rz(2.9981722) q[0];
x q[1];
rz(2.9740779) q[2];
sx q[2];
rz(-1.2021087) q[2];
sx q[2];
rz(0.14376727) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8269315) q[1];
sx q[1];
rz(-2.1912327) q[1];
sx q[1];
rz(0.13038306) q[1];
rz(-pi) q[2];
rz(1.197416) q[3];
sx q[3];
rz(-2.7146517) q[3];
sx q[3];
rz(0.77746848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.20251033) q[2];
sx q[2];
rz(-2.9896917) q[2];
sx q[2];
rz(2.0738156) q[2];
rz(2.0104525) q[3];
sx q[3];
rz(-2.307297) q[3];
sx q[3];
rz(-3.06156) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15602569) q[0];
sx q[0];
rz(-1.9859059) q[0];
sx q[0];
rz(0.40618968) q[0];
rz(-1.4093026) q[1];
sx q[1];
rz(-2.1235762) q[1];
sx q[1];
rz(1.0850151) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99434851) q[0];
sx q[0];
rz(-1.3965333) q[0];
sx q[0];
rz(-0.6263349) q[0];
x q[1];
rz(2.5312838) q[2];
sx q[2];
rz(-1.5670735) q[2];
sx q[2];
rz(-0.10476724) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0250562) q[1];
sx q[1];
rz(-0.79306125) q[1];
sx q[1];
rz(2.039394) q[1];
rz(-2.1187339) q[3];
sx q[3];
rz(-1.0617439) q[3];
sx q[3];
rz(-1.7804543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5857508) q[2];
sx q[2];
rz(-2.5297574) q[2];
sx q[2];
rz(-2.2108868) q[2];
rz(-1.7454923) q[3];
sx q[3];
rz(-1.7788818) q[3];
sx q[3];
rz(3.0220368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4569106) q[0];
sx q[0];
rz(-0.97701183) q[0];
sx q[0];
rz(-1.8344301) q[0];
rz(2.7556509) q[1];
sx q[1];
rz(-1.394505) q[1];
sx q[1];
rz(2.1070259) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32337727) q[0];
sx q[0];
rz(-1.744483) q[0];
sx q[0];
rz(-3.0370569) q[0];
rz(-pi) q[1];
rz(-2.6808591) q[2];
sx q[2];
rz(-1.3855532) q[2];
sx q[2];
rz(-0.68250193) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.17321302) q[1];
sx q[1];
rz(-0.82260859) q[1];
sx q[1];
rz(-2.2065225) q[1];
rz(3.0183082) q[3];
sx q[3];
rz(-0.81466952) q[3];
sx q[3];
rz(-2.2565763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4474386) q[2];
sx q[2];
rz(-2.4093781) q[2];
sx q[2];
rz(0.36699692) q[2];
rz(-2.0224109) q[3];
sx q[3];
rz(-2.7435591) q[3];
sx q[3];
rz(-1.5252349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3601396) q[0];
sx q[0];
rz(-1.9522788) q[0];
sx q[0];
rz(2.0576394) q[0];
rz(3.0837334) q[1];
sx q[1];
rz(-1.8744938) q[1];
sx q[1];
rz(-1.4300463) q[1];
rz(2.9686676) q[2];
sx q[2];
rz(-1.4739707) q[2];
sx q[2];
rz(-2.9116918) q[2];
rz(-0.41856159) q[3];
sx q[3];
rz(-0.7444612) q[3];
sx q[3];
rz(0.85519467) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
