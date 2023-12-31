OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.86413971) q[0];
sx q[0];
rz(-1.5530518) q[0];
sx q[0];
rz(1.6341524) q[0];
rz(1.5965257) q[1];
sx q[1];
rz(2.5453321) q[1];
sx q[1];
rz(8.8095713) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.213744) q[0];
sx q[0];
rz(-0.97457492) q[0];
sx q[0];
rz(-0.56791373) q[0];
rz(-pi) q[1];
rz(2.7841714) q[2];
sx q[2];
rz(-1.7454141) q[2];
sx q[2];
rz(-2.0703966) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8021009) q[1];
sx q[1];
rz(-1.9994945) q[1];
sx q[1];
rz(2.2080253) q[1];
x q[2];
rz(2.8350713) q[3];
sx q[3];
rz(-1.74311) q[3];
sx q[3];
rz(2.4337208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.4404099) q[2];
sx q[2];
rz(-1.5298693) q[2];
sx q[2];
rz(2.8033076) q[2];
rz(1.7017378) q[3];
sx q[3];
rz(-2.2262636) q[3];
sx q[3];
rz(2.2556944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97025362) q[0];
sx q[0];
rz(-0.71115029) q[0];
sx q[0];
rz(-0.030348226) q[0];
rz(-0.066210315) q[1];
sx q[1];
rz(-2.1538484) q[1];
sx q[1];
rz(-1.5240086) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0328007) q[0];
sx q[0];
rz(-1.5691225) q[0];
sx q[0];
rz(1.770442) q[0];
rz(-pi) q[1];
rz(1.5288058) q[2];
sx q[2];
rz(-2.6868372) q[2];
sx q[2];
rz(3.0595879) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0028573) q[1];
sx q[1];
rz(-0.5746952) q[1];
sx q[1];
rz(-1.6981305) q[1];
rz(-pi) q[2];
rz(0.094035427) q[3];
sx q[3];
rz(-1.7574851) q[3];
sx q[3];
rz(1.3446913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7559738) q[2];
sx q[2];
rz(-0.9884584) q[2];
sx q[2];
rz(1.9937817) q[2];
rz(-1.3267481) q[3];
sx q[3];
rz(-1.8170522) q[3];
sx q[3];
rz(2.9045048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0148934) q[0];
sx q[0];
rz(-0.48148695) q[0];
sx q[0];
rz(2.8258064) q[0];
rz(-2.2029927) q[1];
sx q[1];
rz(-1.4626075) q[1];
sx q[1];
rz(-0.25207239) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85900753) q[0];
sx q[0];
rz(-0.95882817) q[0];
sx q[0];
rz(0.99143272) q[0];
rz(-pi) q[1];
x q[1];
rz(0.75906934) q[2];
sx q[2];
rz(-1.0319064) q[2];
sx q[2];
rz(-2.3784504) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.38110106) q[1];
sx q[1];
rz(-0.64844202) q[1];
sx q[1];
rz(1.2566503) q[1];
rz(-2.4113703) q[3];
sx q[3];
rz(-0.50958868) q[3];
sx q[3];
rz(2.2024221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0198274) q[2];
sx q[2];
rz(-0.44744197) q[2];
sx q[2];
rz(-3.1075409) q[2];
rz(-0.017459067) q[3];
sx q[3];
rz(-1.7826467) q[3];
sx q[3];
rz(-2.0461369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
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
rz(1.0871672) q[1];
sx q[1];
rz(-0.68030578) q[1];
sx q[1];
rz(2.4345051) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0079460572) q[0];
sx q[0];
rz(-2.5292853) q[0];
sx q[0];
rz(-0.72511073) q[0];
rz(-pi) q[1];
rz(-0.44666501) q[2];
sx q[2];
rz(-1.6840877) q[2];
sx q[2];
rz(-2.5926673) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.45338079) q[1];
sx q[1];
rz(-2.4773295) q[1];
sx q[1];
rz(-1.3761671) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1726923) q[3];
sx q[3];
rz(-2.7728191) q[3];
sx q[3];
rz(-2.8440648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.84918555) q[2];
sx q[2];
rz(-1.1881928) q[2];
sx q[2];
rz(-0.46009955) q[2];
rz(1.397331) q[3];
sx q[3];
rz(-1.5887235) q[3];
sx q[3];
rz(-0.24266711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3209155) q[0];
sx q[0];
rz(-0.21235947) q[0];
sx q[0];
rz(-1.7472349) q[0];
rz(-1.0955411) q[1];
sx q[1];
rz(-1.6004326) q[1];
sx q[1];
rz(2.8869693) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8481962) q[0];
sx q[0];
rz(-0.080349803) q[0];
sx q[0];
rz(1.7424165) q[0];
x q[1];
rz(-0.014572797) q[2];
sx q[2];
rz(-2.1483148) q[2];
sx q[2];
rz(-2.2968959) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9533206) q[1];
sx q[1];
rz(-0.70536648) q[1];
sx q[1];
rz(-0.377368) q[1];
rz(-1.0765431) q[3];
sx q[3];
rz(-0.77270618) q[3];
sx q[3];
rz(2.5172174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1191117) q[2];
sx q[2];
rz(-0.20038651) q[2];
sx q[2];
rz(1.3767892) q[2];
rz(1.4962176) q[3];
sx q[3];
rz(-1.508537) q[3];
sx q[3];
rz(-2.0549324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6901533) q[0];
sx q[0];
rz(-2.4017161) q[0];
sx q[0];
rz(2.8421463) q[0];
rz(-1.0401789) q[1];
sx q[1];
rz(-1.6957915) q[1];
sx q[1];
rz(0.20656955) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20139192) q[0];
sx q[0];
rz(-0.79027806) q[0];
sx q[0];
rz(-1.9076365) q[0];
rz(0.95025392) q[2];
sx q[2];
rz(-0.62180078) q[2];
sx q[2];
rz(-1.7813462) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8787074) q[1];
sx q[1];
rz(-2.4135114) q[1];
sx q[1];
rz(1.3036149) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9566831) q[3];
sx q[3];
rz(-0.69291249) q[3];
sx q[3];
rz(0.087156765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.841659) q[2];
sx q[2];
rz(-1.724023) q[2];
sx q[2];
rz(-0.28277961) q[2];
rz(0.81280604) q[3];
sx q[3];
rz(-2.7225284) q[3];
sx q[3];
rz(-2.9747484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(0.92641002) q[0];
sx q[0];
rz(-2.0880501) q[0];
sx q[0];
rz(0.38152951) q[0];
rz(0.58386699) q[1];
sx q[1];
rz(-0.54324141) q[1];
sx q[1];
rz(-1.8136224) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0921811) q[0];
sx q[0];
rz(-1.4565399) q[0];
sx q[0];
rz(1.1294424) q[0];
x q[1];
rz(-0.46220025) q[2];
sx q[2];
rz(-1.000058) q[2];
sx q[2];
rz(-2.5525023) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2139637) q[1];
sx q[1];
rz(-2.2447526) q[1];
sx q[1];
rz(2.9272635) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4626571) q[3];
sx q[3];
rz(-0.62641615) q[3];
sx q[3];
rz(-0.39144799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.61838377) q[2];
sx q[2];
rz(-2.0760459) q[2];
sx q[2];
rz(1.7810812) q[2];
rz(1.7112188) q[3];
sx q[3];
rz(-2.1332707) q[3];
sx q[3];
rz(0.84806228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9946063) q[0];
sx q[0];
rz(-1.9817579) q[0];
sx q[0];
rz(-0.18187901) q[0];
rz(-0.47422844) q[1];
sx q[1];
rz(-1.0206181) q[1];
sx q[1];
rz(-2.1906733) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5247105) q[0];
sx q[0];
rz(-1.8572154) q[0];
sx q[0];
rz(-2.8911203) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0792564) q[2];
sx q[2];
rz(-1.8364292) q[2];
sx q[2];
rz(-0.072364256) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.37590313) q[1];
sx q[1];
rz(-1.6357058) q[1];
sx q[1];
rz(0.20903559) q[1];
rz(0.37787921) q[3];
sx q[3];
rz(-1.9482908) q[3];
sx q[3];
rz(-2.8187403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9178847) q[2];
sx q[2];
rz(-2.6612838) q[2];
sx q[2];
rz(3.0656832) q[2];
rz(-0.54801303) q[3];
sx q[3];
rz(-1.8173822) q[3];
sx q[3];
rz(2.5089335) q[3];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52371812) q[0];
sx q[0];
rz(-2.0848367) q[0];
sx q[0];
rz(1.7653718) q[0];
rz(-2.7245522) q[1];
sx q[1];
rz(-1.7224256) q[1];
sx q[1];
rz(2.4818647) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8575681) q[0];
sx q[0];
rz(-2.095247) q[0];
sx q[0];
rz(-2.2542623) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9218036) q[2];
sx q[2];
rz(-2.347749) q[2];
sx q[2];
rz(1.1500051) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8068741) q[1];
sx q[1];
rz(-1.8611307) q[1];
sx q[1];
rz(-0.90805407) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.057007313) q[3];
sx q[3];
rz(-1.032864) q[3];
sx q[3];
rz(1.2008592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.62347162) q[2];
sx q[2];
rz(-2.3770964) q[2];
sx q[2];
rz(2.1155604) q[2];
rz(-0.14885151) q[3];
sx q[3];
rz(-2.1089349) q[3];
sx q[3];
rz(0.025645105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24348564) q[0];
sx q[0];
rz(-0.74111104) q[0];
sx q[0];
rz(1.3056668) q[0];
rz(1.9650412) q[1];
sx q[1];
rz(-1.2780317) q[1];
sx q[1];
rz(-2.1059039) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0730597) q[0];
sx q[0];
rz(-1.8951299) q[0];
sx q[0];
rz(0.43138327) q[0];
x q[1];
rz(-0.91949384) q[2];
sx q[2];
rz(-2.3323625) q[2];
sx q[2];
rz(-2.2724255) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.021691572) q[1];
sx q[1];
rz(-0.38858116) q[1];
sx q[1];
rz(-0.67967023) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3922937) q[3];
sx q[3];
rz(-1.370508) q[3];
sx q[3];
rz(1.1009969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5130561) q[2];
sx q[2];
rz(-2.3302902) q[2];
sx q[2];
rz(-1.9899842) q[2];
rz(-0.51268762) q[3];
sx q[3];
rz(-1.0980462) q[3];
sx q[3];
rz(0.36469665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7619027) q[0];
sx q[0];
rz(-1.7871478) q[0];
sx q[0];
rz(-2.3085069) q[0];
rz(-1.5079386) q[1];
sx q[1];
rz(-2.5588551) q[1];
sx q[1];
rz(-0.48164639) q[1];
rz(-2.7307636) q[2];
sx q[2];
rz(-1.6098235) q[2];
sx q[2];
rz(3.1237684) q[2];
rz(1.5619754) q[3];
sx q[3];
rz(-2.1891441) q[3];
sx q[3];
rz(1.7563663) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
