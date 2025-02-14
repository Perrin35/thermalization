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
rz(-2.6920707) q[0];
sx q[0];
rz(-1.7188526) q[0];
sx q[0];
rz(2.0916405) q[0];
rz(2.6226251) q[1];
sx q[1];
rz(-1.608404) q[1];
sx q[1];
rz(3.0906711) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1919132) q[0];
sx q[0];
rz(-2.5891719) q[0];
sx q[0];
rz(-0.64903736) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.21587431) q[2];
sx q[2];
rz(-1.2435438) q[2];
sx q[2];
rz(0.33282166) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4733153) q[1];
sx q[1];
rz(-0.44679579) q[1];
sx q[1];
rz(2.3832641) q[1];
x q[2];
rz(-2.256278) q[3];
sx q[3];
rz(-1.3847794) q[3];
sx q[3];
rz(1.7486339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8672987) q[2];
sx q[2];
rz(-2.2653502) q[2];
sx q[2];
rz(-1.743861) q[2];
rz(2.6869669) q[3];
sx q[3];
rz(-0.8780829) q[3];
sx q[3];
rz(1.4286058) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4947263) q[0];
sx q[0];
rz(-0.88748256) q[0];
sx q[0];
rz(2.5383762) q[0];
rz(-1.2546722) q[1];
sx q[1];
rz(-2.5277977) q[1];
sx q[1];
rz(2.5616554) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7042687) q[0];
sx q[0];
rz(-0.56341972) q[0];
sx q[0];
rz(-0.31524865) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5712319) q[2];
sx q[2];
rz(-1.1820536) q[2];
sx q[2];
rz(1.1136356) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1973639) q[1];
sx q[1];
rz(-1.6104638) q[1];
sx q[1];
rz(-0.96645379) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5744561) q[3];
sx q[3];
rz(-1.2492164) q[3];
sx q[3];
rz(-1.086832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6569528) q[2];
sx q[2];
rz(-2.3886949) q[2];
sx q[2];
rz(-2.7638655) q[2];
rz(1.4332625) q[3];
sx q[3];
rz(-1.5092311) q[3];
sx q[3];
rz(-1.1908971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3483873) q[0];
sx q[0];
rz(-0.054940104) q[0];
sx q[0];
rz(1.4980263) q[0];
rz(1.161423) q[1];
sx q[1];
rz(-1.8893416) q[1];
sx q[1];
rz(0.84322554) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6093685) q[0];
sx q[0];
rz(-0.62590963) q[0];
sx q[0];
rz(1.468991) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.21768985) q[2];
sx q[2];
rz(-1.9799616) q[2];
sx q[2];
rz(2.7646551) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9595782) q[1];
sx q[1];
rz(-1.6977662) q[1];
sx q[1];
rz(0.89163973) q[1];
x q[2];
rz(0.39210658) q[3];
sx q[3];
rz(-0.74015731) q[3];
sx q[3];
rz(1.8405434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2735542) q[2];
sx q[2];
rz(-2.4783583) q[2];
sx q[2];
rz(-0.794945) q[2];
rz(0.76977175) q[3];
sx q[3];
rz(-0.92307463) q[3];
sx q[3];
rz(2.9845089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.58761) q[0];
sx q[0];
rz(-1.8164604) q[0];
sx q[0];
rz(-2.8532568) q[0];
rz(-0.68710697) q[1];
sx q[1];
rz(-1.4930875) q[1];
sx q[1];
rz(-1.4289325) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2401348) q[0];
sx q[0];
rz(-3.1278538) q[0];
sx q[0];
rz(0.49679784) q[0];
rz(-pi) q[1];
rz(0.40224125) q[2];
sx q[2];
rz(-1.4415485) q[2];
sx q[2];
rz(2.8856173) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1667119) q[1];
sx q[1];
rz(-0.26435164) q[1];
sx q[1];
rz(-2.8893785) q[1];
rz(-2.3219548) q[3];
sx q[3];
rz(-0.82806168) q[3];
sx q[3];
rz(2.6520906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5700506) q[2];
sx q[2];
rz(-2.1752581) q[2];
sx q[2];
rz(-0.16560444) q[2];
rz(-3.0770732) q[3];
sx q[3];
rz(-0.12972984) q[3];
sx q[3];
rz(1.2714413) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22873838) q[0];
sx q[0];
rz(-1.9485291) q[0];
sx q[0];
rz(2.2304529) q[0];
rz(0.97081026) q[1];
sx q[1];
rz(-1.3263005) q[1];
sx q[1];
rz(-2.2784746) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9476833) q[0];
sx q[0];
rz(-1.2811617) q[0];
sx q[0];
rz(-1.0632443) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4865033) q[2];
sx q[2];
rz(-2.8721923) q[2];
sx q[2];
rz(2.0790554) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.48820254) q[1];
sx q[1];
rz(-2.0987256) q[1];
sx q[1];
rz(1.3009682) q[1];
rz(-pi) q[2];
rz(1.7105402) q[3];
sx q[3];
rz(-0.31523963) q[3];
sx q[3];
rz(-1.3073352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0866278) q[2];
sx q[2];
rz(-1.2530155) q[2];
sx q[2];
rz(0.48913726) q[2];
rz(-0.9984115) q[3];
sx q[3];
rz(-2.9676134) q[3];
sx q[3];
rz(-3.0424931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84422207) q[0];
sx q[0];
rz(-2.4956644) q[0];
sx q[0];
rz(0.55322629) q[0];
rz(-0.79752254) q[1];
sx q[1];
rz(-0.99634606) q[1];
sx q[1];
rz(-1.1513938) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40276179) q[0];
sx q[0];
rz(-0.58332764) q[0];
sx q[0];
rz(2.8999694) q[0];
rz(-1.3998447) q[2];
sx q[2];
rz(-2.1151849) q[2];
sx q[2];
rz(0.65786568) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7094671) q[1];
sx q[1];
rz(-1.3015987) q[1];
sx q[1];
rz(1.7749191) q[1];
rz(-pi) q[2];
rz(1.0955515) q[3];
sx q[3];
rz(-2.4184347) q[3];
sx q[3];
rz(2.1840546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3449751) q[2];
sx q[2];
rz(-1.7285708) q[2];
sx q[2];
rz(-0.83149347) q[2];
rz(1.3384532) q[3];
sx q[3];
rz(-1.0748539) q[3];
sx q[3];
rz(-2.2510546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3727342) q[0];
sx q[0];
rz(-1.2059728) q[0];
sx q[0];
rz(0.15705577) q[0];
rz(1.318469) q[1];
sx q[1];
rz(-0.942197) q[1];
sx q[1];
rz(2.7838321) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8028594) q[0];
sx q[0];
rz(-2.0293616) q[0];
sx q[0];
rz(-1.073097) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.90334185) q[2];
sx q[2];
rz(-1.2111665) q[2];
sx q[2];
rz(2.0697167) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.58531144) q[1];
sx q[1];
rz(-2.1211294) q[1];
sx q[1];
rz(0.78720168) q[1];
x q[2];
rz(-2.5671183) q[3];
sx q[3];
rz(-1.7041429) q[3];
sx q[3];
rz(0.80783081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1498108) q[2];
sx q[2];
rz(-1.2573743) q[2];
sx q[2];
rz(3.120976) q[2];
rz(1.8990382) q[3];
sx q[3];
rz(-0.81762448) q[3];
sx q[3];
rz(-2.8500565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6108516) q[0];
sx q[0];
rz(-1.1854956) q[0];
sx q[0];
rz(0.32671842) q[0];
rz(-0.84699455) q[1];
sx q[1];
rz(-2.0920483) q[1];
sx q[1];
rz(-2.0832031) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.055870067) q[0];
sx q[0];
rz(-2.6230109) q[0];
sx q[0];
rz(-1.7118171) q[0];
x q[1];
rz(-2.5024274) q[2];
sx q[2];
rz(-1.7129538) q[2];
sx q[2];
rz(-0.41691142) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2742334) q[1];
sx q[1];
rz(-1.6489442) q[1];
sx q[1];
rz(-1.9539425) q[1];
rz(-pi) q[2];
rz(-1.0048546) q[3];
sx q[3];
rz(-1.3019239) q[3];
sx q[3];
rz(-3.1389126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0224926) q[2];
sx q[2];
rz(-2.7969226) q[2];
sx q[2];
rz(-1.9692839) q[2];
rz(-3.1398224) q[3];
sx q[3];
rz(-0.5032731) q[3];
sx q[3];
rz(0.9790023) q[3];
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
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2718647) q[0];
sx q[0];
rz(-0.80369049) q[0];
sx q[0];
rz(-1.3386238) q[0];
rz(1.9610693) q[1];
sx q[1];
rz(-2.0484643) q[1];
sx q[1];
rz(2.6611633) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.521489) q[0];
sx q[0];
rz(-2.0248981) q[0];
sx q[0];
rz(-0.2257077) q[0];
rz(-2.9715638) q[2];
sx q[2];
rz(-1.5289834) q[2];
sx q[2];
rz(0.60838503) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.94225223) q[1];
sx q[1];
rz(-2.128731) q[1];
sx q[1];
rz(0.11726484) q[1];
rz(-1.3442743) q[3];
sx q[3];
rz(-1.220928) q[3];
sx q[3];
rz(2.1893152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7052445) q[2];
sx q[2];
rz(-2.1465325) q[2];
sx q[2];
rz(-2.056541) q[2];
rz(-2.0294225) q[3];
sx q[3];
rz(-2.2472436) q[3];
sx q[3];
rz(2.3280242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28806624) q[0];
sx q[0];
rz(-2.1948094) q[0];
sx q[0];
rz(2.1296401) q[0];
rz(2.2400253) q[1];
sx q[1];
rz(-2.8755867) q[1];
sx q[1];
rz(2.576135) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5150499) q[0];
sx q[0];
rz(-1.7789655) q[0];
sx q[0];
rz(2.2295537) q[0];
rz(1.9996634) q[2];
sx q[2];
rz(-1.326401) q[2];
sx q[2];
rz(1.7606869) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1484851) q[1];
sx q[1];
rz(-2.4186181) q[1];
sx q[1];
rz(-1.7112247) q[1];
rz(-pi) q[2];
rz(1.2540482) q[3];
sx q[3];
rz(-1.4707139) q[3];
sx q[3];
rz(-2.7730178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3507583) q[2];
sx q[2];
rz(-1.5956722) q[2];
sx q[2];
rz(-1.8809543) q[2];
rz(-2.6993921) q[3];
sx q[3];
rz(-2.111777) q[3];
sx q[3];
rz(1.376576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5653391) q[0];
sx q[0];
rz(-2.2591142) q[0];
sx q[0];
rz(2.2478065) q[0];
rz(1.5671989) q[1];
sx q[1];
rz(-1.6849453) q[1];
sx q[1];
rz(-1.4243855) q[1];
rz(1.8509751) q[2];
sx q[2];
rz(-1.1600672) q[2];
sx q[2];
rz(1.5008012) q[2];
rz(1.9954197) q[3];
sx q[3];
rz(-0.38214798) q[3];
sx q[3];
rz(1.8761763) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
