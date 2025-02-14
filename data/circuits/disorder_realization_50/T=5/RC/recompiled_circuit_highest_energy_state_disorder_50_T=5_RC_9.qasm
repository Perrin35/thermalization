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
rz(-2.4432776) q[0];
sx q[0];
rz(-2.7628216) q[0];
sx q[0];
rz(1.9842499) q[0];
rz(0.78504374) q[1];
sx q[1];
rz(-0.68612376) q[1];
sx q[1];
rz(2.6148028) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.049126712) q[0];
sx q[0];
rz(-1.2727588) q[0];
sx q[0];
rz(-0.98182337) q[0];
rz(-0.85028591) q[2];
sx q[2];
rz(-0.97480259) q[2];
sx q[2];
rz(1.13499) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.59299027) q[1];
sx q[1];
rz(-2.0128662) q[1];
sx q[1];
rz(1.6111239) q[1];
x q[2];
rz(-3.111638) q[3];
sx q[3];
rz(-2.1308793) q[3];
sx q[3];
rz(-0.64489472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.36246768) q[2];
sx q[2];
rz(-2.1802826) q[2];
sx q[2];
rz(-2.989952) q[2];
rz(-2.9996297) q[3];
sx q[3];
rz(-1.6715334) q[3];
sx q[3];
rz(2.0040472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.66459429) q[0];
sx q[0];
rz(-0.73200309) q[0];
sx q[0];
rz(-0.81749302) q[0];
rz(3.0290161) q[1];
sx q[1];
rz(-1.45603) q[1];
sx q[1];
rz(-2.2419825) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0257638) q[0];
sx q[0];
rz(-1.8432105) q[0];
sx q[0];
rz(1.8194356) q[0];
x q[1];
rz(-2.4075872) q[2];
sx q[2];
rz(-0.81799928) q[2];
sx q[2];
rz(0.98246511) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1020722) q[1];
sx q[1];
rz(-1.5145434) q[1];
sx q[1];
rz(-0.7489167) q[1];
x q[2];
rz(2.2475776) q[3];
sx q[3];
rz(-1.1543659) q[3];
sx q[3];
rz(-0.82899603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.67109913) q[2];
sx q[2];
rz(-1.495139) q[2];
sx q[2];
rz(-2.1136843) q[2];
rz(-0.73728621) q[3];
sx q[3];
rz(-1.6643915) q[3];
sx q[3];
rz(-3.1173053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-2.1360433) q[0];
sx q[0];
rz(-0.5558973) q[0];
sx q[0];
rz(-1.1935724) q[0];
rz(1.8434803) q[1];
sx q[1];
rz(-1.6622512) q[1];
sx q[1];
rz(1.6568291) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.68736) q[0];
sx q[0];
rz(-3.0890565) q[0];
sx q[0];
rz(-2.0665069) q[0];
rz(-pi) q[1];
rz(3.0151691) q[2];
sx q[2];
rz(-1.7806541) q[2];
sx q[2];
rz(-2.1239547) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.72276593) q[1];
sx q[1];
rz(-1.9120941) q[1];
sx q[1];
rz(-1.7143634) q[1];
rz(1.9336352) q[3];
sx q[3];
rz(-1.6479744) q[3];
sx q[3];
rz(0.046260351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9803311) q[2];
sx q[2];
rz(-1.9183466) q[2];
sx q[2];
rz(0.70183357) q[2];
rz(3.0724604) q[3];
sx q[3];
rz(-0.88874236) q[3];
sx q[3];
rz(-0.70772901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.055701) q[0];
sx q[0];
rz(-1.1393071) q[0];
sx q[0];
rz(0.62537801) q[0];
rz(-1.0944132) q[1];
sx q[1];
rz(-1.5233327) q[1];
sx q[1];
rz(-1.3166924) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65526456) q[0];
sx q[0];
rz(-1.0287971) q[0];
sx q[0];
rz(2.9744451) q[0];
rz(-pi) q[1];
rz(-1.6676297) q[2];
sx q[2];
rz(-1.2487354) q[2];
sx q[2];
rz(-1.9063973) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3299371) q[1];
sx q[1];
rz(-1.3906758) q[1];
sx q[1];
rz(-0.7094769) q[1];
x q[2];
rz(2.5685124) q[3];
sx q[3];
rz(-2.2551558) q[3];
sx q[3];
rz(2.7671368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7350498) q[2];
sx q[2];
rz(-2.4989765) q[2];
sx q[2];
rz(-3.1114846) q[2];
rz(2.7249469) q[3];
sx q[3];
rz(-1.7130339) q[3];
sx q[3];
rz(3.0863975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77867126) q[0];
sx q[0];
rz(-1.2782949) q[0];
sx q[0];
rz(-1.2177421) q[0];
rz(-2.9171004) q[1];
sx q[1];
rz(-1.4935363) q[1];
sx q[1];
rz(2.7405558) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.359085) q[0];
sx q[0];
rz(-1.9422724) q[0];
sx q[0];
rz(1.1228193) q[0];
x q[1];
rz(-1.941576) q[2];
sx q[2];
rz(-2.1919498) q[2];
sx q[2];
rz(-2.7948684) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2485473) q[1];
sx q[1];
rz(-1.5233938) q[1];
sx q[1];
rz(0.23706146) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4698074) q[3];
sx q[3];
rz(-1.7324395) q[3];
sx q[3];
rz(-2.9532331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2991221) q[2];
sx q[2];
rz(-1.2229908) q[2];
sx q[2];
rz(-2.6341338) q[2];
rz(-2.2211645) q[3];
sx q[3];
rz(-0.43693742) q[3];
sx q[3];
rz(0.18032716) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64666635) q[0];
sx q[0];
rz(-0.05519069) q[0];
sx q[0];
rz(-2.1221509) q[0];
rz(1.9693718) q[1];
sx q[1];
rz(-1.1727138) q[1];
sx q[1];
rz(1.2219465) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2766314) q[0];
sx q[0];
rz(-0.91827938) q[0];
sx q[0];
rz(0.44256532) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.055962) q[2];
sx q[2];
rz(-1.5400572) q[2];
sx q[2];
rz(1.5167502) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1091369) q[1];
sx q[1];
rz(-1.0278406) q[1];
sx q[1];
rz(1.2486267) q[1];
rz(1.1282519) q[3];
sx q[3];
rz(-2.5199515) q[3];
sx q[3];
rz(0.78237247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.4814066) q[2];
sx q[2];
rz(-0.6654827) q[2];
sx q[2];
rz(-2.0484203) q[2];
rz(-0.53705755) q[3];
sx q[3];
rz(-1.4635181) q[3];
sx q[3];
rz(2.4721036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90401232) q[0];
sx q[0];
rz(-2.8222988) q[0];
sx q[0];
rz(-1.4599266) q[0];
rz(1.5096674) q[1];
sx q[1];
rz(-2.5131112) q[1];
sx q[1];
rz(-0.07930886) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8312165) q[0];
sx q[0];
rz(-2.8908911) q[0];
sx q[0];
rz(-1.901907) q[0];
x q[1];
rz(-2.5517188) q[2];
sx q[2];
rz(-0.24352077) q[2];
sx q[2];
rz(-0.68984725) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.39187688) q[1];
sx q[1];
rz(-1.2971216) q[1];
sx q[1];
rz(2.7065212) q[1];
rz(-pi) q[2];
rz(2.4238911) q[3];
sx q[3];
rz(-1.6738179) q[3];
sx q[3];
rz(2.3915402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1377533) q[2];
sx q[2];
rz(-1.1649818) q[2];
sx q[2];
rz(0.4129146) q[2];
rz(1.2952992) q[3];
sx q[3];
rz(-0.92419878) q[3];
sx q[3];
rz(2.7220461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9846648) q[0];
sx q[0];
rz(-2.4913737) q[0];
sx q[0];
rz(0.28537634) q[0];
rz(-3.0889619) q[1];
sx q[1];
rz(-1.4080518) q[1];
sx q[1];
rz(0.1836798) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7516458) q[0];
sx q[0];
rz(-0.11706287) q[0];
sx q[0];
rz(-1.879891) q[0];
rz(-pi) q[1];
rz(-1.7247612) q[2];
sx q[2];
rz(-1.2885258) q[2];
sx q[2];
rz(2.8090614) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.86183) q[1];
sx q[1];
rz(-0.97721902) q[1];
sx q[1];
rz(1.66278) q[1];
rz(-pi) q[2];
rz(-0.37723549) q[3];
sx q[3];
rz(-0.68359112) q[3];
sx q[3];
rz(-2.205276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.24667428) q[2];
sx q[2];
rz(-1.7343212) q[2];
sx q[2];
rz(2.0764009) q[2];
rz(-1.4554321) q[3];
sx q[3];
rz(-1.8543517) q[3];
sx q[3];
rz(-0.58568946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10373779) q[0];
sx q[0];
rz(-2.3268564) q[0];
sx q[0];
rz(1.6023585) q[0];
rz(-2.1922951) q[1];
sx q[1];
rz(-2.0930591) q[1];
sx q[1];
rz(-1.4303713) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.024139013) q[0];
sx q[0];
rz(-1.1012337) q[0];
sx q[0];
rz(-0.59654327) q[0];
rz(-2.2995641) q[2];
sx q[2];
rz(-1.5146189) q[2];
sx q[2];
rz(-2.2420355) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.51586823) q[1];
sx q[1];
rz(-1.6004381) q[1];
sx q[1];
rz(-0.29582204) q[1];
rz(-pi) q[2];
rz(-0.28258459) q[3];
sx q[3];
rz(-1.9842098) q[3];
sx q[3];
rz(1.9081209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1324233) q[2];
sx q[2];
rz(-1.8370266) q[2];
sx q[2];
rz(-3.0150748) q[2];
rz(-2.8070519) q[3];
sx q[3];
rz(-0.25944513) q[3];
sx q[3];
rz(-0.87219605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8406521) q[0];
sx q[0];
rz(-0.88467389) q[0];
sx q[0];
rz(0.51710039) q[0];
rz(-3.0866947) q[1];
sx q[1];
rz(-1.5093191) q[1];
sx q[1];
rz(3.0518234) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15799284) q[0];
sx q[0];
rz(-1.1674034) q[0];
sx q[0];
rz(-0.97121111) q[0];
rz(-pi) q[1];
rz(2.1316809) q[2];
sx q[2];
rz(-1.5176306) q[2];
sx q[2];
rz(-1.9694984) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.59261428) q[1];
sx q[1];
rz(-1.4121118) q[1];
sx q[1];
rz(0.25513809) q[1];
rz(-pi) q[2];
rz(2.91342) q[3];
sx q[3];
rz(-1.2163278) q[3];
sx q[3];
rz(1.4972403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.53139293) q[2];
sx q[2];
rz(-1.1343845) q[2];
sx q[2];
rz(2.4143977) q[2];
rz(-2.3860892) q[3];
sx q[3];
rz(-0.38755363) q[3];
sx q[3];
rz(2.9288647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94201921) q[0];
sx q[0];
rz(-1.6006391) q[0];
sx q[0];
rz(-2.0508456) q[0];
rz(-1.749281) q[1];
sx q[1];
rz(-1.6284457) q[1];
sx q[1];
rz(1.5096691) q[1];
rz(1.4452151) q[2];
sx q[2];
rz(-2.4332252) q[2];
sx q[2];
rz(-0.38548584) q[2];
rz(-0.94115067) q[3];
sx q[3];
rz(-2.0654021) q[3];
sx q[3];
rz(0.13765814) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
