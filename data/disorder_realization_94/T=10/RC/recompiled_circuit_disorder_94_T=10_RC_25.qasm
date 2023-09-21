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
rz(-0.59626055) q[1];
sx q[1];
rz(-2.526386) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0628375) q[0];
sx q[0];
rz(-0.7987928) q[0];
sx q[0];
rz(0.90057217) q[0];
x q[1];
rz(1.7569321) q[2];
sx q[2];
rz(-1.2190483) q[2];
sx q[2];
rz(0.56439161) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7429744) q[1];
sx q[1];
rz(-2.3906039) q[1];
sx q[1];
rz(-2.2258334) q[1];
x q[2];
rz(0.52317046) q[3];
sx q[3];
rz(-0.35029951) q[3];
sx q[3];
rz(-1.7821799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7011828) q[2];
sx q[2];
rz(-1.6117233) q[2];
sx q[2];
rz(2.8033076) q[2];
rz(1.7017378) q[3];
sx q[3];
rz(-2.2262636) q[3];
sx q[3];
rz(-0.88589823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.171339) q[0];
sx q[0];
rz(-2.4304424) q[0];
sx q[0];
rz(-0.030348226) q[0];
rz(0.066210315) q[1];
sx q[1];
rz(-2.1538484) q[1];
sx q[1];
rz(-1.617584) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.671316) q[0];
sx q[0];
rz(-0.19965262) q[0];
sx q[0];
rz(1.5623564) q[0];
rz(1.5288058) q[2];
sx q[2];
rz(-0.4547555) q[2];
sx q[2];
rz(0.08200478) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6807032) q[1];
sx q[1];
rz(-1.5017121) q[1];
sx q[1];
rz(0.99980385) q[1];
x q[2];
rz(2.0321235) q[3];
sx q[3];
rz(-0.20878775) q[3];
sx q[3];
rz(1.8148282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7559738) q[2];
sx q[2];
rz(-0.9884584) q[2];
sx q[2];
rz(1.9937817) q[2];
rz(1.8148445) q[3];
sx q[3];
rz(-1.3245405) q[3];
sx q[3];
rz(-2.9045048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1266992) q[0];
sx q[0];
rz(-0.48148695) q[0];
sx q[0];
rz(0.31578627) q[0];
rz(-0.93859998) q[1];
sx q[1];
rz(-1.4626075) q[1];
sx q[1];
rz(0.25207239) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7092428) q[0];
sx q[0];
rz(-2.3253257) q[0];
sx q[0];
rz(2.4791251) q[0];
rz(-pi) q[1];
rz(2.2601068) q[2];
sx q[2];
rz(-2.2027317) q[2];
sx q[2];
rz(-1.2607247) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4430868) q[1];
sx q[1];
rz(-1.7585187) q[1];
sx q[1];
rz(2.1952573) q[1];
rz(-1.2139981) q[3];
sx q[3];
rz(-1.1988415) q[3];
sx q[3];
rz(-1.4043722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1217653) q[2];
sx q[2];
rz(-2.6941507) q[2];
sx q[2];
rz(0.034051731) q[2];
rz(0.017459067) q[3];
sx q[3];
rz(-1.7826467) q[3];
sx q[3];
rz(2.0461369) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5144192) q[0];
sx q[0];
rz(-1.592941) q[0];
sx q[0];
rz(1.6148286) q[0];
rz(1.0871672) q[1];
sx q[1];
rz(-0.68030578) q[1];
sx q[1];
rz(-0.70708752) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93543816) q[0];
sx q[0];
rz(-1.1797138) q[0];
sx q[0];
rz(-2.6576256) q[0];
rz(-pi) q[1];
rz(2.6949276) q[2];
sx q[2];
rz(-1.6840877) q[2];
sx q[2];
rz(-2.5926673) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1781588) q[1];
sx q[1];
rz(-1.4512832) q[1];
sx q[1];
rz(0.91576373) q[1];
x q[2];
rz(2.9261758) q[3];
sx q[3];
rz(-1.2691174) q[3];
sx q[3];
rz(-0.33723436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.84918555) q[2];
sx q[2];
rz(-1.1881928) q[2];
sx q[2];
rz(2.6814931) q[2];
rz(-1.7442616) q[3];
sx q[3];
rz(-1.5887235) q[3];
sx q[3];
rz(2.8989255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3209155) q[0];
sx q[0];
rz(-2.9292332) q[0];
sx q[0];
rz(-1.3943577) q[0];
rz(-2.0460515) q[1];
sx q[1];
rz(-1.54116) q[1];
sx q[1];
rz(-0.25462338) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8481962) q[0];
sx q[0];
rz(-0.080349803) q[0];
sx q[0];
rz(-1.7424165) q[0];
rz(-pi) q[1];
rz(-1.5484372) q[2];
sx q[2];
rz(-0.5776814) q[2];
sx q[2];
rz(-0.81800848) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8498807) q[1];
sx q[1];
rz(-2.2176718) q[1];
sx q[1];
rz(-1.8748267) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0765431) q[3];
sx q[3];
rz(-0.77270618) q[3];
sx q[3];
rz(-0.62437526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.022481) q[2];
sx q[2];
rz(-0.20038651) q[2];
sx q[2];
rz(1.3767892) q[2];
rz(-1.6453751) q[3];
sx q[3];
rz(-1.6330556) q[3];
sx q[3];
rz(-1.0866603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.45143932) q[0];
sx q[0];
rz(-0.7398766) q[0];
sx q[0];
rz(-0.29944637) q[0];
rz(-2.1014138) q[1];
sx q[1];
rz(-1.4458011) q[1];
sx q[1];
rz(-2.9350231) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9402007) q[0];
sx q[0];
rz(-2.3513146) q[0];
sx q[0];
rz(1.9076365) q[0];
x q[1];
rz(-2.7467696) q[2];
sx q[2];
rz(-1.0770814) q[2];
sx q[2];
rz(2.0815108) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0530015) q[1];
sx q[1];
rz(-2.2676761) q[1];
sx q[1];
rz(-0.23115302) q[1];
x q[2];
rz(-0.30287403) q[3];
sx q[3];
rz(-0.93749638) q[3];
sx q[3];
rz(2.5686222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.841659) q[2];
sx q[2];
rz(-1.724023) q[2];
sx q[2];
rz(0.28277961) q[2];
rz(-0.81280604) q[3];
sx q[3];
rz(-0.41906425) q[3];
sx q[3];
rz(0.16684428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2151826) q[0];
sx q[0];
rz(-2.0880501) q[0];
sx q[0];
rz(2.7600631) q[0];
rz(0.58386699) q[1];
sx q[1];
rz(-0.54324141) q[1];
sx q[1];
rz(-1.8136224) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5751942) q[0];
sx q[0];
rz(-1.1325206) q[0];
sx q[0];
rz(-0.12624329) q[0];
x q[1];
rz(-0.46220025) q[2];
sx q[2];
rz(-2.1415347) q[2];
sx q[2];
rz(2.5525023) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3634062) q[1];
sx q[1];
rz(-1.4038329) q[1];
sx q[1];
rz(2.2560675) q[1];
rz(1.1442723) q[3];
sx q[3];
rz(-1.0970308) q[3];
sx q[3];
rz(1.9667448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5232089) q[2];
sx q[2];
rz(-2.0760459) q[2];
sx q[2];
rz(1.7810812) q[2];
rz(1.7112188) q[3];
sx q[3];
rz(-1.0083219) q[3];
sx q[3];
rz(-0.84806228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9946063) q[0];
sx q[0];
rz(-1.1598347) q[0];
sx q[0];
rz(2.9597136) q[0];
rz(0.47422844) q[1];
sx q[1];
rz(-2.1209746) q[1];
sx q[1];
rz(-2.1906733) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2598341) q[0];
sx q[0];
rz(-1.330733) q[0];
sx q[0];
rz(1.865922) q[0];
rz(-1.7958926) q[2];
sx q[2];
rz(-2.8689119) q[2];
sx q[2];
rz(0.30579145) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2434395) q[1];
sx q[1];
rz(-2.9228518) q[1];
sx q[1];
rz(-2.8380413) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1675646) q[3];
sx q[3];
rz(-1.2207165) q[3];
sx q[3];
rz(-2.0389327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2237079) q[2];
sx q[2];
rz(-0.48030883) q[2];
sx q[2];
rz(3.0656832) q[2];
rz(0.54801303) q[3];
sx q[3];
rz(-1.3242105) q[3];
sx q[3];
rz(-0.63265911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6178745) q[0];
sx q[0];
rz(-1.0567559) q[0];
sx q[0];
rz(1.3762208) q[0];
rz(-2.7245522) q[1];
sx q[1];
rz(-1.4191671) q[1];
sx q[1];
rz(0.65972796) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8575681) q[0];
sx q[0];
rz(-1.0463456) q[0];
sx q[0];
rz(-0.88733034) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3525891) q[2];
sx q[2];
rz(-2.3404684) q[2];
sx q[2];
rz(1.6831236) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8068741) q[1];
sx q[1];
rz(-1.8611307) q[1];
sx q[1];
rz(2.2335386) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0321484) q[3];
sx q[3];
rz(-1.6197455) q[3];
sx q[3];
rz(0.34070542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.518121) q[2];
sx q[2];
rz(-2.3770964) q[2];
sx q[2];
rz(2.1155604) q[2];
rz(-2.9927411) q[3];
sx q[3];
rz(-2.1089349) q[3];
sx q[3];
rz(3.1159475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.898107) q[0];
sx q[0];
rz(-2.4004816) q[0];
sx q[0];
rz(-1.8359258) q[0];
rz(1.9650412) q[1];
sx q[1];
rz(-1.8635609) q[1];
sx q[1];
rz(2.1059039) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0730597) q[0];
sx q[0];
rz(-1.2464628) q[0];
sx q[0];
rz(0.43138327) q[0];
rz(-pi) q[1];
rz(-2.2220988) q[2];
sx q[2];
rz(-2.3323625) q[2];
sx q[2];
rz(2.2724255) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.021691572) q[1];
sx q[1];
rz(-0.38858116) q[1];
sx q[1];
rz(-2.4619224) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3003179) q[3];
sx q[3];
rz(-2.3016553) q[3];
sx q[3];
rz(-2.4887816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5130561) q[2];
sx q[2];
rz(-2.3302902) q[2];
sx q[2];
rz(1.1516085) q[2];
rz(-2.628905) q[3];
sx q[3];
rz(-2.0435464) q[3];
sx q[3];
rz(0.36469665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37968996) q[0];
sx q[0];
rz(-1.3544449) q[0];
sx q[0];
rz(0.83308573) q[0];
rz(-1.6336541) q[1];
sx q[1];
rz(-0.58273756) q[1];
sx q[1];
rz(2.6599463) q[1];
rz(-0.41082906) q[2];
sx q[2];
rz(-1.5317691) q[2];
sx q[2];
rz(-0.017824235) q[2];
rz(-0.012398331) q[3];
sx q[3];
rz(-0.61840246) q[3];
sx q[3];
rz(1.7411504) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
