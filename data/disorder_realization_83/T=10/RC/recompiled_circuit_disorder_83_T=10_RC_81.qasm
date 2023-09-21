OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5542334) q[0];
sx q[0];
rz(-2.1589307) q[0];
sx q[0];
rz(0.76173705) q[0];
rz(0.76454437) q[1];
sx q[1];
rz(4.2188797) q[1];
sx q[1];
rz(10.168434) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0640472) q[0];
sx q[0];
rz(-0.24928688) q[0];
sx q[0];
rz(1.3011978) q[0];
rz(-1.6913937) q[2];
sx q[2];
rz(-2.2199124) q[2];
sx q[2];
rz(-2.803034) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.88120645) q[1];
sx q[1];
rz(-1.2938061) q[1];
sx q[1];
rz(-1.9828412) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0697883) q[3];
sx q[3];
rz(-1.9339438) q[3];
sx q[3];
rz(2.8920409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4347697) q[2];
sx q[2];
rz(-2.7396024) q[2];
sx q[2];
rz(0.10786954) q[2];
rz(2.9989631) q[3];
sx q[3];
rz(-1.4167891) q[3];
sx q[3];
rz(0.67255783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.07664872) q[0];
sx q[0];
rz(-2.3363484) q[0];
sx q[0];
rz(2.9108677) q[0];
rz(-1.327286) q[1];
sx q[1];
rz(-2.4704411) q[1];
sx q[1];
rz(-0.040963106) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7777268) q[0];
sx q[0];
rz(-1.4085318) q[0];
sx q[0];
rz(1.8198387) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3054753) q[2];
sx q[2];
rz(-0.40514075) q[2];
sx q[2];
rz(0.9678313) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.37521024) q[1];
sx q[1];
rz(-0.76645215) q[1];
sx q[1];
rz(0.68739989) q[1];
rz(1.1702234) q[3];
sx q[3];
rz(-1.2188606) q[3];
sx q[3];
rz(-0.34917253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9791947) q[2];
sx q[2];
rz(-1.4731864) q[2];
sx q[2];
rz(-2.5734148) q[2];
rz(0.5125106) q[3];
sx q[3];
rz(-0.5957225) q[3];
sx q[3];
rz(-0.56186831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42999643) q[0];
sx q[0];
rz(-1.5189518) q[0];
sx q[0];
rz(-2.6913753) q[0];
rz(-1.2954767) q[1];
sx q[1];
rz(-1.9772915) q[1];
sx q[1];
rz(-0.67726642) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8836425) q[0];
sx q[0];
rz(-1.6227239) q[0];
sx q[0];
rz(-2.9731263) q[0];
rz(1.5696987) q[2];
sx q[2];
rz(-1.5771616) q[2];
sx q[2];
rz(0.15417834) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.22563572) q[1];
sx q[1];
rz(-1.3924053) q[1];
sx q[1];
rz(0.13326463) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1797818) q[3];
sx q[3];
rz(-2.0174332) q[3];
sx q[3];
rz(3.1004268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.21851097) q[2];
sx q[2];
rz(-1.1542902) q[2];
sx q[2];
rz(-0.046860524) q[2];
rz(-2.3299407) q[3];
sx q[3];
rz(-2.4620158) q[3];
sx q[3];
rz(2.9714382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1470404) q[0];
sx q[0];
rz(-0.11163286) q[0];
sx q[0];
rz(-3.0901093) q[0];
rz(0.4908081) q[1];
sx q[1];
rz(-2.1665116) q[1];
sx q[1];
rz(-1.9690008) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20576142) q[0];
sx q[0];
rz(-1.4371705) q[0];
sx q[0];
rz(0.66098102) q[0];
rz(-pi) q[1];
rz(-3.1191191) q[2];
sx q[2];
rz(-0.60612504) q[2];
sx q[2];
rz(1.4215353) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7395775) q[1];
sx q[1];
rz(-2.6699454) q[1];
sx q[1];
rz(2.2775843) q[1];
rz(-pi) q[2];
rz(-1.7309534) q[3];
sx q[3];
rz(-2.7048116) q[3];
sx q[3];
rz(0.76524759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.88901687) q[2];
sx q[2];
rz(-2.2106407) q[2];
sx q[2];
rz(2.6546997) q[2];
rz(-1.1934818) q[3];
sx q[3];
rz(-0.32961696) q[3];
sx q[3];
rz(-0.025432767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.837773) q[0];
sx q[0];
rz(-1.0205512) q[0];
sx q[0];
rz(1.3254962) q[0];
rz(-0.59108132) q[1];
sx q[1];
rz(-2.4609844) q[1];
sx q[1];
rz(2.4868734) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7703169) q[0];
sx q[0];
rz(-1.2219984) q[0];
sx q[0];
rz(2.6692997) q[0];
x q[1];
rz(-2.7227313) q[2];
sx q[2];
rz(-1.3992116) q[2];
sx q[2];
rz(-2.3562743) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4795585) q[1];
sx q[1];
rz(-2.2231632) q[1];
sx q[1];
rz(-2.0550834) q[1];
rz(-1.4773024) q[3];
sx q[3];
rz(-1.1870541) q[3];
sx q[3];
rz(1.5941217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6902265) q[2];
sx q[2];
rz(-1.2008685) q[2];
sx q[2];
rz(-2.5704685) q[2];
rz(-2.5518104) q[3];
sx q[3];
rz(-0.46001205) q[3];
sx q[3];
rz(-3.0736198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.1968483) q[0];
sx q[0];
rz(-1.9510883) q[0];
sx q[0];
rz(-0.29904547) q[0];
rz(-1.8213182) q[1];
sx q[1];
rz(-2.8838005) q[1];
sx q[1];
rz(1.6437795) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70443557) q[0];
sx q[0];
rz(-2.2309982) q[0];
sx q[0];
rz(3.0755088) q[0];
rz(-0.97700714) q[2];
sx q[2];
rz(-2.0679681) q[2];
sx q[2];
rz(-1.1360816) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1491579) q[1];
sx q[1];
rz(-2.5870393) q[1];
sx q[1];
rz(-2.1402332) q[1];
x q[2];
rz(1.3214345) q[3];
sx q[3];
rz(-1.8421838) q[3];
sx q[3];
rz(0.67941487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6254639) q[2];
sx q[2];
rz(-1.4279782) q[2];
sx q[2];
rz(-3.1141172) q[2];
rz(-2.6190858) q[3];
sx q[3];
rz(-2.3404739) q[3];
sx q[3];
rz(-2.3251422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-0.41247535) q[0];
sx q[0];
rz(-1.1183879) q[0];
sx q[0];
rz(-0.11418848) q[0];
rz(0.97822899) q[1];
sx q[1];
rz(-2.5949635) q[1];
sx q[1];
rz(0.79089975) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.233477) q[0];
sx q[0];
rz(-2.18963) q[0];
sx q[0];
rz(3.0254362) q[0];
rz(-0.36206836) q[2];
sx q[2];
rz(-2.3644591) q[2];
sx q[2];
rz(-2.2287378) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.447688) q[1];
sx q[1];
rz(-1.4920007) q[1];
sx q[1];
rz(-1.9359971) q[1];
rz(-pi) q[2];
rz(-2.2579455) q[3];
sx q[3];
rz(-1.8568608) q[3];
sx q[3];
rz(-1.4253695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2699282) q[2];
sx q[2];
rz(-1.1130788) q[2];
sx q[2];
rz(0.86501914) q[2];
rz(2.4690752) q[3];
sx q[3];
rz(-2.0130242) q[3];
sx q[3];
rz(-0.095656693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6487811) q[0];
sx q[0];
rz(-0.5195986) q[0];
sx q[0];
rz(-2.6742324) q[0];
rz(-0.53721792) q[1];
sx q[1];
rz(-2.1610114) q[1];
sx q[1];
rz(0.25407243) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92914903) q[0];
sx q[0];
rz(-1.7663167) q[0];
sx q[0];
rz(1.50302) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.28283624) q[2];
sx q[2];
rz(-3.1361702) q[2];
sx q[2];
rz(-1.1495513) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.64995631) q[1];
sx q[1];
rz(-1.2600113) q[1];
sx q[1];
rz(-1.4375356) q[1];
x q[2];
rz(-0.11717511) q[3];
sx q[3];
rz(-2.5391038) q[3];
sx q[3];
rz(-2.6206827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.49999985) q[2];
sx q[2];
rz(-2.3708512) q[2];
sx q[2];
rz(-2.8059778) q[2];
rz(-2.8149758) q[3];
sx q[3];
rz(-2.2461522) q[3];
sx q[3];
rz(1.4183104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0582054) q[0];
sx q[0];
rz(-0.21033062) q[0];
sx q[0];
rz(-2.0876419) q[0];
rz(-0.15696934) q[1];
sx q[1];
rz(-1.7228246) q[1];
sx q[1];
rz(-0.98186791) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1218402) q[0];
sx q[0];
rz(-1.4593908) q[0];
sx q[0];
rz(0.41282546) q[0];
x q[1];
rz(0.1083381) q[2];
sx q[2];
rz(-1.0045369) q[2];
sx q[2];
rz(-0.94891753) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8874444) q[1];
sx q[1];
rz(-1.6283298) q[1];
sx q[1];
rz(1.7565586) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.710886) q[3];
sx q[3];
rz(-1.3804187) q[3];
sx q[3];
rz(3.0829317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2173569) q[2];
sx q[2];
rz(-2.0841667) q[2];
sx q[2];
rz(-0.99739933) q[2];
rz(2.635397) q[3];
sx q[3];
rz(-2.180407) q[3];
sx q[3];
rz(-0.034328073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2845594) q[0];
sx q[0];
rz(-2.4301346) q[0];
sx q[0];
rz(0.51666623) q[0];
rz(2.754028) q[1];
sx q[1];
rz(-1.060408) q[1];
sx q[1];
rz(-0.26836747) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3049406) q[0];
sx q[0];
rz(-0.9965082) q[0];
sx q[0];
rz(2.4313297) q[0];
rz(-2.7029413) q[2];
sx q[2];
rz(-1.029656) q[2];
sx q[2];
rz(-1.7739319) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.2942218) q[1];
sx q[1];
rz(-2.1547744) q[1];
sx q[1];
rz(-1.2482615) q[1];
rz(2.1488232) q[3];
sx q[3];
rz(-1.7939995) q[3];
sx q[3];
rz(2.2900667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2051852) q[2];
sx q[2];
rz(-2.3582017) q[2];
sx q[2];
rz(0.56383413) q[2];
rz(-1.9514203) q[3];
sx q[3];
rz(-2.1824013) q[3];
sx q[3];
rz(-0.64175516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47678369) q[0];
sx q[0];
rz(-1.2510779) q[0];
sx q[0];
rz(-1.0673987) q[0];
rz(1.8019567) q[1];
sx q[1];
rz(-1.7032774) q[1];
sx q[1];
rz(1.3443321) q[1];
rz(2.6230326) q[2];
sx q[2];
rz(-1.3052365) q[2];
sx q[2];
rz(-2.4076751) q[2];
rz(2.8292538) q[3];
sx q[3];
rz(-1.3795508) q[3];
sx q[3];
rz(3.0637904) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];