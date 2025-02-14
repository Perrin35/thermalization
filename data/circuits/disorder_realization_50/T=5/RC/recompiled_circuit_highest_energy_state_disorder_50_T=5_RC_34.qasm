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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8136381) q[0];
sx q[0];
rz(-1.0109954) q[0];
sx q[0];
rz(-0.35388057) q[0];
rz(-pi) q[1];
x q[1];
rz(0.73424642) q[2];
sx q[2];
rz(-0.99319211) q[2];
sx q[2];
rz(-2.2478888) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.96054582) q[1];
sx q[1];
rz(-1.6072453) q[1];
sx q[1];
rz(-0.44238449) q[1];
x q[2];
rz(1.5230701) q[3];
sx q[3];
rz(-0.56079799) q[3];
sx q[3];
rz(-2.5530397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.36246768) q[2];
sx q[2];
rz(-0.96131009) q[2];
sx q[2];
rz(2.989952) q[2];
rz(-0.14196299) q[3];
sx q[3];
rz(-1.4700593) q[3];
sx q[3];
rz(2.0040472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4769984) q[0];
sx q[0];
rz(-2.4095896) q[0];
sx q[0];
rz(-0.81749302) q[0];
rz(-0.11257653) q[1];
sx q[1];
rz(-1.45603) q[1];
sx q[1];
rz(-2.2419825) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0257638) q[0];
sx q[0];
rz(-1.2983822) q[0];
sx q[0];
rz(-1.8194356) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4075872) q[2];
sx q[2];
rz(-0.81799928) q[2];
sx q[2];
rz(-2.1591275) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6625329) q[1];
sx q[1];
rz(-2.3182437) q[1];
sx q[1];
rz(1.4940628) q[1];
rz(-pi) q[2];
rz(0.894015) q[3];
sx q[3];
rz(-1.9872267) q[3];
sx q[3];
rz(-0.82899603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.67109913) q[2];
sx q[2];
rz(-1.6464536) q[2];
sx q[2];
rz(1.0279083) q[2];
rz(2.4043064) q[3];
sx q[3];
rz(-1.6643915) q[3];
sx q[3];
rz(0.024287311) q[3];
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
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1360433) q[0];
sx q[0];
rz(-0.5558973) q[0];
sx q[0];
rz(-1.9480202) q[0];
rz(1.8434803) q[1];
sx q[1];
rz(-1.6622512) q[1];
sx q[1];
rz(1.6568291) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95052108) q[0];
sx q[0];
rz(-1.5245887) q[0];
sx q[0];
rz(0.025006983) q[0];
rz(-pi) q[1];
rz(-0.12642352) q[2];
sx q[2];
rz(-1.3609386) q[2];
sx q[2];
rz(-1.017638) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.79968184) q[1];
sx q[1];
rz(-1.7060301) q[1];
sx q[1];
rz(-0.34457259) q[1];
rz(-pi) q[2];
rz(3.0590634) q[3];
sx q[3];
rz(-1.9325053) q[3];
sx q[3];
rz(-1.495273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.16126157) q[2];
sx q[2];
rz(-1.9183466) q[2];
sx q[2];
rz(-0.70183357) q[2];
rz(3.0724604) q[3];
sx q[3];
rz(-2.2528503) q[3];
sx q[3];
rz(0.70772901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.055701) q[0];
sx q[0];
rz(-2.0022855) q[0];
sx q[0];
rz(0.62537801) q[0];
rz(1.0944132) q[1];
sx q[1];
rz(-1.6182599) q[1];
sx q[1];
rz(1.8249003) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65526456) q[0];
sx q[0];
rz(-2.1127955) q[0];
sx q[0];
rz(0.16714759) q[0];
x q[1];
rz(0.32347347) q[2];
sx q[2];
rz(-1.4789561) q[2];
sx q[2];
rz(2.8367281) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.55344501) q[1];
sx q[1];
rz(-2.4134753) q[1];
sx q[1];
rz(-0.27256984) q[1];
x q[2];
rz(0.80009256) q[3];
sx q[3];
rz(-1.1372107) q[3];
sx q[3];
rz(-2.3326002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.40654287) q[2];
sx q[2];
rz(-2.4989765) q[2];
sx q[2];
rz(-3.1114846) q[2];
rz(2.7249469) q[3];
sx q[3];
rz(-1.4285587) q[3];
sx q[3];
rz(0.055195181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(1.9238506) q[0];
rz(-2.9171004) q[1];
sx q[1];
rz(-1.6480564) q[1];
sx q[1];
rz(0.40103689) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.359085) q[0];
sx q[0];
rz(-1.9422724) q[0];
sx q[0];
rz(-1.1228193) q[0];
rz(-pi) q[1];
rz(-0.46868344) q[2];
sx q[2];
rz(-2.4309553) q[2];
sx q[2];
rz(0.24220322) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.89304533) q[1];
sx q[1];
rz(-1.6181989) q[1];
sx q[1];
rz(0.23706146) q[1];
rz(-pi) q[2];
rz(-2.5878536) q[3];
sx q[3];
rz(-0.19036346) q[3];
sx q[3];
rz(0.37356627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2991221) q[2];
sx q[2];
rz(-1.9186019) q[2];
sx q[2];
rz(-2.6341338) q[2];
rz(2.2211645) q[3];
sx q[3];
rz(-2.7046552) q[3];
sx q[3];
rz(-2.9612655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64666635) q[0];
sx q[0];
rz(-0.05519069) q[0];
sx q[0];
rz(-1.0194417) q[0];
rz(-1.1722209) q[1];
sx q[1];
rz(-1.1727138) q[1];
sx q[1];
rz(-1.9196462) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9393765) q[0];
sx q[0];
rz(-0.76986137) q[0];
sx q[0];
rz(-1.059993) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1062791) q[2];
sx q[2];
rz(-2.0853634) q[2];
sx q[2];
rz(0.071431486) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3675569) q[1];
sx q[1];
rz(-1.2962711) q[1];
sx q[1];
rz(-0.56660272) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1453247) q[3];
sx q[3];
rz(-1.3187416) q[3];
sx q[3];
rz(1.9854522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.4814066) q[2];
sx q[2];
rz(-0.6654827) q[2];
sx q[2];
rz(-2.0484203) q[2];
rz(2.6045351) q[3];
sx q[3];
rz(-1.6780746) q[3];
sx q[3];
rz(0.66948906) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90401232) q[0];
sx q[0];
rz(-0.31929382) q[0];
sx q[0];
rz(-1.681666) q[0];
rz(1.5096674) q[1];
sx q[1];
rz(-2.5131112) q[1];
sx q[1];
rz(-0.07930886) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3103762) q[0];
sx q[0];
rz(-0.25070158) q[0];
sx q[0];
rz(1.2396856) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5517188) q[2];
sx q[2];
rz(-0.24352077) q[2];
sx q[2];
rz(0.68984725) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.435959) q[1];
sx q[1];
rz(-0.50928947) q[1];
sx q[1];
rz(-2.5540256) q[1];
rz(-pi) q[2];
rz(-2.4238911) q[3];
sx q[3];
rz(-1.4677748) q[3];
sx q[3];
rz(-0.75005248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1377533) q[2];
sx q[2];
rz(-1.1649818) q[2];
sx q[2];
rz(-0.4129146) q[2];
rz(1.2952992) q[3];
sx q[3];
rz(-2.2173939) q[3];
sx q[3];
rz(0.41954654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9846648) q[0];
sx q[0];
rz(-2.4913737) q[0];
sx q[0];
rz(-0.28537634) q[0];
rz(-3.0889619) q[1];
sx q[1];
rz(-1.7335408) q[1];
sx q[1];
rz(2.9579128) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7010371) q[0];
sx q[0];
rz(-1.6822878) q[0];
sx q[0];
rz(-0.035758408) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6551774) q[2];
sx q[2];
rz(-0.32054311) q[2];
sx q[2];
rz(-2.965791) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7990098) q[1];
sx q[1];
rz(-1.4945806) q[1];
sx q[1];
rz(-0.59554312) q[1];
x q[2];
rz(-2.4934108) q[3];
sx q[3];
rz(-1.8055918) q[3];
sx q[3];
rz(-0.33644331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.24667428) q[2];
sx q[2];
rz(-1.4072714) q[2];
sx q[2];
rz(-1.0651917) q[2];
rz(-1.4554321) q[3];
sx q[3];
rz(-1.287241) q[3];
sx q[3];
rz(-2.5559032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10373779) q[0];
sx q[0];
rz(-2.3268564) q[0];
sx q[0];
rz(-1.6023585) q[0];
rz(2.1922951) q[1];
sx q[1];
rz(-1.0485336) q[1];
sx q[1];
rz(-1.4303713) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8447782) q[0];
sx q[0];
rz(-1.0460633) q[0];
sx q[0];
rz(1.0206166) q[0];
rz(-0.075242356) q[2];
sx q[2];
rz(-0.84343592) q[2];
sx q[2];
rz(0.6211578) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.51586823) q[1];
sx q[1];
rz(-1.5411545) q[1];
sx q[1];
rz(-0.29582204) q[1];
rz(-pi) q[2];
rz(-1.0046047) q[3];
sx q[3];
rz(-2.6454685) q[3];
sx q[3];
rz(-1.8593018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0091693) q[2];
sx q[2];
rz(-1.3045661) q[2];
sx q[2];
rz(-3.0150748) q[2];
rz(-2.8070519) q[3];
sx q[3];
rz(-2.8821475) q[3];
sx q[3];
rz(0.87219605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8406521) q[0];
sx q[0];
rz(-2.2569188) q[0];
sx q[0];
rz(2.6244923) q[0];
rz(0.05489796) q[1];
sx q[1];
rz(-1.6322735) q[1];
sx q[1];
rz(-3.0518234) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.207491) q[0];
sx q[0];
rz(-0.708552) q[0];
sx q[0];
rz(-0.92325489) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1316809) q[2];
sx q[2];
rz(-1.5176306) q[2];
sx q[2];
rz(1.9694984) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.43347142) q[1];
sx q[1];
rz(-2.8420487) q[1];
sx q[1];
rz(0.56510651) q[1];
x q[2];
rz(-1.0221972) q[3];
sx q[3];
rz(-2.7226603) q[3];
sx q[3];
rz(1.054712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.53139293) q[2];
sx q[2];
rz(-1.1343845) q[2];
sx q[2];
rz(-2.4143977) q[2];
rz(-0.7555035) q[3];
sx q[3];
rz(-2.754039) q[3];
sx q[3];
rz(-0.21272794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1995734) q[0];
sx q[0];
rz(-1.5409536) q[0];
sx q[0];
rz(1.090747) q[0];
rz(-1.3923116) q[1];
sx q[1];
rz(-1.5131469) q[1];
sx q[1];
rz(-1.6319235) q[1];
rz(0.86633273) q[2];
sx q[2];
rz(-1.6523747) q[2];
sx q[2];
rz(-2.0518641) q[2];
rz(-2.5531664) q[3];
sx q[3];
rz(-1.0259494) q[3];
sx q[3];
rz(-1.1001724) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
