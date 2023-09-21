OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.8310175) q[0];
sx q[0];
rz(-2.0895045) q[0];
sx q[0];
rz(-1.6488099) q[0];
rz(-2.2812023) q[1];
sx q[1];
rz(-0.69505039) q[1];
sx q[1];
rz(1.0746497) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2225392) q[0];
sx q[0];
rz(-3.1226282) q[0];
sx q[0];
rz(2.8595964) q[0];
rz(-0.78975241) q[2];
sx q[2];
rz(-1.6842004) q[2];
sx q[2];
rz(-1.5577158) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5584471) q[1];
sx q[1];
rz(-1.6513283) q[1];
sx q[1];
rz(-0.12195485) q[1];
rz(2.6082637) q[3];
sx q[3];
rz(-1.1733857) q[3];
sx q[3];
rz(1.0226137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.82912123) q[2];
sx q[2];
rz(-2.1437936) q[2];
sx q[2];
rz(-1.5134229) q[2];
rz(1.1734022) q[3];
sx q[3];
rz(-1.6956001) q[3];
sx q[3];
rz(-2.9287958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56869498) q[0];
sx q[0];
rz(-0.64210367) q[0];
sx q[0];
rz(-1.2233541) q[0];
rz(-0.16076316) q[1];
sx q[1];
rz(-1.5630961) q[1];
sx q[1];
rz(1.6404023) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26186865) q[0];
sx q[0];
rz(-0.12517087) q[0];
sx q[0];
rz(2.9414888) q[0];
rz(0.31200774) q[2];
sx q[2];
rz(-1.1216251) q[2];
sx q[2];
rz(0.95865196) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7565631) q[1];
sx q[1];
rz(-0.90365138) q[1];
sx q[1];
rz(-1.3515616) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1797446) q[3];
sx q[3];
rz(-1.565233) q[3];
sx q[3];
rz(0.55913505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.286065) q[2];
sx q[2];
rz(-0.76670402) q[2];
sx q[2];
rz(1.5141053) q[2];
rz(-1.1668011) q[3];
sx q[3];
rz(-1.5224001) q[3];
sx q[3];
rz(-2.2272026) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0062155) q[0];
sx q[0];
rz(-2.0897946) q[0];
sx q[0];
rz(-2.0626383) q[0];
rz(-1.9619933) q[1];
sx q[1];
rz(-1.0410573) q[1];
sx q[1];
rz(1.6378145) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2730946) q[0];
sx q[0];
rz(-1.4833741) q[0];
sx q[0];
rz(-0.065386535) q[0];
rz(-pi) q[1];
rz(0.9861228) q[2];
sx q[2];
rz(-1.1650656) q[2];
sx q[2];
rz(3.0509146) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1551731) q[1];
sx q[1];
rz(-1.5281614) q[1];
sx q[1];
rz(-0.30100545) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.91434874) q[3];
sx q[3];
rz(-1.4667061) q[3];
sx q[3];
rz(2.5641233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7109795) q[2];
sx q[2];
rz(-0.47982275) q[2];
sx q[2];
rz(-0.7437931) q[2];
rz(-2.8841694) q[3];
sx q[3];
rz(-2.0580723) q[3];
sx q[3];
rz(-2.553489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1733615) q[0];
sx q[0];
rz(-0.1156725) q[0];
sx q[0];
rz(-1.051735) q[0];
rz(0.60032088) q[1];
sx q[1];
rz(-0.61906639) q[1];
sx q[1];
rz(1.1490885) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34654348) q[0];
sx q[0];
rz(-1.330089) q[0];
sx q[0];
rz(-1.4562777) q[0];
rz(-pi) q[1];
rz(-2.7096268) q[2];
sx q[2];
rz(-1.0216121) q[2];
sx q[2];
rz(0.84749046) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2266352) q[1];
sx q[1];
rz(-1.2775704) q[1];
sx q[1];
rz(-0.83421591) q[1];
rz(-pi) q[2];
rz(1.9150919) q[3];
sx q[3];
rz(-0.85840423) q[3];
sx q[3];
rz(2.815849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8950618) q[2];
sx q[2];
rz(-1.4738513) q[2];
sx q[2];
rz(-0.66876283) q[2];
rz(-1.8956005) q[3];
sx q[3];
rz(-0.99530882) q[3];
sx q[3];
rz(0.016383735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.4887061) q[0];
sx q[0];
rz(-1.9911433) q[0];
sx q[0];
rz(-2.0892129) q[0];
rz(1.2106238) q[1];
sx q[1];
rz(-1.4167507) q[1];
sx q[1];
rz(-2.2573684) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.592504) q[0];
sx q[0];
rz(-1.7486524) q[0];
sx q[0];
rz(1.7643719) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0074559) q[2];
sx q[2];
rz(-1.9433125) q[2];
sx q[2];
rz(2.8357752) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.21806949) q[1];
sx q[1];
rz(-1.7150208) q[1];
sx q[1];
rz(-1.5037392) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7719284) q[3];
sx q[3];
rz(-2.0737994) q[3];
sx q[3];
rz(0.30077416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9050682) q[2];
sx q[2];
rz(-0.69513598) q[2];
sx q[2];
rz(2.0489676) q[2];
rz(2.8975899) q[3];
sx q[3];
rz(-2.2677393) q[3];
sx q[3];
rz(-2.3905335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7364175) q[0];
sx q[0];
rz(-0.002679499) q[0];
sx q[0];
rz(-1.8238235) q[0];
rz(-2.4353943) q[1];
sx q[1];
rz(-1.6786007) q[1];
sx q[1];
rz(-0.96907369) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5613865) q[0];
sx q[0];
rz(-1.5302587) q[0];
sx q[0];
rz(1.4617306) q[0];
x q[1];
rz(-1.5029807) q[2];
sx q[2];
rz(-2.2837451) q[2];
sx q[2];
rz(-2.3358047) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.25989306) q[1];
sx q[1];
rz(-2.0310146) q[1];
sx q[1];
rz(-2.0395181) q[1];
rz(-pi) q[2];
x q[2];
rz(8*pi/11) q[3];
sx q[3];
rz(-2.487605) q[3];
sx q[3];
rz(1.6884782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6270854) q[2];
sx q[2];
rz(-0.58834499) q[2];
sx q[2];
rz(0.43760854) q[2];
rz(3.0833516) q[3];
sx q[3];
rz(-0.85844675) q[3];
sx q[3];
rz(-1.5313914) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3095734) q[0];
sx q[0];
rz(-2.11634) q[0];
sx q[0];
rz(1.3597885) q[0];
rz(1.4490022) q[1];
sx q[1];
rz(-2.2429008) q[1];
sx q[1];
rz(-1.4356027) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3754417) q[0];
sx q[0];
rz(-1.0784082) q[0];
sx q[0];
rz(-2.173461) q[0];
rz(-pi) q[1];
x q[1];
rz(0.018354015) q[2];
sx q[2];
rz(-1.5510501) q[2];
sx q[2];
rz(-2.8543618) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2885399) q[1];
sx q[1];
rz(-1.9318046) q[1];
sx q[1];
rz(0.32220528) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7123659) q[3];
sx q[3];
rz(-1.1042522) q[3];
sx q[3];
rz(2.6847117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.74063611) q[2];
sx q[2];
rz(-1.2601968) q[2];
sx q[2];
rz(2.9675972) q[2];
rz(2.1334355) q[3];
sx q[3];
rz(-2.5139136) q[3];
sx q[3];
rz(-0.15667668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.290264) q[0];
sx q[0];
rz(-2.2785318) q[0];
sx q[0];
rz(0.1329578) q[0];
rz(2.6538387) q[1];
sx q[1];
rz(-2.1552591) q[1];
sx q[1];
rz(1.7763604) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93850219) q[0];
sx q[0];
rz(-1.7068958) q[0];
sx q[0];
rz(0.076119856) q[0];
rz(-pi) q[1];
rz(-1.4012785) q[2];
sx q[2];
rz(-3.0119544) q[2];
sx q[2];
rz(-2.4469751) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3268765) q[1];
sx q[1];
rz(-1.6232345) q[1];
sx q[1];
rz(-0.53175064) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.385231) q[3];
sx q[3];
rz(-1.984333) q[3];
sx q[3];
rz(2.6698649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.075295538) q[2];
sx q[2];
rz(-1.8981045) q[2];
sx q[2];
rz(2.7189642) q[2];
rz(-2.1832809) q[3];
sx q[3];
rz(-0.13923968) q[3];
sx q[3];
rz(2.9039834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1388824) q[0];
sx q[0];
rz(-0.30160987) q[0];
sx q[0];
rz(0.31059206) q[0];
rz(-0.25767913) q[1];
sx q[1];
rz(-1.5035166) q[1];
sx q[1];
rz(2.2176567) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3850708) q[0];
sx q[0];
rz(-1.1413045) q[0];
sx q[0];
rz(0.71682741) q[0];
rz(-pi) q[1];
x q[1];
rz(0.02248259) q[2];
sx q[2];
rz(-0.48869952) q[2];
sx q[2];
rz(1.2475916) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1549346) q[1];
sx q[1];
rz(-1.0549874) q[1];
sx q[1];
rz(2.6845279) q[1];
rz(0.71420788) q[3];
sx q[3];
rz(-2.0651157) q[3];
sx q[3];
rz(0.89356542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.49736398) q[2];
sx q[2];
rz(-2.6022537) q[2];
sx q[2];
rz(1.7473934) q[2];
rz(1.1692858) q[3];
sx q[3];
rz(-2.101427) q[3];
sx q[3];
rz(-2.6031301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.37344638) q[0];
sx q[0];
rz(-3.0420619) q[0];
sx q[0];
rz(1.753153) q[0];
rz(-3.1346079) q[1];
sx q[1];
rz(-2.3313637) q[1];
sx q[1];
rz(-3.1034234) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.105135) q[0];
sx q[0];
rz(-1.5441226) q[0];
sx q[0];
rz(2.1795576) q[0];
x q[1];
rz(1.795904) q[2];
sx q[2];
rz(-2.7499866) q[2];
sx q[2];
rz(3.0181146) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0109509) q[1];
sx q[1];
rz(-2.7656733) q[1];
sx q[1];
rz(1.6846659) q[1];
rz(-0.90922728) q[3];
sx q[3];
rz(-1.070676) q[3];
sx q[3];
rz(1.7928894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1756211) q[2];
sx q[2];
rz(-1.9922545) q[2];
sx q[2];
rz(-1.1981298) q[2];
rz(-1.0839869) q[3];
sx q[3];
rz(-2.4626648) q[3];
sx q[3];
rz(-2.9057124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9772298) q[0];
sx q[0];
rz(-1.2156163) q[0];
sx q[0];
rz(-1.6396133) q[0];
rz(1.6434796) q[1];
sx q[1];
rz(-2.5479981) q[1];
sx q[1];
rz(2.610276) q[1];
rz(0.48508587) q[2];
sx q[2];
rz(-3.017364) q[2];
sx q[2];
rz(-0.45807522) q[2];
rz(-2.2723972) q[3];
sx q[3];
rz(-1.8164608) q[3];
sx q[3];
rz(1.7432004) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
