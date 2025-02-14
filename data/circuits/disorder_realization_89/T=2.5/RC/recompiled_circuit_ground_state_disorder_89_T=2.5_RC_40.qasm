OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.9417579) q[0];
sx q[0];
rz(-2.4500442) q[0];
sx q[0];
rz(-2.3681695) q[0];
rz(3.050488) q[1];
sx q[1];
rz(-0.79371047) q[1];
sx q[1];
rz(1.8082126) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1379194) q[0];
sx q[0];
rz(-0.9124476) q[0];
sx q[0];
rz(-1.5821032) q[0];
rz(-pi) q[1];
rz(2.0699507) q[2];
sx q[2];
rz(-1.3678275) q[2];
sx q[2];
rz(2.2405193) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4030093) q[1];
sx q[1];
rz(-0.75323757) q[1];
sx q[1];
rz(-2.1639216) q[1];
rz(0.01706546) q[3];
sx q[3];
rz(-1.6380042) q[3];
sx q[3];
rz(-2.1484571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.29204631) q[2];
sx q[2];
rz(-1.8201733) q[2];
sx q[2];
rz(-2.7267961) q[2];
rz(-1.241812) q[3];
sx q[3];
rz(-2.0262227) q[3];
sx q[3];
rz(-0.19150664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6280129) q[0];
sx q[0];
rz(-0.11592557) q[0];
sx q[0];
rz(-2.7017748) q[0];
rz(-0.10305931) q[1];
sx q[1];
rz(-0.28918806) q[1];
sx q[1];
rz(1.1691079) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70003521) q[0];
sx q[0];
rz(-2.1934783) q[0];
sx q[0];
rz(-0.17130987) q[0];
rz(-pi) q[1];
rz(-1.493426) q[2];
sx q[2];
rz(-2.5499509) q[2];
sx q[2];
rz(-2.735052) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7767765) q[1];
sx q[1];
rz(-2.4668982) q[1];
sx q[1];
rz(-1.8614091) q[1];
rz(1.3430994) q[3];
sx q[3];
rz(-0.72245729) q[3];
sx q[3];
rz(-1.8060773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.64572) q[2];
sx q[2];
rz(-1.7503259) q[2];
sx q[2];
rz(-0.24620852) q[2];
rz(-1.5919033) q[3];
sx q[3];
rz(-0.61384765) q[3];
sx q[3];
rz(1.7483819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.15568) q[0];
sx q[0];
rz(-1.14013) q[0];
sx q[0];
rz(0.21173665) q[0];
rz(-3.1302997) q[1];
sx q[1];
rz(-2.3432422) q[1];
sx q[1];
rz(-1.4439772) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2539719) q[0];
sx q[0];
rz(-1.6059552) q[0];
sx q[0];
rz(1.5707234) q[0];
rz(-pi) q[1];
rz(-1.3668187) q[2];
sx q[2];
rz(-2.8958671) q[2];
sx q[2];
rz(1.3947006) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8577598) q[1];
sx q[1];
rz(-0.74537828) q[1];
sx q[1];
rz(-2.5733295) q[1];
rz(-pi) q[2];
rz(0.45076799) q[3];
sx q[3];
rz(-1.4891461) q[3];
sx q[3];
rz(-1.1686366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6649365) q[2];
sx q[2];
rz(-1.0756476) q[2];
sx q[2];
rz(2.5962489) q[2];
rz(-2.2319345) q[3];
sx q[3];
rz(-2.6792512) q[3];
sx q[3];
rz(-1.2130515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.720541) q[0];
sx q[0];
rz(-0.67890972) q[0];
sx q[0];
rz(0.82756591) q[0];
rz(1.1830117) q[1];
sx q[1];
rz(-1.2748101) q[1];
sx q[1];
rz(-1.2659198) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7050305) q[0];
sx q[0];
rz(-1.4799191) q[0];
sx q[0];
rz(1.4687125) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7782486) q[2];
sx q[2];
rz(-1.5816763) q[2];
sx q[2];
rz(-1.915773) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3490679) q[1];
sx q[1];
rz(-1.5390087) q[1];
sx q[1];
rz(0.94520251) q[1];
rz(-pi) q[2];
rz(-0.87235116) q[3];
sx q[3];
rz(-2.8954801) q[3];
sx q[3];
rz(1.2239561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2563235) q[2];
sx q[2];
rz(-1.4726535) q[2];
sx q[2];
rz(-3.0598158) q[2];
rz(2.8335588) q[3];
sx q[3];
rz(-0.48397288) q[3];
sx q[3];
rz(1.5380194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1643739) q[0];
sx q[0];
rz(-2.866221) q[0];
sx q[0];
rz(-0.070505738) q[0];
rz(2.8071857) q[1];
sx q[1];
rz(-2.6102378) q[1];
sx q[1];
rz(-2.6883584) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95246661) q[0];
sx q[0];
rz(-2.1299908) q[0];
sx q[0];
rz(-2.1711012) q[0];
x q[1];
rz(0.22780995) q[2];
sx q[2];
rz(-1.1880298) q[2];
sx q[2];
rz(1.5973951) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6649225) q[1];
sx q[1];
rz(-1.3241265) q[1];
sx q[1];
rz(1.6849405) q[1];
rz(-2.2599561) q[3];
sx q[3];
rz(-0.83587468) q[3];
sx q[3];
rz(2.3595032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.8411023) q[2];
sx q[2];
rz(-1.6359685) q[2];
sx q[2];
rz(2.9312768) q[2];
rz(0.39120832) q[3];
sx q[3];
rz(-2.936383) q[3];
sx q[3];
rz(0.44262639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6415569) q[0];
sx q[0];
rz(-2.0099202) q[0];
sx q[0];
rz(0.17182194) q[0];
rz(1.7956644) q[1];
sx q[1];
rz(-0.95564061) q[1];
sx q[1];
rz(1.1489493) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1412107) q[0];
sx q[0];
rz(-0.11185574) q[0];
sx q[0];
rz(-0.07312183) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9093018) q[2];
sx q[2];
rz(-1.293876) q[2];
sx q[2];
rz(1.0804707) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2021804) q[1];
sx q[1];
rz(-1.3883852) q[1];
sx q[1];
rz(0.98947452) q[1];
rz(-pi) q[2];
x q[2];
rz(0.20307417) q[3];
sx q[3];
rz(-1.0725029) q[3];
sx q[3];
rz(-2.4808675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.020784) q[2];
sx q[2];
rz(-2.1773931) q[2];
sx q[2];
rz(-0.42050335) q[2];
rz(-1.6546107) q[3];
sx q[3];
rz(-1.9872811) q[3];
sx q[3];
rz(-1.9871064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3625951) q[0];
sx q[0];
rz(-1.915246) q[0];
sx q[0];
rz(-2.763789) q[0];
rz(1.3899577) q[1];
sx q[1];
rz(-1.708834) q[1];
sx q[1];
rz(2.7094944) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8209131) q[0];
sx q[0];
rz(-2.8519676) q[0];
sx q[0];
rz(0.65307133) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.34382208) q[2];
sx q[2];
rz(-0.09843407) q[2];
sx q[2];
rz(-1.7119031) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3687146) q[1];
sx q[1];
rz(-1.1091103) q[1];
sx q[1];
rz(1.6108455) q[1];
x q[2];
rz(2.0666408) q[3];
sx q[3];
rz(-0.73710873) q[3];
sx q[3];
rz(0.51966156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3720588) q[2];
sx q[2];
rz(-0.62166119) q[2];
sx q[2];
rz(0.2335693) q[2];
rz(-1.4522878) q[3];
sx q[3];
rz(-1.4315616) q[3];
sx q[3];
rz(1.5805894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5643144) q[0];
sx q[0];
rz(-2.1121341) q[0];
sx q[0];
rz(-2.8218063) q[0];
rz(2.2794967) q[1];
sx q[1];
rz(-1.2513221) q[1];
sx q[1];
rz(-1.3592985) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0657227) q[0];
sx q[0];
rz(-1.0130873) q[0];
sx q[0];
rz(0.88904492) q[0];
x q[1];
rz(2.1807387) q[2];
sx q[2];
rz(-0.24343389) q[2];
sx q[2];
rz(3.0185901) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.086445623) q[1];
sx q[1];
rz(-2.3999955) q[1];
sx q[1];
rz(-2.8838709) q[1];
rz(-pi) q[2];
x q[2];
rz(0.27761658) q[3];
sx q[3];
rz(-0.95562387) q[3];
sx q[3];
rz(2.6108116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.84997815) q[2];
sx q[2];
rz(-0.41899592) q[2];
sx q[2];
rz(2.6739547) q[2];
rz(2.6575798) q[3];
sx q[3];
rz(-1.368112) q[3];
sx q[3];
rz(1.8638994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9676232) q[0];
sx q[0];
rz(-1.4894217) q[0];
sx q[0];
rz(2.0066579) q[0];
rz(3.0813772) q[1];
sx q[1];
rz(-2.4048012) q[1];
sx q[1];
rz(-1.6103475) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13847199) q[0];
sx q[0];
rz(-1.1116189) q[0];
sx q[0];
rz(-0.91158406) q[0];
rz(-0.61988735) q[2];
sx q[2];
rz(-0.7893749) q[2];
sx q[2];
rz(2.3933586) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.169226) q[1];
sx q[1];
rz(-1.3080096) q[1];
sx q[1];
rz(1.6714833) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6100735) q[3];
sx q[3];
rz(-0.78900064) q[3];
sx q[3];
rz(1.6226744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.81426364) q[2];
sx q[2];
rz(-1.9244104) q[2];
sx q[2];
rz(-2.4129756) q[2];
rz(-0.21555756) q[3];
sx q[3];
rz(-1.7451127) q[3];
sx q[3];
rz(2.047915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88290596) q[0];
sx q[0];
rz(-0.031143324) q[0];
sx q[0];
rz(-2.8236142) q[0];
rz(-1.7440354) q[1];
sx q[1];
rz(-1.8122858) q[1];
sx q[1];
rz(0.2831645) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8737333) q[0];
sx q[0];
rz(-1.443131) q[0];
sx q[0];
rz(-3.0588849) q[0];
rz(0.7541025) q[2];
sx q[2];
rz(-1.6267383) q[2];
sx q[2];
rz(-2.2365582) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6826815) q[1];
sx q[1];
rz(-1.1106756) q[1];
sx q[1];
rz(-2.3691872) q[1];
x q[2];
rz(3.0491676) q[3];
sx q[3];
rz(-1.3636949) q[3];
sx q[3];
rz(2.7566119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2240923) q[2];
sx q[2];
rz(-1.6250236) q[2];
sx q[2];
rz(-3.0523114) q[2];
rz(-1.3034405) q[3];
sx q[3];
rz(-2.3035514) q[3];
sx q[3];
rz(-2.1681521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3598809) q[0];
sx q[0];
rz(-0.64890535) q[0];
sx q[0];
rz(-2.1606408) q[0];
rz(2.5557062) q[1];
sx q[1];
rz(-1.8460907) q[1];
sx q[1];
rz(1.5150217) q[1];
rz(2.9286251) q[2];
sx q[2];
rz(-1.3942547) q[2];
sx q[2];
rz(-0.19993776) q[2];
rz(0.99021749) q[3];
sx q[3];
rz(-0.44689842) q[3];
sx q[3];
rz(-2.3200271) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
