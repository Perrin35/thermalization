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
rz(3.0124445) q[0];
sx q[0];
rz(-1.4718066) q[0];
sx q[0];
rz(2.1762525) q[0];
rz(-3.1211634) q[1];
sx q[1];
rz(-1.6579707) q[1];
sx q[1];
rz(1.3774011) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5063874) q[0];
sx q[0];
rz(-2.0034408) q[0];
sx q[0];
rz(-2.3385027) q[0];
rz(1.6135869) q[2];
sx q[2];
rz(-0.64633179) q[2];
sx q[2];
rz(-0.94852322) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.5038576) q[1];
sx q[1];
rz(-1.4252932) q[1];
sx q[1];
rz(-2.9521431) q[1];
rz(-0.018947424) q[3];
sx q[3];
rz(-2.3403882) q[3];
sx q[3];
rz(-0.5294746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.91780245) q[2];
sx q[2];
rz(-1.4897572) q[2];
sx q[2];
rz(-1.6415143) q[2];
rz(0.7817868) q[3];
sx q[3];
rz(-0.89038554) q[3];
sx q[3];
rz(-2.3894943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.016945275) q[0];
sx q[0];
rz(-2.3638159) q[0];
sx q[0];
rz(0.97050226) q[0];
rz(-1.8151201) q[1];
sx q[1];
rz(-1.7992203) q[1];
sx q[1];
rz(0.91189799) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8489979) q[0];
sx q[0];
rz(-2.1540897) q[0];
sx q[0];
rz(1.769577) q[0];
x q[1];
rz(-1.3746337) q[2];
sx q[2];
rz(-1.4028143) q[2];
sx q[2];
rz(-0.84174918) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.40373793) q[1];
sx q[1];
rz(-2.4399839) q[1];
sx q[1];
rz(1.5075633) q[1];
x q[2];
rz(1.1158783) q[3];
sx q[3];
rz(-1.5568241) q[3];
sx q[3];
rz(-0.25532237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9087002) q[2];
sx q[2];
rz(-1.1129271) q[2];
sx q[2];
rz(-2.1523037) q[2];
rz(-0.42012897) q[3];
sx q[3];
rz(-2.1412854) q[3];
sx q[3];
rz(-0.72107983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.99783889) q[0];
sx q[0];
rz(-1.9793352) q[0];
sx q[0];
rz(1.0379399) q[0];
rz(2.2833917) q[1];
sx q[1];
rz(-2.61519) q[1];
sx q[1];
rz(0.16955489) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35992453) q[0];
sx q[0];
rz(-1.5667652) q[0];
sx q[0];
rz(1.7143334) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.900678) q[2];
sx q[2];
rz(-0.49374106) q[2];
sx q[2];
rz(2.7639219) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.922561) q[1];
sx q[1];
rz(-1.2777849) q[1];
sx q[1];
rz(1.0508363) q[1];
rz(-pi) q[2];
rz(-2.9870728) q[3];
sx q[3];
rz(-1.2121046) q[3];
sx q[3];
rz(1.7254988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1797336) q[2];
sx q[2];
rz(-2.6229975) q[2];
sx q[2];
rz(-0.40599424) q[2];
rz(-2.3640442) q[3];
sx q[3];
rz(-0.33027875) q[3];
sx q[3];
rz(2.3269261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4985713) q[0];
sx q[0];
rz(-1.850147) q[0];
sx q[0];
rz(0.92798573) q[0];
rz(-3.0827177) q[1];
sx q[1];
rz(-1.6935655) q[1];
sx q[1];
rz(1.7108797) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1851981) q[0];
sx q[0];
rz(-1.1692337) q[0];
sx q[0];
rz(-0.08423452) q[0];
rz(-pi) q[1];
rz(-2.7785382) q[2];
sx q[2];
rz(-1.7573234) q[2];
sx q[2];
rz(2.8721953) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2637564) q[1];
sx q[1];
rz(-2.0621513) q[1];
sx q[1];
rz(0.25313913) q[1];
rz(-pi) q[2];
rz(-1.1673314) q[3];
sx q[3];
rz(-0.66934062) q[3];
sx q[3];
rz(0.40237967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6835988) q[2];
sx q[2];
rz(-1.6432089) q[2];
sx q[2];
rz(1.9963473) q[2];
rz(-2.2942719) q[3];
sx q[3];
rz(-1.8073795) q[3];
sx q[3];
rz(1.5020717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1640846) q[0];
sx q[0];
rz(-1.1825528) q[0];
sx q[0];
rz(0.384828) q[0];
rz(-0.79967868) q[1];
sx q[1];
rz(-2.622602) q[1];
sx q[1];
rz(-1.5934561) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1037796) q[0];
sx q[0];
rz(-1.3290958) q[0];
sx q[0];
rz(1.3194618) q[0];
x q[1];
rz(2.9351531) q[2];
sx q[2];
rz(-1.9597133) q[2];
sx q[2];
rz(-2.2909475) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0210798) q[1];
sx q[1];
rz(-1.2215633) q[1];
sx q[1];
rz(1.2362739) q[1];
rz(-pi) q[2];
rz(-1.0043275) q[3];
sx q[3];
rz(-2.5133555) q[3];
sx q[3];
rz(0.35997501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3758214) q[2];
sx q[2];
rz(-2.4479726) q[2];
sx q[2];
rz(0.080246298) q[2];
rz(-2.5806228) q[3];
sx q[3];
rz(-1.4813981) q[3];
sx q[3];
rz(2.5188353) q[3];
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
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65492594) q[0];
sx q[0];
rz(-0.66513649) q[0];
sx q[0];
rz(-1.8810062) q[0];
rz(-0.27443019) q[1];
sx q[1];
rz(-1.471289) q[1];
sx q[1];
rz(2.4226277) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.052375) q[0];
sx q[0];
rz(-1.9196654) q[0];
sx q[0];
rz(-0.87651663) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0098835398) q[2];
sx q[2];
rz(-2.3172969) q[2];
sx q[2];
rz(0.020992779) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5605041) q[1];
sx q[1];
rz(-1.6916654) q[1];
sx q[1];
rz(0.64847364) q[1];
rz(-pi) q[2];
rz(2.9230887) q[3];
sx q[3];
rz(-1.2077304) q[3];
sx q[3];
rz(-0.6930815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6545973) q[2];
sx q[2];
rz(-1.2664814) q[2];
sx q[2];
rz(0.60301644) q[2];
rz(1.6675789) q[3];
sx q[3];
rz(-1.0242198) q[3];
sx q[3];
rz(2.730864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6239768) q[0];
sx q[0];
rz(-1.5331974) q[0];
sx q[0];
rz(0.18727592) q[0];
rz(-2.1693443) q[1];
sx q[1];
rz(-0.15521237) q[1];
sx q[1];
rz(-0.42047277) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12026006) q[0];
sx q[0];
rz(-1.9870523) q[0];
sx q[0];
rz(0.72159213) q[0];
rz(-pi) q[1];
rz(0.37897972) q[2];
sx q[2];
rz(-1.0595269) q[2];
sx q[2];
rz(-0.015470964) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.5198316) q[1];
sx q[1];
rz(-2.6606307) q[1];
sx q[1];
rz(0.93641187) q[1];
x q[2];
rz(1.917508) q[3];
sx q[3];
rz(-0.60255749) q[3];
sx q[3];
rz(1.4977221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.04574) q[2];
sx q[2];
rz(-1.1008215) q[2];
sx q[2];
rz(-2.8908758) q[2];
rz(-1.8935253) q[3];
sx q[3];
rz(-1.3919132) q[3];
sx q[3];
rz(-0.59984508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2774169) q[0];
sx q[0];
rz(-0.59159788) q[0];
sx q[0];
rz(2.199882) q[0];
rz(-0.038453728) q[1];
sx q[1];
rz(-1.7105303) q[1];
sx q[1];
rz(-3.0252735) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7021892) q[0];
sx q[0];
rz(-2.2941414) q[0];
sx q[0];
rz(2.1795616) q[0];
rz(-0.88960008) q[2];
sx q[2];
rz(-1.7187013) q[2];
sx q[2];
rz(-3.0513632) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6623936) q[1];
sx q[1];
rz(-0.31649193) q[1];
sx q[1];
rz(-2.5280158) q[1];
rz(-pi) q[2];
rz(1.6809471) q[3];
sx q[3];
rz(-1.9578736) q[3];
sx q[3];
rz(0.36110669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.8470799) q[2];
sx q[2];
rz(-2.3307762) q[2];
sx q[2];
rz(-1.026356) q[2];
rz(1.2096842) q[3];
sx q[3];
rz(-2.1116202) q[3];
sx q[3];
rz(2.1799555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.472979) q[0];
sx q[0];
rz(-2.0675779) q[0];
sx q[0];
rz(0.37856722) q[0];
rz(2.0319132) q[1];
sx q[1];
rz(-2.8329284) q[1];
sx q[1];
rz(-3.1139156) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33998128) q[0];
sx q[0];
rz(-1.2072354) q[0];
sx q[0];
rz(2.5021863) q[0];
rz(-pi) q[1];
rz(-2.2760324) q[2];
sx q[2];
rz(-2.8606374) q[2];
sx q[2];
rz(-2.11175) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.81144801) q[1];
sx q[1];
rz(-1.7551577) q[1];
sx q[1];
rz(-1.2069279) q[1];
x q[2];
rz(1.1430986) q[3];
sx q[3];
rz(-1.2790247) q[3];
sx q[3];
rz(-1.9960718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3045584) q[2];
sx q[2];
rz(-2.1874032) q[2];
sx q[2];
rz(0.41998106) q[2];
rz(0.84152451) q[3];
sx q[3];
rz(-2.1244815) q[3];
sx q[3];
rz(-0.37874547) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76641744) q[0];
sx q[0];
rz(-0.11688047) q[0];
sx q[0];
rz(0.28512678) q[0];
rz(0.20052234) q[1];
sx q[1];
rz(-1.2031809) q[1];
sx q[1];
rz(0.83555046) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.052976089) q[0];
sx q[0];
rz(-1.1167913) q[0];
sx q[0];
rz(2.3331621) q[0];
rz(-1.5642688) q[2];
sx q[2];
rz(-0.92374974) q[2];
sx q[2];
rz(-1.4359635) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0802059) q[1];
sx q[1];
rz(-2.5701984) q[1];
sx q[1];
rz(2.4953043) q[1];
rz(-pi) q[2];
rz(-0.16160139) q[3];
sx q[3];
rz(-2.9029663) q[3];
sx q[3];
rz(-2.0493226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7150813) q[2];
sx q[2];
rz(-2.0501523) q[2];
sx q[2];
rz(2.4439028) q[2];
rz(0.52608025) q[3];
sx q[3];
rz(-0.79280058) q[3];
sx q[3];
rz(1.5066159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0195011) q[0];
sx q[0];
rz(-2.2631336) q[0];
sx q[0];
rz(2.7418131) q[0];
rz(-1.1493692) q[1];
sx q[1];
rz(-0.72723564) q[1];
sx q[1];
rz(-0.74692187) q[1];
rz(-1.9952075) q[2];
sx q[2];
rz(-0.72327264) q[2];
sx q[2];
rz(-2.8930152) q[2];
rz(-1.2739194) q[3];
sx q[3];
rz(-1.0424725) q[3];
sx q[3];
rz(-2.7155664) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
