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
rz(0.29851222) q[0];
sx q[0];
rz(-1.8520344) q[0];
sx q[0];
rz(0.02331743) q[0];
rz(-2.159637) q[1];
sx q[1];
rz(7.0206849) q[1];
sx q[1];
rz(7.7590396) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6029089) q[0];
sx q[0];
rz(-1.6404466) q[0];
sx q[0];
rz(-3.0979741) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9010004) q[2];
sx q[2];
rz(-1.6088952) q[2];
sx q[2];
rz(0.35534258) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6827794) q[1];
sx q[1];
rz(-2.2662357) q[1];
sx q[1];
rz(2.0942445) q[1];
x q[2];
rz(-2.2470583) q[3];
sx q[3];
rz(-1.1522376) q[3];
sx q[3];
rz(-0.422131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.6797356) q[2];
sx q[2];
rz(-1.6697786) q[2];
sx q[2];
rz(2.8409345) q[2];
rz(-2.4833208) q[3];
sx q[3];
rz(-0.39977795) q[3];
sx q[3];
rz(-2.5532653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7268426) q[0];
sx q[0];
rz(-0.86597935) q[0];
sx q[0];
rz(0.17587371) q[0];
rz(-0.56100065) q[1];
sx q[1];
rz(-2.439552) q[1];
sx q[1];
rz(0.19634136) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4787653) q[0];
sx q[0];
rz(-1.2365554) q[0];
sx q[0];
rz(-2.0408551) q[0];
x q[1];
rz(0.11949338) q[2];
sx q[2];
rz(-1.6867078) q[2];
sx q[2];
rz(1.8279861) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2113116) q[1];
sx q[1];
rz(-1.0553816) q[1];
sx q[1];
rz(-0.0023018259) q[1];
rz(-pi) q[2];
rz(-1.4370055) q[3];
sx q[3];
rz(-1.9113298) q[3];
sx q[3];
rz(1.2029369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0520797) q[2];
sx q[2];
rz(-2.5255346) q[2];
sx q[2];
rz(1.4698131) q[2];
rz(-2.6164264) q[3];
sx q[3];
rz(-1.4273806) q[3];
sx q[3];
rz(2.4770881) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1216275) q[0];
sx q[0];
rz(-2.2522734) q[0];
sx q[0];
rz(-0.63051939) q[0];
rz(1.1446674) q[1];
sx q[1];
rz(-1.8960709) q[1];
sx q[1];
rz(-0.05680457) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6358179) q[0];
sx q[0];
rz(-2.5774101) q[0];
sx q[0];
rz(0.83585541) q[0];
rz(-pi) q[1];
rz(1.0089097) q[2];
sx q[2];
rz(-0.76912921) q[2];
sx q[2];
rz(-3.1388856) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9020128) q[1];
sx q[1];
rz(-2.6551592) q[1];
sx q[1];
rz(2.835266) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1278387) q[3];
sx q[3];
rz(-2.6747799) q[3];
sx q[3];
rz(-1.6396322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8593665) q[2];
sx q[2];
rz(-1.4713918) q[2];
sx q[2];
rz(2.7024506) q[2];
rz(1.3965083) q[3];
sx q[3];
rz(-1.8789995) q[3];
sx q[3];
rz(2.5631574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0466995) q[0];
sx q[0];
rz(-2.4201604) q[0];
sx q[0];
rz(2.1266345) q[0];
rz(3.0789442) q[1];
sx q[1];
rz(-0.95476127) q[1];
sx q[1];
rz(-2.8573724) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13260176) q[0];
sx q[0];
rz(-1.9997445) q[0];
sx q[0];
rz(-1.1852077) q[0];
x q[1];
rz(1.7106871) q[2];
sx q[2];
rz(-1.0770742) q[2];
sx q[2];
rz(1.7715275) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9584459) q[1];
sx q[1];
rz(-1.5588405) q[1];
sx q[1];
rz(-0.35838106) q[1];
rz(1.5169859) q[3];
sx q[3];
rz(-0.25370666) q[3];
sx q[3];
rz(1.0347317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3011498) q[2];
sx q[2];
rz(-2.132685) q[2];
sx q[2];
rz(0.76048771) q[2];
rz(2.8216951) q[3];
sx q[3];
rz(-2.0127998) q[3];
sx q[3];
rz(-1.0285876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8496721) q[0];
sx q[0];
rz(-0.55067647) q[0];
sx q[0];
rz(-0.40529761) q[0];
rz(0.48794508) q[1];
sx q[1];
rz(-1.1064203) q[1];
sx q[1];
rz(-2.778756) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3317446) q[0];
sx q[0];
rz(-0.54143006) q[0];
sx q[0];
rz(2.5141513) q[0];
rz(-pi) q[1];
rz(1.1929197) q[2];
sx q[2];
rz(-2.123017) q[2];
sx q[2];
rz(-0.044980031) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2606973) q[1];
sx q[1];
rz(-0.22628566) q[1];
sx q[1];
rz(-0.21842167) q[1];
rz(-1.8764349) q[3];
sx q[3];
rz(-0.89857093) q[3];
sx q[3];
rz(1.5346988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.001658) q[2];
sx q[2];
rz(-1.1589061) q[2];
sx q[2];
rz(0.72251594) q[2];
rz(-0.51269382) q[3];
sx q[3];
rz(-2.2721458) q[3];
sx q[3];
rz(1.2041913) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8743415) q[0];
sx q[0];
rz(-0.63684547) q[0];
sx q[0];
rz(-2.8668218) q[0];
rz(0.35708669) q[1];
sx q[1];
rz(-1.5029224) q[1];
sx q[1];
rz(1.9836609) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5362893) q[0];
sx q[0];
rz(-0.9035631) q[0];
sx q[0];
rz(-2.8298004) q[0];
x q[1];
rz(0.91041301) q[2];
sx q[2];
rz(-1.0779625) q[2];
sx q[2];
rz(-3.1255258) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.14640644) q[1];
sx q[1];
rz(-0.9764834) q[1];
sx q[1];
rz(-1.3239149) q[1];
rz(2.3401392) q[3];
sx q[3];
rz(-1.0545925) q[3];
sx q[3];
rz(0.046122915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.57924119) q[2];
sx q[2];
rz(-1.3902105) q[2];
sx q[2];
rz(0.58970279) q[2];
rz(1.8288745) q[3];
sx q[3];
rz(-1.4785942) q[3];
sx q[3];
rz(2.5780799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(0.24484672) q[0];
sx q[0];
rz(-2.1488996) q[0];
sx q[0];
rz(2.8712811) q[0];
rz(-2.7145794) q[1];
sx q[1];
rz(-1.7662363) q[1];
sx q[1];
rz(2.7332773) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9664551) q[0];
sx q[0];
rz(-1.2027367) q[0];
sx q[0];
rz(-1.4428663) q[0];
rz(-2.2683581) q[2];
sx q[2];
rz(-1.4931275) q[2];
sx q[2];
rz(2.6447062) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.282849) q[1];
sx q[1];
rz(-2.3925684) q[1];
sx q[1];
rz(-0.94094679) q[1];
rz(-pi) q[2];
rz(0.72260489) q[3];
sx q[3];
rz(-1.4935023) q[3];
sx q[3];
rz(-2.3376436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2030877) q[2];
sx q[2];
rz(-1.0838584) q[2];
sx q[2];
rz(-2.2173524) q[2];
rz(-0.83485323) q[3];
sx q[3];
rz(-0.82930851) q[3];
sx q[3];
rz(2.1448962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.1860109) q[0];
sx q[0];
rz(-2.7614433) q[0];
sx q[0];
rz(-2.7741449) q[0];
rz(2.5732749) q[1];
sx q[1];
rz(-1.808993) q[1];
sx q[1];
rz(2.7265991) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85431722) q[0];
sx q[0];
rz(-2.0234031) q[0];
sx q[0];
rz(0.29083473) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0239564) q[2];
sx q[2];
rz(-2.0342361) q[2];
sx q[2];
rz(1.6187504) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4492682) q[1];
sx q[1];
rz(-2.0776761) q[1];
sx q[1];
rz(-2.436422) q[1];
x q[2];
rz(3.1363963) q[3];
sx q[3];
rz(-1.9165766) q[3];
sx q[3];
rz(0.41692641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7830398) q[2];
sx q[2];
rz(-2.2542605) q[2];
sx q[2];
rz(-2.3084194) q[2];
rz(-0.35946515) q[3];
sx q[3];
rz(-2.0514252) q[3];
sx q[3];
rz(1.5045213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9678765) q[0];
sx q[0];
rz(-2.739527) q[0];
sx q[0];
rz(-0.014884431) q[0];
rz(2.0170276) q[1];
sx q[1];
rz(-0.24528565) q[1];
sx q[1];
rz(-0.10717779) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10725645) q[0];
sx q[0];
rz(-1.4217303) q[0];
sx q[0];
rz(-2.3019522) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6769019) q[2];
sx q[2];
rz(-2.26311) q[2];
sx q[2];
rz(-1.2274881) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9669144) q[1];
sx q[1];
rz(-1.0941097) q[1];
sx q[1];
rz(-1.7382453) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.26132432) q[3];
sx q[3];
rz(-1.3991941) q[3];
sx q[3];
rz(-1.2444519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5643481) q[2];
sx q[2];
rz(-1.3078835) q[2];
sx q[2];
rz(2.2819819) q[2];
rz(2.882242) q[3];
sx q[3];
rz(-1.7813464) q[3];
sx q[3];
rz(-2.7249961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6421826) q[0];
sx q[0];
rz(-0.12633093) q[0];
sx q[0];
rz(1.9835749) q[0];
rz(-0.49531373) q[1];
sx q[1];
rz(-1.1664349) q[1];
sx q[1];
rz(2.8445629) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7659781) q[0];
sx q[0];
rz(-0.31904623) q[0];
sx q[0];
rz(-1.5578397) q[0];
x q[1];
rz(-3.0115602) q[2];
sx q[2];
rz(-2.8388925) q[2];
sx q[2];
rz(-1.5720194) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2386855) q[1];
sx q[1];
rz(-2.8082962) q[1];
sx q[1];
rz(-2.6413647) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2333651) q[3];
sx q[3];
rz(-1.5836827) q[3];
sx q[3];
rz(-2.9290309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.872252) q[2];
sx q[2];
rz(-2.7597235) q[2];
sx q[2];
rz(-0.86088172) q[2];
rz(2.6779029) q[3];
sx q[3];
rz(-0.90353614) q[3];
sx q[3];
rz(1.1369107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.1562445) q[0];
sx q[0];
rz(-1.5679659) q[0];
sx q[0];
rz(1.1563942) q[0];
rz(-1.8078177) q[1];
sx q[1];
rz(-1.6009686) q[1];
sx q[1];
rz(0.70645465) q[1];
rz(-2.1352519) q[2];
sx q[2];
rz(-2.153572) q[2];
sx q[2];
rz(-1.5700271) q[2];
rz(0.64846497) q[3];
sx q[3];
rz(-0.59352661) q[3];
sx q[3];
rz(1.2980579) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
