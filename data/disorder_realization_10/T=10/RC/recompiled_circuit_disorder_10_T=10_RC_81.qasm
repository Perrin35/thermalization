OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7005641) q[0];
sx q[0];
rz(1.1428042) q[0];
sx q[0];
rz(11.354843) q[0];
rz(2.9149574) q[1];
sx q[1];
rz(-1.5645138) q[1];
sx q[1];
rz(-0.29830631) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.053781833) q[0];
sx q[0];
rz(-0.90667533) q[0];
sx q[0];
rz(1.5009297) q[0];
rz(-1.83625) q[2];
sx q[2];
rz(-2.7976118) q[2];
sx q[2];
rz(-1.6814107) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5391985) q[1];
sx q[1];
rz(-1.0326003) q[1];
sx q[1];
rz(-0.56166517) q[1];
rz(-0.084007752) q[3];
sx q[3];
rz(-1.6735895) q[3];
sx q[3];
rz(-1.8889129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6478708) q[2];
sx q[2];
rz(-1.2080668) q[2];
sx q[2];
rz(0.99386627) q[2];
rz(2.1422051) q[3];
sx q[3];
rz(-1.9013654) q[3];
sx q[3];
rz(-2.5527111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64269972) q[0];
sx q[0];
rz(-1.2058586) q[0];
sx q[0];
rz(-3.1233741) q[0];
rz(2.3253564) q[1];
sx q[1];
rz(-1.0304334) q[1];
sx q[1];
rz(0.47168628) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2361006) q[0];
sx q[0];
rz(-0.24370757) q[0];
sx q[0];
rz(-2.7729176) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8505487) q[2];
sx q[2];
rz(-1.0824167) q[2];
sx q[2];
rz(2.5807057) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.43863152) q[1];
sx q[1];
rz(-1.303181) q[1];
sx q[1];
rz(0.57499927) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96305965) q[3];
sx q[3];
rz(-1.4415603) q[3];
sx q[3];
rz(-2.3605763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.54962426) q[2];
sx q[2];
rz(-1.1529808) q[2];
sx q[2];
rz(-2.1726051) q[2];
rz(-0.5747059) q[3];
sx q[3];
rz(-0.55137268) q[3];
sx q[3];
rz(-1.0415174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4330924) q[0];
sx q[0];
rz(-2.0860465) q[0];
sx q[0];
rz(2.95978) q[0];
rz(-2.0388942) q[1];
sx q[1];
rz(-1.6405374) q[1];
sx q[1];
rz(-1.6859432) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11765471) q[0];
sx q[0];
rz(-2.0485989) q[0];
sx q[0];
rz(-2.2295203) q[0];
rz(-pi) q[1];
rz(0.74080148) q[2];
sx q[2];
rz(-1.891279) q[2];
sx q[2];
rz(-1.6522811) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.94565832) q[1];
sx q[1];
rz(-0.12609005) q[1];
sx q[1];
rz(-0.083421589) q[1];
x q[2];
rz(2.5163243) q[3];
sx q[3];
rz(-1.6571952) q[3];
sx q[3];
rz(2.3141935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.87749798) q[2];
sx q[2];
rz(-1.4870746) q[2];
sx q[2];
rz(0.96763119) q[2];
rz(-2.4140221) q[3];
sx q[3];
rz(-1.8811767) q[3];
sx q[3];
rz(2.9038866) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85686344) q[0];
sx q[0];
rz(-2.6155222) q[0];
sx q[0];
rz(2.5033584) q[0];
rz(1.1278641) q[1];
sx q[1];
rz(-0.82740873) q[1];
sx q[1];
rz(-1.2329873) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5467984) q[0];
sx q[0];
rz(-1.5017171) q[0];
sx q[0];
rz(-2.6462206) q[0];
rz(-pi) q[1];
rz(-2.4013176) q[2];
sx q[2];
rz(-0.7358272) q[2];
sx q[2];
rz(-2.1528113) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2265046) q[1];
sx q[1];
rz(-2.9595032) q[1];
sx q[1];
rz(-1.7712797) q[1];
rz(-pi) q[2];
rz(0.87807699) q[3];
sx q[3];
rz(-1.3948166) q[3];
sx q[3];
rz(0.39719492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2234852) q[2];
sx q[2];
rz(-1.2780259) q[2];
sx q[2];
rz(-0.81400648) q[2];
rz(-2.0984086) q[3];
sx q[3];
rz(-0.63101763) q[3];
sx q[3];
rz(-1.1842747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2994613) q[0];
sx q[0];
rz(-1.3548387) q[0];
sx q[0];
rz(0.89170757) q[0];
rz(1.2437598) q[1];
sx q[1];
rz(-1.3777106) q[1];
sx q[1];
rz(-0.2125425) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71056238) q[0];
sx q[0];
rz(-3.115603) q[0];
sx q[0];
rz(2.2463069) q[0];
rz(-pi) q[1];
rz(1.6264621) q[2];
sx q[2];
rz(-2.5362483) q[2];
sx q[2];
rz(0.84673131) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6032431) q[1];
sx q[1];
rz(-1.0298567) q[1];
sx q[1];
rz(-2.4451838) q[1];
rz(-pi) q[2];
rz(0.84380031) q[3];
sx q[3];
rz(-1.2146597) q[3];
sx q[3];
rz(0.86405495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1084958) q[2];
sx q[2];
rz(-0.97110811) q[2];
sx q[2];
rz(-0.5212211) q[2];
rz(1.3850348) q[3];
sx q[3];
rz(-1.228046) q[3];
sx q[3];
rz(-1.8732171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2824771) q[0];
sx q[0];
rz(-0.92882597) q[0];
sx q[0];
rz(-0.22512063) q[0];
rz(-1.3549995) q[1];
sx q[1];
rz(-2.1332108) q[1];
sx q[1];
rz(0.37757847) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7901944) q[0];
sx q[0];
rz(-1.5927918) q[0];
sx q[0];
rz(-1.8368594) q[0];
x q[1];
rz(0.2874561) q[2];
sx q[2];
rz(-1.6171347) q[2];
sx q[2];
rz(3.1099144) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4695417) q[1];
sx q[1];
rz(-0.44476032) q[1];
sx q[1];
rz(1.8163535) q[1];
rz(-pi) q[2];
rz(-2.3466831) q[3];
sx q[3];
rz(-1.5494293) q[3];
sx q[3];
rz(1.4710609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1569415) q[2];
sx q[2];
rz(-0.7876544) q[2];
sx q[2];
rz(2.6605576) q[2];
rz(-0.40361079) q[3];
sx q[3];
rz(-1.0389046) q[3];
sx q[3];
rz(2.8267982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.0890546) q[0];
sx q[0];
rz(-0.60482329) q[0];
sx q[0];
rz(-0.19454923) q[0];
rz(-0.21952195) q[1];
sx q[1];
rz(-1.4621282) q[1];
sx q[1];
rz(2.887168) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2339904) q[0];
sx q[0];
rz(-1.4617209) q[0];
sx q[0];
rz(-1.2489737) q[0];
rz(2.1830325) q[2];
sx q[2];
rz(-1.6974291) q[2];
sx q[2];
rz(-0.74778344) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.88528819) q[1];
sx q[1];
rz(-1.8418152) q[1];
sx q[1];
rz(-0.46896743) q[1];
rz(-pi) q[2];
rz(2.3354704) q[3];
sx q[3];
rz(-2.5297909) q[3];
sx q[3];
rz(2.7968189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.97757942) q[2];
sx q[2];
rz(-1.3914725) q[2];
sx q[2];
rz(-1.8257726) q[2];
rz(-2.2655462) q[3];
sx q[3];
rz(-3.0026569) q[3];
sx q[3];
rz(1.0036489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9309689) q[0];
sx q[0];
rz(-2.7656778) q[0];
sx q[0];
rz(1.4550495) q[0];
rz(-2.3176106) q[1];
sx q[1];
rz(-2.1897557) q[1];
sx q[1];
rz(1.5751858) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9109089) q[0];
sx q[0];
rz(-2.3678603) q[0];
sx q[0];
rz(0.45549972) q[0];
rz(-pi) q[1];
rz(-1.8640395) q[2];
sx q[2];
rz(-0.55985057) q[2];
sx q[2];
rz(0.44550371) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7745061) q[1];
sx q[1];
rz(-1.204406) q[1];
sx q[1];
rz(-0.3831425) q[1];
rz(-pi) q[2];
rz(-1.1493707) q[3];
sx q[3];
rz(-2.0391658) q[3];
sx q[3];
rz(-0.8120265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.87551293) q[2];
sx q[2];
rz(-1.9078887) q[2];
sx q[2];
rz(2.8640462) q[2];
rz(-1.4510441) q[3];
sx q[3];
rz(-0.45193672) q[3];
sx q[3];
rz(2.1267166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16185109) q[0];
sx q[0];
rz(-0.45270544) q[0];
sx q[0];
rz(-1.45654) q[0];
rz(0.62943554) q[1];
sx q[1];
rz(-1.9742191) q[1];
sx q[1];
rz(-1.1368407) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0757383) q[0];
sx q[0];
rz(-0.59959164) q[0];
sx q[0];
rz(1.0435186) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7669737) q[2];
sx q[2];
rz(-1.9667452) q[2];
sx q[2];
rz(-1.7916726) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7965664) q[1];
sx q[1];
rz(-1.6654286) q[1];
sx q[1];
rz(-1.8887397) q[1];
rz(-pi) q[2];
rz(1.9433446) q[3];
sx q[3];
rz(-2.2714943) q[3];
sx q[3];
rz(2.8230599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.64951605) q[2];
sx q[2];
rz(-1.4294383) q[2];
sx q[2];
rz(-2.1006404) q[2];
rz(-3.1395636) q[3];
sx q[3];
rz(-2.7414331) q[3];
sx q[3];
rz(-2.8295529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3025538) q[0];
sx q[0];
rz(-2.9170687) q[0];
sx q[0];
rz(-0.94605207) q[0];
rz(2.229915) q[1];
sx q[1];
rz(-1.9263575) q[1];
sx q[1];
rz(-0.61202234) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6291954) q[0];
sx q[0];
rz(-1.1149659) q[0];
sx q[0];
rz(-0.7645316) q[0];
rz(-pi) q[1];
rz(0.34531784) q[2];
sx q[2];
rz(-1.7974263) q[2];
sx q[2];
rz(-2.741284) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8814197) q[1];
sx q[1];
rz(-2.2248785) q[1];
sx q[1];
rz(-1.7263078) q[1];
rz(-0.58935921) q[3];
sx q[3];
rz(-0.61925626) q[3];
sx q[3];
rz(-0.15976957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0845906) q[2];
sx q[2];
rz(-0.64208639) q[2];
sx q[2];
rz(2.4882312) q[2];
rz(0.35081321) q[3];
sx q[3];
rz(-1.5272798) q[3];
sx q[3];
rz(-2.4408834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6012797) q[0];
sx q[0];
rz(-1.5355587) q[0];
sx q[0];
rz(-3.0328947) q[0];
rz(-0.75469771) q[1];
sx q[1];
rz(-1.3194059) q[1];
sx q[1];
rz(-1.5059765) q[1];
rz(-2.1379708) q[2];
sx q[2];
rz(-0.97273172) q[2];
sx q[2];
rz(1.6956971) q[2];
rz(-1.894542) q[3];
sx q[3];
rz(-1.4907881) q[3];
sx q[3];
rz(-1.1992906) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
