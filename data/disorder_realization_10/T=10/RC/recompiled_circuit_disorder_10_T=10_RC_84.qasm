OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4410285) q[0];
sx q[0];
rz(-1.1428042) q[0];
sx q[0];
rz(-1.2115275) q[0];
rz(2.9149574) q[1];
sx q[1];
rz(-1.5645138) q[1];
sx q[1];
rz(-0.29830631) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.053781833) q[0];
sx q[0];
rz(-2.2349173) q[0];
sx q[0];
rz(1.640663) q[0];
x q[1];
rz(1.9036129) q[2];
sx q[2];
rz(-1.6593854) q[2];
sx q[2];
rz(0.36117902) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.71516192) q[1];
sx q[1];
rz(-2.3843345) q[1];
sx q[1];
rz(0.84233474) q[1];
rz(-1.467642) q[3];
sx q[3];
rz(-1.6543596) q[3];
sx q[3];
rz(-0.32675693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4937218) q[2];
sx q[2];
rz(-1.9335258) q[2];
sx q[2];
rz(-2.1477264) q[2];
rz(2.1422051) q[3];
sx q[3];
rz(-1.9013654) q[3];
sx q[3];
rz(0.58888155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64269972) q[0];
sx q[0];
rz(-1.2058586) q[0];
sx q[0];
rz(-0.018218536) q[0];
rz(-2.3253564) q[1];
sx q[1];
rz(-2.1111592) q[1];
sx q[1];
rz(-2.6699064) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52662151) q[0];
sx q[0];
rz(-1.7978298) q[0];
sx q[0];
rz(-1.660166) q[0];
rz(-pi) q[1];
rz(1.8505487) q[2];
sx q[2];
rz(-1.0824167) q[2];
sx q[2];
rz(-2.5807057) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.43863152) q[1];
sx q[1];
rz(-1.303181) q[1];
sx q[1];
rz(-2.5665934) q[1];
rz(-pi) q[2];
rz(-0.96305965) q[3];
sx q[3];
rz(-1.7000323) q[3];
sx q[3];
rz(-2.3605763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5919684) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4330924) q[0];
sx q[0];
rz(-1.0555462) q[0];
sx q[0];
rz(0.18181268) q[0];
rz(2.0388942) q[1];
sx q[1];
rz(-1.5010553) q[1];
sx q[1];
rz(-1.6859432) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2244814) q[0];
sx q[0];
rz(-2.3492976) q[0];
sx q[0];
rz(2.2729421) q[0];
rz(-pi) q[1];
rz(-2.6845034) q[2];
sx q[2];
rz(-2.3466913) q[2];
sx q[2];
rz(2.7283816) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1118483) q[1];
sx q[1];
rz(-1.4451471) q[1];
sx q[1];
rz(1.5813584) q[1];
x q[2];
rz(2.5163243) q[3];
sx q[3];
rz(-1.6571952) q[3];
sx q[3];
rz(-0.82739917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.87749798) q[2];
sx q[2];
rz(-1.4870746) q[2];
sx q[2];
rz(0.96763119) q[2];
rz(0.72757059) q[3];
sx q[3];
rz(-1.8811767) q[3];
sx q[3];
rz(2.9038866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.2847292) q[0];
sx q[0];
rz(-2.6155222) q[0];
sx q[0];
rz(0.63823429) q[0];
rz(-1.1278641) q[1];
sx q[1];
rz(-0.82740873) q[1];
sx q[1];
rz(1.2329873) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1283135) q[0];
sx q[0];
rz(-1.0767125) q[0];
sx q[0];
rz(-1.6492776) q[0];
x q[1];
rz(-2.5523283) q[2];
sx q[2];
rz(-1.1009842) q[2];
sx q[2];
rz(-1.1772917) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2265046) q[1];
sx q[1];
rz(-2.9595032) q[1];
sx q[1];
rz(-1.7712797) q[1];
x q[2];
rz(-0.22709417) q[3];
sx q[3];
rz(-2.2507651) q[3];
sx q[3];
rz(-1.029315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.91810742) q[2];
sx q[2];
rz(-1.2780259) q[2];
sx q[2];
rz(-2.3275862) q[2];
rz(-2.0984086) q[3];
sx q[3];
rz(-0.63101763) q[3];
sx q[3];
rz(1.957318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84213132) q[0];
sx q[0];
rz(-1.786754) q[0];
sx q[0];
rz(0.89170757) q[0];
rz(-1.8978329) q[1];
sx q[1];
rz(-1.3777106) q[1];
sx q[1];
rz(-0.2125425) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.034886995) q[0];
sx q[0];
rz(-1.5910774) q[0];
sx q[0];
rz(3.1253392) q[0];
rz(-pi) q[1];
rz(0.96617713) q[2];
sx q[2];
rz(-1.6024616) q[2];
sx q[2];
rz(0.67827536) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.58442851) q[1];
sx q[1];
rz(-2.2884528) q[1];
sx q[1];
rz(-2.3889956) q[1];
rz(-pi) q[2];
rz(-0.46194525) q[3];
sx q[3];
rz(-2.2432703) q[3];
sx q[3];
rz(2.7355821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.033096878) q[2];
sx q[2];
rz(-2.1704845) q[2];
sx q[2];
rz(2.6203716) q[2];
rz(-1.7565578) q[3];
sx q[3];
rz(-1.9135467) q[3];
sx q[3];
rz(1.8732171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8591156) q[0];
sx q[0];
rz(-2.2127667) q[0];
sx q[0];
rz(-2.916472) q[0];
rz(-1.3549995) q[1];
sx q[1];
rz(-1.0083818) q[1];
sx q[1];
rz(-0.37757847) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3513983) q[0];
sx q[0];
rz(-1.5927918) q[0];
sx q[0];
rz(-1.8368594) q[0];
x q[1];
rz(-0.16212459) q[2];
sx q[2];
rz(-2.8505278) q[2];
sx q[2];
rz(1.7578917) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0203591) q[1];
sx q[1];
rz(-1.6755783) q[1];
sx q[1];
rz(-2.0038414) q[1];
rz(-0.029929786) q[3];
sx q[3];
rz(-0.79513351) q[3];
sx q[3];
rz(-3.0628169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1569415) q[2];
sx q[2];
rz(-0.7876544) q[2];
sx q[2];
rz(-0.48103508) q[2];
rz(-2.7379819) q[3];
sx q[3];
rz(-2.1026881) q[3];
sx q[3];
rz(2.8267982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(1.0890546) q[0];
sx q[0];
rz(-0.60482329) q[0];
sx q[0];
rz(-2.9470434) q[0];
rz(0.21952195) q[1];
sx q[1];
rz(-1.6794645) q[1];
sx q[1];
rz(2.887168) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90760224) q[0];
sx q[0];
rz(-1.4617209) q[0];
sx q[0];
rz(-1.8926189) q[0];
x q[1];
rz(1.7888072) q[2];
sx q[2];
rz(-0.62354747) q[2];
sx q[2];
rz(2.4965198) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.82032694) q[1];
sx q[1];
rz(-1.1202381) q[1];
sx q[1];
rz(-1.8727559) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0394578) q[3];
sx q[3];
rz(-1.9797167) q[3];
sx q[3];
rz(-1.2498145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1640132) q[2];
sx q[2];
rz(-1.3914725) q[2];
sx q[2];
rz(-1.8257726) q[2];
rz(-2.2655462) q[3];
sx q[3];
rz(-0.13893572) q[3];
sx q[3];
rz(2.1379437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9309689) q[0];
sx q[0];
rz(-0.3759149) q[0];
sx q[0];
rz(-1.4550495) q[0];
rz(2.3176106) q[1];
sx q[1];
rz(-0.95183698) q[1];
sx q[1];
rz(-1.5664068) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5113735) q[0];
sx q[0];
rz(-2.2492118) q[0];
sx q[0];
rz(-1.1648965) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1112061) q[2];
sx q[2];
rz(-1.724913) q[2];
sx q[2];
rz(0.87481462) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7944784) q[1];
sx q[1];
rz(-1.2142668) q[1];
sx q[1];
rz(-1.1785248) q[1];
rz(-pi) q[2];
rz(-0.50623399) q[3];
sx q[3];
rz(-1.1971548) q[3];
sx q[3];
rz(-0.95844275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.87551293) q[2];
sx q[2];
rz(-1.2337039) q[2];
sx q[2];
rz(-2.8640462) q[2];
rz(1.4510441) q[3];
sx q[3];
rz(-0.45193672) q[3];
sx q[3];
rz(1.014876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.9797416) q[0];
sx q[0];
rz(-2.6888872) q[0];
sx q[0];
rz(-1.45654) q[0];
rz(0.62943554) q[1];
sx q[1];
rz(-1.9742191) q[1];
sx q[1];
rz(2.004752) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0946227) q[0];
sx q[0];
rz(-1.8587062) q[0];
sx q[0];
rz(-1.0372435) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.37461899) q[2];
sx q[2];
rz(-1.1748474) q[2];
sx q[2];
rz(-1.7916726) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7965664) q[1];
sx q[1];
rz(-1.6654286) q[1];
sx q[1];
rz(-1.8887397) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7342019) q[3];
sx q[3];
rz(-2.3630777) q[3];
sx q[3];
rz(-0.22637573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.64951605) q[2];
sx q[2];
rz(-1.7121544) q[2];
sx q[2];
rz(2.1006404) q[2];
rz(3.1395636) q[3];
sx q[3];
rz(-0.40015951) q[3];
sx q[3];
rz(-2.8295529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8390389) q[0];
sx q[0];
rz(-0.22452393) q[0];
sx q[0];
rz(-0.94605207) q[0];
rz(2.229915) q[1];
sx q[1];
rz(-1.9263575) q[1];
sx q[1];
rz(2.5295703) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6529022) q[0];
sx q[0];
rz(-2.2757747) q[0];
sx q[0];
rz(-2.5253354) q[0];
rz(1.3304747) q[2];
sx q[2];
rz(-1.2346621) q[2];
sx q[2];
rz(1.0898332) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.40572383) q[1];
sx q[1];
rz(-1.4475665) q[1];
sx q[1];
rz(0.65995364) q[1];
x q[2];
rz(1.1935812) q[3];
sx q[3];
rz(-2.0743138) q[3];
sx q[3];
rz(-2.6138888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0570021) q[2];
sx q[2];
rz(-0.64208639) q[2];
sx q[2];
rz(2.4882312) q[2];
rz(0.35081321) q[3];
sx q[3];
rz(-1.5272798) q[3];
sx q[3];
rz(0.70070926) q[3];
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
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54031298) q[0];
sx q[0];
rz(-1.606034) q[0];
sx q[0];
rz(0.10869797) q[0];
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
