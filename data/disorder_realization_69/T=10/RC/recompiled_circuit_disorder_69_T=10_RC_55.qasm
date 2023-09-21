OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7944613) q[0];
sx q[0];
rz(-2.1262655) q[0];
sx q[0];
rz(2.6740958) q[0];
rz(-0.52019083) q[1];
sx q[1];
rz(-1.3462892) q[1];
sx q[1];
rz(2.4612114) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2517393) q[0];
sx q[0];
rz(-1.7008874) q[0];
sx q[0];
rz(-2.2962909) q[0];
rz(-pi) q[1];
rz(-1.8664076) q[2];
sx q[2];
rz(-1.3829074) q[2];
sx q[2];
rz(-2.0762028) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.80402741) q[1];
sx q[1];
rz(-2.0731888) q[1];
sx q[1];
rz(-2.875209) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3939788) q[3];
sx q[3];
rz(-1.2415981) q[3];
sx q[3];
rz(0.65229177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9900069) q[2];
sx q[2];
rz(-2.1534584) q[2];
sx q[2];
rz(3.0541259) q[2];
rz(0.72922373) q[3];
sx q[3];
rz(-2.7681523) q[3];
sx q[3];
rz(-0.27174404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4366348) q[0];
sx q[0];
rz(-2.3925662) q[0];
sx q[0];
rz(1.0013642) q[0];
rz(2.9691866) q[1];
sx q[1];
rz(-2.0253851) q[1];
sx q[1];
rz(2.6175245) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.090416106) q[0];
sx q[0];
rz(-2.2220082) q[0];
sx q[0];
rz(0.37474664) q[0];
rz(-0.42627724) q[2];
sx q[2];
rz(-2.2239416) q[2];
sx q[2];
rz(0.041235812) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.14428917) q[1];
sx q[1];
rz(-1.803949) q[1];
sx q[1];
rz(-0.38645978) q[1];
rz(-pi) q[2];
rz(-1.4161795) q[3];
sx q[3];
rz(-2.1025476) q[3];
sx q[3];
rz(-2.1173409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9841763) q[2];
sx q[2];
rz(-0.45980644) q[2];
sx q[2];
rz(-1.3400419) q[2];
rz(2.346758) q[3];
sx q[3];
rz(-2.0017616) q[3];
sx q[3];
rz(-0.036858233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3443417) q[0];
sx q[0];
rz(-0.69673711) q[0];
sx q[0];
rz(0.62477338) q[0];
rz(0.96427381) q[1];
sx q[1];
rz(-0.48502973) q[1];
sx q[1];
rz(-0.18951167) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7127842) q[0];
sx q[0];
rz(-1.3002987) q[0];
sx q[0];
rz(2.4012647) q[0];
x q[1];
rz(-2.8787829) q[2];
sx q[2];
rz(-2.7542979) q[2];
sx q[2];
rz(-2.5592321) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.0057919766) q[1];
sx q[1];
rz(-1.6264919) q[1];
sx q[1];
rz(0.23975753) q[1];
rz(-2.965851) q[3];
sx q[3];
rz(-0.97676859) q[3];
sx q[3];
rz(0.16080805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0345962) q[2];
sx q[2];
rz(-2.2980289) q[2];
sx q[2];
rz(-1.3872046) q[2];
rz(0.43131367) q[3];
sx q[3];
rz(-1.8547736) q[3];
sx q[3];
rz(1.0881902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65961924) q[0];
sx q[0];
rz(-0.94835931) q[0];
sx q[0];
rz(1.5959928) q[0];
rz(2.003147) q[1];
sx q[1];
rz(-2.3661416) q[1];
sx q[1];
rz(2.0514354) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2290303) q[0];
sx q[0];
rz(-1.1704485) q[0];
sx q[0];
rz(-1.1926665) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1074739) q[2];
sx q[2];
rz(-0.82674971) q[2];
sx q[2];
rz(2.2817051) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.412147) q[1];
sx q[1];
rz(-2.4033961) q[1];
sx q[1];
rz(-0.63936887) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7146043) q[3];
sx q[3];
rz(-1.7592351) q[3];
sx q[3];
rz(2.6915336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.99889341) q[2];
sx q[2];
rz(-0.45140758) q[2];
sx q[2];
rz(-2.2857655) q[2];
rz(-1.9479729) q[3];
sx q[3];
rz(-1.5193628) q[3];
sx q[3];
rz(2.2284171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49044931) q[0];
sx q[0];
rz(-2.1676846) q[0];
sx q[0];
rz(-0.14973101) q[0];
rz(0.99114746) q[1];
sx q[1];
rz(-1.9350355) q[1];
sx q[1];
rz(1.1700464) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67384185) q[0];
sx q[0];
rz(-1.2497328) q[0];
sx q[0];
rz(-2.2348316) q[0];
rz(-2.2511475) q[2];
sx q[2];
rz(-1.3157985) q[2];
sx q[2];
rz(-1.4444218) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7157117) q[1];
sx q[1];
rz(-0.60957805) q[1];
sx q[1];
rz(0.73519911) q[1];
rz(-pi) q[2];
rz(2.6036161) q[3];
sx q[3];
rz(-1.7378983) q[3];
sx q[3];
rz(-2.3218384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2698007) q[2];
sx q[2];
rz(-1.5912) q[2];
sx q[2];
rz(-3.0043547) q[2];
rz(-1.3373226) q[3];
sx q[3];
rz(-2.5604355) q[3];
sx q[3];
rz(1.2924682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9542434) q[0];
sx q[0];
rz(-1.7631148) q[0];
sx q[0];
rz(1.3504008) q[0];
rz(-0.84287914) q[1];
sx q[1];
rz(-0.73892361) q[1];
sx q[1];
rz(0.67289105) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56129365) q[0];
sx q[0];
rz(-1.6742047) q[0];
sx q[0];
rz(-2.1456477) q[0];
rz(-pi) q[1];
rz(-1.5422836) q[2];
sx q[2];
rz(-2.4270504) q[2];
sx q[2];
rz(-1.5039832) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7676047) q[1];
sx q[1];
rz(-1.7786221) q[1];
sx q[1];
rz(-2.5307104) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4937431) q[3];
sx q[3];
rz(-1.9814081) q[3];
sx q[3];
rz(-1.3239469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.792753) q[2];
sx q[2];
rz(-1.9323843) q[2];
sx q[2];
rz(2.7065281) q[2];
rz(1.3600291) q[3];
sx q[3];
rz(-2.3924148) q[3];
sx q[3];
rz(-2.8939261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.172794) q[0];
sx q[0];
rz(-0.89247576) q[0];
sx q[0];
rz(1.0621747) q[0];
rz(2.0299714) q[1];
sx q[1];
rz(-1.2373135) q[1];
sx q[1];
rz(-1.4020845) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37004334) q[0];
sx q[0];
rz(-1.2983027) q[0];
sx q[0];
rz(-0.38099654) q[0];
rz(-pi) q[1];
rz(-0.24537556) q[2];
sx q[2];
rz(-0.78005314) q[2];
sx q[2];
rz(1.1749554) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3558823) q[1];
sx q[1];
rz(-1.4175698) q[1];
sx q[1];
rz(2.7760837) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13774638) q[3];
sx q[3];
rz(-0.98689729) q[3];
sx q[3];
rz(0.13970845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4222251) q[2];
sx q[2];
rz(-2.5740467) q[2];
sx q[2];
rz(2.4613703) q[2];
rz(0.42823544) q[3];
sx q[3];
rz(-1.8869583) q[3];
sx q[3];
rz(-1.359882) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.632804) q[0];
sx q[0];
rz(-1.7254242) q[0];
sx q[0];
rz(1.592214) q[0];
rz(2.8920065) q[1];
sx q[1];
rz(-1.9892178) q[1];
sx q[1];
rz(-0.54135281) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2490059) q[0];
sx q[0];
rz(-1.5595946) q[0];
sx q[0];
rz(-2.5773994) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2301684) q[2];
sx q[2];
rz(-2.0001786) q[2];
sx q[2];
rz(-0.18061772) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9882422) q[1];
sx q[1];
rz(-1.9820947) q[1];
sx q[1];
rz(-2.3565355) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1369282) q[3];
sx q[3];
rz(-0.38032535) q[3];
sx q[3];
rz(2.5120146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.044518746) q[2];
sx q[2];
rz(-2.7842583) q[2];
sx q[2];
rz(-2.3642335) q[2];
rz(0.85123953) q[3];
sx q[3];
rz(-1.0612396) q[3];
sx q[3];
rz(1.8204934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38354307) q[0];
sx q[0];
rz(-1.8109011) q[0];
sx q[0];
rz(0.0099649075) q[0];
rz(1.0154356) q[1];
sx q[1];
rz(-2.3773057) q[1];
sx q[1];
rz(1.6962956) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2485219) q[0];
sx q[0];
rz(-2.0914441) q[0];
sx q[0];
rz(0.74670665) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.88405774) q[2];
sx q[2];
rz(-2.3381655) q[2];
sx q[2];
rz(-2.6510309) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6362793) q[1];
sx q[1];
rz(-1.5552102) q[1];
sx q[1];
rz(-0.63025766) q[1];
rz(-pi) q[2];
rz(2.2441838) q[3];
sx q[3];
rz(-1.0703147) q[3];
sx q[3];
rz(-0.15411479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4426667) q[2];
sx q[2];
rz(-2.2103504) q[2];
sx q[2];
rz(0.60738579) q[2];
rz(-1.7025042) q[3];
sx q[3];
rz(-1.7581698) q[3];
sx q[3];
rz(-0.79090345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33567515) q[0];
sx q[0];
rz(-1.6573925) q[0];
sx q[0];
rz(2.5277396) q[0];
rz(1.0461668) q[1];
sx q[1];
rz(-0.26509735) q[1];
sx q[1];
rz(2.7493431) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0948254) q[0];
sx q[0];
rz(-1.9691159) q[0];
sx q[0];
rz(2.4368068) q[0];
rz(-0.73108436) q[2];
sx q[2];
rz(-0.77027551) q[2];
sx q[2];
rz(1.1243656) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0872005) q[1];
sx q[1];
rz(-2.0093845) q[1];
sx q[1];
rz(-0.88263504) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5997821) q[3];
sx q[3];
rz(-1.7757123) q[3];
sx q[3];
rz(-0.25587413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0739416) q[2];
sx q[2];
rz(-1.3839046) q[2];
sx q[2];
rz(-2.5411141) q[2];
rz(2.0742119) q[3];
sx q[3];
rz(-1.3200656) q[3];
sx q[3];
rz(-1.0935121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4338715) q[0];
sx q[0];
rz(-0.43294551) q[0];
sx q[0];
rz(-1.659163) q[0];
rz(-2.1451163) q[1];
sx q[1];
rz(-1.4947718) q[1];
sx q[1];
rz(1.6047118) q[1];
rz(-1.9267351) q[2];
sx q[2];
rz(-0.83649737) q[2];
sx q[2];
rz(2.5964824) q[2];
rz(-2.2157833) q[3];
sx q[3];
rz(-0.79173294) q[3];
sx q[3];
rz(-1.790495) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
