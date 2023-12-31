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
rz(-0.46749687) q[0];
rz(-0.52019083) q[1];
sx q[1];
rz(1.7953035) q[1];
sx q[1];
rz(10.105159) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43359783) q[0];
sx q[0];
rz(-2.2888219) q[0];
sx q[0];
rz(2.9684767) q[0];
rz(-pi) q[1];
rz(-2.1490287) q[2];
sx q[2];
rz(-0.34878584) q[2];
sx q[2];
rz(-1.055583) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.28847028) q[1];
sx q[1];
rz(-2.5783357) q[1];
sx q[1];
rz(-1.1239777) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.666113) q[3];
sx q[3];
rz(-2.7694422) q[3];
sx q[3];
rz(-1.1572157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.15158571) q[2];
sx q[2];
rz(-0.9881343) q[2];
sx q[2];
rz(-0.087466784) q[2];
rz(-0.72922373) q[3];
sx q[3];
rz(-0.37344033) q[3];
sx q[3];
rz(2.8698486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7049578) q[0];
sx q[0];
rz(-2.3925662) q[0];
sx q[0];
rz(-2.1402284) q[0];
rz(0.17240605) q[1];
sx q[1];
rz(-2.0253851) q[1];
sx q[1];
rz(0.52406812) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0511765) q[0];
sx q[0];
rz(-0.91958445) q[0];
sx q[0];
rz(0.37474664) q[0];
x q[1];
rz(-0.87191138) q[2];
sx q[2];
rz(-1.9053835) q[2];
sx q[2];
rz(1.3427693) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3327649) q[1];
sx q[1];
rz(-1.9462703) q[1];
sx q[1];
rz(1.8217702) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.53701138) q[3];
sx q[3];
rz(-1.7039263) q[3];
sx q[3];
rz(2.6739124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9841763) q[2];
sx q[2];
rz(-2.6817862) q[2];
sx q[2];
rz(-1.3400419) q[2];
rz(-0.79483461) q[3];
sx q[3];
rz(-2.0017616) q[3];
sx q[3];
rz(3.1047344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79725093) q[0];
sx q[0];
rz(-0.69673711) q[0];
sx q[0];
rz(0.62477338) q[0];
rz(-2.1773188) q[1];
sx q[1];
rz(-0.48502973) q[1];
sx q[1];
rz(2.952081) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0441168) q[0];
sx q[0];
rz(-2.2783845) q[0];
sx q[0];
rz(-1.9301027) q[0];
rz(-1.6763716) q[2];
sx q[2];
rz(-1.9441248) q[2];
sx q[2];
rz(0.8651274) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1358007) q[1];
sx q[1];
rz(-1.5151007) q[1];
sx q[1];
rz(-2.9018351) q[1];
rz(-pi) q[2];
rz(-0.17574163) q[3];
sx q[3];
rz(-2.1648241) q[3];
sx q[3];
rz(0.16080805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1069964) q[2];
sx q[2];
rz(-0.84356374) q[2];
sx q[2];
rz(-1.754388) q[2];
rz(-2.710279) q[3];
sx q[3];
rz(-1.8547736) q[3];
sx q[3];
rz(1.0881902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4819734) q[0];
sx q[0];
rz(-2.1932333) q[0];
sx q[0];
rz(-1.5455998) q[0];
rz(-1.1384456) q[1];
sx q[1];
rz(-2.3661416) q[1];
sx q[1];
rz(2.0514354) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6369612) q[0];
sx q[0];
rz(-1.9177027) q[0];
sx q[0];
rz(-2.7142801) q[0];
rz(1.0341187) q[2];
sx q[2];
rz(-0.82674971) q[2];
sx q[2];
rz(-0.85988753) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7294457) q[1];
sx q[1];
rz(-0.73819654) q[1];
sx q[1];
rz(-0.63936887) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7100536) q[3];
sx q[3];
rz(-2.6772237) q[3];
sx q[3];
rz(-2.4114256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1426992) q[2];
sx q[2];
rz(-2.6901851) q[2];
sx q[2];
rz(2.2857655) q[2];
rz(-1.9479729) q[3];
sx q[3];
rz(-1.6222298) q[3];
sx q[3];
rz(-2.2284171) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6511433) q[0];
sx q[0];
rz(-0.97390807) q[0];
sx q[0];
rz(-2.9918616) q[0];
rz(0.99114746) q[1];
sx q[1];
rz(-1.2065572) q[1];
sx q[1];
rz(-1.1700464) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2802551) q[0];
sx q[0];
rz(-2.4147408) q[0];
sx q[0];
rz(-1.0759541) q[0];
rz(-0.32355002) q[2];
sx q[2];
rz(-2.2252482) q[2];
sx q[2];
rz(-2.8138585) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7331446) q[1];
sx q[1];
rz(-2.0093579) q[1];
sx q[1];
rz(1.1327729) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.53797651) q[3];
sx q[3];
rz(-1.4036944) q[3];
sx q[3];
rz(2.3218384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2698007) q[2];
sx q[2];
rz(-1.5503927) q[2];
sx q[2];
rz(3.0043547) q[2];
rz(1.3373226) q[3];
sx q[3];
rz(-2.5604355) q[3];
sx q[3];
rz(-1.2924682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9542434) q[0];
sx q[0];
rz(-1.7631148) q[0];
sx q[0];
rz(-1.3504008) q[0];
rz(0.84287914) q[1];
sx q[1];
rz(-0.73892361) q[1];
sx q[1];
rz(-0.67289105) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85149375) q[0];
sx q[0];
rz(-2.5585472) q[0];
sx q[0];
rz(1.3821938) q[0];
rz(-1.599309) q[2];
sx q[2];
rz(-0.71454222) q[2];
sx q[2];
rz(1.6376094) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3402965) q[1];
sx q[1];
rz(-0.9749037) q[1];
sx q[1];
rz(-1.8227541) q[1];
rz(-pi) q[2];
rz(-1.6478496) q[3];
sx q[3];
rz(-1.1601845) q[3];
sx q[3];
rz(-1.3239469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.3488397) q[2];
sx q[2];
rz(-1.2092084) q[2];
sx q[2];
rz(-0.43506452) q[2];
rz(1.3600291) q[3];
sx q[3];
rz(-2.3924148) q[3];
sx q[3];
rz(-2.8939261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.172794) q[0];
sx q[0];
rz(-2.2491169) q[0];
sx q[0];
rz(-2.0794179) q[0];
rz(1.1116213) q[1];
sx q[1];
rz(-1.2373135) q[1];
sx q[1];
rz(-1.7395082) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8334478) q[0];
sx q[0];
rz(-1.2045367) q[0];
sx q[0];
rz(1.2783947) q[0];
rz(-1.3349322) q[2];
sx q[2];
rz(-2.3216322) q[2];
sx q[2];
rz(1.6279398) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.40560383) q[1];
sx q[1];
rz(-2.7466008) q[1];
sx q[1];
rz(-2.7337381) q[1];
x q[2];
rz(-3.0038463) q[3];
sx q[3];
rz(-2.1546954) q[3];
sx q[3];
rz(-3.0018842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7193675) q[2];
sx q[2];
rz(-2.5740467) q[2];
sx q[2];
rz(-2.4613703) q[2];
rz(-2.7133572) q[3];
sx q[3];
rz(-1.2546344) q[3];
sx q[3];
rz(-1.7817106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5087886) q[0];
sx q[0];
rz(-1.4161685) q[0];
sx q[0];
rz(1.5493786) q[0];
rz(2.8920065) q[1];
sx q[1];
rz(-1.9892178) q[1];
sx q[1];
rz(2.6002398) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66051018) q[0];
sx q[0];
rz(-0.56429243) q[0];
sx q[0];
rz(-0.020945992) q[0];
rz(0.45221046) q[2];
sx q[2];
rz(-1.2621677) q[2];
sx q[2];
rz(1.897915) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9882422) q[1];
sx q[1];
rz(-1.1594979) q[1];
sx q[1];
rz(-0.78505713) q[1];
x q[2];
rz(0.3803216) q[3];
sx q[3];
rz(-1.5725279) q[3];
sx q[3];
rz(0.94554949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0970739) q[2];
sx q[2];
rz(-2.7842583) q[2];
sx q[2];
rz(2.3642335) q[2];
rz(0.85123953) q[3];
sx q[3];
rz(-1.0612396) q[3];
sx q[3];
rz(-1.3210993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7580496) q[0];
sx q[0];
rz(-1.3306916) q[0];
sx q[0];
rz(0.0099649075) q[0];
rz(-2.126157) q[1];
sx q[1];
rz(-0.76428691) q[1];
sx q[1];
rz(1.4452971) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89307071) q[0];
sx q[0];
rz(-1.0501486) q[0];
sx q[0];
rz(-2.394886) q[0];
rz(-pi) q[1];
rz(-0.89500918) q[2];
sx q[2];
rz(-2.0446606) q[2];
sx q[2];
rz(-0.5627788) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.086844079) q[1];
sx q[1];
rz(-2.5111685) q[1];
sx q[1];
rz(0.026442095) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2441838) q[3];
sx q[3];
rz(-1.0703147) q[3];
sx q[3];
rz(0.15411479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.69892591) q[2];
sx q[2];
rz(-0.93124229) q[2];
sx q[2];
rz(-0.60738579) q[2];
rz(-1.4390885) q[3];
sx q[3];
rz(-1.3834229) q[3];
sx q[3];
rz(2.3506892) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33567515) q[0];
sx q[0];
rz(-1.6573925) q[0];
sx q[0];
rz(2.5277396) q[0];
rz(-2.0954258) q[1];
sx q[1];
rz(-2.8764953) q[1];
sx q[1];
rz(0.39224958) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2989202) q[0];
sx q[0];
rz(-0.93085104) q[0];
sx q[0];
rz(-1.0660893) q[0];
rz(-pi) q[1];
rz(2.1456111) q[2];
sx q[2];
rz(-1.0258342) q[2];
sx q[2];
rz(2.0202707) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.039672893) q[1];
sx q[1];
rz(-0.79636785) q[1];
sx q[1];
rz(2.2069195) q[1];
rz(0.20499968) q[3];
sx q[3];
rz(-1.5424171) q[3];
sx q[3];
rz(-1.8207707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0739416) q[2];
sx q[2];
rz(-1.757688) q[2];
sx q[2];
rz(-2.5411141) q[2];
rz(-1.0673808) q[3];
sx q[3];
rz(-1.821527) q[3];
sx q[3];
rz(-2.0480806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4338715) q[0];
sx q[0];
rz(-2.7086471) q[0];
sx q[0];
rz(1.4824296) q[0];
rz(-2.1451163) q[1];
sx q[1];
rz(-1.4947718) q[1];
sx q[1];
rz(1.6047118) q[1];
rz(0.76657961) q[2];
sx q[2];
rz(-1.3091514) q[2];
sx q[2];
rz(-1.8717629) q[2];
rz(2.5946887) q[3];
sx q[3];
rz(-2.1756267) q[3];
sx q[3];
rz(-0.970943) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
