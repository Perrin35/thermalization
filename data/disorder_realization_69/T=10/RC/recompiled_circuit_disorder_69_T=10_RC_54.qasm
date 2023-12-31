OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.34713137) q[0];
sx q[0];
rz(-1.0153271) q[0];
sx q[0];
rz(0.46749687) q[0];
rz(-0.52019083) q[1];
sx q[1];
rz(-1.3462892) q[1];
sx q[1];
rz(2.4612114) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2517393) q[0];
sx q[0];
rz(-1.7008874) q[0];
sx q[0];
rz(2.2962909) q[0];
x q[1];
rz(0.19619588) q[2];
sx q[2];
rz(-1.8610524) q[2];
sx q[2];
rz(0.44858518) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8531224) q[1];
sx q[1];
rz(-2.5783357) q[1];
sx q[1];
rz(1.1239777) q[1];
x q[2];
rz(-0.33403553) q[3];
sx q[3];
rz(-1.4035657) q[3];
sx q[3];
rz(0.86080307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.15158571) q[2];
sx q[2];
rz(-2.1534584) q[2];
sx q[2];
rz(3.0541259) q[2];
rz(-0.72922373) q[3];
sx q[3];
rz(-2.7681523) q[3];
sx q[3];
rz(0.27174404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7049578) q[0];
sx q[0];
rz(-0.74902642) q[0];
sx q[0];
rz(-2.1402284) q[0];
rz(-0.17240605) q[1];
sx q[1];
rz(-1.1162076) q[1];
sx q[1];
rz(0.52406812) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8952626) q[0];
sx q[0];
rz(-1.2753914) q[0];
sx q[0];
rz(-2.2569879) q[0];
x q[1];
rz(-1.0753724) q[2];
sx q[2];
rz(-0.76250695) q[2];
sx q[2];
rz(-2.5410595) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9429051) q[1];
sx q[1];
rz(-2.6932979) q[1];
sx q[1];
rz(-2.5793736) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53701138) q[3];
sx q[3];
rz(-1.4376663) q[3];
sx q[3];
rz(2.6739124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.15741631) q[2];
sx q[2];
rz(-0.45980644) q[2];
sx q[2];
rz(1.8015507) q[2];
rz(2.346758) q[3];
sx q[3];
rz(-1.1398311) q[3];
sx q[3];
rz(-3.1047344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79725093) q[0];
sx q[0];
rz(-2.4448555) q[0];
sx q[0];
rz(0.62477338) q[0];
rz(-0.96427381) q[1];
sx q[1];
rz(-0.48502973) q[1];
sx q[1];
rz(-2.952081) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7127842) q[0];
sx q[0];
rz(-1.3002987) q[0];
sx q[0];
rz(-2.4012647) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7663642) q[2];
sx q[2];
rz(-1.669075) q[2];
sx q[2];
rz(-0.74429846) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5513969) q[1];
sx q[1];
rz(-1.8101748) q[1];
sx q[1];
rz(1.6281284) q[1];
rz(-pi) q[2];
x q[2];
rz(2.965851) q[3];
sx q[3];
rz(-2.1648241) q[3];
sx q[3];
rz(-2.9807846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1069964) q[2];
sx q[2];
rz(-0.84356374) q[2];
sx q[2];
rz(-1.3872046) q[2];
rz(-2.710279) q[3];
sx q[3];
rz(-1.8547736) q[3];
sx q[3];
rz(1.0881902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4819734) q[0];
sx q[0];
rz(-0.94835931) q[0];
sx q[0];
rz(-1.5455998) q[0];
rz(-1.1384456) q[1];
sx q[1];
rz(-2.3661416) q[1];
sx q[1];
rz(2.0514354) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6369612) q[0];
sx q[0];
rz(-1.22389) q[0];
sx q[0];
rz(2.7142801) q[0];
rz(-pi) q[1];
rz(1.0341187) q[2];
sx q[2];
rz(-0.82674971) q[2];
sx q[2];
rz(-0.85988753) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.62413299) q[1];
sx q[1];
rz(-1.0003261) q[1];
sx q[1];
rz(-1.073451) q[1];
rz(-pi) q[2];
rz(-0.4269883) q[3];
sx q[3];
rz(-1.7592351) q[3];
sx q[3];
rz(-2.6915336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.99889341) q[2];
sx q[2];
rz(-0.45140758) q[2];
sx q[2];
rz(0.85582716) q[2];
rz(-1.1936197) q[3];
sx q[3];
rz(-1.6222298) q[3];
sx q[3];
rz(2.2284171) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6511433) q[0];
sx q[0];
rz(-0.97390807) q[0];
sx q[0];
rz(2.9918616) q[0];
rz(-0.99114746) q[1];
sx q[1];
rz(-1.2065572) q[1];
sx q[1];
rz(-1.9715462) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4677508) q[0];
sx q[0];
rz(-1.8918599) q[0];
sx q[0];
rz(2.2348316) q[0];
x q[1];
rz(0.32355002) q[2];
sx q[2];
rz(-0.91634446) q[2];
sx q[2];
rz(0.3277342) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7331446) q[1];
sx q[1];
rz(-1.1322347) q[1];
sx q[1];
rz(2.0088197) q[1];
rz(-pi) q[2];
rz(-1.3768457) q[3];
sx q[3];
rz(-2.100482) q[3];
sx q[3];
rz(-0.84996163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2698007) q[2];
sx q[2];
rz(-1.5912) q[2];
sx q[2];
rz(0.13723792) q[2];
rz(-1.3373226) q[3];
sx q[3];
rz(-0.58115712) q[3];
sx q[3];
rz(-1.2924682) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9542434) q[0];
sx q[0];
rz(-1.3784778) q[0];
sx q[0];
rz(-1.7911918) q[0];
rz(-0.84287914) q[1];
sx q[1];
rz(-2.402669) q[1];
sx q[1];
rz(-0.67289105) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2900989) q[0];
sx q[0];
rz(-0.58304542) q[0];
sx q[0];
rz(-1.7593988) q[0];
rz(-pi) q[1];
x q[1];
rz(0.024725155) q[2];
sx q[2];
rz(-2.284986) q[2];
sx q[2];
rz(-1.6753472) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8012961) q[1];
sx q[1];
rz(-0.9749037) q[1];
sx q[1];
rz(-1.8227541) q[1];
rz(-2.7298922) q[3];
sx q[3];
rz(-1.6414335) q[3];
sx q[3];
rz(-0.27765805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.3488397) q[2];
sx q[2];
rz(-1.9323843) q[2];
sx q[2];
rz(0.43506452) q[2];
rz(-1.3600291) q[3];
sx q[3];
rz(-2.3924148) q[3];
sx q[3];
rz(-0.24766651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.172794) q[0];
sx q[0];
rz(-2.2491169) q[0];
sx q[0];
rz(-1.0621747) q[0];
rz(1.1116213) q[1];
sx q[1];
rz(-1.2373135) q[1];
sx q[1];
rz(-1.7395082) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60915011) q[0];
sx q[0];
rz(-2.6770868) q[0];
sx q[0];
rz(-2.4971278) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8066605) q[2];
sx q[2];
rz(-0.81996041) q[2];
sx q[2];
rz(1.6279398) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2981616) q[1];
sx q[1];
rz(-1.9318252) q[1];
sx q[1];
rz(1.7346738) q[1];
rz(-pi) q[2];
rz(0.98251179) q[3];
sx q[3];
rz(-1.4559828) q[3];
sx q[3];
rz(1.3548152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7193675) q[2];
sx q[2];
rz(-0.56754595) q[2];
sx q[2];
rz(2.4613703) q[2];
rz(2.7133572) q[3];
sx q[3];
rz(-1.2546344) q[3];
sx q[3];
rz(1.7817106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(0.632804) q[0];
sx q[0];
rz(-1.7254242) q[0];
sx q[0];
rz(1.592214) q[0];
rz(2.8920065) q[1];
sx q[1];
rz(-1.1523749) q[1];
sx q[1];
rz(-2.6002398) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66051018) q[0];
sx q[0];
rz(-2.5773002) q[0];
sx q[0];
rz(-3.1206467) q[0];
x q[1];
rz(0.63033732) q[2];
sx q[2];
rz(-0.54140831) q[2];
sx q[2];
rz(0.88592096) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9882422) q[1];
sx q[1];
rz(-1.9820947) q[1];
sx q[1];
rz(0.78505713) q[1];
rz(0.0046644966) q[3];
sx q[3];
rz(-0.38032535) q[3];
sx q[3];
rz(0.62957803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.044518746) q[2];
sx q[2];
rz(-0.35733435) q[2];
sx q[2];
rz(0.77735916) q[2];
rz(-0.85123953) q[3];
sx q[3];
rz(-2.080353) q[3];
sx q[3];
rz(-1.3210993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7580496) q[0];
sx q[0];
rz(-1.3306916) q[0];
sx q[0];
rz(3.1316277) q[0];
rz(-2.126157) q[1];
sx q[1];
rz(-2.3773057) q[1];
sx q[1];
rz(1.6962956) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2485219) q[0];
sx q[0];
rz(-2.0914441) q[0];
sx q[0];
rz(-0.74670665) q[0];
rz(-pi) q[1];
rz(2.2465835) q[2];
sx q[2];
rz(-1.0969321) q[2];
sx q[2];
rz(0.5627788) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6362793) q[1];
sx q[1];
rz(-1.5863824) q[1];
sx q[1];
rz(0.63025766) q[1];
rz(-pi) q[2];
x q[2];
rz(0.85083665) q[3];
sx q[3];
rz(-0.81504226) q[3];
sx q[3];
rz(-1.9581865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4426667) q[2];
sx q[2];
rz(-2.2103504) q[2];
sx q[2];
rz(2.5342069) q[2];
rz(1.7025042) q[3];
sx q[3];
rz(-1.3834229) q[3];
sx q[3];
rz(2.3506892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33567515) q[0];
sx q[0];
rz(-1.6573925) q[0];
sx q[0];
rz(-2.5277396) q[0];
rz(2.0954258) q[1];
sx q[1];
rz(-0.26509735) q[1];
sx q[1];
rz(-2.7493431) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0948254) q[0];
sx q[0];
rz(-1.1724768) q[0];
sx q[0];
rz(0.70478583) q[0];
rz(-pi) q[1];
rz(2.1456111) q[2];
sx q[2];
rz(-2.1157584) q[2];
sx q[2];
rz(-2.0202707) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2892462) q[1];
sx q[1];
rz(-0.95818555) q[1];
sx q[1];
rz(0.5457408) q[1];
x q[2];
rz(-1.5418105) q[3];
sx q[3];
rz(-1.7757123) q[3];
sx q[3];
rz(2.8857185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.067651) q[2];
sx q[2];
rz(-1.757688) q[2];
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
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.4338715) q[0];
sx q[0];
rz(-2.7086471) q[0];
sx q[0];
rz(1.4824296) q[0];
rz(0.99647635) q[1];
sx q[1];
rz(-1.4947718) q[1];
sx q[1];
rz(1.6047118) q[1];
rz(-0.36841064) q[2];
sx q[2];
rz(-2.340292) q[2];
sx q[2];
rz(3.1030263) q[2];
rz(-0.546904) q[3];
sx q[3];
rz(-2.1756267) q[3];
sx q[3];
rz(-0.970943) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
