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
rz(1.7953035) q[1];
sx q[1];
rz(10.105159) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43359783) q[0];
sx q[0];
rz(-0.85277075) q[0];
sx q[0];
rz(-0.17311592) q[0];
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
sx q[0];
rz(-pi/2) q[0];
rz(2.2441794) q[1];
sx q[1];
rz(-1.3379828) q[1];
sx q[1];
rz(-2.0884872) q[1];
rz(-pi) q[2];
rz(-0.33403553) q[3];
sx q[3];
rz(-1.4035657) q[3];
sx q[3];
rz(0.86080307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.15158571) q[2];
sx q[2];
rz(-0.9881343) q[2];
sx q[2];
rz(3.0541259) q[2];
rz(2.4123689) q[3];
sx q[3];
rz(-0.37344033) q[3];
sx q[3];
rz(2.8698486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4366348) q[0];
sx q[0];
rz(-0.74902642) q[0];
sx q[0];
rz(2.1402284) q[0];
rz(-0.17240605) q[1];
sx q[1];
rz(-1.1162076) q[1];
sx q[1];
rz(0.52406812) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8952626) q[0];
sx q[0];
rz(-1.2753914) q[0];
sx q[0];
rz(2.2569879) q[0];
rz(0.42627724) q[2];
sx q[2];
rz(-2.2239416) q[2];
sx q[2];
rz(3.1003568) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9429051) q[1];
sx q[1];
rz(-2.6932979) q[1];
sx q[1];
rz(2.5793736) q[1];
rz(-1.7254132) q[3];
sx q[3];
rz(-1.039045) q[3];
sx q[3];
rz(-2.1173409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.15741631) q[2];
sx q[2];
rz(-2.6817862) q[2];
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
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3443417) q[0];
sx q[0];
rz(-2.4448555) q[0];
sx q[0];
rz(0.62477338) q[0];
rz(0.96427381) q[1];
sx q[1];
rz(-0.48502973) q[1];
sx q[1];
rz(2.952081) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42650578) q[0];
sx q[0];
rz(-0.77930342) q[0];
sx q[0];
rz(2.7515609) q[0];
rz(2.7663642) q[2];
sx q[2];
rz(-1.669075) q[2];
sx q[2];
rz(-2.3972942) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1358007) q[1];
sx q[1];
rz(-1.5151007) q[1];
sx q[1];
rz(2.9018351) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8241006) q[3];
sx q[3];
rz(-2.5251303) q[3];
sx q[3];
rz(-2.9951819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0345962) q[2];
sx q[2];
rz(-2.2980289) q[2];
sx q[2];
rz(1.3872046) q[2];
rz(2.710279) q[3];
sx q[3];
rz(-1.286819) q[3];
sx q[3];
rz(-2.0534024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(2.4819734) q[0];
sx q[0];
rz(-0.94835931) q[0];
sx q[0];
rz(1.5455998) q[0];
rz(2.003147) q[1];
sx q[1];
rz(-0.7754511) q[1];
sx q[1];
rz(1.0901573) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2290303) q[0];
sx q[0];
rz(-1.1704485) q[0];
sx q[0];
rz(1.9489261) q[0];
rz(-pi) q[1];
rz(2.3218669) q[2];
sx q[2];
rz(-1.956454) q[2];
sx q[2];
rz(1.0939329) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7294457) q[1];
sx q[1];
rz(-2.4033961) q[1];
sx q[1];
rz(2.5022238) q[1];
x q[2];
rz(2.7146043) q[3];
sx q[3];
rz(-1.3823576) q[3];
sx q[3];
rz(2.6915336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.99889341) q[2];
sx q[2];
rz(-0.45140758) q[2];
sx q[2];
rz(-2.2857655) q[2];
rz(-1.9479729) q[3];
sx q[3];
rz(-1.6222298) q[3];
sx q[3];
rz(-2.2284171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6511433) q[0];
sx q[0];
rz(-0.97390807) q[0];
sx q[0];
rz(-0.14973101) q[0];
rz(-2.1504452) q[1];
sx q[1];
rz(-1.9350355) q[1];
sx q[1];
rz(1.1700464) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2802551) q[0];
sx q[0];
rz(-2.4147408) q[0];
sx q[0];
rz(2.0656385) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8180426) q[2];
sx q[2];
rz(-0.91634446) q[2];
sx q[2];
rz(2.8138585) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7331446) q[1];
sx q[1];
rz(-2.0093579) q[1];
sx q[1];
rz(-1.1327729) q[1];
rz(1.7647469) q[3];
sx q[3];
rz(-2.100482) q[3];
sx q[3];
rz(2.291631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.871792) q[2];
sx q[2];
rz(-1.5503927) q[2];
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
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18734922) q[0];
sx q[0];
rz(-1.7631148) q[0];
sx q[0];
rz(-1.3504008) q[0];
rz(-2.2987135) q[1];
sx q[1];
rz(-0.73892361) q[1];
sx q[1];
rz(-0.67289105) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.580299) q[0];
sx q[0];
rz(-1.6742047) q[0];
sx q[0];
rz(-2.1456477) q[0];
rz(0.85645533) q[2];
sx q[2];
rz(-1.5894784) q[2];
sx q[2];
rz(3.0532388) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.91025483) q[1];
sx q[1];
rz(-0.64097039) q[1];
sx q[1];
rz(2.7892968) q[1];
rz(-2.9665885) q[3];
sx q[3];
rz(-2.724218) q[3];
sx q[3];
rz(-2.008703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.172794) q[0];
sx q[0];
rz(-0.89247576) q[0];
sx q[0];
rz(2.0794179) q[0];
rz(2.0299714) q[1];
sx q[1];
rz(-1.2373135) q[1];
sx q[1];
rz(-1.4020845) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5324425) q[0];
sx q[0];
rz(-2.6770868) q[0];
sx q[0];
rz(2.4971278) q[0];
rz(1.3349322) q[2];
sx q[2];
rz(-2.3216322) q[2];
sx q[2];
rz(1.5136528) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.78571037) q[1];
sx q[1];
rz(-1.4175698) q[1];
sx q[1];
rz(2.7760837) q[1];
rz(-3.0038463) q[3];
sx q[3];
rz(-2.1546954) q[3];
sx q[3];
rz(-3.0018842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7193675) q[2];
sx q[2];
rz(-2.5740467) q[2];
sx q[2];
rz(-2.4613703) q[2];
rz(0.42823544) q[3];
sx q[3];
rz(-1.8869583) q[3];
sx q[3];
rz(1.7817106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.632804) q[0];
sx q[0];
rz(-1.7254242) q[0];
sx q[0];
rz(-1.592214) q[0];
rz(2.8920065) q[1];
sx q[1];
rz(-1.9892178) q[1];
sx q[1];
rz(-0.54135281) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66051018) q[0];
sx q[0];
rz(-0.56429243) q[0];
sx q[0];
rz(3.1206467) q[0];
rz(-pi) q[1];
rz(-1.2301684) q[2];
sx q[2];
rz(-2.0001786) q[2];
sx q[2];
rz(-2.9609749) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7975446) q[1];
sx q[1];
rz(-2.27564) q[1];
sx q[1];
rz(1.0182347) q[1];
x q[2];
rz(-0.3803216) q[3];
sx q[3];
rz(-1.5725279) q[3];
sx q[3];
rz(-0.94554949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.044518746) q[2];
sx q[2];
rz(-2.7842583) q[2];
sx q[2];
rz(0.77735916) q[2];
rz(2.2903531) q[3];
sx q[3];
rz(-2.080353) q[3];
sx q[3];
rz(-1.3210993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7580496) q[0];
sx q[0];
rz(-1.3306916) q[0];
sx q[0];
rz(0.0099649075) q[0];
rz(2.126157) q[1];
sx q[1];
rz(-0.76428691) q[1];
sx q[1];
rz(1.6962956) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9706668) q[0];
sx q[0];
rz(-0.88060856) q[0];
sx q[0];
rz(2.4404581) q[0];
rz(-pi) q[1];
rz(-2.2465835) q[2];
sx q[2];
rz(-1.0969321) q[2];
sx q[2];
rz(2.5788139) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.086844079) q[1];
sx q[1];
rz(-0.63042414) q[1];
sx q[1];
rz(0.026442095) q[1];
rz(-pi) q[2];
rz(-0.85083665) q[3];
sx q[3];
rz(-2.3265504) q[3];
sx q[3];
rz(1.1834061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.69892591) q[2];
sx q[2];
rz(-2.2103504) q[2];
sx q[2];
rz(-0.60738579) q[2];
rz(1.7025042) q[3];
sx q[3];
rz(-1.7581698) q[3];
sx q[3];
rz(-2.3506892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8059175) q[0];
sx q[0];
rz(-1.4842002) q[0];
sx q[0];
rz(2.5277396) q[0];
rz(1.0461668) q[1];
sx q[1];
rz(-0.26509735) q[1];
sx q[1];
rz(2.7493431) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0948254) q[0];
sx q[0];
rz(-1.1724768) q[0];
sx q[0];
rz(-2.4368068) q[0];
rz(-pi) q[1];
x q[1];
rz(0.62551542) q[2];
sx q[2];
rz(-1.0872456) q[2];
sx q[2];
rz(3.0160883) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0872005) q[1];
sx q[1];
rz(-2.0093845) q[1];
sx q[1];
rz(2.2589576) q[1];
x q[2];
rz(-0.13855374) q[3];
sx q[3];
rz(-2.9346653) q[3];
sx q[3];
rz(3.0272527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.067651) q[2];
sx q[2];
rz(-1.3839046) q[2];
sx q[2];
rz(-0.60047853) q[2];
rz(-1.0673808) q[3];
sx q[3];
rz(-1.821527) q[3];
sx q[3];
rz(-2.0480806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7077211) q[0];
sx q[0];
rz(-2.7086471) q[0];
sx q[0];
rz(1.4824296) q[0];
rz(2.1451163) q[1];
sx q[1];
rz(-1.6468208) q[1];
sx q[1];
rz(-1.5368808) q[1];
rz(2.773182) q[2];
sx q[2];
rz(-2.340292) q[2];
sx q[2];
rz(3.1030263) q[2];
rz(2.2157833) q[3];
sx q[3];
rz(-2.3498597) q[3];
sx q[3];
rz(1.3510977) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
