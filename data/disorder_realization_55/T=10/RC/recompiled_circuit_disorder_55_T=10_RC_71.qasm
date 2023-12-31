OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.54685932) q[0];
sx q[0];
rz(-1.62513) q[0];
sx q[0];
rz(-0.2642785) q[0];
rz(-0.9737941) q[1];
sx q[1];
rz(5.073054) q[1];
sx q[1];
rz(10.160025) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41266325) q[0];
sx q[0];
rz(-0.67617765) q[0];
sx q[0];
rz(2.9039608) q[0];
rz(-pi) q[1];
rz(-2.5392883) q[2];
sx q[2];
rz(-2.3659083) q[2];
sx q[2];
rz(1.3210981) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3056065) q[1];
sx q[1];
rz(-1.4164093) q[1];
sx q[1];
rz(-1.2812213) q[1];
x q[2];
rz(0.04282184) q[3];
sx q[3];
rz(-2.560727) q[3];
sx q[3];
rz(-2.7659741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.66951093) q[2];
sx q[2];
rz(-1.8410204) q[2];
sx q[2];
rz(-1.1038587) q[2];
rz(1.8707229) q[3];
sx q[3];
rz(-1.2277675) q[3];
sx q[3];
rz(0.27403533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1141777) q[0];
sx q[0];
rz(-0.52863055) q[0];
sx q[0];
rz(-0.43637481) q[0];
rz(-2.6787058) q[1];
sx q[1];
rz(-1.0375689) q[1];
sx q[1];
rz(-0.26611051) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9929745) q[0];
sx q[0];
rz(-2.2668112) q[0];
sx q[0];
rz(0.50177411) q[0];
rz(0.99533178) q[2];
sx q[2];
rz(-1.7082214) q[2];
sx q[2];
rz(-0.95057887) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4437342) q[1];
sx q[1];
rz(-2.2599972) q[1];
sx q[1];
rz(-0.84390784) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8310043) q[3];
sx q[3];
rz(-2.5187413) q[3];
sx q[3];
rz(0.42850307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.46488547) q[2];
sx q[2];
rz(-1.2767982) q[2];
sx q[2];
rz(-2.6300988) q[2];
rz(2.3320847) q[3];
sx q[3];
rz(-1.6102689) q[3];
sx q[3];
rz(2.8539343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
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
rz(2.4354316) q[0];
sx q[0];
rz(-1.5478739) q[0];
sx q[0];
rz(-0.92873746) q[0];
rz(1.4061032) q[1];
sx q[1];
rz(-2.4415253) q[1];
sx q[1];
rz(1.3471289) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72235332) q[0];
sx q[0];
rz(-1.1142715) q[0];
sx q[0];
rz(1.2786352) q[0];
rz(1.0727097) q[2];
sx q[2];
rz(-2.8194397) q[2];
sx q[2];
rz(-2.1307532) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9007064) q[1];
sx q[1];
rz(-0.97929231) q[1];
sx q[1];
rz(0.35066168) q[1];
rz(0.8637572) q[3];
sx q[3];
rz(-1.6861526) q[3];
sx q[3];
rz(-0.74187169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1469664) q[2];
sx q[2];
rz(-1.3866084) q[2];
sx q[2];
rz(1.6195126) q[2];
rz(0.26432031) q[3];
sx q[3];
rz(-1.0364573) q[3];
sx q[3];
rz(0.49595293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79214823) q[0];
sx q[0];
rz(-1.148372) q[0];
sx q[0];
rz(2.175892) q[0];
rz(-0.72215885) q[1];
sx q[1];
rz(-1.637371) q[1];
sx q[1];
rz(2.5818363) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1346261) q[0];
sx q[0];
rz(-0.43381938) q[0];
sx q[0];
rz(-2.9475648) q[0];
rz(-pi) q[1];
rz(2.431805) q[2];
sx q[2];
rz(-0.94883942) q[2];
sx q[2];
rz(-0.70993916) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.30448118) q[1];
sx q[1];
rz(-1.6164391) q[1];
sx q[1];
rz(-0.95506217) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3523931) q[3];
sx q[3];
rz(-1.7390828) q[3];
sx q[3];
rz(-2.5326953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.42797783) q[2];
sx q[2];
rz(-1.6327991) q[2];
sx q[2];
rz(-1.4245865) q[2];
rz(-0.26040855) q[3];
sx q[3];
rz(-1.3924761) q[3];
sx q[3];
rz(2.7105455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99073064) q[0];
sx q[0];
rz(-1.6435511) q[0];
sx q[0];
rz(0.36636233) q[0];
rz(-1.5461961) q[1];
sx q[1];
rz(-0.55391824) q[1];
sx q[1];
rz(-2.7979134) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3643091) q[0];
sx q[0];
rz(-1.0316348) q[0];
sx q[0];
rz(2.3794412) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5630066) q[2];
sx q[2];
rz(-2.1861976) q[2];
sx q[2];
rz(1.7631284) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1829454) q[1];
sx q[1];
rz(-1.4970386) q[1];
sx q[1];
rz(1.0720836) q[1];
rz(-0.021275612) q[3];
sx q[3];
rz(-2.2145503) q[3];
sx q[3];
rz(-1.4847886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7498103) q[2];
sx q[2];
rz(-0.85931531) q[2];
sx q[2];
rz(3.0878477) q[2];
rz(1.404445) q[3];
sx q[3];
rz(-0.497538) q[3];
sx q[3];
rz(0.29156175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4869726) q[0];
sx q[0];
rz(-0.64990652) q[0];
sx q[0];
rz(2.3858331) q[0];
rz(-3.1164363) q[1];
sx q[1];
rz(-2.2143366) q[1];
sx q[1];
rz(2.8818534) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.040859) q[0];
sx q[0];
rz(-1.20964) q[0];
sx q[0];
rz(0.33613899) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9070712) q[2];
sx q[2];
rz(-1.5894801) q[2];
sx q[2];
rz(-2.3980354) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7139587) q[1];
sx q[1];
rz(-2.7379588) q[1];
sx q[1];
rz(2.1778818) q[1];
rz(-pi) q[2];
rz(1.5335347) q[3];
sx q[3];
rz(-0.44964368) q[3];
sx q[3];
rz(2.8763308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6529237) q[2];
sx q[2];
rz(-2.3355464) q[2];
sx q[2];
rz(0.10061131) q[2];
rz(-0.18209022) q[3];
sx q[3];
rz(-2.263335) q[3];
sx q[3];
rz(-1.3476868) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9942193) q[0];
sx q[0];
rz(-1.4202776) q[0];
sx q[0];
rz(-0.31016645) q[0];
rz(-0.50225964) q[1];
sx q[1];
rz(-2.6175833) q[1];
sx q[1];
rz(-0.60595864) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92354846) q[0];
sx q[0];
rz(-0.96203066) q[0];
sx q[0];
rz(-1.2952842) q[0];
x q[1];
rz(0.29166834) q[2];
sx q[2];
rz(-1.3783611) q[2];
sx q[2];
rz(-2.8067436) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0026605) q[1];
sx q[1];
rz(-1.2298889) q[1];
sx q[1];
rz(-1.0374116) q[1];
rz(-2.3378387) q[3];
sx q[3];
rz(-0.70138068) q[3];
sx q[3];
rz(2.3826117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.24017748) q[2];
sx q[2];
rz(-1.1837974) q[2];
sx q[2];
rz(2.288738) q[2];
rz(-1.7715706) q[3];
sx q[3];
rz(-1.6829237) q[3];
sx q[3];
rz(0.20496932) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9119499) q[0];
sx q[0];
rz(-2.5456173) q[0];
sx q[0];
rz(1.6802616) q[0];
rz(-1.4029067) q[1];
sx q[1];
rz(-0.97424126) q[1];
sx q[1];
rz(0.064037474) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2647576) q[0];
sx q[0];
rz(-1.7581853) q[0];
sx q[0];
rz(-0.95215709) q[0];
x q[1];
rz(-2.8352751) q[2];
sx q[2];
rz(-2.0567354) q[2];
sx q[2];
rz(1.475032) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.21676) q[1];
sx q[1];
rz(-1.2158582) q[1];
sx q[1];
rz(-1.6096398) q[1];
x q[2];
rz(-1.2588345) q[3];
sx q[3];
rz(-1.6083816) q[3];
sx q[3];
rz(3.104044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.091207592) q[2];
sx q[2];
rz(-2.51077) q[2];
sx q[2];
rz(1.5861661) q[2];
rz(-2.2533916) q[3];
sx q[3];
rz(-1.9017838) q[3];
sx q[3];
rz(2.2122038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14426194) q[0];
sx q[0];
rz(-2.0781131) q[0];
sx q[0];
rz(-1.4319179) q[0];
rz(-0.56888467) q[1];
sx q[1];
rz(-0.535393) q[1];
sx q[1];
rz(-2.0137537) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.226798) q[0];
sx q[0];
rz(-0.1495805) q[0];
sx q[0];
rz(3.0339255) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.88840975) q[2];
sx q[2];
rz(-0.75196224) q[2];
sx q[2];
rz(-1.0351406) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4081501) q[1];
sx q[1];
rz(-0.24194716) q[1];
sx q[1];
rz(-1.6192295) q[1];
x q[2];
rz(1.0206251) q[3];
sx q[3];
rz(-2.5084825) q[3];
sx q[3];
rz(1.2596631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.52788064) q[2];
sx q[2];
rz(-0.85313672) q[2];
sx q[2];
rz(2.9679427) q[2];
rz(2.8052143) q[3];
sx q[3];
rz(-1.222638) q[3];
sx q[3];
rz(-2.7500847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7062475) q[0];
sx q[0];
rz(-0.58723891) q[0];
sx q[0];
rz(1.4655112) q[0];
rz(-2.3174875) q[1];
sx q[1];
rz(-1.6128287) q[1];
sx q[1];
rz(-2.5691659) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3819645) q[0];
sx q[0];
rz(-2.1158764) q[0];
sx q[0];
rz(0.99960021) q[0];
rz(-pi) q[1];
rz(-2.0881537) q[2];
sx q[2];
rz(-0.51689076) q[2];
sx q[2];
rz(-3.0992532) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.010667) q[1];
sx q[1];
rz(-1.3548684) q[1];
sx q[1];
rz(-0.24433498) q[1];
rz(-pi) q[2];
rz(0.33705538) q[3];
sx q[3];
rz(-0.73738499) q[3];
sx q[3];
rz(0.22102236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.77999014) q[2];
sx q[2];
rz(-2.6670691) q[2];
sx q[2];
rz(2.6043716) q[2];
rz(2.0843263) q[3];
sx q[3];
rz(-0.89151645) q[3];
sx q[3];
rz(0.69361544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2873516) q[0];
sx q[0];
rz(-1.1653405) q[0];
sx q[0];
rz(-1.5821138) q[0];
rz(-2.6782425) q[1];
sx q[1];
rz(-2.2644823) q[1];
sx q[1];
rz(1.5092441) q[1];
rz(0.91607416) q[2];
sx q[2];
rz(-1.6486042) q[2];
sx q[2];
rz(-0.83124607) q[2];
rz(3.0795931) q[3];
sx q[3];
rz(-1.0835032) q[3];
sx q[3];
rz(-0.54855357) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
