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
rz(1.4690118) q[0];
sx q[0];
rz(-0.91977683) q[0];
sx q[0];
rz(0.65499175) q[0];
rz(1.6545777) q[1];
sx q[1];
rz(4.5643212) q[1];
sx q[1];
rz(8.1017452) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1974063) q[0];
sx q[0];
rz(-1.8776769) q[0];
sx q[0];
rz(2.5635864) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8048067) q[2];
sx q[2];
rz(-1.26097) q[2];
sx q[2];
rz(1.7676304) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2046584) q[1];
sx q[1];
rz(-1.8465252) q[1];
sx q[1];
rz(3.0380121) q[1];
rz(-2.4134656) q[3];
sx q[3];
rz(-1.1525197) q[3];
sx q[3];
rz(-2.3969216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.46372867) q[2];
sx q[2];
rz(-2.4552796) q[2];
sx q[2];
rz(0.65832552) q[2];
rz(2.5074734) q[3];
sx q[3];
rz(-1.442797) q[3];
sx q[3];
rz(-1.3707976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(3.051801) q[0];
sx q[0];
rz(-0.3312411) q[0];
sx q[0];
rz(-0.077089699) q[0];
rz(-0.58473051) q[1];
sx q[1];
rz(-1.3405776) q[1];
sx q[1];
rz(-1.9416521) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8611885) q[0];
sx q[0];
rz(-1.1139835) q[0];
sx q[0];
rz(-0.21976075) q[0];
x q[1];
rz(3.1121677) q[2];
sx q[2];
rz(-1.4064854) q[2];
sx q[2];
rz(1.2677594) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5189831) q[1];
sx q[1];
rz(-0.72777343) q[1];
sx q[1];
rz(-0.98811291) q[1];
rz(-1.6317815) q[3];
sx q[3];
rz(-1.012523) q[3];
sx q[3];
rz(-1.3577485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.2526907) q[2];
sx q[2];
rz(-0.98068792) q[2];
sx q[2];
rz(-2.4661031) q[2];
rz(-0.0096970079) q[3];
sx q[3];
rz(-0.24295013) q[3];
sx q[3];
rz(-0.01827904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6460687) q[0];
sx q[0];
rz(-1.6697474) q[0];
sx q[0];
rz(-2.38696) q[0];
rz(-0.010628788) q[1];
sx q[1];
rz(-1.7517215) q[1];
sx q[1];
rz(2.1307814) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6957798) q[0];
sx q[0];
rz(-0.58923972) q[0];
sx q[0];
rz(2.1408098) q[0];
rz(-0.88917261) q[2];
sx q[2];
rz(-1.5487567) q[2];
sx q[2];
rz(0.88543788) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.52151187) q[1];
sx q[1];
rz(-2.3151868) q[1];
sx q[1];
rz(1.6537731) q[1];
x q[2];
rz(-0.63186462) q[3];
sx q[3];
rz(-0.87040983) q[3];
sx q[3];
rz(-2.4720342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2406771) q[2];
sx q[2];
rz(-2.7132576) q[2];
sx q[2];
rz(0.50673103) q[2];
rz(1.8690551) q[3];
sx q[3];
rz(-1.3576018) q[3];
sx q[3];
rz(0.27819628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5946567) q[0];
sx q[0];
rz(-1.2242067) q[0];
sx q[0];
rz(-1.1154255) q[0];
rz(-1.2660654) q[1];
sx q[1];
rz(-1.5645809) q[1];
sx q[1];
rz(2.26684) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0698358) q[0];
sx q[0];
rz(-1.6594145) q[0];
sx q[0];
rz(-1.8231234) q[0];
x q[1];
rz(-0.018847887) q[2];
sx q[2];
rz(-2.224378) q[2];
sx q[2];
rz(-3.1342521) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.03434788) q[1];
sx q[1];
rz(-2.6470091) q[1];
sx q[1];
rz(-2.5550708) q[1];
rz(-pi) q[2];
x q[2];
rz(2.277209) q[3];
sx q[3];
rz(-1.7076075) q[3];
sx q[3];
rz(-2.9446326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0227585) q[2];
sx q[2];
rz(-2.6648882) q[2];
sx q[2];
rz(1.9258707) q[2];
rz(1.7440965) q[3];
sx q[3];
rz(-1.5236866) q[3];
sx q[3];
rz(-1.4309179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8656411) q[0];
sx q[0];
rz(-0.10198274) q[0];
sx q[0];
rz(1.4707461) q[0];
rz(-1.7631081) q[1];
sx q[1];
rz(-0.78796092) q[1];
sx q[1];
rz(1.7313622) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5263162) q[0];
sx q[0];
rz(-1.7782167) q[0];
sx q[0];
rz(-2.0735969) q[0];
rz(2.8431434) q[2];
sx q[2];
rz(-2.4458234) q[2];
sx q[2];
rz(-1.4571822) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8961337) q[1];
sx q[1];
rz(-1.8832379) q[1];
sx q[1];
rz(-1.5691084) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.35378176) q[3];
sx q[3];
rz(-1.980034) q[3];
sx q[3];
rz(1.080846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0146279) q[2];
sx q[2];
rz(-0.85863272) q[2];
sx q[2];
rz(1.3237759) q[2];
rz(1.3592367) q[3];
sx q[3];
rz(-1.2991354) q[3];
sx q[3];
rz(-2.7544379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-1.0303665) q[0];
sx q[0];
rz(-1.68196) q[0];
sx q[0];
rz(0.02221814) q[0];
rz(-0.70023316) q[1];
sx q[1];
rz(-1.8902238) q[1];
sx q[1];
rz(0.95692316) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6895804) q[0];
sx q[0];
rz(-0.71419026) q[0];
sx q[0];
rz(0.22269188) q[0];
rz(-1.3889203) q[2];
sx q[2];
rz(-1.3564912) q[2];
sx q[2];
rz(-0.3265115) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.14084223) q[1];
sx q[1];
rz(-2.6472808) q[1];
sx q[1];
rz(-2.1782081) q[1];
rz(-pi) q[2];
rz(0.33958667) q[3];
sx q[3];
rz(-1.1533779) q[3];
sx q[3];
rz(-2.0805032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2172829) q[2];
sx q[2];
rz(-2.0192912) q[2];
sx q[2];
rz(-1.4420606) q[2];
rz(-1.1066655) q[3];
sx q[3];
rz(-2.536074) q[3];
sx q[3];
rz(1.1520011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7149413) q[0];
sx q[0];
rz(-1.8639257) q[0];
sx q[0];
rz(-2.6746124) q[0];
rz(1.3106208) q[1];
sx q[1];
rz(-0.95342031) q[1];
sx q[1];
rz(0.39042815) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4425499) q[0];
sx q[0];
rz(-0.14420284) q[0];
sx q[0];
rz(-0.22871916) q[0];
x q[1];
rz(1.2762345) q[2];
sx q[2];
rz(-1.6737564) q[2];
sx q[2];
rz(-1.1338794) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5891287) q[1];
sx q[1];
rz(-0.93163449) q[1];
sx q[1];
rz(-2.5659849) q[1];
x q[2];
rz(-2.9938145) q[3];
sx q[3];
rz(-0.96333233) q[3];
sx q[3];
rz(-1.0318508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9499669) q[2];
sx q[2];
rz(-1.9248631) q[2];
sx q[2];
rz(1.0170271) q[2];
rz(0.93938604) q[3];
sx q[3];
rz(-0.87299577) q[3];
sx q[3];
rz(0.92923195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7153213) q[0];
sx q[0];
rz(-2.4763698) q[0];
sx q[0];
rz(-1.2128879) q[0];
rz(-0.60316482) q[1];
sx q[1];
rz(-2.0096571) q[1];
sx q[1];
rz(-0.42339465) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.903666) q[0];
sx q[0];
rz(-2.4366424) q[0];
sx q[0];
rz(-1.4506819) q[0];
rz(-pi) q[1];
rz(3.064761) q[2];
sx q[2];
rz(-1.72121) q[2];
sx q[2];
rz(0.39932775) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0871157) q[1];
sx q[1];
rz(-1.7810139) q[1];
sx q[1];
rz(1.3940548) q[1];
rz(-2.9717507) q[3];
sx q[3];
rz(-1.563493) q[3];
sx q[3];
rz(-0.23582349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6508871) q[2];
sx q[2];
rz(-0.31279534) q[2];
sx q[2];
rz(2.1631964) q[2];
rz(-2.4466416) q[3];
sx q[3];
rz(-1.9415104) q[3];
sx q[3];
rz(-1.8470701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.097718358) q[0];
sx q[0];
rz(-0.23625034) q[0];
sx q[0];
rz(3.0978715) q[0];
rz(-1.1737191) q[1];
sx q[1];
rz(-0.81467384) q[1];
sx q[1];
rz(2.3497605) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7515562) q[0];
sx q[0];
rz(-1.9851713) q[0];
sx q[0];
rz(-0.66585559) q[0];
rz(-pi) q[1];
rz(0.20251198) q[2];
sx q[2];
rz(-1.2767451) q[2];
sx q[2];
rz(0.45627764) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.83306317) q[1];
sx q[1];
rz(-1.2385784) q[1];
sx q[1];
rz(1.1221207) q[1];
rz(-pi) q[2];
x q[2];
rz(0.65483396) q[3];
sx q[3];
rz(-2.5635701) q[3];
sx q[3];
rz(0.27930799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9110079) q[2];
sx q[2];
rz(-2.8318996) q[2];
sx q[2];
rz(-1.2370375) q[2];
rz(0.48219314) q[3];
sx q[3];
rz(-2.2174406) q[3];
sx q[3];
rz(-2.6528416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1173387) q[0];
sx q[0];
rz(-2.2145705) q[0];
sx q[0];
rz(1.8035969) q[0];
rz(2.2894739) q[1];
sx q[1];
rz(-1.2282649) q[1];
sx q[1];
rz(-1.9911912) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8088732) q[0];
sx q[0];
rz(-2.3925663) q[0];
sx q[0];
rz(-2.3179884) q[0];
rz(-pi) q[1];
rz(-1.2537635) q[2];
sx q[2];
rz(-1.6765127) q[2];
sx q[2];
rz(0.52717613) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.89429606) q[1];
sx q[1];
rz(-1.1793696) q[1];
sx q[1];
rz(-2.9904234) q[1];
x q[2];
rz(-0.97148599) q[3];
sx q[3];
rz(-1.2462292) q[3];
sx q[3];
rz(-2.6443986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.77356768) q[2];
sx q[2];
rz(-2.1358392) q[2];
sx q[2];
rz(2.8311938) q[2];
rz(-0.6271022) q[3];
sx q[3];
rz(-0.98395625) q[3];
sx q[3];
rz(1.1712801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8695759) q[0];
sx q[0];
rz(-1.5740812) q[0];
sx q[0];
rz(1.5480702) q[0];
rz(1.2177474) q[1];
sx q[1];
rz(-1.6883313) q[1];
sx q[1];
rz(2.2399166) q[1];
rz(-2.9318342) q[2];
sx q[2];
rz(-2.0611613) q[2];
sx q[2];
rz(0.47917889) q[2];
rz(-1.6321833) q[3];
sx q[3];
rz(-1.2493549) q[3];
sx q[3];
rz(-3.1102104) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
