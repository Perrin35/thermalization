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
rz(0.02778223) q[0];
sx q[0];
rz(-1.2454998) q[0];
sx q[0];
rz(1.0863289) q[0];
rz(2.0431986) q[1];
sx q[1];
rz(-1.9280704) q[1];
sx q[1];
rz(1.9953802) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6667696) q[0];
sx q[0];
rz(-0.36628977) q[0];
sx q[0];
rz(-1.1046871) q[0];
rz(1.97922) q[2];
sx q[2];
rz(-2.1222489) q[2];
sx q[2];
rz(1.1609032) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0806345) q[1];
sx q[1];
rz(-2.3976711) q[1];
sx q[1];
rz(0.77772909) q[1];
rz(-pi) q[2];
x q[2];
rz(0.2262874) q[3];
sx q[3];
rz(-2.0620966) q[3];
sx q[3];
rz(2.0313583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9132797) q[2];
sx q[2];
rz(-1.5756807) q[2];
sx q[2];
rz(2.7794465) q[2];
rz(0.95237887) q[3];
sx q[3];
rz(-0.55914545) q[3];
sx q[3];
rz(0.70498103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1700965) q[0];
sx q[0];
rz(-0.2146475) q[0];
sx q[0];
rz(-1.6098518) q[0];
rz(1.8929298) q[1];
sx q[1];
rz(-1.5208533) q[1];
sx q[1];
rz(-2.2889287) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8448317) q[0];
sx q[0];
rz(-1.0687901) q[0];
sx q[0];
rz(0.39433483) q[0];
x q[1];
rz(-0.12262362) q[2];
sx q[2];
rz(-1.2957591) q[2];
sx q[2];
rz(3.1240535) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.98140796) q[1];
sx q[1];
rz(-1.7716737) q[1];
sx q[1];
rz(-2.1698879) q[1];
rz(-pi) q[2];
x q[2];
rz(2.144935) q[3];
sx q[3];
rz(-0.64145815) q[3];
sx q[3];
rz(1.5456089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2109005) q[2];
sx q[2];
rz(-2.1810668) q[2];
sx q[2];
rz(1.1962013) q[2];
rz(-3.1148552) q[3];
sx q[3];
rz(-0.49632159) q[3];
sx q[3];
rz(1.668821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6579599) q[0];
sx q[0];
rz(-0.25068972) q[0];
sx q[0];
rz(0.87580097) q[0];
rz(-0.35975131) q[1];
sx q[1];
rz(-1.7774372) q[1];
sx q[1];
rz(-2.13805) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8784766) q[0];
sx q[0];
rz(-2.5777393) q[0];
sx q[0];
rz(3.0082358) q[0];
rz(0.92873145) q[2];
sx q[2];
rz(-1.8936529) q[2];
sx q[2];
rz(-1.0869458) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7454804) q[1];
sx q[1];
rz(-2.2428233) q[1];
sx q[1];
rz(-1.7538096) q[1];
x q[2];
rz(0.67387132) q[3];
sx q[3];
rz(-2.0099075) q[3];
sx q[3];
rz(-1.278745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.2063736) q[2];
sx q[2];
rz(-0.71722764) q[2];
sx q[2];
rz(-0.75500542) q[2];
rz(0.098091789) q[3];
sx q[3];
rz(-2.5710929) q[3];
sx q[3];
rz(0.3977972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0574684) q[0];
sx q[0];
rz(-1.5115154) q[0];
sx q[0];
rz(-0.64495069) q[0];
rz(1.0954674) q[1];
sx q[1];
rz(-2.7332833) q[1];
sx q[1];
rz(2.0300949) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43954014) q[0];
sx q[0];
rz(-1.6344392) q[0];
sx q[0];
rz(-1.511277) q[0];
rz(-1.4884645) q[2];
sx q[2];
rz(-1.4042743) q[2];
sx q[2];
rz(2.4660267) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7202338) q[1];
sx q[1];
rz(-1.3610657) q[1];
sx q[1];
rz(1.7214107) q[1];
x q[2];
rz(2.901279) q[3];
sx q[3];
rz(-2.3860768) q[3];
sx q[3];
rz(0.66080581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5081386) q[2];
sx q[2];
rz(-2.3544669) q[2];
sx q[2];
rz(-2.656929) q[2];
rz(-0.44114068) q[3];
sx q[3];
rz(-1.7287247) q[3];
sx q[3];
rz(1.2539366) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6345374) q[0];
sx q[0];
rz(-2.2978954) q[0];
sx q[0];
rz(-1.7150568) q[0];
rz(-0.85975319) q[1];
sx q[1];
rz(-2.8883002) q[1];
sx q[1];
rz(-2.3111129) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7380421) q[0];
sx q[0];
rz(-1.0071686) q[0];
sx q[0];
rz(0.72569816) q[0];
rz(-2.203622) q[2];
sx q[2];
rz(-1.5245078) q[2];
sx q[2];
rz(0.37424311) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3871206) q[1];
sx q[1];
rz(-0.74229147) q[1];
sx q[1];
rz(0.20739389) q[1];
rz(-1.4111341) q[3];
sx q[3];
rz(-1.214802) q[3];
sx q[3];
rz(1.0245263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6972203) q[2];
sx q[2];
rz(-1.0151007) q[2];
sx q[2];
rz(-2.3946136) q[2];
rz(-2.8367481) q[3];
sx q[3];
rz(-1.3957142) q[3];
sx q[3];
rz(1.9029118) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9656669) q[0];
sx q[0];
rz(-1.8400064) q[0];
sx q[0];
rz(0.23342215) q[0];
rz(0.4625136) q[1];
sx q[1];
rz(-1.1902483) q[1];
sx q[1];
rz(-2.7375284) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8916805) q[0];
sx q[0];
rz(-0.95522987) q[0];
sx q[0];
rz(2.8717442) q[0];
x q[1];
rz(-2.7163804) q[2];
sx q[2];
rz(-1.325144) q[2];
sx q[2];
rz(0.52580357) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7884454) q[1];
sx q[1];
rz(-0.57203469) q[1];
sx q[1];
rz(2.9414909) q[1];
rz(-0.63066354) q[3];
sx q[3];
rz(-1.3913097) q[3];
sx q[3];
rz(1.6398182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3290951) q[2];
sx q[2];
rz(-1.3500682) q[2];
sx q[2];
rz(2.4883032) q[2];
rz(-1.6434044) q[3];
sx q[3];
rz(-2.523962) q[3];
sx q[3];
rz(0.60630715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6752328) q[0];
sx q[0];
rz(-1.9283858) q[0];
sx q[0];
rz(-2.8208222) q[0];
rz(2.4877211) q[1];
sx q[1];
rz(-2.8120698) q[1];
sx q[1];
rz(-0.060506495) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4757547) q[0];
sx q[0];
rz(-1.3011969) q[0];
sx q[0];
rz(-0.17335391) q[0];
rz(-pi) q[1];
rz(0.67530151) q[2];
sx q[2];
rz(-0.85682288) q[2];
sx q[2];
rz(-1.7973822) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6937032) q[1];
sx q[1];
rz(-0.9606908) q[1];
sx q[1];
rz(-2.542716) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2809847) q[3];
sx q[3];
rz(-2.6250556) q[3];
sx q[3];
rz(2.9252657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.74762809) q[2];
sx q[2];
rz(-1.4614033) q[2];
sx q[2];
rz(0.72236577) q[2];
rz(0.28313053) q[3];
sx q[3];
rz(-0.88909283) q[3];
sx q[3];
rz(-2.6319671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79001456) q[0];
sx q[0];
rz(-0.55792648) q[0];
sx q[0];
rz(-2.956692) q[0];
rz(2.6376873) q[1];
sx q[1];
rz(-1.7920707) q[1];
sx q[1];
rz(2.6769743) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9757742) q[0];
sx q[0];
rz(-1.4041472) q[0];
sx q[0];
rz(-1.7383132) q[0];
rz(-pi) q[1];
rz(1.4746998) q[2];
sx q[2];
rz(-1.155174) q[2];
sx q[2];
rz(-1.9643652) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.020430001) q[1];
sx q[1];
rz(-1.1926225) q[1];
sx q[1];
rz(2.6932823) q[1];
rz(1.7691602) q[3];
sx q[3];
rz(-0.86097211) q[3];
sx q[3];
rz(1.7539303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9087312) q[2];
sx q[2];
rz(-1.6249012) q[2];
sx q[2];
rz(2.0880584) q[2];
rz(0.72707027) q[3];
sx q[3];
rz(-1.2289685) q[3];
sx q[3];
rz(-0.54522902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9497249) q[0];
sx q[0];
rz(-1.9504915) q[0];
sx q[0];
rz(-0.18549347) q[0];
rz(2.5718451) q[1];
sx q[1];
rz(-0.92420095) q[1];
sx q[1];
rz(1.1572908) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13979736) q[0];
sx q[0];
rz(-1.9855864) q[0];
sx q[0];
rz(-1.7263361) q[0];
x q[1];
rz(-3.0915843) q[2];
sx q[2];
rz(-2.2029624) q[2];
sx q[2];
rz(0.061246733) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.696582) q[1];
sx q[1];
rz(-0.60234501) q[1];
sx q[1];
rz(0.38378832) q[1];
rz(-pi) q[2];
x q[2];
rz(0.80482153) q[3];
sx q[3];
rz(-1.8880883) q[3];
sx q[3];
rz(-2.1498093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7374604) q[2];
sx q[2];
rz(-2.9170333) q[2];
sx q[2];
rz(-1.4004716) q[2];
rz(0.36462668) q[3];
sx q[3];
rz(-2.3783763) q[3];
sx q[3];
rz(-2.3188685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51065651) q[0];
sx q[0];
rz(-0.83844227) q[0];
sx q[0];
rz(3.1070218) q[0];
rz(2.8055387) q[1];
sx q[1];
rz(-2.2269109) q[1];
sx q[1];
rz(-0.58759442) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4776992) q[0];
sx q[0];
rz(-1.4283533) q[0];
sx q[0];
rz(2.7201128) q[0];
rz(0.73070902) q[2];
sx q[2];
rz(-1.3255505) q[2];
sx q[2];
rz(-2.2718549) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.42660824) q[1];
sx q[1];
rz(-2.9054925) q[1];
sx q[1];
rz(-1.3900422) q[1];
rz(2.4617599) q[3];
sx q[3];
rz(-1.203591) q[3];
sx q[3];
rz(2.1757371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8079638) q[2];
sx q[2];
rz(-0.89666349) q[2];
sx q[2];
rz(2.4133852) q[2];
rz(-2.7306469) q[3];
sx q[3];
rz(-0.22308895) q[3];
sx q[3];
rz(-1.2724916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67644607) q[0];
sx q[0];
rz(-1.5560173) q[0];
sx q[0];
rz(-1.7898855) q[0];
rz(0.43988718) q[1];
sx q[1];
rz(-0.37275795) q[1];
sx q[1];
rz(2.8511924) q[1];
rz(-0.34392233) q[2];
sx q[2];
rz(-1.2871598) q[2];
sx q[2];
rz(-3.1393307) q[2];
rz(-0.10671814) q[3];
sx q[3];
rz(-1.5244457) q[3];
sx q[3];
rz(2.5929858) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
