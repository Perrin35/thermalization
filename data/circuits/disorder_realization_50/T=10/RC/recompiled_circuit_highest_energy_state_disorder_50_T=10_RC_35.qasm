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
rz(0.4515689) q[0];
sx q[0];
rz(-1.7234001) q[0];
sx q[0];
rz(0.42887846) q[0];
rz(-0.14313993) q[1];
sx q[1];
rz(-2.0506471) q[1];
sx q[1];
rz(0.080088869) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0319666) q[0];
sx q[0];
rz(-1.4916149) q[0];
sx q[0];
rz(-1.5575141) q[0];
x q[1];
rz(-1.2099138) q[2];
sx q[2];
rz(-1.6671125) q[2];
sx q[2];
rz(0.99260867) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7987207) q[1];
sx q[1];
rz(-1.3075383) q[1];
sx q[1];
rz(0.85204551) q[1];
rz(-pi) q[2];
rz(-0.24899068) q[3];
sx q[3];
rz(-2.5496229) q[3];
sx q[3];
rz(-2.0651061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0668138) q[2];
sx q[2];
rz(-2.0902233) q[2];
sx q[2];
rz(1.4220413) q[2];
rz(-1.7285796) q[3];
sx q[3];
rz(-1.541411) q[3];
sx q[3];
rz(-1.6933256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0195352) q[0];
sx q[0];
rz(-1.3966565) q[0];
sx q[0];
rz(0.33357093) q[0];
rz(0.16593274) q[1];
sx q[1];
rz(-0.3438147) q[1];
sx q[1];
rz(-0.75210345) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1307169) q[0];
sx q[0];
rz(-0.97144043) q[0];
sx q[0];
rz(-0.54979558) q[0];
x q[1];
rz(0.11801783) q[2];
sx q[2];
rz(-0.22023295) q[2];
sx q[2];
rz(1.756246) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1811407) q[1];
sx q[1];
rz(-1.6799095) q[1];
sx q[1];
rz(-1.0427598) q[1];
rz(-2.4346967) q[3];
sx q[3];
rz(-0.33875468) q[3];
sx q[3];
rz(-1.930463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0636474) q[2];
sx q[2];
rz(-2.2293978) q[2];
sx q[2];
rz(0.65703854) q[2];
rz(0.71197236) q[3];
sx q[3];
rz(-0.65198055) q[3];
sx q[3];
rz(2.3967192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7065358) q[0];
sx q[0];
rz(-2.5223795) q[0];
sx q[0];
rz(1.0414498) q[0];
rz(2.7667747) q[1];
sx q[1];
rz(-1.5033787) q[1];
sx q[1];
rz(-0.55694881) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3303012) q[0];
sx q[0];
rz(-2.3200703) q[0];
sx q[0];
rz(-2.6514451) q[0];
rz(-pi) q[1];
rz(-1.1191423) q[2];
sx q[2];
rz(-1.5030393) q[2];
sx q[2];
rz(-2.6299143) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8136602) q[1];
sx q[1];
rz(-0.90444649) q[1];
sx q[1];
rz(-1.7009237) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4905373) q[3];
sx q[3];
rz(-1.3811033) q[3];
sx q[3];
rz(-0.22471353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.53691429) q[2];
sx q[2];
rz(-0.92427173) q[2];
sx q[2];
rz(-0.43517819) q[2];
rz(-2.6168881) q[3];
sx q[3];
rz(-2.8080431) q[3];
sx q[3];
rz(-0.13478002) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0404469) q[0];
sx q[0];
rz(-2.0772159) q[0];
sx q[0];
rz(-1.5144298) q[0];
rz(2.0928404) q[1];
sx q[1];
rz(-0.92275134) q[1];
sx q[1];
rz(1.705816) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7483317) q[0];
sx q[0];
rz(-2.9611641) q[0];
sx q[0];
rz(2.9004651) q[0];
rz(-1.9210774) q[2];
sx q[2];
rz(-2.0868013) q[2];
sx q[2];
rz(-2.8523741) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9658815) q[1];
sx q[1];
rz(-2.615965) q[1];
sx q[1];
rz(1.4429343) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0346927) q[3];
sx q[3];
rz(-1.8827065) q[3];
sx q[3];
rz(2.8297765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.40020308) q[2];
sx q[2];
rz(-1.783739) q[2];
sx q[2];
rz(-3.1164361) q[2];
rz(1.4870421) q[3];
sx q[3];
rz(-0.39792037) q[3];
sx q[3];
rz(-2.3253843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
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
rz(0.095489278) q[0];
sx q[0];
rz(-0.64952055) q[0];
sx q[0];
rz(-0.15897861) q[0];
rz(-2.036463) q[1];
sx q[1];
rz(-0.72560328) q[1];
sx q[1];
rz(-0.22430688) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0353349) q[0];
sx q[0];
rz(-0.70959896) q[0];
sx q[0];
rz(-0.64340083) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0105466) q[2];
sx q[2];
rz(-1.59861) q[2];
sx q[2];
rz(-1.5459488) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.080706656) q[1];
sx q[1];
rz(-1.6861486) q[1];
sx q[1];
rz(-1.0447698) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1929694) q[3];
sx q[3];
rz(-2.5558839) q[3];
sx q[3];
rz(-1.0568413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.20532456) q[2];
sx q[2];
rz(-1.3621796) q[2];
sx q[2];
rz(-0.74860191) q[2];
rz(0.11048206) q[3];
sx q[3];
rz(-1.2275077) q[3];
sx q[3];
rz(2.7301679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1614138) q[0];
sx q[0];
rz(-0.36853376) q[0];
sx q[0];
rz(-0.73032105) q[0];
rz(1.5018564) q[1];
sx q[1];
rz(-1.7514936) q[1];
sx q[1];
rz(2.9072442) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6662009) q[0];
sx q[0];
rz(-2.0115543) q[0];
sx q[0];
rz(-2.4077364) q[0];
rz(-pi) q[1];
rz(-2.2147398) q[2];
sx q[2];
rz(-2.0704464) q[2];
sx q[2];
rz(2.814803) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1533656) q[1];
sx q[1];
rz(-1.4082137) q[1];
sx q[1];
rz(-1.6189727) q[1];
rz(-1.430351) q[3];
sx q[3];
rz(-1.5101119) q[3];
sx q[3];
rz(1.7588601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.948287) q[2];
sx q[2];
rz(-0.81393465) q[2];
sx q[2];
rz(1.8602271) q[2];
rz(-2.5365601) q[3];
sx q[3];
rz(-1.0993967) q[3];
sx q[3];
rz(-1.9090778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.989711) q[0];
sx q[0];
rz(-1.879377) q[0];
sx q[0];
rz(-2.5285517) q[0];
rz(-2.4609861) q[1];
sx q[1];
rz(-1.1111518) q[1];
sx q[1];
rz(-1.3253164) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8794358) q[0];
sx q[0];
rz(-0.61816321) q[0];
sx q[0];
rz(-2.6690041) q[0];
rz(-1.6540131) q[2];
sx q[2];
rz(-1.7832547) q[2];
sx q[2];
rz(2.0709289) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.054106709) q[1];
sx q[1];
rz(-2.1313868) q[1];
sx q[1];
rz(-1.7203548) q[1];
x q[2];
rz(0.57215055) q[3];
sx q[3];
rz(-0.66416262) q[3];
sx q[3];
rz(-0.45156878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.94720542) q[2];
sx q[2];
rz(-1.2014061) q[2];
sx q[2];
rz(0.80805937) q[2];
rz(-2.4957073) q[3];
sx q[3];
rz(-0.26791993) q[3];
sx q[3];
rz(2.8382235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5757669) q[0];
sx q[0];
rz(-0.72661) q[0];
sx q[0];
rz(-2.6988244) q[0];
rz(0.59204656) q[1];
sx q[1];
rz(-2.4202085) q[1];
sx q[1];
rz(-2.9659042) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4018049) q[0];
sx q[0];
rz(-1.2953838) q[0];
sx q[0];
rz(-2.7729183) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6106583) q[2];
sx q[2];
rz(-2.022036) q[2];
sx q[2];
rz(0.892837) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.86398839) q[1];
sx q[1];
rz(-0.93011198) q[1];
sx q[1];
rz(-1.6246269) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2288228) q[3];
sx q[3];
rz(-2.4173749) q[3];
sx q[3];
rz(-0.66373827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5438133) q[2];
sx q[2];
rz(-2.3976349) q[2];
sx q[2];
rz(-2.234484) q[2];
rz(-2.234327) q[3];
sx q[3];
rz(-1.7472569) q[3];
sx q[3];
rz(0.11769122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1842781) q[0];
sx q[0];
rz(-2.2397581) q[0];
sx q[0];
rz(-2.9515475) q[0];
rz(-1.5589145) q[1];
sx q[1];
rz(-1.546944) q[1];
sx q[1];
rz(1.83439) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0810159) q[0];
sx q[0];
rz(-1.9344442) q[0];
sx q[0];
rz(2.9878997) q[0];
rz(0.3385915) q[2];
sx q[2];
rz(-2.2420852) q[2];
sx q[2];
rz(0.012635144) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1550473) q[1];
sx q[1];
rz(-1.5380171) q[1];
sx q[1];
rz(1.4630586) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8966859) q[3];
sx q[3];
rz(-0.4462113) q[3];
sx q[3];
rz(-3.0065766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0731395) q[2];
sx q[2];
rz(-1.0871202) q[2];
sx q[2];
rz(-2.9404409) q[2];
rz(1.0226095) q[3];
sx q[3];
rz(-2.4590731) q[3];
sx q[3];
rz(-1.6020017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96711838) q[0];
sx q[0];
rz(-2.0956464) q[0];
sx q[0];
rz(2.2651267) q[0];
rz(2.5190952) q[1];
sx q[1];
rz(-0.21177706) q[1];
sx q[1];
rz(-2.7988787) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4349608) q[0];
sx q[0];
rz(-1.5049269) q[0];
sx q[0];
rz(0.94191282) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0194655) q[2];
sx q[2];
rz(-1.2363648) q[2];
sx q[2];
rz(1.8129406) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5610159) q[1];
sx q[1];
rz(-0.80272148) q[1];
sx q[1];
rz(-0.64508857) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8233772) q[3];
sx q[3];
rz(-1.1942721) q[3];
sx q[3];
rz(1.2600217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1891302) q[2];
sx q[2];
rz(-1.6996982) q[2];
sx q[2];
rz(-0.44378898) q[2];
rz(-2.022838) q[3];
sx q[3];
rz(-1.5951364) q[3];
sx q[3];
rz(-0.55383033) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2832058) q[0];
sx q[0];
rz(-1.5571742) q[0];
sx q[0];
rz(-1.5207186) q[0];
rz(-0.063477909) q[1];
sx q[1];
rz(-1.7964446) q[1];
sx q[1];
rz(-1.7367015) q[1];
rz(1.9012801) q[2];
sx q[2];
rz(-1.5122432) q[2];
sx q[2];
rz(1.6846364) q[2];
rz(2.1370112) q[3];
sx q[3];
rz(-1.3642642) q[3];
sx q[3];
rz(3.0929079) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
