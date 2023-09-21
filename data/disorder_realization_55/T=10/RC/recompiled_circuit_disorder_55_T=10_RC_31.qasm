OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5947333) q[0];
sx q[0];
rz(-1.5164627) q[0];
sx q[0];
rz(-2.8773142) q[0];
rz(2.1677986) q[1];
sx q[1];
rz(-1.9314613) q[1];
sx q[1];
rz(2.4063453) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7289294) q[0];
sx q[0];
rz(-0.67617765) q[0];
sx q[0];
rz(2.9039608) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0779607) q[2];
sx q[2];
rz(-0.9557561) q[2];
sx q[2];
rz(1.0539436) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3056065) q[1];
sx q[1];
rz(-1.4164093) q[1];
sx q[1];
rz(1.2812213) q[1];
x q[2];
rz(0.04282184) q[3];
sx q[3];
rz(-0.58086568) q[3];
sx q[3];
rz(2.7659741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4720817) q[2];
sx q[2];
rz(-1.8410204) q[2];
sx q[2];
rz(-1.1038587) q[2];
rz(1.2708698) q[3];
sx q[3];
rz(-1.9138252) q[3];
sx q[3];
rz(-2.8675573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1141777) q[0];
sx q[0];
rz(-0.52863055) q[0];
sx q[0];
rz(2.7052178) q[0];
rz(0.46288681) q[1];
sx q[1];
rz(-1.0375689) q[1];
sx q[1];
rz(-0.26611051) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3811831) q[0];
sx q[0];
rz(-1.1927483) q[0];
sx q[0];
rz(-0.80947431) q[0];
rz(-0.16337784) q[2];
sx q[2];
rz(-2.1401569) q[2];
sx q[2];
rz(2.6100104) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6978584) q[1];
sx q[1];
rz(-0.88159544) q[1];
sx q[1];
rz(-0.84390784) q[1];
rz(-pi) q[2];
rz(-2.9588685) q[3];
sx q[3];
rz(-0.97191873) q[3];
sx q[3];
rz(3.0298508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6767072) q[2];
sx q[2];
rz(-1.2767982) q[2];
sx q[2];
rz(2.6300988) q[2];
rz(-2.3320847) q[3];
sx q[3];
rz(-1.5313238) q[3];
sx q[3];
rz(2.8539343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4354316) q[0];
sx q[0];
rz(-1.5937188) q[0];
sx q[0];
rz(-2.2128552) q[0];
rz(1.4061032) q[1];
sx q[1];
rz(-2.4415253) q[1];
sx q[1];
rz(1.3471289) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0179694) q[0];
sx q[0];
rz(-0.53640134) q[0];
sx q[0];
rz(2.6111952) q[0];
rz(1.2855661) q[2];
sx q[2];
rz(-1.4189548) q[2];
sx q[2];
rz(1.0361995) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.24088629) q[1];
sx q[1];
rz(-2.1623003) q[1];
sx q[1];
rz(-2.790931) q[1];
rz(1.3942765) q[3];
sx q[3];
rz(-0.71478292) q[3];
sx q[3];
rz(-0.96283462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9946263) q[2];
sx q[2];
rz(-1.7549843) q[2];
sx q[2];
rz(-1.6195126) q[2];
rz(-0.26432031) q[3];
sx q[3];
rz(-2.1051354) q[3];
sx q[3];
rz(-2.6456397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79214823) q[0];
sx q[0];
rz(-1.9932207) q[0];
sx q[0];
rz(-2.175892) q[0];
rz(-2.4194338) q[1];
sx q[1];
rz(-1.5042217) q[1];
sx q[1];
rz(2.5818363) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9352919) q[0];
sx q[0];
rz(-1.1456523) q[0];
sx q[0];
rz(1.6598808) q[0];
x q[1];
rz(2.3279834) q[2];
sx q[2];
rz(-1.0126197) q[2];
sx q[2];
rz(1.3249601) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.30448118) q[1];
sx q[1];
rz(-1.5251535) q[1];
sx q[1];
rz(2.1865305) q[1];
rz(-pi) q[2];
rz(1.3341321) q[3];
sx q[3];
rz(-2.3458614) q[3];
sx q[3];
rz(-0.79470293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7136148) q[2];
sx q[2];
rz(-1.6327991) q[2];
sx q[2];
rz(1.4245865) q[2];
rz(-0.26040855) q[3];
sx q[3];
rz(-1.3924761) q[3];
sx q[3];
rz(2.7105455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99073064) q[0];
sx q[0];
rz(-1.6435511) q[0];
sx q[0];
rz(-2.7752303) q[0];
rz(1.5461961) q[1];
sx q[1];
rz(-2.5876744) q[1];
sx q[1];
rz(-2.7979134) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8415547) q[0];
sx q[0];
rz(-0.90111387) q[0];
sx q[0];
rz(2.4276053) q[0];
rz(-pi) q[1];
x q[1];
rz(0.86954388) q[2];
sx q[2];
rz(-1.1079259) q[2];
sx q[2];
rz(-0.16823828) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3949563) q[1];
sx q[1];
rz(-0.50368217) q[1];
sx q[1];
rz(1.4175182) q[1];
rz(-0.92693365) q[3];
sx q[3];
rz(-1.5537795) q[3];
sx q[3];
rz(-0.073236853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7498103) q[2];
sx q[2];
rz(-0.85931531) q[2];
sx q[2];
rz(0.053744944) q[2];
rz(1.7371477) q[3];
sx q[3];
rz(-0.497538) q[3];
sx q[3];
rz(-0.29156175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
rz(-0.4869726) q[0];
sx q[0];
rz(-2.4916861) q[0];
sx q[0];
rz(-0.75575954) q[0];
rz(0.02515633) q[1];
sx q[1];
rz(-0.92725602) q[1];
sx q[1];
rz(0.25973928) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.320967) q[0];
sx q[0];
rz(-2.6532986) q[0];
sx q[0];
rz(0.85296209) q[0];
x q[1];
rz(0.01979205) q[2];
sx q[2];
rz(-1.9070101) q[2];
sx q[2];
rz(-0.83376955) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.427634) q[1];
sx q[1];
rz(-2.7379588) q[1];
sx q[1];
rz(0.9637109) q[1];
rz(1.1214244) q[3];
sx q[3];
rz(-1.5546038) q[3];
sx q[3];
rz(-1.8024973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6529237) q[2];
sx q[2];
rz(-2.3355464) q[2];
sx q[2];
rz(0.10061131) q[2];
rz(2.9595024) q[3];
sx q[3];
rz(-0.87825769) q[3];
sx q[3];
rz(1.3476868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9942193) q[0];
sx q[0];
rz(-1.4202776) q[0];
sx q[0];
rz(-0.31016645) q[0];
rz(-2.639333) q[1];
sx q[1];
rz(-2.6175833) q[1];
sx q[1];
rz(-2.535634) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92354846) q[0];
sx q[0];
rz(-2.179562) q[0];
sx q[0];
rz(1.8463085) q[0];
rz(2.546054) q[2];
sx q[2];
rz(-2.7936802) q[2];
sx q[2];
rz(-0.66875848) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1389321) q[1];
sx q[1];
rz(-1.9117038) q[1];
sx q[1];
rz(1.0374116) q[1];
rz(-2.6113854) q[3];
sx q[3];
rz(-1.0876417) q[3];
sx q[3];
rz(0.14164856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.24017748) q[2];
sx q[2];
rz(-1.9577953) q[2];
sx q[2];
rz(0.85285464) q[2];
rz(-1.7715706) q[3];
sx q[3];
rz(-1.4586689) q[3];
sx q[3];
rz(2.9366233) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2296427) q[0];
sx q[0];
rz(-2.5456173) q[0];
sx q[0];
rz(1.6802616) q[0];
rz(1.4029067) q[1];
sx q[1];
rz(-2.1673514) q[1];
sx q[1];
rz(-3.0775552) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2647576) q[0];
sx q[0];
rz(-1.7581853) q[0];
sx q[0];
rz(2.1894356) q[0];
rz(-1.0520347) q[2];
sx q[2];
rz(-2.5737737) q[2];
sx q[2];
rz(-2.0702814) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9248326) q[1];
sx q[1];
rz(-1.9257345) q[1];
sx q[1];
rz(-1.5319529) q[1];
rz(-3.1021032) q[3];
sx q[3];
rz(-1.8825304) q[3];
sx q[3];
rz(1.5962275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.091207592) q[2];
sx q[2];
rz(-0.63082266) q[2];
sx q[2];
rz(1.5861661) q[2];
rz(2.2533916) q[3];
sx q[3];
rz(-1.9017838) q[3];
sx q[3];
rz(0.92938882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-2.6061997) q[1];
sx q[1];
rz(-1.127839) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9040684) q[0];
sx q[0];
rz(-1.5868109) q[0];
sx q[0];
rz(-0.1487271) q[0];
x q[1];
rz(-0.94294195) q[2];
sx q[2];
rz(-1.1254416) q[2];
sx q[2];
rz(3.359059e-05) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3582663) q[1];
sx q[1];
rz(-1.3291385) q[1];
sx q[1];
rz(-0.011947167) q[1];
x q[2];
rz(0.3663775) q[3];
sx q[3];
rz(-1.0421703) q[3];
sx q[3];
rz(-0.60929326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.613712) q[2];
sx q[2];
rz(-2.2884559) q[2];
sx q[2];
rz(-0.17364994) q[2];
rz(-2.8052143) q[3];
sx q[3];
rz(-1.9189546) q[3];
sx q[3];
rz(0.39150795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7062475) q[0];
sx q[0];
rz(-2.5543537) q[0];
sx q[0];
rz(-1.4655112) q[0];
rz(2.3174875) q[1];
sx q[1];
rz(-1.5287639) q[1];
sx q[1];
rz(0.5724268) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3819645) q[0];
sx q[0];
rz(-1.0257162) q[0];
sx q[0];
rz(-2.1419924) q[0];
rz(-2.0881537) q[2];
sx q[2];
rz(-2.6247019) q[2];
sx q[2];
rz(-0.042339485) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.13092566) q[1];
sx q[1];
rz(-1.7867242) q[1];
sx q[1];
rz(-0.24433498) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.33705538) q[3];
sx q[3];
rz(-0.73738499) q[3];
sx q[3];
rz(2.9205703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.77999014) q[2];
sx q[2];
rz(-0.47452351) q[2];
sx q[2];
rz(2.6043716) q[2];
rz(1.0572664) q[3];
sx q[3];
rz(-0.89151645) q[3];
sx q[3];
rz(-0.69361544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2873516) q[0];
sx q[0];
rz(-1.9762522) q[0];
sx q[0];
rz(1.5594788) q[0];
rz(0.46335012) q[1];
sx q[1];
rz(-2.2644823) q[1];
sx q[1];
rz(1.5092441) q[1];
rz(-1.6981381) q[2];
sx q[2];
rz(-2.4829395) q[2];
sx q[2];
rz(0.84045835) q[2];
rz(-1.6871917) q[3];
sx q[3];
rz(-2.6506861) q[3];
sx q[3];
rz(-0.41674137) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];