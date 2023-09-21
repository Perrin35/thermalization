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
rz(4.6580553) q[0];
sx q[0];
rz(9.1604995) q[0];
rz(2.1677986) q[1];
sx q[1];
rz(-1.9314613) q[1];
sx q[1];
rz(-0.73524737) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7289294) q[0];
sx q[0];
rz(-2.465415) q[0];
sx q[0];
rz(-0.23763188) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0636319) q[2];
sx q[2];
rz(-2.1858366) q[2];
sx q[2];
rz(-1.0539436) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.83598614) q[1];
sx q[1];
rz(-1.7251833) q[1];
sx q[1];
rz(1.2812213) q[1];
rz(2.5611476) q[3];
sx q[3];
rz(-1.547303) q[3];
sx q[3];
rz(-1.1593727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4720817) q[2];
sx q[2];
rz(-1.8410204) q[2];
sx q[2];
rz(-2.0377339) q[2];
rz(-1.2708698) q[3];
sx q[3];
rz(-1.2277675) q[3];
sx q[3];
rz(0.27403533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1141777) q[0];
sx q[0];
rz(-0.52863055) q[0];
sx q[0];
rz(-2.7052178) q[0];
rz(2.6787058) q[1];
sx q[1];
rz(-2.1040237) q[1];
sx q[1];
rz(-0.26611051) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14861815) q[0];
sx q[0];
rz(-2.2668112) q[0];
sx q[0];
rz(-0.50177411) q[0];
rz(-2.9782148) q[2];
sx q[2];
rz(-2.1401569) q[2];
sx q[2];
rz(0.53158224) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4437342) q[1];
sx q[1];
rz(-0.88159544) q[1];
sx q[1];
rz(-2.2976848) q[1];
rz(0.9640785) q[3];
sx q[3];
rz(-1.7214516) q[3];
sx q[3];
rz(-1.3552624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.46488547) q[2];
sx q[2];
rz(-1.2767982) q[2];
sx q[2];
rz(-0.51149386) q[2];
rz(2.3320847) q[3];
sx q[3];
rz(-1.5313238) q[3];
sx q[3];
rz(-2.8539343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-2.4354316) q[0];
sx q[0];
rz(-1.5937188) q[0];
sx q[0];
rz(-0.92873746) q[0];
rz(-1.4061032) q[1];
sx q[1];
rz(-2.4415253) q[1];
sx q[1];
rz(-1.3471289) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72235332) q[0];
sx q[0];
rz(-1.1142715) q[0];
sx q[0];
rz(1.8629575) q[0];
rz(-pi) q[1];
rz(-0.15813078) q[2];
sx q[2];
rz(-1.2889382) q[2];
sx q[2];
rz(-2.6513197) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9007064) q[1];
sx q[1];
rz(-0.97929231) q[1];
sx q[1];
rz(2.790931) q[1];
rz(-1.3942765) q[3];
sx q[3];
rz(-0.71478292) q[3];
sx q[3];
rz(-2.178758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9946263) q[2];
sx q[2];
rz(-1.3866084) q[2];
sx q[2];
rz(1.6195126) q[2];
rz(-2.8772723) q[3];
sx q[3];
rz(-1.0364573) q[3];
sx q[3];
rz(0.49595293) q[3];
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
rz(-2.3494444) q[0];
sx q[0];
rz(-1.148372) q[0];
sx q[0];
rz(0.96570063) q[0];
rz(-2.4194338) q[1];
sx q[1];
rz(-1.5042217) q[1];
sx q[1];
rz(2.5818363) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0069665466) q[0];
sx q[0];
rz(-0.43381938) q[0];
sx q[0];
rz(-0.19402786) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3279834) q[2];
sx q[2];
rz(-1.0126197) q[2];
sx q[2];
rz(1.8166325) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.30448118) q[1];
sx q[1];
rz(-1.5251535) q[1];
sx q[1];
rz(0.95506217) q[1];
rz(-1.3341321) q[3];
sx q[3];
rz(-0.79573123) q[3];
sx q[3];
rz(2.3468897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7136148) q[2];
sx q[2];
rz(-1.5087936) q[2];
sx q[2];
rz(-1.7170061) q[2];
rz(2.8811841) q[3];
sx q[3];
rz(-1.7491165) q[3];
sx q[3];
rz(-2.7105455) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99073064) q[0];
sx q[0];
rz(-1.6435511) q[0];
sx q[0];
rz(-0.36636233) q[0];
rz(1.5953966) q[1];
sx q[1];
rz(-0.55391824) q[1];
sx q[1];
rz(-2.7979134) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30003798) q[0];
sx q[0];
rz(-0.90111387) q[0];
sx q[0];
rz(0.71398736) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2720488) q[2];
sx q[2];
rz(-2.0336667) q[2];
sx q[2];
rz(0.16823828) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3949563) q[1];
sx q[1];
rz(-0.50368217) q[1];
sx q[1];
rz(-1.4175182) q[1];
rz(-pi) q[2];
rz(-0.92693365) q[3];
sx q[3];
rz(-1.5878131) q[3];
sx q[3];
rz(0.073236853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3917824) q[2];
sx q[2];
rz(-2.2822773) q[2];
sx q[2];
rz(3.0878477) q[2];
rz(-1.7371477) q[3];
sx q[3];
rz(-2.6440547) q[3];
sx q[3];
rz(-0.29156175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6546201) q[0];
sx q[0];
rz(-0.64990652) q[0];
sx q[0];
rz(0.75575954) q[0];
rz(0.02515633) q[1];
sx q[1];
rz(-0.92725602) q[1];
sx q[1];
rz(-2.8818534) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4070968) q[0];
sx q[0];
rz(-1.8844814) q[0];
sx q[0];
rz(-1.19019) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.01979205) q[2];
sx q[2];
rz(-1.2345825) q[2];
sx q[2];
rz(-0.83376955) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7163135) q[1];
sx q[1];
rz(-1.7967766) q[1];
sx q[1];
rz(1.2334361) q[1];
x q[2];
rz(1.1214244) q[3];
sx q[3];
rz(-1.5546038) q[3];
sx q[3];
rz(-1.8024973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6529237) q[2];
sx q[2];
rz(-2.3355464) q[2];
sx q[2];
rz(3.0409813) q[2];
rz(2.9595024) q[3];
sx q[3];
rz(-2.263335) q[3];
sx q[3];
rz(-1.3476868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-1.1473734) q[0];
sx q[0];
rz(-1.7213151) q[0];
sx q[0];
rz(-0.31016645) q[0];
rz(0.50225964) q[1];
sx q[1];
rz(-2.6175833) q[1];
sx q[1];
rz(0.60595864) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92354846) q[0];
sx q[0];
rz(-0.96203066) q[0];
sx q[0];
rz(1.2952842) q[0];
rz(-2.546054) q[2];
sx q[2];
rz(-0.34791246) q[2];
sx q[2];
rz(2.4728342) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1389321) q[1];
sx q[1];
rz(-1.9117038) q[1];
sx q[1];
rz(1.0374116) q[1];
x q[2];
rz(-0.80375399) q[3];
sx q[3];
rz(-2.440212) q[3];
sx q[3];
rz(-0.758981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9014152) q[2];
sx q[2];
rz(-1.9577953) q[2];
sx q[2];
rz(-0.85285464) q[2];
rz(1.3700221) q[3];
sx q[3];
rz(-1.6829237) q[3];
sx q[3];
rz(0.20496932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2296427) q[0];
sx q[0];
rz(-0.59597534) q[0];
sx q[0];
rz(1.461331) q[0];
rz(1.7386859) q[1];
sx q[1];
rz(-2.1673514) q[1];
sx q[1];
rz(-0.064037474) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2647576) q[0];
sx q[0];
rz(-1.3834073) q[0];
sx q[0];
rz(-2.1894356) q[0];
rz(-pi) q[1];
rz(1.0520347) q[2];
sx q[2];
rz(-0.56781893) q[2];
sx q[2];
rz(1.0713112) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9248326) q[1];
sx q[1];
rz(-1.2158582) q[1];
sx q[1];
rz(1.6096398) q[1];
rz(-1.2588345) q[3];
sx q[3];
rz(-1.6083816) q[3];
sx q[3];
rz(3.104044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0503851) q[2];
sx q[2];
rz(-0.63082266) q[2];
sx q[2];
rz(1.5554265) q[2];
rz(0.88820109) q[3];
sx q[3];
rz(-1.2398088) q[3];
sx q[3];
rz(0.92938882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.9973307) q[0];
sx q[0];
rz(-1.0634796) q[0];
sx q[0];
rz(-1.7096747) q[0];
rz(2.572708) q[1];
sx q[1];
rz(-0.535393) q[1];
sx q[1];
rz(1.127839) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.226798) q[0];
sx q[0];
rz(-0.1495805) q[0];
sx q[0];
rz(0.10766715) q[0];
x q[1];
rz(-2.2531829) q[2];
sx q[2];
rz(-0.75196224) q[2];
sx q[2];
rz(1.0351406) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.20967083) q[1];
sx q[1];
rz(-1.5591963) q[1];
sx q[1];
rz(-1.329122) q[1];
x q[2];
rz(-2.1209675) q[3];
sx q[3];
rz(-0.63311011) q[3];
sx q[3];
rz(1.8819295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.613712) q[2];
sx q[2];
rz(-2.2884559) q[2];
sx q[2];
rz(-2.9679427) q[2];
rz(-2.8052143) q[3];
sx q[3];
rz(-1.222638) q[3];
sx q[3];
rz(2.7500847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4353452) q[0];
sx q[0];
rz(-0.58723891) q[0];
sx q[0];
rz(-1.6760814) q[0];
rz(0.82410518) q[1];
sx q[1];
rz(-1.5287639) q[1];
sx q[1];
rz(2.5691659) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49003285) q[0];
sx q[0];
rz(-0.76793725) q[0];
sx q[0];
rz(0.72816531) q[0];
rz(-1.053439) q[2];
sx q[2];
rz(-0.51689076) q[2];
sx q[2];
rz(3.0992532) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.13092566) q[1];
sx q[1];
rz(-1.7867242) q[1];
sx q[1];
rz(-2.8972577) q[1];
rz(-pi) q[2];
rz(-2.4329348) q[3];
sx q[3];
rz(-1.346568) q[3];
sx q[3];
rz(1.538016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3616025) q[2];
sx q[2];
rz(-0.47452351) q[2];
sx q[2];
rz(-0.53722107) q[2];
rz(2.0843263) q[3];
sx q[3];
rz(-0.89151645) q[3];
sx q[3];
rz(0.69361544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.2873516) q[0];
sx q[0];
rz(-1.1653405) q[0];
sx q[0];
rz(-1.5821138) q[0];
rz(-0.46335012) q[1];
sx q[1];
rz(-0.87711038) q[1];
sx q[1];
rz(-1.6323485) q[1];
rz(1.4434545) q[2];
sx q[2];
rz(-2.4829395) q[2];
sx q[2];
rz(0.84045835) q[2];
rz(1.4544009) q[3];
sx q[3];
rz(-2.6506861) q[3];
sx q[3];
rz(-0.41674137) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
