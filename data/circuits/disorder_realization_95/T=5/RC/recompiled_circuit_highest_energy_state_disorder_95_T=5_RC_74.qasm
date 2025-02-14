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
rz(1.3746102) q[0];
sx q[0];
rz(-2.1687825) q[0];
sx q[0];
rz(1.2306124) q[0];
rz(-1.542701) q[1];
sx q[1];
rz(-0.35294947) q[1];
sx q[1];
rz(0.43500873) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4458076) q[0];
sx q[0];
rz(-1.2381993) q[0];
sx q[0];
rz(3.0363054) q[0];
x q[1];
rz(-1.4392008) q[2];
sx q[2];
rz(-1.7965266) q[2];
sx q[2];
rz(0.44039681) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8486464) q[1];
sx q[1];
rz(-0.82445626) q[1];
sx q[1];
rz(2.7566977) q[1];
rz(-1.1633662) q[3];
sx q[3];
rz(-3.1071783) q[3];
sx q[3];
rz(-0.84190166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1970485) q[2];
sx q[2];
rz(-2.7272482) q[2];
sx q[2];
rz(-2.22331) q[2];
rz(2.7783172) q[3];
sx q[3];
rz(-1.2493635) q[3];
sx q[3];
rz(-1.8173119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0488224) q[0];
sx q[0];
rz(-0.73960441) q[0];
sx q[0];
rz(-2.0059465) q[0];
rz(-2.8593235) q[1];
sx q[1];
rz(-1.6259364) q[1];
sx q[1];
rz(2.938882) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36798635) q[0];
sx q[0];
rz(-1.752089) q[0];
sx q[0];
rz(1.7530935) q[0];
rz(1.545867) q[2];
sx q[2];
rz(-0.48183107) q[2];
sx q[2];
rz(2.3584751) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7248279) q[1];
sx q[1];
rz(-1.6931122) q[1];
sx q[1];
rz(0.45682795) q[1];
rz(1.560503) q[3];
sx q[3];
rz(-1.5655616) q[3];
sx q[3];
rz(0.43983398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2477766) q[2];
sx q[2];
rz(-0.42422366) q[2];
sx q[2];
rz(-3.0730548) q[2];
rz(2.2872772) q[3];
sx q[3];
rz(-1.8770542) q[3];
sx q[3];
rz(-3.1356649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(1.051556) q[0];
sx q[0];
rz(-2.5904901) q[0];
sx q[0];
rz(-2.0447482) q[0];
rz(2.5966273) q[1];
sx q[1];
rz(-0.68383354) q[1];
sx q[1];
rz(-2.9473238) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0939464) q[0];
sx q[0];
rz(-1.1585278) q[0];
sx q[0];
rz(-2.7591989) q[0];
rz(-1.0164073) q[2];
sx q[2];
rz(-1.28456) q[2];
sx q[2];
rz(-0.58974671) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4020445) q[1];
sx q[1];
rz(-0.6303936) q[1];
sx q[1];
rz(2.5044548) q[1];
rz(-0.13432947) q[3];
sx q[3];
rz(-1.9212011) q[3];
sx q[3];
rz(-2.2131395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.37476173) q[2];
sx q[2];
rz(-1.3730048) q[2];
sx q[2];
rz(0.72804803) q[2];
rz(-2.9211365) q[3];
sx q[3];
rz(-2.9900592) q[3];
sx q[3];
rz(2.5098586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24882889) q[0];
sx q[0];
rz(-1.6680102) q[0];
sx q[0];
rz(-2.072075) q[0];
rz(-2.2546774) q[1];
sx q[1];
rz(-2.6353757) q[1];
sx q[1];
rz(1.3217529) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2489172) q[0];
sx q[0];
rz(-1.6098611) q[0];
sx q[0];
rz(-0.43401735) q[0];
rz(-pi) q[1];
rz(0.23897929) q[2];
sx q[2];
rz(-1.5859446) q[2];
sx q[2];
rz(-0.8123443) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.29577434) q[1];
sx q[1];
rz(-2.008518) q[1];
sx q[1];
rz(0.079016642) q[1];
rz(-pi) q[2];
rz(-2.1683919) q[3];
sx q[3];
rz(-2.0102484) q[3];
sx q[3];
rz(3.1126461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.47361031) q[2];
sx q[2];
rz(-1.0429635) q[2];
sx q[2];
rz(-1.2874416) q[2];
rz(-2.1227664) q[3];
sx q[3];
rz(-2.5835218) q[3];
sx q[3];
rz(0.67600018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92893112) q[0];
sx q[0];
rz(-0.23290817) q[0];
sx q[0];
rz(2.2610597) q[0];
rz(-3.0048043) q[1];
sx q[1];
rz(-2.9158178) q[1];
sx q[1];
rz(0.61019623) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.033025893) q[0];
sx q[0];
rz(-1.5582425) q[0];
sx q[0];
rz(-0.92856284) q[0];
rz(-pi) q[1];
x q[1];
rz(0.43164189) q[2];
sx q[2];
rz(-1.4315617) q[2];
sx q[2];
rz(2.8406258) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.62214855) q[1];
sx q[1];
rz(-0.3560027) q[1];
sx q[1];
rz(0.58575969) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7022252) q[3];
sx q[3];
rz(-1.4476848) q[3];
sx q[3];
rz(-1.1151552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.60787624) q[2];
sx q[2];
rz(-0.72489649) q[2];
sx q[2];
rz(1.2302715) q[2];
rz(2.2599334) q[3];
sx q[3];
rz(-1.6917546) q[3];
sx q[3];
rz(1.8264044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6864258) q[0];
sx q[0];
rz(-1.5979586) q[0];
sx q[0];
rz(2.0879188) q[0];
rz(1.3764489) q[1];
sx q[1];
rz(-2.0339298) q[1];
sx q[1];
rz(0.68127662) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.412144) q[0];
sx q[0];
rz(-0.870847) q[0];
sx q[0];
rz(-1.5922209) q[0];
rz(2.7679484) q[2];
sx q[2];
rz(-2.3609997) q[2];
sx q[2];
rz(-3.0958423) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.14398814) q[1];
sx q[1];
rz(-2.0625084) q[1];
sx q[1];
rz(2.8365447) q[1];
rz(-pi) q[2];
rz(2.2627955) q[3];
sx q[3];
rz(-1.5419772) q[3];
sx q[3];
rz(-1.4895213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.072558746) q[2];
sx q[2];
rz(-2.7636038) q[2];
sx q[2];
rz(-2.6663713) q[2];
rz(1.2230988) q[3];
sx q[3];
rz(-0.94448543) q[3];
sx q[3];
rz(-0.16330115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2326736) q[0];
sx q[0];
rz(-2.1894046) q[0];
sx q[0];
rz(-0.21937823) q[0];
rz(2.1351922) q[1];
sx q[1];
rz(-1.4234797) q[1];
sx q[1];
rz(1.9361608) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69898283) q[0];
sx q[0];
rz(-0.97215334) q[0];
sx q[0];
rz(-3.1151616) q[0];
rz(-pi) q[1];
rz(1.8379301) q[2];
sx q[2];
rz(-1.1877737) q[2];
sx q[2];
rz(2.2522425) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0675065) q[1];
sx q[1];
rz(-1.5536904) q[1];
sx q[1];
rz(-0.19931227) q[1];
rz(-pi) q[2];
rz(-3.0053776) q[3];
sx q[3];
rz(-2.0117617) q[3];
sx q[3];
rz(-1.7441526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.22126234) q[2];
sx q[2];
rz(-0.39322501) q[2];
sx q[2];
rz(-0.66366759) q[2];
rz(2.4237733) q[3];
sx q[3];
rz(-0.65687537) q[3];
sx q[3];
rz(-2.9653911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32888907) q[0];
sx q[0];
rz(-2.1825574) q[0];
sx q[0];
rz(-1.9004199) q[0];
rz(-0.42789704) q[1];
sx q[1];
rz(-1.5792081) q[1];
sx q[1];
rz(1.4235628) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5964029) q[0];
sx q[0];
rz(-1.3529142) q[0];
sx q[0];
rz(-1.6473739) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.15870416) q[2];
sx q[2];
rz(-0.85904143) q[2];
sx q[2];
rz(1.9090599) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2124205) q[1];
sx q[1];
rz(-2.4732686) q[1];
sx q[1];
rz(1.6051488) q[1];
rz(-pi) q[2];
rz(-1.6070494) q[3];
sx q[3];
rz(-1.9055618) q[3];
sx q[3];
rz(-0.70543766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9972035) q[2];
sx q[2];
rz(-2.453697) q[2];
sx q[2];
rz(2.8456861) q[2];
rz(1.0812673) q[3];
sx q[3];
rz(-1.6176977) q[3];
sx q[3];
rz(-0.14779873) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5645771) q[0];
sx q[0];
rz(-2.5595589) q[0];
sx q[0];
rz(-2.4356595) q[0];
rz(2.9544746) q[1];
sx q[1];
rz(-2.3046604) q[1];
sx q[1];
rz(0.44609889) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4574544) q[0];
sx q[0];
rz(-1.7230526) q[0];
sx q[0];
rz(1.5726552) q[0];
rz(-pi) q[1];
rz(-0.31766625) q[2];
sx q[2];
rz(-2.0216172) q[2];
sx q[2];
rz(2.7433155) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0072807) q[1];
sx q[1];
rz(-2.4680302) q[1];
sx q[1];
rz(-2.1271655) q[1];
x q[2];
rz(1.3996077) q[3];
sx q[3];
rz(-2.1406731) q[3];
sx q[3];
rz(0.2967473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.17911653) q[2];
sx q[2];
rz(-2.1160782) q[2];
sx q[2];
rz(1.4013438) q[2];
rz(0.60351795) q[3];
sx q[3];
rz(-2.9602435) q[3];
sx q[3];
rz(-3.0028499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.03183455) q[0];
sx q[0];
rz(-1.4048445) q[0];
sx q[0];
rz(3.0314714) q[0];
rz(-1.3212063) q[1];
sx q[1];
rz(-0.92480129) q[1];
sx q[1];
rz(-2.9168911) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6455655) q[0];
sx q[0];
rz(-1.9651881) q[0];
sx q[0];
rz(2.140122) q[0];
rz(-0.16188936) q[2];
sx q[2];
rz(-1.8542552) q[2];
sx q[2];
rz(3.0155011) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4367392) q[1];
sx q[1];
rz(-0.74932832) q[1];
sx q[1];
rz(2.8981461) q[1];
rz(1.7431926) q[3];
sx q[3];
rz(-1.5972923) q[3];
sx q[3];
rz(0.66988984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3118185) q[2];
sx q[2];
rz(-0.31121397) q[2];
sx q[2];
rz(-2.6149926) q[2];
rz(2.7569568) q[3];
sx q[3];
rz(-1.2472109) q[3];
sx q[3];
rz(-0.23469901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93152355) q[0];
sx q[0];
rz(-0.43545224) q[0];
sx q[0];
rz(-0.42238105) q[0];
rz(2.3474563) q[1];
sx q[1];
rz(-1.7105449) q[1];
sx q[1];
rz(1.7078043) q[1];
rz(2.2591545) q[2];
sx q[2];
rz(-2.0954162) q[2];
sx q[2];
rz(1.4388234) q[2];
rz(1.2023114) q[3];
sx q[3];
rz(-1.9574584) q[3];
sx q[3];
rz(-3.0493469) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
