OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.24703439) q[0];
sx q[0];
rz(4.3972754) q[0];
sx q[0];
rz(9.7527405) q[0];
rz(2.9070931) q[1];
sx q[1];
rz(-0.20107888) q[1];
sx q[1];
rz(0.091436401) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8592005) q[0];
sx q[0];
rz(-0.99635591) q[0];
sx q[0];
rz(-1.0919071) q[0];
rz(-pi) q[1];
rz(-1.7056866) q[2];
sx q[2];
rz(-1.3139259) q[2];
sx q[2];
rz(-0.97193064) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.92802231) q[1];
sx q[1];
rz(-1.3334647) q[1];
sx q[1];
rz(-0.084753239) q[1];
rz(-2.923008) q[3];
sx q[3];
rz(-2.1544475) q[3];
sx q[3];
rz(-2.1368795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4771007) q[2];
sx q[2];
rz(-2.1680919) q[2];
sx q[2];
rz(-1.1260024) q[2];
rz(0.27515718) q[3];
sx q[3];
rz(-2.5313009) q[3];
sx q[3];
rz(2.2385105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0098410957) q[0];
sx q[0];
rz(-0.69350243) q[0];
sx q[0];
rz(-2.4480208) q[0];
rz(1.0961078) q[1];
sx q[1];
rz(-2.1577436) q[1];
sx q[1];
rz(-0.19031659) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099146518) q[0];
sx q[0];
rz(-2.6876246) q[0];
sx q[0];
rz(2.0367665) q[0];
rz(-2.5311567) q[2];
sx q[2];
rz(-2.5154841) q[2];
sx q[2];
rz(2.6420643) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7479334) q[1];
sx q[1];
rz(-0.84888443) q[1];
sx q[1];
rz(3.1118803) q[1];
x q[2];
rz(2.4137647) q[3];
sx q[3];
rz(-1.2065294) q[3];
sx q[3];
rz(1.4621853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6341614) q[2];
sx q[2];
rz(-2.4298411) q[2];
sx q[2];
rz(-1.2197536) q[2];
rz(0.1427342) q[3];
sx q[3];
rz(-2.1295363) q[3];
sx q[3];
rz(2.232961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74293566) q[0];
sx q[0];
rz(-2.0516899) q[0];
sx q[0];
rz(2.3213342) q[0];
rz(0.29207686) q[1];
sx q[1];
rz(-2.0675817) q[1];
sx q[1];
rz(1.8935727) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4688063) q[0];
sx q[0];
rz(-2.5350223) q[0];
sx q[0];
rz(2.9787105) q[0];
rz(-1.5065932) q[2];
sx q[2];
rz(-1.777613) q[2];
sx q[2];
rz(-1.9931672) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6229912) q[1];
sx q[1];
rz(-2.2366183) q[1];
sx q[1];
rz(-1.4856505) q[1];
rz(-1.1850584) q[3];
sx q[3];
rz(-0.68813656) q[3];
sx q[3];
rz(3.0582173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.32039207) q[2];
sx q[2];
rz(-1.0504861) q[2];
sx q[2];
rz(-1.2134264) q[2];
rz(2.976867) q[3];
sx q[3];
rz(-2.2475593) q[3];
sx q[3];
rz(-0.081610672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40959013) q[0];
sx q[0];
rz(-1.859917) q[0];
sx q[0];
rz(0.18606342) q[0];
rz(-0.20446725) q[1];
sx q[1];
rz(-2.6696413) q[1];
sx q[1];
rz(1.8444555) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0766749) q[0];
sx q[0];
rz(-3.0041822) q[0];
sx q[0];
rz(-1.8924367) q[0];
rz(-pi) q[1];
rz(2.9739431) q[2];
sx q[2];
rz(-1.285458) q[2];
sx q[2];
rz(0.97380762) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.90721005) q[1];
sx q[1];
rz(-2.130059) q[1];
sx q[1];
rz(0.84749605) q[1];
rz(-pi) q[2];
rz(-1.4947055) q[3];
sx q[3];
rz(-2.9615059) q[3];
sx q[3];
rz(-2.0538581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.90594784) q[2];
sx q[2];
rz(-2.3489958) q[2];
sx q[2];
rz(-2.2223991) q[2];
rz(-2.8202608) q[3];
sx q[3];
rz(-1.0682169) q[3];
sx q[3];
rz(1.8937768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17386757) q[0];
sx q[0];
rz(-2.6350832) q[0];
sx q[0];
rz(2.2633973) q[0];
rz(-1.8163266) q[1];
sx q[1];
rz(-1.4122496) q[1];
sx q[1];
rz(-1.7153046) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4195376) q[0];
sx q[0];
rz(-2.4800322) q[0];
sx q[0];
rz(1.4369591) q[0];
rz(-1.8833141) q[2];
sx q[2];
rz(-1.3442355) q[2];
sx q[2];
rz(2.8001919) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9702455) q[1];
sx q[1];
rz(-2.4023224) q[1];
sx q[1];
rz(1.7319748) q[1];
rz(-1.250508) q[3];
sx q[3];
rz(-0.8257782) q[3];
sx q[3];
rz(0.68009963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2720126) q[2];
sx q[2];
rz(-11*pi/13) q[2];
sx q[2];
rz(0.041794725) q[2];
rz(-0.061491866) q[3];
sx q[3];
rz(-2.0525335) q[3];
sx q[3];
rz(2.7048236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3549266) q[0];
sx q[0];
rz(-0.085952856) q[0];
sx q[0];
rz(-2.1110995) q[0];
rz(-2.4018535) q[1];
sx q[1];
rz(-1.5286427) q[1];
sx q[1];
rz(2.5700263) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3540928) q[0];
sx q[0];
rz(-1.0618853) q[0];
sx q[0];
rz(2.4718168) q[0];
rz(-0.90494855) q[2];
sx q[2];
rz(-1.4706503) q[2];
sx q[2];
rz(3.1040994) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7476269) q[1];
sx q[1];
rz(-2.6193301) q[1];
sx q[1];
rz(1.0737435) q[1];
rz(-pi) q[2];
rz(0.81482859) q[3];
sx q[3];
rz(-1.0199254) q[3];
sx q[3];
rz(-1.8967472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.17343865) q[2];
sx q[2];
rz(-2.3996694) q[2];
sx q[2];
rz(-2.5773933) q[2];
rz(-0.12600222) q[3];
sx q[3];
rz(-1.6832451) q[3];
sx q[3];
rz(-2.8469767) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0512222) q[0];
sx q[0];
rz(-2.9416961) q[0];
sx q[0];
rz(2.4293161) q[0];
rz(2.6157216) q[1];
sx q[1];
rz(-0.41627517) q[1];
sx q[1];
rz(-0.66551048) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77008477) q[0];
sx q[0];
rz(-1.5718282) q[0];
sx q[0];
rz(-1.8130215) q[0];
rz(-pi) q[1];
x q[1];
rz(2.268157) q[2];
sx q[2];
rz(-1.4566112) q[2];
sx q[2];
rz(0.21542491) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4588786) q[1];
sx q[1];
rz(-0.51358089) q[1];
sx q[1];
rz(1.1125803) q[1];
rz(2.9969278) q[3];
sx q[3];
rz(-1.33107) q[3];
sx q[3];
rz(3.1160115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.15726382) q[2];
sx q[2];
rz(-2.4815346) q[2];
sx q[2];
rz(1.7187913) q[2];
rz(0.11519365) q[3];
sx q[3];
rz(-2.6264103) q[3];
sx q[3];
rz(2.9490024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0916864) q[0];
sx q[0];
rz(-2.005907) q[0];
sx q[0];
rz(2.9507622) q[0];
rz(0.62675369) q[1];
sx q[1];
rz(-2.125506) q[1];
sx q[1];
rz(-2.802882) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8341634) q[0];
sx q[0];
rz(-1.2273664) q[0];
sx q[0];
rz(-0.15983454) q[0];
x q[1];
rz(-2.68967) q[2];
sx q[2];
rz(-1.7799313) q[2];
sx q[2];
rz(2.506633) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.415886) q[1];
sx q[1];
rz(-0.48735122) q[1];
sx q[1];
rz(-1.332167) q[1];
x q[2];
rz(2.8645105) q[3];
sx q[3];
rz(-1.8859366) q[3];
sx q[3];
rz(1.806123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4954341) q[2];
sx q[2];
rz(-0.63445264) q[2];
sx q[2];
rz(2.4411566) q[2];
rz(0.8979848) q[3];
sx q[3];
rz(-1.286641) q[3];
sx q[3];
rz(-0.17351304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.609628) q[0];
sx q[0];
rz(-0.48159972) q[0];
sx q[0];
rz(-2.4832446) q[0];
rz(-0.61093962) q[1];
sx q[1];
rz(-1.8860203) q[1];
sx q[1];
rz(3.0019965) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0567719) q[0];
sx q[0];
rz(-3.0773101) q[0];
sx q[0];
rz(-1.8875185) q[0];
x q[1];
rz(-0.058768674) q[2];
sx q[2];
rz(-0.81548703) q[2];
sx q[2];
rz(-1.8975443) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.82397205) q[1];
sx q[1];
rz(-0.52597731) q[1];
sx q[1];
rz(-0.13336639) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6406052) q[3];
sx q[3];
rz(-0.6456635) q[3];
sx q[3];
rz(2.0696236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5247941) q[2];
sx q[2];
rz(-1.0107661) q[2];
sx q[2];
rz(2.810478) q[2];
rz(2.3838499) q[3];
sx q[3];
rz(-2.7527633) q[3];
sx q[3];
rz(-0.087879114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8452633) q[0];
sx q[0];
rz(-2.1305278) q[0];
sx q[0];
rz(0.18558003) q[0];
rz(-1.0962076) q[1];
sx q[1];
rz(-0.21462333) q[1];
sx q[1];
rz(-1.4846444) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.164924) q[0];
sx q[0];
rz(-0.52545588) q[0];
sx q[0];
rz(2.1303961) q[0];
x q[1];
rz(-1.2946285) q[2];
sx q[2];
rz(-1.3535) q[2];
sx q[2];
rz(-1.5664139) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.85877307) q[1];
sx q[1];
rz(-2.7532137) q[1];
sx q[1];
rz(-1.8358843) q[1];
x q[2];
rz(-1.5649892) q[3];
sx q[3];
rz(-0.70492893) q[3];
sx q[3];
rz(-0.73202902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.93402702) q[2];
sx q[2];
rz(-0.88946122) q[2];
sx q[2];
rz(-0.55220848) q[2];
rz(-0.77783716) q[3];
sx q[3];
rz(-0.85770291) q[3];
sx q[3];
rz(-2.623693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8778397) q[0];
sx q[0];
rz(-1.5120266) q[0];
sx q[0];
rz(1.7396447) q[0];
rz(-2.9539625) q[1];
sx q[1];
rz(-1.7788806) q[1];
sx q[1];
rz(2.3685041) q[1];
rz(0.27992237) q[2];
sx q[2];
rz(-1.4301849) q[2];
sx q[2];
rz(2.4591597) q[2];
rz(-0.88541661) q[3];
sx q[3];
rz(-1.0049184) q[3];
sx q[3];
rz(0.64941209) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
