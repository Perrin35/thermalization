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
rz(0.61252874) q[0];
sx q[0];
rz(-2.5084485) q[0];
sx q[0];
rz(0.36344114) q[0];
rz(-0.58703047) q[1];
sx q[1];
rz(3.670985) q[1];
sx q[1];
rz(10.932473) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80565174) q[0];
sx q[0];
rz(-1.664289) q[0];
sx q[0];
rz(3.0911618) q[0];
rz(2.5560651) q[2];
sx q[2];
rz(-2.3351933) q[2];
sx q[2];
rz(1.6438649) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9132293) q[1];
sx q[1];
rz(-2.1268561) q[1];
sx q[1];
rz(-0.25198752) q[1];
x q[2];
rz(-2.1430339) q[3];
sx q[3];
rz(-0.84484085) q[3];
sx q[3];
rz(2.69119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.868136) q[2];
sx q[2];
rz(-0.09082219) q[2];
sx q[2];
rz(2.691213) q[2];
rz(3.0716589) q[3];
sx q[3];
rz(-2.0262597) q[3];
sx q[3];
rz(0.44225881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.062334199) q[0];
sx q[0];
rz(-0.82759768) q[0];
sx q[0];
rz(-2.8424679) q[0];
rz(-0.66080719) q[1];
sx q[1];
rz(-2.0313163) q[1];
sx q[1];
rz(2.2181236) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7851334) q[0];
sx q[0];
rz(-1.3970427) q[0];
sx q[0];
rz(1.2453402) q[0];
rz(-pi) q[1];
rz(-2.3908565) q[2];
sx q[2];
rz(-0.3171277) q[2];
sx q[2];
rz(0.37130585) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8750117) q[1];
sx q[1];
rz(-1.204334) q[1];
sx q[1];
rz(-1.5024356) q[1];
x q[2];
rz(2.1577182) q[3];
sx q[3];
rz(-1.8070584) q[3];
sx q[3];
rz(2.7191702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6580761) q[2];
sx q[2];
rz(-0.4008171) q[2];
sx q[2];
rz(0.083902396) q[2];
rz(-2.1462601) q[3];
sx q[3];
rz(-1.9845125) q[3];
sx q[3];
rz(-1.6307433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29626572) q[0];
sx q[0];
rz(-1.9060598) q[0];
sx q[0];
rz(-0.96257019) q[0];
rz(-2.8341017) q[1];
sx q[1];
rz(-2.2566785) q[1];
sx q[1];
rz(-2.9226551) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5318659) q[0];
sx q[0];
rz(-1.5012739) q[0];
sx q[0];
rz(2.7951434) q[0];
rz(-pi) q[1];
rz(-0.23425183) q[2];
sx q[2];
rz(-0.90595379) q[2];
sx q[2];
rz(-1.6981036) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.57812095) q[1];
sx q[1];
rz(-1.5679661) q[1];
sx q[1];
rz(1.1712537) q[1];
rz(-pi) q[2];
rz(0.89618613) q[3];
sx q[3];
rz(-2.0749669) q[3];
sx q[3];
rz(-1.4880796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4585877) q[2];
sx q[2];
rz(-1.747749) q[2];
sx q[2];
rz(1.5031507) q[2];
rz(2.891244) q[3];
sx q[3];
rz(-2.4177347) q[3];
sx q[3];
rz(-2.943128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64971626) q[0];
sx q[0];
rz(-2.8174077) q[0];
sx q[0];
rz(0.31679308) q[0];
rz(-0.016238796) q[1];
sx q[1];
rz(-0.79252487) q[1];
sx q[1];
rz(2.8932755) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1127437) q[0];
sx q[0];
rz(-2.1857977) q[0];
sx q[0];
rz(-0.93407913) q[0];
rz(-0.14713471) q[2];
sx q[2];
rz(-2.6700041) q[2];
sx q[2];
rz(1.4855282) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5491396) q[1];
sx q[1];
rz(-0.35129476) q[1];
sx q[1];
rz(-0.37254803) q[1];
rz(-pi) q[2];
rz(1.789757) q[3];
sx q[3];
rz(-1.3830502) q[3];
sx q[3];
rz(-0.93175542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2636286) q[2];
sx q[2];
rz(-0.10136494) q[2];
sx q[2];
rz(-1.144009) q[2];
rz(2.3794203) q[3];
sx q[3];
rz(-2.6165369) q[3];
sx q[3];
rz(-0.99153668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2852729) q[0];
sx q[0];
rz(-3.0506436) q[0];
sx q[0];
rz(2.0891304) q[0];
rz(-0.4902803) q[1];
sx q[1];
rz(-1.0888381) q[1];
sx q[1];
rz(-3.0163684) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76716048) q[0];
sx q[0];
rz(-1.1781791) q[0];
sx q[0];
rz(-0.36628337) q[0];
rz(-1.8502153) q[2];
sx q[2];
rz(-2.1012857) q[2];
sx q[2];
rz(0.16044469) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.4983212) q[1];
sx q[1];
rz(-1.3210728) q[1];
sx q[1];
rz(3.0672464) q[1];
x q[2];
rz(1.8242314) q[3];
sx q[3];
rz(-0.7281239) q[3];
sx q[3];
rz(0.65092939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.48309717) q[2];
sx q[2];
rz(-2.6039637) q[2];
sx q[2];
rz(-1.5527234) q[2];
rz(3.0850278) q[3];
sx q[3];
rz(-1.5141234) q[3];
sx q[3];
rz(0.27730832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5618943) q[0];
sx q[0];
rz(-2.6772006) q[0];
sx q[0];
rz(-2.7138746) q[0];
rz(-0.5411717) q[1];
sx q[1];
rz(-0.41053694) q[1];
sx q[1];
rz(-1.3637095) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3078474) q[0];
sx q[0];
rz(-2.5302093) q[0];
sx q[0];
rz(1.3788542) q[0];
x q[1];
rz(1.5132163) q[2];
sx q[2];
rz(-1.1442746) q[2];
sx q[2];
rz(0.7698537) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.8783166) q[1];
sx q[1];
rz(-1.8748302) q[1];
sx q[1];
rz(-0.85460198) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1306222) q[3];
sx q[3];
rz(-1.1580165) q[3];
sx q[3];
rz(0.00020325155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.042818729) q[2];
sx q[2];
rz(-2.6218178) q[2];
sx q[2];
rz(2.7322312) q[2];
rz(2.9218819) q[3];
sx q[3];
rz(-1.5493871) q[3];
sx q[3];
rz(1.6519914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4396502) q[0];
sx q[0];
rz(-0.42124978) q[0];
sx q[0];
rz(-0.14081328) q[0];
rz(0.5367865) q[1];
sx q[1];
rz(-1.355492) q[1];
sx q[1];
rz(2.7649194) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0842782) q[0];
sx q[0];
rz(-1.8207826) q[0];
sx q[0];
rz(2.9765073) q[0];
x q[1];
rz(-0.48483168) q[2];
sx q[2];
rz(-0.73532447) q[2];
sx q[2];
rz(1.5905274) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6766549) q[1];
sx q[1];
rz(-0.79494093) q[1];
sx q[1];
rz(1.388768) q[1];
x q[2];
rz(0.41531947) q[3];
sx q[3];
rz(-1.4276081) q[3];
sx q[3];
rz(-0.48986942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0385527) q[2];
sx q[2];
rz(-1.1662177) q[2];
sx q[2];
rz(-0.80123365) q[2];
rz(1.0771105) q[3];
sx q[3];
rz(-0.25682768) q[3];
sx q[3];
rz(-0.1424772) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.899026) q[0];
sx q[0];
rz(-2.9809256) q[0];
sx q[0];
rz(2.2005079) q[0];
rz(-2.2451027) q[1];
sx q[1];
rz(-1.4184378) q[1];
sx q[1];
rz(2.3268708) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77234948) q[0];
sx q[0];
rz(-2.6008478) q[0];
sx q[0];
rz(-0.7665828) q[0];
rz(-pi) q[1];
rz(-2.8553634) q[2];
sx q[2];
rz(-0.17795086) q[2];
sx q[2];
rz(-2.9432757) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5960434) q[1];
sx q[1];
rz(-1.909797) q[1];
sx q[1];
rz(0.91722644) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0031189) q[3];
sx q[3];
rz(-1.8879297) q[3];
sx q[3];
rz(0.88256146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5251289) q[2];
sx q[2];
rz(-2.5069671) q[2];
sx q[2];
rz(2.5647707) q[2];
rz(-0.47354627) q[3];
sx q[3];
rz(-1.1487995) q[3];
sx q[3];
rz(-0.043341652) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43805495) q[0];
sx q[0];
rz(-2.3860768) q[0];
sx q[0];
rz(-0.087906539) q[0];
rz(1.5535304) q[1];
sx q[1];
rz(-0.5843662) q[1];
sx q[1];
rz(2.8839674) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9842897) q[0];
sx q[0];
rz(-2.7098795) q[0];
sx q[0];
rz(1.8583511) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.94341523) q[2];
sx q[2];
rz(-1.5134654) q[2];
sx q[2];
rz(1.0939404) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.39752326) q[1];
sx q[1];
rz(-1.9998614) q[1];
sx q[1];
rz(1.3123162) q[1];
rz(0.78225953) q[3];
sx q[3];
rz(-2.3965947) q[3];
sx q[3];
rz(-1.2040602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0660925) q[2];
sx q[2];
rz(-1.4614212) q[2];
sx q[2];
rz(-3.063805) q[2];
rz(2.0870239) q[3];
sx q[3];
rz(-2.4013897) q[3];
sx q[3];
rz(2.4906702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8973549) q[0];
sx q[0];
rz(-3.1081508) q[0];
sx q[0];
rz(-2.7058386) q[0];
rz(-0.032595366) q[1];
sx q[1];
rz(-2.4041912) q[1];
sx q[1];
rz(1.7936515) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12440269) q[0];
sx q[0];
rz(-1.5778871) q[0];
sx q[0];
rz(0.75049863) q[0];
x q[1];
rz(0.52705898) q[2];
sx q[2];
rz(-1.7816451) q[2];
sx q[2];
rz(-1.859889) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.15534066) q[1];
sx q[1];
rz(-1.5568451) q[1];
sx q[1];
rz(1.6068684) q[1];
rz(-pi) q[2];
rz(-0.3147821) q[3];
sx q[3];
rz(-1.7932995) q[3];
sx q[3];
rz(0.088929847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9255479) q[2];
sx q[2];
rz(-1.0341158) q[2];
sx q[2];
rz(0.75630581) q[2];
rz(-3.0594399) q[3];
sx q[3];
rz(-0.3923471) q[3];
sx q[3];
rz(2.4503585) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3638851) q[0];
sx q[0];
rz(-2.4019699) q[0];
sx q[0];
rz(-0.85516047) q[0];
rz(1.1126407) q[1];
sx q[1];
rz(-1.1323977) q[1];
sx q[1];
rz(-0.36066537) q[1];
rz(-2.9654824) q[2];
sx q[2];
rz(-1.817504) q[2];
sx q[2];
rz(0.70155406) q[2];
rz(-2.1300456) q[3];
sx q[3];
rz(-1.9725298) q[3];
sx q[3];
rz(-1.6172668) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
