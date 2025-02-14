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
rz(3.7747369) q[0];
sx q[0];
rz(9.7882191) q[0];
rz(-0.58703047) q[1];
sx q[1];
rz(3.670985) q[1];
sx q[1];
rz(10.932473) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3717361) q[0];
sx q[0];
rz(-1.6210067) q[0];
sx q[0];
rz(-1.4771853) q[0];
rz(-pi) q[1];
rz(-2.5560651) q[2];
sx q[2];
rz(-2.3351933) q[2];
sx q[2];
rz(-1.6438649) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.22836338) q[1];
sx q[1];
rz(-2.1268561) q[1];
sx q[1];
rz(-2.8896051) q[1];
rz(-0.54777439) q[3];
sx q[3];
rz(-2.2507077) q[3];
sx q[3];
rz(0.31991968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2734566) q[2];
sx q[2];
rz(-0.09082219) q[2];
sx q[2];
rz(-0.4503797) q[2];
rz(3.0716589) q[3];
sx q[3];
rz(-2.0262597) q[3];
sx q[3];
rz(-2.6993338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.062334199) q[0];
sx q[0];
rz(-0.82759768) q[0];
sx q[0];
rz(-2.8424679) q[0];
rz(2.4807855) q[1];
sx q[1];
rz(-2.0313163) q[1];
sx q[1];
rz(-0.92346907) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4006996) q[0];
sx q[0];
rz(-2.7741196) q[0];
sx q[0];
rz(-1.0687554) q[0];
rz(-1.3505351) q[2];
sx q[2];
rz(-1.3407605) q[2];
sx q[2];
rz(-2.736614) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8128494) q[1];
sx q[1];
rz(-1.6346115) q[1];
sx q[1];
rz(2.7743474) q[1];
rz(-pi) q[2];
rz(-0.98387444) q[3];
sx q[3];
rz(-1.8070584) q[3];
sx q[3];
rz(2.7191702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6580761) q[2];
sx q[2];
rz(-0.4008171) q[2];
sx q[2];
rz(-3.0576903) q[2];
rz(0.99533254) q[3];
sx q[3];
rz(-1.1570802) q[3];
sx q[3];
rz(1.6307433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29626572) q[0];
sx q[0];
rz(-1.9060598) q[0];
sx q[0];
rz(0.96257019) q[0];
rz(-0.30749097) q[1];
sx q[1];
rz(-2.2566785) q[1];
sx q[1];
rz(2.9226551) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22902952) q[0];
sx q[0];
rz(-2.7885127) q[0];
sx q[0];
rz(-2.9393239) q[0];
x q[1];
rz(2.9073408) q[2];
sx q[2];
rz(-2.2356389) q[2];
sx q[2];
rz(1.6981036) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5634717) q[1];
sx q[1];
rz(-1.5736266) q[1];
sx q[1];
rz(-1.970339) q[1];
rz(-pi) q[2];
rz(2.2454065) q[3];
sx q[3];
rz(-1.0666258) q[3];
sx q[3];
rz(-1.4880796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4585877) q[2];
sx q[2];
rz(-1.747749) q[2];
sx q[2];
rz(-1.5031507) q[2];
rz(2.891244) q[3];
sx q[3];
rz(-2.4177347) q[3];
sx q[3];
rz(-2.943128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4918764) q[0];
sx q[0];
rz(-2.8174077) q[0];
sx q[0];
rz(2.8247996) q[0];
rz(-0.016238796) q[1];
sx q[1];
rz(-2.3490678) q[1];
sx q[1];
rz(-2.8932755) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.028848948) q[0];
sx q[0];
rz(-0.95579493) q[0];
sx q[0];
rz(-2.2075135) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6454208) q[2];
sx q[2];
rz(-2.0368825) q[2];
sx q[2];
rz(1.8209195) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5491396) q[1];
sx q[1];
rz(-0.35129476) q[1];
sx q[1];
rz(-2.7690446) q[1];
rz(-0.85217441) q[3];
sx q[3];
rz(-2.8541454) q[3];
sx q[3];
rz(-1.3369657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2636286) q[2];
sx q[2];
rz(-3.0402277) q[2];
sx q[2];
rz(1.144009) q[2];
rz(0.7621724) q[3];
sx q[3];
rz(-0.52505571) q[3];
sx q[3];
rz(2.150056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85631973) q[0];
sx q[0];
rz(-0.090949051) q[0];
sx q[0];
rz(-1.0524622) q[0];
rz(-2.6513124) q[1];
sx q[1];
rz(-2.0527546) q[1];
sx q[1];
rz(0.12522423) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1220684) q[0];
sx q[0];
rz(-0.53036371) q[0];
sx q[0];
rz(2.2838462) q[0];
rz(-0.54791656) q[2];
sx q[2];
rz(-1.8109908) q[2];
sx q[2];
rz(-1.5870767) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.4983212) q[1];
sx q[1];
rz(-1.8205199) q[1];
sx q[1];
rz(3.0672464) q[1];
rz(2.282827) q[3];
sx q[3];
rz(-1.4031583) q[3];
sx q[3];
rz(2.0307547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6584955) q[2];
sx q[2];
rz(-0.53762895) q[2];
sx q[2];
rz(-1.5888692) q[2];
rz(3.0850278) q[3];
sx q[3];
rz(-1.6274692) q[3];
sx q[3];
rz(-0.27730832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5618943) q[0];
sx q[0];
rz(-2.6772006) q[0];
sx q[0];
rz(-0.42771801) q[0];
rz(-0.5411717) q[1];
sx q[1];
rz(-0.41053694) q[1];
sx q[1];
rz(-1.3637095) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60073131) q[0];
sx q[0];
rz(-0.97222881) q[0];
sx q[0];
rz(-3.0086584) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.42714675) q[2];
sx q[2];
rz(-1.5183798) q[2];
sx q[2];
rz(-2.3644931) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.8783166) q[1];
sx q[1];
rz(-1.8748302) q[1];
sx q[1];
rz(0.85460198) q[1];
x q[2];
rz(-1.545752) q[3];
sx q[3];
rz(-0.4129172) q[3];
sx q[3];
rz(-3.1144547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0987739) q[2];
sx q[2];
rz(-2.6218178) q[2];
sx q[2];
rz(2.7322312) q[2];
rz(0.21971075) q[3];
sx q[3];
rz(-1.5493871) q[3];
sx q[3];
rz(1.4896013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4396502) q[0];
sx q[0];
rz(-2.7203429) q[0];
sx q[0];
rz(0.14081328) q[0];
rz(0.5367865) q[1];
sx q[1];
rz(-1.355492) q[1];
sx q[1];
rz(-0.37667325) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46464065) q[0];
sx q[0];
rz(-0.29862216) q[0];
sx q[0];
rz(-0.99891164) q[0];
rz(-1.9697625) q[2];
sx q[2];
rz(-2.2062183) q[2];
sx q[2];
rz(-0.93346006) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.72195327) q[1];
sx q[1];
rz(-0.7925539) q[1];
sx q[1];
rz(-0.18246095) q[1];
rz(0.34318681) q[3];
sx q[3];
rz(-2.7036441) q[3];
sx q[3];
rz(1.7476976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0385527) q[2];
sx q[2];
rz(-1.1662177) q[2];
sx q[2];
rz(2.340359) q[2];
rz(2.0644821) q[3];
sx q[3];
rz(-2.884765) q[3];
sx q[3];
rz(-0.1424772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.899026) q[0];
sx q[0];
rz(-2.9809256) q[0];
sx q[0];
rz(-0.94108474) q[0];
rz(2.2451027) q[1];
sx q[1];
rz(-1.4184378) q[1];
sx q[1];
rz(-2.3268708) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3692432) q[0];
sx q[0];
rz(-0.54074484) q[0];
sx q[0];
rz(2.3750099) q[0];
rz(-pi) q[1];
rz(-2.9707387) q[2];
sx q[2];
rz(-1.5207982) q[2];
sx q[2];
rz(1.4871666) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.22413177) q[1];
sx q[1];
rz(-0.96012174) q[1];
sx q[1];
rz(0.4179722) q[1];
rz(-pi) q[2];
rz(-0.37121935) q[3];
sx q[3];
rz(-2.1069848) q[3];
sx q[3];
rz(-0.49193383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5251289) q[2];
sx q[2];
rz(-0.63462555) q[2];
sx q[2];
rz(-0.57682192) q[2];
rz(-0.47354627) q[3];
sx q[3];
rz(-1.9927931) q[3];
sx q[3];
rz(-3.098251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7035377) q[0];
sx q[0];
rz(-2.3860768) q[0];
sx q[0];
rz(0.087906539) q[0];
rz(1.5880623) q[1];
sx q[1];
rz(-0.5843662) q[1];
sx q[1];
rz(-2.8839674) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6759153) q[0];
sx q[0];
rz(-1.6897461) q[0];
sx q[0];
rz(1.9867935) q[0];
x q[1];
rz(-2.1981774) q[2];
sx q[2];
rz(-1.5134654) q[2];
sx q[2];
rz(2.0476523) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.39752326) q[1];
sx q[1];
rz(-1.1417312) q[1];
sx q[1];
rz(1.8292765) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5622665) q[3];
sx q[3];
rz(-2.0690479) q[3];
sx q[3];
rz(-2.1439596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.075500175) q[2];
sx q[2];
rz(-1.6801715) q[2];
sx q[2];
rz(-0.077787682) q[2];
rz(2.0870239) q[3];
sx q[3];
rz(-0.74020296) q[3];
sx q[3];
rz(-2.4906702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
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
rz(1.2442378) q[0];
sx q[0];
rz(-0.033441823) q[0];
sx q[0];
rz(-0.43575409) q[0];
rz(0.032595366) q[1];
sx q[1];
rz(-0.73740149) q[1];
sx q[1];
rz(1.7936515) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6885869) q[0];
sx q[0];
rz(-0.82032114) q[0];
sx q[0];
rz(1.561101) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6145337) q[2];
sx q[2];
rz(-1.3599476) q[2];
sx q[2];
rz(-1.859889) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7256335) q[1];
sx q[1];
rz(-1.5347278) q[1];
sx q[1];
rz(-0.013960268) q[1];
rz(-0.63107154) q[3];
sx q[3];
rz(-2.7582599) q[3];
sx q[3];
rz(2.255343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9255479) q[2];
sx q[2];
rz(-2.1074769) q[2];
sx q[2];
rz(2.3852868) q[2];
rz(0.082152724) q[3];
sx q[3];
rz(-0.3923471) q[3];
sx q[3];
rz(2.4503585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77770752) q[0];
sx q[0];
rz(-2.4019699) q[0];
sx q[0];
rz(-0.85516047) q[0];
rz(-1.1126407) q[1];
sx q[1];
rz(-2.009195) q[1];
sx q[1];
rz(2.7809273) q[1];
rz(2.9654824) q[2];
sx q[2];
rz(-1.3240887) q[2];
sx q[2];
rz(-2.4400386) q[2];
rz(2.2459946) q[3];
sx q[3];
rz(-0.67586312) q[3];
sx q[3];
rz(-2.6296658) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
