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
rz(-3.2847326) q[1];
sx q[1];
rz(-1.0909456) q[1];
sx q[1];
rz(9.3446891) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53777909) q[0];
sx q[0];
rz(-1.5575557) q[0];
sx q[0];
rz(0.079188345) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9316788) q[2];
sx q[2];
rz(-1.4744802) q[2];
sx q[2];
rz(-2.148984) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.061191843) q[1];
sx q[1];
rz(-2.3843003) q[1];
sx q[1];
rz(1.1822834) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4065885) q[3];
sx q[3];
rz(-0.99942151) q[3];
sx q[3];
rz(-2.3624275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.074778883) q[2];
sx q[2];
rz(-2.0902233) q[2];
sx q[2];
rz(1.7195513) q[2];
rz(-1.7285796) q[3];
sx q[3];
rz(-1.541411) q[3];
sx q[3];
rz(-1.6933256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(1.0195352) q[0];
sx q[0];
rz(-1.3966565) q[0];
sx q[0];
rz(-2.8080217) q[0];
rz(-0.16593274) q[1];
sx q[1];
rz(-2.797778) q[1];
sx q[1];
rz(2.3894892) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2488326) q[0];
sx q[0];
rz(-2.016883) q[0];
sx q[0];
rz(-2.2463138) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0235748) q[2];
sx q[2];
rz(-2.9213597) q[2];
sx q[2];
rz(-1.756246) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.57422728) q[1];
sx q[1];
rz(-0.53814473) q[1];
sx q[1];
rz(1.7848915) q[1];
rz(-pi) q[2];
rz(-2.8798298) q[3];
sx q[3];
rz(-1.7883375) q[3];
sx q[3];
rz(1.0378154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0636474) q[2];
sx q[2];
rz(-0.91219488) q[2];
sx q[2];
rz(-0.65703854) q[2];
rz(-0.71197236) q[3];
sx q[3];
rz(-2.4896121) q[3];
sx q[3];
rz(2.3967192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4350568) q[0];
sx q[0];
rz(-2.5223795) q[0];
sx q[0];
rz(2.1001429) q[0];
rz(0.374818) q[1];
sx q[1];
rz(-1.638214) q[1];
sx q[1];
rz(2.5846438) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81129148) q[0];
sx q[0];
rz(-2.3200703) q[0];
sx q[0];
rz(2.6514451) q[0];
rz(-1.416549) q[2];
sx q[2];
rz(-2.6852312) q[2];
sx q[2];
rz(-2.2211494) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8136602) q[1];
sx q[1];
rz(-0.90444649) q[1];
sx q[1];
rz(-1.440669) q[1];
rz(2.8347765) q[3];
sx q[3];
rz(-2.4673438) q[3];
sx q[3];
rz(-1.1034654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6046784) q[2];
sx q[2];
rz(-2.2173209) q[2];
sx q[2];
rz(-0.43517819) q[2];
rz(-0.52470454) q[3];
sx q[3];
rz(-2.8080431) q[3];
sx q[3];
rz(-3.0068126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10114577) q[0];
sx q[0];
rz(-2.0772159) q[0];
sx q[0];
rz(1.6271628) q[0];
rz(-1.0487522) q[1];
sx q[1];
rz(-0.92275134) q[1];
sx q[1];
rz(1.705816) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.993282) q[0];
sx q[0];
rz(-1.7459501) q[0];
sx q[0];
rz(-1.5272642) q[0];
x q[1];
rz(0.54404144) q[2];
sx q[2];
rz(-2.5269507) q[2];
sx q[2];
rz(-0.34814542) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8183349) q[1];
sx q[1];
rz(-1.0498974) q[1];
sx q[1];
rz(-3.0677615) q[1];
rz(-pi) q[2];
x q[2];
rz(0.34603186) q[3];
sx q[3];
rz(-1.1308987) q[3];
sx q[3];
rz(1.4113246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.40020308) q[2];
sx q[2];
rz(-1.3578537) q[2];
sx q[2];
rz(3.1164361) q[2];
rz(1.4870421) q[3];
sx q[3];
rz(-2.7436723) q[3];
sx q[3];
rz(-0.81620836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0461034) q[0];
sx q[0];
rz(-2.4920721) q[0];
sx q[0];
rz(0.15897861) q[0];
rz(-2.036463) q[1];
sx q[1];
rz(-0.72560328) q[1];
sx q[1];
rz(-0.22430688) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.255757) q[0];
sx q[0];
rz(-1.0224662) q[0];
sx q[0];
rz(-2.0465536) q[0];
rz(-1.6231027) q[2];
sx q[2];
rz(-2.5807267) q[2];
sx q[2];
rz(-3.0724409) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6857924) q[1];
sx q[1];
rz(-2.6042301) q[1];
sx q[1];
rz(-1.7975897) q[1];
rz(-0.57711011) q[3];
sx q[3];
rz(-1.464586) q[3];
sx q[3];
rz(0.67536992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.20532456) q[2];
sx q[2];
rz(-1.3621796) q[2];
sx q[2];
rz(0.74860191) q[2];
rz(0.11048206) q[3];
sx q[3];
rz(-1.2275077) q[3];
sx q[3];
rz(-0.41142472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1614138) q[0];
sx q[0];
rz(-2.7730589) q[0];
sx q[0];
rz(-0.73032105) q[0];
rz(-1.5018564) q[1];
sx q[1];
rz(-1.7514936) q[1];
sx q[1];
rz(-2.9072442) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6662009) q[0];
sx q[0];
rz(-2.0115543) q[0];
sx q[0];
rz(-0.73385629) q[0];
rz(2.2147398) q[2];
sx q[2];
rz(-1.0711462) q[2];
sx q[2];
rz(-0.32678963) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6987353) q[1];
sx q[1];
rz(-0.1695098) q[1];
sx q[1];
rz(-0.28556602) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1612915) q[3];
sx q[3];
rz(-2.9886768) q[3];
sx q[3];
rz(-2.5482938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.19330567) q[2];
sx q[2];
rz(-0.81393465) q[2];
sx q[2];
rz(-1.8602271) q[2];
rz(-2.5365601) q[3];
sx q[3];
rz(-2.0421959) q[3];
sx q[3];
rz(1.9090778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.989711) q[0];
sx q[0];
rz(-1.2622156) q[0];
sx q[0];
rz(-0.61304098) q[0];
rz(2.4609861) q[1];
sx q[1];
rz(-1.1111518) q[1];
sx q[1];
rz(-1.8162762) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2277057) q[0];
sx q[0];
rz(-1.3038346) q[0];
sx q[0];
rz(0.56446979) q[0];
x q[1];
rz(-2.9284186) q[2];
sx q[2];
rz(-1.4894549) q[2];
sx q[2];
rz(2.6590462) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.054106709) q[1];
sx q[1];
rz(-2.1313868) q[1];
sx q[1];
rz(-1.4212379) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5595253) q[3];
sx q[3];
rz(-1.9110693) q[3];
sx q[3];
rz(-0.64990625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.94720542) q[2];
sx q[2];
rz(-1.9401865) q[2];
sx q[2];
rz(2.3335333) q[2];
rz(2.4957073) q[3];
sx q[3];
rz(-2.8736727) q[3];
sx q[3];
rz(2.8382235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5658257) q[0];
sx q[0];
rz(-0.72661) q[0];
sx q[0];
rz(-2.6988244) q[0];
rz(2.5495461) q[1];
sx q[1];
rz(-2.4202085) q[1];
sx q[1];
rz(-0.17568849) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21766996) q[0];
sx q[0];
rz(-2.6852073) q[0];
sx q[0];
rz(0.66500591) q[0];
x q[1];
rz(-1.5309344) q[2];
sx q[2];
rz(-2.022036) q[2];
sx q[2];
rz(-2.2487557) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4669816) q[1];
sx q[1];
rz(-1.5276485) q[1];
sx q[1];
rz(-0.64137913) q[1];
rz(-pi) q[2];
rz(-0.91276987) q[3];
sx q[3];
rz(-2.4173749) q[3];
sx q[3];
rz(-2.4778544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.59777933) q[2];
sx q[2];
rz(-0.74395776) q[2];
sx q[2];
rz(-0.90710863) q[2];
rz(-0.90726566) q[3];
sx q[3];
rz(-1.7472569) q[3];
sx q[3];
rz(3.0239014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1842781) q[0];
sx q[0];
rz(-2.2397581) q[0];
sx q[0];
rz(-0.19004518) q[0];
rz(1.5589145) q[1];
sx q[1];
rz(-1.546944) q[1];
sx q[1];
rz(1.3072026) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0810159) q[0];
sx q[0];
rz(-1.9344442) q[0];
sx q[0];
rz(2.9878997) q[0];
rz(-pi) q[1];
rz(-2.2707269) q[2];
sx q[2];
rz(-1.8339089) q[2];
sx q[2];
rz(-1.342529) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9865454) q[1];
sx q[1];
rz(-1.6035756) q[1];
sx q[1];
rz(1.4630586) q[1];
rz(-0.43453026) q[3];
sx q[3];
rz(-1.4659681) q[3];
sx q[3];
rz(-1.2140345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.068453161) q[2];
sx q[2];
rz(-2.0544724) q[2];
sx q[2];
rz(-2.9404409) q[2];
rz(1.0226095) q[3];
sx q[3];
rz(-0.68251959) q[3];
sx q[3];
rz(1.6020017) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96711838) q[0];
sx q[0];
rz(-1.0459463) q[0];
sx q[0];
rz(0.87646595) q[0];
rz(-0.62249741) q[1];
sx q[1];
rz(-0.21177706) q[1];
sx q[1];
rz(0.34271398) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7066318) q[0];
sx q[0];
rz(-1.6366658) q[0];
sx q[0];
rz(-2.1996798) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2335974) q[2];
sx q[2];
rz(-0.35524455) q[2];
sx q[2];
rz(-0.97081414) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6500998) q[1];
sx q[1];
rz(-2.0180185) q[1];
sx q[1];
rz(-2.4504721) q[1];
rz(1.8233772) q[3];
sx q[3];
rz(-1.9473206) q[3];
sx q[3];
rz(-1.8815709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9524625) q[2];
sx q[2];
rz(-1.4418944) q[2];
sx q[2];
rz(2.6978037) q[2];
rz(1.1187547) q[3];
sx q[3];
rz(-1.5951364) q[3];
sx q[3];
rz(2.5877623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8583869) q[0];
sx q[0];
rz(-1.5844185) q[0];
sx q[0];
rz(1.6208741) q[0];
rz(-0.063477909) q[1];
sx q[1];
rz(-1.7964446) q[1];
sx q[1];
rz(-1.7367015) q[1];
rz(-3.079698) q[2];
sx q[2];
rz(-1.9006922) q[2];
sx q[2];
rz(-3.047826) q[2];
rz(-1.0045814) q[3];
sx q[3];
rz(-1.3642642) q[3];
sx q[3];
rz(3.0929079) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
