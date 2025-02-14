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
rz(-2.6900238) q[0];
sx q[0];
rz(-1.4181925) q[0];
sx q[0];
rz(-0.42887846) q[0];
rz(-3.2847326) q[1];
sx q[1];
rz(-1.0909456) q[1];
sx q[1];
rz(9.3446891) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0319666) q[0];
sx q[0];
rz(-1.4916149) q[0];
sx q[0];
rz(1.5840785) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8378788) q[2];
sx q[2];
rz(-0.3729698) q[2];
sx q[2];
rz(0.8276865) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.061191843) q[1];
sx q[1];
rz(-0.75729232) q[1];
sx q[1];
rz(-1.9593092) q[1];
rz(-pi) q[2];
x q[2];
rz(2.892602) q[3];
sx q[3];
rz(-0.59196977) q[3];
sx q[3];
rz(-1.0764866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.074778883) q[2];
sx q[2];
rz(-2.0902233) q[2];
sx q[2];
rz(1.7195513) q[2];
rz(-1.413013) q[3];
sx q[3];
rz(-1.6001817) q[3];
sx q[3];
rz(-1.6933256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1220575) q[0];
sx q[0];
rz(-1.3966565) q[0];
sx q[0];
rz(-0.33357093) q[0];
rz(0.16593274) q[1];
sx q[1];
rz(-2.797778) q[1];
sx q[1];
rz(-2.3894892) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.89276) q[0];
sx q[0];
rz(-1.1247096) q[0];
sx q[0];
rz(2.2463138) q[0];
rz(-pi) q[1];
rz(1.5444438) q[2];
sx q[2];
rz(-1.7894723) q[2];
sx q[2];
rz(1.5062576) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5673654) q[1];
sx q[1];
rz(-0.53814473) q[1];
sx q[1];
rz(-1.3567011) q[1];
rz(2.8798298) q[3];
sx q[3];
rz(-1.3532551) q[3];
sx q[3];
rz(1.0378154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.077945262) q[2];
sx q[2];
rz(-2.2293978) q[2];
sx q[2];
rz(-2.4845541) q[2];
rz(-0.71197236) q[3];
sx q[3];
rz(-0.65198055) q[3];
sx q[3];
rz(0.74487346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4350568) q[0];
sx q[0];
rz(-0.61921316) q[0];
sx q[0];
rz(2.1001429) q[0];
rz(-2.7667747) q[1];
sx q[1];
rz(-1.638214) q[1];
sx q[1];
rz(-0.55694881) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81129148) q[0];
sx q[0];
rz(-0.82152237) q[0];
sx q[0];
rz(-2.6514451) q[0];
x q[1];
rz(-3.0663112) q[2];
sx q[2];
rz(-1.1202552) q[2];
sx q[2];
rz(-2.0496429) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8136602) q[1];
sx q[1];
rz(-2.2371462) q[1];
sx q[1];
rz(-1.440669) q[1];
rz(-pi) q[2];
rz(-1.8076423) q[3];
sx q[3];
rz(-0.93333731) q[3];
sx q[3];
rz(1.6528296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.53691429) q[2];
sx q[2];
rz(-2.2173209) q[2];
sx q[2];
rz(2.7064145) q[2];
rz(2.6168881) q[3];
sx q[3];
rz(-0.33354959) q[3];
sx q[3];
rz(3.0068126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10114577) q[0];
sx q[0];
rz(-2.0772159) q[0];
sx q[0];
rz(1.5144298) q[0];
rz(1.0487522) q[1];
sx q[1];
rz(-2.2188413) q[1];
sx q[1];
rz(-1.4357766) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.993282) q[0];
sx q[0];
rz(-1.3956426) q[0];
sx q[0];
rz(-1.5272642) q[0];
rz(-pi) q[1];
rz(2.5975512) q[2];
sx q[2];
rz(-0.61464192) q[2];
sx q[2];
rz(2.7934472) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.857261) q[1];
sx q[1];
rz(-1.5067717) q[1];
sx q[1];
rz(2.0928736) q[1];
x q[2];
rz(1.1069) q[3];
sx q[3];
rz(-1.2588862) q[3];
sx q[3];
rz(2.8297765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.40020308) q[2];
sx q[2];
rz(-1.3578537) q[2];
sx q[2];
rz(3.1164361) q[2];
rz(-1.6545506) q[3];
sx q[3];
rz(-2.7436723) q[3];
sx q[3];
rz(-0.81620836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0461034) q[0];
sx q[0];
rz(-0.64952055) q[0];
sx q[0];
rz(-0.15897861) q[0];
rz(-1.1051296) q[1];
sx q[1];
rz(-0.72560328) q[1];
sx q[1];
rz(0.22430688) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8858356) q[0];
sx q[0];
rz(-2.1191264) q[0];
sx q[0];
rz(-1.095039) q[0];
rz(-pi) q[1];
rz(2.131046) q[2];
sx q[2];
rz(-1.59861) q[2];
sx q[2];
rz(1.5459488) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7182279) q[1];
sx q[1];
rz(-2.0929687) q[1];
sx q[1];
rz(-3.0084064) q[1];
rz(-pi) q[2];
rz(-1.4442611) q[3];
sx q[3];
rz(-0.99735051) q[3];
sx q[3];
rz(-0.8265178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9362681) q[2];
sx q[2];
rz(-1.779413) q[2];
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
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9801789) q[0];
sx q[0];
rz(-2.7730589) q[0];
sx q[0];
rz(-0.73032105) q[0];
rz(-1.6397363) q[1];
sx q[1];
rz(-1.390099) q[1];
sx q[1];
rz(-2.9072442) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6788819) q[0];
sx q[0];
rz(-2.2214799) q[0];
sx q[0];
rz(-1.0048871) q[0];
x q[1];
rz(2.5426834) q[2];
sx q[2];
rz(-2.1258168) q[2];
sx q[2];
rz(2.2428494) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.40962666) q[1];
sx q[1];
rz(-1.5232558) q[1];
sx q[1];
rz(2.9788245) q[1];
rz(3.0803063) q[3];
sx q[3];
rz(-1.7109814) q[3];
sx q[3];
rz(2.9621027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.19330567) q[2];
sx q[2];
rz(-2.327658) q[2];
sx q[2];
rz(1.2813655) q[2];
rz(2.5365601) q[3];
sx q[3];
rz(-1.0993967) q[3];
sx q[3];
rz(1.9090778) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.989711) q[0];
sx q[0];
rz(-1.2622156) q[0];
sx q[0];
rz(2.5285517) q[0];
rz(2.4609861) q[1];
sx q[1];
rz(-2.0304408) q[1];
sx q[1];
rz(1.8162762) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2277057) q[0];
sx q[0];
rz(-1.3038346) q[0];
sx q[0];
rz(2.5771229) q[0];
rz(-1.4875796) q[2];
sx q[2];
rz(-1.3583379) q[2];
sx q[2];
rz(2.0709289) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.22205849) q[1];
sx q[1];
rz(-0.57812837) q[1];
sx q[1];
rz(2.9085552) q[1];
x q[2];
rz(2.5595253) q[3];
sx q[3];
rz(-1.2305233) q[3];
sx q[3];
rz(2.4916864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1943872) q[2];
sx q[2];
rz(-1.9401865) q[2];
sx q[2];
rz(0.80805937) q[2];
rz(-0.64588532) q[3];
sx q[3];
rz(-0.26791993) q[3];
sx q[3];
rz(0.3033692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
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
rz(1.5658257) q[0];
sx q[0];
rz(-0.72661) q[0];
sx q[0];
rz(0.44276825) q[0];
rz(-0.59204656) q[1];
sx q[1];
rz(-2.4202085) q[1];
sx q[1];
rz(2.9659042) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93568703) q[0];
sx q[0];
rz(-1.924951) q[0];
sx q[0];
rz(-1.8649578) q[0];
x q[1];
rz(-1.6106583) q[2];
sx q[2];
rz(-2.022036) q[2];
sx q[2];
rz(-0.892837) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.95388874) q[1];
sx q[1];
rz(-0.64262455) q[1];
sx q[1];
rz(-3.0695555) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6457186) q[3];
sx q[3];
rz(-2.1225833) q[3];
sx q[3];
rz(1.4648579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.59777933) q[2];
sx q[2];
rz(-0.74395776) q[2];
sx q[2];
rz(2.234484) q[2];
rz(2.234327) q[3];
sx q[3];
rz(-1.3943358) q[3];
sx q[3];
rz(-3.0239014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-1.1842781) q[0];
sx q[0];
rz(-0.90183455) q[0];
sx q[0];
rz(-0.19004518) q[0];
rz(1.5589145) q[1];
sx q[1];
rz(-1.5946486) q[1];
sx q[1];
rz(1.83439) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.060576749) q[0];
sx q[0];
rz(-1.2071484) q[0];
sx q[0];
rz(2.9878997) q[0];
rz(1.1747423) q[2];
sx q[2];
rz(-2.4017334) q[2];
sx q[2];
rz(-2.6137874) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9865454) q[1];
sx q[1];
rz(-1.6035756) q[1];
sx q[1];
rz(-1.4630586) q[1];
x q[2];
rz(-1.4553189) q[3];
sx q[3];
rz(-1.1388123) q[3];
sx q[3];
rz(0.40528709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.068453161) q[2];
sx q[2];
rz(-1.0871202) q[2];
sx q[2];
rz(0.20115176) q[2];
rz(1.0226095) q[3];
sx q[3];
rz(-0.68251959) q[3];
sx q[3];
rz(-1.539591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96711838) q[0];
sx q[0];
rz(-1.0459463) q[0];
sx q[0];
rz(2.2651267) q[0];
rz(-0.62249741) q[1];
sx q[1];
rz(-0.21177706) q[1];
sx q[1];
rz(0.34271398) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4349608) q[0];
sx q[0];
rz(-1.6366658) q[0];
sx q[0];
rz(-0.94191282) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9079952) q[2];
sx q[2];
rz(-0.35524455) q[2];
sx q[2];
rz(0.97081414) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4058447) q[1];
sx q[1];
rz(-2.1830507) q[1];
sx q[1];
rz(2.1275747) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7539135) q[3];
sx q[3];
rz(-1.3362543) q[3];
sx q[3];
rz(-2.925433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9524625) q[2];
sx q[2];
rz(-1.6996982) q[2];
sx q[2];
rz(-0.44378898) q[2];
rz(2.022838) q[3];
sx q[3];
rz(-1.5464562) q[3];
sx q[3];
rz(-0.55383033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2832058) q[0];
sx q[0];
rz(-1.5571742) q[0];
sx q[0];
rz(-1.5207186) q[0];
rz(3.0781147) q[1];
sx q[1];
rz(-1.7964446) q[1];
sx q[1];
rz(-1.7367015) q[1];
rz(3.079698) q[2];
sx q[2];
rz(-1.2409004) q[2];
sx q[2];
rz(0.093766669) q[2];
rz(1.9431498) q[3];
sx q[3];
rz(-0.59881864) q[3];
sx q[3];
rz(-1.3074085) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
