OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.71894574) q[0];
sx q[0];
rz(-2.5526241) q[0];
sx q[0];
rz(-0.11864057) q[0];
rz(2.150382) q[1];
sx q[1];
rz(-1.041643) q[1];
sx q[1];
rz(-2.6278194) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2048747) q[0];
sx q[0];
rz(-2.4002909) q[0];
sx q[0];
rz(-2.3573849) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5085601) q[2];
sx q[2];
rz(-2.0449491) q[2];
sx q[2];
rz(-1.7218423) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2456677) q[1];
sx q[1];
rz(-1.8918719) q[1];
sx q[1];
rz(1.3447869) q[1];
rz(-pi) q[2];
rz(-0.59207497) q[3];
sx q[3];
rz(-1.5993803) q[3];
sx q[3];
rz(1.1118654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.674268) q[2];
sx q[2];
rz(-1.5200204) q[2];
sx q[2];
rz(-0.33217126) q[2];
rz(0.25003555) q[3];
sx q[3];
rz(-1.2269521) q[3];
sx q[3];
rz(1.2927607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8914723) q[0];
sx q[0];
rz(-1.253506) q[0];
sx q[0];
rz(2.6943595) q[0];
rz(-2.1121292) q[1];
sx q[1];
rz(-0.58164683) q[1];
sx q[1];
rz(0.10118016) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44199884) q[0];
sx q[0];
rz(-1.9706107) q[0];
sx q[0];
rz(2.6599822) q[0];
rz(-pi) q[1];
rz(-2.9415628) q[2];
sx q[2];
rz(-1.0009655) q[2];
sx q[2];
rz(0.20117682) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.90535962) q[1];
sx q[1];
rz(-1.3049076) q[1];
sx q[1];
rz(-2.7120525) q[1];
rz(-1.7216136) q[3];
sx q[3];
rz(-0.64435092) q[3];
sx q[3];
rz(2.2805813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0589361) q[2];
sx q[2];
rz(-0.75581789) q[2];
sx q[2];
rz(-1.5400881) q[2];
rz(-0.61795175) q[3];
sx q[3];
rz(-0.58303419) q[3];
sx q[3];
rz(-1.9790953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5471632) q[0];
sx q[0];
rz(-2.1367441) q[0];
sx q[0];
rz(-0.51602236) q[0];
rz(1.0524606) q[1];
sx q[1];
rz(-2.2155589) q[1];
sx q[1];
rz(1.9693536) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2207551) q[0];
sx q[0];
rz(-1.9948378) q[0];
sx q[0];
rz(1.5829257) q[0];
rz(-3.0686106) q[2];
sx q[2];
rz(-2.2799645) q[2];
sx q[2];
rz(-2.2617509) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0111085) q[1];
sx q[1];
rz(-1.3846894) q[1];
sx q[1];
rz(-2.0687201) q[1];
rz(-3.1276305) q[3];
sx q[3];
rz(-0.80028557) q[3];
sx q[3];
rz(1.6228907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.27511328) q[2];
sx q[2];
rz(-1.2629513) q[2];
sx q[2];
rz(2.4647554) q[2];
rz(0.46946851) q[3];
sx q[3];
rz(-2.2468086) q[3];
sx q[3];
rz(1.2137871) q[3];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6770099) q[0];
sx q[0];
rz(-2.1362169) q[0];
sx q[0];
rz(-1.1997892) q[0];
rz(0.59263539) q[1];
sx q[1];
rz(-1.0721782) q[1];
sx q[1];
rz(1.6823654) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71818411) q[0];
sx q[0];
rz(-2.1482824) q[0];
sx q[0];
rz(-1.8052438) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8267385) q[2];
sx q[2];
rz(-0.78849604) q[2];
sx q[2];
rz(-2.7389604) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.063733405) q[1];
sx q[1];
rz(-2.1850536) q[1];
sx q[1];
rz(-0.31645223) q[1];
rz(-pi) q[2];
rz(-2.4577971) q[3];
sx q[3];
rz(-1.2714362) q[3];
sx q[3];
rz(-1.9366858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.49372855) q[2];
sx q[2];
rz(-1.0961327) q[2];
sx q[2];
rz(0.68191051) q[2];
rz(2.72825) q[3];
sx q[3];
rz(-1.5324493) q[3];
sx q[3];
rz(1.1558862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2979564) q[0];
sx q[0];
rz(-1.4763259) q[0];
sx q[0];
rz(-2.8635039) q[0];
rz(0.9043215) q[1];
sx q[1];
rz(-1.6878637) q[1];
sx q[1];
rz(-2.1542737) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.322418) q[0];
sx q[0];
rz(-1.3331945) q[0];
sx q[0];
rz(2.1850719) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9402307) q[2];
sx q[2];
rz(-3.0018531) q[2];
sx q[2];
rz(-1.0248549) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.029475538) q[1];
sx q[1];
rz(-1.5095081) q[1];
sx q[1];
rz(1.1809096) q[1];
x q[2];
rz(-2.1073527) q[3];
sx q[3];
rz(-1.8066386) q[3];
sx q[3];
rz(-0.72808108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.014293369) q[2];
sx q[2];
rz(-2.1743446) q[2];
sx q[2];
rz(-0.75931749) q[2];
rz(1.6589818) q[3];
sx q[3];
rz(-1.9448152) q[3];
sx q[3];
rz(-0.35879859) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9229014) q[0];
sx q[0];
rz(-2.3221115) q[0];
sx q[0];
rz(-0.6915834) q[0];
rz(0.85451952) q[1];
sx q[1];
rz(-2.4311192) q[1];
sx q[1];
rz(0.42133322) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74752676) q[0];
sx q[0];
rz(-1.1390862) q[0];
sx q[0];
rz(-3.0661446) q[0];
rz(-pi) q[1];
rz(2.8745804) q[2];
sx q[2];
rz(-1.2894372) q[2];
sx q[2];
rz(-2.6002778) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5616963) q[1];
sx q[1];
rz(-1.6547799) q[1];
sx q[1];
rz(2.7799699) q[1];
x q[2];
rz(1.297703) q[3];
sx q[3];
rz(-2.1583302) q[3];
sx q[3];
rz(-1.3466101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2104346) q[2];
sx q[2];
rz(-1.6633818) q[2];
sx q[2];
rz(0.20273905) q[2];
rz(-1.5984795) q[3];
sx q[3];
rz(-2.8212382) q[3];
sx q[3];
rz(-1.4690442) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8572674) q[0];
sx q[0];
rz(-1.7131282) q[0];
sx q[0];
rz(-0.39886928) q[0];
rz(1.9786037) q[1];
sx q[1];
rz(-0.82610026) q[1];
sx q[1];
rz(-2.861048) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6847072) q[0];
sx q[0];
rz(-2.5734757) q[0];
sx q[0];
rz(-0.76596188) q[0];
rz(-pi) q[1];
rz(-2.248204) q[2];
sx q[2];
rz(-1.2268434) q[2];
sx q[2];
rz(-2.266573) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.39574341) q[1];
sx q[1];
rz(-1.1458995) q[1];
sx q[1];
rz(0.33558553) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.79645282) q[3];
sx q[3];
rz(-2.011125) q[3];
sx q[3];
rz(-1.274282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3205388) q[2];
sx q[2];
rz(-2.1619449) q[2];
sx q[2];
rz(-1.4873571) q[2];
rz(-2.1982543) q[3];
sx q[3];
rz(-0.92883674) q[3];
sx q[3];
rz(1.1613891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5281552) q[0];
sx q[0];
rz(-1.4035839) q[0];
sx q[0];
rz(2.6392537) q[0];
rz(0.54652864) q[1];
sx q[1];
rz(-1.2556475) q[1];
sx q[1];
rz(2.6796403) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7936458) q[0];
sx q[0];
rz(-2.0610621) q[0];
sx q[0];
rz(-0.33167524) q[0];
x q[1];
rz(-1.6335604) q[2];
sx q[2];
rz(-0.6522921) q[2];
sx q[2];
rz(1.4029274) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2425804) q[1];
sx q[1];
rz(-1.1625152) q[1];
sx q[1];
rz(-3.0228843) q[1];
rz(-pi) q[2];
rz(0.53047385) q[3];
sx q[3];
rz(-1.5598205) q[3];
sx q[3];
rz(-2.5297942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.168557) q[2];
sx q[2];
rz(-2.904197) q[2];
sx q[2];
rz(-2.9711704) q[2];
rz(0.45525822) q[3];
sx q[3];
rz(-1.5452496) q[3];
sx q[3];
rz(0.71793238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97487226) q[0];
sx q[0];
rz(-1.5019187) q[0];
sx q[0];
rz(2.9611452) q[0];
rz(-0.51668733) q[1];
sx q[1];
rz(-1.5645212) q[1];
sx q[1];
rz(0.43389854) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0297697) q[0];
sx q[0];
rz(-1.0445458) q[0];
sx q[0];
rz(-0.49825251) q[0];
rz(-1.1230311) q[2];
sx q[2];
rz(-1.2553213) q[2];
sx q[2];
rz(-1.4184784) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5499159) q[1];
sx q[1];
rz(-2.2001451) q[1];
sx q[1];
rz(-2.9931328) q[1];
rz(-pi) q[2];
rz(-0.29124864) q[3];
sx q[3];
rz(-0.90121239) q[3];
sx q[3];
rz(-0.88676363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5608998) q[2];
sx q[2];
rz(-2.7036724) q[2];
sx q[2];
rz(-1.8221347) q[2];
rz(-0.27160078) q[3];
sx q[3];
rz(-1.3180132) q[3];
sx q[3];
rz(-1.1118838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.44529799) q[0];
sx q[0];
rz(-0.91067186) q[0];
sx q[0];
rz(-1.4060422) q[0];
rz(1.7156853) q[1];
sx q[1];
rz(-2.0366663) q[1];
sx q[1];
rz(-1.2474733) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1734764) q[0];
sx q[0];
rz(-1.0884459) q[0];
sx q[0];
rz(-0.29598693) q[0];
rz(1.1389473) q[2];
sx q[2];
rz(-0.12933918) q[2];
sx q[2];
rz(0.51233385) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.48161246) q[1];
sx q[1];
rz(-1.9126737) q[1];
sx q[1];
rz(1.6141763) q[1];
rz(2.3214809) q[3];
sx q[3];
rz(-2.3305428) q[3];
sx q[3];
rz(0.91203989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4282816) q[2];
sx q[2];
rz(-2.9094628) q[2];
sx q[2];
rz(-0.10078079) q[2];
rz(3.1320324) q[3];
sx q[3];
rz(-1.5848426) q[3];
sx q[3];
rz(-1.5338219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.99209256) q[0];
sx q[0];
rz(-1.7763573) q[0];
sx q[0];
rz(-1.1471164) q[0];
rz(2.5485582) q[1];
sx q[1];
rz(-1.4321764) q[1];
sx q[1];
rz(2.166116) q[1];
rz(-1.4040074) q[2];
sx q[2];
rz(-0.40344147) q[2];
sx q[2];
rz(-0.16735195) q[2];
rz(-1.5884052) q[3];
sx q[3];
rz(-1.7918158) q[3];
sx q[3];
rz(0.50504897) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
