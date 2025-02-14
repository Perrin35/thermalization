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
rz(0.51377327) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.936718) q[0];
sx q[0];
rz(-0.74130171) q[0];
sx q[0];
rz(2.3573849) q[0];
rz(0.47494048) q[2];
sx q[2];
rz(-1.6261592) q[2];
sx q[2];
rz(2.9621015) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.39756718) q[1];
sx q[1];
rz(-1.356522) q[1];
sx q[1];
rz(0.32886966) q[1];
x q[2];
rz(1.6052395) q[3];
sx q[3];
rz(-2.1625966) q[3];
sx q[3];
rz(-0.47815048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.674268) q[2];
sx q[2];
rz(-1.5200204) q[2];
sx q[2];
rz(0.33217126) q[2];
rz(-2.8915571) q[3];
sx q[3];
rz(-1.2269521) q[3];
sx q[3];
rz(-1.8488319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8914723) q[0];
sx q[0];
rz(-1.253506) q[0];
sx q[0];
rz(-2.6943595) q[0];
rz(-1.0294634) q[1];
sx q[1];
rz(-0.58164683) q[1];
sx q[1];
rz(-0.10118016) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48866901) q[0];
sx q[0];
rz(-2.5258668) q[0];
sx q[0];
rz(2.4020345) q[0];
x q[1];
rz(-2.9415628) q[2];
sx q[2];
rz(-2.1406271) q[2];
sx q[2];
rz(2.9404158) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9969719) q[1];
sx q[1];
rz(-0.50077866) q[1];
sx q[1];
rz(2.5624496) q[1];
rz(-pi) q[2];
rz(0.11240837) q[3];
sx q[3];
rz(-2.2066457) q[3];
sx q[3];
rz(-2.4684255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0589361) q[2];
sx q[2];
rz(-2.3857748) q[2];
sx q[2];
rz(-1.6015046) q[2];
rz(-2.5236409) q[3];
sx q[3];
rz(-0.58303419) q[3];
sx q[3];
rz(1.9790953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59442941) q[0];
sx q[0];
rz(-1.0048486) q[0];
sx q[0];
rz(2.6255703) q[0];
rz(1.0524606) q[1];
sx q[1];
rz(-0.92603374) q[1];
sx q[1];
rz(1.1722391) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2207551) q[0];
sx q[0];
rz(-1.1467548) q[0];
sx q[0];
rz(1.5829257) q[0];
x q[1];
rz(-1.486023) q[2];
sx q[2];
rz(-0.71226487) q[2];
sx q[2];
rz(-2.149947) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7683923) q[1];
sx q[1];
rz(-0.52881634) q[1];
sx q[1];
rz(1.1952728) q[1];
x q[2];
rz(-0.013962176) q[3];
sx q[3];
rz(-0.80028557) q[3];
sx q[3];
rz(1.518702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8664794) q[2];
sx q[2];
rz(-1.2629513) q[2];
sx q[2];
rz(-0.67683721) q[2];
rz(2.6721241) q[3];
sx q[3];
rz(-0.89478409) q[3];
sx q[3];
rz(-1.9278056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4645828) q[0];
sx q[0];
rz(-1.0053758) q[0];
sx q[0];
rz(1.9418035) q[0];
rz(-2.5489573) q[1];
sx q[1];
rz(-2.0694144) q[1];
sx q[1];
rz(1.4592272) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71818411) q[0];
sx q[0];
rz(-0.99331021) q[0];
sx q[0];
rz(1.3363488) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8921669) q[2];
sx q[2];
rz(-0.81455961) q[2];
sx q[2];
rz(-0.047175353) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.693599) q[1];
sx q[1];
rz(-1.3136615) q[1];
sx q[1];
rz(-2.2092381) q[1];
rz(-1.9497037) q[3];
sx q[3];
rz(-2.2188596) q[3];
sx q[3];
rz(0.60175446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6478641) q[2];
sx q[2];
rz(-1.0961327) q[2];
sx q[2];
rz(-2.4596821) q[2];
rz(-2.72825) q[3];
sx q[3];
rz(-1.6091434) q[3];
sx q[3];
rz(1.1558862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-0.84363627) q[0];
sx q[0];
rz(-1.6652668) q[0];
sx q[0];
rz(0.27808878) q[0];
rz(2.2372712) q[1];
sx q[1];
rz(-1.4537289) q[1];
sx q[1];
rz(-2.1542737) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0678735) q[0];
sx q[0];
rz(-0.6530531) q[0];
sx q[0];
rz(-1.1730173) q[0];
rz(-2.9402307) q[2];
sx q[2];
rz(-3.0018531) q[2];
sx q[2];
rz(1.0248549) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1121171) q[1];
sx q[1];
rz(-1.6320845) q[1];
sx q[1];
rz(-1.1809096) q[1];
x q[2];
rz(-2.868949) q[3];
sx q[3];
rz(-1.050625) q[3];
sx q[3];
rz(2.4369654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1272993) q[2];
sx q[2];
rz(-2.1743446) q[2];
sx q[2];
rz(-0.75931749) q[2];
rz(-1.6589818) q[3];
sx q[3];
rz(-1.9448152) q[3];
sx q[3];
rz(-2.7827941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21869126) q[0];
sx q[0];
rz(-0.8194812) q[0];
sx q[0];
rz(2.4500093) q[0];
rz(0.85451952) q[1];
sx q[1];
rz(-0.71047345) q[1];
sx q[1];
rz(2.7202594) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3499419) q[0];
sx q[0];
rz(-1.6393108) q[0];
sx q[0];
rz(1.1380026) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3106903) q[2];
sx q[2];
rz(-2.75616) q[2];
sx q[2];
rz(1.3192434) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9142575) q[1];
sx q[1];
rz(-2.7707639) q[1];
sx q[1];
rz(2.9079958) q[1];
x q[2];
rz(-0.60507994) q[3];
sx q[3];
rz(-1.7972094) q[3];
sx q[3];
rz(-0.070158557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.93115807) q[2];
sx q[2];
rz(-1.6633818) q[2];
sx q[2];
rz(-0.20273905) q[2];
rz(1.5984795) q[3];
sx q[3];
rz(-0.32035443) q[3];
sx q[3];
rz(-1.4690442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(0.28432524) q[0];
sx q[0];
rz(-1.4284644) q[0];
sx q[0];
rz(-0.39886928) q[0];
rz(1.162989) q[1];
sx q[1];
rz(-0.82610026) q[1];
sx q[1];
rz(2.861048) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6847072) q[0];
sx q[0];
rz(-2.5734757) q[0];
sx q[0];
rz(0.76596188) q[0];
rz(-pi) q[1];
rz(-2.0899827) q[2];
sx q[2];
rz(-0.74724846) q[2];
sx q[2];
rz(2.0488536) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.39574341) q[1];
sx q[1];
rz(-1.9956931) q[1];
sx q[1];
rz(-2.8060071) q[1];
x q[2];
rz(-0.97784247) q[3];
sx q[3];
rz(-0.8675608) q[3];
sx q[3];
rz(-3.0271526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3205388) q[2];
sx q[2];
rz(-2.1619449) q[2];
sx q[2];
rz(1.4873571) q[2];
rz(2.1982543) q[3];
sx q[3];
rz(-0.92883674) q[3];
sx q[3];
rz(1.9802035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5281552) q[0];
sx q[0];
rz(-1.4035839) q[0];
sx q[0];
rz(2.6392537) q[0];
rz(2.595064) q[1];
sx q[1];
rz(-1.2556475) q[1];
sx q[1];
rz(0.46195236) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1621409) q[0];
sx q[0];
rz(-2.5573423) q[0];
sx q[0];
rz(1.0229848) q[0];
rz(0.91945474) q[2];
sx q[2];
rz(-1.5327138) q[2];
sx q[2];
rz(3.0236261) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2425804) q[1];
sx q[1];
rz(-1.1625152) q[1];
sx q[1];
rz(-3.0228843) q[1];
rz(-pi) q[2];
rz(-0.53047385) q[3];
sx q[3];
rz(-1.5598205) q[3];
sx q[3];
rz(-0.61179841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.97303566) q[2];
sx q[2];
rz(-2.904197) q[2];
sx q[2];
rz(2.9711704) q[2];
rz(2.6863344) q[3];
sx q[3];
rz(-1.5452496) q[3];
sx q[3];
rz(2.4236603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(0.51668733) q[1];
sx q[1];
rz(-1.5645212) q[1];
sx q[1];
rz(2.7076941) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11182298) q[0];
sx q[0];
rz(-1.0445458) q[0];
sx q[0];
rz(0.49825251) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0185616) q[2];
sx q[2];
rz(-1.2553213) q[2];
sx q[2];
rz(1.4184784) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.89131195) q[1];
sx q[1];
rz(-1.6906594) q[1];
sx q[1];
rz(2.2054172) q[1];
rz(-2.2614165) q[3];
sx q[3];
rz(-1.7978906) q[3];
sx q[3];
rz(-2.2736062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.58069289) q[2];
sx q[2];
rz(-0.43792024) q[2];
sx q[2];
rz(-1.8221347) q[2];
rz(-0.27160078) q[3];
sx q[3];
rz(-1.8235794) q[3];
sx q[3];
rz(-2.0297089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.6962947) q[0];
sx q[0];
rz(-2.2309208) q[0];
sx q[0];
rz(1.7355504) q[0];
rz(1.4259074) q[1];
sx q[1];
rz(-1.1049263) q[1];
sx q[1];
rz(-1.2474733) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46215993) q[0];
sx q[0];
rz(-1.3094256) q[0];
sx q[0];
rz(2.0716459) q[0];
x q[1];
rz(-0.054385238) q[2];
sx q[2];
rz(-1.453389) q[2];
sx q[2];
rz(2.1942127) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0746324) q[1];
sx q[1];
rz(-1.6116643) q[1];
sx q[1];
rz(2.7994179) q[1];
x q[2];
rz(2.5188451) q[3];
sx q[3];
rz(-2.1295665) q[3];
sx q[3];
rz(-1.2947242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4282816) q[2];
sx q[2];
rz(-2.9094628) q[2];
sx q[2];
rz(-0.10078079) q[2];
rz(3.1320324) q[3];
sx q[3];
rz(-1.5848426) q[3];
sx q[3];
rz(1.6077707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99209256) q[0];
sx q[0];
rz(-1.7763573) q[0];
sx q[0];
rz(-1.1471164) q[0];
rz(-0.59303444) q[1];
sx q[1];
rz(-1.4321764) q[1];
sx q[1];
rz(2.166116) q[1];
rz(1.1723761) q[2];
sx q[2];
rz(-1.6360184) q[2];
sx q[2];
rz(1.2498275) q[2];
rz(1.5531874) q[3];
sx q[3];
rz(-1.7918158) q[3];
sx q[3];
rz(0.50504897) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
