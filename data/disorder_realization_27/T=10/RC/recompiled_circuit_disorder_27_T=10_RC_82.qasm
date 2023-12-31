OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.1448016) q[0];
sx q[0];
rz(0.15455833) q[0];
sx q[0];
rz(6.9757087) q[0];
rz(1.9321631) q[1];
sx q[1];
rz(-1.2485319) q[1];
sx q[1];
rz(-1.385153) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4614852) q[0];
sx q[0];
rz(-0.84851096) q[0];
sx q[0];
rz(-2.0018342) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.295747) q[2];
sx q[2];
rz(-2.1056471) q[2];
sx q[2];
rz(1.5799074) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7117118) q[1];
sx q[1];
rz(-1.8036588) q[1];
sx q[1];
rz(1.4038605) q[1];
rz(2.4888943) q[3];
sx q[3];
rz(-1.9737118) q[3];
sx q[3];
rz(-2.9626915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8866855) q[2];
sx q[2];
rz(-2.343785) q[2];
sx q[2];
rz(2.936426) q[2];
rz(2.3702879) q[3];
sx q[3];
rz(-2.3588534) q[3];
sx q[3];
rz(-1.1024968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(2.7339864) q[0];
sx q[0];
rz(-2.3953231) q[0];
sx q[0];
rz(0.45390391) q[0];
rz(2.1167963) q[1];
sx q[1];
rz(-0.4075993) q[1];
sx q[1];
rz(1.227238) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5204029) q[0];
sx q[0];
rz(-1.6463917) q[0];
sx q[0];
rz(-1.7002506) q[0];
rz(-pi) q[1];
rz(-1.111755) q[2];
sx q[2];
rz(-1.0978205) q[2];
sx q[2];
rz(-0.74726653) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.36453516) q[1];
sx q[1];
rz(-0.39452663) q[1];
sx q[1];
rz(-2.415034) q[1];
x q[2];
rz(0.20083986) q[3];
sx q[3];
rz(-1.7240932) q[3];
sx q[3];
rz(-0.64985819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.041302117) q[2];
sx q[2];
rz(-1.9561448) q[2];
sx q[2];
rz(-2.5741637) q[2];
rz(0.36519095) q[3];
sx q[3];
rz(-1.7285715) q[3];
sx q[3];
rz(-2.173483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.658618) q[0];
sx q[0];
rz(-2.5768319) q[0];
sx q[0];
rz(-2.2429402) q[0];
rz(2.1458416) q[1];
sx q[1];
rz(-1.5581222) q[1];
sx q[1];
rz(2.8083037) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8329187) q[0];
sx q[0];
rz(-1.1090288) q[0];
sx q[0];
rz(-0.69899107) q[0];
x q[1];
rz(3.0436685) q[2];
sx q[2];
rz(-1.3946748) q[2];
sx q[2];
rz(-2.2224561) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8771613) q[1];
sx q[1];
rz(-2.1794771) q[1];
sx q[1];
rz(2.9560637) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1786147) q[3];
sx q[3];
rz(-0.64219785) q[3];
sx q[3];
rz(1.9574788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4553392) q[2];
sx q[2];
rz(-1.3505961) q[2];
sx q[2];
rz(-1.9906445) q[2];
rz(-0.84093705) q[3];
sx q[3];
rz(-1.1154113) q[3];
sx q[3];
rz(-1.9870728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6999917) q[0];
sx q[0];
rz(-1.561152) q[0];
sx q[0];
rz(2.4568795) q[0];
rz(-1.0355863) q[1];
sx q[1];
rz(-0.5077478) q[1];
sx q[1];
rz(1.9365786) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9158463) q[0];
sx q[0];
rz(-1.6041184) q[0];
sx q[0];
rz(-0.87945625) q[0];
x q[1];
rz(-2.9377405) q[2];
sx q[2];
rz(-0.9365558) q[2];
sx q[2];
rz(-0.39433345) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4971784) q[1];
sx q[1];
rz(-0.89343151) q[1];
sx q[1];
rz(2.7888984) q[1];
x q[2];
rz(1.6550001) q[3];
sx q[3];
rz(-2.4349672) q[3];
sx q[3];
rz(-0.018761793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8923607) q[2];
sx q[2];
rz(-1.426733) q[2];
sx q[2];
rz(-2.7704346) q[2];
rz(-1.4012198) q[3];
sx q[3];
rz(-2.4818082) q[3];
sx q[3];
rz(-1.1192809) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.086833) q[0];
sx q[0];
rz(-0.78614569) q[0];
sx q[0];
rz(-3.0084685) q[0];
rz(-2.1482824) q[1];
sx q[1];
rz(-1.7555833) q[1];
sx q[1];
rz(0.55508074) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13984891) q[0];
sx q[0];
rz(-1.2554597) q[0];
sx q[0];
rz(-3.128201) q[0];
x q[1];
rz(3.073408) q[2];
sx q[2];
rz(-1.1567332) q[2];
sx q[2];
rz(-0.33472862) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.811378) q[1];
sx q[1];
rz(-1.1256309) q[1];
sx q[1];
rz(-0.021275714) q[1];
rz(0.080901905) q[3];
sx q[3];
rz(-1.6228075) q[3];
sx q[3];
rz(-0.73976433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.83539) q[2];
sx q[2];
rz(-1.0058879) q[2];
sx q[2];
rz(0.13892697) q[2];
rz(2.1991918) q[3];
sx q[3];
rz(-1.5706294) q[3];
sx q[3];
rz(-2.5901103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54365629) q[0];
sx q[0];
rz(-0.56977001) q[0];
sx q[0];
rz(-0.57975769) q[0];
rz(3.014091) q[1];
sx q[1];
rz(-1.9516877) q[1];
sx q[1];
rz(-1.6019843) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67028763) q[0];
sx q[0];
rz(-1.7738713) q[0];
sx q[0];
rz(0.85104403) q[0];
x q[1];
rz(-0.13055735) q[2];
sx q[2];
rz(-1.5255543) q[2];
sx q[2];
rz(-3.0322078) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0582038) q[1];
sx q[1];
rz(-1.997588) q[1];
sx q[1];
rz(-2.8298488) q[1];
rz(-0.29069889) q[3];
sx q[3];
rz(-2.2346367) q[3];
sx q[3];
rz(1.3568527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5876028) q[2];
sx q[2];
rz(-0.25045276) q[2];
sx q[2];
rz(-2.8721151) q[2];
rz(0.23412165) q[3];
sx q[3];
rz(-2.6168489) q[3];
sx q[3];
rz(0.060119303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6948029) q[0];
sx q[0];
rz(-0.61820522) q[0];
sx q[0];
rz(-0.68429464) q[0];
rz(-3.0220095) q[1];
sx q[1];
rz(-1.2922829) q[1];
sx q[1];
rz(0.51876846) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2001901) q[0];
sx q[0];
rz(-1.4369643) q[0];
sx q[0];
rz(2.3076513) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3709626) q[2];
sx q[2];
rz(-1.4118886) q[2];
sx q[2];
rz(-0.7030013) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.59473945) q[1];
sx q[1];
rz(-0.93187983) q[1];
sx q[1];
rz(-2.1584595) q[1];
rz(-pi) q[2];
rz(1.5090844) q[3];
sx q[3];
rz(-2.7402174) q[3];
sx q[3];
rz(2.9369831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2542904) q[2];
sx q[2];
rz(-1.3672978) q[2];
sx q[2];
rz(-0.56345144) q[2];
rz(-0.051579483) q[3];
sx q[3];
rz(-2.0023465) q[3];
sx q[3];
rz(0.95190597) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1812487) q[0];
sx q[0];
rz(-2.612817) q[0];
sx q[0];
rz(-1.7425591) q[0];
rz(-2.3545806) q[1];
sx q[1];
rz(-1.1287289) q[1];
sx q[1];
rz(0.74434892) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23119152) q[0];
sx q[0];
rz(-1.4230799) q[0];
sx q[0];
rz(-1.5157248) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3049576) q[2];
sx q[2];
rz(-0.53005866) q[2];
sx q[2];
rz(0.30345464) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5163755) q[1];
sx q[1];
rz(-2.5853734) q[1];
sx q[1];
rz(0.43362995) q[1];
rz(1.5943424) q[3];
sx q[3];
rz(-1.1163201) q[3];
sx q[3];
rz(0.2621032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4259592) q[2];
sx q[2];
rz(-1.7461494) q[2];
sx q[2];
rz(-0.60950935) q[2];
rz(0.65731796) q[3];
sx q[3];
rz(-0.64353639) q[3];
sx q[3];
rz(-0.26143423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9534) q[0];
sx q[0];
rz(-3.0472026) q[0];
sx q[0];
rz(1.5040065) q[0];
rz(1.2414744) q[1];
sx q[1];
rz(-1.1499317) q[1];
sx q[1];
rz(0.77493587) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79561728) q[0];
sx q[0];
rz(-2.2240763) q[0];
sx q[0];
rz(2.6443291) q[0];
rz(-pi) q[1];
rz(2.9853285) q[2];
sx q[2];
rz(-1.1322349) q[2];
sx q[2];
rz(1.4807448) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7271125) q[1];
sx q[1];
rz(-0.4948805) q[1];
sx q[1];
rz(-2.7369569) q[1];
rz(-1.7899412) q[3];
sx q[3];
rz(-1.7100167) q[3];
sx q[3];
rz(-0.62300357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9541786) q[2];
sx q[2];
rz(-0.22970197) q[2];
sx q[2];
rz(2.9837218) q[2];
rz(-1.9291417) q[3];
sx q[3];
rz(-1.0586497) q[3];
sx q[3];
rz(-1.8635748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0697486) q[0];
sx q[0];
rz(-2.1691515) q[0];
sx q[0];
rz(-2.9272595) q[0];
rz(0.65746039) q[1];
sx q[1];
rz(-2.9174556) q[1];
sx q[1];
rz(1.0459895) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25046529) q[0];
sx q[0];
rz(-1.9292826) q[0];
sx q[0];
rz(-1.4573218) q[0];
rz(-pi) q[1];
rz(-0.043141456) q[2];
sx q[2];
rz(-2.5501745) q[2];
sx q[2];
rz(0.45515781) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6007081) q[1];
sx q[1];
rz(-2.4661459) q[1];
sx q[1];
rz(-2.1727968) q[1];
x q[2];
rz(1.5973741) q[3];
sx q[3];
rz(-1.1655032) q[3];
sx q[3];
rz(1.6534896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7252698) q[2];
sx q[2];
rz(-3.0815093) q[2];
sx q[2];
rz(1.0894758) q[2];
rz(-1.5754835) q[3];
sx q[3];
rz(-1.8982866) q[3];
sx q[3];
rz(2.4889448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8626704) q[0];
sx q[0];
rz(-0.60386064) q[0];
sx q[0];
rz(0.84558564) q[0];
rz(1.6090341) q[1];
sx q[1];
rz(-1.4668203) q[1];
sx q[1];
rz(-1.1046881) q[1];
rz(-1.2621461) q[2];
sx q[2];
rz(-1.9625004) q[2];
sx q[2];
rz(0.99448268) q[2];
rz(-3.0388721) q[3];
sx q[3];
rz(-0.552388) q[3];
sx q[3];
rz(-1.296476) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
