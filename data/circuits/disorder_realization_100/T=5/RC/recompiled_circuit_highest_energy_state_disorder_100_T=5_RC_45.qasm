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
rz(1.963653) q[0];
sx q[0];
rz(-0.093129245) q[0];
sx q[0];
rz(7.0153305) q[0];
rz(-1.5692476) q[1];
sx q[1];
rz(-0.88580004) q[1];
sx q[1];
rz(-2.7028309) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8092361) q[0];
sx q[0];
rz(-1.2389823) q[0];
sx q[0];
rz(-3.0829264) q[0];
x q[1];
rz(-0.26023999) q[2];
sx q[2];
rz(-1.8450882) q[2];
sx q[2];
rz(-0.14014527) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6631515) q[1];
sx q[1];
rz(-1.3125129) q[1];
sx q[1];
rz(-1.5228935) q[1];
rz(-pi) q[2];
rz(-2.6832163) q[3];
sx q[3];
rz(-1.5273558) q[3];
sx q[3];
rz(0.14940748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0438805) q[2];
sx q[2];
rz(-3.1011797) q[2];
sx q[2];
rz(2.2229693) q[2];
rz(1.0527323) q[3];
sx q[3];
rz(-2.2248) q[3];
sx q[3];
rz(2.2241433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.265428) q[0];
sx q[0];
rz(-0.81887236) q[0];
sx q[0];
rz(2.2692666) q[0];
rz(-1.0126975) q[1];
sx q[1];
rz(-1.6302949) q[1];
sx q[1];
rz(0.33908078) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1992897) q[0];
sx q[0];
rz(-0.45594076) q[0];
sx q[0];
rz(2.845167) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6800687) q[2];
sx q[2];
rz(-2.0465896) q[2];
sx q[2];
rz(-0.78936374) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.87552947) q[1];
sx q[1];
rz(-1.5632946) q[1];
sx q[1];
rz(3.0999567) q[1];
rz(-pi) q[2];
rz(-1.9267155) q[3];
sx q[3];
rz(-2.7971075) q[3];
sx q[3];
rz(3.0657299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.33112153) q[2];
sx q[2];
rz(-1.1473742) q[2];
sx q[2];
rz(-0.60532153) q[2];
rz(2.807054) q[3];
sx q[3];
rz(-2.7693222) q[3];
sx q[3];
rz(0.076586671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89638585) q[0];
sx q[0];
rz(-2.4672282) q[0];
sx q[0];
rz(-1.6687923) q[0];
rz(0.32707602) q[1];
sx q[1];
rz(-1.3027124) q[1];
sx q[1];
rz(-2.8715141) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6657406) q[0];
sx q[0];
rz(-1.0531296) q[0];
sx q[0];
rz(-1.7980018) q[0];
rz(-1.142092) q[2];
sx q[2];
rz(-1.639591) q[2];
sx q[2];
rz(0.69359932) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1535608) q[1];
sx q[1];
rz(-0.84756339) q[1];
sx q[1];
rz(0.60232298) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5669426) q[3];
sx q[3];
rz(-1.4906297) q[3];
sx q[3];
rz(-1.4269258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.104287) q[2];
sx q[2];
rz(-1.9992) q[2];
sx q[2];
rz(-1.8541065) q[2];
rz(-1.2152952) q[3];
sx q[3];
rz(-1.7561965) q[3];
sx q[3];
rz(2.9638885) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6215068) q[0];
sx q[0];
rz(-0.44317133) q[0];
sx q[0];
rz(2.8218414) q[0];
rz(-0.41140914) q[1];
sx q[1];
rz(-1.5965867) q[1];
sx q[1];
rz(0.080332669) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7112996) q[0];
sx q[0];
rz(-1.3724494) q[0];
sx q[0];
rz(1.8848778) q[0];
x q[1];
rz(2.3351203) q[2];
sx q[2];
rz(-1.6739168) q[2];
sx q[2];
rz(-1.1191302) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6957533) q[1];
sx q[1];
rz(-1.8421122) q[1];
sx q[1];
rz(-2.325227) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5307054) q[3];
sx q[3];
rz(-1.3502933) q[3];
sx q[3];
rz(-2.084651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8722998) q[2];
sx q[2];
rz(-0.23739693) q[2];
sx q[2];
rz(3.0966975) q[2];
rz(2.6063555) q[3];
sx q[3];
rz(-1.5072631) q[3];
sx q[3];
rz(3.0285192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.991268) q[0];
sx q[0];
rz(-2.4327705) q[0];
sx q[0];
rz(-2.2659361) q[0];
rz(-0.21370299) q[1];
sx q[1];
rz(-0.49476606) q[1];
sx q[1];
rz(1.5692086) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.101709) q[0];
sx q[0];
rz(-1.3613762) q[0];
sx q[0];
rz(1.8906192) q[0];
rz(-pi) q[1];
rz(-1.3415229) q[2];
sx q[2];
rz(-1.1439644) q[2];
sx q[2];
rz(-0.62139702) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9006648) q[1];
sx q[1];
rz(-0.87605175) q[1];
sx q[1];
rz(1.4347632) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6650136) q[3];
sx q[3];
rz(-1.251128) q[3];
sx q[3];
rz(-1.5227384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.18465061) q[2];
sx q[2];
rz(-0.1447548) q[2];
sx q[2];
rz(1.0682028) q[2];
rz(-0.60643658) q[3];
sx q[3];
rz(-1.9304099) q[3];
sx q[3];
rz(3.1312805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(0.40474263) q[0];
sx q[0];
rz(-0.24979845) q[0];
sx q[0];
rz(-1.3909719) q[0];
rz(-0.5237611) q[1];
sx q[1];
rz(-0.57558376) q[1];
sx q[1];
rz(1.9578804) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42126781) q[0];
sx q[0];
rz(-2.8479332) q[0];
sx q[0];
rz(1.7633506) q[0];
rz(-pi) q[1];
rz(1.6464116) q[2];
sx q[2];
rz(-1.5174149) q[2];
sx q[2];
rz(-0.23301197) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1866926) q[1];
sx q[1];
rz(-2.1789393) q[1];
sx q[1];
rz(-0.18540877) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2565485) q[3];
sx q[3];
rz(-2.9687299) q[3];
sx q[3];
rz(0.96641738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.70048731) q[2];
sx q[2];
rz(-2.5429071) q[2];
sx q[2];
rz(-1.0318476) q[2];
rz(1.4619689) q[3];
sx q[3];
rz(-2.4470191) q[3];
sx q[3];
rz(1.7803378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.4782891) q[0];
sx q[0];
rz(-2.7155868) q[0];
sx q[0];
rz(1.4248832) q[0];
rz(2.5583963) q[1];
sx q[1];
rz(-1.9164663) q[1];
sx q[1];
rz(-3.084175) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9641477) q[0];
sx q[0];
rz(-2.3741407) q[0];
sx q[0];
rz(-1.5657809) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.67441794) q[2];
sx q[2];
rz(-2.4092374) q[2];
sx q[2];
rz(-0.44841097) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.26099953) q[1];
sx q[1];
rz(-0.94124244) q[1];
sx q[1];
rz(1.7884607) q[1];
x q[2];
rz(-1.1693277) q[3];
sx q[3];
rz(-0.9801995) q[3];
sx q[3];
rz(0.99311738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3164369) q[2];
sx q[2];
rz(-1.8437443) q[2];
sx q[2];
rz(1.8334897) q[2];
rz(-1.3206652) q[3];
sx q[3];
rz(-2.2949009) q[3];
sx q[3];
rz(2.2804885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2945781) q[0];
sx q[0];
rz(-2.4569643) q[0];
sx q[0];
rz(1.8735877) q[0];
rz(2.0199203) q[1];
sx q[1];
rz(-1.2753762) q[1];
sx q[1];
rz(-2.4915288) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92597678) q[0];
sx q[0];
rz(-2.4575666) q[0];
sx q[0];
rz(-2.2355272) q[0];
x q[1];
rz(2.4576538) q[2];
sx q[2];
rz(-1.9994447) q[2];
sx q[2];
rz(-0.49698439) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2686058) q[1];
sx q[1];
rz(-2.112231) q[1];
sx q[1];
rz(1.1100015) q[1];
rz(2.2560589) q[3];
sx q[3];
rz(-0.69708744) q[3];
sx q[3];
rz(-3.0892885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.140363) q[2];
sx q[2];
rz(-1.4823806) q[2];
sx q[2];
rz(2.7492827) q[2];
rz(-2.9366734) q[3];
sx q[3];
rz(-0.50060087) q[3];
sx q[3];
rz(-0.71995455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0023163) q[0];
sx q[0];
rz(-1.5520232) q[0];
sx q[0];
rz(-1.4578777) q[0];
rz(2.814759) q[1];
sx q[1];
rz(-1.9953597) q[1];
sx q[1];
rz(-1.0439509) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21699587) q[0];
sx q[0];
rz(-1.1490273) q[0];
sx q[0];
rz(2.934458) q[0];
x q[1];
rz(-1.7108261) q[2];
sx q[2];
rz(-1.1159889) q[2];
sx q[2];
rz(-2.8190913) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.20919101) q[1];
sx q[1];
rz(-1.4889917) q[1];
sx q[1];
rz(1.9664759) q[1];
rz(-pi) q[2];
rz(-2.2208728) q[3];
sx q[3];
rz(-1.3082005) q[3];
sx q[3];
rz(-2.289734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8331208) q[2];
sx q[2];
rz(-2.1127508) q[2];
sx q[2];
rz(-2.2157045) q[2];
rz(-1.067767) q[3];
sx q[3];
rz(-2.0177149) q[3];
sx q[3];
rz(-1.7759751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2866216) q[0];
sx q[0];
rz(-1.6678896) q[0];
sx q[0];
rz(0.58706748) q[0];
rz(2.558737) q[1];
sx q[1];
rz(-0.28935495) q[1];
sx q[1];
rz(0.48097441) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5339587) q[0];
sx q[0];
rz(-1.1635499) q[0];
sx q[0];
rz(-1.7986222) q[0];
x q[1];
rz(-1.0985259) q[2];
sx q[2];
rz(-2.741217) q[2];
sx q[2];
rz(0.085299678) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9945337) q[1];
sx q[1];
rz(-1.9939341) q[1];
sx q[1];
rz(-0.59696859) q[1];
rz(2.7417861) q[3];
sx q[3];
rz(-1.3862063) q[3];
sx q[3];
rz(2.1124482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.19758548) q[2];
sx q[2];
rz(-0.53637594) q[2];
sx q[2];
rz(1.3295004) q[2];
rz(-0.78209376) q[3];
sx q[3];
rz(-0.32556459) q[3];
sx q[3];
rz(-1.9935002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44163497) q[0];
sx q[0];
rz(-1.7229719) q[0];
sx q[0];
rz(-1.4650719) q[0];
rz(1.5326473) q[1];
sx q[1];
rz(-1.9736704) q[1];
sx q[1];
rz(-0.71221487) q[1];
rz(-2.4423626) q[2];
sx q[2];
rz(-2.2742827) q[2];
sx q[2];
rz(-0.90894528) q[2];
rz(-2.6431177) q[3];
sx q[3];
rz(-1.7117906) q[3];
sx q[3];
rz(-1.6556647) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
