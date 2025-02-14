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
rz(0.024700392) q[0];
sx q[0];
rz(8.165897) q[0];
sx q[0];
rz(11.293647) q[0];
rz(1.3920353) q[1];
sx q[1];
rz(-1.4104383) q[1];
sx q[1];
rz(2.2478204) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0451242) q[0];
sx q[0];
rz(-2.22825) q[0];
sx q[0];
rz(2.3923001) q[0];
x q[1];
rz(1.3145163) q[2];
sx q[2];
rz(-1.8602784) q[2];
sx q[2];
rz(0.50965259) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8402579) q[1];
sx q[1];
rz(-0.71941676) q[1];
sx q[1];
rz(-2.4459228) q[1];
rz(-0.046272101) q[3];
sx q[3];
rz(-0.36847575) q[3];
sx q[3];
rz(-0.61532332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.36665234) q[2];
sx q[2];
rz(-1.4619991) q[2];
sx q[2];
rz(-0.53202638) q[2];
rz(0.73412791) q[3];
sx q[3];
rz(-2.8668154) q[3];
sx q[3];
rz(2.0348569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72320139) q[0];
sx q[0];
rz(-0.98576236) q[0];
sx q[0];
rz(0.080408737) q[0];
rz(0.53781992) q[1];
sx q[1];
rz(-1.0686914) q[1];
sx q[1];
rz(-0.72371662) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2576669) q[0];
sx q[0];
rz(-0.91086713) q[0];
sx q[0];
rz(-1.0038478) q[0];
rz(-pi) q[1];
rz(0.89214708) q[2];
sx q[2];
rz(-2.3513459) q[2];
sx q[2];
rz(-2.9228766) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.557505) q[1];
sx q[1];
rz(-0.98230108) q[1];
sx q[1];
rz(-3.1397538) q[1];
rz(0.51809727) q[3];
sx q[3];
rz(-2.2541917) q[3];
sx q[3];
rz(0.70113126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2878652) q[2];
sx q[2];
rz(-2.008805) q[2];
sx q[2];
rz(2.2354324) q[2];
rz(-0.3624889) q[3];
sx q[3];
rz(-1.5972842) q[3];
sx q[3];
rz(-0.046886142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2694117) q[0];
sx q[0];
rz(-1.325664) q[0];
sx q[0];
rz(1.0915225) q[0];
rz(-0.12256924) q[1];
sx q[1];
rz(-1.7054319) q[1];
sx q[1];
rz(-1.3465808) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99341644) q[0];
sx q[0];
rz(-1.5412586) q[0];
sx q[0];
rz(1.9269153) q[0];
rz(-pi) q[1];
rz(-2.7287219) q[2];
sx q[2];
rz(-0.70253583) q[2];
sx q[2];
rz(1.2324049) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.31981644) q[1];
sx q[1];
rz(-0.21142658) q[1];
sx q[1];
rz(1.3216622) q[1];
rz(-2.9503959) q[3];
sx q[3];
rz(-0.53732291) q[3];
sx q[3];
rz(2.4464698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.71318212) q[2];
sx q[2];
rz(-1.1298263) q[2];
sx q[2];
rz(0.23207363) q[2];
rz(-1.9706767) q[3];
sx q[3];
rz(-2.5210095) q[3];
sx q[3];
rz(0.28461972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.93037) q[0];
sx q[0];
rz(-0.79293293) q[0];
sx q[0];
rz(-1.5027745) q[0];
rz(-0.59208313) q[1];
sx q[1];
rz(-1.9860257) q[1];
sx q[1];
rz(2.2519462) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98418035) q[0];
sx q[0];
rz(-2.0286467) q[0];
sx q[0];
rz(-2.1628987) q[0];
rz(-pi) q[1];
rz(1.1355033) q[2];
sx q[2];
rz(-2.6084628) q[2];
sx q[2];
rz(1.0077623) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0977323) q[1];
sx q[1];
rz(-1.7508888) q[1];
sx q[1];
rz(-1.8921683) q[1];
rz(-pi) q[2];
rz(2.3208404) q[3];
sx q[3];
rz(-0.53583691) q[3];
sx q[3];
rz(-1.395063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4952116) q[2];
sx q[2];
rz(-1.6782328) q[2];
sx q[2];
rz(0.67592534) q[2];
rz(1.7050788) q[3];
sx q[3];
rz(-2.8016475) q[3];
sx q[3];
rz(0.16726141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9306358) q[0];
sx q[0];
rz(-2.8093331) q[0];
sx q[0];
rz(-0.19530547) q[0];
rz(3.0209814) q[1];
sx q[1];
rz(-0.21574012) q[1];
sx q[1];
rz(-0.19048555) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57225655) q[0];
sx q[0];
rz(-1.4829206) q[0];
sx q[0];
rz(-1.6356619) q[0];
rz(-pi) q[1];
rz(0.67660013) q[2];
sx q[2];
rz(-1.6381421) q[2];
sx q[2];
rz(1.2852033) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.021564158) q[1];
sx q[1];
rz(-3.0282674) q[1];
sx q[1];
rz(-0.47685195) q[1];
rz(-pi) q[2];
rz(-0.56511527) q[3];
sx q[3];
rz(-2.4451974) q[3];
sx q[3];
rz(-0.90000641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6327989) q[2];
sx q[2];
rz(-1.6570647) q[2];
sx q[2];
rz(1.18139) q[2];
rz(-1.7975636) q[3];
sx q[3];
rz(-2.400178) q[3];
sx q[3];
rz(0.77176362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33264273) q[0];
sx q[0];
rz(-1.3430261) q[0];
sx q[0];
rz(2.0020265) q[0];
rz(-0.66967669) q[1];
sx q[1];
rz(-0.91595903) q[1];
sx q[1];
rz(-1.12961) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78167262) q[0];
sx q[0];
rz(-1.2914713) q[0];
sx q[0];
rz(1.2646227) q[0];
rz(0.30252075) q[2];
sx q[2];
rz(-2.4008958) q[2];
sx q[2];
rz(0.64756913) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.46907963) q[1];
sx q[1];
rz(-1.7653078) q[1];
sx q[1];
rz(-0.37176337) q[1];
x q[2];
rz(-1.4127168) q[3];
sx q[3];
rz(-1.7493093) q[3];
sx q[3];
rz(-1.3016537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.42664042) q[2];
sx q[2];
rz(-1.4098097) q[2];
sx q[2];
rz(-0.40531522) q[2];
rz(0.82695588) q[3];
sx q[3];
rz(-0.6260286) q[3];
sx q[3];
rz(-0.85957447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3655132) q[0];
sx q[0];
rz(-0.022553355) q[0];
sx q[0];
rz(0.93589163) q[0];
rz(-2.2731958) q[1];
sx q[1];
rz(-2.6135542) q[1];
sx q[1];
rz(-1.7220928) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0763797) q[0];
sx q[0];
rz(-1.2166436) q[0];
sx q[0];
rz(-0.54535086) q[0];
rz(-pi) q[1];
rz(-1.0638191) q[2];
sx q[2];
rz(-1.3666743) q[2];
sx q[2];
rz(-2.7046471) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6554113) q[1];
sx q[1];
rz(-2.2194194) q[1];
sx q[1];
rz(2.5470665) q[1];
rz(-pi) q[2];
rz(-0.21162947) q[3];
sx q[3];
rz(-3.1397925) q[3];
sx q[3];
rz(-1.3789919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.3203656) q[2];
sx q[2];
rz(-1.3996404) q[2];
sx q[2];
rz(1.0199176) q[2];
rz(0.50926456) q[3];
sx q[3];
rz(-1.441322) q[3];
sx q[3];
rz(0.078484623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.7996063) q[0];
sx q[0];
rz(-1.6260363) q[0];
sx q[0];
rz(1.1588143) q[0];
rz(-0.47053567) q[1];
sx q[1];
rz(-1.7262986) q[1];
sx q[1];
rz(1.4656969) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3780791) q[0];
sx q[0];
rz(-1.5944885) q[0];
sx q[0];
rz(0.050934249) q[0];
rz(2.7297425) q[2];
sx q[2];
rz(-0.73397103) q[2];
sx q[2];
rz(0.07871544) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1701057) q[1];
sx q[1];
rz(-2.4772812) q[1];
sx q[1];
rz(-3.0168135) q[1];
x q[2];
rz(1.2837778) q[3];
sx q[3];
rz(-0.94285184) q[3];
sx q[3];
rz(-2.7046741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1450242) q[2];
sx q[2];
rz(-1.327876) q[2];
sx q[2];
rz(0.81929755) q[2];
rz(-0.79646349) q[3];
sx q[3];
rz(-0.4070681) q[3];
sx q[3];
rz(1.4888633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9374989) q[0];
sx q[0];
rz(-3.0452073) q[0];
sx q[0];
rz(-2.3199484) q[0];
rz(-0.67283982) q[1];
sx q[1];
rz(-2.5655589) q[1];
sx q[1];
rz(1.0427262) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13623986) q[0];
sx q[0];
rz(-1.7757105) q[0];
sx q[0];
rz(1.4968275) q[0];
rz(-pi) q[1];
x q[1];
rz(0.77206467) q[2];
sx q[2];
rz(-2.0520718) q[2];
sx q[2];
rz(1.5779863) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1112177) q[1];
sx q[1];
rz(-1.0008604) q[1];
sx q[1];
rz(-2.6968677) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25298499) q[3];
sx q[3];
rz(-1.4232985) q[3];
sx q[3];
rz(-1.9852888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.78160703) q[2];
sx q[2];
rz(-1.8693417) q[2];
sx q[2];
rz(0.84235111) q[2];
rz(0.53884566) q[3];
sx q[3];
rz(-1.6539961) q[3];
sx q[3];
rz(2.6174788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46776029) q[0];
sx q[0];
rz(-0.14586511) q[0];
sx q[0];
rz(-1.9061331) q[0];
rz(0.94888672) q[1];
sx q[1];
rz(-2.3088539) q[1];
sx q[1];
rz(-1.6611151) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7296071) q[0];
sx q[0];
rz(-2.3959037) q[0];
sx q[0];
rz(-1.9870158) q[0];
rz(0.5625426) q[2];
sx q[2];
rz(-2.1058561) q[2];
sx q[2];
rz(-1.8059274) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.93714156) q[1];
sx q[1];
rz(-1.4960663) q[1];
sx q[1];
rz(0.16507574) q[1];
x q[2];
rz(-1.3289024) q[3];
sx q[3];
rz(-1.2513565) q[3];
sx q[3];
rz(-2.0177359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4296055) q[2];
sx q[2];
rz(-1.7179855) q[2];
sx q[2];
rz(-0.31022662) q[2];
rz(-0.45983908) q[3];
sx q[3];
rz(-0.55791563) q[3];
sx q[3];
rz(1.7673813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1772169) q[0];
sx q[0];
rz(-1.1342659) q[0];
sx q[0];
rz(-2.076617) q[0];
rz(0.89827697) q[1];
sx q[1];
rz(-2.5986462) q[1];
sx q[1];
rz(2.6959261) q[1];
rz(-2.2353371) q[2];
sx q[2];
rz(-1.2438085) q[2];
sx q[2];
rz(-0.96240733) q[2];
rz(0.38705714) q[3];
sx q[3];
rz(-2.195812) q[3];
sx q[3];
rz(0.14433911) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
