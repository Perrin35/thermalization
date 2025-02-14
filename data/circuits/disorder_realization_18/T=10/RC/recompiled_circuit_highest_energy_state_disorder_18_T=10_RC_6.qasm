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
rz(-1.350116) q[0];
sx q[0];
rz(3.8173563) q[0];
sx q[0];
rz(9.4326333) q[0];
rz(-2.3738742) q[1];
sx q[1];
rz(-1.6176728) q[1];
sx q[1];
rz(-2.89892) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15930106) q[0];
sx q[0];
rz(-0.86573273) q[0];
sx q[0];
rz(-0.16601913) q[0];
rz(-pi) q[1];
rz(-0.45969339) q[2];
sx q[2];
rz(-1.2412583) q[2];
sx q[2];
rz(-0.38063637) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.85181273) q[1];
sx q[1];
rz(-2.7875126) q[1];
sx q[1];
rz(3.0733527) q[1];
x q[2];
rz(-2.4759329) q[3];
sx q[3];
rz(-1.7659617) q[3];
sx q[3];
rz(1.5100556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1050038) q[2];
sx q[2];
rz(-1.236311) q[2];
sx q[2];
rz(-0.71199065) q[2];
rz(-2.8365734) q[3];
sx q[3];
rz(-2.9192393) q[3];
sx q[3];
rz(3.0960848) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5104093) q[0];
sx q[0];
rz(-1.864186) q[0];
sx q[0];
rz(-1.3674059) q[0];
rz(-3.101688) q[1];
sx q[1];
rz(-0.80723643) q[1];
sx q[1];
rz(2.5423999) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78135787) q[0];
sx q[0];
rz(-1.0971945) q[0];
sx q[0];
rz(-1.4417157) q[0];
rz(-pi) q[1];
rz(-0.81920287) q[2];
sx q[2];
rz(-1.1408198) q[2];
sx q[2];
rz(2.5366304) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6015905) q[1];
sx q[1];
rz(-0.7379325) q[1];
sx q[1];
rz(0.5088536) q[1];
rz(-pi) q[2];
rz(2.8017524) q[3];
sx q[3];
rz(-2.426894) q[3];
sx q[3];
rz(3.0106737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0590608) q[2];
sx q[2];
rz(-2.579687) q[2];
sx q[2];
rz(1.3085636) q[2];
rz(-3.1176873) q[3];
sx q[3];
rz(-1.5289565) q[3];
sx q[3];
rz(1.5628975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6858653) q[0];
sx q[0];
rz(-0.66615921) q[0];
sx q[0];
rz(2.300793) q[0];
rz(-1.0288382) q[1];
sx q[1];
rz(-0.40887555) q[1];
sx q[1];
rz(-1.4604481) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.235229) q[0];
sx q[0];
rz(-1.1113941) q[0];
sx q[0];
rz(1.7306585) q[0];
x q[1];
rz(1.1991791) q[2];
sx q[2];
rz(-2.016168) q[2];
sx q[2];
rz(-2.1567313) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9509995) q[1];
sx q[1];
rz(-1.237932) q[1];
sx q[1];
rz(-2.4458363) q[1];
rz(-1.5140686) q[3];
sx q[3];
rz(-0.9642082) q[3];
sx q[3];
rz(1.1067672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.2568405) q[2];
sx q[2];
rz(-1.931087) q[2];
sx q[2];
rz(-2.9849226) q[2];
rz(-1.6414075) q[3];
sx q[3];
rz(-1.898396) q[3];
sx q[3];
rz(-0.18178864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0665862) q[0];
sx q[0];
rz(-1.8445419) q[0];
sx q[0];
rz(3.060925) q[0];
rz(-2.0196041) q[1];
sx q[1];
rz(-2.5009218) q[1];
sx q[1];
rz(1.5871619) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64645221) q[0];
sx q[0];
rz(-2.3123154) q[0];
sx q[0];
rz(-0.089228169) q[0];
rz(0.3360668) q[2];
sx q[2];
rz(-3.0185351) q[2];
sx q[2];
rz(2.6390136) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5836559) q[1];
sx q[1];
rz(-0.93975818) q[1];
sx q[1];
rz(0.93871745) q[1];
rz(-pi) q[2];
rz(1.2359592) q[3];
sx q[3];
rz(-0.99042884) q[3];
sx q[3];
rz(2.133647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2598205) q[2];
sx q[2];
rz(-1.0910923) q[2];
sx q[2];
rz(-2.1176977) q[2];
rz(-2.3648868) q[3];
sx q[3];
rz(-1.1506162) q[3];
sx q[3];
rz(1.0030494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-0.40989947) q[0];
sx q[0];
rz(-2.2875146) q[0];
sx q[0];
rz(-0.62752974) q[0];
rz(0.69951406) q[1];
sx q[1];
rz(-1.0001837) q[1];
sx q[1];
rz(-0.90739179) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32226792) q[0];
sx q[0];
rz(-1.5307679) q[0];
sx q[0];
rz(0.75314523) q[0];
x q[1];
rz(0.19030119) q[2];
sx q[2];
rz(-1.4138599) q[2];
sx q[2];
rz(-0.67853329) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6945362) q[1];
sx q[1];
rz(-1.9354104) q[1];
sx q[1];
rz(-2.1237808) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6697407) q[3];
sx q[3];
rz(-0.53831783) q[3];
sx q[3];
rz(-2.7145999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0044535) q[2];
sx q[2];
rz(-1.2848102) q[2];
sx q[2];
rz(-0.84490204) q[2];
rz(2.4431303) q[3];
sx q[3];
rz(-2.4013077) q[3];
sx q[3];
rz(2.3499878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38248211) q[0];
sx q[0];
rz(-1.4610721) q[0];
sx q[0];
rz(1.8735029) q[0];
rz(-1.6146487) q[1];
sx q[1];
rz(-1.0918278) q[1];
sx q[1];
rz(0.39438927) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6170664) q[0];
sx q[0];
rz(-1.6163278) q[0];
sx q[0];
rz(-2.2152053) q[0];
rz(-pi) q[1];
rz(-3.0494681) q[2];
sx q[2];
rz(-2.4045304) q[2];
sx q[2];
rz(-1.9305522) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.39526597) q[1];
sx q[1];
rz(-1.526064) q[1];
sx q[1];
rz(2.4886063) q[1];
rz(-pi) q[2];
rz(-1.4930356) q[3];
sx q[3];
rz(-1.2111411) q[3];
sx q[3];
rz(0.80870562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7093198) q[2];
sx q[2];
rz(-3.0797112) q[2];
sx q[2];
rz(-0.56149948) q[2];
rz(-1.1310486) q[3];
sx q[3];
rz(-2.3384422) q[3];
sx q[3];
rz(2.6530182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5973709) q[0];
sx q[0];
rz(-0.18246305) q[0];
sx q[0];
rz(-2.9083948) q[0];
rz(3.0274262) q[1];
sx q[1];
rz(-2.3515067) q[1];
sx q[1];
rz(-0.17359576) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9812935) q[0];
sx q[0];
rz(-0.7795142) q[0];
sx q[0];
rz(-2.1008089) q[0];
rz(-pi) q[1];
rz(2.2701413) q[2];
sx q[2];
rz(-2.5837499) q[2];
sx q[2];
rz(2.7047472) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6379695) q[1];
sx q[1];
rz(-1.733128) q[1];
sx q[1];
rz(-1.7945402) q[1];
rz(-3.0775142) q[3];
sx q[3];
rz(-0.18905583) q[3];
sx q[3];
rz(-1.0156251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0357828) q[2];
sx q[2];
rz(-2.1828987) q[2];
sx q[2];
rz(-2.2868273) q[2];
rz(-0.0828951) q[3];
sx q[3];
rz(-1.9708743) q[3];
sx q[3];
rz(3.0875201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4891124) q[0];
sx q[0];
rz(-1.2615477) q[0];
sx q[0];
rz(-1.1789119) q[0];
rz(-1.0014125) q[1];
sx q[1];
rz(-1.8779571) q[1];
sx q[1];
rz(-0.15403919) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0309054) q[0];
sx q[0];
rz(-2.2461235) q[0];
sx q[0];
rz(-1.9934325) q[0];
rz(-2.270584) q[2];
sx q[2];
rz(-0.99620512) q[2];
sx q[2];
rz(2.1432997) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5402962) q[1];
sx q[1];
rz(-2.5754693) q[1];
sx q[1];
rz(-0.98928063) q[1];
x q[2];
rz(0.2897183) q[3];
sx q[3];
rz(-1.3179835) q[3];
sx q[3];
rz(3.1186274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.68086326) q[2];
sx q[2];
rz(-1.0793842) q[2];
sx q[2];
rz(2.4737849) q[2];
rz(-1.4051416) q[3];
sx q[3];
rz(-1.3677771) q[3];
sx q[3];
rz(-0.63647979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3619096) q[0];
sx q[0];
rz(-1.6661665) q[0];
sx q[0];
rz(0.59378004) q[0];
rz(1.8608015) q[1];
sx q[1];
rz(-2.180876) q[1];
sx q[1];
rz(-1.8437754) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5820044) q[0];
sx q[0];
rz(-0.98752092) q[0];
sx q[0];
rz(-1.6571664) q[0];
x q[1];
rz(-0.4164575) q[2];
sx q[2];
rz(-0.67085941) q[2];
sx q[2];
rz(2.351298) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5173856) q[1];
sx q[1];
rz(-0.38100699) q[1];
sx q[1];
rz(-0.34492774) q[1];
rz(-2.3343349) q[3];
sx q[3];
rz(-1.4516677) q[3];
sx q[3];
rz(-1.8594683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.39591509) q[2];
sx q[2];
rz(-3.119645) q[2];
sx q[2];
rz(-1.6321261) q[2];
rz(-0.11219003) q[3];
sx q[3];
rz(-1.0271007) q[3];
sx q[3];
rz(-1.3794911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47141948) q[0];
sx q[0];
rz(-1.0059953) q[0];
sx q[0];
rz(-2.8759586) q[0];
rz(-2.1521125) q[1];
sx q[1];
rz(-1.2629291) q[1];
sx q[1];
rz(-0.33214733) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6393747) q[0];
sx q[0];
rz(-1.9042339) q[0];
sx q[0];
rz(0.0040472814) q[0];
rz(-pi) q[1];
x q[1];
rz(0.99224706) q[2];
sx q[2];
rz(-1.8162811) q[2];
sx q[2];
rz(0.13260376) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.89236081) q[1];
sx q[1];
rz(-2.7977825) q[1];
sx q[1];
rz(1.4268141) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0933558) q[3];
sx q[3];
rz(-1.3833481) q[3];
sx q[3];
rz(-1.6316044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0997448) q[2];
sx q[2];
rz(-2.0691278) q[2];
sx q[2];
rz(-2.9912046) q[2];
rz(-2.2930875) q[3];
sx q[3];
rz(-0.34981194) q[3];
sx q[3];
rz(-1.8457886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083241845) q[0];
sx q[0];
rz(-1.8235089) q[0];
sx q[0];
rz(1.9705809) q[0];
rz(1.3311483) q[1];
sx q[1];
rz(-0.7970627) q[1];
sx q[1];
rz(2.0373559) q[1];
rz(1.4990357) q[2];
sx q[2];
rz(-0.74081479) q[2];
sx q[2];
rz(-0.011394636) q[2];
rz(3.1114586) q[3];
sx q[3];
rz(-2.4345955) q[3];
sx q[3];
rz(1.101936) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
