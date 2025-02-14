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
rz(0.76771843) q[1];
sx q[1];
rz(-1.5239198) q[1];
sx q[1];
rz(2.89892) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8382671) q[0];
sx q[0];
rz(-1.6969862) q[0];
sx q[0];
rz(-0.85889205) q[0];
rz(-pi) q[1];
rz(-0.45969339) q[2];
sx q[2];
rz(-1.9003344) q[2];
sx q[2];
rz(-2.7609563) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3585904) q[1];
sx q[1];
rz(-1.5471518) q[1];
sx q[1];
rz(0.35332291) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.30978608) q[3];
sx q[3];
rz(-0.68947809) q[3];
sx q[3];
rz(2.8387217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0365888) q[2];
sx q[2];
rz(-1.236311) q[2];
sx q[2];
rz(0.71199065) q[2];
rz(0.30501929) q[3];
sx q[3];
rz(-0.22235338) q[3];
sx q[3];
rz(0.045507889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(0.73030534) q[0];
sx q[0];
rz(-1.6856025) q[0];
sx q[0];
rz(-2.6645917) q[0];
rz(-pi) q[1];
rz(-0.56053253) q[2];
sx q[2];
rz(-2.2403702) q[2];
sx q[2];
rz(2.5472484) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.89394648) q[1];
sx q[1];
rz(-2.1987913) q[1];
sx q[1];
rz(-1.9878073) q[1];
rz(1.2892339) q[3];
sx q[3];
rz(-0.90471876) q[3];
sx q[3];
rz(-0.56872845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0590608) q[2];
sx q[2];
rz(-0.56190562) q[2];
sx q[2];
rz(-1.833029) q[2];
rz(0.023905309) q[3];
sx q[3];
rz(-1.5289565) q[3];
sx q[3];
rz(1.5628975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4557274) q[0];
sx q[0];
rz(-2.4754334) q[0];
sx q[0];
rz(2.300793) q[0];
rz(-2.1127545) q[1];
sx q[1];
rz(-2.7327171) q[1];
sx q[1];
rz(-1.4604481) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59305916) q[0];
sx q[0];
rz(-1.7139627) q[0];
sx q[0];
rz(-2.6770704) q[0];
rz(-pi) q[1];
rz(0.65030789) q[2];
sx q[2];
rz(-0.57194369) q[2];
sx q[2];
rz(-0.24947333) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1138304) q[1];
sx q[1];
rz(-2.2215469) q[1];
sx q[1];
rz(1.9939994) q[1];
rz(-pi) q[2];
rz(-0.60734235) q[3];
sx q[3];
rz(-1.6173956) q[3];
sx q[3];
rz(0.49639116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.2568405) q[2];
sx q[2];
rz(-1.931087) q[2];
sx q[2];
rz(0.15667008) q[2];
rz(1.6414075) q[3];
sx q[3];
rz(-1.2431966) q[3];
sx q[3];
rz(-0.18178864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0750065) q[0];
sx q[0];
rz(-1.8445419) q[0];
sx q[0];
rz(-0.080667607) q[0];
rz(1.1219885) q[1];
sx q[1];
rz(-2.5009218) q[1];
sx q[1];
rz(-1.5544308) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.156896) q[0];
sx q[0];
rz(-1.6365572) q[0];
sx q[0];
rz(-0.82729152) q[0];
rz(-pi) q[1];
rz(1.5300314) q[2];
sx q[2];
rz(-1.6869378) q[2];
sx q[2];
rz(-0.16414205) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4505541) q[1];
sx q[1];
rz(-2.2802556) q[1];
sx q[1];
rz(-0.67994173) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5347363) q[3];
sx q[3];
rz(-1.8491866) q[3];
sx q[3];
rz(-0.75137072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.88177219) q[2];
sx q[2];
rz(-1.0910923) q[2];
sx q[2];
rz(-2.1176977) q[2];
rz(2.3648868) q[3];
sx q[3];
rz(-1.9909765) q[3];
sx q[3];
rz(-2.1385433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40989947) q[0];
sx q[0];
rz(-2.2875146) q[0];
sx q[0];
rz(-2.5140629) q[0];
rz(-0.69951406) q[1];
sx q[1];
rz(-1.0001837) q[1];
sx q[1];
rz(0.90739179) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8555657) q[0];
sx q[0];
rz(-2.3231909) q[0];
sx q[0];
rz(-1.5159541) q[0];
rz(-2.9512915) q[2];
sx q[2];
rz(-1.7277328) q[2];
sx q[2];
rz(0.67853329) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7414455) q[1];
sx q[1];
rz(-0.65175599) q[1];
sx q[1];
rz(-2.1992127) q[1];
rz(2.1069585) q[3];
sx q[3];
rz(-1.6214633) q[3];
sx q[3];
rz(1.2288273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.13713914) q[2];
sx q[2];
rz(-1.2848102) q[2];
sx q[2];
rz(-2.2966906) q[2];
rz(0.6984624) q[3];
sx q[3];
rz(-0.74028492) q[3];
sx q[3];
rz(-0.79160488) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7591105) q[0];
sx q[0];
rz(-1.4610721) q[0];
sx q[0];
rz(1.8735029) q[0];
rz(-1.5269439) q[1];
sx q[1];
rz(-1.0918278) q[1];
sx q[1];
rz(-0.39438927) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1295107) q[0];
sx q[0];
rz(-2.2144268) q[0];
sx q[0];
rz(3.0846616) q[0];
rz(1.6541078) q[2];
sx q[2];
rz(-2.3040161) q[2];
sx q[2];
rz(1.8064229) q[2];
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
rz(-0.6529863) q[1];
rz(-pi) q[2];
rz(0.36065336) q[3];
sx q[3];
rz(-1.4980199) q[3];
sx q[3];
rz(-2.3520855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7093198) q[2];
sx q[2];
rz(-3.0797112) q[2];
sx q[2];
rz(2.5800932) q[2];
rz(-2.0105441) q[3];
sx q[3];
rz(-2.3384422) q[3];
sx q[3];
rz(0.48857442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5442218) q[0];
sx q[0];
rz(-2.9591296) q[0];
sx q[0];
rz(2.9083948) q[0];
rz(-3.0274262) q[1];
sx q[1];
rz(-2.3515067) q[1];
sx q[1];
rz(-2.9679969) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1602992) q[0];
sx q[0];
rz(-2.3620785) q[0];
sx q[0];
rz(-1.0407838) q[0];
rz(-pi) q[1];
x q[1];
rz(0.87145135) q[2];
sx q[2];
rz(-2.5837499) q[2];
sx q[2];
rz(0.43684549) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6919903) q[1];
sx q[1];
rz(-2.8659645) q[1];
sx q[1];
rz(-0.93494934) q[1];
rz(-3.0775142) q[3];
sx q[3];
rz(-2.9525368) q[3];
sx q[3];
rz(-2.1259675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1058098) q[2];
sx q[2];
rz(-0.95869392) q[2];
sx q[2];
rz(-0.85476533) q[2];
rz(-3.0586976) q[3];
sx q[3];
rz(-1.1707183) q[3];
sx q[3];
rz(-0.054072592) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6524803) q[0];
sx q[0];
rz(-1.2615477) q[0];
sx q[0];
rz(1.9626807) q[0];
rz(-2.1401801) q[1];
sx q[1];
rz(-1.2636355) q[1];
sx q[1];
rz(-0.15403919) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1106873) q[0];
sx q[0];
rz(-2.2461235) q[0];
sx q[0];
rz(1.1481601) q[0];
rz(0.78275795) q[2];
sx q[2];
rz(-2.2680757) q[2];
sx q[2];
rz(-3.1408666) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5402962) q[1];
sx q[1];
rz(-0.5661234) q[1];
sx q[1];
rz(2.152312) q[1];
rz(-pi) q[2];
rz(2.8518744) q[3];
sx q[3];
rz(-1.8236092) q[3];
sx q[3];
rz(3.1186274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4607294) q[2];
sx q[2];
rz(-2.0622084) q[2];
sx q[2];
rz(-2.4737849) q[2];
rz(1.7364511) q[3];
sx q[3];
rz(-1.3677771) q[3];
sx q[3];
rz(-0.63647979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3619096) q[0];
sx q[0];
rz(-1.6661665) q[0];
sx q[0];
rz(-2.5478126) q[0];
rz(-1.8608015) q[1];
sx q[1];
rz(-2.180876) q[1];
sx q[1];
rz(-1.2978172) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4036538) q[0];
sx q[0];
rz(-0.58890051) q[0];
sx q[0];
rz(3.011601) q[0];
rz(1.2601398) q[2];
sx q[2];
rz(-2.1755078) q[2];
sx q[2];
rz(-0.2762281) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.62420708) q[1];
sx q[1];
rz(-0.38100699) q[1];
sx q[1];
rz(-0.34492774) q[1];
rz(-2.3343349) q[3];
sx q[3];
rz(-1.689925) q[3];
sx q[3];
rz(1.8594683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.39591509) q[2];
sx q[2];
rz(-0.021947689) q[2];
sx q[2];
rz(-1.6321261) q[2];
rz(-3.0294026) q[3];
sx q[3];
rz(-1.0271007) q[3];
sx q[3];
rz(1.3794911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6701732) q[0];
sx q[0];
rz(-1.0059953) q[0];
sx q[0];
rz(2.8759586) q[0];
rz(-0.98948014) q[1];
sx q[1];
rz(-1.2629291) q[1];
sx q[1];
rz(-2.8094453) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0743389) q[0];
sx q[0];
rz(-1.5746207) q[0];
sx q[0];
rz(1.2373562) q[0];
rz(-pi) q[1];
x q[1];
rz(0.29075622) q[2];
sx q[2];
rz(-2.1298879) q[2];
sx q[2];
rz(1.5459932) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0451584) q[1];
sx q[1];
rz(-1.2306884) q[1];
sx q[1];
rz(0.051326871) q[1];
rz(-0.048236851) q[3];
sx q[3];
rz(-1.7582446) q[3];
sx q[3];
rz(-1.6316044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0997448) q[2];
sx q[2];
rz(-2.0691278) q[2];
sx q[2];
rz(-0.15038807) q[2];
rz(-2.2930875) q[3];
sx q[3];
rz(-0.34981194) q[3];
sx q[3];
rz(-1.8457886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083241845) q[0];
sx q[0];
rz(-1.8235089) q[0];
sx q[0];
rz(1.9705809) q[0];
rz(1.8104443) q[1];
sx q[1];
rz(-2.34453) q[1];
sx q[1];
rz(-1.1042368) q[1];
rz(0.065481117) q[2];
sx q[2];
rz(-0.83233287) q[2];
sx q[2];
rz(-3.0558791) q[2];
rz(1.5450618) q[3];
sx q[3];
rz(-0.86418695) q[3];
sx q[3];
rz(1.0623111) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
