OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.22566158) q[0];
sx q[0];
rz(-2.2731279) q[0];
sx q[0];
rz(-2.948569) q[0];
rz(-1.9999737) q[1];
sx q[1];
rz(-2.7116099) q[1];
sx q[1];
rz(-2.4584682) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8773168) q[0];
sx q[0];
rz(-1.6289992) q[0];
sx q[0];
rz(0.067406128) q[0];
rz(2.9002951) q[2];
sx q[2];
rz(-1.9120875) q[2];
sx q[2];
rz(1.2933921) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9477387) q[1];
sx q[1];
rz(-2.1734936) q[1];
sx q[1];
rz(-1.0268474) q[1];
x q[2];
rz(1.0115511) q[3];
sx q[3];
rz(-2.4992001) q[3];
sx q[3];
rz(-2.0410283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0156988) q[2];
sx q[2];
rz(-1.761972) q[2];
sx q[2];
rz(-2.0430298) q[2];
rz(1.0788318) q[3];
sx q[3];
rz(-2.1964985) q[3];
sx q[3];
rz(-1.1014972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9803479) q[0];
sx q[0];
rz(-2.8256567) q[0];
sx q[0];
rz(2.9336477) q[0];
rz(2.5646599) q[1];
sx q[1];
rz(-0.88795841) q[1];
sx q[1];
rz(-1.6764486) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6786878) q[0];
sx q[0];
rz(-1.5683187) q[0];
sx q[0];
rz(-2.4173827) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0589774) q[2];
sx q[2];
rz(-1.3435875) q[2];
sx q[2];
rz(0.28085923) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.84491731) q[1];
sx q[1];
rz(-2.1293318) q[1];
sx q[1];
rz(-1.0122453) q[1];
rz(-pi) q[2];
x q[2];
rz(0.37104718) q[3];
sx q[3];
rz(-0.78073946) q[3];
sx q[3];
rz(-0.77663976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0186105) q[2];
sx q[2];
rz(-0.98615042) q[2];
sx q[2];
rz(0.1097651) q[2];
rz(-0.62260735) q[3];
sx q[3];
rz(-2.770335) q[3];
sx q[3];
rz(1.6842779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1784172) q[0];
sx q[0];
rz(-2.2816179) q[0];
sx q[0];
rz(2.6254568) q[0];
rz(0.57488817) q[1];
sx q[1];
rz(-2.2153885) q[1];
sx q[1];
rz(0.80054545) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.752906) q[0];
sx q[0];
rz(-1.9770925) q[0];
sx q[0];
rz(2.9857062) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3163047) q[2];
sx q[2];
rz(-2.2824259) q[2];
sx q[2];
rz(-0.7691783) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2785981) q[1];
sx q[1];
rz(-1.1578214) q[1];
sx q[1];
rz(-1.5763271) q[1];
x q[2];
rz(0.86977264) q[3];
sx q[3];
rz(-1.8043892) q[3];
sx q[3];
rz(-1.1886532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4642554) q[2];
sx q[2];
rz(-0.3192454) q[2];
sx q[2];
rz(-1.3595954) q[2];
rz(-2.1740186) q[3];
sx q[3];
rz(-1.8680957) q[3];
sx q[3];
rz(-1.7165855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99252218) q[0];
sx q[0];
rz(-1.8795805) q[0];
sx q[0];
rz(2.676679) q[0];
rz(-2.7930296) q[1];
sx q[1];
rz(-2.8788853) q[1];
sx q[1];
rz(2.0565313) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59180075) q[0];
sx q[0];
rz(-1.340197) q[0];
sx q[0];
rz(0.53996284) q[0];
rz(1.9374574) q[2];
sx q[2];
rz(-1.8030093) q[2];
sx q[2];
rz(-2.1249078) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3738721) q[1];
sx q[1];
rz(-1.671119) q[1];
sx q[1];
rz(-0.14212455) q[1];
rz(-pi) q[2];
rz(-2.5892341) q[3];
sx q[3];
rz(-0.68867749) q[3];
sx q[3];
rz(-2.3104582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.00502914) q[2];
sx q[2];
rz(-2.0932784) q[2];
sx q[2];
rz(-0.68112779) q[2];
rz(2.629771) q[3];
sx q[3];
rz(-2.8183283) q[3];
sx q[3];
rz(2.8779023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.27914771) q[0];
sx q[0];
rz(-2.6036766) q[0];
sx q[0];
rz(1.408668) q[0];
rz(-2.7092343) q[1];
sx q[1];
rz(-0.84195781) q[1];
sx q[1];
rz(-2.1599105) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78445804) q[0];
sx q[0];
rz(-1.1103837) q[0];
sx q[0];
rz(0.61607342) q[0];
rz(-pi) q[1];
rz(1.8129187) q[2];
sx q[2];
rz(-1.7244581) q[2];
sx q[2];
rz(-0.15649934) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.423646) q[1];
sx q[1];
rz(-1.7182933) q[1];
sx q[1];
rz(-2.6731078) q[1];
rz(-pi) q[2];
rz(1.4622299) q[3];
sx q[3];
rz(-2.541399) q[3];
sx q[3];
rz(1.1443646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0354054) q[2];
sx q[2];
rz(-1.6550487) q[2];
sx q[2];
rz(0.48941082) q[2];
rz(1.0148467) q[3];
sx q[3];
rz(-2.615052) q[3];
sx q[3];
rz(1.813252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6569825) q[0];
sx q[0];
rz(-2.9841612) q[0];
sx q[0];
rz(2.7691675) q[0];
rz(1.3308446) q[1];
sx q[1];
rz(-2.1061888) q[1];
sx q[1];
rz(0.16528027) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80711354) q[0];
sx q[0];
rz(-1.4757336) q[0];
sx q[0];
rz(-1.6001742) q[0];
rz(-2.4539102) q[2];
sx q[2];
rz(-1.624794) q[2];
sx q[2];
rz(2.8962367) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2416934) q[1];
sx q[1];
rz(-2.4447933) q[1];
sx q[1];
rz(2.4844869) q[1];
rz(-1.1092471) q[3];
sx q[3];
rz(-2.647532) q[3];
sx q[3];
rz(-0.1915313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2281987) q[2];
sx q[2];
rz(-2.1134351) q[2];
sx q[2];
rz(1.6983263) q[2];
rz(2.7741487) q[3];
sx q[3];
rz(-1.2747217) q[3];
sx q[3];
rz(-2.929556) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7845602) q[0];
sx q[0];
rz(-2.050188) q[0];
sx q[0];
rz(0.34926397) q[0];
rz(0.7473942) q[1];
sx q[1];
rz(-0.29577297) q[1];
sx q[1];
rz(0.73648891) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9518785) q[0];
sx q[0];
rz(-1.6197546) q[0];
sx q[0];
rz(-3.1244856) q[0];
rz(-pi) q[1];
rz(0.94524224) q[2];
sx q[2];
rz(-2.7316299) q[2];
sx q[2];
rz(-1.9117102) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6474825) q[1];
sx q[1];
rz(-2.3586015) q[1];
sx q[1];
rz(0.25581911) q[1];
x q[2];
rz(-2.001708) q[3];
sx q[3];
rz(-2.0322554) q[3];
sx q[3];
rz(1.0678837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.13654576) q[2];
sx q[2];
rz(-1.2282635) q[2];
sx q[2];
rz(0.53156701) q[2];
rz(-0.68938869) q[3];
sx q[3];
rz(-1.4644943) q[3];
sx q[3];
rz(-1.2954856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0625967) q[0];
sx q[0];
rz(-1.9364708) q[0];
sx q[0];
rz(1.2980365) q[0];
rz(-0.80728665) q[1];
sx q[1];
rz(-1.9629982) q[1];
sx q[1];
rz(2.2198026) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67547922) q[0];
sx q[0];
rz(-2.4371494) q[0];
sx q[0];
rz(-2.9324313) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2577031) q[2];
sx q[2];
rz(-0.76258341) q[2];
sx q[2];
rz(1.8956172) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7505155) q[1];
sx q[1];
rz(-0.32954307) q[1];
sx q[1];
rz(-2.4604172) q[1];
rz(2.962035) q[3];
sx q[3];
rz(-0.87053821) q[3];
sx q[3];
rz(1.2683271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.90199295) q[2];
sx q[2];
rz(-0.56100503) q[2];
sx q[2];
rz(1.9699338) q[2];
rz(-1.7840067) q[3];
sx q[3];
rz(-1.4566908) q[3];
sx q[3];
rz(-1.4484891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25306025) q[0];
sx q[0];
rz(-1.1760412) q[0];
sx q[0];
rz(-1.6660447) q[0];
rz(-1.5400003) q[1];
sx q[1];
rz(-1.4614636) q[1];
sx q[1];
rz(2.4618861) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98201671) q[0];
sx q[0];
rz(-1.3681108) q[0];
sx q[0];
rz(-0.52635877) q[0];
rz(-pi) q[1];
rz(-2.26608) q[2];
sx q[2];
rz(-1.7098223) q[2];
sx q[2];
rz(-1.8347486) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.45021536) q[1];
sx q[1];
rz(-2.2274349) q[1];
sx q[1];
rz(-2.0037829) q[1];
x q[2];
rz(0.67919517) q[3];
sx q[3];
rz(-1.067357) q[3];
sx q[3];
rz(-0.60590832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0150962) q[2];
sx q[2];
rz(-1.482778) q[2];
sx q[2];
rz(-0.91040197) q[2];
rz(2.4662468) q[3];
sx q[3];
rz(-0.93674913) q[3];
sx q[3];
rz(-2.2909686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.05242059) q[0];
sx q[0];
rz(-1.9071254) q[0];
sx q[0];
rz(-0.7146548) q[0];
rz(-0.71406281) q[1];
sx q[1];
rz(-2.1866182) q[1];
sx q[1];
rz(1.9649327) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6274174) q[0];
sx q[0];
rz(-1.5222933) q[0];
sx q[0];
rz(1.4072627) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5942957) q[2];
sx q[2];
rz(-2.1314869) q[2];
sx q[2];
rz(-2.1292357) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1365876) q[1];
sx q[1];
rz(-1.0141546) q[1];
sx q[1];
rz(-3.0535166) q[1];
x q[2];
rz(2.9858399) q[3];
sx q[3];
rz(-0.51625801) q[3];
sx q[3];
rz(1.932365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.24361336) q[2];
sx q[2];
rz(-1.4695797) q[2];
sx q[2];
rz(-1.127355) q[2];
rz(0.35774287) q[3];
sx q[3];
rz(-2.2704411) q[3];
sx q[3];
rz(-2.0991142) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42416278) q[0];
sx q[0];
rz(-1.9308199) q[0];
sx q[0];
rz(0.45817026) q[0];
rz(-0.39623109) q[1];
sx q[1];
rz(-3.1165262) q[1];
sx q[1];
rz(-2.810626) q[1];
rz(-2.836543) q[2];
sx q[2];
rz(-2.3813644) q[2];
sx q[2];
rz(2.7347953) q[2];
rz(-1.5789923) q[3];
sx q[3];
rz(-1.2200439) q[3];
sx q[3];
rz(-0.72190819) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];