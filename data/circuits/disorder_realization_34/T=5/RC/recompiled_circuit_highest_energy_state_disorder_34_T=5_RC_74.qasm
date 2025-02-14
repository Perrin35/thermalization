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
rz(1.3833157) q[0];
sx q[0];
rz(-1.436469) q[0];
sx q[0];
rz(-2.1767148) q[0];
rz(-0.75196737) q[1];
sx q[1];
rz(-0.42958346) q[1];
sx q[1];
rz(2.092195) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9193662) q[0];
sx q[0];
rz(-2.1458186) q[0];
sx q[0];
rz(2.7336043) q[0];
rz(-pi) q[1];
rz(1.9732835) q[2];
sx q[2];
rz(-0.3503198) q[2];
sx q[2];
rz(-0.27327368) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5677898) q[1];
sx q[1];
rz(-1.192523) q[1];
sx q[1];
rz(2.2371464) q[1];
x q[2];
rz(-0.14127381) q[3];
sx q[3];
rz(-2.5577099) q[3];
sx q[3];
rz(-0.48030765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6875978) q[2];
sx q[2];
rz(-1.5291841) q[2];
sx q[2];
rz(-3.122984) q[2];
rz(0.54443693) q[3];
sx q[3];
rz(-0.33014044) q[3];
sx q[3];
rz(-2.4936567) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7740771) q[0];
sx q[0];
rz(-1.4544961) q[0];
sx q[0];
rz(-2.6065705) q[0];
rz(-2.6016443) q[1];
sx q[1];
rz(-2.5060563) q[1];
sx q[1];
rz(1.3998869) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2311418) q[0];
sx q[0];
rz(-2.5325091) q[0];
sx q[0];
rz(2.2876491) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.82615699) q[2];
sx q[2];
rz(-0.75010502) q[2];
sx q[2];
rz(-0.78924417) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.78838879) q[1];
sx q[1];
rz(-1.8192768) q[1];
sx q[1];
rz(0.39239359) q[1];
rz(-pi) q[2];
rz(-2.4055491) q[3];
sx q[3];
rz(-0.76220817) q[3];
sx q[3];
rz(2.8413642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0743559) q[2];
sx q[2];
rz(-1.3398193) q[2];
sx q[2];
rz(-2.2155217) q[2];
rz(0.98008424) q[3];
sx q[3];
rz(-0.96433774) q[3];
sx q[3];
rz(0.66942352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3493018) q[0];
sx q[0];
rz(-1.2783569) q[0];
sx q[0];
rz(1.6287623) q[0];
rz(0.28383645) q[1];
sx q[1];
rz(-0.92548871) q[1];
sx q[1];
rz(-1.1121174) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4130198) q[0];
sx q[0];
rz(-1.5778825) q[0];
sx q[0];
rz(1.5636233) q[0];
x q[1];
rz(-0.9588963) q[2];
sx q[2];
rz(-3.009353) q[2];
sx q[2];
rz(-2.7041777) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.49431153) q[1];
sx q[1];
rz(-2.0485281) q[1];
sx q[1];
rz(-1.7812438) q[1];
rz(-0.25896163) q[3];
sx q[3];
rz(-1.9230712) q[3];
sx q[3];
rz(-1.6083628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6185559) q[2];
sx q[2];
rz(-1.3448389) q[2];
sx q[2];
rz(1.1019361) q[2];
rz(-0.88895041) q[3];
sx q[3];
rz(-0.44822732) q[3];
sx q[3];
rz(-2.2811208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76982826) q[0];
sx q[0];
rz(-0.17426057) q[0];
sx q[0];
rz(-0.68921047) q[0];
rz(2.8269732) q[1];
sx q[1];
rz(-0.92051053) q[1];
sx q[1];
rz(-1.4124195) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5310881) q[0];
sx q[0];
rz(-0.83685368) q[0];
sx q[0];
rz(2.713301) q[0];
rz(-pi) q[1];
rz(0.41823776) q[2];
sx q[2];
rz(-2.9091638) q[2];
sx q[2];
rz(-1.4168036) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9758983) q[1];
sx q[1];
rz(-1.1912575) q[1];
sx q[1];
rz(-1.5089773) q[1];
rz(1.9288428) q[3];
sx q[3];
rz(-0.51135495) q[3];
sx q[3];
rz(2.8062888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7249001) q[2];
sx q[2];
rz(-2.576454) q[2];
sx q[2];
rz(-2.5856384) q[2];
rz(-1.2712831) q[3];
sx q[3];
rz(-1.8794182) q[3];
sx q[3];
rz(-2.8506193) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3274662) q[0];
sx q[0];
rz(-3.0004582) q[0];
sx q[0];
rz(0.59447527) q[0];
rz(2.5796083) q[1];
sx q[1];
rz(-0.85834208) q[1];
sx q[1];
rz(-2.1655703) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41244477) q[0];
sx q[0];
rz(-1.2093822) q[0];
sx q[0];
rz(-1.1114208) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.030841737) q[2];
sx q[2];
rz(-2.387306) q[2];
sx q[2];
rz(-3.0109143) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1072717) q[1];
sx q[1];
rz(-1.8684505) q[1];
sx q[1];
rz(-0.68796449) q[1];
rz(-1.7877949) q[3];
sx q[3];
rz(-1.2760538) q[3];
sx q[3];
rz(0.47302305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.48656616) q[2];
sx q[2];
rz(-1.8883294) q[2];
sx q[2];
rz(-0.61069926) q[2];
rz(2.8357909) q[3];
sx q[3];
rz(-2.1761201) q[3];
sx q[3];
rz(1.3293728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82764757) q[0];
sx q[0];
rz(-0.018445404) q[0];
sx q[0];
rz(-1.6754643) q[0];
rz(2.3620391) q[1];
sx q[1];
rz(-1.9773217) q[1];
sx q[1];
rz(0.73878845) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9290498) q[0];
sx q[0];
rz(-1.3669786) q[0];
sx q[0];
rz(-1.8575559) q[0];
rz(3.1081057) q[2];
sx q[2];
rz(-2.8346363) q[2];
sx q[2];
rz(2.592776) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.87118426) q[1];
sx q[1];
rz(-2.1201907) q[1];
sx q[1];
rz(2.710538) q[1];
rz(-pi) q[2];
rz(-0.8035369) q[3];
sx q[3];
rz(-0.49852405) q[3];
sx q[3];
rz(-0.62832181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5256727) q[2];
sx q[2];
rz(-1.2391261) q[2];
sx q[2];
rz(2.2789148) q[2];
rz(-0.45977965) q[3];
sx q[3];
rz(-1.6254057) q[3];
sx q[3];
rz(-1.9988683) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90726844) q[0];
sx q[0];
rz(-2.7643804) q[0];
sx q[0];
rz(1.020485) q[0];
rz(3.0842969) q[1];
sx q[1];
rz(-1.4833996) q[1];
sx q[1];
rz(-0.78757706) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5334398) q[0];
sx q[0];
rz(-2.449547) q[0];
sx q[0];
rz(0.37779053) q[0];
rz(-0.25813132) q[2];
sx q[2];
rz(-1.6399334) q[2];
sx q[2];
rz(-3.0237232) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.139442) q[1];
sx q[1];
rz(-2.6137335) q[1];
sx q[1];
rz(-1.2107641) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1586886) q[3];
sx q[3];
rz(-1.5421151) q[3];
sx q[3];
rz(0.3405638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.82137498) q[2];
sx q[2];
rz(-1.8193974) q[2];
sx q[2];
rz(-2.5569432) q[2];
rz(-0.65230495) q[3];
sx q[3];
rz(-1.9125331) q[3];
sx q[3];
rz(1.8023796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.898191) q[0];
sx q[0];
rz(-1.3404055) q[0];
sx q[0];
rz(-2.8133494) q[0];
rz(-0.39168656) q[1];
sx q[1];
rz(-1.1513386) q[1];
sx q[1];
rz(1.6627056) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9098556) q[0];
sx q[0];
rz(-1.3074991) q[0];
sx q[0];
rz(1.3366827) q[0];
x q[1];
rz(2.4166862) q[2];
sx q[2];
rz(-2.4838243) q[2];
sx q[2];
rz(2.7810514) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.45962557) q[1];
sx q[1];
rz(-2.8468968) q[1];
sx q[1];
rz(1.9035643) q[1];
rz(2.9313179) q[3];
sx q[3];
rz(-1.9326903) q[3];
sx q[3];
rz(2.0023605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.11999232) q[2];
sx q[2];
rz(-1.3730717) q[2];
sx q[2];
rz(-0.78869406) q[2];
rz(-2.9005519) q[3];
sx q[3];
rz(-0.92625109) q[3];
sx q[3];
rz(0.81361667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64860827) q[0];
sx q[0];
rz(-2.6032175) q[0];
sx q[0];
rz(-2.1797144) q[0];
rz(2.9073763) q[1];
sx q[1];
rz(-2.0030256) q[1];
sx q[1];
rz(-0.73807565) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80958145) q[0];
sx q[0];
rz(-1.9652848) q[0];
sx q[0];
rz(1.8277192) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0162651) q[2];
sx q[2];
rz(-0.98410749) q[2];
sx q[2];
rz(-2.4054804) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9540087) q[1];
sx q[1];
rz(-2.6212647) q[1];
sx q[1];
rz(0.95611683) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8515737) q[3];
sx q[3];
rz(-0.92377907) q[3];
sx q[3];
rz(-3.0386277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.48478475) q[2];
sx q[2];
rz(-1.022555) q[2];
sx q[2];
rz(-1.2602932) q[2];
rz(1.8079181) q[3];
sx q[3];
rz(-1.0013872) q[3];
sx q[3];
rz(2.8792152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-1.242908) q[0];
sx q[0];
rz(-0.81166357) q[0];
sx q[0];
rz(0.86724487) q[0];
rz(1.0137089) q[1];
sx q[1];
rz(-1.8127245) q[1];
sx q[1];
rz(2.0713461) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7263171) q[0];
sx q[0];
rz(-1.4196102) q[0];
sx q[0];
rz(3.0180406) q[0];
rz(-pi) q[1];
x q[1];
rz(2.321601) q[2];
sx q[2];
rz(-0.89897663) q[2];
sx q[2];
rz(-1.9209678) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.874843) q[1];
sx q[1];
rz(-1.0581746) q[1];
sx q[1];
rz(-1.1260518) q[1];
rz(2.4730014) q[3];
sx q[3];
rz(-2.9644659) q[3];
sx q[3];
rz(1.9884381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.16537198) q[2];
sx q[2];
rz(-2.3207211) q[2];
sx q[2];
rz(-1.9099859) q[2];
rz(-1.6407137) q[3];
sx q[3];
rz(-2.9314163) q[3];
sx q[3];
rz(-0.73178449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-0.83723849) q[0];
sx q[0];
rz(-1.1501034) q[0];
sx q[0];
rz(1.3788086) q[0];
rz(2.7267743) q[1];
sx q[1];
rz(-1.3777614) q[1];
sx q[1];
rz(-1.6090964) q[1];
rz(1.6258705) q[2];
sx q[2];
rz(-1.1255506) q[2];
sx q[2];
rz(-0.19565565) q[2];
rz(-0.61607331) q[3];
sx q[3];
rz(-2.2884634) q[3];
sx q[3];
rz(-0.4175755) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
