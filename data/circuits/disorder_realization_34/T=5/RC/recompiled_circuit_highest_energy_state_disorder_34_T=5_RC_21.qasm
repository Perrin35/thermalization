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
rz(0.96487784) q[0];
rz(2.3896253) q[1];
sx q[1];
rz(-2.7120092) q[1];
sx q[1];
rz(1.0493976) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2222264) q[0];
sx q[0];
rz(-0.99577409) q[0];
sx q[0];
rz(-0.40798835) q[0];
rz(2.9994316) q[2];
sx q[2];
rz(-1.8920676) q[2];
sx q[2];
rz(0.69883332) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.2856372) q[1];
sx q[1];
rz(-2.1826943) q[1];
sx q[1];
rz(-2.6735071) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4780294) q[3];
sx q[3];
rz(-2.1481107) q[3];
sx q[3];
rz(-0.64914671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.45399484) q[2];
sx q[2];
rz(-1.6124085) q[2];
sx q[2];
rz(0.018608658) q[2];
rz(-2.5971557) q[3];
sx q[3];
rz(-2.8114522) q[3];
sx q[3];
rz(-0.64793599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3675156) q[0];
sx q[0];
rz(-1.4544961) q[0];
sx q[0];
rz(2.6065705) q[0];
rz(-2.6016443) q[1];
sx q[1];
rz(-2.5060563) q[1];
sx q[1];
rz(-1.7417057) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28093058) q[0];
sx q[0];
rz(-1.1854404) q[0];
sx q[0];
rz(-1.0866665) q[0];
rz(-pi) q[1];
x q[1];
rz(0.97008791) q[2];
sx q[2];
rz(-2.0510489) q[2];
sx q[2];
rz(0.18829543) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3185841) q[1];
sx q[1];
rz(-2.6806147) q[1];
sx q[1];
rz(-2.5557842) q[1];
rz(-pi) q[2];
rz(2.4055491) q[3];
sx q[3];
rz(-0.76220817) q[3];
sx q[3];
rz(-2.8413642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0743559) q[2];
sx q[2];
rz(-1.8017733) q[2];
sx q[2];
rz(0.92607099) q[2];
rz(-0.98008424) q[3];
sx q[3];
rz(-0.96433774) q[3];
sx q[3];
rz(-0.66942352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3493018) q[0];
sx q[0];
rz(-1.2783569) q[0];
sx q[0];
rz(-1.5128304) q[0];
rz(-0.28383645) q[1];
sx q[1];
rz(-2.2161039) q[1];
sx q[1];
rz(-1.1121174) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4130198) q[0];
sx q[0];
rz(-1.5637102) q[0];
sx q[0];
rz(-1.5636233) q[0];
rz(1.4623423) q[2];
sx q[2];
rz(-1.4949833) q[2];
sx q[2];
rz(2.6159942) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.92924547) q[1];
sx q[1];
rz(-2.622859) q[1];
sx q[1];
rz(-0.38350819) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2073969) q[3];
sx q[3];
rz(-1.3280686) q[3];
sx q[3];
rz(3.0128765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5230368) q[2];
sx q[2];
rz(-1.7967537) q[2];
sx q[2];
rz(1.1019361) q[2];
rz(2.2526422) q[3];
sx q[3];
rz(-2.6933653) q[3];
sx q[3];
rz(2.2811208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76982826) q[0];
sx q[0];
rz(-0.17426057) q[0];
sx q[0];
rz(2.4523822) q[0];
rz(0.31461942) q[1];
sx q[1];
rz(-0.92051053) q[1];
sx q[1];
rz(1.4124195) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1293761) q[0];
sx q[0];
rz(-2.3123992) q[0];
sx q[0];
rz(-1.1393121) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9285703) q[2];
sx q[2];
rz(-1.4771059) q[2];
sx q[2];
rz(2.5793864) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3312348) q[1];
sx q[1];
rz(-2.7572933) q[1];
sx q[1];
rz(-0.15366252) q[1];
rz(-0.19417089) q[3];
sx q[3];
rz(-1.0946858) q[3];
sx q[3];
rz(-0.74060696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.41669258) q[2];
sx q[2];
rz(-0.56513864) q[2];
sx q[2];
rz(0.55595428) q[2];
rz(-1.2712831) q[3];
sx q[3];
rz(-1.2621745) q[3];
sx q[3];
rz(2.8506193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3274662) q[0];
sx q[0];
rz(-0.1411345) q[0];
sx q[0];
rz(2.5471174) q[0];
rz(-0.56198436) q[1];
sx q[1];
rz(-2.2832506) q[1];
sx q[1];
rz(-0.97602239) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41244477) q[0];
sx q[0];
rz(-1.9322104) q[0];
sx q[0];
rz(1.1114208) q[0];
rz(-pi) q[1];
rz(-2.3875434) q[2];
sx q[2];
rz(-1.5496786) q[2];
sx q[2];
rz(-1.4625975) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.70006338) q[1];
sx q[1];
rz(-2.2231327) q[1];
sx q[1];
rz(-1.1928012) q[1];
rz(-pi) q[2];
rz(-0.30140437) q[3];
sx q[3];
rz(-1.3632953) q[3];
sx q[3];
rz(1.0338155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6550265) q[2];
sx q[2];
rz(-1.2532633) q[2];
sx q[2];
rz(2.5308934) q[2];
rz(0.30580172) q[3];
sx q[3];
rz(-0.96547258) q[3];
sx q[3];
rz(1.3293728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.3139451) q[0];
sx q[0];
rz(-3.1231472) q[0];
sx q[0];
rz(-1.6754643) q[0];
rz(-0.77955359) q[1];
sx q[1];
rz(-1.164271) q[1];
sx q[1];
rz(2.4028042) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9290498) q[0];
sx q[0];
rz(-1.3669786) q[0];
sx q[0];
rz(1.8575559) q[0];
rz(-pi) q[1];
rz(-1.5814085) q[2];
sx q[2];
rz(-1.2640177) q[2];
sx q[2];
rz(2.5576484) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2704084) q[1];
sx q[1];
rz(-1.021402) q[1];
sx q[1];
rz(2.710538) q[1];
rz(-pi) q[2];
rz(0.36130623) q[3];
sx q[3];
rz(-1.2194467) q[3];
sx q[3];
rz(1.460399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5256727) q[2];
sx q[2];
rz(-1.2391261) q[2];
sx q[2];
rz(0.8626779) q[2];
rz(-2.681813) q[3];
sx q[3];
rz(-1.516187) q[3];
sx q[3];
rz(-1.9988683) q[3];
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
rz(-2.2343242) q[0];
sx q[0];
rz(-2.7643804) q[0];
sx q[0];
rz(-1.020485) q[0];
rz(-0.057295784) q[1];
sx q[1];
rz(-1.4833996) q[1];
sx q[1];
rz(2.3540156) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5334398) q[0];
sx q[0];
rz(-2.449547) q[0];
sx q[0];
rz(-0.37779053) q[0];
x q[1];
rz(2.8834613) q[2];
sx q[2];
rz(-1.6399334) q[2];
sx q[2];
rz(0.11786945) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0021506) q[1];
sx q[1];
rz(-2.6137335) q[1];
sx q[1];
rz(1.2107641) q[1];
rz(-pi) q[2];
rz(0.98290409) q[3];
sx q[3];
rz(-1.5421151) q[3];
sx q[3];
rz(-0.3405638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.82137498) q[2];
sx q[2];
rz(-1.3221952) q[2];
sx q[2];
rz(2.5569432) q[2];
rz(-2.4892877) q[3];
sx q[3];
rz(-1.9125331) q[3];
sx q[3];
rz(-1.8023796) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24340165) q[0];
sx q[0];
rz(-1.8011872) q[0];
sx q[0];
rz(-0.32824326) q[0];
rz(-0.39168656) q[1];
sx q[1];
rz(-1.990254) q[1];
sx q[1];
rz(1.4788871) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23173702) q[0];
sx q[0];
rz(-1.3074991) q[0];
sx q[0];
rz(-1.80491) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4166862) q[2];
sx q[2];
rz(-2.4838243) q[2];
sx q[2];
rz(-2.7810514) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4305684) q[1];
sx q[1];
rz(-1.4757753) q[1];
sx q[1];
rz(-1.8501758) q[1];
x q[2];
rz(2.9313179) q[3];
sx q[3];
rz(-1.9326903) q[3];
sx q[3];
rz(2.0023605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.11999232) q[2];
sx q[2];
rz(-1.3730717) q[2];
sx q[2];
rz(2.3528986) q[2];
rz(-2.9005519) q[3];
sx q[3];
rz(-0.92625109) q[3];
sx q[3];
rz(0.81361667) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4929844) q[0];
sx q[0];
rz(-0.5383752) q[0];
sx q[0];
rz(0.96187821) q[0];
rz(-2.9073763) q[1];
sx q[1];
rz(-1.1385671) q[1];
sx q[1];
rz(2.403517) q[1];
rz(pi/2) q[2];
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
x q[1];
rz(-0.57508399) q[2];
sx q[2];
rz(-0.7204537) q[2];
sx q[2];
rz(-1.4478113) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8705504) q[1];
sx q[1];
rz(-1.9890474) q[1];
sx q[1];
rz(2.8224432) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.66633114) q[3];
sx q[3];
rz(-1.347858) q[3];
sx q[3];
rz(1.2957038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.48478475) q[2];
sx q[2];
rz(-1.022555) q[2];
sx q[2];
rz(1.8812995) q[2];
rz(-1.8079181) q[3];
sx q[3];
rz(-1.0013872) q[3];
sx q[3];
rz(0.26237747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8986847) q[0];
sx q[0];
rz(-2.3299291) q[0];
sx q[0];
rz(-0.86724487) q[0];
rz(-1.0137089) q[1];
sx q[1];
rz(-1.8127245) q[1];
sx q[1];
rz(-2.0713461) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7263171) q[0];
sx q[0];
rz(-1.4196102) q[0];
sx q[0];
rz(-0.12355208) q[0];
rz(-pi) q[1];
rz(0.70906822) q[2];
sx q[2];
rz(-0.96154172) q[2];
sx q[2];
rz(2.9035062) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5336736) q[1];
sx q[1];
rz(-1.9551139) q[1];
sx q[1];
rz(2.584105) q[1];
rz(-1.6813047) q[3];
sx q[3];
rz(-1.7095057) q[3];
sx q[3];
rz(-1.3121625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9762207) q[2];
sx q[2];
rz(-2.3207211) q[2];
sx q[2];
rz(-1.9099859) q[2];
rz(-1.5008789) q[3];
sx q[3];
rz(-0.21017635) q[3];
sx q[3];
rz(2.4098082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3043542) q[0];
sx q[0];
rz(-1.9914892) q[0];
sx q[0];
rz(-1.7627841) q[0];
rz(2.7267743) q[1];
sx q[1];
rz(-1.3777614) q[1];
sx q[1];
rz(-1.6090964) q[1];
rz(3.0267486) q[2];
sx q[2];
rz(-0.44841246) q[2];
sx q[2];
rz(2.8186225) q[2];
rz(2.5255193) q[3];
sx q[3];
rz(-2.2884634) q[3];
sx q[3];
rz(-0.4175755) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
