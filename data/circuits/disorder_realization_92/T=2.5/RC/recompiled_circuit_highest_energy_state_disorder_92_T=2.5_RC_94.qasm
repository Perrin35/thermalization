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
rz(0.25843698) q[0];
sx q[0];
rz(2.8539477) q[0];
sx q[0];
rz(10.532425) q[0];
rz(-1.8183174) q[1];
sx q[1];
rz(-0.32185093) q[1];
sx q[1];
rz(0.62825656) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2520066) q[0];
sx q[0];
rz(-1.4528251) q[0];
sx q[0];
rz(1.749176) q[0];
rz(-0.43810644) q[2];
sx q[2];
rz(-2.3692961) q[2];
sx q[2];
rz(-0.4060678) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.19753708) q[1];
sx q[1];
rz(-1.6070131) q[1];
sx q[1];
rz(-2.9128068) q[1];
x q[2];
rz(0.44094687) q[3];
sx q[3];
rz(-2.0618871) q[3];
sx q[3];
rz(-0.90346149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2884752) q[2];
sx q[2];
rz(-1.709781) q[2];
sx q[2];
rz(0.34118578) q[2];
rz(-0.8902542) q[3];
sx q[3];
rz(-1.8743926) q[3];
sx q[3];
rz(-2.6064579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61966908) q[0];
sx q[0];
rz(-2.4578019) q[0];
sx q[0];
rz(-2.4271915) q[0];
rz(-1.5996541) q[1];
sx q[1];
rz(-2.6456412) q[1];
sx q[1];
rz(3.0468859) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9908087) q[0];
sx q[0];
rz(-1.2512733) q[0];
sx q[0];
rz(3.0421769) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.96927283) q[2];
sx q[2];
rz(-1.2899295) q[2];
sx q[2];
rz(1.2796677) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4058806) q[1];
sx q[1];
rz(-2.7196214) q[1];
sx q[1];
rz(-0.60776819) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.27908947) q[3];
sx q[3];
rz(-2.6153339) q[3];
sx q[3];
rz(0.56015362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.16227214) q[2];
sx q[2];
rz(-1.8956192) q[2];
sx q[2];
rz(0.47404131) q[2];
rz(0.76215172) q[3];
sx q[3];
rz(-2.5049329) q[3];
sx q[3];
rz(0.4711841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8629465) q[0];
sx q[0];
rz(-3.0248248) q[0];
sx q[0];
rz(-2.2886724) q[0];
rz(0.22311738) q[1];
sx q[1];
rz(-0.47839636) q[1];
sx q[1];
rz(0.51582897) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86403041) q[0];
sx q[0];
rz(-1.6245961) q[0];
sx q[0];
rz(0.050431099) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2081096) q[2];
sx q[2];
rz(-1.7778983) q[2];
sx q[2];
rz(1.1400956) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.307439) q[1];
sx q[1];
rz(-2.1527228) q[1];
sx q[1];
rz(-0.80640275) q[1];
rz(-pi) q[2];
x q[2];
rz(0.34308691) q[3];
sx q[3];
rz(-0.90215397) q[3];
sx q[3];
rz(-0.17681387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8082661) q[2];
sx q[2];
rz(-2.1648679) q[2];
sx q[2];
rz(3.0466363) q[2];
rz(0.44148463) q[3];
sx q[3];
rz(-1.7850826) q[3];
sx q[3];
rz(-1.9879742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3917291) q[0];
sx q[0];
rz(-0.57498217) q[0];
sx q[0];
rz(1.0561426) q[0];
rz(-0.9437584) q[1];
sx q[1];
rz(-0.25139233) q[1];
sx q[1];
rz(-1.6897078) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7565931) q[0];
sx q[0];
rz(-2.2384522) q[0];
sx q[0];
rz(0.7636015) q[0];
rz(0.46197666) q[2];
sx q[2];
rz(-1.2179795) q[2];
sx q[2];
rz(1.8227673) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6744212) q[1];
sx q[1];
rz(-1.8113664) q[1];
sx q[1];
rz(-1.6075385) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5783674) q[3];
sx q[3];
rz(-1.6814582) q[3];
sx q[3];
rz(-0.31294926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4132495) q[2];
sx q[2];
rz(-0.80060935) q[2];
sx q[2];
rz(-0.70017868) q[2];
rz(-0.87107301) q[3];
sx q[3];
rz(-1.0659404) q[3];
sx q[3];
rz(-1.357249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0102608) q[0];
sx q[0];
rz(-2.6483674) q[0];
sx q[0];
rz(-2.6569271) q[0];
rz(-2.2635745) q[1];
sx q[1];
rz(-1.1596102) q[1];
sx q[1];
rz(2.7427618) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6054886) q[0];
sx q[0];
rz(-1.7057422) q[0];
sx q[0];
rz(-0.55079726) q[0];
rz(-pi) q[1];
rz(2.7473248) q[2];
sx q[2];
rz(-0.74299251) q[2];
sx q[2];
rz(-2.791648) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4927401) q[1];
sx q[1];
rz(-1.2034186) q[1];
sx q[1];
rz(0.87441625) q[1];
rz(0.54581235) q[3];
sx q[3];
rz(-0.94792507) q[3];
sx q[3];
rz(-2.2198077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6681119) q[2];
sx q[2];
rz(-1.6010189) q[2];
sx q[2];
rz(-2.2637746) q[2];
rz(-0.38464883) q[3];
sx q[3];
rz(-2.2967702) q[3];
sx q[3];
rz(1.1788684) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28410742) q[0];
sx q[0];
rz(-2.6326038) q[0];
sx q[0];
rz(-0.2704764) q[0];
rz(0.30666223) q[1];
sx q[1];
rz(-2.6276734) q[1];
sx q[1];
rz(-2.564863) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3424042) q[0];
sx q[0];
rz(-1.4831838) q[0];
sx q[0];
rz(-1.3848928) q[0];
rz(-0.88097303) q[2];
sx q[2];
rz(-1.2861797) q[2];
sx q[2];
rz(1.6135482) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2564078) q[1];
sx q[1];
rz(-1.1756304) q[1];
sx q[1];
rz(-2.1705519) q[1];
rz(-pi) q[2];
rz(-1.5848198) q[3];
sx q[3];
rz(-1.2360667) q[3];
sx q[3];
rz(2.1317792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7755255) q[2];
sx q[2];
rz(-1.0696573) q[2];
sx q[2];
rz(-0.97037399) q[2];
rz(2.7905285) q[3];
sx q[3];
rz(-2.909436) q[3];
sx q[3];
rz(-0.27950132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2358667) q[0];
sx q[0];
rz(-0.38385639) q[0];
sx q[0];
rz(2.6406777) q[0];
rz(-0.72055912) q[1];
sx q[1];
rz(-0.84545285) q[1];
sx q[1];
rz(0.93773425) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7739627) q[0];
sx q[0];
rz(-1.367462) q[0];
sx q[0];
rz(-3.0995661) q[0];
x q[1];
rz(0.65848041) q[2];
sx q[2];
rz(-0.57504932) q[2];
sx q[2];
rz(2.5031896) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5590458) q[1];
sx q[1];
rz(-1.2101047) q[1];
sx q[1];
rz(-0.8949032) q[1];
rz(-pi) q[2];
rz(3.1295007) q[3];
sx q[3];
rz(-2.4547057) q[3];
sx q[3];
rz(1.505132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.131989) q[2];
sx q[2];
rz(-0.14543532) q[2];
sx q[2];
rz(-0.84849882) q[2];
rz(-0.84544539) q[3];
sx q[3];
rz(-2.4852821) q[3];
sx q[3];
rz(1.9158624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41035143) q[0];
sx q[0];
rz(-2.1827965) q[0];
sx q[0];
rz(0.89286667) q[0];
rz(-0.061575312) q[1];
sx q[1];
rz(-0.67191809) q[1];
sx q[1];
rz(-0.16199131) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91770691) q[0];
sx q[0];
rz(-0.10334238) q[0];
sx q[0];
rz(1.4204106) q[0];
x q[1];
rz(1.4561557) q[2];
sx q[2];
rz(-0.87213665) q[2];
sx q[2];
rz(1.4574775) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2714241) q[1];
sx q[1];
rz(-1.4335911) q[1];
sx q[1];
rz(-0.31183621) q[1];
rz(-pi) q[2];
rz(1.4668767) q[3];
sx q[3];
rz(-1.6641518) q[3];
sx q[3];
rz(0.68478497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.050345) q[2];
sx q[2];
rz(-0.42751905) q[2];
sx q[2];
rz(2.4166935) q[2];
rz(2.8643518) q[3];
sx q[3];
rz(-0.96479779) q[3];
sx q[3];
rz(-3.0550756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8898833) q[0];
sx q[0];
rz(-2.4885663) q[0];
sx q[0];
rz(3.1157893) q[0];
rz(-1.39894) q[1];
sx q[1];
rz(-1.6394697) q[1];
sx q[1];
rz(3.0339411) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7079332) q[0];
sx q[0];
rz(-1.7042394) q[0];
sx q[0];
rz(-0.99876499) q[0];
rz(2.0094821) q[2];
sx q[2];
rz(-2.0033547) q[2];
sx q[2];
rz(0.74299845) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.15756179) q[1];
sx q[1];
rz(-1.5347297) q[1];
sx q[1];
rz(-0.47337516) q[1];
rz(-pi) q[2];
rz(1.4571545) q[3];
sx q[3];
rz(-1.5530619) q[3];
sx q[3];
rz(-2.4271698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4699576) q[2];
sx q[2];
rz(-1.2592955) q[2];
sx q[2];
rz(2.0476445) q[2];
rz(-0.57540244) q[3];
sx q[3];
rz(-2.3016774) q[3];
sx q[3];
rz(-1.1168787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.6082918) q[0];
sx q[0];
rz(-2.5550483) q[0];
sx q[0];
rz(2.4285512) q[0];
rz(2.9167922) q[1];
sx q[1];
rz(-0.59200042) q[1];
sx q[1];
rz(-1.0932659) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2042052) q[0];
sx q[0];
rz(-2.2474216) q[0];
sx q[0];
rz(-0.94955523) q[0];
rz(-pi) q[1];
x q[1];
rz(0.20168882) q[2];
sx q[2];
rz(-1.8400157) q[2];
sx q[2];
rz(1.9112916) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8902407) q[1];
sx q[1];
rz(-1.4595493) q[1];
sx q[1];
rz(-0.45808582) q[1];
rz(1.9841927) q[3];
sx q[3];
rz(-1.6781312) q[3];
sx q[3];
rz(1.3232376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.0010058086) q[2];
sx q[2];
rz(-0.69585496) q[2];
sx q[2];
rz(2.801753) q[2];
rz(-0.042424399) q[3];
sx q[3];
rz(-1.2845311) q[3];
sx q[3];
rz(2.5137918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9545659) q[0];
sx q[0];
rz(-1.5196336) q[0];
sx q[0];
rz(-1.2484311) q[0];
rz(-0.78202248) q[1];
sx q[1];
rz(-0.93688688) q[1];
sx q[1];
rz(-0.72601906) q[1];
rz(-2.4154631) q[2];
sx q[2];
rz(-2.6604322) q[2];
sx q[2];
rz(-0.12728035) q[2];
rz(2.0841712) q[3];
sx q[3];
rz(-0.48013018) q[3];
sx q[3];
rz(-2.6063802) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
