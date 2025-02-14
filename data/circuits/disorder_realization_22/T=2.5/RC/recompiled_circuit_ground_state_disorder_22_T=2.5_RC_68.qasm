OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.9353256) q[0];
sx q[0];
rz(4.7630035) q[0];
sx q[0];
rz(8.0061316) q[0];
rz(-3.0980134) q[1];
sx q[1];
rz(-1.8585304) q[1];
sx q[1];
rz(2.9820774) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7107726) q[0];
sx q[0];
rz(-2.9732467) q[0];
sx q[0];
rz(-1.7414067) q[0];
x q[1];
rz(2.7829917) q[2];
sx q[2];
rz(-2.498543) q[2];
sx q[2];
rz(-1.3728764) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2251896) q[1];
sx q[1];
rz(-1.708751) q[1];
sx q[1];
rz(-2.0504584) q[1];
x q[2];
rz(2.2920424) q[3];
sx q[3];
rz(-1.8462984) q[3];
sx q[3];
rz(1.4412137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.688545) q[2];
sx q[2];
rz(-1.9828372) q[2];
sx q[2];
rz(0.090252074) q[2];
rz(3.0669751) q[3];
sx q[3];
rz(-1.6736232) q[3];
sx q[3];
rz(2.727865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2523044) q[0];
sx q[0];
rz(-1.1587208) q[0];
sx q[0];
rz(1.6878376) q[0];
rz(2.3989035) q[1];
sx q[1];
rz(-1.7748666) q[1];
sx q[1];
rz(2.3765366) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1677645) q[0];
sx q[0];
rz(-0.59685464) q[0];
sx q[0];
rz(0.11886998) q[0];
rz(-pi) q[1];
rz(2.2751669) q[2];
sx q[2];
rz(-1.0110628) q[2];
sx q[2];
rz(1.0407966) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2606352) q[1];
sx q[1];
rz(-0.85399125) q[1];
sx q[1];
rz(-0.44181255) q[1];
rz(-pi) q[2];
x q[2];
rz(0.54799796) q[3];
sx q[3];
rz(-1.573085) q[3];
sx q[3];
rz(1.8016004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.226977) q[2];
sx q[2];
rz(-1.0116297) q[2];
sx q[2];
rz(-1.4545308) q[2];
rz(-2.4948273) q[3];
sx q[3];
rz(-0.55964595) q[3];
sx q[3];
rz(-1.2796992) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0890559) q[0];
sx q[0];
rz(-1.4707969) q[0];
sx q[0];
rz(-3.0964417) q[0];
rz(-1.2308944) q[1];
sx q[1];
rz(-1.8321593) q[1];
sx q[1];
rz(1.0567573) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089356747) q[0];
sx q[0];
rz(-2.1187021) q[0];
sx q[0];
rz(-2.7578081) q[0];
rz(2.0665523) q[2];
sx q[2];
rz(-0.98598749) q[2];
sx q[2];
rz(-0.57384795) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5505728) q[1];
sx q[1];
rz(-1.3168087) q[1];
sx q[1];
rz(-2.2588782) q[1];
rz(0.23356861) q[3];
sx q[3];
rz(-2.2456777) q[3];
sx q[3];
rz(-2.6949331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7354273) q[2];
sx q[2];
rz(-0.96144095) q[2];
sx q[2];
rz(-0.65262922) q[2];
rz(-1.8485707) q[3];
sx q[3];
rz(-0.76904622) q[3];
sx q[3];
rz(0.098171083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3339313) q[0];
sx q[0];
rz(-0.46285358) q[0];
sx q[0];
rz(1.3941143) q[0];
rz(1.3035424) q[1];
sx q[1];
rz(-1.1640254) q[1];
sx q[1];
rz(-2.0533144) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9435642) q[0];
sx q[0];
rz(-1.8450415) q[0];
sx q[0];
rz(-2.345577) q[0];
rz(-pi) q[1];
rz(0.91371961) q[2];
sx q[2];
rz(-0.83763257) q[2];
sx q[2];
rz(-0.73071551) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6034607) q[1];
sx q[1];
rz(-0.96718057) q[1];
sx q[1];
rz(0.97410283) q[1];
rz(-pi) q[2];
rz(2.6422464) q[3];
sx q[3];
rz(-1.5245228) q[3];
sx q[3];
rz(-0.83016992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3950222) q[2];
sx q[2];
rz(-1.3777379) q[2];
sx q[2];
rz(-2.6611888) q[2];
rz(0.26614842) q[3];
sx q[3];
rz(-1.1689309) q[3];
sx q[3];
rz(2.2973785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-3.0323623) q[0];
sx q[0];
rz(-2.5205595) q[0];
sx q[0];
rz(-1.9367628) q[0];
rz(-3.0989528) q[1];
sx q[1];
rz(-1.7667021) q[1];
sx q[1];
rz(0.94243324) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5980127) q[0];
sx q[0];
rz(-2.6930598) q[0];
sx q[0];
rz(-1.2537812) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2082926) q[2];
sx q[2];
rz(-1.0663053) q[2];
sx q[2];
rz(-1.8613033) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.79980274) q[1];
sx q[1];
rz(-1.1936578) q[1];
sx q[1];
rz(1.0339292) q[1];
rz(-pi) q[2];
rz(-1.4616555) q[3];
sx q[3];
rz(-0.93459826) q[3];
sx q[3];
rz(-0.075625758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.11613906) q[2];
sx q[2];
rz(-2.4567273) q[2];
sx q[2];
rz(-1.9174513) q[2];
rz(-2.4141566) q[3];
sx q[3];
rz(-1.6614611) q[3];
sx q[3];
rz(-3.091403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.3392451) q[0];
sx q[0];
rz(-2.6363063) q[0];
sx q[0];
rz(-2.8832866) q[0];
rz(-0.91336617) q[1];
sx q[1];
rz(-0.86562997) q[1];
sx q[1];
rz(2.1736653) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8081849) q[0];
sx q[0];
rz(-0.70751941) q[0];
sx q[0];
rz(0.61715166) q[0];
x q[1];
rz(-2.2931678) q[2];
sx q[2];
rz(-1.8471187) q[2];
sx q[2];
rz(1.6538617) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.591394) q[1];
sx q[1];
rz(-2.0481143) q[1];
sx q[1];
rz(1.5309019) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0276035) q[3];
sx q[3];
rz(-0.99811998) q[3];
sx q[3];
rz(-3.0610994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4919058) q[2];
sx q[2];
rz(-0.68444362) q[2];
sx q[2];
rz(0.41927949) q[2];
rz(0.78420091) q[3];
sx q[3];
rz(-2.1214285) q[3];
sx q[3];
rz(1.3721344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.066561617) q[0];
sx q[0];
rz(-0.12396585) q[0];
sx q[0];
rz(-0.01874622) q[0];
rz(-3.0491414) q[1];
sx q[1];
rz(-1.3321313) q[1];
sx q[1];
rz(-0.094955347) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92871794) q[0];
sx q[0];
rz(-1.8935907) q[0];
sx q[0];
rz(-2.5514437) q[0];
rz(-pi) q[1];
x q[1];
rz(0.59573458) q[2];
sx q[2];
rz(-1.2168665) q[2];
sx q[2];
rz(-2.6806493) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5508037) q[1];
sx q[1];
rz(-1.992464) q[1];
sx q[1];
rz(-1.8487537) q[1];
x q[2];
rz(1.9390887) q[3];
sx q[3];
rz(-0.42516685) q[3];
sx q[3];
rz(1.396338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2350601) q[2];
sx q[2];
rz(-1.1842714) q[2];
sx q[2];
rz(2.2825784) q[2];
rz(-0.61339316) q[3];
sx q[3];
rz(-0.20039138) q[3];
sx q[3];
rz(1.8092417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(2.5253946) q[0];
sx q[0];
rz(-1.2735676) q[0];
sx q[0];
rz(1.6599576) q[0];
rz(2.0741277) q[1];
sx q[1];
rz(-0.72507247) q[1];
sx q[1];
rz(-1.8225089) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6292012) q[0];
sx q[0];
rz(-0.90921697) q[0];
sx q[0];
rz(-0.11920905) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5379792) q[2];
sx q[2];
rz(-0.60342646) q[2];
sx q[2];
rz(-0.44925424) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.47506491) q[1];
sx q[1];
rz(-1.804978) q[1];
sx q[1];
rz(0.98824309) q[1];
rz(-pi) q[2];
rz(0.35679166) q[3];
sx q[3];
rz(-1.4513375) q[3];
sx q[3];
rz(0.91537634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8760425) q[2];
sx q[2];
rz(-0.37028131) q[2];
sx q[2];
rz(-3.1004356) q[2];
rz(1.1228849) q[3];
sx q[3];
rz(-1.4833781) q[3];
sx q[3];
rz(2.6179204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48162833) q[0];
sx q[0];
rz(-2.1631503) q[0];
sx q[0];
rz(1.6690669) q[0];
rz(-0.72498471) q[1];
sx q[1];
rz(-1.1712733) q[1];
sx q[1];
rz(-0.57458007) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6642484) q[0];
sx q[0];
rz(-1.9084404) q[0];
sx q[0];
rz(1.6214423) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9926316) q[2];
sx q[2];
rz(-1.8293194) q[2];
sx q[2];
rz(-2.0542893) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9758399) q[1];
sx q[1];
rz(-2.4226592) q[1];
sx q[1];
rz(0.045054212) q[1];
rz(-pi) q[2];
rz(0.88317849) q[3];
sx q[3];
rz(-1.6856582) q[3];
sx q[3];
rz(-2.3695947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.64728111) q[2];
sx q[2];
rz(-2.8936671) q[2];
sx q[2];
rz(2.1809123) q[2];
rz(-2.658355) q[3];
sx q[3];
rz(-1.5471231) q[3];
sx q[3];
rz(-1.5794992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1305337) q[0];
sx q[0];
rz(-2.647825) q[0];
sx q[0];
rz(2.6997987) q[0];
rz(0.68253016) q[1];
sx q[1];
rz(-1.6923169) q[1];
sx q[1];
rz(-1.0535047) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58013142) q[0];
sx q[0];
rz(-2.2765358) q[0];
sx q[0];
rz(-0.34039008) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0526673) q[2];
sx q[2];
rz(-1.94238) q[2];
sx q[2];
rz(-0.099612923) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.23020959) q[1];
sx q[1];
rz(-0.95048302) q[1];
sx q[1];
rz(2.7432454) q[1];
x q[2];
rz(-2.9656762) q[3];
sx q[3];
rz(-0.39947594) q[3];
sx q[3];
rz(-2.7586058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0423476) q[2];
sx q[2];
rz(-1.8727563) q[2];
sx q[2];
rz(-2.919) q[2];
rz(1.7706722) q[3];
sx q[3];
rz(-1.418908) q[3];
sx q[3];
rz(-0.62209779) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2400146) q[0];
sx q[0];
rz(-1.0116901) q[0];
sx q[0];
rz(-2.1160175) q[0];
rz(1.1657794) q[1];
sx q[1];
rz(-0.60973254) q[1];
sx q[1];
rz(2.9235074) q[1];
rz(0.72543421) q[2];
sx q[2];
rz(-0.67580961) q[2];
sx q[2];
rz(0.11563822) q[2];
rz(3.141116) q[3];
sx q[3];
rz(-1.2491741) q[3];
sx q[3];
rz(0.76330976) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
