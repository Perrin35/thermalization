OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.75582957) q[0];
sx q[0];
rz(-1.4094149) q[0];
sx q[0];
rz(-0.29456079) q[0];
rz(0.28490588) q[1];
sx q[1];
rz(-0.51061881) q[1];
sx q[1];
rz(-2.7198305) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.784214) q[0];
sx q[0];
rz(-1.5445909) q[0];
sx q[0];
rz(-1.5062917) q[0];
rz(-pi) q[1];
rz(1.6099036) q[2];
sx q[2];
rz(-0.88458672) q[2];
sx q[2];
rz(1.4456911) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.64695839) q[1];
sx q[1];
rz(-1.1263493) q[1];
sx q[1];
rz(1.7723945) q[1];
x q[2];
rz(1.4310322) q[3];
sx q[3];
rz(-2.1602109) q[3];
sx q[3];
rz(-2.0624269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1224147) q[2];
sx q[2];
rz(-0.77029595) q[2];
sx q[2];
rz(2.1315234) q[2];
rz(1.4953556) q[3];
sx q[3];
rz(-1.8027179) q[3];
sx q[3];
rz(8*pi/11) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0533957) q[0];
sx q[0];
rz(-1.2258376) q[0];
sx q[0];
rz(-2.2825867) q[0];
rz(-2.7711218) q[1];
sx q[1];
rz(-1.5444688) q[1];
sx q[1];
rz(1.6765615) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1169491) q[0];
sx q[0];
rz(-1.878371) q[0];
sx q[0];
rz(-1.7988388) q[0];
x q[1];
rz(-2.6071762) q[2];
sx q[2];
rz(-0.60306163) q[2];
sx q[2];
rz(0.39819983) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.49711455) q[1];
sx q[1];
rz(-1.9299572) q[1];
sx q[1];
rz(2.4317846) q[1];
rz(-pi) q[2];
rz(0.52821918) q[3];
sx q[3];
rz(-1.8209407) q[3];
sx q[3];
rz(1.6268829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7039965) q[2];
sx q[2];
rz(-1.9251172) q[2];
sx q[2];
rz(-2.583288) q[2];
rz(1.0393418) q[3];
sx q[3];
rz(-1.2149518) q[3];
sx q[3];
rz(0.59282747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6168183) q[0];
sx q[0];
rz(-2.1891948) q[0];
sx q[0];
rz(-2.9779789) q[0];
rz(0.71584654) q[1];
sx q[1];
rz(-0.84638458) q[1];
sx q[1];
rz(1.3735501) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4931902) q[0];
sx q[0];
rz(-1.545854) q[0];
sx q[0];
rz(2.6148952) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9203556) q[2];
sx q[2];
rz(-0.74410838) q[2];
sx q[2];
rz(2.0958054) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.19792412) q[1];
sx q[1];
rz(-2.1197882) q[1];
sx q[1];
rz(-2.3180069) q[1];
rz(-pi) q[2];
rz(1.340629) q[3];
sx q[3];
rz(-1.7880511) q[3];
sx q[3];
rz(2.1118856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0107161) q[2];
sx q[2];
rz(-0.54543442) q[2];
sx q[2];
rz(1.019657) q[2];
rz(-0.9807469) q[3];
sx q[3];
rz(-2.5648983) q[3];
sx q[3];
rz(0.61292928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90137988) q[0];
sx q[0];
rz(-0.69832435) q[0];
sx q[0];
rz(-2.3024978) q[0];
rz(1.6038731) q[1];
sx q[1];
rz(-1.537375) q[1];
sx q[1];
rz(0.61027169) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82204098) q[0];
sx q[0];
rz(-0.10986957) q[0];
sx q[0];
rz(-1.4676453) q[0];
rz(-pi) q[1];
rz(0.82579124) q[2];
sx q[2];
rz(-0.50160393) q[2];
sx q[2];
rz(0.85130168) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.8205386) q[1];
sx q[1];
rz(-0.74756261) q[1];
sx q[1];
rz(1.9588406) q[1];
rz(-0.88704349) q[3];
sx q[3];
rz(-0.700637) q[3];
sx q[3];
rz(-2.5168583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.443976) q[2];
sx q[2];
rz(-0.51091754) q[2];
sx q[2];
rz(-0.237341) q[2];
rz(-2.234263) q[3];
sx q[3];
rz(-1.5932339) q[3];
sx q[3];
rz(1.2835519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5089371) q[0];
sx q[0];
rz(-3.0314358) q[0];
sx q[0];
rz(1.5326112) q[0];
rz(-2.4081047) q[1];
sx q[1];
rz(-1.2669022) q[1];
sx q[1];
rz(-1.8283432) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6812692) q[0];
sx q[0];
rz(-1.8955505) q[0];
sx q[0];
rz(-0.89694174) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0079185) q[2];
sx q[2];
rz(-3.0745227) q[2];
sx q[2];
rz(-1.1941393) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5161799) q[1];
sx q[1];
rz(-1.9429632) q[1];
sx q[1];
rz(-0.401293) q[1];
x q[2];
rz(2.2458514) q[3];
sx q[3];
rz(-0.81479077) q[3];
sx q[3];
rz(0.18273396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0017172) q[2];
sx q[2];
rz(-1.2728609) q[2];
sx q[2];
rz(-0.66429794) q[2];
rz(-2.4364566) q[3];
sx q[3];
rz(-0.64283979) q[3];
sx q[3];
rz(-0.10372182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.076684549) q[0];
sx q[0];
rz(-0.10184558) q[0];
sx q[0];
rz(3.0867807) q[0];
rz(-1.0143771) q[1];
sx q[1];
rz(-1.6631815) q[1];
sx q[1];
rz(0.48318133) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38835634) q[0];
sx q[0];
rz(-1.6325132) q[0];
sx q[0];
rz(1.7524936) q[0];
x q[1];
rz(-2.0334843) q[2];
sx q[2];
rz(-1.6001284) q[2];
sx q[2];
rz(-2.6372452) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.74181611) q[1];
sx q[1];
rz(-2.2175334) q[1];
sx q[1];
rz(2.497655) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11404927) q[3];
sx q[3];
rz(-1.2403135) q[3];
sx q[3];
rz(1.198311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9873535) q[2];
sx q[2];
rz(-2.8581212) q[2];
sx q[2];
rz(-3.0563291) q[2];
rz(-1.2003468) q[3];
sx q[3];
rz(-1.6253358) q[3];
sx q[3];
rz(-2.9912662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7744556) q[0];
sx q[0];
rz(-0.13639233) q[0];
sx q[0];
rz(2.1869587) q[0];
rz(0.57149354) q[1];
sx q[1];
rz(-1.1936455) q[1];
sx q[1];
rz(0.25209299) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82089059) q[0];
sx q[0];
rz(-1.989869) q[0];
sx q[0];
rz(1.4120031) q[0];
rz(-pi) q[1];
x q[1];
rz(0.78105314) q[2];
sx q[2];
rz(-0.99596802) q[2];
sx q[2];
rz(-2.1512254) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.39602173) q[1];
sx q[1];
rz(-2.3704297) q[1];
sx q[1];
rz(1.039617) q[1];
rz(-0.89594706) q[3];
sx q[3];
rz(-0.5657256) q[3];
sx q[3];
rz(1.5783527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.38348848) q[2];
sx q[2];
rz(-2.1540116) q[2];
sx q[2];
rz(-2.2371116) q[2];
rz(1.0036184) q[3];
sx q[3];
rz(-0.86528722) q[3];
sx q[3];
rz(-2.482567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43338183) q[0];
sx q[0];
rz(-3.1383585) q[0];
sx q[0];
rz(-1.7277539) q[0];
rz(-0.66043234) q[1];
sx q[1];
rz(-1.3884037) q[1];
sx q[1];
rz(1.7699014) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28384128) q[0];
sx q[0];
rz(-2.5159266) q[0];
sx q[0];
rz(-1.4197423) q[0];
rz(-pi) q[1];
rz(1.781109) q[2];
sx q[2];
rz(-2.3206707) q[2];
sx q[2];
rz(0.49441499) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0121507) q[1];
sx q[1];
rz(-0.65629849) q[1];
sx q[1];
rz(0.49883573) q[1];
rz(-pi) q[2];
rz(-1.7226108) q[3];
sx q[3];
rz(-1.9278798) q[3];
sx q[3];
rz(-0.15541542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.82683212) q[2];
sx q[2];
rz(-0.62425745) q[2];
sx q[2];
rz(0.39789847) q[2];
rz(2.6509616) q[3];
sx q[3];
rz(-1.5450954) q[3];
sx q[3];
rz(1.3354966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9234377) q[0];
sx q[0];
rz(-1.2175918) q[0];
sx q[0];
rz(-0.23432215) q[0];
rz(-1.461347) q[1];
sx q[1];
rz(-2.318858) q[1];
sx q[1];
rz(-0.27059069) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.035707) q[0];
sx q[0];
rz(-2.8286655) q[0];
sx q[0];
rz(-0.77625113) q[0];
rz(-pi) q[1];
rz(1.2241237) q[2];
sx q[2];
rz(-1.3895831) q[2];
sx q[2];
rz(-2.5172362) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.27442676) q[1];
sx q[1];
rz(-1.5741036) q[1];
sx q[1];
rz(0.13722384) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0398265) q[3];
sx q[3];
rz(-1.093285) q[3];
sx q[3];
rz(-2.4448642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8018735) q[2];
sx q[2];
rz(-0.084771307) q[2];
sx q[2];
rz(-1.6516997) q[2];
rz(-0.70288944) q[3];
sx q[3];
rz(-1.1542164) q[3];
sx q[3];
rz(-1.1606914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(0.40437317) q[0];
sx q[0];
rz(-1.5727366) q[0];
sx q[0];
rz(-0.15429601) q[0];
rz(-2.1986296) q[1];
sx q[1];
rz(-2.0097201) q[1];
sx q[1];
rz(0.74238366) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48001227) q[0];
sx q[0];
rz(-2.1614657) q[0];
sx q[0];
rz(2.3175879) q[0];
x q[1];
rz(-2.7102091) q[2];
sx q[2];
rz(-2.1200392) q[2];
sx q[2];
rz(1.0647578) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.66354499) q[1];
sx q[1];
rz(-1.5596584) q[1];
sx q[1];
rz(2.7958109) q[1];
rz(-1.5998245) q[3];
sx q[3];
rz(-1.4283984) q[3];
sx q[3];
rz(0.27424973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2877038) q[2];
sx q[2];
rz(-2.0159857) q[2];
sx q[2];
rz(1.5596191) q[2];
rz(-2.93086) q[3];
sx q[3];
rz(-0.78453523) q[3];
sx q[3];
rz(2.6081086) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29006526) q[0];
sx q[0];
rz(-1.0354488) q[0];
sx q[0];
rz(-2.5356472) q[0];
rz(1.3416946) q[1];
sx q[1];
rz(-1.1732027) q[1];
sx q[1];
rz(-1.8684594) q[1];
rz(-2.8200091) q[2];
sx q[2];
rz(-2.2939773) q[2];
sx q[2];
rz(1.1800223) q[2];
rz(0.7704173) q[3];
sx q[3];
rz(-0.19376783) q[3];
sx q[3];
rz(2.7289058) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
