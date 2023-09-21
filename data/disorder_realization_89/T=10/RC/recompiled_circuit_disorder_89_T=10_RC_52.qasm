OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3857631) q[0];
sx q[0];
rz(-1.7321777) q[0];
sx q[0];
rz(-2.8470319) q[0];
rz(0.28490588) q[1];
sx q[1];
rz(-0.51061881) q[1];
sx q[1];
rz(0.4217622) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3573787) q[0];
sx q[0];
rz(-1.5970018) q[0];
sx q[0];
rz(1.635301) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6099036) q[2];
sx q[2];
rz(-2.2570059) q[2];
sx q[2];
rz(1.4456911) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4946343) q[1];
sx q[1];
rz(-1.1263493) q[1];
sx q[1];
rz(-1.7723945) q[1];
rz(-1.4310322) q[3];
sx q[3];
rz(-2.1602109) q[3];
sx q[3];
rz(-1.0791657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0191779) q[2];
sx q[2];
rz(-0.77029595) q[2];
sx q[2];
rz(-1.0100693) q[2];
rz(-1.4953556) q[3];
sx q[3];
rz(-1.8027179) q[3];
sx q[3];
rz(-8*pi/11) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0533957) q[0];
sx q[0];
rz(-1.915755) q[0];
sx q[0];
rz(-2.2825867) q[0];
rz(2.7711218) q[1];
sx q[1];
rz(-1.5971239) q[1];
sx q[1];
rz(-1.4650311) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61629907) q[0];
sx q[0];
rz(-1.787961) q[0];
sx q[0];
rz(0.3152245) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2334521) q[2];
sx q[2];
rz(-2.0806899) q[2];
sx q[2];
rz(-2.9166729) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.49711455) q[1];
sx q[1];
rz(-1.2116355) q[1];
sx q[1];
rz(2.4317846) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6133735) q[3];
sx q[3];
rz(-1.8209407) q[3];
sx q[3];
rz(-1.6268829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4375962) q[2];
sx q[2];
rz(-1.9251172) q[2];
sx q[2];
rz(-0.55830467) q[2];
rz(2.1022508) q[3];
sx q[3];
rz(-1.9266409) q[3];
sx q[3];
rz(-2.5487652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6168183) q[0];
sx q[0];
rz(-0.95239788) q[0];
sx q[0];
rz(-2.9779789) q[0];
rz(-2.4257461) q[1];
sx q[1];
rz(-0.84638458) q[1];
sx q[1];
rz(-1.7680426) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87953075) q[0];
sx q[0];
rz(-2.6143605) q[0];
sx q[0];
rz(0.049588163) q[0];
x q[1];
rz(-1.221237) q[2];
sx q[2];
rz(-0.74410838) q[2];
sx q[2];
rz(-1.0457872) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.19792412) q[1];
sx q[1];
rz(-1.0218044) q[1];
sx q[1];
rz(2.3180069) q[1];
rz(-pi) q[2];
rz(-0.80188607) q[3];
sx q[3];
rz(-2.826414) q[3];
sx q[3];
rz(-1.2847881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.13087656) q[2];
sx q[2];
rz(-0.54543442) q[2];
sx q[2];
rz(2.1219357) q[2];
rz(-2.1608458) q[3];
sx q[3];
rz(-2.5648983) q[3];
sx q[3];
rz(2.5286634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2402128) q[0];
sx q[0];
rz(-2.4432683) q[0];
sx q[0];
rz(-2.3024978) q[0];
rz(-1.6038731) q[1];
sx q[1];
rz(-1.537375) q[1];
sx q[1];
rz(2.531321) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85128879) q[0];
sx q[0];
rz(-1.5595058) q[0];
sx q[0];
rz(-1.6800866) q[0];
rz(-1.9539815) q[2];
sx q[2];
rz(-1.9028579) q[2];
sx q[2];
rz(-3.1021169) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3210541) q[1];
sx q[1];
rz(-2.39403) q[1];
sx q[1];
rz(-1.9588406) q[1];
rz(2.2545492) q[3];
sx q[3];
rz(-2.4409557) q[3];
sx q[3];
rz(-0.62473434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6976167) q[2];
sx q[2];
rz(-2.6306751) q[2];
sx q[2];
rz(-2.9042517) q[2];
rz(2.234263) q[3];
sx q[3];
rz(-1.5932339) q[3];
sx q[3];
rz(1.8580407) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5089371) q[0];
sx q[0];
rz(-0.11015686) q[0];
sx q[0];
rz(-1.5326112) q[0];
rz(-2.4081047) q[1];
sx q[1];
rz(-1.8746904) q[1];
sx q[1];
rz(-1.3132494) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7816313) q[0];
sx q[0];
rz(-2.2035723) q[0];
sx q[0];
rz(-0.40681337) q[0];
rz(-pi) q[1];
rz(1.6275431) q[2];
sx q[2];
rz(-1.6065671) q[2];
sx q[2];
rz(-2.2030731) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.65391958) q[1];
sx q[1];
rz(-2.6012585) q[1];
sx q[1];
rz(2.3565156) q[1];
x q[2];
rz(2.26236) q[3];
sx q[3];
rz(-2.0428265) q[3];
sx q[3];
rz(1.8902682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1398754) q[2];
sx q[2];
rz(-1.2728609) q[2];
sx q[2];
rz(0.66429794) q[2];
rz(2.4364566) q[3];
sx q[3];
rz(-2.4987529) q[3];
sx q[3];
rz(3.0378708) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0649081) q[0];
sx q[0];
rz(-3.0397471) q[0];
sx q[0];
rz(-0.054811906) q[0];
rz(1.0143771) q[1];
sx q[1];
rz(-1.4784112) q[1];
sx q[1];
rz(-2.6584113) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2830402) q[0];
sx q[0];
rz(-2.9498093) q[0];
sx q[0];
rz(-1.2412846) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0334843) q[2];
sx q[2];
rz(-1.5414642) q[2];
sx q[2];
rz(-0.5043475) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.74181611) q[1];
sx q[1];
rz(-2.2175334) q[1];
sx q[1];
rz(0.6439376) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2504962) q[3];
sx q[3];
rz(-2.7926676) q[3];
sx q[3];
rz(-0.85897173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9873535) q[2];
sx q[2];
rz(-2.8581212) q[2];
sx q[2];
rz(-3.0563291) q[2];
rz(1.2003468) q[3];
sx q[3];
rz(-1.6253358) q[3];
sx q[3];
rz(-0.15032642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7744556) q[0];
sx q[0];
rz(-0.13639233) q[0];
sx q[0];
rz(0.95463395) q[0];
rz(-0.57149354) q[1];
sx q[1];
rz(-1.9479472) q[1];
sx q[1];
rz(0.25209299) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82089059) q[0];
sx q[0];
rz(-1.1517236) q[0];
sx q[0];
rz(-1.7295895) q[0];
rz(-0.7438296) q[2];
sx q[2];
rz(-0.93223909) q[2];
sx q[2];
rz(-3.0628052) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8512307) q[1];
sx q[1];
rz(-2.2154659) q[1];
sx q[1];
rz(-2.6840997) q[1];
x q[2];
rz(0.37766405) q[3];
sx q[3];
rz(-2.0026243) q[3];
sx q[3];
rz(0.8197195) q[3];
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
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7082108) q[0];
sx q[0];
rz(-0.0032341783) q[0];
sx q[0];
rz(-1.4138387) q[0];
rz(-0.66043234) q[1];
sx q[1];
rz(-1.753189) q[1];
sx q[1];
rz(-1.7699014) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9773974) q[0];
sx q[0];
rz(-1.482555) q[0];
sx q[0];
rz(0.95055687) q[0];
rz(-1.781109) q[2];
sx q[2];
rz(-0.82092199) q[2];
sx q[2];
rz(-2.6471777) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.129442) q[1];
sx q[1];
rz(-0.65629849) q[1];
sx q[1];
rz(-0.49883573) q[1];
x q[2];
rz(0.38512226) q[3];
sx q[3];
rz(-2.754854) q[3];
sx q[3];
rz(2.573607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3147605) q[2];
sx q[2];
rz(-2.5173352) q[2];
sx q[2];
rz(-2.7436942) q[2];
rz(2.6509616) q[3];
sx q[3];
rz(-1.5964973) q[3];
sx q[3];
rz(1.8060961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(0.21815498) q[0];
sx q[0];
rz(-1.2175918) q[0];
sx q[0];
rz(0.23432215) q[0];
rz(-1.461347) q[1];
sx q[1];
rz(-0.82273465) q[1];
sx q[1];
rz(0.27059069) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10588564) q[0];
sx q[0];
rz(-2.8286655) q[0];
sx q[0];
rz(2.3653415) q[0];
rz(-2.0653535) q[2];
sx q[2];
rz(-2.7521172) q[2];
sx q[2];
rz(1.4091834) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2724243) q[1];
sx q[1];
rz(-0.13726343) q[1];
sx q[1];
rz(3.1174201) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.10176615) q[3];
sx q[3];
rz(-1.093285) q[3];
sx q[3];
rz(-2.4448642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.33971912) q[2];
sx q[2];
rz(-3.0568213) q[2];
sx q[2];
rz(-1.6516997) q[2];
rz(-0.70288944) q[3];
sx q[3];
rz(-1.9873762) q[3];
sx q[3];
rz(-1.9809013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7372195) q[0];
sx q[0];
rz(-1.5727366) q[0];
sx q[0];
rz(0.15429601) q[0];
rz(-0.94296304) q[1];
sx q[1];
rz(-1.1318726) q[1];
sx q[1];
rz(-2.399209) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5267747) q[0];
sx q[0];
rz(-2.1702538) q[0];
sx q[0];
rz(0.74032289) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1701199) q[2];
sx q[2];
rz(-0.68442548) q[2];
sx q[2];
rz(-1.3542092) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.87634516) q[1];
sx q[1];
rz(-2.7956388) q[1];
sx q[1];
rz(0.032851263) q[1];
x q[2];
rz(2.9418482) q[3];
sx q[3];
rz(-2.9962857) q[3];
sx q[3];
rz(-2.6655281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2877038) q[2];
sx q[2];
rz(-2.0159857) q[2];
sx q[2];
rz(-1.5596191) q[2];
rz(0.21073267) q[3];
sx q[3];
rz(-2.3570574) q[3];
sx q[3];
rz(0.5334841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8515274) q[0];
sx q[0];
rz(-2.1061438) q[0];
sx q[0];
rz(0.60594546) q[0];
rz(1.3416946) q[1];
sx q[1];
rz(-1.1732027) q[1];
sx q[1];
rz(-1.8684594) q[1];
rz(2.3201597) q[2];
sx q[2];
rz(-1.3315622) q[2];
sx q[2];
rz(-0.17377725) q[2];
rz(-3.001694) q[3];
sx q[3];
rz(-1.4362873) q[3];
sx q[3];
rz(0.39713058) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
