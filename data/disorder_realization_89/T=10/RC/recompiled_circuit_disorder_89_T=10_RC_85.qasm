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
rz(-2.7198305) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.598782) q[0];
sx q[0];
rz(-3.071975) q[0];
sx q[0];
rz(1.9570062) q[0];
rz(-pi) q[1];
x q[1];
rz(0.68658457) q[2];
sx q[2];
rz(-1.6010487) q[2];
sx q[2];
rz(-0.1498915) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.64695839) q[1];
sx q[1];
rz(-1.1263493) q[1];
sx q[1];
rz(1.3691982) q[1];
rz(-pi) q[2];
rz(-1.7105605) q[3];
sx q[3];
rz(-2.1602109) q[3];
sx q[3];
rz(-2.0624269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0191779) q[2];
sx q[2];
rz(-2.3712967) q[2];
sx q[2];
rz(-1.0100693) q[2];
rz(1.4953556) q[3];
sx q[3];
rz(-1.3388747) q[3];
sx q[3];
rz(3*pi/11) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0881969) q[0];
sx q[0];
rz(-1.2258376) q[0];
sx q[0];
rz(2.2825867) q[0];
rz(-2.7711218) q[1];
sx q[1];
rz(-1.5444688) q[1];
sx q[1];
rz(1.6765615) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5252936) q[0];
sx q[0];
rz(-1.787961) q[0];
sx q[0];
rz(2.8263682) q[0];
x q[1];
rz(-2.6066166) q[2];
sx q[2];
rz(-1.2777395) q[2];
sx q[2];
rz(1.5154293) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6795923) q[1];
sx q[1];
rz(-2.3604217) q[1];
sx q[1];
rz(-0.52266927) q[1];
rz(-pi) q[2];
rz(0.46918418) q[3];
sx q[3];
rz(-2.5622534) q[3];
sx q[3];
rz(2.7964696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7039965) q[2];
sx q[2];
rz(-1.9251172) q[2];
sx q[2];
rz(0.55830467) q[2];
rz(-2.1022508) q[3];
sx q[3];
rz(-1.2149518) q[3];
sx q[3];
rz(0.59282747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52477437) q[0];
sx q[0];
rz(-2.1891948) q[0];
sx q[0];
rz(0.16361374) q[0];
rz(-2.4257461) q[1];
sx q[1];
rz(-2.2952081) q[1];
sx q[1];
rz(-1.3735501) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87953075) q[0];
sx q[0];
rz(-0.52723215) q[0];
sx q[0];
rz(0.049588163) q[0];
rz(-1.9203556) q[2];
sx q[2];
rz(-0.74410838) q[2];
sx q[2];
rz(1.0457872) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.19792412) q[1];
sx q[1];
rz(-1.0218044) q[1];
sx q[1];
rz(-2.3180069) q[1];
rz(-pi) q[2];
rz(2.3397066) q[3];
sx q[3];
rz(-0.31517866) q[3];
sx q[3];
rz(1.2847881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.13087656) q[2];
sx q[2];
rz(-0.54543442) q[2];
sx q[2];
rz(-1.019657) q[2];
rz(2.1608458) q[3];
sx q[3];
rz(-0.5766944) q[3];
sx q[3];
rz(2.5286634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(2.2402128) q[0];
sx q[0];
rz(-2.4432683) q[0];
sx q[0];
rz(0.83909488) q[0];
rz(1.6038731) q[1];
sx q[1];
rz(-1.6042177) q[1];
sx q[1];
rz(2.531321) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3195517) q[0];
sx q[0];
rz(-0.10986957) q[0];
sx q[0];
rz(1.4676453) q[0];
rz(-pi) q[1];
rz(2.3158014) q[2];
sx q[2];
rz(-2.6399887) q[2];
sx q[2];
rz(0.85130168) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3210541) q[1];
sx q[1];
rz(-2.39403) q[1];
sx q[1];
rz(1.9588406) q[1];
x q[2];
rz(-2.1498333) q[3];
sx q[3];
rz(-1.1513396) q[3];
sx q[3];
rz(-1.6384244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.443976) q[2];
sx q[2];
rz(-0.51091754) q[2];
sx q[2];
rz(2.9042517) q[2];
rz(-0.90732968) q[3];
sx q[3];
rz(-1.5483587) q[3];
sx q[3];
rz(1.2835519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5089371) q[0];
sx q[0];
rz(-3.0314358) q[0];
sx q[0];
rz(-1.5326112) q[0];
rz(-0.73348796) q[1];
sx q[1];
rz(-1.8746904) q[1];
sx q[1];
rz(-1.8283432) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8719296) q[0];
sx q[0];
rz(-0.73686826) q[0];
sx q[0];
rz(2.065573) q[0];
rz(-pi) q[1];
rz(-0.03582844) q[2];
sx q[2];
rz(-1.6275068) q[2];
sx q[2];
rz(2.5113475) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9338867) q[1];
sx q[1];
rz(-1.943214) q[1];
sx q[1];
rz(1.1697342) q[1];
x q[2];
rz(0.87923268) q[3];
sx q[3];
rz(-2.0428265) q[3];
sx q[3];
rz(-1.8902682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0017172) q[2];
sx q[2];
rz(-1.2728609) q[2];
sx q[2];
rz(-0.66429794) q[2];
rz(-2.4364566) q[3];
sx q[3];
rz(-2.4987529) q[3];
sx q[3];
rz(0.10372182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.076684549) q[0];
sx q[0];
rz(-0.10184558) q[0];
sx q[0];
rz(3.0867807) q[0];
rz(2.1272155) q[1];
sx q[1];
rz(-1.6631815) q[1];
sx q[1];
rz(-2.6584113) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85855243) q[0];
sx q[0];
rz(-0.1917834) q[0];
sx q[0];
rz(-1.900308) q[0];
rz(1.6364355) q[2];
sx q[2];
rz(-0.46354957) q[2];
sx q[2];
rz(2.133873) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.74181611) q[1];
sx q[1];
rz(-2.2175334) q[1];
sx q[1];
rz(0.6439376) q[1];
x q[2];
rz(1.8910965) q[3];
sx q[3];
rz(-0.34892504) q[3];
sx q[3];
rz(2.2826209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.15423916) q[2];
sx q[2];
rz(-2.8581212) q[2];
sx q[2];
rz(-0.085263578) q[2];
rz(1.2003468) q[3];
sx q[3];
rz(-1.6253358) q[3];
sx q[3];
rz(-0.15032642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7744556) q[0];
sx q[0];
rz(-3.0052003) q[0];
sx q[0];
rz(-2.1869587) q[0];
rz(0.57149354) q[1];
sx q[1];
rz(-1.9479472) q[1];
sx q[1];
rz(2.8894997) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82089059) q[0];
sx q[0];
rz(-1.989869) q[0];
sx q[0];
rz(-1.7295895) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3977631) q[2];
sx q[2];
rz(-2.2093536) q[2];
sx q[2];
rz(-0.078787412) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8512307) q[1];
sx q[1];
rz(-2.2154659) q[1];
sx q[1];
rz(2.6840997) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0310568) q[3];
sx q[3];
rz(-1.9122951) q[3];
sx q[3];
rz(2.5550492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.38348848) q[2];
sx q[2];
rz(-2.1540116) q[2];
sx q[2];
rz(2.2371116) q[2];
rz(-2.1379743) q[3];
sx q[3];
rz(-0.86528722) q[3];
sx q[3];
rz(-2.482567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43338183) q[0];
sx q[0];
rz(-0.0032341783) q[0];
sx q[0];
rz(1.4138387) q[0];
rz(-2.4811603) q[1];
sx q[1];
rz(-1.3884037) q[1];
sx q[1];
rz(1.3716912) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8577514) q[0];
sx q[0];
rz(-2.5159266) q[0];
sx q[0];
rz(-1.7218504) q[0];
rz(-pi) q[1];
x q[1];
rz(0.22050627) q[2];
sx q[2];
rz(-0.77324152) q[2];
sx q[2];
rz(-2.3436433) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9904069) q[1];
sx q[1];
rz(-1.2745665) q[1];
sx q[1];
rz(0.59466655) q[1];
rz(2.7564704) q[3];
sx q[3];
rz(-2.754854) q[3];
sx q[3];
rz(-2.573607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.82683212) q[2];
sx q[2];
rz(-2.5173352) q[2];
sx q[2];
rz(-2.7436942) q[2];
rz(-2.6509616) q[3];
sx q[3];
rz(-1.5450954) q[3];
sx q[3];
rz(1.8060961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21815498) q[0];
sx q[0];
rz(-1.9240009) q[0];
sx q[0];
rz(-0.23432215) q[0];
rz(-1.6802457) q[1];
sx q[1];
rz(-0.82273465) q[1];
sx q[1];
rz(-0.27059069) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4280677) q[0];
sx q[0];
rz(-1.7881835) q[0];
sx q[0];
rz(2.9146951) q[0];
x q[1];
rz(-0.19240304) q[2];
sx q[2];
rz(-1.9115598) q[2];
sx q[2];
rz(-0.88142384) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2968263) q[1];
sx q[1];
rz(-1.7080194) q[1];
sx q[1];
rz(1.574135) q[1];
rz(-pi) q[2];
rz(1.3769334) q[3];
sx q[3];
rz(-0.48741515) q[3];
sx q[3];
rz(2.2262115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8018735) q[2];
sx q[2];
rz(-0.084771307) q[2];
sx q[2];
rz(1.489893) q[2];
rz(0.70288944) q[3];
sx q[3];
rz(-1.9873762) q[3];
sx q[3];
rz(-1.1606914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7372195) q[0];
sx q[0];
rz(-1.5727366) q[0];
sx q[0];
rz(-0.15429601) q[0];
rz(2.1986296) q[1];
sx q[1];
rz(-1.1318726) q[1];
sx q[1];
rz(0.74238366) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61481793) q[0];
sx q[0];
rz(-0.97133884) q[0];
sx q[0];
rz(0.74032289) q[0];
rz(2.1701199) q[2];
sx q[2];
rz(-2.4571672) q[2];
sx q[2];
rz(1.3542092) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2303289) q[1];
sx q[1];
rz(-1.2250369) q[1];
sx q[1];
rz(-1.5589578) q[1];
x q[2];
rz(1.5417682) q[3];
sx q[3];
rz(-1.4283984) q[3];
sx q[3];
rz(-2.8673429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8538889) q[2];
sx q[2];
rz(-1.1256069) q[2];
sx q[2];
rz(1.5819736) q[2];
rz(-2.93086) q[3];
sx q[3];
rz(-2.3570574) q[3];
sx q[3];
rz(0.5334841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8515274) q[0];
sx q[0];
rz(-1.0354488) q[0];
sx q[0];
rz(-2.5356472) q[0];
rz(1.3416946) q[1];
sx q[1];
rz(-1.1732027) q[1];
sx q[1];
rz(-1.8684594) q[1];
rz(-0.32158357) q[2];
sx q[2];
rz(-0.84761534) q[2];
sx q[2];
rz(-1.9615704) q[2];
rz(2.3711754) q[3];
sx q[3];
rz(-2.9478248) q[3];
sx q[3];
rz(-0.41268681) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
