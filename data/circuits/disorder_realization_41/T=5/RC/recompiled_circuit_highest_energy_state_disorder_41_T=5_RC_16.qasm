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
rz(2.3693585) q[0];
sx q[0];
rz(-1.9263664) q[0];
sx q[0];
rz(1.1507432) q[0];
rz(2.8590705) q[1];
sx q[1];
rz(-1.1257659) q[1];
sx q[1];
rz(-1.1123302) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38508666) q[0];
sx q[0];
rz(-1.6876564) q[0];
sx q[0];
rz(0.26330848) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0289828) q[2];
sx q[2];
rz(-0.99163429) q[2];
sx q[2];
rz(0.65161588) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2124651) q[1];
sx q[1];
rz(-1.0914789) q[1];
sx q[1];
rz(-2.1838837) q[1];
x q[2];
rz(-0.69975515) q[3];
sx q[3];
rz(-2.2198921) q[3];
sx q[3];
rz(3.0720667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.34899601) q[2];
sx q[2];
rz(-2.0670321) q[2];
sx q[2];
rz(1.7462771) q[2];
rz(-0.058517728) q[3];
sx q[3];
rz(-1.4165712) q[3];
sx q[3];
rz(-1.3099366) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7175452) q[0];
sx q[0];
rz(-0.40322867) q[0];
sx q[0];
rz(0.44977093) q[0];
rz(-2.6768118) q[1];
sx q[1];
rz(-1.8315146) q[1];
sx q[1];
rz(-0.25258499) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9802482) q[0];
sx q[0];
rz(-0.35156116) q[0];
sx q[0];
rz(1.6173961) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4812883) q[2];
sx q[2];
rz(-1.2129158) q[2];
sx q[2];
rz(-1.344327) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.77923939) q[1];
sx q[1];
rz(-1.5736516) q[1];
sx q[1];
rz(1.8523907) q[1];
x q[2];
rz(-2.1869867) q[3];
sx q[3];
rz(-0.86087117) q[3];
sx q[3];
rz(1.4675642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.8978591) q[2];
sx q[2];
rz(-0.89577883) q[2];
sx q[2];
rz(-0.016156999) q[2];
rz(0.92205087) q[3];
sx q[3];
rz(-2.1478896) q[3];
sx q[3];
rz(3.0097358) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.09318) q[0];
sx q[0];
rz(-1.484363) q[0];
sx q[0];
rz(-2.4600001) q[0];
rz(2.5438578) q[1];
sx q[1];
rz(-1.455541) q[1];
sx q[1];
rz(2.8173503) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63211381) q[0];
sx q[0];
rz(-1.1369708) q[0];
sx q[0];
rz(1.9471517) q[0];
rz(-1.6382255) q[2];
sx q[2];
rz(-1.0649657) q[2];
sx q[2];
rz(2.1492599) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5432918) q[1];
sx q[1];
rz(-0.17379119) q[1];
sx q[1];
rz(-1.8164746) q[1];
x q[2];
rz(-2.1799654) q[3];
sx q[3];
rz(-1.7097843) q[3];
sx q[3];
rz(-2.9442996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2250259) q[2];
sx q[2];
rz(-0.39174199) q[2];
sx q[2];
rz(-1.0507091) q[2];
rz(0.74058908) q[3];
sx q[3];
rz(-1.8480443) q[3];
sx q[3];
rz(0.061323969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4968313) q[0];
sx q[0];
rz(-2.3065541) q[0];
sx q[0];
rz(2.187425) q[0];
rz(-2.0046115) q[1];
sx q[1];
rz(-1.6258207) q[1];
sx q[1];
rz(2.383393) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6369259) q[0];
sx q[0];
rz(-1.530679) q[0];
sx q[0];
rz(1.841943) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7654645) q[2];
sx q[2];
rz(-0.83651453) q[2];
sx q[2];
rz(1.0502953) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0688732) q[1];
sx q[1];
rz(-2.2617635) q[1];
sx q[1];
rz(3.0145538) q[1];
rz(-pi) q[2];
x q[2];
rz(0.011298374) q[3];
sx q[3];
rz(-0.73888981) q[3];
sx q[3];
rz(-2.9113462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.053146426) q[2];
sx q[2];
rz(-0.55067486) q[2];
sx q[2];
rz(1.3147563) q[2];
rz(2.8407319) q[3];
sx q[3];
rz(-1.8010537) q[3];
sx q[3];
rz(-0.68527591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9984197) q[0];
sx q[0];
rz(-1.5851861) q[0];
sx q[0];
rz(-1.6402624) q[0];
rz(-2.7550664) q[1];
sx q[1];
rz(-2.0730348) q[1];
sx q[1];
rz(1.5466669) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3980557) q[0];
sx q[0];
rz(-0.72303444) q[0];
sx q[0];
rz(2.8788752) q[0];
rz(-0.32471409) q[2];
sx q[2];
rz(-2.1544059) q[2];
sx q[2];
rz(-1.8167855) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3749126) q[1];
sx q[1];
rz(-1.1385575) q[1];
sx q[1];
rz(-1.47827) q[1];
x q[2];
rz(-0.57379471) q[3];
sx q[3];
rz(-1.6125896) q[3];
sx q[3];
rz(-2.6190663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3950562) q[2];
sx q[2];
rz(-0.48607963) q[2];
sx q[2];
rz(-2.0716095) q[2];
rz(-0.88548958) q[3];
sx q[3];
rz(-2.077379) q[3];
sx q[3];
rz(-1.989367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92538658) q[0];
sx q[0];
rz(-2.8031271) q[0];
sx q[0];
rz(2.0881407) q[0];
rz(-2.9523051) q[1];
sx q[1];
rz(-2.3129251) q[1];
sx q[1];
rz(1.4922356) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5439592) q[0];
sx q[0];
rz(-1.2043556) q[0];
sx q[0];
rz(2.9935915) q[0];
rz(-pi) q[1];
rz(2.9290861) q[2];
sx q[2];
rz(-1.2709317) q[2];
sx q[2];
rz(2.1455255) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.63000667) q[1];
sx q[1];
rz(-1.6701856) q[1];
sx q[1];
rz(-2.6376827) q[1];
rz(-pi) q[2];
rz(-1.6118649) q[3];
sx q[3];
rz(-2.3784935) q[3];
sx q[3];
rz(-1.7880609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.80387992) q[2];
sx q[2];
rz(-2.1264075) q[2];
sx q[2];
rz(2.0886776) q[2];
rz(3.0247011) q[3];
sx q[3];
rz(-2.0785073) q[3];
sx q[3];
rz(-1.6803928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18148024) q[0];
sx q[0];
rz(-2.5616665) q[0];
sx q[0];
rz(-0.54452288) q[0];
rz(-2.8879884) q[1];
sx q[1];
rz(-0.2519775) q[1];
sx q[1];
rz(2.8869218) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7527251) q[0];
sx q[0];
rz(-1.6427186) q[0];
sx q[0];
rz(3.1353463) q[0];
rz(-pi) q[1];
rz(-0.20607941) q[2];
sx q[2];
rz(-0.40073088) q[2];
sx q[2];
rz(1.8054659) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0905625) q[1];
sx q[1];
rz(-1.5080308) q[1];
sx q[1];
rz(2.323137) q[1];
x q[2];
rz(0.56470498) q[3];
sx q[3];
rz(-2.2484591) q[3];
sx q[3];
rz(-1.1260179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9165245) q[2];
sx q[2];
rz(-2.2038286) q[2];
sx q[2];
rz(2.7226105) q[2];
rz(0.032912832) q[3];
sx q[3];
rz(-1.2337647) q[3];
sx q[3];
rz(-2.029443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5527749) q[0];
sx q[0];
rz(-2.3989615) q[0];
sx q[0];
rz(1.3277998) q[0];
rz(1.0325507) q[1];
sx q[1];
rz(-1.3464709) q[1];
sx q[1];
rz(0.58951497) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96846554) q[0];
sx q[0];
rz(-1.0778946) q[0];
sx q[0];
rz(-2.1987134) q[0];
rz(-pi) q[1];
rz(2.8500798) q[2];
sx q[2];
rz(-1.8290724) q[2];
sx q[2];
rz(-1.2746537) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2736328) q[1];
sx q[1];
rz(-0.74791779) q[1];
sx q[1];
rz(-0.91435097) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0716075) q[3];
sx q[3];
rz(-1.2195865) q[3];
sx q[3];
rz(-2.2413188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0699658) q[2];
sx q[2];
rz(-0.67487851) q[2];
sx q[2];
rz(-0.91342941) q[2];
rz(2.6895798) q[3];
sx q[3];
rz(-1.858523) q[3];
sx q[3];
rz(1.7970596) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87208676) q[0];
sx q[0];
rz(-2.1156613) q[0];
sx q[0];
rz(1.6312697) q[0];
rz(-1.3144846) q[1];
sx q[1];
rz(-2.1035106) q[1];
sx q[1];
rz(-2.7297535) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5740032) q[0];
sx q[0];
rz(-1.8471944) q[0];
sx q[0];
rz(0.92692356) q[0];
rz(-pi) q[1];
rz(2.9963295) q[2];
sx q[2];
rz(-2.8194339) q[2];
sx q[2];
rz(-1.3346498) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.50166124) q[1];
sx q[1];
rz(-1.2218214) q[1];
sx q[1];
rz(-2.5759012) q[1];
x q[2];
rz(-1.8985073) q[3];
sx q[3];
rz(-1.9372889) q[3];
sx q[3];
rz(-1.9043363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1371548) q[2];
sx q[2];
rz(-1.4919446) q[2];
sx q[2];
rz(-0.23492661) q[2];
rz(-2.2233985) q[3];
sx q[3];
rz(-2.4166959) q[3];
sx q[3];
rz(1.0880067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.039577) q[0];
sx q[0];
rz(-0.35440847) q[0];
sx q[0];
rz(1.8812195) q[0];
rz(-1.4541939) q[1];
sx q[1];
rz(-1.4581542) q[1];
sx q[1];
rz(-1.4448602) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4871976) q[0];
sx q[0];
rz(-0.62614589) q[0];
sx q[0];
rz(2.3930413) q[0];
x q[1];
rz(-2.6131975) q[2];
sx q[2];
rz(-0.95167347) q[2];
sx q[2];
rz(-2.2592659) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6548938) q[1];
sx q[1];
rz(-1.4461262) q[1];
sx q[1];
rz(-0.71790865) q[1];
rz(-pi) q[2];
rz(-2.4244196) q[3];
sx q[3];
rz(-2.6568718) q[3];
sx q[3];
rz(-0.53418272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5789648) q[2];
sx q[2];
rz(-2.8426888) q[2];
sx q[2];
rz(0.12358269) q[2];
rz(-0.28892162) q[3];
sx q[3];
rz(-1.2762504) q[3];
sx q[3];
rz(-2.6821274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73744437) q[0];
sx q[0];
rz(-1.7973719) q[0];
sx q[0];
rz(1.6112882) q[0];
rz(0.96724802) q[1];
sx q[1];
rz(-2.1687242) q[1];
sx q[1];
rz(0.8868934) q[1];
rz(1.5161472) q[2];
sx q[2];
rz(-2.1693632) q[2];
sx q[2];
rz(-2.7827435) q[2];
rz(2.2933949) q[3];
sx q[3];
rz(-2.316939) q[3];
sx q[3];
rz(1.4311781) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
