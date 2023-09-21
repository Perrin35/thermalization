OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.80823094) q[0];
sx q[0];
rz(-2.1699177) q[0];
sx q[0];
rz(1.4810286) q[0];
rz(-0.12227585) q[1];
sx q[1];
rz(3.0552157) q[1];
sx q[1];
rz(9.4019158) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.011933) q[0];
sx q[0];
rz(-1.4225905) q[0];
sx q[0];
rz(-3.0236142) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0055662) q[2];
sx q[2];
rz(-2.044894) q[2];
sx q[2];
rz(-0.20027645) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3373128) q[1];
sx q[1];
rz(-2.7645281) q[1];
sx q[1];
rz(1.1537329) q[1];
rz(-pi) q[2];
rz(1.919235) q[3];
sx q[3];
rz(-0.64823965) q[3];
sx q[3];
rz(-3.0188308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7227398) q[2];
sx q[2];
rz(-0.8272233) q[2];
sx q[2];
rz(-0.76618761) q[2];
rz(-0.17151672) q[3];
sx q[3];
rz(-0.73232108) q[3];
sx q[3];
rz(-2.5550301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3294753) q[0];
sx q[0];
rz(-2.7024039) q[0];
sx q[0];
rz(-3.0549333) q[0];
rz(-0.86241972) q[1];
sx q[1];
rz(-2.5984867) q[1];
sx q[1];
rz(3.0564953) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9239685) q[0];
sx q[0];
rz(-0.94857615) q[0];
sx q[0];
rz(-1.5909965) q[0];
x q[1];
rz(2.000196) q[2];
sx q[2];
rz(-2.4944802) q[2];
sx q[2];
rz(2.5246758) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9566006) q[1];
sx q[1];
rz(-1.6829832) q[1];
sx q[1];
rz(0.086602028) q[1];
rz(-0.3195023) q[3];
sx q[3];
rz(-2.0204244) q[3];
sx q[3];
rz(-3.086123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8791085) q[2];
sx q[2];
rz(-1.5905249) q[2];
sx q[2];
rz(1.4651728) q[2];
rz(1.9654467) q[3];
sx q[3];
rz(-2.2361103) q[3];
sx q[3];
rz(2.8951077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48433205) q[0];
sx q[0];
rz(-3.0069139) q[0];
sx q[0];
rz(-2.235967) q[0];
rz(-2.8088645) q[1];
sx q[1];
rz(-0.85488027) q[1];
sx q[1];
rz(1.7571626) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15554409) q[0];
sx q[0];
rz(-1.4882003) q[0];
sx q[0];
rz(-2.6655469) q[0];
rz(-pi) q[1];
rz(-2.7483814) q[2];
sx q[2];
rz(-2.0438497) q[2];
sx q[2];
rz(0.87393239) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.78450173) q[1];
sx q[1];
rz(-0.75575268) q[1];
sx q[1];
rz(-3.109039) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0043886) q[3];
sx q[3];
rz(-2.0319788) q[3];
sx q[3];
rz(-3.0571292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.24421346) q[2];
sx q[2];
rz(-0.34003568) q[2];
sx q[2];
rz(1.3282233) q[2];
rz(-1.3978847) q[3];
sx q[3];
rz(-0.54026794) q[3];
sx q[3];
rz(0.11169294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.06552799) q[0];
sx q[0];
rz(-2.9632443) q[0];
sx q[0];
rz(-0.67101014) q[0];
rz(-1.3440075) q[1];
sx q[1];
rz(-1.9799415) q[1];
sx q[1];
rz(2.9342594) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3509388) q[0];
sx q[0];
rz(-1.6121943) q[0];
sx q[0];
rz(-0.008965094) q[0];
rz(-1.7447128) q[2];
sx q[2];
rz(-1.6007649) q[2];
sx q[2];
rz(2.7001691) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3934717) q[1];
sx q[1];
rz(-1.3337743) q[1];
sx q[1];
rz(2.9841828) q[1];
rz(3.1158833) q[3];
sx q[3];
rz(-2.6298012) q[3];
sx q[3];
rz(-2.5185086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8308782) q[2];
sx q[2];
rz(-2.0534616) q[2];
sx q[2];
rz(-0.80835289) q[2];
rz(-2.3338142) q[3];
sx q[3];
rz(-2.419796) q[3];
sx q[3];
rz(2.9479153) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59072524) q[0];
sx q[0];
rz(-2.4243675) q[0];
sx q[0];
rz(0.89548683) q[0];
rz(2.8764309) q[1];
sx q[1];
rz(-2.3064955) q[1];
sx q[1];
rz(-0.53363824) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69913188) q[0];
sx q[0];
rz(-1.2958382) q[0];
sx q[0];
rz(1.86637) q[0];
rz(-0.76782121) q[2];
sx q[2];
rz(-1.7349093) q[2];
sx q[2];
rz(0.62190157) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.071760885) q[1];
sx q[1];
rz(-1.3511718) q[1];
sx q[1];
rz(-1.4211618) q[1];
rz(1.8495464) q[3];
sx q[3];
rz(-2.4933443) q[3];
sx q[3];
rz(1.3384782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2395893) q[2];
sx q[2];
rz(-1.5633554) q[2];
sx q[2];
rz(-0.56838244) q[2];
rz(1.2166294) q[3];
sx q[3];
rz(-0.47074461) q[3];
sx q[3];
rz(0.34354982) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6333273) q[0];
sx q[0];
rz(-0.87843043) q[0];
sx q[0];
rz(1.1761965) q[0];
rz(1.5856702) q[1];
sx q[1];
rz(-2.1891179) q[1];
sx q[1];
rz(-2.1369381) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4906824) q[0];
sx q[0];
rz(-2.6255529) q[0];
sx q[0];
rz(-0.74923058) q[0];
x q[1];
rz(0.19095687) q[2];
sx q[2];
rz(-2.7381884) q[2];
sx q[2];
rz(-1.0952589) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.1130044) q[1];
sx q[1];
rz(-0.60739809) q[1];
sx q[1];
rz(2.7007156) q[1];
x q[2];
rz(-1.5435013) q[3];
sx q[3];
rz(-0.41397646) q[3];
sx q[3];
rz(0.15347813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7375609) q[2];
sx q[2];
rz(-0.99016756) q[2];
sx q[2];
rz(-2.693434) q[2];
rz(0.29799497) q[3];
sx q[3];
rz(-0.69152504) q[3];
sx q[3];
rz(-0.34860778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77666831) q[0];
sx q[0];
rz(-1.3013327) q[0];
sx q[0];
rz(-2.4834494) q[0];
rz(1.2843885) q[1];
sx q[1];
rz(-0.45416608) q[1];
sx q[1];
rz(2.4694494) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1561688) q[0];
sx q[0];
rz(-2.9754313) q[0];
sx q[0];
rz(-1.6155433) q[0];
x q[1];
rz(-2.0124112) q[2];
sx q[2];
rz(-1.1814983) q[2];
sx q[2];
rz(2.0814975) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.044683177) q[1];
sx q[1];
rz(-1.9367366) q[1];
sx q[1];
rz(-2.6198322) q[1];
rz(-pi) q[2];
rz(2.4403205) q[3];
sx q[3];
rz(-1.4605195) q[3];
sx q[3];
rz(2.9993204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.60823524) q[2];
sx q[2];
rz(-0.63295263) q[2];
sx q[2];
rz(2.7427924) q[2];
rz(2.3170025) q[3];
sx q[3];
rz(-1.7493068) q[3];
sx q[3];
rz(2.5430172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-2.7010715) q[0];
sx q[0];
rz(-0.43061391) q[0];
sx q[0];
rz(-0.15821247) q[0];
rz(1.0868866) q[1];
sx q[1];
rz(-2.465076) q[1];
sx q[1];
rz(2.3349082) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0772484) q[0];
sx q[0];
rz(-1.4532538) q[0];
sx q[0];
rz(-2.006152) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4629455) q[2];
sx q[2];
rz(-1.3662405) q[2];
sx q[2];
rz(-2.4754935) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9629509) q[1];
sx q[1];
rz(-2.2370403) q[1];
sx q[1];
rz(-1.9402177) q[1];
rz(1.2699548) q[3];
sx q[3];
rz(-1.5402208) q[3];
sx q[3];
rz(-1.2310864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.57758254) q[2];
sx q[2];
rz(-1.6489886) q[2];
sx q[2];
rz(2.1419443) q[2];
rz(-3.0380761) q[3];
sx q[3];
rz(-3.0045356) q[3];
sx q[3];
rz(-2.0172393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0209811) q[0];
sx q[0];
rz(-0.7779026) q[0];
sx q[0];
rz(-2.4840684) q[0];
rz(0.28655562) q[1];
sx q[1];
rz(-0.92266881) q[1];
sx q[1];
rz(-0.34067571) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1239615) q[0];
sx q[0];
rz(-0.88273662) q[0];
sx q[0];
rz(2.7244363) q[0];
rz(2.1358143) q[2];
sx q[2];
rz(-1.7592332) q[2];
sx q[2];
rz(-0.68067683) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.162519) q[1];
sx q[1];
rz(-2.5538429) q[1];
sx q[1];
rz(1.2390922) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4021923) q[3];
sx q[3];
rz(-1.3753034) q[3];
sx q[3];
rz(-0.80381264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.75491536) q[2];
sx q[2];
rz(-1.9415386) q[2];
sx q[2];
rz(-0.49794751) q[2];
rz(-2.8399816) q[3];
sx q[3];
rz(-2.7409654) q[3];
sx q[3];
rz(0.51914674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4122445) q[0];
sx q[0];
rz(-3.0817741) q[0];
sx q[0];
rz(2.8630032) q[0];
rz(-0.57922286) q[1];
sx q[1];
rz(-0.93943739) q[1];
sx q[1];
rz(0.07671193) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.840832) q[0];
sx q[0];
rz(-1.4474086) q[0];
sx q[0];
rz(1.2542017) q[0];
rz(-pi) q[1];
rz(2.9701091) q[2];
sx q[2];
rz(-2.137261) q[2];
sx q[2];
rz(-0.4085853) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0753239) q[1];
sx q[1];
rz(-1.9362209) q[1];
sx q[1];
rz(-1.5931904) q[1];
x q[2];
rz(1.4823227) q[3];
sx q[3];
rz(-1.3946597) q[3];
sx q[3];
rz(2.4000771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3672093) q[2];
sx q[2];
rz(-2.4343906) q[2];
sx q[2];
rz(-0.45551604) q[2];
rz(-2.7101743) q[3];
sx q[3];
rz(-3.016267) q[3];
sx q[3];
rz(-3.07807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9615622) q[0];
sx q[0];
rz(-1.7407692) q[0];
sx q[0];
rz(-2.2041007) q[0];
rz(0.58615276) q[1];
sx q[1];
rz(-1.502232) q[1];
sx q[1];
rz(-1.4486817) q[1];
rz(0.2858328) q[2];
sx q[2];
rz(-1.1163859) q[2];
sx q[2];
rz(-3.1039539) q[2];
rz(-1.205411) q[3];
sx q[3];
rz(-1.1689417) q[3];
sx q[3];
rz(0.4369215) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];