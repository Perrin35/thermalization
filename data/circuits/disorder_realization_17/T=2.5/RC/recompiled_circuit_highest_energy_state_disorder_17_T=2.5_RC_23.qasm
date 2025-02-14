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
rz(-2.1276346) q[0];
sx q[0];
rz(-1.4528217) q[0];
sx q[0];
rz(0.59803522) q[0];
rz(1.1818089) q[1];
sx q[1];
rz(2.476517) q[1];
sx q[1];
rz(9.2718931) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6939815) q[0];
sx q[0];
rz(-1.6504399) q[0];
sx q[0];
rz(-1.4586186) q[0];
rz(-0.48798765) q[2];
sx q[2];
rz(-1.6225015) q[2];
sx q[2];
rz(2.0590797) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8064639) q[1];
sx q[1];
rz(-1.4452308) q[1];
sx q[1];
rz(-1.7837202) q[1];
rz(-0.045957248) q[3];
sx q[3];
rz(-1.5678239) q[3];
sx q[3];
rz(-1.7645022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.87024706) q[2];
sx q[2];
rz(-2.953244) q[2];
sx q[2];
rz(-1.1790454) q[2];
rz(-2.7315268) q[3];
sx q[3];
rz(-1.054801) q[3];
sx q[3];
rz(0.89455354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9533185) q[0];
sx q[0];
rz(-0.66903791) q[0];
sx q[0];
rz(-2.8642995) q[0];
rz(2.9412728) q[1];
sx q[1];
rz(-0.89626139) q[1];
sx q[1];
rz(-0.083981363) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0344211) q[0];
sx q[0];
rz(-2.0382529) q[0];
sx q[0];
rz(0.22308992) q[0];
x q[1];
rz(-0.72767392) q[2];
sx q[2];
rz(-1.4598338) q[2];
sx q[2];
rz(2.3425255) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0711956) q[1];
sx q[1];
rz(-2.4807811) q[1];
sx q[1];
rz(-0.7699716) q[1];
rz(-pi) q[2];
rz(0.53392729) q[3];
sx q[3];
rz(-1.3442497) q[3];
sx q[3];
rz(1.2372897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1818485) q[2];
sx q[2];
rz(-0.66755787) q[2];
sx q[2];
rz(-0.7461156) q[2];
rz(2.5203846) q[3];
sx q[3];
rz(-0.64380232) q[3];
sx q[3];
rz(-1.3889036) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8212432) q[0];
sx q[0];
rz(-0.72750434) q[0];
sx q[0];
rz(-1.9387091) q[0];
rz(-2.7299643) q[1];
sx q[1];
rz(-1.8609214) q[1];
sx q[1];
rz(-1.113755) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8804563) q[0];
sx q[0];
rz(-2.2703982) q[0];
sx q[0];
rz(-2.1300132) q[0];
rz(0.44373893) q[2];
sx q[2];
rz(-2.1390171) q[2];
sx q[2];
rz(-2.2310471) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.71761638) q[1];
sx q[1];
rz(-0.90002093) q[1];
sx q[1];
rz(2.9714279) q[1];
rz(-pi) q[2];
rz(2.5630234) q[3];
sx q[3];
rz(-2.1820543) q[3];
sx q[3];
rz(-1.3276154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0456298) q[2];
sx q[2];
rz(-0.7926422) q[2];
sx q[2];
rz(1.1620129) q[2];
rz(-0.80471936) q[3];
sx q[3];
rz(-1.8685124) q[3];
sx q[3];
rz(1.2812251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3941037) q[0];
sx q[0];
rz(-1.2825613) q[0];
sx q[0];
rz(1.5149186) q[0];
rz(2.6331242) q[1];
sx q[1];
rz(-2.4912806) q[1];
sx q[1];
rz(-1.5079087) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8430458) q[0];
sx q[0];
rz(-2.9905479) q[0];
sx q[0];
rz(1.5848978) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.14105395) q[2];
sx q[2];
rz(-1.3562437) q[2];
sx q[2];
rz(0.44520865) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0993938) q[1];
sx q[1];
rz(-1.458199) q[1];
sx q[1];
rz(1.2641175) q[1];
rz(-pi) q[2];
rz(0.91458264) q[3];
sx q[3];
rz(-1.9176502) q[3];
sx q[3];
rz(-0.28430249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7278829) q[2];
sx q[2];
rz(-1.9472313) q[2];
sx q[2];
rz(0.74753648) q[2];
rz(1.9109292) q[3];
sx q[3];
rz(-2.3321798) q[3];
sx q[3];
rz(-0.49720732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.704945) q[0];
sx q[0];
rz(-1.1282938) q[0];
sx q[0];
rz(1.6823912) q[0];
rz(-2.7875426) q[1];
sx q[1];
rz(-0.67146462) q[1];
sx q[1];
rz(-2.5669602) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8588206) q[0];
sx q[0];
rz(-0.078411113) q[0];
sx q[0];
rz(-0.17228244) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76894297) q[2];
sx q[2];
rz(-1.7684325) q[2];
sx q[2];
rz(-0.96425024) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.84031096) q[1];
sx q[1];
rz(-1.6224098) q[1];
sx q[1];
rz(-0.2272931) q[1];
rz(-pi) q[2];
rz(-0.1149807) q[3];
sx q[3];
rz(-1.3674187) q[3];
sx q[3];
rz(-2.9733244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6358801) q[2];
sx q[2];
rz(-1.1735703) q[2];
sx q[2];
rz(2.9393348) q[2];
rz(2.108719) q[3];
sx q[3];
rz(-0.60690108) q[3];
sx q[3];
rz(-3.0752944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.031484) q[0];
sx q[0];
rz(-0.38723543) q[0];
sx q[0];
rz(-0.37145823) q[0];
rz(0.29608852) q[1];
sx q[1];
rz(-1.5541872) q[1];
sx q[1];
rz(2.9782226) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5510493) q[0];
sx q[0];
rz(-1.3871017) q[0];
sx q[0];
rz(0.50409533) q[0];
rz(-0.54941515) q[2];
sx q[2];
rz(-1.9684458) q[2];
sx q[2];
rz(2.0009918) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0381546) q[1];
sx q[1];
rz(-1.1327599) q[1];
sx q[1];
rz(1.7704727) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9904731) q[3];
sx q[3];
rz(-2.2607779) q[3];
sx q[3];
rz(2.2595764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0801733) q[2];
sx q[2];
rz(-2.8714608) q[2];
sx q[2];
rz(0.6529676) q[2];
rz(0.4717007) q[3];
sx q[3];
rz(-0.74840778) q[3];
sx q[3];
rz(-2.4145943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2111874) q[0];
sx q[0];
rz(-2.0624332) q[0];
sx q[0];
rz(2.1042673) q[0];
rz(2.6264181) q[1];
sx q[1];
rz(-1.8637135) q[1];
sx q[1];
rz(1.3264664) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68149306) q[0];
sx q[0];
rz(-1.1837479) q[0];
sx q[0];
rz(-2.1429064) q[0];
rz(-pi) q[1];
rz(0.35197978) q[2];
sx q[2];
rz(-1.5328141) q[2];
sx q[2];
rz(-2.4427736) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4660037) q[1];
sx q[1];
rz(-1.6674622) q[1];
sx q[1];
rz(-2.3748891) q[1];
x q[2];
rz(0.063422261) q[3];
sx q[3];
rz(-1.8759955) q[3];
sx q[3];
rz(0.88808192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6933763) q[2];
sx q[2];
rz(-2.6770112) q[2];
sx q[2];
rz(0.32619897) q[2];
rz(0.97801963) q[3];
sx q[3];
rz(-1.2468485) q[3];
sx q[3];
rz(0.090156468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4433032) q[0];
sx q[0];
rz(-0.25743085) q[0];
sx q[0];
rz(2.8522016) q[0];
rz(-3.1293213) q[1];
sx q[1];
rz(-0.27996501) q[1];
sx q[1];
rz(1.4853005) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2606364) q[0];
sx q[0];
rz(-2.0800965) q[0];
sx q[0];
rz(-0.80672046) q[0];
rz(-pi) q[1];
rz(0.22356914) q[2];
sx q[2];
rz(-1.5180598) q[2];
sx q[2];
rz(-1.2745672) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0694794) q[1];
sx q[1];
rz(-2.8918307) q[1];
sx q[1];
rz(-0.62432557) q[1];
rz(-0.9910219) q[3];
sx q[3];
rz(-2.5977547) q[3];
sx q[3];
rz(2.4286634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.54899186) q[2];
sx q[2];
rz(-1.5718549) q[2];
sx q[2];
rz(-0.17295095) q[2];
rz(-1.7671827) q[3];
sx q[3];
rz(-1.7632615) q[3];
sx q[3];
rz(-1.4779199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7462815) q[0];
sx q[0];
rz(-2.6938541) q[0];
sx q[0];
rz(-0.81004274) q[0];
rz(0.41656247) q[1];
sx q[1];
rz(-0.77605334) q[1];
sx q[1];
rz(0.3449482) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5826915) q[0];
sx q[0];
rz(-1.4479637) q[0];
sx q[0];
rz(2.4305811) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9760804) q[2];
sx q[2];
rz(-1.8461746) q[2];
sx q[2];
rz(-2.2294758) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4629412) q[1];
sx q[1];
rz(-0.57616568) q[1];
sx q[1];
rz(0.94959308) q[1];
rz(-2.3446001) q[3];
sx q[3];
rz(-0.16599645) q[3];
sx q[3];
rz(2.5696563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.42744669) q[2];
sx q[2];
rz(-2.2553208) q[2];
sx q[2];
rz(-3.0575338) q[2];
rz(-2.2345624) q[3];
sx q[3];
rz(-1.7232938) q[3];
sx q[3];
rz(-0.9575873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8091938) q[0];
sx q[0];
rz(-0.087662307) q[0];
sx q[0];
rz(-2.8293389) q[0];
rz(-0.23278438) q[1];
sx q[1];
rz(-2.7504031) q[1];
sx q[1];
rz(1.8544474) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2256266) q[0];
sx q[0];
rz(-2.4082528) q[0];
sx q[0];
rz(-0.49247964) q[0];
rz(-pi) q[1];
rz(-0.63913362) q[2];
sx q[2];
rz(-1.5331563) q[2];
sx q[2];
rz(-0.89307484) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.171265) q[1];
sx q[1];
rz(-1.1210151) q[1];
sx q[1];
rz(0.27189092) q[1];
rz(-pi) q[2];
rz(1.200807) q[3];
sx q[3];
rz(-0.97856802) q[3];
sx q[3];
rz(-1.5371479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1708258) q[2];
sx q[2];
rz(-1.2831186) q[2];
sx q[2];
rz(1.6440561) q[2];
rz(1.0981285) q[3];
sx q[3];
rz(-0.69881717) q[3];
sx q[3];
rz(0.48925492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.004414) q[0];
sx q[0];
rz(-1.351384) q[0];
sx q[0];
rz(2.2030892) q[0];
rz(1.8347523) q[1];
sx q[1];
rz(-0.88941457) q[1];
sx q[1];
rz(-0.79759146) q[1];
rz(0.14456476) q[2];
sx q[2];
rz(-2.2546386) q[2];
sx q[2];
rz(-2.4972514) q[2];
rz(1.3908006) q[3];
sx q[3];
rz(-0.67435657) q[3];
sx q[3];
rz(-1.6410905) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
