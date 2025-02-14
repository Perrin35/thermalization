OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.0410864) q[0];
sx q[0];
rz(-0.20809986) q[0];
sx q[0];
rz(2.7383374) q[0];
rz(2.9085605) q[1];
sx q[1];
rz(-1.7014528) q[1];
sx q[1];
rz(-0.22388248) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3849523) q[0];
sx q[0];
rz(-2.1813574) q[0];
sx q[0];
rz(-1.8288307) q[0];
rz(-0.86647948) q[2];
sx q[2];
rz(-1.6687376) q[2];
sx q[2];
rz(-2.3549454) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.510281) q[1];
sx q[1];
rz(-1.3461539) q[1];
sx q[1];
rz(2.2366877) q[1];
rz(2.3305064) q[3];
sx q[3];
rz(-1.038365) q[3];
sx q[3];
rz(-2.6528751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.69748059) q[2];
sx q[2];
rz(-0.31284249) q[2];
sx q[2];
rz(-3.1057788) q[2];
rz(-0.73016417) q[3];
sx q[3];
rz(-1.9883479) q[3];
sx q[3];
rz(2.9353976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25074729) q[0];
sx q[0];
rz(-2.146281) q[0];
sx q[0];
rz(-2.9440951) q[0];
rz(-0.47710553) q[1];
sx q[1];
rz(-0.66480607) q[1];
sx q[1];
rz(2.3359931) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6960736) q[0];
sx q[0];
rz(-1.3851926) q[0];
sx q[0];
rz(2.473102) q[0];
rz(-pi) q[1];
rz(-2.7609772) q[2];
sx q[2];
rz(-1.632649) q[2];
sx q[2];
rz(0.46910367) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9220613) q[1];
sx q[1];
rz(-0.94486713) q[1];
sx q[1];
rz(0.74447592) q[1];
rz(-1.7116551) q[3];
sx q[3];
rz(-1.4211486) q[3];
sx q[3];
rz(0.7310673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8311367) q[2];
sx q[2];
rz(-1.6418991) q[2];
sx q[2];
rz(1.4031225) q[2];
rz(2.5095615) q[3];
sx q[3];
rz(-2.8681614) q[3];
sx q[3];
rz(3.0145751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3056575) q[0];
sx q[0];
rz(-0.62863612) q[0];
sx q[0];
rz(2.054731) q[0];
rz(0.19733363) q[1];
sx q[1];
rz(-1.5060164) q[1];
sx q[1];
rz(0.33263439) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63276382) q[0];
sx q[0];
rz(-2.6869505) q[0];
sx q[0];
rz(1.8931382) q[0];
x q[1];
rz(-0.83586043) q[2];
sx q[2];
rz(-1.5243013) q[2];
sx q[2];
rz(0.64620668) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8197934) q[1];
sx q[1];
rz(-1.395257) q[1];
sx q[1];
rz(-1.790253) q[1];
rz(-pi) q[2];
rz(-2.4147291) q[3];
sx q[3];
rz(-1.7439902) q[3];
sx q[3];
rz(1.7693335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7670333) q[2];
sx q[2];
rz(-1.550753) q[2];
sx q[2];
rz(-1.4683051) q[2];
rz(-3.1324006) q[3];
sx q[3];
rz(-0.46949783) q[3];
sx q[3];
rz(2.3495242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.511635) q[0];
sx q[0];
rz(-2.503643) q[0];
sx q[0];
rz(1.6988423) q[0];
rz(2.1171782) q[1];
sx q[1];
rz(-1.9099648) q[1];
sx q[1];
rz(2.3146497) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31931811) q[0];
sx q[0];
rz(-1.7546904) q[0];
sx q[0];
rz(-0.075550373) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9873136) q[2];
sx q[2];
rz(-0.88866975) q[2];
sx q[2];
rz(0.23274225) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0961541) q[1];
sx q[1];
rz(-1.3609582) q[1];
sx q[1];
rz(2.7010598) q[1];
x q[2];
rz(0.12110658) q[3];
sx q[3];
rz(-0.89205974) q[3];
sx q[3];
rz(0.27287441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3247165) q[2];
sx q[2];
rz(-2.3574895) q[2];
sx q[2];
rz(-1.8640222) q[2];
rz(2.4534524) q[3];
sx q[3];
rz(-1.0253996) q[3];
sx q[3];
rz(-0.96243206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0012896) q[0];
sx q[0];
rz(-1.4411417) q[0];
sx q[0];
rz(-0.21743123) q[0];
rz(0.28542074) q[1];
sx q[1];
rz(-2.5358584) q[1];
sx q[1];
rz(2.6434456) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.065154) q[0];
sx q[0];
rz(-2.3992043) q[0];
sx q[0];
rz(2.0641293) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1310546) q[2];
sx q[2];
rz(-0.56098191) q[2];
sx q[2];
rz(-1.9064685) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0464669) q[1];
sx q[1];
rz(-0.53494638) q[1];
sx q[1];
rz(-1.9880268) q[1];
rz(-1.823447) q[3];
sx q[3];
rz(-1.9985191) q[3];
sx q[3];
rz(1.778217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.51776) q[2];
sx q[2];
rz(-2.0443003) q[2];
sx q[2];
rz(-0.1499873) q[2];
rz(-0.013966694) q[3];
sx q[3];
rz(-1.59168) q[3];
sx q[3];
rz(-2.6278031) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49195313) q[0];
sx q[0];
rz(-2.5998901) q[0];
sx q[0];
rz(-0.81277043) q[0];
rz(1.4179519) q[1];
sx q[1];
rz(-1.5128472) q[1];
sx q[1];
rz(-2.0779804) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3115754) q[0];
sx q[0];
rz(-0.63236134) q[0];
sx q[0];
rz(2.564179) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9779149) q[2];
sx q[2];
rz(-1.5874169) q[2];
sx q[2];
rz(1.7134242) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8001682) q[1];
sx q[1];
rz(-1.582195) q[1];
sx q[1];
rz(-2.5428548) q[1];
rz(1.8251347) q[3];
sx q[3];
rz(-1.7539548) q[3];
sx q[3];
rz(1.1381799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.873988) q[2];
sx q[2];
rz(-2.5817817) q[2];
sx q[2];
rz(1.7248636) q[2];
rz(3.0806165) q[3];
sx q[3];
rz(-2.2305326) q[3];
sx q[3];
rz(-1.4643668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7199719) q[0];
sx q[0];
rz(-0.82141972) q[0];
sx q[0];
rz(0.11909568) q[0];
rz(-2.6630317) q[1];
sx q[1];
rz(-2.3275972) q[1];
sx q[1];
rz(0.68971577) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0499303) q[0];
sx q[0];
rz(-0.99794594) q[0];
sx q[0];
rz(2.2747103) q[0];
rz(-pi) q[1];
rz(-1.1999667) q[2];
sx q[2];
rz(-1.1111819) q[2];
sx q[2];
rz(2.4209765) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3773681) q[1];
sx q[1];
rz(-0.86478327) q[1];
sx q[1];
rz(3.0824667) q[1];
rz(1.4004565) q[3];
sx q[3];
rz(-0.719845) q[3];
sx q[3];
rz(3.0468536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0488284) q[2];
sx q[2];
rz(-3.0807639) q[2];
sx q[2];
rz(2.0737341) q[2];
rz(2.974406) q[3];
sx q[3];
rz(-1.7696295) q[3];
sx q[3];
rz(-0.13944496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6086513) q[0];
sx q[0];
rz(-0.81357384) q[0];
sx q[0];
rz(2.2331878) q[0];
rz(-0.99336973) q[1];
sx q[1];
rz(-0.16255957) q[1];
sx q[1];
rz(-2.2672674) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0197948) q[0];
sx q[0];
rz(-3.0082729) q[0];
sx q[0];
rz(-0.13102417) q[0];
rz(3.1221603) q[2];
sx q[2];
rz(-1.8284869) q[2];
sx q[2];
rz(-1.9382221) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.711124) q[1];
sx q[1];
rz(-1.5114771) q[1];
sx q[1];
rz(-1.435168) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.120058) q[3];
sx q[3];
rz(-2.5172699) q[3];
sx q[3];
rz(0.84824991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5270762) q[2];
sx q[2];
rz(-2.9824342) q[2];
sx q[2];
rz(2.2873774) q[2];
rz(-2.6863875) q[3];
sx q[3];
rz(-2.5762317) q[3];
sx q[3];
rz(2.8860886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.663317) q[0];
sx q[0];
rz(-1.9238967) q[0];
sx q[0];
rz(-1.664337) q[0];
rz(1.5746337) q[1];
sx q[1];
rz(-2.4576371) q[1];
sx q[1];
rz(0.73653594) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5631396) q[0];
sx q[0];
rz(-1.3361276) q[0];
sx q[0];
rz(-1.0893954) q[0];
x q[1];
rz(1.6746503) q[2];
sx q[2];
rz(-1.9496634) q[2];
sx q[2];
rz(-1.1844289) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8443201) q[1];
sx q[1];
rz(-1.8836918) q[1];
sx q[1];
rz(-1.827153) q[1];
x q[2];
rz(0.035321354) q[3];
sx q[3];
rz(-2.1675054) q[3];
sx q[3];
rz(-0.91340706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3972724) q[2];
sx q[2];
rz(-2.263415) q[2];
sx q[2];
rz(1.6142023) q[2];
rz(-2.7496036) q[3];
sx q[3];
rz(-1.9422928) q[3];
sx q[3];
rz(-0.049662445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.024260661) q[0];
sx q[0];
rz(-1.4735104) q[0];
sx q[0];
rz(2.9397553) q[0];
rz(0.2233389) q[1];
sx q[1];
rz(-0.52450648) q[1];
sx q[1];
rz(0.39252678) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.057923) q[0];
sx q[0];
rz(-1.6763078) q[0];
sx q[0];
rz(-1.1062463) q[0];
x q[1];
rz(3.0891339) q[2];
sx q[2];
rz(-0.83485583) q[2];
sx q[2];
rz(-2.3273766) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0333657) q[1];
sx q[1];
rz(-1.1443597) q[1];
sx q[1];
rz(0.99797499) q[1];
x q[2];
rz(-2.659117) q[3];
sx q[3];
rz(-1.8419208) q[3];
sx q[3];
rz(1.1678054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7212123) q[2];
sx q[2];
rz(-2.1537697) q[2];
sx q[2];
rz(-2.5610899) q[2];
rz(-2.765559) q[3];
sx q[3];
rz(-0.97085634) q[3];
sx q[3];
rz(-1.5413126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15039438) q[0];
sx q[0];
rz(-1.5491485) q[0];
sx q[0];
rz(-1.3843672) q[0];
rz(-1.6509542) q[1];
sx q[1];
rz(-1.2532267) q[1];
sx q[1];
rz(0.86029235) q[1];
rz(2.3357794) q[2];
sx q[2];
rz(-1.0849107) q[2];
sx q[2];
rz(2.8203865) q[2];
rz(-3.0720465) q[3];
sx q[3];
rz(-0.31252091) q[3];
sx q[3];
rz(-1.9557709) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
