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
rz(-2.1751997) q[0];
sx q[0];
rz(-1.3858495) q[0];
sx q[0];
rz(0.59544271) q[0];
rz(-1.8499941) q[1];
sx q[1];
rz(2.5505677) q[1];
sx q[1];
rz(12.443065) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5714471) q[0];
sx q[0];
rz(-0.71200221) q[0];
sx q[0];
rz(0.90306905) q[0];
rz(-pi) q[1];
rz(-2.3583636) q[2];
sx q[2];
rz(-2.1821664) q[2];
sx q[2];
rz(0.91125823) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9682546) q[1];
sx q[1];
rz(-0.23431331) q[1];
sx q[1];
rz(-1.987275) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7756157) q[3];
sx q[3];
rz(-1.4282188) q[3];
sx q[3];
rz(-0.78998427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3143602) q[2];
sx q[2];
rz(-2.0388956) q[2];
sx q[2];
rz(-0.11191351) q[2];
rz(1.3917475) q[3];
sx q[3];
rz(-1.9839957) q[3];
sx q[3];
rz(0.30645034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4157442) q[0];
sx q[0];
rz(-2.2968676) q[0];
sx q[0];
rz(0.14990212) q[0];
rz(0.28597486) q[1];
sx q[1];
rz(-1.7466702) q[1];
sx q[1];
rz(2.012595) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6647346) q[0];
sx q[0];
rz(-1.6010512) q[0];
sx q[0];
rz(-1.2824351) q[0];
x q[1];
rz(-1.5111792) q[2];
sx q[2];
rz(-1.5562415) q[2];
sx q[2];
rz(2.3090135) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.93176881) q[1];
sx q[1];
rz(-1.2126414) q[1];
sx q[1];
rz(1.7478554) q[1];
rz(-pi) q[2];
x q[2];
rz(0.58993323) q[3];
sx q[3];
rz(-2.1141756) q[3];
sx q[3];
rz(-0.50266279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4631606) q[2];
sx q[2];
rz(-2.2191935) q[2];
sx q[2];
rz(1.9909667) q[2];
rz(-1.9567418) q[3];
sx q[3];
rz(-1.7246282) q[3];
sx q[3];
rz(-3.1209893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6612369) q[0];
sx q[0];
rz(-2.895597) q[0];
sx q[0];
rz(-0.38159698) q[0];
rz(-1.2941788) q[1];
sx q[1];
rz(-1.1062063) q[1];
sx q[1];
rz(2.9433184) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37368942) q[0];
sx q[0];
rz(-1.3438017) q[0];
sx q[0];
rz(0.17902184) q[0];
x q[1];
rz(-3.0610721) q[2];
sx q[2];
rz(-2.2386907) q[2];
sx q[2];
rz(-2.5494247) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1041996) q[1];
sx q[1];
rz(-1.3379119) q[1];
sx q[1];
rz(-1.3523471) q[1];
rz(3.1214633) q[3];
sx q[3];
rz(-1.745589) q[3];
sx q[3];
rz(2.54769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3759489) q[2];
sx q[2];
rz(-2.9889034) q[2];
sx q[2];
rz(-2.622733) q[2];
rz(1.2505924) q[3];
sx q[3];
rz(-1.0873245) q[3];
sx q[3];
rz(2.5605104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64312235) q[0];
sx q[0];
rz(-2.9415218) q[0];
sx q[0];
rz(-0.44922391) q[0];
rz(-1.6953702) q[1];
sx q[1];
rz(-0.64165533) q[1];
sx q[1];
rz(0.62072388) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8410438) q[0];
sx q[0];
rz(-2.4660206) q[0];
sx q[0];
rz(2.9094957) q[0];
x q[1];
rz(-0.16731842) q[2];
sx q[2];
rz(-0.64019055) q[2];
sx q[2];
rz(-0.13969914) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4808958) q[1];
sx q[1];
rz(-2.6891853) q[1];
sx q[1];
rz(0.47641944) q[1];
x q[2];
rz(-0.053106282) q[3];
sx q[3];
rz(-1.0328919) q[3];
sx q[3];
rz(0.90467473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0351403) q[2];
sx q[2];
rz(-1.8702714) q[2];
sx q[2];
rz(-2.980496) q[2];
rz(-2.7437239) q[3];
sx q[3];
rz(-1.6631923) q[3];
sx q[3];
rz(0.14959344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.38254) q[0];
sx q[0];
rz(-1.7790786) q[0];
sx q[0];
rz(-2.8845442) q[0];
rz(2.2611179) q[1];
sx q[1];
rz(-2.5011261) q[1];
sx q[1];
rz(-1.2535198) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.391025) q[0];
sx q[0];
rz(-1.5526958) q[0];
sx q[0];
rz(-1.5347634) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3316111) q[2];
sx q[2];
rz(-1.4000386) q[2];
sx q[2];
rz(2.8851938) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.67318007) q[1];
sx q[1];
rz(-0.77408964) q[1];
sx q[1];
rz(-1.6353398) q[1];
rz(-pi) q[2];
rz(-1.569031) q[3];
sx q[3];
rz(-1.2236508) q[3];
sx q[3];
rz(2.7070759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4132061) q[2];
sx q[2];
rz(-1.9325958) q[2];
sx q[2];
rz(-1.8750635) q[2];
rz(2.5201216) q[3];
sx q[3];
rz(-0.60756835) q[3];
sx q[3];
rz(2.6004041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7414339) q[0];
sx q[0];
rz(-2.3024237) q[0];
sx q[0];
rz(-3.1165282) q[0];
rz(1.9301682) q[1];
sx q[1];
rz(-1.8823267) q[1];
sx q[1];
rz(-2.8866344) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4716041) q[0];
sx q[0];
rz(-1.5758262) q[0];
sx q[0];
rz(0.049335376) q[0];
rz(-pi) q[1];
rz(2.5549501) q[2];
sx q[2];
rz(-1.1544747) q[2];
sx q[2];
rz(3.0468009) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.43422302) q[1];
sx q[1];
rz(-2.1277782) q[1];
sx q[1];
rz(-2.2051478) q[1];
x q[2];
rz(1.8942157) q[3];
sx q[3];
rz(-2.5445523) q[3];
sx q[3];
rz(-0.14123329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8956464) q[2];
sx q[2];
rz(-2.0993555) q[2];
sx q[2];
rz(0.45905217) q[2];
rz(-1.4970655) q[3];
sx q[3];
rz(-0.11681695) q[3];
sx q[3];
rz(0.12115255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6727305) q[0];
sx q[0];
rz(-0.48968306) q[0];
sx q[0];
rz(-0.086061867) q[0];
rz(-0.91880265) q[1];
sx q[1];
rz(-2.0252392) q[1];
sx q[1];
rz(1.2109717) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0491517) q[0];
sx q[0];
rz(-1.7730129) q[0];
sx q[0];
rz(-2.8114955) q[0];
rz(2.2077256) q[2];
sx q[2];
rz(-2.4709765) q[2];
sx q[2];
rz(0.071337168) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0414435) q[1];
sx q[1];
rz(-0.38310928) q[1];
sx q[1];
rz(-0.99975296) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3567621) q[3];
sx q[3];
rz(-0.49477029) q[3];
sx q[3];
rz(0.030232375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3226037) q[2];
sx q[2];
rz(-0.49590597) q[2];
sx q[2];
rz(-0.71713478) q[2];
rz(1.0273733) q[3];
sx q[3];
rz(-1.7337948) q[3];
sx q[3];
rz(-2.7023442) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1502007) q[0];
sx q[0];
rz(-2.1450277) q[0];
sx q[0];
rz(-3.126934) q[0];
rz(-0.75153366) q[1];
sx q[1];
rz(-1.9207759) q[1];
sx q[1];
rz(1.6709447) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.178407) q[0];
sx q[0];
rz(-1.7589594) q[0];
sx q[0];
rz(1.7117731) q[0];
rz(-2.884955) q[2];
sx q[2];
rz(-1.6483288) q[2];
sx q[2];
rz(0.77285779) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4336509) q[1];
sx q[1];
rz(-1.4580863) q[1];
sx q[1];
rz(1.1658843) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5812381) q[3];
sx q[3];
rz(-2.7964948) q[3];
sx q[3];
rz(2.0562248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8810001) q[2];
sx q[2];
rz(-1.1129817) q[2];
sx q[2];
rz(2.8680958) q[2];
rz(1.1547487) q[3];
sx q[3];
rz(-2.4879849) q[3];
sx q[3];
rz(-0.65034136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.038789373) q[0];
sx q[0];
rz(-2.0162855) q[0];
sx q[0];
rz(-2.0661085) q[0];
rz(3.0134046) q[1];
sx q[1];
rz(-0.85103858) q[1];
sx q[1];
rz(2.7959965) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49508383) q[0];
sx q[0];
rz(-2.2301607) q[0];
sx q[0];
rz(-0.20276265) q[0];
x q[1];
rz(-0.32175154) q[2];
sx q[2];
rz(-1.0257799) q[2];
sx q[2];
rz(2.7779752) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.83452713) q[1];
sx q[1];
rz(-0.5448927) q[1];
sx q[1];
rz(0.12172555) q[1];
rz(-pi) q[2];
rz(-2.54353) q[3];
sx q[3];
rz(-1.7081714) q[3];
sx q[3];
rz(-3.1118903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0547611) q[2];
sx q[2];
rz(-2.0589477) q[2];
sx q[2];
rz(1.6205988) q[2];
rz(-1.3711551) q[3];
sx q[3];
rz(-0.75826472) q[3];
sx q[3];
rz(2.7267406) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3138251) q[0];
sx q[0];
rz(-2.1791552) q[0];
sx q[0];
rz(-2.9058822) q[0];
rz(-0.57304263) q[1];
sx q[1];
rz(-2.133281) q[1];
sx q[1];
rz(1.3689573) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9873857) q[0];
sx q[0];
rz(-1.1416178) q[0];
sx q[0];
rz(0.81855075) q[0];
rz(0.70765453) q[2];
sx q[2];
rz(-0.81205149) q[2];
sx q[2];
rz(-2.6132513) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.072242529) q[1];
sx q[1];
rz(-1.4003039) q[1];
sx q[1];
rz(-2.1428698) q[1];
rz(-pi) q[2];
rz(0.24457358) q[3];
sx q[3];
rz(-0.37062708) q[3];
sx q[3];
rz(1.2115492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3706751) q[2];
sx q[2];
rz(-1.8149899) q[2];
sx q[2];
rz(-3.0885922) q[2];
rz(-2.3897589) q[3];
sx q[3];
rz(-1.0337318) q[3];
sx q[3];
rz(0.92541614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5175405) q[0];
sx q[0];
rz(-1.0160099) q[0];
sx q[0];
rz(2.1852063) q[0];
rz(-1.7908295) q[1];
sx q[1];
rz(-0.39719926) q[1];
sx q[1];
rz(-0.15695708) q[1];
rz(1.7012637) q[2];
sx q[2];
rz(-0.87776504) q[2];
sx q[2];
rz(-0.19765111) q[2];
rz(-2.1696321) q[3];
sx q[3];
rz(-2.2836015) q[3];
sx q[3];
rz(-2.9290269) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
