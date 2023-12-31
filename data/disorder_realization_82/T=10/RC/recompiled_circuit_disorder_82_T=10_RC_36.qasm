OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.80469552) q[0];
sx q[0];
rz(-1.0372294) q[0];
sx q[0];
rz(0.35559911) q[0];
rz(2.8388677) q[1];
sx q[1];
rz(5.2390715) q[1];
sx q[1];
rz(10.704438) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.049392603) q[0];
sx q[0];
rz(-0.28486262) q[0];
sx q[0];
rz(-2.0869135) q[0];
rz(-pi) q[1];
rz(2.7202397) q[2];
sx q[2];
rz(-1.3660649) q[2];
sx q[2];
rz(-0.28819627) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.31561786) q[1];
sx q[1];
rz(-2.3992535) q[1];
sx q[1];
rz(-0.62645341) q[1];
x q[2];
rz(-2.9079307) q[3];
sx q[3];
rz(-1.7710925) q[3];
sx q[3];
rz(1.83028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1352284) q[2];
sx q[2];
rz(-0.93868119) q[2];
sx q[2];
rz(-2.485086) q[2];
rz(0.73900977) q[3];
sx q[3];
rz(-2.6813172) q[3];
sx q[3];
rz(2.7242993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1332557) q[0];
sx q[0];
rz(-2.3454741) q[0];
sx q[0];
rz(0.59536368) q[0];
rz(0.061925109) q[1];
sx q[1];
rz(-1.9134816) q[1];
sx q[1];
rz(0.48746902) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66747626) q[0];
sx q[0];
rz(-1.5063018) q[0];
sx q[0];
rz(-1.4897896) q[0];
rz(-pi) q[1];
rz(-1.1335282) q[2];
sx q[2];
rz(-1.184706) q[2];
sx q[2];
rz(-1.3646477) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4844955) q[1];
sx q[1];
rz(-0.24689455) q[1];
sx q[1];
rz(2.2730278) q[1];
x q[2];
rz(1.8295733) q[3];
sx q[3];
rz(-2.8142455) q[3];
sx q[3];
rz(2.5395288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9618824) q[2];
sx q[2];
rz(-1.8790745) q[2];
sx q[2];
rz(0.3953735) q[2];
rz(-2.0987434) q[3];
sx q[3];
rz(-0.53111774) q[3];
sx q[3];
rz(-0.038671967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19668002) q[0];
sx q[0];
rz(-2.0356464) q[0];
sx q[0];
rz(-0.19038598) q[0];
rz(3.0186675) q[1];
sx q[1];
rz(-2.7540837) q[1];
sx q[1];
rz(0.22274676) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5644154) q[0];
sx q[0];
rz(-2.5860997) q[0];
sx q[0];
rz(2.3781858) q[0];
x q[1];
rz(-0.22521714) q[2];
sx q[2];
rz(-1.9334963) q[2];
sx q[2];
rz(-0.97937102) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.45914868) q[1];
sx q[1];
rz(-2.5249081) q[1];
sx q[1];
rz(0.1174121) q[1];
rz(-1.2336897) q[3];
sx q[3];
rz(-1.8030231) q[3];
sx q[3];
rz(0.043134886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8824076) q[2];
sx q[2];
rz(-1.2732482) q[2];
sx q[2];
rz(-1.9474691) q[2];
rz(2.0866701) q[3];
sx q[3];
rz(-1.4849562) q[3];
sx q[3];
rz(0.96364337) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38917437) q[0];
sx q[0];
rz(-2.8635633) q[0];
sx q[0];
rz(1.7383204) q[0];
rz(-1.4933043) q[1];
sx q[1];
rz(-1.5999258) q[1];
sx q[1];
rz(-0.9202252) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7677778) q[0];
sx q[0];
rz(-2.2245363) q[0];
sx q[0];
rz(0.52315229) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4843066) q[2];
sx q[2];
rz(-1.930634) q[2];
sx q[2];
rz(0.53777018) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8060382) q[1];
sx q[1];
rz(-1.1127377) q[1];
sx q[1];
rz(0.6946509) q[1];
rz(-2.782015) q[3];
sx q[3];
rz(-1.7633071) q[3];
sx q[3];
rz(-2.0457091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9437287) q[2];
sx q[2];
rz(-1.4131763) q[2];
sx q[2];
rz(-1.2801923) q[2];
rz(2.3222893) q[3];
sx q[3];
rz(-1.8236482) q[3];
sx q[3];
rz(-1.357648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6827877) q[0];
sx q[0];
rz(-1.7651432) q[0];
sx q[0];
rz(-1.4047594) q[0];
rz(0.76830307) q[1];
sx q[1];
rz(-2.4914425) q[1];
sx q[1];
rz(-2.6884902) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26950726) q[0];
sx q[0];
rz(-3.0591024) q[0];
sx q[0];
rz(-1.1016269) q[0];
x q[1];
rz(-1.808851) q[2];
sx q[2];
rz(-2.9946054) q[2];
sx q[2];
rz(2.5631995) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6140155) q[1];
sx q[1];
rz(-1.7909389) q[1];
sx q[1];
rz(-0.084053587) q[1];
rz(2.1702607) q[3];
sx q[3];
rz(-1.1820275) q[3];
sx q[3];
rz(-2.142981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0129464) q[2];
sx q[2];
rz(-2.3036028) q[2];
sx q[2];
rz(-1.6298693) q[2];
rz(0.72367469) q[3];
sx q[3];
rz(-1.2440163) q[3];
sx q[3];
rz(-2.2912912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.033427514) q[0];
sx q[0];
rz(-1.7280248) q[0];
sx q[0];
rz(-0.23813716) q[0];
rz(2.761633) q[1];
sx q[1];
rz(-1.0521051) q[1];
sx q[1];
rz(-1.628081) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3753525) q[0];
sx q[0];
rz(-2.0413114) q[0];
sx q[0];
rz(1.2067632) q[0];
rz(-2.020535) q[2];
sx q[2];
rz(-1.626484) q[2];
sx q[2];
rz(0.72292098) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5634798) q[1];
sx q[1];
rz(-2.0260749) q[1];
sx q[1];
rz(-0.77225765) q[1];
rz(-pi) q[2];
rz(-1.120818) q[3];
sx q[3];
rz(-1.2202154) q[3];
sx q[3];
rz(-1.9812802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0825519) q[2];
sx q[2];
rz(-2.5566792) q[2];
sx q[2];
rz(-1.9980105) q[2];
rz(-3.011076) q[3];
sx q[3];
rz(-1.7139707) q[3];
sx q[3];
rz(2.9706764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0156353) q[0];
sx q[0];
rz(-0.97706777) q[0];
sx q[0];
rz(2.7440199) q[0];
rz(1.3145087) q[1];
sx q[1];
rz(-0.54388261) q[1];
sx q[1];
rz(1.1869173) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8966973) q[0];
sx q[0];
rz(-0.58422409) q[0];
sx q[0];
rz(-1.0446192) q[0];
rz(-1.5415807) q[2];
sx q[2];
rz(-1.511682) q[2];
sx q[2];
rz(2.7616449) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.742332) q[1];
sx q[1];
rz(-0.84335828) q[1];
sx q[1];
rz(2.9620693) q[1];
x q[2];
rz(-1.2790473) q[3];
sx q[3];
rz(-0.96127931) q[3];
sx q[3];
rz(0.17385829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.26178965) q[2];
sx q[2];
rz(-1.3061085) q[2];
sx q[2];
rz(-1.7857893) q[2];
rz(-1.5073744) q[3];
sx q[3];
rz(-0.58871388) q[3];
sx q[3];
rz(0.9986977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81925201) q[0];
sx q[0];
rz(-1.1871908) q[0];
sx q[0];
rz(2.2858455) q[0];
rz(-0.021727173) q[1];
sx q[1];
rz(-1.0236434) q[1];
sx q[1];
rz(-1.0303248) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.94901) q[0];
sx q[0];
rz(-2.2826676) q[0];
sx q[0];
rz(-2.683995) q[0];
rz(0.048642283) q[2];
sx q[2];
rz(-2.0970793) q[2];
sx q[2];
rz(2.4056048) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.47146637) q[1];
sx q[1];
rz(-1.2738436) q[1];
sx q[1];
rz(2.9832277) q[1];
rz(-pi) q[2];
x q[2];
rz(0.86352591) q[3];
sx q[3];
rz(-2.415495) q[3];
sx q[3];
rz(-1.8973779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8490303) q[2];
sx q[2];
rz(-1.2529255) q[2];
sx q[2];
rz(3.0714152) q[2];
rz(2.3146546) q[3];
sx q[3];
rz(-1.4708054) q[3];
sx q[3];
rz(2.5551445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.4940015) q[0];
sx q[0];
rz(-1.7140472) q[0];
sx q[0];
rz(-0.31627396) q[0];
rz(-1.0519741) q[1];
sx q[1];
rz(-2.4119222) q[1];
sx q[1];
rz(2.0064328) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90056706) q[0];
sx q[0];
rz(-1.7591488) q[0];
sx q[0];
rz(-2.7450949) q[0];
rz(-pi) q[1];
rz(-3.1124434) q[2];
sx q[2];
rz(-1.2505194) q[2];
sx q[2];
rz(-1.0708035) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.89716298) q[1];
sx q[1];
rz(-1.0883372) q[1];
sx q[1];
rz(3.0753115) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9477356) q[3];
sx q[3];
rz(-1.590168) q[3];
sx q[3];
rz(1.3458348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.51745522) q[2];
sx q[2];
rz(-0.74666658) q[2];
sx q[2];
rz(-1.5448145) q[2];
rz(-0.67772135) q[3];
sx q[3];
rz(-0.8422519) q[3];
sx q[3];
rz(-0.28963447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.233376) q[0];
sx q[0];
rz(-0.69013086) q[0];
sx q[0];
rz(0.35183516) q[0];
rz(2.8219163) q[1];
sx q[1];
rz(-2.7668178) q[1];
sx q[1];
rz(0.19616729) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0344714) q[0];
sx q[0];
rz(-2.2959318) q[0];
sx q[0];
rz(-0.79191533) q[0];
x q[1];
rz(3.0912193) q[2];
sx q[2];
rz(-1.7059687) q[2];
sx q[2];
rz(-0.86063517) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0077121) q[1];
sx q[1];
rz(-0.73426437) q[1];
sx q[1];
rz(-1.5478565) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2833061) q[3];
sx q[3];
rz(-2.2755816) q[3];
sx q[3];
rz(-1.5376877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7508042) q[2];
sx q[2];
rz(-1.5020341) q[2];
sx q[2];
rz(-3.0604559) q[2];
rz(1.1817415) q[3];
sx q[3];
rz(-2.2406082) q[3];
sx q[3];
rz(-2.0666163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.022973013) q[0];
sx q[0];
rz(-1.5587627) q[0];
sx q[0];
rz(-1.8297304) q[0];
rz(-1.3148057) q[1];
sx q[1];
rz(-0.79827764) q[1];
sx q[1];
rz(-1.6006443) q[1];
rz(-1.5224456) q[2];
sx q[2];
rz(-1.3370677) q[2];
sx q[2];
rz(1.7481902) q[2];
rz(1.374791) q[3];
sx q[3];
rz(-2.9075665) q[3];
sx q[3];
rz(-2.8961765) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
