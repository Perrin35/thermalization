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
rz(1.7293575) q[0];
sx q[0];
rz(-2.5492302) q[0];
sx q[0];
rz(1.2047729) q[0];
rz(-2.0545948) q[1];
sx q[1];
rz(-0.49538651) q[1];
sx q[1];
rz(-1.187721) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43150768) q[0];
sx q[0];
rz(-2.1174701) q[0];
sx q[0];
rz(-1.7492848) q[0];
rz(0.56894033) q[2];
sx q[2];
rz(-1.1747588) q[2];
sx q[2];
rz(-0.64096156) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.097259911) q[1];
sx q[1];
rz(-9/(16*pi)) q[1];
sx q[1];
rz(-1.3023085) q[1];
rz(-1.2037362) q[3];
sx q[3];
rz(-1.4407183) q[3];
sx q[3];
rz(0.50103044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.54186934) q[2];
sx q[2];
rz(-0.60167998) q[2];
sx q[2];
rz(-2.0987233) q[2];
rz(-0.045182191) q[3];
sx q[3];
rz(-2.9729645) q[3];
sx q[3];
rz(-1.6056304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0586108) q[0];
sx q[0];
rz(-3.0284212) q[0];
sx q[0];
rz(0.97852069) q[0];
rz(1.1631896) q[1];
sx q[1];
rz(-0.99449831) q[1];
sx q[1];
rz(1.8189583) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2319671) q[0];
sx q[0];
rz(-0.84950209) q[0];
sx q[0];
rz(-0.90497525) q[0];
x q[1];
rz(1.00707) q[2];
sx q[2];
rz(-0.33180922) q[2];
sx q[2];
rz(0.72124583) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.757431) q[1];
sx q[1];
rz(-2.4520284) q[1];
sx q[1];
rz(1.9138676) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.69422428) q[3];
sx q[3];
rz(-1.4560008) q[3];
sx q[3];
rz(-0.59370422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.80403745) q[2];
sx q[2];
rz(-2.940371) q[2];
sx q[2];
rz(0.03841722) q[2];
rz(1.5086959) q[3];
sx q[3];
rz(-1.8149866) q[3];
sx q[3];
rz(1.4940777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9134193) q[0];
sx q[0];
rz(-0.67263022) q[0];
sx q[0];
rz(-0.2359373) q[0];
rz(1.1783696) q[1];
sx q[1];
rz(-1.447568) q[1];
sx q[1];
rz(-0.49740121) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.936393) q[0];
sx q[0];
rz(-0.25435624) q[0];
sx q[0];
rz(-2.3533871) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3264333) q[2];
sx q[2];
rz(-2.5962156) q[2];
sx q[2];
rz(0.27622414) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5617852) q[1];
sx q[1];
rz(-0.39424636) q[1];
sx q[1];
rz(-3.0841139) q[1];
rz(-pi) q[2];
rz(1.8286934) q[3];
sx q[3];
rz(-1.329943) q[3];
sx q[3];
rz(-1.1782714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.95911372) q[2];
sx q[2];
rz(-1.4664058) q[2];
sx q[2];
rz(1.2314931) q[2];
rz(1.4500827) q[3];
sx q[3];
rz(-1.3030038) q[3];
sx q[3];
rz(2.2836397) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0602144) q[0];
sx q[0];
rz(-2.2980818) q[0];
sx q[0];
rz(-2.9486935) q[0];
rz(-1.7861722) q[1];
sx q[1];
rz(-1.5602292) q[1];
sx q[1];
rz(0.17098175) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9381704) q[0];
sx q[0];
rz(-1.6294384) q[0];
sx q[0];
rz(-1.0542271) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58329432) q[2];
sx q[2];
rz(-0.82628822) q[2];
sx q[2];
rz(2.9982931) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2747171) q[1];
sx q[1];
rz(-1.1422542) q[1];
sx q[1];
rz(3.0264167) q[1];
x q[2];
rz(1.2661352) q[3];
sx q[3];
rz(-0.6864635) q[3];
sx q[3];
rz(-1.0047508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1811447) q[2];
sx q[2];
rz(-1.170271) q[2];
sx q[2];
rz(0.32285264) q[2];
rz(1.7668096) q[3];
sx q[3];
rz(-1.0389453) q[3];
sx q[3];
rz(-0.84097451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5508995) q[0];
sx q[0];
rz(-1.0223848) q[0];
sx q[0];
rz(0.92556959) q[0];
rz(0.12807056) q[1];
sx q[1];
rz(-1.6431199) q[1];
sx q[1];
rz(3.0421323) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29716897) q[0];
sx q[0];
rz(-0.71650973) q[0];
sx q[0];
rz(1.5707209) q[0];
rz(-pi) q[1];
rz(2.3771466) q[2];
sx q[2];
rz(-0.4808772) q[2];
sx q[2];
rz(-0.99221992) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.22597081) q[1];
sx q[1];
rz(-2.5720937) q[1];
sx q[1];
rz(-1.2036588) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8811536) q[3];
sx q[3];
rz(-2.884293) q[3];
sx q[3];
rz(-2.8043487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9596935) q[2];
sx q[2];
rz(-1.5965896) q[2];
sx q[2];
rz(-2.7679475) q[2];
rz(2.2419808) q[3];
sx q[3];
rz(-2.3989232) q[3];
sx q[3];
rz(2.0067748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9099092) q[0];
sx q[0];
rz(-0.21655701) q[0];
sx q[0];
rz(2.5496971) q[0];
rz(-2.9369211) q[1];
sx q[1];
rz(-2.1494631) q[1];
sx q[1];
rz(-0.051102292) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0979733) q[0];
sx q[0];
rz(-1.7228925) q[0];
sx q[0];
rz(0.19466227) q[0];
x q[1];
rz(1.2418397) q[2];
sx q[2];
rz(-2.4566413) q[2];
sx q[2];
rz(-2.9732413) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.49858677) q[1];
sx q[1];
rz(-2.6664147) q[1];
sx q[1];
rz(2.2982486) q[1];
rz(-pi) q[2];
rz(-2.6079518) q[3];
sx q[3];
rz(-2.009005) q[3];
sx q[3];
rz(-2.4166188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0245612) q[2];
sx q[2];
rz(-2.0270963) q[2];
sx q[2];
rz(1.0943817) q[2];
rz(-0.38703212) q[3];
sx q[3];
rz(-0.47457591) q[3];
sx q[3];
rz(0.3064557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7759906) q[0];
sx q[0];
rz(-2.2633573) q[0];
sx q[0];
rz(2.8940417) q[0];
rz(-1.4546825) q[1];
sx q[1];
rz(-1.8984112) q[1];
sx q[1];
rz(-0.036458485) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.582583) q[0];
sx q[0];
rz(-1.2866308) q[0];
sx q[0];
rz(2.3733632) q[0];
rz(2.4385284) q[2];
sx q[2];
rz(-1.7181239) q[2];
sx q[2];
rz(-2.4670403) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9769638) q[1];
sx q[1];
rz(-1.1067302) q[1];
sx q[1];
rz(1.0010757) q[1];
x q[2];
rz(-0.8730252) q[3];
sx q[3];
rz(-3.0254078) q[3];
sx q[3];
rz(0.98770638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0869007) q[2];
sx q[2];
rz(-0.96618235) q[2];
sx q[2];
rz(1.3521693) q[2];
rz(1.2482721) q[3];
sx q[3];
rz(-1.5636445) q[3];
sx q[3];
rz(-1.1917535) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9098814) q[0];
sx q[0];
rz(-2.8928962) q[0];
sx q[0];
rz(1.6491718) q[0];
rz(-2.9630648) q[1];
sx q[1];
rz(-1.8793722) q[1];
sx q[1];
rz(-0.55007225) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.290925) q[0];
sx q[0];
rz(-1.5907453) q[0];
sx q[0];
rz(1.602142) q[0];
x q[1];
rz(-1.8185316) q[2];
sx q[2];
rz(-1.5093409) q[2];
sx q[2];
rz(-0.8601735) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4143205) q[1];
sx q[1];
rz(-1.4106693) q[1];
sx q[1];
rz(-0.54267197) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2794203) q[3];
sx q[3];
rz(-2.5579427) q[3];
sx q[3];
rz(0.72008163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.34755808) q[2];
sx q[2];
rz(-1.7358235) q[2];
sx q[2];
rz(0.75322914) q[2];
rz(-2.8473162) q[3];
sx q[3];
rz(-0.080065057) q[3];
sx q[3];
rz(0.068597138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8681965) q[0];
sx q[0];
rz(-0.28565872) q[0];
sx q[0];
rz(1.4519325) q[0];
rz(0.1952576) q[1];
sx q[1];
rz(-1.0804907) q[1];
sx q[1];
rz(-1.6498227) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0367835) q[0];
sx q[0];
rz(-1.097479) q[0];
sx q[0];
rz(1.2261054) q[0];
x q[1];
rz(-0.765018) q[2];
sx q[2];
rz(-1.6510626) q[2];
sx q[2];
rz(1.8945872) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0780219) q[1];
sx q[1];
rz(-1.5063022) q[1];
sx q[1];
rz(-1.4149354) q[1];
rz(-3.1141485) q[3];
sx q[3];
rz(-1.5774836) q[3];
sx q[3];
rz(-0.96379507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0407654) q[2];
sx q[2];
rz(-0.46697524) q[2];
sx q[2];
rz(2.2838498) q[2];
rz(-1.4486676) q[3];
sx q[3];
rz(-1.492615) q[3];
sx q[3];
rz(-2.9878476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.0070852) q[0];
sx q[0];
rz(-2.9869098) q[0];
sx q[0];
rz(0.78306985) q[0];
rz(-2.018441) q[1];
sx q[1];
rz(-1.6936561) q[1];
sx q[1];
rz(3.0812841) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9911847) q[0];
sx q[0];
rz(-1.4772319) q[0];
sx q[0];
rz(1.2148598) q[0];
x q[1];
rz(2.1570656) q[2];
sx q[2];
rz(-1.6590377) q[2];
sx q[2];
rz(2.4172104) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0114501) q[1];
sx q[1];
rz(-0.9025652) q[1];
sx q[1];
rz(-1.6570309) q[1];
rz(-pi) q[2];
rz(2.7785886) q[3];
sx q[3];
rz(-1.6771183) q[3];
sx q[3];
rz(0.3029938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.15410885) q[2];
sx q[2];
rz(-2.2769112) q[2];
sx q[2];
rz(1.8314499) q[2];
rz(1.3335258) q[3];
sx q[3];
rz(-1.798505) q[3];
sx q[3];
rz(-1.8346627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46351984) q[0];
sx q[0];
rz(-1.1165883) q[0];
sx q[0];
rz(-0.22620871) q[0];
rz(0.75640596) q[1];
sx q[1];
rz(-2.6466128) q[1];
sx q[1];
rz(-0.80642798) q[1];
rz(-0.020515223) q[2];
sx q[2];
rz(-0.43045577) q[2];
sx q[2];
rz(0.53878709) q[2];
rz(1.1011878) q[3];
sx q[3];
rz(-2.7536177) q[3];
sx q[3];
rz(-2.9256647) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
