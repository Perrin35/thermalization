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
rz(0.6495629) q[0];
sx q[0];
rz(-0.49512884) q[0];
sx q[0];
rz(-2.8861217) q[0];
rz(-0.35222346) q[1];
sx q[1];
rz(-1.7429587) q[1];
sx q[1];
rz(1.7359098) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86780016) q[0];
sx q[0];
rz(-1.5711014) q[0];
sx q[0];
rz(0.0034694151) q[0];
rz(-pi) q[1];
rz(2.954835) q[2];
sx q[2];
rz(-2.6180912) q[2];
sx q[2];
rz(1.7331327) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5454332) q[1];
sx q[1];
rz(-1.1876629) q[1];
sx q[1];
rz(0.0036762357) q[1];
rz(1.6035523) q[3];
sx q[3];
rz(-2.1385953) q[3];
sx q[3];
rz(-0.43772426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.63456717) q[2];
sx q[2];
rz(-0.031012379) q[2];
sx q[2];
rz(-2.3090889) q[2];
rz(-2.3333874) q[3];
sx q[3];
rz(-0.015232239) q[3];
sx q[3];
rz(-0.25850779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38043624) q[0];
sx q[0];
rz(-2.7988837) q[0];
sx q[0];
rz(0.1880745) q[0];
rz(3.0694118) q[1];
sx q[1];
rz(-2.1071823) q[1];
sx q[1];
rz(1.6292876) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8919024) q[0];
sx q[0];
rz(-1.238126) q[0];
sx q[0];
rz(-1.4494677) q[0];
rz(-pi) q[1];
rz(0.030628344) q[2];
sx q[2];
rz(-1.616744) q[2];
sx q[2];
rz(2.6776938) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9413599) q[1];
sx q[1];
rz(-1.5116475) q[1];
sx q[1];
rz(-3.0422387) q[1];
rz(-pi) q[2];
x q[2];
rz(2.333669) q[3];
sx q[3];
rz(-1.1594541) q[3];
sx q[3];
rz(-1.5613256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9709116) q[2];
sx q[2];
rz(-2.4965681) q[2];
sx q[2];
rz(1.8485273) q[2];
rz(-2.7944148) q[3];
sx q[3];
rz(-0.2694338) q[3];
sx q[3];
rz(3.0612883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3229356) q[0];
sx q[0];
rz(-1.289239) q[0];
sx q[0];
rz(0.47054189) q[0];
rz(-1.200354) q[1];
sx q[1];
rz(-2.410694) q[1];
sx q[1];
rz(1.2568731) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6171489) q[0];
sx q[0];
rz(-1.2192508) q[0];
sx q[0];
rz(-1.0008414) q[0];
x q[1];
rz(2.1684709) q[2];
sx q[2];
rz(-2.9785756) q[2];
sx q[2];
rz(-2.4578641) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5611492) q[1];
sx q[1];
rz(-0.014583556) q[1];
sx q[1];
rz(1.4766978) q[1];
rz(-pi) q[2];
rz(1.1670349) q[3];
sx q[3];
rz(-0.65126538) q[3];
sx q[3];
rz(0.854597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5990126) q[2];
sx q[2];
rz(-2.161669) q[2];
sx q[2];
rz(-0.92213321) q[2];
rz(1.3433836) q[3];
sx q[3];
rz(-0.94815367) q[3];
sx q[3];
rz(-1.62742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2682997) q[0];
sx q[0];
rz(-1.0404077) q[0];
sx q[0];
rz(-2.701395) q[0];
rz(1.6104376) q[1];
sx q[1];
rz(-1.4819769) q[1];
sx q[1];
rz(0.24756113) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22641121) q[0];
sx q[0];
rz(-2.0005715) q[0];
sx q[0];
rz(2.129401) q[0];
x q[1];
rz(0.48657067) q[2];
sx q[2];
rz(-0.19127327) q[2];
sx q[2];
rz(-2.4007806) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9369071) q[1];
sx q[1];
rz(-1.2720894) q[1];
sx q[1];
rz(-0.077110962) q[1];
rz(-3.1010755) q[3];
sx q[3];
rz(-2.4926293) q[3];
sx q[3];
rz(1.9891091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2441248) q[2];
sx q[2];
rz(-2.1110057) q[2];
sx q[2];
rz(-1.6991276) q[2];
rz(-2.9407732) q[3];
sx q[3];
rz(-1.2186058) q[3];
sx q[3];
rz(2.0012746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9852801) q[0];
sx q[0];
rz(-2.9887178) q[0];
sx q[0];
rz(2.5961764) q[0];
rz(-0.044005752) q[1];
sx q[1];
rz(-0.017887201) q[1];
sx q[1];
rz(2.6652179) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26075073) q[0];
sx q[0];
rz(-1.6549661) q[0];
sx q[0];
rz(-2.8974786) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.40619855) q[2];
sx q[2];
rz(-1.9425689) q[2];
sx q[2];
rz(2.4570454) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0477476) q[1];
sx q[1];
rz(-2.412891) q[1];
sx q[1];
rz(-1.8456259) q[1];
rz(-pi) q[2];
rz(3.0073193) q[3];
sx q[3];
rz(-2.0760787) q[3];
sx q[3];
rz(-3.0485632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5215317) q[2];
sx q[2];
rz(-2.4510577) q[2];
sx q[2];
rz(-0.36038348) q[2];
rz(-2.9195869) q[3];
sx q[3];
rz(-1.6111776) q[3];
sx q[3];
rz(-1.7300026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5905269) q[0];
sx q[0];
rz(-1.5883625) q[0];
sx q[0];
rz(-1.5850413) q[0];
rz(-3.0146154) q[1];
sx q[1];
rz(-1.304909) q[1];
sx q[1];
rz(3.051905) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68192775) q[0];
sx q[0];
rz(-1.0583504) q[0];
sx q[0];
rz(-2.1905023) q[0];
rz(0.78533919) q[2];
sx q[2];
rz(-2.0302823) q[2];
sx q[2];
rz(1.6096514) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9759102) q[1];
sx q[1];
rz(-1.0780562) q[1];
sx q[1];
rz(1.3085646) q[1];
rz(-pi) q[2];
rz(-0.60097127) q[3];
sx q[3];
rz(-1.8716164) q[3];
sx q[3];
rz(1.5301289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0330607) q[2];
sx q[2];
rz(-2.8046799) q[2];
sx q[2];
rz(2.8899371) q[2];
rz(0.039994914) q[3];
sx q[3];
rz(-0.34188855) q[3];
sx q[3];
rz(2.6778636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43117487) q[0];
sx q[0];
rz(-0.091910563) q[0];
sx q[0];
rz(-0.41496667) q[0];
rz(1.6522495) q[1];
sx q[1];
rz(-0.057561189) q[1];
sx q[1];
rz(0.32589486) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7975066) q[0];
sx q[0];
rz(-2.5523253) q[0];
sx q[0];
rz(-2.2641568) q[0];
rz(-pi) q[1];
rz(-3.0517111) q[2];
sx q[2];
rz(-1.0542068) q[2];
sx q[2];
rz(-0.84026779) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4517026) q[1];
sx q[1];
rz(-2.5696515) q[1];
sx q[1];
rz(-2.6533935) q[1];
rz(0.075906673) q[3];
sx q[3];
rz(-1.2682639) q[3];
sx q[3];
rz(-0.74750604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9140279) q[2];
sx q[2];
rz(-1.808337) q[2];
sx q[2];
rz(-1.0464767) q[2];
rz(0.21876167) q[3];
sx q[3];
rz(-2.1121787) q[3];
sx q[3];
rz(-1.7441162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47531146) q[0];
sx q[0];
rz(-3.1240211) q[0];
sx q[0];
rz(-0.46324357) q[0];
rz(-2.8766368) q[1];
sx q[1];
rz(-0.0015365096) q[1];
sx q[1];
rz(-1.6418246) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2210082) q[0];
sx q[0];
rz(-2.7084368) q[0];
sx q[0];
rz(1.6154352) q[0];
rz(-pi) q[1];
rz(1.6258662) q[2];
sx q[2];
rz(-2.0158421) q[2];
sx q[2];
rz(2.9663756) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3140351) q[1];
sx q[1];
rz(-1.5416939) q[1];
sx q[1];
rz(-1.3489789) q[1];
rz(1.8231976) q[3];
sx q[3];
rz(-1.5034672) q[3];
sx q[3];
rz(-0.52378262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2396592) q[2];
sx q[2];
rz(-2.689211) q[2];
sx q[2];
rz(-1.8233914) q[2];
rz(0.0031331172) q[3];
sx q[3];
rz(-1.2001218) q[3];
sx q[3];
rz(-1.9759294) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5779293) q[0];
sx q[0];
rz(-0.85275537) q[0];
sx q[0];
rz(-0.71459115) q[0];
rz(-2.928012) q[1];
sx q[1];
rz(-0.029284632) q[1];
sx q[1];
rz(-1.2108796) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1909669) q[0];
sx q[0];
rz(-1.4651148) q[0];
sx q[0];
rz(-1.2638541) q[0];
rz(-pi) q[1];
rz(-0.40990746) q[2];
sx q[2];
rz(-2.4033467) q[2];
sx q[2];
rz(-1.8869836) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9309733) q[1];
sx q[1];
rz(-0.60928357) q[1];
sx q[1];
rz(2.0740725) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2641422) q[3];
sx q[3];
rz(-2.8971379) q[3];
sx q[3];
rz(-0.86737421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.96656281) q[2];
sx q[2];
rz(-3.1201456) q[2];
sx q[2];
rz(0.19801298) q[2];
rz(-2.7909732) q[3];
sx q[3];
rz(-1.5703166) q[3];
sx q[3];
rz(-1.143379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7570033) q[0];
sx q[0];
rz(-2.2191255) q[0];
sx q[0];
rz(-2.5272227) q[0];
rz(-0.059582926) q[1];
sx q[1];
rz(-0.091484286) q[1];
sx q[1];
rz(1.4245859) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88286582) q[0];
sx q[0];
rz(-0.75260163) q[0];
sx q[0];
rz(-1.2105788) q[0];
rz(1.5727881) q[2];
sx q[2];
rz(-1.702076) q[2];
sx q[2];
rz(2.7715193) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1150306) q[1];
sx q[1];
rz(-1.4938338) q[1];
sx q[1];
rz(-1.413024) q[1];
x q[2];
rz(-1.102081) q[3];
sx q[3];
rz(-2.6913683) q[3];
sx q[3];
rz(2.1463822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0397348) q[2];
sx q[2];
rz(-0.27145824) q[2];
sx q[2];
rz(-0.52379918) q[2];
rz(1.0456746) q[3];
sx q[3];
rz(-1.8664675) q[3];
sx q[3];
rz(-0.4642134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6687748) q[0];
sx q[0];
rz(-1.5325118) q[0];
sx q[0];
rz(-1.4788628) q[0];
rz(-1.4115903) q[1];
sx q[1];
rz(-0.23163207) q[1];
sx q[1];
rz(-3.0805265) q[1];
rz(-0.642943) q[2];
sx q[2];
rz(-1.3284773) q[2];
sx q[2];
rz(0.20198573) q[2];
rz(0.46661507) q[3];
sx q[3];
rz(-1.7468921) q[3];
sx q[3];
rz(0.56260059) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
