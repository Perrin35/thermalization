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
rz(-1.4122352) q[0];
sx q[0];
rz(-0.5923624) q[0];
sx q[0];
rz(1.9368197) q[0];
rz(-2.0545948) q[1];
sx q[1];
rz(-0.49538651) q[1];
sx q[1];
rz(-1.187721) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43150768) q[0];
sx q[0];
rz(-2.1174701) q[0];
sx q[0];
rz(-1.3923079) q[0];
rz(-2.4815791) q[2];
sx q[2];
rz(-0.68049351) q[2];
sx q[2];
rz(1.4726435) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0443327) q[1];
sx q[1];
rz(-9/(16*pi)) q[1];
sx q[1];
rz(-1.3023085) q[1];
x q[2];
rz(-1.9203482) q[3];
sx q[3];
rz(-0.38844019) q[3];
sx q[3];
rz(-2.3972008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5997233) q[2];
sx q[2];
rz(-0.60167998) q[2];
sx q[2];
rz(-2.0987233) q[2];
rz(0.045182191) q[3];
sx q[3];
rz(-0.16862814) q[3];
sx q[3];
rz(-1.6056304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.0829818) q[0];
sx q[0];
rz(-3.0284212) q[0];
sx q[0];
rz(-2.163072) q[0];
rz(-1.1631896) q[1];
sx q[1];
rz(-2.1470943) q[1];
sx q[1];
rz(-1.3226343) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.103687) q[0];
sx q[0];
rz(-0.9390489) q[0];
sx q[0];
rz(-0.61236713) q[0];
x q[1];
rz(1.00707) q[2];
sx q[2];
rz(-2.8097834) q[2];
sx q[2];
rz(2.4203468) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.757431) q[1];
sx q[1];
rz(-2.4520284) q[1];
sx q[1];
rz(1.227725) q[1];
rz(-2.9632872) q[3];
sx q[3];
rz(-2.4394991) q[3];
sx q[3];
rz(0.84035471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3375552) q[2];
sx q[2];
rz(-2.940371) q[2];
sx q[2];
rz(0.03841722) q[2];
rz(-1.6328968) q[3];
sx q[3];
rz(-1.8149866) q[3];
sx q[3];
rz(-1.647515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.9134193) q[0];
sx q[0];
rz(-0.67263022) q[0];
sx q[0];
rz(-2.9056554) q[0];
rz(1.963223) q[1];
sx q[1];
rz(-1.447568) q[1];
sx q[1];
rz(-2.6441914) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20519964) q[0];
sx q[0];
rz(-2.8872364) q[0];
sx q[0];
rz(-0.78820552) q[0];
rz(-pi) q[1];
rz(-2.3264333) q[2];
sx q[2];
rz(-0.54537702) q[2];
sx q[2];
rz(-0.27622414) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.062089109) q[1];
sx q[1];
rz(-1.5928643) q[1];
sx q[1];
rz(0.39366054) q[1];
rz(-0.24876066) q[3];
sx q[3];
rz(-1.3205055) q[3];
sx q[3];
rz(-2.8119025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1824789) q[2];
sx q[2];
rz(-1.6751869) q[2];
sx q[2];
rz(1.2314931) q[2];
rz(1.4500827) q[3];
sx q[3];
rz(-1.8385889) q[3];
sx q[3];
rz(-2.2836397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-2.0813783) q[0];
sx q[0];
rz(-0.84351081) q[0];
sx q[0];
rz(-0.1928992) q[0];
rz(1.3554205) q[1];
sx q[1];
rz(-1.5813634) q[1];
sx q[1];
rz(2.9706109) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6714013) q[0];
sx q[0];
rz(-0.51958749) q[0];
sx q[0];
rz(1.4524824) q[0];
rz(-0.58329432) q[2];
sx q[2];
rz(-2.3153044) q[2];
sx q[2];
rz(0.14329958) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2747171) q[1];
sx q[1];
rz(-1.1422542) q[1];
sx q[1];
rz(0.11517597) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.2410197) q[3];
sx q[3];
rz(-2.2200826) q[3];
sx q[3];
rz(-0.61862448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1811447) q[2];
sx q[2];
rz(-1.170271) q[2];
sx q[2];
rz(-0.32285264) q[2];
rz(-1.3747831) q[3];
sx q[3];
rz(-1.0389453) q[3];
sx q[3];
rz(2.3006181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(1.5508995) q[0];
sx q[0];
rz(-2.1192079) q[0];
sx q[0];
rz(-0.92556959) q[0];
rz(-0.12807056) q[1];
sx q[1];
rz(-1.6431199) q[1];
sx q[1];
rz(0.099460348) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8443237) q[0];
sx q[0];
rz(-0.8542866) q[0];
sx q[0];
rz(6.5690551e-05) q[0];
x q[1];
rz(0.76444605) q[2];
sx q[2];
rz(-2.6607155) q[2];
sx q[2];
rz(2.1493727) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4872553) q[1];
sx q[1];
rz(-1.0434101) q[1];
sx q[1];
rz(-2.9156963) q[1];
rz(-pi) q[2];
rz(0.080187967) q[3];
sx q[3];
rz(-1.3260475) q[3];
sx q[3];
rz(0.65746869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9596935) q[2];
sx q[2];
rz(-1.5965896) q[2];
sx q[2];
rz(2.7679475) q[2];
rz(-0.89961189) q[3];
sx q[3];
rz(-0.74266946) q[3];
sx q[3];
rz(1.1348178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.9099092) q[0];
sx q[0];
rz(-2.9250356) q[0];
sx q[0];
rz(2.5496971) q[0];
rz(-0.20467155) q[1];
sx q[1];
rz(-0.99212956) q[1];
sx q[1];
rz(-0.051102292) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6986324) q[0];
sx q[0];
rz(-1.3784096) q[0];
sx q[0];
rz(1.7257742) q[0];
rz(1.2418397) q[2];
sx q[2];
rz(-2.4566413) q[2];
sx q[2];
rz(-2.9732413) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7390427) q[1];
sx q[1];
rz(-1.2616757) q[1];
sx q[1];
rz(1.9376631) q[1];
rz(-pi) q[2];
rz(-2.0692286) q[3];
sx q[3];
rz(-1.0921156) q[3];
sx q[3];
rz(-2.5414027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0245612) q[2];
sx q[2];
rz(-1.1144964) q[2];
sx q[2];
rz(1.0943817) q[2];
rz(-2.7545605) q[3];
sx q[3];
rz(-2.6670167) q[3];
sx q[3];
rz(0.3064557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36560202) q[0];
sx q[0];
rz(-0.87823534) q[0];
sx q[0];
rz(-2.8940417) q[0];
rz(-1.4546825) q[1];
sx q[1];
rz(-1.8984112) q[1];
sx q[1];
rz(-0.036458485) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74724417) q[0];
sx q[0];
rz(-2.3009662) q[0];
sx q[0];
rz(-1.1850112) q[0];
rz(-pi) q[1];
rz(1.7629303) q[2];
sx q[2];
rz(-0.87687826) q[2];
sx q[2];
rz(1.0200227) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3452975) q[1];
sx q[1];
rz(-2.4234001) q[1];
sx q[1];
rz(-2.3188306) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8730252) q[3];
sx q[3];
rz(-3.0254078) q[3];
sx q[3];
rz(2.1538863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0869007) q[2];
sx q[2];
rz(-2.1754103) q[2];
sx q[2];
rz(-1.3521693) q[2];
rz(-1.2482721) q[3];
sx q[3];
rz(-1.5636445) q[3];
sx q[3];
rz(-1.9498391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2317113) q[0];
sx q[0];
rz(-2.8928962) q[0];
sx q[0];
rz(-1.4924208) q[0];
rz(-0.17852783) q[1];
sx q[1];
rz(-1.2622204) q[1];
sx q[1];
rz(-0.55007225) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.290925) q[0];
sx q[0];
rz(-1.5508473) q[0];
sx q[0];
rz(-1.602142) q[0];
rz(-1.8185316) q[2];
sx q[2];
rz(-1.6322517) q[2];
sx q[2];
rz(0.8601735) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.25233634) q[1];
sx q[1];
rz(-1.0358216) q[1];
sx q[1];
rz(1.757213) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7356803) q[3];
sx q[3];
rz(-2.0024869) q[3];
sx q[3];
rz(1.5188876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.34755808) q[2];
sx q[2];
rz(-1.4057691) q[2];
sx q[2];
rz(0.75322914) q[2];
rz(-2.8473162) q[3];
sx q[3];
rz(-3.0615276) q[3];
sx q[3];
rz(-0.068597138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2733961) q[0];
sx q[0];
rz(-0.28565872) q[0];
sx q[0];
rz(1.4519325) q[0];
rz(0.1952576) q[1];
sx q[1];
rz(-1.0804907) q[1];
sx q[1];
rz(1.4917699) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8378252) q[0];
sx q[0];
rz(-1.8762824) q[0];
sx q[0];
rz(-0.49834337) q[0];
rz(-pi) q[1];
x q[1];
rz(0.765018) q[2];
sx q[2];
rz(-1.6510626) q[2];
sx q[2];
rz(-1.8945872) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6386913) q[1];
sx q[1];
rz(-1.4152621) q[1];
sx q[1];
rz(-3.0763094) q[1];
rz(0.027444124) q[3];
sx q[3];
rz(-1.5774836) q[3];
sx q[3];
rz(2.1777976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0407654) q[2];
sx q[2];
rz(-0.46697524) q[2];
sx q[2];
rz(0.85774285) q[2];
rz(1.6929251) q[3];
sx q[3];
rz(-1.6489776) q[3];
sx q[3];
rz(-0.15374507) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1345074) q[0];
sx q[0];
rz(-0.15468287) q[0];
sx q[0];
rz(2.3585228) q[0];
rz(-2.018441) q[1];
sx q[1];
rz(-1.6936561) q[1];
sx q[1];
rz(-0.060308594) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6864844) q[0];
sx q[0];
rz(-1.2164854) q[0];
sx q[0];
rz(3.0418116) q[0];
rz(3.0357828) q[2];
sx q[2];
rz(-2.1544837) q[2];
sx q[2];
rz(-2.2367144) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.26877182) q[1];
sx q[1];
rz(-2.4686681) q[1];
sx q[1];
rz(-3.0329143) q[1];
rz(-pi) q[2];
x q[2];
rz(0.36300404) q[3];
sx q[3];
rz(-1.6771183) q[3];
sx q[3];
rz(2.8385988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9874838) q[2];
sx q[2];
rz(-0.86468148) q[2];
sx q[2];
rz(1.3101428) q[2];
rz(1.8080669) q[3];
sx q[3];
rz(-1.798505) q[3];
sx q[3];
rz(1.8346627) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6780728) q[0];
sx q[0];
rz(-2.0250043) q[0];
sx q[0];
rz(2.9153839) q[0];
rz(-2.3851867) q[1];
sx q[1];
rz(-2.6466128) q[1];
sx q[1];
rz(-0.80642798) q[1];
rz(0.020515223) q[2];
sx q[2];
rz(-2.7111369) q[2];
sx q[2];
rz(-2.6028056) q[2];
rz(1.2213094) q[3];
sx q[3];
rz(-1.7428453) q[3];
sx q[3];
rz(1.3475628) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
