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
rz(-1.2047729) q[0];
rz(1.0869979) q[1];
sx q[1];
rz(-2.6462061) q[1];
sx q[1];
rz(-1.9538716) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76554322) q[0];
sx q[0];
rz(-0.57224579) q[0];
sx q[0];
rz(2.8577198) q[0];
x q[1];
rz(2.0314902) q[2];
sx q[2];
rz(-1.0506127) q[2];
sx q[2];
rz(-0.68797639) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7379308) q[1];
sx q[1];
rz(-1.523535) q[1];
sx q[1];
rz(1.3980327) q[1];
x q[2];
rz(1.9203482) q[3];
sx q[3];
rz(-0.38844019) q[3];
sx q[3];
rz(-0.74439186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5997233) q[2];
sx q[2];
rz(-0.60167998) q[2];
sx q[2];
rz(1.0428693) q[2];
rz(-0.045182191) q[3];
sx q[3];
rz(-2.9729645) q[3];
sx q[3];
rz(1.5359623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0586108) q[0];
sx q[0];
rz(-3.0284212) q[0];
sx q[0];
rz(0.97852069) q[0];
rz(-1.9784031) q[1];
sx q[1];
rz(-0.99449831) q[1];
sx q[1];
rz(1.8189583) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2319671) q[0];
sx q[0];
rz(-0.84950209) q[0];
sx q[0];
rz(2.2366174) q[0];
rz(-1.2873994) q[2];
sx q[2];
rz(-1.7457508) q[2];
sx q[2];
rz(-0.31084479) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1910227) q[1];
sx q[1];
rz(-2.2132067) q[1];
sx q[1];
rz(0.2705785) q[1];
rz(-pi) q[2];
x q[2];
rz(0.69422428) q[3];
sx q[3];
rz(-1.6855918) q[3];
sx q[3];
rz(2.5478884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3375552) q[2];
sx q[2];
rz(-0.20122169) q[2];
sx q[2];
rz(-0.03841722) q[2];
rz(1.5086959) q[3];
sx q[3];
rz(-1.8149866) q[3];
sx q[3];
rz(1.4940777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9134193) q[0];
sx q[0];
rz(-2.4689624) q[0];
sx q[0];
rz(-2.9056554) q[0];
rz(-1.963223) q[1];
sx q[1];
rz(-1.447568) q[1];
sx q[1];
rz(2.6441914) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.936393) q[0];
sx q[0];
rz(-0.25435624) q[0];
sx q[0];
rz(-0.78820552) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1549306) q[2];
sx q[2];
rz(-1.207104) q[2];
sx q[2];
rz(-1.9726582) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6420501) q[1];
sx q[1];
rz(-1.1772369) q[1];
sx q[1];
rz(1.5469013) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8286934) q[3];
sx q[3];
rz(-1.8116496) q[3];
sx q[3];
rz(1.9633213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1824789) q[2];
sx q[2];
rz(-1.6751869) q[2];
sx q[2];
rz(-1.9100995) q[2];
rz(1.4500827) q[3];
sx q[3];
rz(-1.8385889) q[3];
sx q[3];
rz(0.85795295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0813783) q[0];
sx q[0];
rz(-2.2980818) q[0];
sx q[0];
rz(-2.9486935) q[0];
rz(-1.3554205) q[1];
sx q[1];
rz(-1.5813634) q[1];
sx q[1];
rz(0.17098175) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9381704) q[0];
sx q[0];
rz(-1.6294384) q[0];
sx q[0];
rz(-1.0542271) q[0];
rz(-pi) q[1];
rz(-0.73604446) q[2];
sx q[2];
rz(-1.9878329) q[2];
sx q[2];
rz(1.2936178) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1384004) q[1];
sx q[1];
rz(-2.698771) q[1];
sx q[1];
rz(-1.8172311) q[1];
x q[2];
rz(-0.2410197) q[3];
sx q[3];
rz(-0.92151002) q[3];
sx q[3];
rz(0.61862448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.96044797) q[2];
sx q[2];
rz(-1.170271) q[2];
sx q[2];
rz(0.32285264) q[2];
rz(1.3747831) q[3];
sx q[3];
rz(-2.1026473) q[3];
sx q[3];
rz(-0.84097451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5508995) q[0];
sx q[0];
rz(-2.1192079) q[0];
sx q[0];
rz(0.92556959) q[0];
rz(3.0135221) q[1];
sx q[1];
rz(-1.6431199) q[1];
sx q[1];
rz(-3.0421323) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8444237) q[0];
sx q[0];
rz(-2.4250829) q[0];
sx q[0];
rz(-1.5708718) q[0];
x q[1];
rz(-2.7814513) q[2];
sx q[2];
rz(-1.2449054) q[2];
sx q[2];
rz(-1.2831836) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.22597081) q[1];
sx q[1];
rz(-2.5720937) q[1];
sx q[1];
rz(1.2036588) q[1];
rz(3.0614047) q[3];
sx q[3];
rz(-1.8155451) q[3];
sx q[3];
rz(-2.484124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1818992) q[2];
sx q[2];
rz(-1.5965896) q[2];
sx q[2];
rz(0.37364513) q[2];
rz(-0.89961189) q[3];
sx q[3];
rz(-0.74266946) q[3];
sx q[3];
rz(1.1348178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9099092) q[0];
sx q[0];
rz(-2.9250356) q[0];
sx q[0];
rz(0.59189558) q[0];
rz(2.9369211) q[1];
sx q[1];
rz(-2.1494631) q[1];
sx q[1];
rz(0.051102292) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44296023) q[0];
sx q[0];
rz(-1.763183) q[0];
sx q[0];
rz(1.4158184) q[0];
rz(-pi) q[1];
rz(2.883576) q[2];
sx q[2];
rz(-2.2126901) q[2];
sx q[2];
rz(-2.8947865) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2846191) q[1];
sx q[1];
rz(-1.2220807) q[1];
sx q[1];
rz(-0.32964175) q[1];
rz(2.3971812) q[3];
sx q[3];
rz(-2.4648414) q[3];
sx q[3];
rz(1.4686506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.11703141) q[2];
sx q[2];
rz(-1.1144964) q[2];
sx q[2];
rz(2.047211) q[2];
rz(0.38703212) q[3];
sx q[3];
rz(-2.6670167) q[3];
sx q[3];
rz(-2.835137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36560202) q[0];
sx q[0];
rz(-2.2633573) q[0];
sx q[0];
rz(-0.24755092) q[0];
rz(1.6869102) q[1];
sx q[1];
rz(-1.2431815) q[1];
sx q[1];
rz(0.036458485) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.582583) q[0];
sx q[0];
rz(-1.2866308) q[0];
sx q[0];
rz(-2.3733632) q[0];
x q[1];
rz(2.9159732) q[2];
sx q[2];
rz(-2.4258412) q[2];
sx q[2];
rz(2.4168487) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4562021) q[1];
sx q[1];
rz(-2.074132) q[1];
sx q[1];
rz(-0.53629843) q[1];
x q[2];
rz(2.2685675) q[3];
sx q[3];
rz(-3.0254078) q[3];
sx q[3];
rz(-2.1538863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0869007) q[2];
sx q[2];
rz(-0.96618235) q[2];
sx q[2];
rz(-1.7894233) q[2];
rz(-1.8933206) q[3];
sx q[3];
rz(-1.5636445) q[3];
sx q[3];
rz(1.9498391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
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
rz(2.5915204) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4208385) q[0];
sx q[0];
rz(-1.6021358) q[0];
sx q[0];
rz(-0.019958812) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8185316) q[2];
sx q[2];
rz(-1.5093409) q[2];
sx q[2];
rz(2.2814192) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8892563) q[1];
sx q[1];
rz(-2.1057711) q[1];
sx q[1];
rz(-1.757213) q[1];
rz(-pi) q[2];
rz(0.86217238) q[3];
sx q[3];
rz(-0.58364999) q[3];
sx q[3];
rz(0.72008163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.34755808) q[2];
sx q[2];
rz(-1.7358235) q[2];
sx q[2];
rz(-0.75322914) q[2];
rz(0.29427648) q[3];
sx q[3];
rz(-3.0615276) q[3];
sx q[3];
rz(3.0729955) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2733961) q[0];
sx q[0];
rz(-0.28565872) q[0];
sx q[0];
rz(1.4519325) q[0];
rz(2.946335) q[1];
sx q[1];
rz(-1.0804907) q[1];
sx q[1];
rz(1.6498227) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8378252) q[0];
sx q[0];
rz(-1.2653102) q[0];
sx q[0];
rz(-0.49834337) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3765747) q[2];
sx q[2];
rz(-1.6510626) q[2];
sx q[2];
rz(1.2470055) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.10354708) q[1];
sx q[1];
rz(-2.9730151) q[1];
sx q[1];
rz(-1.9650616) q[1];
rz(-pi) q[2];
rz(-1.5641065) q[3];
sx q[3];
rz(-1.5982398) q[3];
sx q[3];
rz(-0.60718482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1008272) q[2];
sx q[2];
rz(-0.46697524) q[2];
sx q[2];
rz(2.2838498) q[2];
rz(-1.6929251) q[3];
sx q[3];
rz(-1.492615) q[3];
sx q[3];
rz(2.9878476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1345074) q[0];
sx q[0];
rz(-0.15468287) q[0];
sx q[0];
rz(0.78306985) q[0];
rz(-2.018441) q[1];
sx q[1];
rz(-1.4479366) q[1];
sx q[1];
rz(0.060308594) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4551082) q[0];
sx q[0];
rz(-1.9251072) q[0];
sx q[0];
rz(0.099781009) q[0];
x q[1];
rz(-2.1570656) q[2];
sx q[2];
rz(-1.6590377) q[2];
sx q[2];
rz(0.72438223) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7544514) q[1];
sx q[1];
rz(-1.6384513) q[1];
sx q[1];
rz(2.4715502) q[1];
rz(2.8496212) q[3];
sx q[3];
rz(-2.7640016) q[3];
sx q[3];
rz(-0.99536125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.15410885) q[2];
sx q[2];
rz(-2.2769112) q[2];
sx q[2];
rz(-1.8314499) q[2];
rz(1.3335258) q[3];
sx q[3];
rz(-1.798505) q[3];
sx q[3];
rz(-1.8346627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(1.9202833) q[3];
sx q[3];
rz(-1.3987473) q[3];
sx q[3];
rz(-1.7940298) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
