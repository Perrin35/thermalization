OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2965887) q[0];
sx q[0];
rz(-2.4175329) q[0];
sx q[0];
rz(1.6568503) q[0];
rz(-1.9384664) q[1];
sx q[1];
rz(-2.6180747) q[1];
sx q[1];
rz(0.88820052) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9125497) q[0];
sx q[0];
rz(-2.9648844) q[0];
sx q[0];
rz(-2.8852709) q[0];
rz(-pi) q[1];
rz(-0.16562478) q[2];
sx q[2];
rz(-2.1016444) q[2];
sx q[2];
rz(0.80425516) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8656617) q[1];
sx q[1];
rz(-1.4038424) q[1];
sx q[1];
rz(-1.5702412) q[1];
rz(-pi) q[2];
rz(-0.039460823) q[3];
sx q[3];
rz(-0.69898116) q[3];
sx q[3];
rz(0.82074245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.858294) q[2];
sx q[2];
rz(-0.41559872) q[2];
sx q[2];
rz(-1.1179914) q[2];
rz(0.14532146) q[3];
sx q[3];
rz(-1.558692) q[3];
sx q[3];
rz(0.045923559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5685101) q[0];
sx q[0];
rz(-1.8772323) q[0];
sx q[0];
rz(0.65482393) q[0];
rz(1.2163935) q[1];
sx q[1];
rz(-1.9768068) q[1];
sx q[1];
rz(0.12589802) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.038371041) q[0];
sx q[0];
rz(-1.9403337) q[0];
sx q[0];
rz(-2.080337) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5198031) q[2];
sx q[2];
rz(-0.99025531) q[2];
sx q[2];
rz(-0.34678005) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.241908) q[1];
sx q[1];
rz(-3.0113314) q[1];
sx q[1];
rz(-1.8036519) q[1];
rz(-pi) q[2];
rz(-0.010766518) q[3];
sx q[3];
rz(-2.2942703) q[3];
sx q[3];
rz(2.7303498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0931603) q[2];
sx q[2];
rz(-1.9530714) q[2];
sx q[2];
rz(0.38802567) q[2];
rz(-1.4240501) q[3];
sx q[3];
rz(-2.5035796) q[3];
sx q[3];
rz(0.58732906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5557142) q[0];
sx q[0];
rz(-2.3590187) q[0];
sx q[0];
rz(-0.079332381) q[0];
rz(3.0575867) q[1];
sx q[1];
rz(-2.3386798) q[1];
sx q[1];
rz(1.1598587) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0214329) q[0];
sx q[0];
rz(-2.624357) q[0];
sx q[0];
rz(-1.73818) q[0];
rz(-pi) q[1];
rz(-0.91684219) q[2];
sx q[2];
rz(-0.7407623) q[2];
sx q[2];
rz(-3.0286718) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.68418903) q[1];
sx q[1];
rz(-1.1154798) q[1];
sx q[1];
rz(-1.9940358) q[1];
x q[2];
rz(1.0252762) q[3];
sx q[3];
rz(-1.6189515) q[3];
sx q[3];
rz(-0.84850509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7362061) q[2];
sx q[2];
rz(-1.3233041) q[2];
sx q[2];
rz(-0.1082871) q[2];
rz(0.64374271) q[3];
sx q[3];
rz(-1.0780004) q[3];
sx q[3];
rz(2.5260177) q[3];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9765587) q[0];
sx q[0];
rz(-1.7465916) q[0];
sx q[0];
rz(-1.6595586) q[0];
rz(2.4404793) q[1];
sx q[1];
rz(-1.8891524) q[1];
sx q[1];
rz(2.8569417) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50243044) q[0];
sx q[0];
rz(-0.48361719) q[0];
sx q[0];
rz(1.4101656) q[0];
rz(0.20547159) q[2];
sx q[2];
rz(-1.2887508) q[2];
sx q[2];
rz(0.13946433) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.804867) q[1];
sx q[1];
rz(-1.6972099) q[1];
sx q[1];
rz(3.0625507) q[1];
x q[2];
rz(0.90162006) q[3];
sx q[3];
rz(-1.9603143) q[3];
sx q[3];
rz(-2.8297686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9099137) q[2];
sx q[2];
rz(-0.96857962) q[2];
sx q[2];
rz(-2.5740734) q[2];
rz(-0.41401687) q[3];
sx q[3];
rz(-2.0040138) q[3];
sx q[3];
rz(-2.0295985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48859566) q[0];
sx q[0];
rz(-2.1814006) q[0];
sx q[0];
rz(1.3265142) q[0];
rz(-1.1524221) q[1];
sx q[1];
rz(-1.7633341) q[1];
sx q[1];
rz(-0.93793905) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65813488) q[0];
sx q[0];
rz(-1.5926653) q[0];
sx q[0];
rz(1.4786426) q[0];
rz(-2.6229834) q[2];
sx q[2];
rz(-2.0687639) q[2];
sx q[2];
rz(1.2155611) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2069619) q[1];
sx q[1];
rz(-3.0983371) q[1];
sx q[1];
rz(-2.2224777) q[1];
x q[2];
rz(2.8030843) q[3];
sx q[3];
rz(-2.7673628) q[3];
sx q[3];
rz(0.81774536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9662629) q[2];
sx q[2];
rz(-0.18925174) q[2];
sx q[2];
rz(-1.0413292) q[2];
rz(-1.3025618) q[3];
sx q[3];
rz(-2.0077191) q[3];
sx q[3];
rz(-1.3180102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.702521) q[0];
sx q[0];
rz(-1.1477926) q[0];
sx q[0];
rz(-0.55737108) q[0];
rz(0.56466651) q[1];
sx q[1];
rz(-0.70960418) q[1];
sx q[1];
rz(2.5851137) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5508931) q[0];
sx q[0];
rz(-1.7771582) q[0];
sx q[0];
rz(-2.7569689) q[0];
rz(-pi) q[1];
rz(-0.85530497) q[2];
sx q[2];
rz(-1.4391293) q[2];
sx q[2];
rz(-1.5065187) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6528875) q[1];
sx q[1];
rz(-0.31213752) q[1];
sx q[1];
rz(2.1991792) q[1];
rz(-pi) q[2];
rz(-2.3658386) q[3];
sx q[3];
rz(-1.6072825) q[3];
sx q[3];
rz(-0.072582399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.53211987) q[2];
sx q[2];
rz(-0.76081053) q[2];
sx q[2];
rz(-1.2247941) q[2];
rz(-1.6312284) q[3];
sx q[3];
rz(-1.8211726) q[3];
sx q[3];
rz(-1.640655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0625793) q[0];
sx q[0];
rz(-1.8824848) q[0];
sx q[0];
rz(0.042908948) q[0];
rz(-0.91730109) q[1];
sx q[1];
rz(-0.62364548) q[1];
sx q[1];
rz(2.6409805) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0450889) q[0];
sx q[0];
rz(-1.6052264) q[0];
sx q[0];
rz(1.7486497) q[0];
x q[1];
rz(2.0666276) q[2];
sx q[2];
rz(-0.60056409) q[2];
sx q[2];
rz(1.1662837) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.73270479) q[1];
sx q[1];
rz(-1.9703456) q[1];
sx q[1];
rz(-1.4399745) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9686437) q[3];
sx q[3];
rz(-0.85507353) q[3];
sx q[3];
rz(-2.0778823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6033972) q[2];
sx q[2];
rz(-1.0791225) q[2];
sx q[2];
rz(-0.67561692) q[2];
rz(2.7006941) q[3];
sx q[3];
rz(-1.7023804) q[3];
sx q[3];
rz(0.52136695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.065141) q[0];
sx q[0];
rz(-0.32442176) q[0];
sx q[0];
rz(2.0741529) q[0];
rz(-0.43287977) q[1];
sx q[1];
rz(-1.6378816) q[1];
sx q[1];
rz(1.6759466) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4439125) q[0];
sx q[0];
rz(-1.7051823) q[0];
sx q[0];
rz(2.0356376) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9550214) q[2];
sx q[2];
rz(-1.4636883) q[2];
sx q[2];
rz(-2.4466116) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9626179) q[1];
sx q[1];
rz(-1.0777506) q[1];
sx q[1];
rz(2.7589873) q[1];
rz(-0.64126863) q[3];
sx q[3];
rz(-1.0675758) q[3];
sx q[3];
rz(-0.36676952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.33891588) q[2];
sx q[2];
rz(-2.2854476) q[2];
sx q[2];
rz(1.6652997) q[2];
rz(0.87351292) q[3];
sx q[3];
rz(-0.87681186) q[3];
sx q[3];
rz(2.6749558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2575689) q[0];
sx q[0];
rz(-0.57618657) q[0];
sx q[0];
rz(2.4043758) q[0];
rz(0.018741477) q[1];
sx q[1];
rz(-2.811921) q[1];
sx q[1];
rz(2.2163056) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5697445) q[0];
sx q[0];
rz(-1.0751372) q[0];
sx q[0];
rz(-2.3832688) q[0];
rz(0.66321744) q[2];
sx q[2];
rz(-1.6753917) q[2];
sx q[2];
rz(1.516972) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9097594) q[1];
sx q[1];
rz(-1.4405466) q[1];
sx q[1];
rz(-2.9568683) q[1];
x q[2];
rz(-2.7221189) q[3];
sx q[3];
rz(-0.89722108) q[3];
sx q[3];
rz(0.89299612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.76901889) q[2];
sx q[2];
rz(-1.5344658) q[2];
sx q[2];
rz(2.1378689) q[2];
rz(0.090099661) q[3];
sx q[3];
rz(-0.026244791) q[3];
sx q[3];
rz(1.1130921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1019679) q[0];
sx q[0];
rz(-1.5682546) q[0];
sx q[0];
rz(-1.8819303) q[0];
rz(-2.8758077) q[1];
sx q[1];
rz(-0.62756413) q[1];
sx q[1];
rz(0.75751799) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.950338) q[0];
sx q[0];
rz(-1.3779103) q[0];
sx q[0];
rz(0.14001503) q[0];
x q[1];
rz(-1.3390113) q[2];
sx q[2];
rz(-2.5387562) q[2];
sx q[2];
rz(-1.1834809) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4841008) q[1];
sx q[1];
rz(-1.3490281) q[1];
sx q[1];
rz(-2.0055254) q[1];
rz(3.1086139) q[3];
sx q[3];
rz(-0.71027256) q[3];
sx q[3];
rz(0.23746333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1511128) q[2];
sx q[2];
rz(-1.8258784) q[2];
sx q[2];
rz(0.79375664) q[2];
rz(2.5027067) q[3];
sx q[3];
rz(-2.7975438) q[3];
sx q[3];
rz(-1.0206153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4929852) q[0];
sx q[0];
rz(-2.6633371) q[0];
sx q[0];
rz(2.2289842) q[0];
rz(-1.5007301) q[1];
sx q[1];
rz(-2.224557) q[1];
sx q[1];
rz(1.7932737) q[1];
rz(-1.9008255) q[2];
sx q[2];
rz(-1.1541661) q[2];
sx q[2];
rz(3.0164568) q[2];
rz(1.100148) q[3];
sx q[3];
rz(-0.63334076) q[3];
sx q[3];
rz(0.28950194) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
