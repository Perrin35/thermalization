OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4410285) q[0];
sx q[0];
rz(-1.1428042) q[0];
sx q[0];
rz(-1.2115275) q[0];
rz(-0.22663528) q[1];
sx q[1];
rz(-1.5770788) q[1];
sx q[1];
rz(0.29830631) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.053781833) q[0];
sx q[0];
rz(-0.90667533) q[0];
sx q[0];
rz(-1.5009297) q[0];
rz(-pi) q[1];
rz(3.0478893) q[2];
sx q[2];
rz(-1.9022577) q[2];
sx q[2];
rz(1.9625488) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4264307) q[1];
sx q[1];
rz(-2.3843345) q[1];
sx q[1];
rz(-2.2992579) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.88793036) q[3];
sx q[3];
rz(-0.13266064) q[3];
sx q[3];
rz(2.5761029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4937218) q[2];
sx q[2];
rz(-1.9335258) q[2];
sx q[2];
rz(-2.1477264) q[2];
rz(-2.1422051) q[3];
sx q[3];
rz(-1.2402273) q[3];
sx q[3];
rz(-2.5527111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64269972) q[0];
sx q[0];
rz(-1.9357341) q[0];
sx q[0];
rz(3.1233741) q[0];
rz(2.3253564) q[1];
sx q[1];
rz(-2.1111592) q[1];
sx q[1];
rz(2.6699064) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6149711) q[0];
sx q[0];
rz(-1.7978298) q[0];
sx q[0];
rz(-1.4814266) q[0];
rz(-pi) q[1];
rz(0.50498982) q[2];
sx q[2];
rz(-1.8171176) q[2];
sx q[2];
rz(1.9976975) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.43863152) q[1];
sx q[1];
rz(-1.303181) q[1];
sx q[1];
rz(-0.57499927) q[1];
rz(-pi) q[2];
rz(-0.96305965) q[3];
sx q[3];
rz(-1.4415603) q[3];
sx q[3];
rz(2.3605763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5919684) q[2];
sx q[2];
rz(-1.9886118) q[2];
sx q[2];
rz(0.96898752) q[2];
rz(0.5747059) q[3];
sx q[3];
rz(-2.59022) q[3];
sx q[3];
rz(-1.0415174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4330924) q[0];
sx q[0];
rz(-1.0555462) q[0];
sx q[0];
rz(2.95978) q[0];
rz(1.1026985) q[1];
sx q[1];
rz(-1.5010553) q[1];
sx q[1];
rz(1.6859432) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11765471) q[0];
sx q[0];
rz(-2.0485989) q[0];
sx q[0];
rz(-0.91207232) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.45708926) q[2];
sx q[2];
rz(-2.3466913) q[2];
sx q[2];
rz(-2.7283816) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1118483) q[1];
sx q[1];
rz(-1.4451471) q[1];
sx q[1];
rz(-1.5813584) q[1];
rz(0.62526838) q[3];
sx q[3];
rz(-1.6571952) q[3];
sx q[3];
rz(-2.3141935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2640947) q[2];
sx q[2];
rz(-1.4870746) q[2];
sx q[2];
rz(2.1739615) q[2];
rz(2.4140221) q[3];
sx q[3];
rz(-1.8811767) q[3];
sx q[3];
rz(-2.9038866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2847292) q[0];
sx q[0];
rz(-0.52607042) q[0];
sx q[0];
rz(2.5033584) q[0];
rz(-1.1278641) q[1];
sx q[1];
rz(-0.82740873) q[1];
sx q[1];
rz(1.2329873) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5947943) q[0];
sx q[0];
rz(-1.6398755) q[0];
sx q[0];
rz(-2.6462206) q[0];
rz(2.1190676) q[2];
sx q[2];
rz(-1.0523445) q[2];
sx q[2];
rz(-0.099629121) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9150881) q[1];
sx q[1];
rz(-2.9595032) q[1];
sx q[1];
rz(-1.7712797) q[1];
rz(-1.2992371) q[3];
sx q[3];
rz(-2.4304667) q[3];
sx q[3];
rz(1.7600876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2234852) q[2];
sx q[2];
rz(-1.2780259) q[2];
sx q[2];
rz(-2.3275862) q[2];
rz(1.043184) q[3];
sx q[3];
rz(-0.63101763) q[3];
sx q[3];
rz(-1.1842747) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2994613) q[0];
sx q[0];
rz(-1.3548387) q[0];
sx q[0];
rz(0.89170757) q[0];
rz(-1.2437598) q[1];
sx q[1];
rz(-1.7638821) q[1];
sx q[1];
rz(2.9290501) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71056238) q[0];
sx q[0];
rz(-3.115603) q[0];
sx q[0];
rz(-0.89528577) q[0];
rz(1.6264621) q[2];
sx q[2];
rz(-2.5362483) q[2];
sx q[2];
rz(-2.2948613) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7674539) q[1];
sx q[1];
rz(-0.98857388) q[1];
sx q[1];
rz(0.9064845) q[1];
x q[2];
rz(2.6796474) q[3];
sx q[3];
rz(-2.2432703) q[3];
sx q[3];
rz(2.7355821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.033096878) q[2];
sx q[2];
rz(-2.1704845) q[2];
sx q[2];
rz(2.6203716) q[2];
rz(1.3850348) q[3];
sx q[3];
rz(-1.228046) q[3];
sx q[3];
rz(-1.8732171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2824771) q[0];
sx q[0];
rz(-2.2127667) q[0];
sx q[0];
rz(-2.916472) q[0];
rz(1.3549995) q[1];
sx q[1];
rz(-1.0083818) q[1];
sx q[1];
rz(-2.7640142) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3513983) q[0];
sx q[0];
rz(-1.5488008) q[0];
sx q[0];
rz(-1.8368594) q[0];
rz(-pi) q[1];
rz(2.8541366) q[2];
sx q[2];
rz(-1.524458) q[2];
sx q[2];
rz(3.1099144) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4695417) q[1];
sx q[1];
rz(-2.6968323) q[1];
sx q[1];
rz(-1.3252392) q[1];
x q[2];
rz(-0.029929786) q[3];
sx q[3];
rz(-2.3464591) q[3];
sx q[3];
rz(-0.078775725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.98465115) q[2];
sx q[2];
rz(-0.7876544) q[2];
sx q[2];
rz(0.48103508) q[2];
rz(-0.40361079) q[3];
sx q[3];
rz(-1.0389046) q[3];
sx q[3];
rz(-0.31479442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0890546) q[0];
sx q[0];
rz(-2.5367694) q[0];
sx q[0];
rz(2.9470434) q[0];
rz(-0.21952195) q[1];
sx q[1];
rz(-1.6794645) q[1];
sx q[1];
rz(-2.887168) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97878153) q[0];
sx q[0];
rz(-0.33919507) q[0];
sx q[0];
rz(1.2374872) q[0];
x q[1];
rz(-2.1830325) q[2];
sx q[2];
rz(-1.4441635) q[2];
sx q[2];
rz(-0.74778344) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3212657) q[1];
sx q[1];
rz(-1.1202381) q[1];
sx q[1];
rz(1.2688368) q[1];
rz(-pi) q[2];
rz(0.80612225) q[3];
sx q[3];
rz(-2.5297909) q[3];
sx q[3];
rz(-2.7968189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1640132) q[2];
sx q[2];
rz(-1.7501202) q[2];
sx q[2];
rz(1.3158201) q[2];
rz(-0.87604648) q[3];
sx q[3];
rz(-3.0026569) q[3];
sx q[3];
rz(2.1379437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2106237) q[0];
sx q[0];
rz(-2.7656778) q[0];
sx q[0];
rz(1.4550495) q[0];
rz(-0.82398206) q[1];
sx q[1];
rz(-2.1897557) q[1];
sx q[1];
rz(-1.5751858) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5113735) q[0];
sx q[0];
rz(-2.2492118) q[0];
sx q[0];
rz(1.9766962) q[0];
rz(-pi) q[1];
rz(0.17922108) q[2];
sx q[2];
rz(-2.1041098) q[2];
sx q[2];
rz(-0.78782493) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6187001) q[1];
sx q[1];
rz(-0.52378264) q[1];
sx q[1];
rz(2.3431542) q[1];
x q[2];
rz(2.6353587) q[3];
sx q[3];
rz(-1.1971548) q[3];
sx q[3];
rz(2.1831499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.87551293) q[2];
sx q[2];
rz(-1.9078887) q[2];
sx q[2];
rz(-2.8640462) q[2];
rz(1.4510441) q[3];
sx q[3];
rz(-2.6896559) q[3];
sx q[3];
rz(-1.014876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16185109) q[0];
sx q[0];
rz(-0.45270544) q[0];
sx q[0];
rz(1.6850527) q[0];
rz(-0.62943554) q[1];
sx q[1];
rz(-1.9742191) q[1];
sx q[1];
rz(1.1368407) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0757383) q[0];
sx q[0];
rz(-0.59959164) q[0];
sx q[0];
rz(2.0980741) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2898229) q[2];
sx q[2];
rz(-0.53817828) q[2];
sx q[2];
rz(-2.1449508) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9469103) q[1];
sx q[1];
rz(-1.2543251) q[1];
sx q[1];
rz(3.0419993) q[1];
rz(2.4056899) q[3];
sx q[3];
rz(-1.2888442) q[3];
sx q[3];
rz(-1.499093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4920766) q[2];
sx q[2];
rz(-1.4294383) q[2];
sx q[2];
rz(-2.1006404) q[2];
rz(3.1395636) q[3];
sx q[3];
rz(-0.40015951) q[3];
sx q[3];
rz(0.31203976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3025538) q[0];
sx q[0];
rz(-2.9170687) q[0];
sx q[0];
rz(-0.94605207) q[0];
rz(0.91167766) q[1];
sx q[1];
rz(-1.2152351) q[1];
sx q[1];
rz(2.5295703) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5123972) q[0];
sx q[0];
rz(-2.0266268) q[0];
sx q[0];
rz(-2.3770611) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3304747) q[2];
sx q[2];
rz(-1.2346621) q[2];
sx q[2];
rz(2.0517595) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0079841) q[1];
sx q[1];
rz(-0.66966479) q[1];
sx q[1];
rz(-2.942251) q[1];
rz(-pi) q[2];
rz(-1.9480115) q[3];
sx q[3];
rz(-1.0672788) q[3];
sx q[3];
rz(2.6138888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0570021) q[2];
sx q[2];
rz(-2.4995063) q[2];
sx q[2];
rz(2.4882312) q[2];
rz(-0.35081321) q[3];
sx q[3];
rz(-1.6143129) q[3];
sx q[3];
rz(-2.4408834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6012797) q[0];
sx q[0];
rz(-1.5355587) q[0];
sx q[0];
rz(-3.0328947) q[0];
rz(-0.75469771) q[1];
sx q[1];
rz(-1.3194059) q[1];
sx q[1];
rz(-1.5059765) q[1];
rz(1.0036219) q[2];
sx q[2];
rz(-0.97273172) q[2];
sx q[2];
rz(1.6956971) q[2];
rz(-1.3238974) q[3];
sx q[3];
rz(-2.8084451) q[3];
sx q[3];
rz(0.13767903) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
