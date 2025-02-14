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
rz(2.7163765) q[0];
sx q[0];
rz(-1.334231) q[0];
sx q[0];
rz(0.38036007) q[0];
rz(1.7827787) q[1];
sx q[1];
rz(-0.18667297) q[1];
sx q[1];
rz(0.92460728) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3596852) q[0];
sx q[0];
rz(-1.9115907) q[0];
sx q[0];
rz(-2.4016892) q[0];
rz(-1.2819498) q[2];
sx q[2];
rz(-0.53448662) q[2];
sx q[2];
rz(1.6030902) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.29942611) q[1];
sx q[1];
rz(-1.5314845) q[1];
sx q[1];
rz(-0.8419324) q[1];
rz(-pi) q[2];
rz(0.12792952) q[3];
sx q[3];
rz(-1.5797714) q[3];
sx q[3];
rz(2.8913446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1251462) q[2];
sx q[2];
rz(-0.95637286) q[2];
sx q[2];
rz(2.1535786) q[2];
rz(0.2068578) q[3];
sx q[3];
rz(-1.7346953) q[3];
sx q[3];
rz(-0.44001165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67966953) q[0];
sx q[0];
rz(-1.672687) q[0];
sx q[0];
rz(0.2521421) q[0];
rz(-0.02154669) q[1];
sx q[1];
rz(-0.69308678) q[1];
sx q[1];
rz(-2.2861939) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1148672) q[0];
sx q[0];
rz(-1.572319) q[0];
sx q[0];
rz(-1.5978088) q[0];
rz(2.5874675) q[2];
sx q[2];
rz(-1.9337855) q[2];
sx q[2];
rz(-2.3100694) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.336074) q[1];
sx q[1];
rz(-1.0401498) q[1];
sx q[1];
rz(-0.25605274) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7383201) q[3];
sx q[3];
rz(-2.4174815) q[3];
sx q[3];
rz(1.4693174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1713193) q[2];
sx q[2];
rz(-3.1372941) q[2];
sx q[2];
rz(2.8805736) q[2];
rz(-0.039693443) q[3];
sx q[3];
rz(-1.7429765) q[3];
sx q[3];
rz(-0.54350129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4259341) q[0];
sx q[0];
rz(-1.6682699) q[0];
sx q[0];
rz(0.69450992) q[0];
rz(-2.014324) q[1];
sx q[1];
rz(-2.290461) q[1];
sx q[1];
rz(2.1515501) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9309826) q[0];
sx q[0];
rz(-1.6525804) q[0];
sx q[0];
rz(-2.7922996) q[0];
x q[1];
rz(0.98857356) q[2];
sx q[2];
rz(-1.4043839) q[2];
sx q[2];
rz(-0.63936448) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7496365) q[1];
sx q[1];
rz(-2.5534494) q[1];
sx q[1];
rz(2.4571477) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6022791) q[3];
sx q[3];
rz(-2.1997872) q[3];
sx q[3];
rz(1.9696389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8433044) q[2];
sx q[2];
rz(-0.97323209) q[2];
sx q[2];
rz(-2.6514371) q[2];
rz(0.35689029) q[3];
sx q[3];
rz(-0.39339742) q[3];
sx q[3];
rz(-1.2662158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.056034293) q[0];
sx q[0];
rz(-1.9558676) q[0];
sx q[0];
rz(-2.2991142) q[0];
rz(-0.8017686) q[1];
sx q[1];
rz(-0.27920488) q[1];
sx q[1];
rz(0.78559771) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88670896) q[0];
sx q[0];
rz(-0.87866106) q[0];
sx q[0];
rz(1.2413625) q[0];
rz(-pi) q[1];
rz(0.17274348) q[2];
sx q[2];
rz(-0.72451353) q[2];
sx q[2];
rz(-0.24562626) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3093058) q[1];
sx q[1];
rz(-0.63696276) q[1];
sx q[1];
rz(-2.6061406) q[1];
rz(2.3219452) q[3];
sx q[3];
rz(-1.420212) q[3];
sx q[3];
rz(-0.25170799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.12513146) q[2];
sx q[2];
rz(-0.74840122) q[2];
sx q[2];
rz(2.1330736) q[2];
rz(2.290001) q[3];
sx q[3];
rz(-2.3840756) q[3];
sx q[3];
rz(-1.0803224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(1.9652902) q[0];
sx q[0];
rz(-0.98133123) q[0];
sx q[0];
rz(1.0705795) q[0];
rz(-0.77752441) q[1];
sx q[1];
rz(-2.2309525) q[1];
sx q[1];
rz(2.1308965) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.010653822) q[0];
sx q[0];
rz(-1.9294327) q[0];
sx q[0];
rz(2.6333195) q[0];
x q[1];
rz(-1.2161184) q[2];
sx q[2];
rz(-0.18624072) q[2];
sx q[2];
rz(0.29224381) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.1221558) q[1];
sx q[1];
rz(-2.2194462) q[1];
sx q[1];
rz(2.1557356) q[1];
rz(-pi) q[2];
rz(-1.6701974) q[3];
sx q[3];
rz(-2.1306921) q[3];
sx q[3];
rz(-0.82370629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8396478) q[2];
sx q[2];
rz(-3.0366615) q[2];
sx q[2];
rz(1.5555752) q[2];
rz(-1.0157061) q[3];
sx q[3];
rz(-2.000122) q[3];
sx q[3];
rz(-1.7293845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79065901) q[0];
sx q[0];
rz(-0.18336329) q[0];
sx q[0];
rz(-2.3405128) q[0];
rz(-3.0084897) q[1];
sx q[1];
rz(-1.8201273) q[1];
sx q[1];
rz(-0.4020234) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4426081) q[0];
sx q[0];
rz(-2.3056718) q[0];
sx q[0];
rz(-2.0280272) q[0];
x q[1];
rz(1.5779675) q[2];
sx q[2];
rz(-0.29080393) q[2];
sx q[2];
rz(-2.4343642) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.133144) q[1];
sx q[1];
rz(-2.5767972) q[1];
sx q[1];
rz(-1.0393591) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4834255) q[3];
sx q[3];
rz(-1.6699381) q[3];
sx q[3];
rz(-2.1226573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5037527) q[2];
sx q[2];
rz(-1.736234) q[2];
sx q[2];
rz(2.3657738) q[2];
rz(-1.6860298) q[3];
sx q[3];
rz(-1.9676696) q[3];
sx q[3];
rz(2.1337401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.044947226) q[0];
sx q[0];
rz(-0.89770397) q[0];
sx q[0];
rz(-2.4626379) q[0];
rz(-1.2848162) q[1];
sx q[1];
rz(-2.4285451) q[1];
sx q[1];
rz(1.3791893) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5011936) q[0];
sx q[0];
rz(-1.7441347) q[0];
sx q[0];
rz(2.3480183) q[0];
rz(-pi) q[1];
rz(-0.083456434) q[2];
sx q[2];
rz(-1.2922568) q[2];
sx q[2];
rz(-2.4065774) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0028982) q[1];
sx q[1];
rz(-2.0636673) q[1];
sx q[1];
rz(-1.0652131) q[1];
rz(-1.7709183) q[3];
sx q[3];
rz(-1.516023) q[3];
sx q[3];
rz(1.9211662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3369559) q[2];
sx q[2];
rz(-2.3865484) q[2];
sx q[2];
rz(-1.2389368) q[2];
rz(1.8094481) q[3];
sx q[3];
rz(-0.76056162) q[3];
sx q[3];
rz(-0.24470394) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7234583) q[0];
sx q[0];
rz(-1.9754388) q[0];
sx q[0];
rz(0.13846692) q[0];
rz(0.18374099) q[1];
sx q[1];
rz(-2.467149) q[1];
sx q[1];
rz(0.26652452) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4999102) q[0];
sx q[0];
rz(-1.8962093) q[0];
sx q[0];
rz(2.9703559) q[0];
rz(-pi) q[1];
rz(-1.1600003) q[2];
sx q[2];
rz(-2.5345384) q[2];
sx q[2];
rz(1.0714873) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.16921088) q[1];
sx q[1];
rz(-1.5145497) q[1];
sx q[1];
rz(1.1509658) q[1];
x q[2];
rz(-0.72690225) q[3];
sx q[3];
rz(-1.57734) q[3];
sx q[3];
rz(2.9605799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0900241) q[2];
sx q[2];
rz(-1.8847621) q[2];
sx q[2];
rz(2.9070692) q[2];
rz(1.4759493) q[3];
sx q[3];
rz(-1.539295) q[3];
sx q[3];
rz(-0.75638151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.804857) q[0];
sx q[0];
rz(-0.8661626) q[0];
sx q[0];
rz(0.79972237) q[0];
rz(-2.5526478) q[1];
sx q[1];
rz(-0.84732333) q[1];
sx q[1];
rz(-0.2074997) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8556535) q[0];
sx q[0];
rz(-2.7253599) q[0];
sx q[0];
rz(2.8134384) q[0];
x q[1];
rz(-0.3090119) q[2];
sx q[2];
rz(-2.2591619) q[2];
sx q[2];
rz(-2.7613044) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5415685) q[1];
sx q[1];
rz(-2.3316521) q[1];
sx q[1];
rz(-1.0066973) q[1];
rz(-pi) q[2];
rz(3.112821) q[3];
sx q[3];
rz(-2.040401) q[3];
sx q[3];
rz(2.4503277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.80692446) q[2];
sx q[2];
rz(-0.88403264) q[2];
sx q[2];
rz(0.36726382) q[2];
rz(1.4372829) q[3];
sx q[3];
rz(-1.1742914) q[3];
sx q[3];
rz(-1.1737163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2631898) q[0];
sx q[0];
rz(-2.1387565) q[0];
sx q[0];
rz(-2.3314085) q[0];
rz(-1.1148249) q[1];
sx q[1];
rz(-2.7572542) q[1];
sx q[1];
rz(0.55327639) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.081757717) q[0];
sx q[0];
rz(-2.5011269) q[0];
sx q[0];
rz(-1.4778796) q[0];
rz(2.3478339) q[2];
sx q[2];
rz(-1.1053893) q[2];
sx q[2];
rz(2.2528354) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7825494) q[1];
sx q[1];
rz(-1.3842674) q[1];
sx q[1];
rz(-1.8576417) q[1];
rz(-pi) q[2];
rz(-1.9235885) q[3];
sx q[3];
rz(-0.81656352) q[3];
sx q[3];
rz(1.0366131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1222003) q[2];
sx q[2];
rz(-0.70851749) q[2];
sx q[2];
rz(-0.23966399) q[2];
rz(1.3278809) q[3];
sx q[3];
rz(-2.0769104) q[3];
sx q[3];
rz(0.46752587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0620621) q[0];
sx q[0];
rz(-1.7876328) q[0];
sx q[0];
rz(-2.7476516) q[0];
rz(2.486034) q[1];
sx q[1];
rz(-1.9656904) q[1];
sx q[1];
rz(0.94737731) q[1];
rz(1.9026359) q[2];
sx q[2];
rz(-1.9658026) q[2];
sx q[2];
rz(2.6323242) q[2];
rz(-0.37674698) q[3];
sx q[3];
rz(-2.527311) q[3];
sx q[3];
rz(-0.75025076) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
