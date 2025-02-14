OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.48489269) q[0];
sx q[0];
rz(-3.0527053) q[0];
sx q[0];
rz(0.77699295) q[0];
rz(0.69500336) q[1];
sx q[1];
rz(-3.0106795) q[1];
sx q[1];
rz(-2.793269) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94960864) q[0];
sx q[0];
rz(-2.3460547) q[0];
sx q[0];
rz(-0.27792161) q[0];
x q[1];
rz(2.8856002) q[2];
sx q[2];
rz(-0.29689327) q[2];
sx q[2];
rz(0.96742899) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4433044) q[1];
sx q[1];
rz(-1.0545679) q[1];
sx q[1];
rz(2.1702532) q[1];
rz(-pi) q[2];
rz(1.6197085) q[3];
sx q[3];
rz(-1.1756611) q[3];
sx q[3];
rz(-1.6312158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6473039) q[2];
sx q[2];
rz(-2.5755197) q[2];
sx q[2];
rz(1.0757793) q[2];
rz(0.11476573) q[3];
sx q[3];
rz(-1.4817295) q[3];
sx q[3];
rz(-2.4858294) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39605108) q[0];
sx q[0];
rz(-0.93463722) q[0];
sx q[0];
rz(2.4387687) q[0];
rz(1.5463691) q[1];
sx q[1];
rz(-2.399235) q[1];
sx q[1];
rz(0.57483086) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.069734978) q[0];
sx q[0];
rz(-1.1848579) q[0];
sx q[0];
rz(-1.925408) q[0];
rz(-pi) q[1];
rz(2.5288519) q[2];
sx q[2];
rz(-2.0144955) q[2];
sx q[2];
rz(-0.24335621) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0965754) q[1];
sx q[1];
rz(-0.63961601) q[1];
sx q[1];
rz(2.9208697) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4739875) q[3];
sx q[3];
rz(-0.75112334) q[3];
sx q[3];
rz(2.3977007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.90291643) q[2];
sx q[2];
rz(-1.6409589) q[2];
sx q[2];
rz(2.404876) q[2];
rz(-0.77218562) q[3];
sx q[3];
rz(-0.18852791) q[3];
sx q[3];
rz(2.7948715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18948983) q[0];
sx q[0];
rz(-1.075407) q[0];
sx q[0];
rz(-1.1601675) q[0];
rz(0.53933764) q[1];
sx q[1];
rz(-0.62349206) q[1];
sx q[1];
rz(-3.0711874) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7382875) q[0];
sx q[0];
rz(-2.0576446) q[0];
sx q[0];
rz(-2.3281718) q[0];
x q[1];
rz(-2.3023476) q[2];
sx q[2];
rz(-2.9021724) q[2];
sx q[2];
rz(0.024062238) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.16500948) q[1];
sx q[1];
rz(-1.9515459) q[1];
sx q[1];
rz(-1.6409207) q[1];
rz(-pi) q[2];
rz(2.539864) q[3];
sx q[3];
rz(-0.7454819) q[3];
sx q[3];
rz(-1.1423542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3714441) q[2];
sx q[2];
rz(-1.4263209) q[2];
sx q[2];
rz(2.9798129) q[2];
rz(-2.3700304) q[3];
sx q[3];
rz(-3.0050889) q[3];
sx q[3];
rz(-2.0989044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5914087) q[0];
sx q[0];
rz(-2.3593481) q[0];
sx q[0];
rz(-2.8915306) q[0];
rz(-1.0391883) q[1];
sx q[1];
rz(-2.3999374) q[1];
sx q[1];
rz(-0.80113062) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.259819) q[0];
sx q[0];
rz(-1.7823813) q[0];
sx q[0];
rz(1.5863281) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3336557) q[2];
sx q[2];
rz(-1.3190499) q[2];
sx q[2];
rz(-0.27790127) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.114157) q[1];
sx q[1];
rz(-2.0832169) q[1];
sx q[1];
rz(1.4911265) q[1];
rz(1.2053554) q[3];
sx q[3];
rz(-1.9840999) q[3];
sx q[3];
rz(-0.45729056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.73418009) q[2];
sx q[2];
rz(-1.4449747) q[2];
sx q[2];
rz(-0.59047353) q[2];
rz(-2.7967795) q[3];
sx q[3];
rz(-2.7761288) q[3];
sx q[3];
rz(-1.4709681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.19584) q[0];
sx q[0];
rz(-1.7029637) q[0];
sx q[0];
rz(1.6666743) q[0];
rz(-1.6189812) q[1];
sx q[1];
rz(-1.5983862) q[1];
sx q[1];
rz(-0.408907) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5602138) q[0];
sx q[0];
rz(-1.5239851) q[0];
sx q[0];
rz(1.5623832) q[0];
rz(-pi) q[1];
rz(-2.2719745) q[2];
sx q[2];
rz(-1.8843576) q[2];
sx q[2];
rz(-2.8558848) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.58432128) q[1];
sx q[1];
rz(-1.9964594) q[1];
sx q[1];
rz(-1.0866648) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0237883) q[3];
sx q[3];
rz(-2.4642053) q[3];
sx q[3];
rz(2.7927674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5976379) q[2];
sx q[2];
rz(-3.0089617) q[2];
sx q[2];
rz(2.4908716) q[2];
rz(1.6808274) q[3];
sx q[3];
rz(-1.4131578) q[3];
sx q[3];
rz(-2.7248603) q[3];
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
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8553218) q[0];
sx q[0];
rz(-2.0725508) q[0];
sx q[0];
rz(1.2438114) q[0];
rz(-2.6142201) q[1];
sx q[1];
rz(-1.4445883) q[1];
sx q[1];
rz(-1.3810371) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0422538) q[0];
sx q[0];
rz(-1.5586133) q[0];
sx q[0];
rz(2.8996917) q[0];
x q[1];
rz(-0.63746286) q[2];
sx q[2];
rz(-0.12345498) q[2];
sx q[2];
rz(0.055094624) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.11070793) q[1];
sx q[1];
rz(-2.0763016) q[1];
sx q[1];
rz(2.5896789) q[1];
rz(0.88532577) q[3];
sx q[3];
rz(-2.4866085) q[3];
sx q[3];
rz(1.8670435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1018299) q[2];
sx q[2];
rz(-1.7888174) q[2];
sx q[2];
rz(1.9977894) q[2];
rz(2.3666798) q[3];
sx q[3];
rz(-1.1760005) q[3];
sx q[3];
rz(0.2253783) q[3];
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
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5483122) q[0];
sx q[0];
rz(-0.53638023) q[0];
sx q[0];
rz(0.29658741) q[0];
rz(2.2071154) q[1];
sx q[1];
rz(-0.88948932) q[1];
sx q[1];
rz(2.8737822) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64235425) q[0];
sx q[0];
rz(-1.3569512) q[0];
sx q[0];
rz(-1.7876028) q[0];
rz(-pi) q[1];
rz(2.8765062) q[2];
sx q[2];
rz(-0.62787442) q[2];
sx q[2];
rz(-1.5292847) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6567475) q[1];
sx q[1];
rz(-2.2094927) q[1];
sx q[1];
rz(-1.188422) q[1];
x q[2];
rz(-0.15800627) q[3];
sx q[3];
rz(-1.4337162) q[3];
sx q[3];
rz(2.7003845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3725793) q[2];
sx q[2];
rz(-1.0472426) q[2];
sx q[2];
rz(-0.11865842) q[2];
rz(-0.81698051) q[3];
sx q[3];
rz(-1.9575565) q[3];
sx q[3];
rz(2.584804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7939222) q[0];
sx q[0];
rz(-2.6642647) q[0];
sx q[0];
rz(2.0678066) q[0];
rz(-2.1444164) q[1];
sx q[1];
rz(-1.8938096) q[1];
sx q[1];
rz(1.006385) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0050765) q[0];
sx q[0];
rz(-3.0126326) q[0];
sx q[0];
rz(1.5004679) q[0];
rz(-pi) q[1];
rz(0.7147737) q[2];
sx q[2];
rz(-0.48764418) q[2];
sx q[2];
rz(-0.89370103) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.83676618) q[1];
sx q[1];
rz(-1.8988892) q[1];
sx q[1];
rz(1.7734709) q[1];
rz(-pi) q[2];
rz(2.2131165) q[3];
sx q[3];
rz(-1.9657502) q[3];
sx q[3];
rz(-2.5028281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5239079) q[2];
sx q[2];
rz(-1.074581) q[2];
sx q[2];
rz(-0.96134821) q[2];
rz(-2.5210467) q[3];
sx q[3];
rz(-0.70182645) q[3];
sx q[3];
rz(2.3294241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0945702) q[0];
sx q[0];
rz(-2.4235239) q[0];
sx q[0];
rz(-2.4156003) q[0];
rz(-2.831931) q[1];
sx q[1];
rz(-1.0331215) q[1];
sx q[1];
rz(-2.9594701) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34946467) q[0];
sx q[0];
rz(-0.82577151) q[0];
sx q[0];
rz(-1.4057805) q[0];
rz(-pi) q[1];
rz(-0.51515963) q[2];
sx q[2];
rz(-2.023989) q[2];
sx q[2];
rz(0.31441385) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.88090522) q[1];
sx q[1];
rz(-1.1662848) q[1];
sx q[1];
rz(-1.4244367) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1624231) q[3];
sx q[3];
rz(-0.97399903) q[3];
sx q[3];
rz(-2.2788871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9228268) q[2];
sx q[2];
rz(-1.4309692) q[2];
sx q[2];
rz(-1.9975086) q[2];
rz(0.73793441) q[3];
sx q[3];
rz(-2.4149371) q[3];
sx q[3];
rz(2.0209594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4436307) q[0];
sx q[0];
rz(-1.8101036) q[0];
sx q[0];
rz(-0.6829845) q[0];
rz(1.5890652) q[1];
sx q[1];
rz(-2.2582159) q[1];
sx q[1];
rz(-0.89673269) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67499298) q[0];
sx q[0];
rz(-1.7651801) q[0];
sx q[0];
rz(2.7036576) q[0];
rz(-pi) q[1];
rz(2.1188154) q[2];
sx q[2];
rz(-1.6923133) q[2];
sx q[2];
rz(-2.1650328) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0005714) q[1];
sx q[1];
rz(-1.440417) q[1];
sx q[1];
rz(-1.6597553) q[1];
rz(-0.50161485) q[3];
sx q[3];
rz(-1.1715495) q[3];
sx q[3];
rz(0.58540067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.54301846) q[2];
sx q[2];
rz(-1.2335346) q[2];
sx q[2];
rz(-2.7198071) q[2];
rz(0.56454033) q[3];
sx q[3];
rz(-2.3081503) q[3];
sx q[3];
rz(-0.88258266) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17557872) q[0];
sx q[0];
rz(-1.5400374) q[0];
sx q[0];
rz(-1.2651428) q[0];
rz(1.9059146) q[1];
sx q[1];
rz(-1.0503678) q[1];
sx q[1];
rz(-2.8138524) q[1];
rz(-0.50516456) q[2];
sx q[2];
rz(-0.41675962) q[2];
sx q[2];
rz(2.6909697) q[2];
rz(2.6648585) q[3];
sx q[3];
rz(-1.1179964) q[3];
sx q[3];
rz(-2.6159058) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
