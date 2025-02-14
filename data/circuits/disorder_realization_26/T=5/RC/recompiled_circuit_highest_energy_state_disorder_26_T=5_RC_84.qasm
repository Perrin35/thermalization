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
rz(-2.969279) q[0];
sx q[0];
rz(-0.20792374) q[0];
sx q[0];
rz(-1.2907668) q[0];
rz(-2.9895904) q[1];
sx q[1];
rz(-0.64270371) q[1];
sx q[1];
rz(2.7132577) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59574678) q[0];
sx q[0];
rz(-1.4662678) q[0];
sx q[0];
rz(-1.680611) q[0];
rz(2.7784155) q[2];
sx q[2];
rz(-2.4360949) q[2];
sx q[2];
rz(0.35340912) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5338057) q[1];
sx q[1];
rz(-1.6139107) q[1];
sx q[1];
rz(2.4030466) q[1];
rz(0.21875225) q[3];
sx q[3];
rz(-1.3748589) q[3];
sx q[3];
rz(-2.4234555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9316445) q[2];
sx q[2];
rz(-2.1817709) q[2];
sx q[2];
rz(1.2098562) q[2];
rz(0.079553902) q[3];
sx q[3];
rz(-0.39593655) q[3];
sx q[3];
rz(-2.615926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(-2.2570268) q[0];
sx q[0];
rz(-0.50978065) q[0];
sx q[0];
rz(-0.38401815) q[0];
rz(0.89029038) q[1];
sx q[1];
rz(-1.9046116) q[1];
sx q[1];
rz(-0.30337897) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6174406) q[0];
sx q[0];
rz(-0.57311237) q[0];
sx q[0];
rz(-1.9277431) q[0];
rz(-pi) q[1];
x q[1];
rz(1.134642) q[2];
sx q[2];
rz(-0.56098962) q[2];
sx q[2];
rz(-0.84535384) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.56306707) q[1];
sx q[1];
rz(-1.6614698) q[1];
sx q[1];
rz(-0.062460047) q[1];
rz(-pi) q[2];
rz(2.0201265) q[3];
sx q[3];
rz(-0.96689618) q[3];
sx q[3];
rz(-2.3788798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.38637912) q[2];
sx q[2];
rz(-2.0342125) q[2];
sx q[2];
rz(-0.73491043) q[2];
rz(0.84878659) q[3];
sx q[3];
rz(-2.3990192) q[3];
sx q[3];
rz(-0.36062226) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81247771) q[0];
sx q[0];
rz(-0.37549967) q[0];
sx q[0];
rz(3.062881) q[0];
rz(0.82376897) q[1];
sx q[1];
rz(-2.0562101) q[1];
sx q[1];
rz(-1.2944006) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9322008) q[0];
sx q[0];
rz(-3.0552312) q[0];
sx q[0];
rz(-2.3461653) q[0];
rz(-pi) q[1];
rz(-2.2357777) q[2];
sx q[2];
rz(-0.73630263) q[2];
sx q[2];
rz(-1.6396963) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.30101206) q[1];
sx q[1];
rz(-2.8361179) q[1];
sx q[1];
rz(2.7423647) q[1];
rz(2.2749316) q[3];
sx q[3];
rz(-1.7539644) q[3];
sx q[3];
rz(-2.4742692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2417629) q[2];
sx q[2];
rz(-2.7624942) q[2];
sx q[2];
rz(2.3251593) q[2];
rz(1.624931) q[3];
sx q[3];
rz(-0.95976019) q[3];
sx q[3];
rz(0.70037705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4888332) q[0];
sx q[0];
rz(-2.3214898) q[0];
sx q[0];
rz(2.5732727) q[0];
rz(-1.9991416) q[1];
sx q[1];
rz(-1.3520974) q[1];
sx q[1];
rz(-1.4521339) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7441918) q[0];
sx q[0];
rz(-2.4699587) q[0];
sx q[0];
rz(0.36104843) q[0];
x q[1];
rz(1.6644888) q[2];
sx q[2];
rz(-1.7498651) q[2];
sx q[2];
rz(-0.25824091) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8206222) q[1];
sx q[1];
rz(-1.8838716) q[1];
sx q[1];
rz(2.2105107) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9630757) q[3];
sx q[3];
rz(-1.2367289) q[3];
sx q[3];
rz(-2.8802383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6292754) q[2];
sx q[2];
rz(-2.6111111) q[2];
sx q[2];
rz(3.0647035) q[2];
rz(2.7422089) q[3];
sx q[3];
rz(-2.2279584) q[3];
sx q[3];
rz(-2.5975749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9870616) q[0];
sx q[0];
rz(-0.35683826) q[0];
sx q[0];
rz(-2.8416908) q[0];
rz(1.1497644) q[1];
sx q[1];
rz(-1.0118142) q[1];
sx q[1];
rz(2.6531175) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4048674) q[0];
sx q[0];
rz(-2.953385) q[0];
sx q[0];
rz(-1.3948649) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.68814338) q[2];
sx q[2];
rz(-2.4240085) q[2];
sx q[2];
rz(-0.46844278) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9087968) q[1];
sx q[1];
rz(-2.1543062) q[1];
sx q[1];
rz(-1.5934492) q[1];
rz(-0.64170047) q[3];
sx q[3];
rz(-1.3535168) q[3];
sx q[3];
rz(-3.1171796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.38146314) q[2];
sx q[2];
rz(-0.13801408) q[2];
sx q[2];
rz(1.6628954) q[2];
rz(1.1100769) q[3];
sx q[3];
rz(-2.1434982) q[3];
sx q[3];
rz(-2.4983675) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4019302) q[0];
sx q[0];
rz(-1.1818385) q[0];
sx q[0];
rz(2.8231743) q[0];
rz(0.034612522) q[1];
sx q[1];
rz(-0.56762677) q[1];
sx q[1];
rz(2.6908223) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.634003) q[0];
sx q[0];
rz(-1.5566155) q[0];
sx q[0];
rz(-2.0812278) q[0];
rz(-pi) q[1];
rz(-2.2658453) q[2];
sx q[2];
rz(-1.0330457) q[2];
sx q[2];
rz(1.328383) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.83799441) q[1];
sx q[1];
rz(-0.7571836) q[1];
sx q[1];
rz(-2.1794469) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1799906) q[3];
sx q[3];
rz(-0.79422659) q[3];
sx q[3];
rz(-0.022217928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.047711756) q[2];
sx q[2];
rz(-2.9517951) q[2];
sx q[2];
rz(3.0799358) q[2];
rz(-0.098585248) q[3];
sx q[3];
rz(-2.3969789) q[3];
sx q[3];
rz(-2.4390167) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3848569) q[0];
sx q[0];
rz(-2.0282133) q[0];
sx q[0];
rz(3.1016438) q[0];
rz(0.7705676) q[1];
sx q[1];
rz(-2.4295085) q[1];
sx q[1];
rz(-2.9133453) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21439274) q[0];
sx q[0];
rz(-1.6480371) q[0];
sx q[0];
rz(1.4820251) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1399916) q[2];
sx q[2];
rz(-1.2646156) q[2];
sx q[2];
rz(-1.209335) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3747663) q[1];
sx q[1];
rz(-3.0438381) q[1];
sx q[1];
rz(-2.5663816) q[1];
rz(-0.93668117) q[3];
sx q[3];
rz(-1.2290658) q[3];
sx q[3];
rz(0.069119819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0780636) q[2];
sx q[2];
rz(-2.0827796) q[2];
sx q[2];
rz(-2.883319) q[2];
rz(2.9410948) q[3];
sx q[3];
rz(-2.343488) q[3];
sx q[3];
rz(0.18331461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16266009) q[0];
sx q[0];
rz(-1.7698092) q[0];
sx q[0];
rz(0.3072511) q[0];
rz(1.2112674) q[1];
sx q[1];
rz(-0.4937506) q[1];
sx q[1];
rz(2.5603851) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4330313) q[0];
sx q[0];
rz(-1.975346) q[0];
sx q[0];
rz(-3.1076215) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.18980726) q[2];
sx q[2];
rz(-2.4229089) q[2];
sx q[2];
rz(-0.41658066) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.30347428) q[1];
sx q[1];
rz(-1.7707157) q[1];
sx q[1];
rz(2.3568627) q[1];
x q[2];
rz(-0.56116207) q[3];
sx q[3];
rz(-0.61120874) q[3];
sx q[3];
rz(1.2964013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6792949) q[2];
sx q[2];
rz(-2.0079948) q[2];
sx q[2];
rz(3.1104258) q[2];
rz(-2.9160685) q[3];
sx q[3];
rz(-1.7799107) q[3];
sx q[3];
rz(-1.0303191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.088260055) q[0];
sx q[0];
rz(-3.0971165) q[0];
sx q[0];
rz(-2.4326676) q[0];
rz(-2.9290579) q[1];
sx q[1];
rz(-2.1268851) q[1];
sx q[1];
rz(2.2841891) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1171061) q[0];
sx q[0];
rz(-3.0117234) q[0];
sx q[0];
rz(1.5316567) q[0];
rz(-pi) q[1];
rz(1.2456263) q[2];
sx q[2];
rz(-2.0813) q[2];
sx q[2];
rz(0.73869642) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.45891201) q[1];
sx q[1];
rz(-2.4574728) q[1];
sx q[1];
rz(-2.4980809) q[1];
x q[2];
rz(-0.25257381) q[3];
sx q[3];
rz(-2.4285304) q[3];
sx q[3];
rz(0.79914645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7015486) q[2];
sx q[2];
rz(-0.35228071) q[2];
sx q[2];
rz(-2.3360543) q[2];
rz(-2.7696179) q[3];
sx q[3];
rz(-1.493908) q[3];
sx q[3];
rz(-0.39723799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1086248) q[0];
sx q[0];
rz(-1.2297577) q[0];
sx q[0];
rz(-2.1627872) q[0];
rz(0.72273123) q[1];
sx q[1];
rz(-1.1284072) q[1];
sx q[1];
rz(-0.61789787) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4198087) q[0];
sx q[0];
rz(-0.81994826) q[0];
sx q[0];
rz(0.3799812) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5064836) q[2];
sx q[2];
rz(-0.76580566) q[2];
sx q[2];
rz(-0.85585153) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4010716) q[1];
sx q[1];
rz(-0.24722543) q[1];
sx q[1];
rz(2.3257491) q[1];
rz(0.79094751) q[3];
sx q[3];
rz(-0.11868782) q[3];
sx q[3];
rz(0.044767901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2021947) q[2];
sx q[2];
rz(-1.3556182) q[2];
sx q[2];
rz(0.00016577684) q[2];
rz(0.58445066) q[3];
sx q[3];
rz(-1.0060468) q[3];
sx q[3];
rz(-0.59529006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7547739) q[0];
sx q[0];
rz(-1.5703572) q[0];
sx q[0];
rz(1.5686709) q[0];
rz(1.3407002) q[1];
sx q[1];
rz(-2.0410213) q[1];
sx q[1];
rz(-1.6165728) q[1];
rz(2.2095815) q[2];
sx q[2];
rz(-1.1975653) q[2];
sx q[2];
rz(-1.4690659) q[2];
rz(-0.44092785) q[3];
sx q[3];
rz(-1.2444513) q[3];
sx q[3];
rz(1.2969482) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
