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
rz(-2.2351216) q[0];
sx q[0];
rz(-1.6989166) q[0];
sx q[0];
rz(-0.28191167) q[0];
rz(-2.6126722) q[1];
sx q[1];
rz(-1.6493874) q[1];
sx q[1];
rz(1.563501) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53858763) q[0];
sx q[0];
rz(-1.9360124) q[0];
sx q[0];
rz(0.19293789) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.27484244) q[2];
sx q[2];
rz(-1.5759528) q[2];
sx q[2];
rz(2.8066563) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3696182) q[1];
sx q[1];
rz(-1.6965515) q[1];
sx q[1];
rz(-0.63196504) q[1];
rz(-1.2658843) q[3];
sx q[3];
rz(-0.90581363) q[3];
sx q[3];
rz(-0.18892442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6866744) q[2];
sx q[2];
rz(-0.29759559) q[2];
sx q[2];
rz(0.61304027) q[2];
rz(-2.6679299) q[3];
sx q[3];
rz(-1.9434171) q[3];
sx q[3];
rz(1.4325498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4050201) q[0];
sx q[0];
rz(-0.18017811) q[0];
sx q[0];
rz(0.75102425) q[0];
rz(2.660102) q[1];
sx q[1];
rz(-1.0844237) q[1];
sx q[1];
rz(-0.96985936) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78208047) q[0];
sx q[0];
rz(-2.1954698) q[0];
sx q[0];
rz(2.0332574) q[0];
rz(-pi) q[1];
rz(-0.2158176) q[2];
sx q[2];
rz(-1.6105798) q[2];
sx q[2];
rz(-1.424274) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0326421) q[1];
sx q[1];
rz(-1.3130929) q[1];
sx q[1];
rz(0.9089246) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0131575) q[3];
sx q[3];
rz(-2.1214888) q[3];
sx q[3];
rz(-0.1503508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2227309) q[2];
sx q[2];
rz(-1.4758045) q[2];
sx q[2];
rz(-0.53885031) q[2];
rz(-0.017596267) q[3];
sx q[3];
rz(-0.16418695) q[3];
sx q[3];
rz(2.2331451) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4410412) q[0];
sx q[0];
rz(-2.7563162) q[0];
sx q[0];
rz(0.061263099) q[0];
rz(1.6368658) q[1];
sx q[1];
rz(-0.69925362) q[1];
sx q[1];
rz(-1.1963074) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6590779) q[0];
sx q[0];
rz(-1.8097623) q[0];
sx q[0];
rz(2.2963797) q[0];
rz(-pi) q[1];
rz(1.5012791) q[2];
sx q[2];
rz(-1.3090773) q[2];
sx q[2];
rz(-1.4756605) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.34431008) q[1];
sx q[1];
rz(-0.6629262) q[1];
sx q[1];
rz(1.0271038) q[1];
rz(-pi) q[2];
x q[2];
rz(2.263271) q[3];
sx q[3];
rz(-1.3008144) q[3];
sx q[3];
rz(-2.6019271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4521744) q[2];
sx q[2];
rz(-1.3070561) q[2];
sx q[2];
rz(-0.43251953) q[2];
rz(-1.0906667) q[3];
sx q[3];
rz(-2.5011823) q[3];
sx q[3];
rz(2.045491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5890305) q[0];
sx q[0];
rz(-2.2046389) q[0];
sx q[0];
rz(0.20137782) q[0];
rz(2.6248113) q[1];
sx q[1];
rz(-0.38025451) q[1];
sx q[1];
rz(-0.94863272) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0364591) q[0];
sx q[0];
rz(-0.74315846) q[0];
sx q[0];
rz(2.0997597) q[0];
rz(2.8490813) q[2];
sx q[2];
rz(-0.88380948) q[2];
sx q[2];
rz(-0.77921898) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3740179) q[1];
sx q[1];
rz(-1.5841055) q[1];
sx q[1];
rz(-0.50434169) q[1];
rz(-0.16818436) q[3];
sx q[3];
rz(-2.6144989) q[3];
sx q[3];
rz(-2.9236562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2739233) q[2];
sx q[2];
rz(-0.40931585) q[2];
sx q[2];
rz(2.5841827) q[2];
rz(3.0607767) q[3];
sx q[3];
rz(-1.2405453) q[3];
sx q[3];
rz(1.935299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(1.7252561) q[0];
sx q[0];
rz(-1.6660322) q[0];
sx q[0];
rz(2.1112554) q[0];
rz(2.985785) q[1];
sx q[1];
rz(-1.067433) q[1];
sx q[1];
rz(-0.20419289) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4725114) q[0];
sx q[0];
rz(-1.8127499) q[0];
sx q[0];
rz(1.6415963) q[0];
rz(-pi) q[1];
rz(-1.7014808) q[2];
sx q[2];
rz(-1.7820662) q[2];
sx q[2];
rz(-1.4071161) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9928007) q[1];
sx q[1];
rz(-2.6084503) q[1];
sx q[1];
rz(0.41457446) q[1];
x q[2];
rz(2.89661) q[3];
sx q[3];
rz(-0.7881654) q[3];
sx q[3];
rz(0.80313166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.92945176) q[2];
sx q[2];
rz(-1.6484304) q[2];
sx q[2];
rz(-2.7280651) q[2];
rz(2.2996969) q[3];
sx q[3];
rz(-1.3827518) q[3];
sx q[3];
rz(2.1690185) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.095116422) q[0];
sx q[0];
rz(-1.0914509) q[0];
sx q[0];
rz(3.1360151) q[0];
rz(-1.3439517) q[1];
sx q[1];
rz(-0.19897142) q[1];
sx q[1];
rz(-0.70095789) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98788059) q[0];
sx q[0];
rz(-1.0544262) q[0];
sx q[0];
rz(-0.51306458) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2118819) q[2];
sx q[2];
rz(-1.2107009) q[2];
sx q[2];
rz(1.4907229) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8474104) q[1];
sx q[1];
rz(-0.27833101) q[1];
sx q[1];
rz(-2.9918482) q[1];
x q[2];
rz(2.088955) q[3];
sx q[3];
rz(-1.0598044) q[3];
sx q[3];
rz(2.4571153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8274716) q[2];
sx q[2];
rz(-0.33385971) q[2];
sx q[2];
rz(-1.1673048) q[2];
rz(-0.88268924) q[3];
sx q[3];
rz(-0.76789951) q[3];
sx q[3];
rz(0.32514969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90671396) q[0];
sx q[0];
rz(-2.0218847) q[0];
sx q[0];
rz(-3.0562905) q[0];
rz(0.75434297) q[1];
sx q[1];
rz(-0.5032379) q[1];
sx q[1];
rz(-2.1139961) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7639973) q[0];
sx q[0];
rz(-0.21276424) q[0];
sx q[0];
rz(-3.0932674) q[0];
rz(-pi) q[1];
rz(2.5449488) q[2];
sx q[2];
rz(-1.5663356) q[2];
sx q[2];
rz(-0.73847929) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.8366521) q[1];
sx q[1];
rz(-1.2479559) q[1];
sx q[1];
rz(0.12860997) q[1];
x q[2];
rz(0.87966921) q[3];
sx q[3];
rz(-2.4604049) q[3];
sx q[3];
rz(-0.78024125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.1520285) q[2];
sx q[2];
rz(-0.98735183) q[2];
sx q[2];
rz(2.7222471) q[2];
rz(2.5069405) q[3];
sx q[3];
rz(-2.3366163) q[3];
sx q[3];
rz(-1.1254719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80970508) q[0];
sx q[0];
rz(-0.87600791) q[0];
sx q[0];
rz(0.69212717) q[0];
rz(-0.57811111) q[1];
sx q[1];
rz(-2.7222241) q[1];
sx q[1];
rz(-1.7587761) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22631881) q[0];
sx q[0];
rz(-1.5513486) q[0];
sx q[0];
rz(1.6383199) q[0];
rz(0.67202576) q[2];
sx q[2];
rz(-2.1161072) q[2];
sx q[2];
rz(-1.3915075) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9005712) q[1];
sx q[1];
rz(-1.3374778) q[1];
sx q[1];
rz(-2.1395348) q[1];
rz(-pi) q[2];
x q[2];
rz(0.84660952) q[3];
sx q[3];
rz(-1.203152) q[3];
sx q[3];
rz(-1.6366307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.77070037) q[2];
sx q[2];
rz(-0.53052491) q[2];
sx q[2];
rz(-1.716506) q[2];
rz(-0.51618451) q[3];
sx q[3];
rz(-2.1631212) q[3];
sx q[3];
rz(1.2396575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6462964) q[0];
sx q[0];
rz(-1.7072059) q[0];
sx q[0];
rz(2.7700951) q[0];
rz(-0.83483541) q[1];
sx q[1];
rz(-2.3284262) q[1];
sx q[1];
rz(-0.017597839) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71401063) q[0];
sx q[0];
rz(-1.0508063) q[0];
sx q[0];
rz(-2.8218357) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8210635) q[2];
sx q[2];
rz(-2.2907718) q[2];
sx q[2];
rz(-1.7257476) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.086603348) q[1];
sx q[1];
rz(-1.2003683) q[1];
sx q[1];
rz(-2.5309674) q[1];
x q[2];
rz(-0.60815717) q[3];
sx q[3];
rz(-1.8377234) q[3];
sx q[3];
rz(1.1692695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1060433) q[2];
sx q[2];
rz(-2.273166) q[2];
sx q[2];
rz(3.0412728) q[2];
rz(0.22942461) q[3];
sx q[3];
rz(-2.094163) q[3];
sx q[3];
rz(2.6917698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4556274) q[0];
sx q[0];
rz(-0.28814155) q[0];
sx q[0];
rz(-0.57383865) q[0];
rz(-2.0876743) q[1];
sx q[1];
rz(-1.046448) q[1];
sx q[1];
rz(2.4651249) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5528111) q[0];
sx q[0];
rz(-1.8491321) q[0];
sx q[0];
rz(1.8235444) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0208548) q[2];
sx q[2];
rz(-0.46560198) q[2];
sx q[2];
rz(0.6731205) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.45470787) q[1];
sx q[1];
rz(-1.0446207) q[1];
sx q[1];
rz(1.4684907) q[1];
rz(-pi) q[2];
rz(-1.2419534) q[3];
sx q[3];
rz(-2.5848205) q[3];
sx q[3];
rz(2.013735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.68710589) q[2];
sx q[2];
rz(-2.3837619) q[2];
sx q[2];
rz(1.494361) q[2];
rz(-2.2729661) q[3];
sx q[3];
rz(-2.4354911) q[3];
sx q[3];
rz(-2.6715265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(0.070521991) q[0];
sx q[0];
rz(-1.949953) q[0];
sx q[0];
rz(-1.9198887) q[0];
rz(-1.8078049) q[1];
sx q[1];
rz(-1.318327) q[1];
sx q[1];
rz(-0.64073906) q[1];
rz(-3.0467544) q[2];
sx q[2];
rz(-0.76764501) q[2];
sx q[2];
rz(-2.1697247) q[2];
rz(-3.0972539) q[3];
sx q[3];
rz(-1.7462742) q[3];
sx q[3];
rz(-2.3901226) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
