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
rz(-1.6406887) q[0];
sx q[0];
rz(-2.357643) q[0];
sx q[0];
rz(0.10070237) q[0];
rz(-1.4251703) q[1];
sx q[1];
rz(-0.66882381) q[1];
sx q[1];
rz(1.5688904) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4948499) q[0];
sx q[0];
rz(-1.6892489) q[0];
sx q[0];
rz(-0.54991566) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4670829) q[2];
sx q[2];
rz(-1.3000037) q[2];
sx q[2];
rz(2.873604) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0364337) q[1];
sx q[1];
rz(-1.729125) q[1];
sx q[1];
rz(1.4033844) q[1];
rz(-2.334758) q[3];
sx q[3];
rz(-0.63043046) q[3];
sx q[3];
rz(0.85496074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0147741) q[2];
sx q[2];
rz(-2.7744881) q[2];
sx q[2];
rz(-1.611562) q[2];
rz(1.4989467) q[3];
sx q[3];
rz(-0.25224125) q[3];
sx q[3];
rz(0.69935316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.67795578) q[0];
sx q[0];
rz(-1.3709443) q[0];
sx q[0];
rz(0.26588765) q[0];
rz(-1.5222585) q[1];
sx q[1];
rz(-1.3109173) q[1];
sx q[1];
rz(-1.6552077) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9135654) q[0];
sx q[0];
rz(-1.5633213) q[0];
sx q[0];
rz(-2.3662503) q[0];
x q[1];
rz(2.5780074) q[2];
sx q[2];
rz(-1.3662837) q[2];
sx q[2];
rz(-0.16789189) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.34936541) q[1];
sx q[1];
rz(-2.0369895) q[1];
sx q[1];
rz(1.8833877) q[1];
rz(0.19274917) q[3];
sx q[3];
rz(-1.1067353) q[3];
sx q[3];
rz(1.0100067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8365525) q[2];
sx q[2];
rz(-2.063648) q[2];
sx q[2];
rz(-3.0744699) q[2];
rz(2.6302795) q[3];
sx q[3];
rz(-1.3498243) q[3];
sx q[3];
rz(2.2093723) q[3];
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
rz(-0.2387282) q[0];
sx q[0];
rz(-1.6599864) q[0];
sx q[0];
rz(0.1486775) q[0];
rz(0.45064926) q[1];
sx q[1];
rz(-1.5842178) q[1];
sx q[1];
rz(-2.2221036) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4212358) q[0];
sx q[0];
rz(-1.3034058) q[0];
sx q[0];
rz(-1.8084099) q[0];
x q[1];
rz(2.2581159) q[2];
sx q[2];
rz(-1.7224489) q[2];
sx q[2];
rz(0.54017457) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1053523) q[1];
sx q[1];
rz(-0.72004623) q[1];
sx q[1];
rz(-0.68036637) q[1];
rz(2.6387003) q[3];
sx q[3];
rz(-0.38030312) q[3];
sx q[3];
rz(1.8675493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.86842361) q[2];
sx q[2];
rz(-1.0852852) q[2];
sx q[2];
rz(-0.85624179) q[2];
rz(2.26561) q[3];
sx q[3];
rz(-0.3951422) q[3];
sx q[3];
rz(-1.9877079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7885389) q[0];
sx q[0];
rz(-1.5748698) q[0];
sx q[0];
rz(0.62623155) q[0];
rz(-1.7010472) q[1];
sx q[1];
rz(-2.3150621) q[1];
sx q[1];
rz(0.70762077) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0514446) q[0];
sx q[0];
rz(-2.2934545) q[0];
sx q[0];
rz(3.0279006) q[0];
x q[1];
rz(-1.2265251) q[2];
sx q[2];
rz(-0.81721899) q[2];
sx q[2];
rz(-0.57932779) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.172239) q[1];
sx q[1];
rz(-1.8392977) q[1];
sx q[1];
rz(1.4866465) q[1];
x q[2];
rz(1.510101) q[3];
sx q[3];
rz(-0.48572054) q[3];
sx q[3];
rz(2.751089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6669199) q[2];
sx q[2];
rz(-1.8644574) q[2];
sx q[2];
rz(-0.39371583) q[2];
rz(0.79931021) q[3];
sx q[3];
rz(-0.36472133) q[3];
sx q[3];
rz(1.5289615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9202775) q[0];
sx q[0];
rz(-2.3096313) q[0];
sx q[0];
rz(3.063391) q[0];
rz(0.68866628) q[1];
sx q[1];
rz(-0.93910256) q[1];
sx q[1];
rz(0.92855531) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9720042) q[0];
sx q[0];
rz(-1.0079591) q[0];
sx q[0];
rz(1.2006951) q[0];
rz(0.86805328) q[2];
sx q[2];
rz(-0.45449025) q[2];
sx q[2];
rz(2.8502685) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4457561) q[1];
sx q[1];
rz(-1.1923881) q[1];
sx q[1];
rz(1.3354785) q[1];
rz(-pi) q[2];
rz(2.6070091) q[3];
sx q[3];
rz(-2.4494736) q[3];
sx q[3];
rz(0.87833709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3277305) q[2];
sx q[2];
rz(-2.6321415) q[2];
sx q[2];
rz(0.65208411) q[2];
rz(-0.70513519) q[3];
sx q[3];
rz(-1.4950246) q[3];
sx q[3];
rz(0.42858538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9014277) q[0];
sx q[0];
rz(-2.9225898) q[0];
sx q[0];
rz(-1.3448311) q[0];
rz(0.52974686) q[1];
sx q[1];
rz(-1.2836645) q[1];
sx q[1];
rz(-1.8240428) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3968788) q[0];
sx q[0];
rz(-1.4900123) q[0];
sx q[0];
rz(2.1470444) q[0];
x q[1];
rz(-2.8805783) q[2];
sx q[2];
rz(-1.3512638) q[2];
sx q[2];
rz(-0.21071502) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.81281137) q[1];
sx q[1];
rz(-1.9095632) q[1];
sx q[1];
rz(0.12943204) q[1];
rz(-1.2413792) q[3];
sx q[3];
rz(-0.53295153) q[3];
sx q[3];
rz(-0.50607133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.162447) q[2];
sx q[2];
rz(-1.8628758) q[2];
sx q[2];
rz(-2.0549959) q[2];
rz(-3.047591) q[3];
sx q[3];
rz(-2.0407929) q[3];
sx q[3];
rz(-2.6619958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78733665) q[0];
sx q[0];
rz(-0.84699637) q[0];
sx q[0];
rz(0.081548668) q[0];
rz(1.5230007) q[1];
sx q[1];
rz(-1.0146419) q[1];
sx q[1];
rz(2.5340705) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1475246) q[0];
sx q[0];
rz(-0.77558625) q[0];
sx q[0];
rz(1.8344384) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6742769) q[2];
sx q[2];
rz(-1.226034) q[2];
sx q[2];
rz(-1.4490102) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5931892) q[1];
sx q[1];
rz(-1.4680937) q[1];
sx q[1];
rz(0.013506933) q[1];
rz(-pi) q[2];
x q[2];
rz(0.37346249) q[3];
sx q[3];
rz(-1.499553) q[3];
sx q[3];
rz(0.072865818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.088323204) q[2];
sx q[2];
rz(-0.20603267) q[2];
sx q[2];
rz(2.6982809) q[2];
rz(-0.045470227) q[3];
sx q[3];
rz(-1.1683522) q[3];
sx q[3];
rz(2.5209691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4745673) q[0];
sx q[0];
rz(-0.98282951) q[0];
sx q[0];
rz(2.5734651) q[0];
rz(1.2688295) q[1];
sx q[1];
rz(-1.9867691) q[1];
sx q[1];
rz(2.5297129) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6997212) q[0];
sx q[0];
rz(-1.4825025) q[0];
sx q[0];
rz(-2.7631212) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.9857432) q[2];
sx q[2];
rz(-0.32801706) q[2];
sx q[2];
rz(0.32250052) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5962741) q[1];
sx q[1];
rz(-2.5891212) q[1];
sx q[1];
rz(1.0956863) q[1];
rz(-pi) q[2];
rz(-1.7692206) q[3];
sx q[3];
rz(-1.7412609) q[3];
sx q[3];
rz(-0.19399079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2511217) q[2];
sx q[2];
rz(-1.2476363) q[2];
sx q[2];
rz(-0.20784155) q[2];
rz(-3.1244997) q[3];
sx q[3];
rz(-2.1199675) q[3];
sx q[3];
rz(-1.0926854) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6297778) q[0];
sx q[0];
rz(-0.20447525) q[0];
sx q[0];
rz(-1.7275607) q[0];
rz(3.0746225) q[1];
sx q[1];
rz(-1.3221062) q[1];
sx q[1];
rz(-3.0019143) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4614452) q[0];
sx q[0];
rz(-1.2530039) q[0];
sx q[0];
rz(-1.2387973) q[0];
rz(-pi) q[1];
rz(-1.1316179) q[2];
sx q[2];
rz(-2.1358523) q[2];
sx q[2];
rz(-1.4840839) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5433257) q[1];
sx q[1];
rz(-1.3451335) q[1];
sx q[1];
rz(-1.3545827) q[1];
x q[2];
rz(-3.090276) q[3];
sx q[3];
rz(-1.4017377) q[3];
sx q[3];
rz(-0.07758197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.70860538) q[2];
sx q[2];
rz(-1.2689509) q[2];
sx q[2];
rz(-2.404786) q[2];
rz(1.6144729) q[3];
sx q[3];
rz(-1.0090642) q[3];
sx q[3];
rz(1.8710776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0438185) q[0];
sx q[0];
rz(-2.173824) q[0];
sx q[0];
rz(-1.1678102) q[0];
rz(2.4763926) q[1];
sx q[1];
rz(-1.8635112) q[1];
sx q[1];
rz(2.1591689) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5031122) q[0];
sx q[0];
rz(-1.3706722) q[0];
sx q[0];
rz(-1.0699231) q[0];
rz(0.63913156) q[2];
sx q[2];
rz(-1.35372) q[2];
sx q[2];
rz(0.21737305) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3150627) q[1];
sx q[1];
rz(-0.86931397) q[1];
sx q[1];
rz(1.6531273) q[1];
rz(-1.7299557) q[3];
sx q[3];
rz(-2.7388262) q[3];
sx q[3];
rz(-1.7006766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3490225) q[2];
sx q[2];
rz(-2.4945365) q[2];
sx q[2];
rz(1.6046074) q[2];
rz(-1.0200621) q[3];
sx q[3];
rz(-1.5217109) q[3];
sx q[3];
rz(0.1608688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0237324) q[0];
sx q[0];
rz(-1.5845789) q[0];
sx q[0];
rz(-1.3316863) q[0];
rz(-2.3125519) q[1];
sx q[1];
rz(-1.4844898) q[1];
sx q[1];
rz(-0.97113562) q[1];
rz(1.2185417) q[2];
sx q[2];
rz(-2.2178219) q[2];
sx q[2];
rz(-2.1597663) q[2];
rz(2.0547413) q[3];
sx q[3];
rz(-1.9515358) q[3];
sx q[3];
rz(0.87072434) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
