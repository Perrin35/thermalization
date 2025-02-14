OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.174515) q[0];
sx q[0];
rz(-1.6835901) q[0];
sx q[0];
rz(-0.30012497) q[0];
rz(-1.9441654) q[1];
sx q[1];
rz(-1.5265042) q[1];
sx q[1];
rz(-1.6405029) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1688582) q[0];
sx q[0];
rz(-1.5669113) q[0];
sx q[0];
rz(-0.015104276) q[0];
rz(-pi) q[1];
rz(0.80279843) q[2];
sx q[2];
rz(-1.8539696) q[2];
sx q[2];
rz(2.6263023) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.38084322) q[1];
sx q[1];
rz(-0.92807209) q[1];
sx q[1];
rz(1.2973143) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5759104) q[3];
sx q[3];
rz(-0.76200125) q[3];
sx q[3];
rz(-0.68540547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7370721) q[2];
sx q[2];
rz(-1.9998735) q[2];
sx q[2];
rz(0.70297757) q[2];
rz(0.56973714) q[3];
sx q[3];
rz(-0.8786141) q[3];
sx q[3];
rz(-1.5841293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0074145929) q[0];
sx q[0];
rz(-0.61892048) q[0];
sx q[0];
rz(-1.9084357) q[0];
rz(0.59457072) q[1];
sx q[1];
rz(-2.1023127) q[1];
sx q[1];
rz(2.0436683) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1933408) q[0];
sx q[0];
rz(-1.7714714) q[0];
sx q[0];
rz(3.0737123) q[0];
rz(0.64867257) q[2];
sx q[2];
rz(-0.59174109) q[2];
sx q[2];
rz(1.6206738) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6095396) q[1];
sx q[1];
rz(-1.5455478) q[1];
sx q[1];
rz(-1.6482501) q[1];
rz(-pi) q[2];
rz(-3.0497562) q[3];
sx q[3];
rz(-2.3478234) q[3];
sx q[3];
rz(-0.61749279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7132831) q[2];
sx q[2];
rz(-1.2998591) q[2];
sx q[2];
rz(1.606288) q[2];
rz(0.71896583) q[3];
sx q[3];
rz(-0.9459559) q[3];
sx q[3];
rz(-2.9276221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(2.3035901) q[0];
sx q[0];
rz(-0.8135697) q[0];
sx q[0];
rz(2.549262) q[0];
rz(1.2447641) q[1];
sx q[1];
rz(-2.0015621) q[1];
sx q[1];
rz(-0.14561428) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1942822) q[0];
sx q[0];
rz(-0.5049476) q[0];
sx q[0];
rz(2.8584252) q[0];
rz(-pi) q[1];
rz(-2.2592741) q[2];
sx q[2];
rz(-1.6711418) q[2];
sx q[2];
rz(-1.060263) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.56359953) q[1];
sx q[1];
rz(-2.1195998) q[1];
sx q[1];
rz(0.33434681) q[1];
rz(2.0508025) q[3];
sx q[3];
rz(-2.601805) q[3];
sx q[3];
rz(1.9095498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.738872) q[2];
sx q[2];
rz(-2.1108184) q[2];
sx q[2];
rz(-1.2325475) q[2];
rz(-0.23972073) q[3];
sx q[3];
rz(-2.0870049) q[3];
sx q[3];
rz(-2.6565552) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0204912) q[0];
sx q[0];
rz(-0.70206577) q[0];
sx q[0];
rz(-3.1109911) q[0];
rz(-2.5864511) q[1];
sx q[1];
rz(-1.4786913) q[1];
sx q[1];
rz(1.0928924) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8776833) q[0];
sx q[0];
rz(-2.110024) q[0];
sx q[0];
rz(-2.3506021) q[0];
rz(-pi) q[1];
rz(3.080064) q[2];
sx q[2];
rz(-0.76632753) q[2];
sx q[2];
rz(1.7249223) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2177201) q[1];
sx q[1];
rz(-1.3663379) q[1];
sx q[1];
rz(-2.2108498) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3044358) q[3];
sx q[3];
rz(-2.2686696) q[3];
sx q[3];
rz(0.18750377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0838202) q[2];
sx q[2];
rz(-1.7744935) q[2];
sx q[2];
rz(0.27285451) q[2];
rz(-1.7504292) q[3];
sx q[3];
rz(-2.302156) q[3];
sx q[3];
rz(2.1896037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6735753) q[0];
sx q[0];
rz(-1.8033569) q[0];
sx q[0];
rz(1.8433628) q[0];
rz(-0.6908373) q[1];
sx q[1];
rz(-1.9419443) q[1];
sx q[1];
rz(-1.615049) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95193255) q[0];
sx q[0];
rz(-2.2554923) q[0];
sx q[0];
rz(-0.71345274) q[0];
rz(-pi) q[1];
rz(-1.2882907) q[2];
sx q[2];
rz(-1.2307248) q[2];
sx q[2];
rz(-1.1929389) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1880875) q[1];
sx q[1];
rz(-1.5449448) q[1];
sx q[1];
rz(0.8698747) q[1];
rz(-0.68256179) q[3];
sx q[3];
rz(-2.6877211) q[3];
sx q[3];
rz(1.0938494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4096628) q[2];
sx q[2];
rz(-0.80594984) q[2];
sx q[2];
rz(0.16732495) q[2];
rz(-1.3903728) q[3];
sx q[3];
rz(-2.6087587) q[3];
sx q[3];
rz(0.65756857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9996027) q[0];
sx q[0];
rz(-2.7775192) q[0];
sx q[0];
rz(-1.2784736) q[0];
rz(-1.7059884) q[1];
sx q[1];
rz(-1.6537138) q[1];
sx q[1];
rz(-2.7659168) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73248752) q[0];
sx q[0];
rz(-1.300749) q[0];
sx q[0];
rz(0.23464111) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2907545) q[2];
sx q[2];
rz(-2.9172643) q[2];
sx q[2];
rz(-1.4590603) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7573812) q[1];
sx q[1];
rz(-0.39296341) q[1];
sx q[1];
rz(-2.1860366) q[1];
x q[2];
rz(-2.385684) q[3];
sx q[3];
rz(-1.7497853) q[3];
sx q[3];
rz(0.32839113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0659236) q[2];
sx q[2];
rz(-2.0772987) q[2];
sx q[2];
rz(2.2806878) q[2];
rz(1.9979477) q[3];
sx q[3];
rz(-0.30503169) q[3];
sx q[3];
rz(0.27826571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.025573108) q[0];
sx q[0];
rz(-1.3207859) q[0];
sx q[0];
rz(-2.386911) q[0];
rz(-2.2017551) q[1];
sx q[1];
rz(-0.87516963) q[1];
sx q[1];
rz(1.8584724) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96961731) q[0];
sx q[0];
rz(-1.8775058) q[0];
sx q[0];
rz(0.0071010751) q[0];
rz(-pi) q[1];
rz(-0.91249429) q[2];
sx q[2];
rz(-2.5844838) q[2];
sx q[2];
rz(-2.0832286) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2055473) q[1];
sx q[1];
rz(-0.48921555) q[1];
sx q[1];
rz(-0.65566109) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96638443) q[3];
sx q[3];
rz(-2.3570286) q[3];
sx q[3];
rz(-1.094629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.54520404) q[2];
sx q[2];
rz(-2.0241006) q[2];
sx q[2];
rz(1.0835353) q[2];
rz(-0.74294535) q[3];
sx q[3];
rz(-1.6094306) q[3];
sx q[3];
rz(-1.4325745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
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
rz(2.1525986) q[0];
sx q[0];
rz(-2.3567663) q[0];
sx q[0];
rz(-2.3866744) q[0];
rz(-2.6424291) q[1];
sx q[1];
rz(-2.1776336) q[1];
sx q[1];
rz(-2.4430433) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7344727) q[0];
sx q[0];
rz(-1.2615347) q[0];
sx q[0];
rz(-2.5109108) q[0];
x q[1];
rz(-0.28251799) q[2];
sx q[2];
rz(-0.48136863) q[2];
sx q[2];
rz(-0.11061874) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0341728) q[1];
sx q[1];
rz(-2.1028958) q[1];
sx q[1];
rz(3.1183395) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5732865) q[3];
sx q[3];
rz(-0.41908136) q[3];
sx q[3];
rz(-0.97657874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2373206) q[2];
sx q[2];
rz(-0.99232435) q[2];
sx q[2];
rz(-0.80238706) q[2];
rz(2.5896416) q[3];
sx q[3];
rz(-0.80058432) q[3];
sx q[3];
rz(2.5210181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17860831) q[0];
sx q[0];
rz(-1.7010138) q[0];
sx q[0];
rz(-1.1075903) q[0];
rz(-1.7604609) q[1];
sx q[1];
rz(-1.3637204) q[1];
sx q[1];
rz(1.6345056) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7716146) q[0];
sx q[0];
rz(-0.42015892) q[0];
sx q[0];
rz(2.7922947) q[0];
rz(-pi) q[1];
x q[1];
rz(0.31839026) q[2];
sx q[2];
rz(-1.9620775) q[2];
sx q[2];
rz(1.1423781) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9025165) q[1];
sx q[1];
rz(-2.7488144) q[1];
sx q[1];
rz(0.3879438) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.082066925) q[3];
sx q[3];
rz(-0.82857271) q[3];
sx q[3];
rz(0.16211331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0845118) q[2];
sx q[2];
rz(-0.19597404) q[2];
sx q[2];
rz(1.2723119) q[2];
rz(-1.48014) q[3];
sx q[3];
rz(-2.0062165) q[3];
sx q[3];
rz(-0.60012668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24499527) q[0];
sx q[0];
rz(-2.5801881) q[0];
sx q[0];
rz(1.2835314) q[0];
rz(-0.01783477) q[1];
sx q[1];
rz(-0.63648883) q[1];
sx q[1];
rz(-3.0160115) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3572306) q[0];
sx q[0];
rz(-1.0166369) q[0];
sx q[0];
rz(-3.1249376) q[0];
rz(-pi) q[1];
x q[1];
rz(0.56350817) q[2];
sx q[2];
rz(-2.3803408) q[2];
sx q[2];
rz(-0.13680563) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4976313) q[1];
sx q[1];
rz(-2.5217274) q[1];
sx q[1];
rz(-2.3513112) q[1];
rz(-1.3036895) q[3];
sx q[3];
rz(-1.7648762) q[3];
sx q[3];
rz(-1.1052856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3632726) q[2];
sx q[2];
rz(-1.8724226) q[2];
sx q[2];
rz(-0.8030836) q[2];
rz(-2.965029) q[3];
sx q[3];
rz(-1.3253788) q[3];
sx q[3];
rz(2.363502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17978996) q[0];
sx q[0];
rz(-0.50146865) q[0];
sx q[0];
rz(1.5487221) q[0];
rz(1.322562) q[1];
sx q[1];
rz(-1.3643199) q[1];
sx q[1];
rz(1.551569) q[1];
rz(0.24257913) q[2];
sx q[2];
rz(-1.6318113) q[2];
sx q[2];
rz(3.0616374) q[2];
rz(1.5104891) q[3];
sx q[3];
rz(-1.2338439) q[3];
sx q[3];
rz(-3.0046786) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
