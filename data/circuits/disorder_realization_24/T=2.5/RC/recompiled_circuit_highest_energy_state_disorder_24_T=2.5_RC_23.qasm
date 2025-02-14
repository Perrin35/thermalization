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
rz(-0.51198045) q[0];
sx q[0];
rz(-2.5706302) q[0];
sx q[0];
rz(-2.9538739) q[0];
rz(-2.3277148) q[1];
sx q[1];
rz(-1.4717646) q[1];
sx q[1];
rz(1.8522813) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1027689) q[0];
sx q[0];
rz(-1.7784405) q[0];
sx q[0];
rz(-0.24216346) q[0];
x q[1];
rz(1.600483) q[2];
sx q[2];
rz(-2.1412537) q[2];
sx q[2];
rz(-0.96889673) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6481406) q[1];
sx q[1];
rz(-1.3228184) q[1];
sx q[1];
rz(-0.47391717) q[1];
rz(-pi) q[2];
rz(1.9240407) q[3];
sx q[3];
rz(-1.3526256) q[3];
sx q[3];
rz(2.6835203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2671555) q[2];
sx q[2];
rz(-1.1031373) q[2];
sx q[2];
rz(-0.77679408) q[2];
rz(1.8393501) q[3];
sx q[3];
rz(-1.4812508) q[3];
sx q[3];
rz(1.9384025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5341107) q[0];
sx q[0];
rz(-0.48180875) q[0];
sx q[0];
rz(-0.99824655) q[0];
rz(-0.36969319) q[1];
sx q[1];
rz(-1.0414711) q[1];
sx q[1];
rz(1.1179771) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8000172) q[0];
sx q[0];
rz(-1.1693926) q[0];
sx q[0];
rz(0.97824162) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4684899) q[2];
sx q[2];
rz(-2.3323698) q[2];
sx q[2];
rz(0.26462091) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3433132) q[1];
sx q[1];
rz(-1.717853) q[1];
sx q[1];
rz(-1.1215854) q[1];
rz(2.8727628) q[3];
sx q[3];
rz(-2.1850366) q[3];
sx q[3];
rz(-1.5996831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.86576858) q[2];
sx q[2];
rz(-1.7663225) q[2];
sx q[2];
rz(2.7562874) q[2];
rz(-3.1028808) q[3];
sx q[3];
rz(-2.7888515) q[3];
sx q[3];
rz(-2.501131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5019219) q[0];
sx q[0];
rz(-2.6646035) q[0];
sx q[0];
rz(1.3013526) q[0];
rz(-0.057706984) q[1];
sx q[1];
rz(-0.4464018) q[1];
sx q[1];
rz(-1.3267964) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1088789) q[0];
sx q[0];
rz(-2.0427867) q[0];
sx q[0];
rz(2.1817939) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4313372) q[2];
sx q[2];
rz(-1.4118952) q[2];
sx q[2];
rz(0.1316084) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.88736278) q[1];
sx q[1];
rz(-2.0944203) q[1];
sx q[1];
rz(1.2382522) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7884198) q[3];
sx q[3];
rz(-0.88017094) q[3];
sx q[3];
rz(1.6817125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3441299) q[2];
sx q[2];
rz(-0.8351438) q[2];
sx q[2];
rz(0.75378913) q[2];
rz(1.3695184) q[3];
sx q[3];
rz(-1.4621719) q[3];
sx q[3];
rz(-2.4863825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9264483) q[0];
sx q[0];
rz(-1.0855874) q[0];
sx q[0];
rz(1.7468859) q[0];
rz(-1.9649547) q[1];
sx q[1];
rz(-1.5086915) q[1];
sx q[1];
rz(-0.36697695) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.079355571) q[0];
sx q[0];
rz(-1.5141762) q[0];
sx q[0];
rz(1.6376154) q[0];
rz(0.33432482) q[2];
sx q[2];
rz(-1.8211289) q[2];
sx q[2];
rz(1.804753) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5287197) q[1];
sx q[1];
rz(-2.7320128) q[1];
sx q[1];
rz(0.50480857) q[1];
x q[2];
rz(0.22032471) q[3];
sx q[3];
rz(-2.7661565) q[3];
sx q[3];
rz(-2.8130949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.41161141) q[2];
sx q[2];
rz(-1.9133762) q[2];
sx q[2];
rz(-1.270594) q[2];
rz(2.1805084) q[3];
sx q[3];
rz(-1.7226487) q[3];
sx q[3];
rz(-3.0911176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5807895) q[0];
sx q[0];
rz(-2.2152948) q[0];
sx q[0];
rz(1.8286937) q[0];
rz(-1.987223) q[1];
sx q[1];
rz(-1.1416953) q[1];
sx q[1];
rz(-2.5837574) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2656888) q[0];
sx q[0];
rz(-2.3703528) q[0];
sx q[0];
rz(-0.79637595) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6959582) q[2];
sx q[2];
rz(-2.9639396) q[2];
sx q[2];
rz(2.7522617) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3250998) q[1];
sx q[1];
rz(-1.268728) q[1];
sx q[1];
rz(-2.3343071) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6799503) q[3];
sx q[3];
rz(-1.7064693) q[3];
sx q[3];
rz(-1.5380579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4981726) q[2];
sx q[2];
rz(-1.8489445) q[2];
sx q[2];
rz(-0.84929973) q[2];
rz(0.15527209) q[3];
sx q[3];
rz(-0.97584358) q[3];
sx q[3];
rz(0.37080216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3575386) q[0];
sx q[0];
rz(-0.58670601) q[0];
sx q[0];
rz(2.6724755) q[0];
rz(-0.47438374) q[1];
sx q[1];
rz(-0.7862888) q[1];
sx q[1];
rz(-0.75538409) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2952399) q[0];
sx q[0];
rz(-1.8420514) q[0];
sx q[0];
rz(-3.1370107) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7003661) q[2];
sx q[2];
rz(-2.4640814) q[2];
sx q[2];
rz(-2.4063619) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5893822) q[1];
sx q[1];
rz(-0.24866074) q[1];
sx q[1];
rz(-1.0316487) q[1];
rz(-pi) q[2];
x q[2];
rz(0.83306082) q[3];
sx q[3];
rz(-1.2445881) q[3];
sx q[3];
rz(-1.2556094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0834121) q[2];
sx q[2];
rz(-1.8184793) q[2];
sx q[2];
rz(-2.9537436) q[2];
rz(1.8958873) q[3];
sx q[3];
rz(-1.3789504) q[3];
sx q[3];
rz(1.0904306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2096527) q[0];
sx q[0];
rz(-1.2670452) q[0];
sx q[0];
rz(0.39091045) q[0];
rz(0.93005013) q[1];
sx q[1];
rz(-1.7101733) q[1];
sx q[1];
rz(-1.3040868) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.025607312) q[0];
sx q[0];
rz(-1.5918333) q[0];
sx q[0];
rz(-1.2198835) q[0];
x q[1];
rz(-0.2371929) q[2];
sx q[2];
rz(-0.35065864) q[2];
sx q[2];
rz(-1.2192977) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.75445777) q[1];
sx q[1];
rz(-1.7384733) q[1];
sx q[1];
rz(-2.0281591) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5518858) q[3];
sx q[3];
rz(-1.867037) q[3];
sx q[3];
rz(-0.66413524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2685252) q[2];
sx q[2];
rz(-0.26669845) q[2];
sx q[2];
rz(0.11338691) q[2];
rz(3.125627) q[3];
sx q[3];
rz(-1.4957875) q[3];
sx q[3];
rz(-0.16460831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3575344) q[0];
sx q[0];
rz(-1.4736195) q[0];
sx q[0];
rz(-0.12406021) q[0];
rz(3.0198174) q[1];
sx q[1];
rz(-0.75141326) q[1];
sx q[1];
rz(1.442499) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0826113) q[0];
sx q[0];
rz(-1.9537203) q[0];
sx q[0];
rz(0.10848606) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8496672) q[2];
sx q[2];
rz(-1.5291844) q[2];
sx q[2];
rz(-2.6440563) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5296764) q[1];
sx q[1];
rz(-1.9135336) q[1];
sx q[1];
rz(3.0004098) q[1];
rz(-pi) q[2];
rz(-2.4115218) q[3];
sx q[3];
rz(-1.3833117) q[3];
sx q[3];
rz(1.2227675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.96559912) q[2];
sx q[2];
rz(-1.3600574) q[2];
sx q[2];
rz(2.3967801) q[2];
rz(-0.92721573) q[3];
sx q[3];
rz(-1.5085446) q[3];
sx q[3];
rz(-2.0492699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6532779) q[0];
sx q[0];
rz(-1.7019685) q[0];
sx q[0];
rz(-0.22536817) q[0];
rz(1.1124181) q[1];
sx q[1];
rz(-0.77443361) q[1];
sx q[1];
rz(-2.2427799) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32466896) q[0];
sx q[0];
rz(-0.82624871) q[0];
sx q[0];
rz(2.3284495) q[0];
rz(-1.1583011) q[2];
sx q[2];
rz(-1.0935244) q[2];
sx q[2];
rz(-2.0606747) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.73279335) q[1];
sx q[1];
rz(-0.91161722) q[1];
sx q[1];
rz(0.2562457) q[1];
x q[2];
rz(1.6055272) q[3];
sx q[3];
rz(-0.29303778) q[3];
sx q[3];
rz(1.422783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0972458) q[2];
sx q[2];
rz(-1.711742) q[2];
sx q[2];
rz(0.35935768) q[2];
rz(2.8906631) q[3];
sx q[3];
rz(-1.072262) q[3];
sx q[3];
rz(0.76771626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74442416) q[0];
sx q[0];
rz(-3.011062) q[0];
sx q[0];
rz(-2.2644444) q[0];
rz(-1.5769222) q[1];
sx q[1];
rz(-1.501187) q[1];
sx q[1];
rz(-2.3559779) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6495384) q[0];
sx q[0];
rz(-1.3527332) q[0];
sx q[0];
rz(-1.8137003) q[0];
rz(-0.53391407) q[2];
sx q[2];
rz(-1.1964239) q[2];
sx q[2];
rz(1.716734) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.12745107) q[1];
sx q[1];
rz(-1.2268865) q[1];
sx q[1];
rz(1.6456855) q[1];
rz(-pi) q[2];
rz(-0.50650017) q[3];
sx q[3];
rz(-1.7025547) q[3];
sx q[3];
rz(1.7817287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2239573) q[2];
sx q[2];
rz(-2.3515297) q[2];
sx q[2];
rz(0.7381953) q[2];
rz(-0.92292845) q[3];
sx q[3];
rz(-1.8332558) q[3];
sx q[3];
rz(-2.8452528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34687635) q[0];
sx q[0];
rz(-1.7807757) q[0];
sx q[0];
rz(-2.1697252) q[0];
rz(-2.234266) q[1];
sx q[1];
rz(-1.0846039) q[1];
sx q[1];
rz(-0.32122282) q[1];
rz(-2.1512866) q[2];
sx q[2];
rz(-1.9038426) q[2];
sx q[2];
rz(2.7833084) q[2];
rz(2.3774556) q[3];
sx q[3];
rz(-2.4485179) q[3];
sx q[3];
rz(-2.5701523) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
