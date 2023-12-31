OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5053951) q[0];
sx q[0];
rz(-2.8656821) q[0];
sx q[0];
rz(-1.3077868) q[0];
rz(-2.0055327) q[1];
sx q[1];
rz(4.0772822) q[1];
sx q[1];
rz(4.7128591) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9309064) q[0];
sx q[0];
rz(-1.2331729) q[0];
sx q[0];
rz(0.36436413) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.76531305) q[2];
sx q[2];
rz(-1.0368477) q[2];
sx q[2];
rz(0.050616654) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.03304122) q[1];
sx q[1];
rz(-1.812495) q[1];
sx q[1];
rz(-0.34717314) q[1];
rz(0.19574638) q[3];
sx q[3];
rz(-1.2139075) q[3];
sx q[3];
rz(1.2008047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2661665) q[2];
sx q[2];
rz(-0.29310075) q[2];
sx q[2];
rz(-1.1323294) q[2];
rz(-1.4663565) q[3];
sx q[3];
rz(-1.8050067) q[3];
sx q[3];
rz(-1.0124538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.19673008) q[0];
sx q[0];
rz(-0.20962993) q[0];
sx q[0];
rz(-2.9557513) q[0];
rz(-0.56022412) q[1];
sx q[1];
rz(-1.8461684) q[1];
sx q[1];
rz(-0.21683189) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6235979) q[0];
sx q[0];
rz(-1.2835842) q[0];
sx q[0];
rz(0.90201305) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47302834) q[2];
sx q[2];
rz(-1.066726) q[2];
sx q[2];
rz(-2.8216528) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9456957) q[1];
sx q[1];
rz(-2.2356114) q[1];
sx q[1];
rz(2.871454) q[1];
x q[2];
rz(-2.300755) q[3];
sx q[3];
rz(-0.84078046) q[3];
sx q[3];
rz(-2.9527612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.310114) q[2];
sx q[2];
rz(-0.82565132) q[2];
sx q[2];
rz(1.2878093) q[2];
rz(-2.3790322) q[3];
sx q[3];
rz(-1.972714) q[3];
sx q[3];
rz(-0.30502239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4644311) q[0];
sx q[0];
rz(-0.34496775) q[0];
sx q[0];
rz(0.60423869) q[0];
rz(1.3263946) q[1];
sx q[1];
rz(-1.3605958) q[1];
sx q[1];
rz(-2.2089829) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7639097) q[0];
sx q[0];
rz(-1.8881067) q[0];
sx q[0];
rz(-3.0850287) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.89838018) q[2];
sx q[2];
rz(-2.8542238) q[2];
sx q[2];
rz(-2.3786366) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8285117) q[1];
sx q[1];
rz(-2.5701437) q[1];
sx q[1];
rz(-2.9213195) q[1];
x q[2];
rz(-1.014939) q[3];
sx q[3];
rz(-1.1851289) q[3];
sx q[3];
rz(2.8876497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9937667) q[2];
sx q[2];
rz(-2.0596762) q[2];
sx q[2];
rz(-1.0926584) q[2];
rz(2.5993733) q[3];
sx q[3];
rz(-2.0565624) q[3];
sx q[3];
rz(-2.1742163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3595235) q[0];
sx q[0];
rz(-0.096465915) q[0];
sx q[0];
rz(0.50022593) q[0];
rz(-0.80530986) q[1];
sx q[1];
rz(-1.9814682) q[1];
sx q[1];
rz(-1.4979699) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.995979) q[0];
sx q[0];
rz(-2.0659475) q[0];
sx q[0];
rz(-0.33546319) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3802714) q[2];
sx q[2];
rz(-1.5152144) q[2];
sx q[2];
rz(-1.7855102) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.30724635) q[1];
sx q[1];
rz(-1.8028959) q[1];
sx q[1];
rz(-2.8744065) q[1];
x q[2];
rz(0.46521503) q[3];
sx q[3];
rz(-1.9808931) q[3];
sx q[3];
rz(1.0614392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3952289) q[2];
sx q[2];
rz(-2.5791898) q[2];
sx q[2];
rz(-0.70181075) q[2];
rz(-2.3102405) q[3];
sx q[3];
rz(-2.1777007) q[3];
sx q[3];
rz(2.519616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9005301) q[0];
sx q[0];
rz(-2.5456972) q[0];
sx q[0];
rz(0.81533122) q[0];
rz(-1.5218081) q[1];
sx q[1];
rz(-0.83414572) q[1];
sx q[1];
rz(-1.048208) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5883023) q[0];
sx q[0];
rz(-2.7748845) q[0];
sx q[0];
rz(-2.328863) q[0];
rz(-pi) q[1];
rz(0.016096073) q[2];
sx q[2];
rz(-0.65158366) q[2];
sx q[2];
rz(-2.1799257) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.198846) q[1];
sx q[1];
rz(-1.8837351) q[1];
sx q[1];
rz(1.320977) q[1];
x q[2];
rz(-2.7311677) q[3];
sx q[3];
rz(-2.1459208) q[3];
sx q[3];
rz(-2.9068973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.52577019) q[2];
sx q[2];
rz(-0.56695357) q[2];
sx q[2];
rz(-1.099951) q[2];
rz(-2.3163017) q[3];
sx q[3];
rz(-2.0992978) q[3];
sx q[3];
rz(-0.88551372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0734171) q[0];
sx q[0];
rz(-2.5475579) q[0];
sx q[0];
rz(0.90240479) q[0];
rz(-2.1249318) q[1];
sx q[1];
rz(-1.0598176) q[1];
sx q[1];
rz(-0.12983233) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79486217) q[0];
sx q[0];
rz(-2.4299893) q[0];
sx q[0];
rz(-0.58563389) q[0];
rz(-2.1075222) q[2];
sx q[2];
rz(-0.64646361) q[2];
sx q[2];
rz(-1.5922286) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8772014) q[1];
sx q[1];
rz(-1.0069205) q[1];
sx q[1];
rz(-1.1370204) q[1];
rz(-pi) q[2];
rz(0.31452175) q[3];
sx q[3];
rz(-0.57146996) q[3];
sx q[3];
rz(-0.86626569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8292024) q[2];
sx q[2];
rz(-2.1924993) q[2];
sx q[2];
rz(-2.9373346) q[2];
rz(1.9355109) q[3];
sx q[3];
rz(-1.6198502) q[3];
sx q[3];
rz(-0.23541418) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7234574) q[0];
sx q[0];
rz(-1.8122939) q[0];
sx q[0];
rz(1.6947421) q[0];
rz(-1.8824668) q[1];
sx q[1];
rz(-2.1513758) q[1];
sx q[1];
rz(-0.68626219) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9335564) q[0];
sx q[0];
rz(-2.4718923) q[0];
sx q[0];
rz(2.3963388) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3115736) q[2];
sx q[2];
rz(-2.0247211) q[2];
sx q[2];
rz(1.1656851) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2974907) q[1];
sx q[1];
rz(-0.19980783) q[1];
sx q[1];
rz(1.1872477) q[1];
x q[2];
rz(1.8920184) q[3];
sx q[3];
rz(-1.9739082) q[3];
sx q[3];
rz(1.5954799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.69616047) q[2];
sx q[2];
rz(-1.7636718) q[2];
sx q[2];
rz(-3.1398204) q[2];
rz(2.5799675) q[3];
sx q[3];
rz(-0.91149819) q[3];
sx q[3];
rz(-1.6368438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5381662) q[0];
sx q[0];
rz(-0.68646938) q[0];
sx q[0];
rz(1.4461393) q[0];
rz(-0.7810477) q[1];
sx q[1];
rz(-1.8361517) q[1];
sx q[1];
rz(-1.6400281) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7088889) q[0];
sx q[0];
rz(-0.49312691) q[0];
sx q[0];
rz(-1.6119484) q[0];
rz(-pi) q[1];
x q[1];
rz(3.115032) q[2];
sx q[2];
rz(-1.6325258) q[2];
sx q[2];
rz(0.82690566) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1135243) q[1];
sx q[1];
rz(-1.6599732) q[1];
sx q[1];
rz(2.813617) q[1];
x q[2];
rz(-2.5128965) q[3];
sx q[3];
rz(-2.0572212) q[3];
sx q[3];
rz(1.4592255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.35187307) q[2];
sx q[2];
rz(-1.3616273) q[2];
sx q[2];
rz(1.3191351) q[2];
rz(-1.2119279) q[3];
sx q[3];
rz(-1.2865678) q[3];
sx q[3];
rz(-0.31931988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8050352) q[0];
sx q[0];
rz(-2.5890077) q[0];
sx q[0];
rz(1.2040899) q[0];
rz(-2.7583292) q[1];
sx q[1];
rz(-2.6158694) q[1];
sx q[1];
rz(2.7899172) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29031819) q[0];
sx q[0];
rz(-1.1997249) q[0];
sx q[0];
rz(-2.0593658) q[0];
rz(-pi) q[1];
rz(-2.7792764) q[2];
sx q[2];
rz(-3*pi/13) q[2];
sx q[2];
rz(2.97646) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0723567) q[1];
sx q[1];
rz(-1.5178536) q[1];
sx q[1];
rz(-0.80979053) q[1];
x q[2];
rz(-1.6230691) q[3];
sx q[3];
rz(-1.4613323) q[3];
sx q[3];
rz(-2.0930406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7982771) q[2];
sx q[2];
rz(-1.1077935) q[2];
sx q[2];
rz(-1.8593672) q[2];
rz(-1.4964237) q[3];
sx q[3];
rz(-1.5346425) q[3];
sx q[3];
rz(-2.0857247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6431817) q[0];
sx q[0];
rz(-1.8739941) q[0];
sx q[0];
rz(2.9472651) q[0];
rz(2.1037897) q[1];
sx q[1];
rz(-0.56832814) q[1];
sx q[1];
rz(1.0338354) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72458306) q[0];
sx q[0];
rz(-1.4405182) q[0];
sx q[0];
rz(-0.91086046) q[0];
rz(-pi) q[1];
x q[1];
rz(0.98722234) q[2];
sx q[2];
rz(-2.3505031) q[2];
sx q[2];
rz(0.9466048) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2227877) q[1];
sx q[1];
rz(-0.58090392) q[1];
sx q[1];
rz(1.7564299) q[1];
rz(1.3528353) q[3];
sx q[3];
rz(-2.2721604) q[3];
sx q[3];
rz(-2.4234114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0620492) q[2];
sx q[2];
rz(-2.1958308) q[2];
sx q[2];
rz(2.5058084) q[2];
rz(0.27030269) q[3];
sx q[3];
rz(-2.342194) q[3];
sx q[3];
rz(1.5283782) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6939659) q[0];
sx q[0];
rz(-1.8287369) q[0];
sx q[0];
rz(1.0736314) q[0];
rz(-1.7059965) q[1];
sx q[1];
rz(-1.5789079) q[1];
sx q[1];
rz(0.78067738) q[1];
rz(0.96799093) q[2];
sx q[2];
rz(-1.5445166) q[2];
sx q[2];
rz(-1.9894285) q[2];
rz(-0.070449645) q[3];
sx q[3];
rz(-1.2415213) q[3];
sx q[3];
rz(0.51125676) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
