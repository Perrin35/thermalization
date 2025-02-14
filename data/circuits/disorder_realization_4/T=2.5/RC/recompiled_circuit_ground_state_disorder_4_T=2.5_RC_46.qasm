OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5623915) q[0];
sx q[0];
rz(-2.3810823) q[0];
sx q[0];
rz(2.1265246) q[0];
rz(-1.1633582) q[1];
sx q[1];
rz(3.562869) q[1];
sx q[1];
rz(8.1864551) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2531812) q[0];
sx q[0];
rz(-2.0813353) q[0];
sx q[0];
rz(2.2967413) q[0];
rz(-pi) q[1];
rz(2.2767645) q[2];
sx q[2];
rz(-0.53542407) q[2];
sx q[2];
rz(-0.12685093) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1972292) q[1];
sx q[1];
rz(-0.60603332) q[1];
sx q[1];
rz(0.89971772) q[1];
rz(-pi) q[2];
x q[2];
rz(0.56391509) q[3];
sx q[3];
rz(-1.7132572) q[3];
sx q[3];
rz(0.13232732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.51332659) q[2];
sx q[2];
rz(-1.790975) q[2];
sx q[2];
rz(-0.70293054) q[2];
rz(-0.85567307) q[3];
sx q[3];
rz(-0.84508768) q[3];
sx q[3];
rz(-1.8446911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6913476) q[0];
sx q[0];
rz(-2.0424728) q[0];
sx q[0];
rz(-2.3968089) q[0];
rz(-1.567747) q[1];
sx q[1];
rz(-0.66190043) q[1];
sx q[1];
rz(-2.1994798) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73977913) q[0];
sx q[0];
rz(-1.7727858) q[0];
sx q[0];
rz(-0.76633472) q[0];
x q[1];
rz(0.079374119) q[2];
sx q[2];
rz(-2.3217776) q[2];
sx q[2];
rz(0.83460966) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2971645) q[1];
sx q[1];
rz(-1.7065863) q[1];
sx q[1];
rz(-2.335615) q[1];
x q[2];
rz(-1.510624) q[3];
sx q[3];
rz(-2.1561557) q[3];
sx q[3];
rz(-0.28096052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.88910237) q[2];
sx q[2];
rz(-2.1847051) q[2];
sx q[2];
rz(-2.0460184) q[2];
rz(-1.6628294) q[3];
sx q[3];
rz(-0.79083276) q[3];
sx q[3];
rz(1.3862632) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0304994) q[0];
sx q[0];
rz(-1.4249304) q[0];
sx q[0];
rz(1.4053364) q[0];
rz(1.3966712) q[1];
sx q[1];
rz(-1.3537355) q[1];
sx q[1];
rz(-0.90000802) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42381313) q[0];
sx q[0];
rz(-0.85067816) q[0];
sx q[0];
rz(2.6507342) q[0];
rz(0.58061231) q[2];
sx q[2];
rz(-2.0502649) q[2];
sx q[2];
rz(-1.978385) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.69086466) q[1];
sx q[1];
rz(-2.303586) q[1];
sx q[1];
rz(1.3370675) q[1];
rz(-pi) q[2];
rz(2.4803512) q[3];
sx q[3];
rz(-2.0783278) q[3];
sx q[3];
rz(0.37933168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6806543) q[2];
sx q[2];
rz(-2.9260981) q[2];
sx q[2];
rz(0.94949618) q[2];
rz(0.62120581) q[3];
sx q[3];
rz(-2.216279) q[3];
sx q[3];
rz(1.9882103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7200274) q[0];
sx q[0];
rz(-1.8864487) q[0];
sx q[0];
rz(-3.0928639) q[0];
rz(0.85982927) q[1];
sx q[1];
rz(-2.8099334) q[1];
sx q[1];
rz(2.8299832) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45845073) q[0];
sx q[0];
rz(-1.299017) q[0];
sx q[0];
rz(-2.0608451) q[0];
x q[1];
rz(2.0987058) q[2];
sx q[2];
rz(-1.2706869) q[2];
sx q[2];
rz(-1.2712511) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.58285634) q[1];
sx q[1];
rz(-1.8355825) q[1];
sx q[1];
rz(-2.2469421) q[1];
rz(0.29693691) q[3];
sx q[3];
rz(-1.9949759) q[3];
sx q[3];
rz(2.3564828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9638046) q[2];
sx q[2];
rz(-0.21169855) q[2];
sx q[2];
rz(-1.6507899) q[2];
rz(-0.35495159) q[3];
sx q[3];
rz(-1.7345411) q[3];
sx q[3];
rz(1.2995592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-2.0519003) q[0];
sx q[0];
rz(-2.7841452) q[0];
sx q[0];
rz(1.8748913) q[0];
rz(-0.1162687) q[1];
sx q[1];
rz(-2.1513042) q[1];
sx q[1];
rz(1.5142534) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26402347) q[0];
sx q[0];
rz(-2.2248587) q[0];
sx q[0];
rz(-2.7436849) q[0];
x q[1];
rz(1.9562938) q[2];
sx q[2];
rz(-0.23699871) q[2];
sx q[2];
rz(-2.1182107) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0872495) q[1];
sx q[1];
rz(-1.5197481) q[1];
sx q[1];
rz(1.0264978) q[1];
rz(-0.037260696) q[3];
sx q[3];
rz(-2.2989591) q[3];
sx q[3];
rz(-1.3986349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9194455) q[2];
sx q[2];
rz(-2.6667892) q[2];
sx q[2];
rz(-1.1754645) q[2];
rz(0.90421024) q[3];
sx q[3];
rz(-1.2491106) q[3];
sx q[3];
rz(0.97833943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0745875) q[0];
sx q[0];
rz(-0.33127221) q[0];
sx q[0];
rz(-0.40147716) q[0];
rz(1.1270771) q[1];
sx q[1];
rz(-1.8311484) q[1];
sx q[1];
rz(0.81261596) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3881587) q[0];
sx q[0];
rz(-1.5749314) q[0];
sx q[0];
rz(-2.2049516) q[0];
rz(-1.8046384) q[2];
sx q[2];
rz(-1.8385781) q[2];
sx q[2];
rz(3.0674231) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9795831) q[1];
sx q[1];
rz(-0.97921959) q[1];
sx q[1];
rz(-1.3237557) q[1];
x q[2];
rz(-0.74515588) q[3];
sx q[3];
rz(-1.1155012) q[3];
sx q[3];
rz(0.32768341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9408985) q[2];
sx q[2];
rz(-2.5914067) q[2];
sx q[2];
rz(1.5563439) q[2];
rz(-2.9511792) q[3];
sx q[3];
rz(-0.75473458) q[3];
sx q[3];
rz(1.93369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3170526) q[0];
sx q[0];
rz(-0.38924488) q[0];
sx q[0];
rz(-0.12271605) q[0];
rz(1.9484733) q[1];
sx q[1];
rz(-0.90843186) q[1];
sx q[1];
rz(-0.76748031) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7372663) q[0];
sx q[0];
rz(-2.4929059) q[0];
sx q[0];
rz(1.4859597) q[0];
rz(-pi) q[1];
rz(1.6558455) q[2];
sx q[2];
rz(-1.1901996) q[2];
sx q[2];
rz(-2.2589661) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.40015107) q[1];
sx q[1];
rz(-1.5269244) q[1];
sx q[1];
rz(-0.9167826) q[1];
rz(3.1311099) q[3];
sx q[3];
rz(-0.40626486) q[3];
sx q[3];
rz(-0.83028136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3976589) q[2];
sx q[2];
rz(-3.1199516) q[2];
sx q[2];
rz(-1.5519315) q[2];
rz(-1.4895561) q[3];
sx q[3];
rz(-1.4068312) q[3];
sx q[3];
rz(-2.0261197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3234696) q[0];
sx q[0];
rz(-0.36190811) q[0];
sx q[0];
rz(2.7662011) q[0];
rz(-0.88343945) q[1];
sx q[1];
rz(-1.6049623) q[1];
sx q[1];
rz(2.4129131) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.323384) q[0];
sx q[0];
rz(-2.2012156) q[0];
sx q[0];
rz(2.4805043) q[0];
x q[1];
rz(-1.402114) q[2];
sx q[2];
rz(-1.8034435) q[2];
sx q[2];
rz(0.53984914) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9364215) q[1];
sx q[1];
rz(-1.9269619) q[1];
sx q[1];
rz(-2.2486399) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.65547319) q[3];
sx q[3];
rz(-2.270916) q[3];
sx q[3];
rz(1.164013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4484619) q[2];
sx q[2];
rz(-1.5631661) q[2];
sx q[2];
rz(2.3942088) q[2];
rz(1.4061617) q[3];
sx q[3];
rz(-1.2179255) q[3];
sx q[3];
rz(1.5860484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9091699) q[0];
sx q[0];
rz(-2.2012043) q[0];
sx q[0];
rz(0.7255834) q[0];
rz(-1.3369417) q[1];
sx q[1];
rz(-1.2533816) q[1];
sx q[1];
rz(-2.5673089) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5919898) q[0];
sx q[0];
rz(-1.4583602) q[0];
sx q[0];
rz(1.7312584) q[0];
x q[1];
rz(-2.7321283) q[2];
sx q[2];
rz(-1.0396233) q[2];
sx q[2];
rz(2.766618) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0813539) q[1];
sx q[1];
rz(-1.4574086) q[1];
sx q[1];
rz(0.84073034) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1508997) q[3];
sx q[3];
rz(-2.266006) q[3];
sx q[3];
rz(0.08838544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1672704) q[2];
sx q[2];
rz(-1.7631301) q[2];
sx q[2];
rz(-0.66217011) q[2];
rz(-2.3769489) q[3];
sx q[3];
rz(-1.5352826) q[3];
sx q[3];
rz(-1.3242599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26226703) q[0];
sx q[0];
rz(-2.5536394) q[0];
sx q[0];
rz(1.3978488) q[0];
rz(1.0700048) q[1];
sx q[1];
rz(-1.1281697) q[1];
sx q[1];
rz(-1.8096583) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2764212) q[0];
sx q[0];
rz(-1.7214473) q[0];
sx q[0];
rz(-0.030929203) q[0];
x q[1];
rz(1.9851763) q[2];
sx q[2];
rz(-1.7935026) q[2];
sx q[2];
rz(3.0321995) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8269013) q[1];
sx q[1];
rz(-1.3413635) q[1];
sx q[1];
rz(1.3342085) q[1];
rz(-pi) q[2];
rz(2.2666107) q[3];
sx q[3];
rz(-2.0379645) q[3];
sx q[3];
rz(2.9206729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8269044) q[2];
sx q[2];
rz(-1.1004227) q[2];
sx q[2];
rz(-2.9617214) q[2];
rz(1.7449069) q[3];
sx q[3];
rz(-1.7161918) q[3];
sx q[3];
rz(-2.9250308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3708645) q[0];
sx q[0];
rz(-2.6612119) q[0];
sx q[0];
rz(-2.2330855) q[0];
rz(0.52275672) q[1];
sx q[1];
rz(-2.0261384) q[1];
sx q[1];
rz(-1.1631858) q[1];
rz(1.9942453) q[2];
sx q[2];
rz(-2.1622787) q[2];
sx q[2];
rz(2.5522243) q[2];
rz(-1.9369851) q[3];
sx q[3];
rz(-1.1259176) q[3];
sx q[3];
rz(2.068145) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
