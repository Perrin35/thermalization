OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.22566158) q[0];
sx q[0];
rz(-2.2731279) q[0];
sx q[0];
rz(0.19302364) q[0];
rz(1.141619) q[1];
sx q[1];
rz(-0.42998278) q[1];
sx q[1];
rz(2.4584682) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1238414) q[0];
sx q[0];
rz(-3.0525644) q[0];
sx q[0];
rz(2.428399) q[0];
x q[1];
rz(2.9002951) q[2];
sx q[2];
rz(-1.2295051) q[2];
sx q[2];
rz(-1.2933921) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7072304) q[1];
sx q[1];
rz(-1.1303567) q[1];
sx q[1];
rz(-0.67727725) q[1];
x q[2];
rz(-1.0115511) q[3];
sx q[3];
rz(-2.4992001) q[3];
sx q[3];
rz(-1.1005644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1258939) q[2];
sx q[2];
rz(-1.761972) q[2];
sx q[2];
rz(1.0985628) q[2];
rz(2.0627608) q[3];
sx q[3];
rz(-0.94509411) q[3];
sx q[3];
rz(-1.1014972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9803479) q[0];
sx q[0];
rz(-2.8256567) q[0];
sx q[0];
rz(-2.9336477) q[0];
rz(2.5646599) q[1];
sx q[1];
rz(-0.88795841) q[1];
sx q[1];
rz(1.4651441) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0365021) q[0];
sx q[0];
rz(-2.4173792) q[0];
sx q[0];
rz(0.0037395517) q[0];
rz(-pi) q[1];
rz(-1.0826153) q[2];
sx q[2];
rz(-1.7980051) q[2];
sx q[2];
rz(2.8607334) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1190471) q[1];
sx q[1];
rz(-2.3733768) q[1];
sx q[1];
rz(-0.70336282) q[1];
rz(-0.37104718) q[3];
sx q[3];
rz(-0.78073946) q[3];
sx q[3];
rz(0.77663976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0186105) q[2];
sx q[2];
rz(-2.1554422) q[2];
sx q[2];
rz(-3.0318276) q[2];
rz(0.62260735) q[3];
sx q[3];
rz(-2.770335) q[3];
sx q[3];
rz(1.4573147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96317545) q[0];
sx q[0];
rz(-2.2816179) q[0];
sx q[0];
rz(-0.51613581) q[0];
rz(0.57488817) q[1];
sx q[1];
rz(-2.2153885) q[1];
sx q[1];
rz(0.80054545) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8974509) q[0];
sx q[0];
rz(-1.4276917) q[0];
sx q[0];
rz(1.1600526) q[0];
rz(0.66652253) q[2];
sx q[2];
rz(-0.98072532) q[2];
sx q[2];
rz(-2.9556264) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2785981) q[1];
sx q[1];
rz(-1.9837712) q[1];
sx q[1];
rz(-1.5763271) q[1];
rz(-1.21739) q[3];
sx q[3];
rz(-0.73261515) q[3];
sx q[3];
rz(0.64980799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4642554) q[2];
sx q[2];
rz(-0.3192454) q[2];
sx q[2];
rz(-1.7819972) q[2];
rz(0.96757403) q[3];
sx q[3];
rz(-1.273497) q[3];
sx q[3];
rz(1.7165855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1490705) q[0];
sx q[0];
rz(-1.2620121) q[0];
sx q[0];
rz(0.46491369) q[0];
rz(0.34856302) q[1];
sx q[1];
rz(-2.8788853) q[1];
sx q[1];
rz(-1.0850614) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5497919) q[0];
sx q[0];
rz(-1.8013957) q[0];
sx q[0];
rz(0.53996284) q[0];
x q[1];
rz(-0.24809804) q[2];
sx q[2];
rz(-1.9271701) q[2];
sx q[2];
rz(-0.46596371) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.76772056) q[1];
sx q[1];
rz(-1.671119) q[1];
sx q[1];
rz(-2.9994681) q[1];
x q[2];
rz(2.5303909) q[3];
sx q[3];
rz(-1.9107606) q[3];
sx q[3];
rz(1.9577648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.00502914) q[2];
sx q[2];
rz(-2.0932784) q[2];
sx q[2];
rz(0.68112779) q[2];
rz(2.629771) q[3];
sx q[3];
rz(-2.8183283) q[3];
sx q[3];
rz(-0.26369035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(2.8624449) q[0];
sx q[0];
rz(-0.537916) q[0];
sx q[0];
rz(-1.7329247) q[0];
rz(-2.7092343) q[1];
sx q[1];
rz(-2.2996348) q[1];
sx q[1];
rz(-0.98168215) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0911134) q[0];
sx q[0];
rz(-2.1149153) q[0];
sx q[0];
rz(1.0247466) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8129187) q[2];
sx q[2];
rz(-1.4171346) q[2];
sx q[2];
rz(-0.15649934) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.86451929) q[1];
sx q[1];
rz(-2.652087) q[1];
sx q[1];
rz(-2.8237052) q[1];
rz(-pi) q[2];
rz(-1.4622299) q[3];
sx q[3];
rz(-0.6001937) q[3];
sx q[3];
rz(-1.997228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1061873) q[2];
sx q[2];
rz(-1.4865439) q[2];
sx q[2];
rz(-0.48941082) q[2];
rz(1.0148467) q[3];
sx q[3];
rz(-2.615052) q[3];
sx q[3];
rz(1.813252) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4846102) q[0];
sx q[0];
rz(-2.9841612) q[0];
sx q[0];
rz(0.37242517) q[0];
rz(-1.8107481) q[1];
sx q[1];
rz(-1.0354038) q[1];
sx q[1];
rz(-0.16528027) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80711354) q[0];
sx q[0];
rz(-1.4757336) q[0];
sx q[0];
rz(1.5414184) q[0];
rz(-pi) q[1];
rz(1.500962) q[2];
sx q[2];
rz(-2.2572821) q[2];
sx q[2];
rz(1.3697461) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2416934) q[1];
sx q[1];
rz(-2.4447933) q[1];
sx q[1];
rz(-2.4844869) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23541707) q[3];
sx q[3];
rz(-1.132292) q[3];
sx q[3];
rz(-0.32270839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.91339397) q[2];
sx q[2];
rz(-2.1134351) q[2];
sx q[2];
rz(1.6983263) q[2];
rz(0.36744395) q[3];
sx q[3];
rz(-1.2747217) q[3];
sx q[3];
rz(2.929556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3570324) q[0];
sx q[0];
rz(-1.0914047) q[0];
sx q[0];
rz(2.7923287) q[0];
rz(-2.3941984) q[1];
sx q[1];
rz(-2.8458197) q[1];
sx q[1];
rz(2.4051037) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38024494) q[0];
sx q[0];
rz(-1.5537098) q[0];
sx q[0];
rz(-1.5218309) q[0];
rz(-pi) q[1];
rz(1.2320802) q[2];
sx q[2];
rz(-1.806353) q[2];
sx q[2];
rz(-0.92600694) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4941102) q[1];
sx q[1];
rz(-2.3586015) q[1];
sx q[1];
rz(2.8857735) q[1];
rz(-pi) q[2];
rz(-2.6408259) q[3];
sx q[3];
rz(-1.9541249) q[3];
sx q[3];
rz(2.4367743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.13654576) q[2];
sx q[2];
rz(-1.9133291) q[2];
sx q[2];
rz(0.53156701) q[2];
rz(2.452204) q[3];
sx q[3];
rz(-1.4644943) q[3];
sx q[3];
rz(1.846107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0625967) q[0];
sx q[0];
rz(-1.9364708) q[0];
sx q[0];
rz(-1.8435562) q[0];
rz(0.80728665) q[1];
sx q[1];
rz(-1.1785945) q[1];
sx q[1];
rz(2.2198026) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4661134) q[0];
sx q[0];
rz(-2.4371494) q[0];
sx q[0];
rz(0.2091614) q[0];
rz(-pi) q[1];
rz(-1.2577031) q[2];
sx q[2];
rz(-0.76258341) q[2];
sx q[2];
rz(-1.2459754) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3075879) q[1];
sx q[1];
rz(-1.7760135) q[1];
sx q[1];
rz(0.25968857) q[1];
rz(-1.3619625) q[3];
sx q[3];
rz(-0.71912557) q[3];
sx q[3];
rz(1.5428839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.90199295) q[2];
sx q[2];
rz(-2.5805876) q[2];
sx q[2];
rz(1.1716589) q[2];
rz(-1.7840067) q[3];
sx q[3];
rz(-1.4566908) q[3];
sx q[3];
rz(1.6931036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25306025) q[0];
sx q[0];
rz(-1.9655515) q[0];
sx q[0];
rz(1.6660447) q[0];
rz(1.5400003) q[1];
sx q[1];
rz(-1.4614636) q[1];
sx q[1];
rz(0.67970651) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98201671) q[0];
sx q[0];
rz(-1.3681108) q[0];
sx q[0];
rz(2.6152339) q[0];
x q[1];
rz(0.87551261) q[2];
sx q[2];
rz(-1.7098223) q[2];
sx q[2];
rz(1.306844) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.745984) q[1];
sx q[1];
rz(-1.9095699) q[1];
sx q[1];
rz(0.70396522) q[1];
rz(2.4623975) q[3];
sx q[3];
rz(-1.067357) q[3];
sx q[3];
rz(-2.5356843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1264964) q[2];
sx q[2];
rz(-1.6588147) q[2];
sx q[2];
rz(0.91040197) q[2];
rz(0.67534584) q[3];
sx q[3];
rz(-0.93674913) q[3];
sx q[3];
rz(2.2909686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0891721) q[0];
sx q[0];
rz(-1.9071254) q[0];
sx q[0];
rz(2.4269379) q[0];
rz(-0.71406281) q[1];
sx q[1];
rz(-0.95497447) q[1];
sx q[1];
rz(-1.9649327) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.048621) q[0];
sx q[0];
rz(-1.4074568) q[0];
sx q[0];
rz(3.0924348) q[0];
rz(-pi) q[1];
rz(2.2628225) q[2];
sx q[2];
rz(-2.3792017) q[2];
sx q[2];
rz(-2.9825485) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9709819) q[1];
sx q[1];
rz(-0.5628399) q[1];
sx q[1];
rz(-1.7112205) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6305466) q[3];
sx q[3];
rz(-1.647445) q[3];
sx q[3];
rz(-2.6443036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8979793) q[2];
sx q[2];
rz(-1.4695797) q[2];
sx q[2];
rz(2.0142377) q[2];
rz(0.35774287) q[3];
sx q[3];
rz(-2.2704411) q[3];
sx q[3];
rz(-2.0991142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7174299) q[0];
sx q[0];
rz(-1.9308199) q[0];
sx q[0];
rz(0.45817026) q[0];
rz(2.7453616) q[1];
sx q[1];
rz(-3.1165262) q[1];
sx q[1];
rz(-2.810626) q[1];
rz(1.8489807) q[2];
sx q[2];
rz(-0.85360151) q[2];
sx q[2];
rz(-3.1384946) q[2];
rz(-0.35076326) q[3];
sx q[3];
rz(-1.5784932) q[3];
sx q[3];
rz(-2.2898883) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
