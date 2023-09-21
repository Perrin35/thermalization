OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.84500399) q[0];
sx q[0];
rz(-0.72405976) q[0];
sx q[0];
rz(-1.6568503) q[0];
rz(1.2031263) q[1];
sx q[1];
rz(-0.523518) q[1];
sx q[1];
rz(2.2533921) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6523525) q[0];
sx q[0];
rz(-1.3999192) q[0];
sx q[0];
rz(-1.5255552) q[0];
rz(-pi) q[1];
rz(2.9759679) q[2];
sx q[2];
rz(-1.0399482) q[2];
sx q[2];
rz(-0.80425516) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8468195) q[1];
sx q[1];
rz(-1.5702489) q[1];
sx q[1];
rz(-0.16695395) q[1];
rz(-pi) q[2];
rz(-1.5376484) q[3];
sx q[3];
rz(-0.87246934) q[3];
sx q[3];
rz(-2.2693199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.28329864) q[2];
sx q[2];
rz(-2.7259939) q[2];
sx q[2];
rz(-1.1179914) q[2];
rz(0.14532146) q[3];
sx q[3];
rz(-1.5829007) q[3];
sx q[3];
rz(-0.045923559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5730826) q[0];
sx q[0];
rz(-1.8772323) q[0];
sx q[0];
rz(-0.65482393) q[0];
rz(1.9251992) q[1];
sx q[1];
rz(-1.9768068) q[1];
sx q[1];
rz(3.0156946) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.038371041) q[0];
sx q[0];
rz(-1.201259) q[0];
sx q[0];
rz(-1.0612556) q[0];
x q[1];
rz(0.6217896) q[2];
sx q[2];
rz(-0.99025531) q[2];
sx q[2];
rz(0.34678005) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1344578) q[1];
sx q[1];
rz(-1.6975228) q[1];
sx q[1];
rz(0.030220672) q[1];
rz(-pi) q[2];
rz(1.5586073) q[3];
sx q[3];
rz(-2.4180531) q[3];
sx q[3];
rz(-2.7140868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.048432365) q[2];
sx q[2];
rz(-1.9530714) q[2];
sx q[2];
rz(0.38802567) q[2];
rz(1.4240501) q[3];
sx q[3];
rz(-2.5035796) q[3];
sx q[3];
rz(2.5542636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5557142) q[0];
sx q[0];
rz(-2.3590187) q[0];
sx q[0];
rz(-0.079332381) q[0];
rz(-0.084005984) q[1];
sx q[1];
rz(-2.3386798) q[1];
sx q[1];
rz(1.1598587) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0214329) q[0];
sx q[0];
rz(-2.624357) q[0];
sx q[0];
rz(-1.73818) q[0];
rz(-2.2247505) q[2];
sx q[2];
rz(-0.7407623) q[2];
sx q[2];
rz(-0.11292085) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.11324727) q[1];
sx q[1];
rz(-2.5302437) q[1];
sx q[1];
rz(-2.4436185) q[1];
rz(-pi) q[2];
rz(1.4781811) q[3];
sx q[3];
rz(-2.5941656) q[3];
sx q[3];
rz(0.80143354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7362061) q[2];
sx q[2];
rz(-1.3233041) q[2];
sx q[2];
rz(3.0333056) q[2];
rz(0.64374271) q[3];
sx q[3];
rz(-1.0780004) q[3];
sx q[3];
rz(2.5260177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9765587) q[0];
sx q[0];
rz(-1.3950011) q[0];
sx q[0];
rz(1.4820341) q[0];
rz(-0.7011134) q[1];
sx q[1];
rz(-1.8891524) q[1];
sx q[1];
rz(2.8569417) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92589256) q[0];
sx q[0];
rz(-1.6452351) q[0];
sx q[0];
rz(2.0490993) q[0];
rz(-2.9361211) q[2];
sx q[2];
rz(-1.2887508) q[2];
sx q[2];
rz(-3.0021283) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.33672562) q[1];
sx q[1];
rz(-1.4443828) q[1];
sx q[1];
rz(-3.0625507) q[1];
rz(-pi) q[2];
rz(0.986226) q[3];
sx q[3];
rz(-0.75891906) q[3];
sx q[3];
rz(0.81134568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9099137) q[2];
sx q[2];
rz(-0.96857962) q[2];
sx q[2];
rz(0.56751928) q[2];
rz(0.41401687) q[3];
sx q[3];
rz(-1.1375789) q[3];
sx q[3];
rz(-2.0295985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48859566) q[0];
sx q[0];
rz(-2.1814006) q[0];
sx q[0];
rz(1.8150785) q[0];
rz(1.1524221) q[1];
sx q[1];
rz(-1.3782586) q[1];
sx q[1];
rz(-0.93793905) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4834578) q[0];
sx q[0];
rz(-1.5926653) q[0];
sx q[0];
rz(1.6629501) q[0];
rz(-1.0114952) q[2];
sx q[2];
rz(-2.0213631) q[2];
sx q[2];
rz(0.089103854) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8541969) q[1];
sx q[1];
rz(-1.5970267) q[1];
sx q[1];
rz(-1.5363974) q[1];
rz(0.33850833) q[3];
sx q[3];
rz(-0.37422985) q[3];
sx q[3];
rz(-2.3238473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1753297) q[2];
sx q[2];
rz(-0.18925174) q[2];
sx q[2];
rz(-1.0413292) q[2];
rz(1.8390309) q[3];
sx q[3];
rz(-2.0077191) q[3];
sx q[3];
rz(1.8235824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.702521) q[0];
sx q[0];
rz(-1.9938001) q[0];
sx q[0];
rz(2.5842216) q[0];
rz(-2.5769261) q[1];
sx q[1];
rz(-0.70960418) q[1];
sx q[1];
rz(-0.55647892) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59069955) q[0];
sx q[0];
rz(-1.3644344) q[0];
sx q[0];
rz(-0.38462374) q[0];
rz(-1.7700023) q[2];
sx q[2];
rz(-0.7253941) q[2];
sx q[2];
rz(-0.2142011) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.52292697) q[1];
sx q[1];
rz(-1.7523125) q[1];
sx q[1];
rz(1.8261441) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5197093) q[3];
sx q[3];
rz(-0.79569492) q[3];
sx q[3];
rz(-1.6076128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6094728) q[2];
sx q[2];
rz(-0.76081053) q[2];
sx q[2];
rz(-1.2247941) q[2];
rz(1.5103643) q[3];
sx q[3];
rz(-1.8211726) q[3];
sx q[3];
rz(-1.640655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0790134) q[0];
sx q[0];
rz(-1.8824848) q[0];
sx q[0];
rz(-3.0986837) q[0];
rz(0.91730109) q[1];
sx q[1];
rz(-0.62364548) q[1];
sx q[1];
rz(-2.6409805) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6635054) q[0];
sx q[0];
rz(-0.18112077) q[0];
sx q[0];
rz(-1.3785133) q[0];
x q[1];
rz(1.0285573) q[2];
sx q[2];
rz(-1.2985897) q[2];
sx q[2];
rz(0.82424639) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7350406) q[1];
sx q[1];
rz(-2.7222689) q[1];
sx q[1];
rz(-0.29962824) q[1];
rz(1.9686437) q[3];
sx q[3];
rz(-0.85507353) q[3];
sx q[3];
rz(1.0637103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6033972) q[2];
sx q[2];
rz(-1.0791225) q[2];
sx q[2];
rz(-0.67561692) q[2];
rz(2.7006941) q[3];
sx q[3];
rz(-1.4392122) q[3];
sx q[3];
rz(2.6202257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.065141) q[0];
sx q[0];
rz(-2.8171709) q[0];
sx q[0];
rz(-2.0741529) q[0];
rz(0.43287977) q[1];
sx q[1];
rz(-1.6378816) q[1];
sx q[1];
rz(1.4656461) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.081799) q[0];
sx q[0];
rz(-2.0311211) q[0];
sx q[0];
rz(-0.1501118) q[0];
rz(-0.52532105) q[2];
sx q[2];
rz(-0.21481951) q[2];
sx q[2];
rz(-2.7810682) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.88361909) q[1];
sx q[1];
rz(-0.61425629) q[1];
sx q[1];
rz(-2.1780464) q[1];
rz(0.74387868) q[3];
sx q[3];
rz(-2.3489967) q[3];
sx q[3];
rz(-2.511123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8026768) q[2];
sx q[2];
rz(-0.85614506) q[2];
sx q[2];
rz(-1.476293) q[2];
rz(0.87351292) q[3];
sx q[3];
rz(-0.87681186) q[3];
sx q[3];
rz(2.6749558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.8840238) q[0];
sx q[0];
rz(-2.5654061) q[0];
sx q[0];
rz(2.4043758) q[0];
rz(3.1228512) q[1];
sx q[1];
rz(-2.811921) q[1];
sx q[1];
rz(-2.2163056) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6773274) q[0];
sx q[0];
rz(-0.87809169) q[0];
sx q[0];
rz(2.4753184) q[0];
rz(-pi) q[1];
rz(1.4383573) q[2];
sx q[2];
rz(-0.91184154) q[2];
sx q[2];
rz(-0.027564136) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.195203) q[1];
sx q[1];
rz(-2.9159947) q[1];
sx q[1];
rz(2.5220847) q[1];
rz(-pi) q[2];
rz(-2.7221189) q[3];
sx q[3];
rz(-0.89722108) q[3];
sx q[3];
rz(-2.2485965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3725738) q[2];
sx q[2];
rz(-1.5344658) q[2];
sx q[2];
rz(2.1378689) q[2];
rz(-3.051493) q[3];
sx q[3];
rz(-3.1153479) q[3];
sx q[3];
rz(-1.1130921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1019679) q[0];
sx q[0];
rz(-1.5682546) q[0];
sx q[0];
rz(1.8819303) q[0];
rz(2.8758077) q[1];
sx q[1];
rz(-0.62756413) q[1];
sx q[1];
rz(2.3840747) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82523358) q[0];
sx q[0];
rz(-2.9037583) q[0];
sx q[0];
rz(2.1912078) q[0];
rz(-1.8025814) q[2];
sx q[2];
rz(-2.5387562) q[2];
sx q[2];
rz(-1.9581118) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6706558) q[1];
sx q[1];
rz(-2.6568012) q[1];
sx q[1];
rz(-2.0623341) q[1];
rz(1.5424472) q[3];
sx q[3];
rz(-2.2806014) q[3];
sx q[3];
rz(0.2809487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1511128) q[2];
sx q[2];
rz(-1.8258784) q[2];
sx q[2];
rz(2.347836) q[2];
rz(-0.63888597) q[3];
sx q[3];
rz(-2.7975438) q[3];
sx q[3];
rz(-1.0206153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4929852) q[0];
sx q[0];
rz(-2.6633371) q[0];
sx q[0];
rz(2.2289842) q[0];
rz(-1.6408625) q[1];
sx q[1];
rz(-0.91703569) q[1];
sx q[1];
rz(-1.348319) q[1];
rz(0.43754229) q[2];
sx q[2];
rz(-1.2699288) q[2];
sx q[2];
rz(-1.5581836) q[2];
rz(0.32140857) q[3];
sx q[3];
rz(-1.0151498) q[3];
sx q[3];
rz(-2.2890454) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];