OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.9159311) q[0];
sx q[0];
rz(-0.8684648) q[0];
sx q[0];
rz(2.948569) q[0];
rz(-1.9999737) q[1];
sx q[1];
rz(3.5715754) q[1];
sx q[1];
rz(6.9663098) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8389991) q[0];
sx q[0];
rz(-1.6380881) q[0];
sx q[0];
rz(1.6291314) q[0];
x q[1];
rz(2.9002951) q[2];
sx q[2];
rz(-1.9120875) q[2];
sx q[2];
rz(1.2933921) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.7636823) q[3];
sx q[3];
rz(-2.1034735) q[3];
sx q[3];
rz(-2.7045254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0156988) q[2];
sx q[2];
rz(-1.3796207) q[2];
sx q[2];
rz(2.0430298) q[2];
rz(-2.0627608) q[3];
sx q[3];
rz(-0.94509411) q[3];
sx q[3];
rz(1.1014972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1612448) q[0];
sx q[0];
rz(-0.315936) q[0];
sx q[0];
rz(-0.20794491) q[0];
rz(-2.5646599) q[1];
sx q[1];
rz(-2.2536342) q[1];
sx q[1];
rz(-1.6764486) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11008308) q[0];
sx q[0];
rz(-2.2950036) q[0];
sx q[0];
rz(-1.5674885) q[0];
rz(-2.8855578) q[2];
sx q[2];
rz(-1.0962152) q[2];
sx q[2];
rz(1.408996) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2966753) q[1];
sx q[1];
rz(-1.0122609) q[1];
sx q[1];
rz(-2.1293473) q[1];
rz(-pi) q[2];
rz(2.39605) q[3];
sx q[3];
rz(-1.8288444) q[3];
sx q[3];
rz(-2.6170956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0186105) q[2];
sx q[2];
rz(-2.1554422) q[2];
sx q[2];
rz(-0.1097651) q[2];
rz(2.5189853) q[3];
sx q[3];
rz(-0.37125769) q[3];
sx q[3];
rz(-1.6842779) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1784172) q[0];
sx q[0];
rz(-2.2816179) q[0];
sx q[0];
rz(2.6254568) q[0];
rz(-0.57488817) q[1];
sx q[1];
rz(-2.2153885) q[1];
sx q[1];
rz(-0.80054545) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0101937) q[0];
sx q[0];
rz(-0.43361615) q[0];
sx q[0];
rz(1.9171159) q[0];
rz(0.66652253) q[2];
sx q[2];
rz(-0.98072532) q[2];
sx q[2];
rz(0.18596622) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2785981) q[1];
sx q[1];
rz(-1.1578214) q[1];
sx q[1];
rz(-1.5763271) q[1];
rz(-2.839746) q[3];
sx q[3];
rz(-0.89248025) q[3];
sx q[3];
rz(2.9523926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.67733726) q[2];
sx q[2];
rz(-2.8223473) q[2];
sx q[2];
rz(1.3595954) q[2];
rz(-0.96757403) q[3];
sx q[3];
rz(-1.8680957) q[3];
sx q[3];
rz(-1.4250071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1490705) q[0];
sx q[0];
rz(-1.2620121) q[0];
sx q[0];
rz(0.46491369) q[0];
rz(2.7930296) q[1];
sx q[1];
rz(-0.26270738) q[1];
sx q[1];
rz(2.0565313) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1151428) q[0];
sx q[0];
rz(-1.0466252) q[0];
sx q[0];
rz(1.8379704) q[0];
x q[1];
rz(2.8934946) q[2];
sx q[2];
rz(-1.2144226) q[2];
sx q[2];
rz(0.46596371) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3241868) q[1];
sx q[1];
rz(-1.7122014) q[1];
sx q[1];
rz(1.6721339) q[1];
x q[2];
rz(0.55235858) q[3];
sx q[3];
rz(-2.4529152) q[3];
sx q[3];
rz(-0.83113447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1365635) q[2];
sx q[2];
rz(-1.0483142) q[2];
sx q[2];
rz(2.4604649) q[2];
rz(0.51182169) q[3];
sx q[3];
rz(-2.8183283) q[3];
sx q[3];
rz(-2.8779023) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8624449) q[0];
sx q[0];
rz(-0.537916) q[0];
sx q[0];
rz(1.7329247) q[0];
rz(-0.43235835) q[1];
sx q[1];
rz(-2.2996348) q[1];
sx q[1];
rz(0.98168215) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0504793) q[0];
sx q[0];
rz(-2.1149153) q[0];
sx q[0];
rz(-2.1168461) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.9972516) q[2];
sx q[2];
rz(-2.855636) q[2];
sx q[2];
rz(-2.2821102) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.423646) q[1];
sx q[1];
rz(-1.7182933) q[1];
sx q[1];
rz(2.6731078) q[1];
x q[2];
rz(1.4622299) q[3];
sx q[3];
rz(-0.6001937) q[3];
sx q[3];
rz(-1.1443646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0354054) q[2];
sx q[2];
rz(-1.6550487) q[2];
sx q[2];
rz(0.48941082) q[2];
rz(2.126746) q[3];
sx q[3];
rz(-2.615052) q[3];
sx q[3];
rz(-1.813252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.6569825) q[0];
sx q[0];
rz(-2.9841612) q[0];
sx q[0];
rz(-0.37242517) q[0];
rz(-1.8107481) q[1];
sx q[1];
rz(-2.1061888) q[1];
sx q[1];
rz(0.16528027) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76089345) q[0];
sx q[0];
rz(-1.5415511) q[0];
sx q[0];
rz(0.0951035) q[0];
rz(-3.0566453) q[2];
sx q[2];
rz(-0.68945486) q[2];
sx q[2];
rz(-1.2598318) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2051516) q[1];
sx q[1];
rz(-1.1679822) q[1];
sx q[1];
rz(-0.58516296) q[1];
rz(2.9061756) q[3];
sx q[3];
rz(-2.0093007) q[3];
sx q[3];
rz(-2.8188843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2281987) q[2];
sx q[2];
rz(-2.1134351) q[2];
sx q[2];
rz(1.4432663) q[2];
rz(-0.36744395) q[3];
sx q[3];
rz(-1.2747217) q[3];
sx q[3];
rz(-2.929556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3570324) q[0];
sx q[0];
rz(-2.050188) q[0];
sx q[0];
rz(-2.7923287) q[0];
rz(-2.3941984) q[1];
sx q[1];
rz(-2.8458197) q[1];
sx q[1];
rz(-0.73648891) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6155646) q[0];
sx q[0];
rz(-0.051858735) q[0];
sx q[0];
rz(-1.9066914) q[0];
rz(-pi) q[1];
rz(0.94524224) q[2];
sx q[2];
rz(-2.7316299) q[2];
sx q[2];
rz(-1.9117102) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0349717) q[1];
sx q[1];
rz(-1.391341) q[1];
sx q[1];
rz(-2.3751395) q[1];
x q[2];
rz(-0.69865366) q[3];
sx q[3];
rz(-0.62056382) q[3];
sx q[3];
rz(0.26649775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.13654576) q[2];
sx q[2];
rz(-1.9133291) q[2];
sx q[2];
rz(-0.53156701) q[2];
rz(-0.68938869) q[3];
sx q[3];
rz(-1.6770984) q[3];
sx q[3];
rz(-1.846107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0625967) q[0];
sx q[0];
rz(-1.2051219) q[0];
sx q[0];
rz(1.2980365) q[0];
rz(-0.80728665) q[1];
sx q[1];
rz(-1.9629982) q[1];
sx q[1];
rz(2.2198026) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4066276) q[0];
sx q[0];
rz(-1.4359183) q[0];
sx q[0];
rz(-0.69358967) q[0];
x q[1];
rz(-2.8554101) q[2];
sx q[2];
rz(-2.2879062) q[2];
sx q[2];
rz(-1.66695) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8340048) q[1];
sx q[1];
rz(-1.7760135) q[1];
sx q[1];
rz(-0.25968857) q[1];
x q[2];
rz(2.962035) q[3];
sx q[3];
rz(-0.87053821) q[3];
sx q[3];
rz(1.2683271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.90199295) q[2];
sx q[2];
rz(-0.56100503) q[2];
sx q[2];
rz(1.1716589) q[2];
rz(1.3575859) q[3];
sx q[3];
rz(-1.4566908) q[3];
sx q[3];
rz(1.6931036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25306025) q[0];
sx q[0];
rz(-1.9655515) q[0];
sx q[0];
rz(-1.4755479) q[0];
rz(-1.6015923) q[1];
sx q[1];
rz(-1.4614636) q[1];
sx q[1];
rz(-2.4618861) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4363791) q[0];
sx q[0];
rz(-1.0562911) q[0];
sx q[0];
rz(-1.3374469) q[0];
rz(-1.3557415) q[2];
sx q[2];
rz(-0.70677033) q[2];
sx q[2];
rz(-0.099345318) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.745984) q[1];
sx q[1];
rz(-1.2320227) q[1];
sx q[1];
rz(0.70396522) q[1];
x q[2];
rz(2.4623975) q[3];
sx q[3];
rz(-2.0742356) q[3];
sx q[3];
rz(-0.60590832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0150962) q[2];
sx q[2];
rz(-1.482778) q[2];
sx q[2];
rz(-0.91040197) q[2];
rz(-0.67534584) q[3];
sx q[3];
rz(-0.93674913) q[3];
sx q[3];
rz(-2.2909686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0891721) q[0];
sx q[0];
rz(-1.2344673) q[0];
sx q[0];
rz(-2.4269379) q[0];
rz(0.71406281) q[1];
sx q[1];
rz(-2.1866182) q[1];
sx q[1];
rz(-1.9649327) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0929716) q[0];
sx q[0];
rz(-1.4074568) q[0];
sx q[0];
rz(0.049157814) q[0];
rz(-pi) q[1];
rz(2.2048336) q[2];
sx q[2];
rz(-1.1144131) q[2];
sx q[2];
rz(-2.2697743) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1706108) q[1];
sx q[1];
rz(-2.5787528) q[1];
sx q[1];
rz(-1.7112205) q[1];
x q[2];
rz(-2.9858399) q[3];
sx q[3];
rz(-0.51625801) q[3];
sx q[3];
rz(1.2092276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.24361336) q[2];
sx q[2];
rz(-1.4695797) q[2];
sx q[2];
rz(-1.127355) q[2];
rz(-0.35774287) q[3];
sx q[3];
rz(-0.8711516) q[3];
sx q[3];
rz(1.0424785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
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
rz(-0.73666112) q[2];
sx q[2];
rz(-1.3623289) q[2];
sx q[2];
rz(1.3883432) q[2];
rz(-1.5626004) q[3];
sx q[3];
rz(-1.9215487) q[3];
sx q[3];
rz(2.4196845) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];