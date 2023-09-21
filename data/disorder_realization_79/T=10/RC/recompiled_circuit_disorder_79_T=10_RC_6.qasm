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
rz(7.1516501) q[0];
sx q[0];
rz(9.2317543) q[0];
rz(-1.9999737) q[1];
sx q[1];
rz(-2.7116099) q[1];
sx q[1];
rz(-2.4584682) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0177512) q[0];
sx q[0];
rz(-0.089028247) q[0];
sx q[0];
rz(-2.428399) q[0];
x q[1];
rz(-1.2201266) q[2];
sx q[2];
rz(-1.3436683) q[2];
sx q[2];
rz(-2.7820058) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.19385399) q[1];
sx q[1];
rz(-2.1734936) q[1];
sx q[1];
rz(2.1147453) q[1];
rz(-2.7636823) q[3];
sx q[3];
rz(-2.1034735) q[3];
sx q[3];
rz(-2.7045254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1258939) q[2];
sx q[2];
rz(-1.3796207) q[2];
sx q[2];
rz(-1.0985628) q[2];
rz(-2.0627608) q[3];
sx q[3];
rz(-2.1964985) q[3];
sx q[3];
rz(-1.1014972) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1612448) q[0];
sx q[0];
rz(-0.315936) q[0];
sx q[0];
rz(-2.9336477) q[0];
rz(0.57693276) q[1];
sx q[1];
rz(-2.2536342) q[1];
sx q[1];
rz(1.4651441) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0365021) q[0];
sx q[0];
rz(-0.72421342) q[0];
sx q[0];
rz(-0.0037395517) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0287839) q[2];
sx q[2];
rz(-2.6070242) q[2];
sx q[2];
rz(0.88876681) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1190471) q[1];
sx q[1];
rz(-0.76821583) q[1];
sx q[1];
rz(2.4382298) q[1];
rz(-pi) q[2];
rz(-0.74554262) q[3];
sx q[3];
rz(-1.8288444) q[3];
sx q[3];
rz(-2.6170956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1229822) q[2];
sx q[2];
rz(-0.98615042) q[2];
sx q[2];
rz(-3.0318276) q[2];
rz(2.5189853) q[3];
sx q[3];
rz(-0.37125769) q[3];
sx q[3];
rz(-1.6842779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96317545) q[0];
sx q[0];
rz(-0.85997471) q[0];
sx q[0];
rz(2.6254568) q[0];
rz(2.5667045) q[1];
sx q[1];
rz(-0.92620414) q[1];
sx q[1];
rz(-2.3410472) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.752906) q[0];
sx q[0];
rz(-1.1645002) q[0];
sx q[0];
rz(-2.9857062) q[0];
rz(-2.4750701) q[2];
sx q[2];
rz(-0.98072532) q[2];
sx q[2];
rz(0.18596622) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.87677466) q[1];
sx q[1];
rz(-2.7285828) q[1];
sx q[1];
rz(-3.1289711) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.3018467) q[3];
sx q[3];
rz(-2.2491124) q[3];
sx q[3];
rz(2.9523926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.67733726) q[2];
sx q[2];
rz(-2.8223473) q[2];
sx q[2];
rz(-1.3595954) q[2];
rz(-0.96757403) q[3];
sx q[3];
rz(-1.8680957) q[3];
sx q[3];
rz(-1.4250071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99252218) q[0];
sx q[0];
rz(-1.8795805) q[0];
sx q[0];
rz(2.676679) q[0];
rz(-2.7930296) q[1];
sx q[1];
rz(-0.26270738) q[1];
sx q[1];
rz(1.0850614) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5497919) q[0];
sx q[0];
rz(-1.340197) q[0];
sx q[0];
rz(2.6016298) q[0];
x q[1];
rz(2.8934946) q[2];
sx q[2];
rz(-1.9271701) q[2];
sx q[2];
rz(-0.46596371) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3738721) q[1];
sx q[1];
rz(-1.4704736) q[1];
sx q[1];
rz(-2.9994681) q[1];
rz(-1.9784847) q[3];
sx q[3];
rz(-0.99916047) q[3];
sx q[3];
rz(0.15743263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.00502914) q[2];
sx q[2];
rz(-2.0932784) q[2];
sx q[2];
rz(-0.68112779) q[2];
rz(-0.51182169) q[3];
sx q[3];
rz(-0.32326439) q[3];
sx q[3];
rz(-2.8779023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8624449) q[0];
sx q[0];
rz(-2.6036766) q[0];
sx q[0];
rz(1.408668) q[0];
rz(-2.7092343) q[1];
sx q[1];
rz(-0.84195781) q[1];
sx q[1];
rz(0.98168215) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3571346) q[0];
sx q[0];
rz(-2.031209) q[0];
sx q[0];
rz(-0.61607342) q[0];
x q[1];
rz(1.8129187) q[2];
sx q[2];
rz(-1.7244581) q[2];
sx q[2];
rz(2.9850933) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7179467) q[1];
sx q[1];
rz(-1.7182933) q[1];
sx q[1];
rz(-0.46848483) q[1];
rz(-3.0675689) q[3];
sx q[3];
rz(-2.166966) q[3];
sx q[3];
rz(-2.1285469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1061873) q[2];
sx q[2];
rz(-1.4865439) q[2];
sx q[2];
rz(2.6521818) q[2];
rz(-2.126746) q[3];
sx q[3];
rz(-0.5265407) q[3];
sx q[3];
rz(-1.813252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6569825) q[0];
sx q[0];
rz(-0.15743142) q[0];
sx q[0];
rz(-0.37242517) q[0];
rz(-1.8107481) q[1];
sx q[1];
rz(-2.1061888) q[1];
sx q[1];
rz(0.16528027) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0342456) q[0];
sx q[0];
rz(-0.099485569) q[0];
sx q[0];
rz(0.29883595) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0566453) q[2];
sx q[2];
rz(-0.68945486) q[2];
sx q[2];
rz(1.2598318) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2416934) q[1];
sx q[1];
rz(-2.4447933) q[1];
sx q[1];
rz(0.65710575) q[1];
x q[2];
rz(-0.23541707) q[3];
sx q[3];
rz(-1.132292) q[3];
sx q[3];
rz(-0.32270839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.91339397) q[2];
sx q[2];
rz(-1.0281576) q[2];
sx q[2];
rz(1.4432663) q[2];
rz(0.36744395) q[3];
sx q[3];
rz(-1.8668709) q[3];
sx q[3];
rz(0.2120367) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3570324) q[0];
sx q[0];
rz(-1.0914047) q[0];
sx q[0];
rz(2.7923287) q[0];
rz(-0.7473942) q[1];
sx q[1];
rz(-0.29577297) q[1];
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
rz(1.6197617) q[0];
x q[1];
rz(-2.1963504) q[2];
sx q[2];
rz(-0.40996273) q[2];
sx q[2];
rz(-1.2298825) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0349717) q[1];
sx q[1];
rz(-1.7502516) q[1];
sx q[1];
rz(-0.76645318) q[1];
x q[2];
rz(1.1398846) q[3];
sx q[3];
rz(-2.0322554) q[3];
sx q[3];
rz(-2.0737089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0050469) q[2];
sx q[2];
rz(-1.9133291) q[2];
sx q[2];
rz(-2.6100256) q[2];
rz(-0.68938869) q[3];
sx q[3];
rz(-1.4644943) q[3];
sx q[3];
rz(-1.2954856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
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
rz(-0.73496504) q[0];
sx q[0];
rz(-1.7056744) q[0];
sx q[0];
rz(0.69358967) q[0];
x q[1];
rz(1.8838896) q[2];
sx q[2];
rz(-0.76258341) q[2];
sx q[2];
rz(1.8956172) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8242952) q[1];
sx q[1];
rz(-1.824914) q[1];
sx q[1];
rz(1.7829249) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7796302) q[3];
sx q[3];
rz(-0.71912557) q[3];
sx q[3];
rz(-1.5428839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2395997) q[2];
sx q[2];
rz(-0.56100503) q[2];
sx q[2];
rz(1.9699338) q[2];
rz(1.7840067) q[3];
sx q[3];
rz(-1.4566908) q[3];
sx q[3];
rz(1.4484891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25306025) q[0];
sx q[0];
rz(-1.9655515) q[0];
sx q[0];
rz(1.6660447) q[0];
rz(1.6015923) q[1];
sx q[1];
rz(-1.4614636) q[1];
sx q[1];
rz(-0.67970651) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4363791) q[0];
sx q[0];
rz(-1.0562911) q[0];
sx q[0];
rz(-1.8041457) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9613413) q[2];
sx q[2];
rz(-2.2580574) q[2];
sx q[2];
rz(-0.3790516) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9437127) q[1];
sx q[1];
rz(-2.3730952) q[1];
sx q[1];
rz(-0.49853034) q[1];
x q[2];
rz(-0.9548095) q[3];
sx q[3];
rz(-2.1534352) q[3];
sx q[3];
rz(0.59347502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1264964) q[2];
sx q[2];
rz(-1.6588147) q[2];
sx q[2];
rz(-2.2311907) q[2];
rz(-0.67534584) q[3];
sx q[3];
rz(-2.2048435) q[3];
sx q[3];
rz(-0.85062406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-3.0891721) q[0];
sx q[0];
rz(-1.2344673) q[0];
sx q[0];
rz(-2.4269379) q[0];
rz(-2.4275298) q[1];
sx q[1];
rz(-2.1866182) q[1];
sx q[1];
rz(1.1766599) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51417527) q[0];
sx q[0];
rz(-1.5222933) q[0];
sx q[0];
rz(-1.7343299) q[0];
x q[1];
rz(-2.5942957) q[2];
sx q[2];
rz(-2.1314869) q[2];
sx q[2];
rz(-1.0123569) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0050051) q[1];
sx q[1];
rz(-2.1274381) q[1];
sx q[1];
rz(-0.088076061) q[1];
rz(-pi) q[2];
rz(-0.51104607) q[3];
sx q[3];
rz(-1.4941477) q[3];
sx q[3];
rz(0.49728909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8979793) q[2];
sx q[2];
rz(-1.4695797) q[2];
sx q[2];
rz(-1.127355) q[2];
rz(0.35774287) q[3];
sx q[3];
rz(-2.2704411) q[3];
sx q[3];
rz(1.0424785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42416278) q[0];
sx q[0];
rz(-1.2107727) q[0];
sx q[0];
rz(-2.6834224) q[0];
rz(-0.39623109) q[1];
sx q[1];
rz(-3.1165262) q[1];
sx q[1];
rz(-2.810626) q[1];
rz(-1.8489807) q[2];
sx q[2];
rz(-2.2879911) q[2];
sx q[2];
rz(0.0030980274) q[2];
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
