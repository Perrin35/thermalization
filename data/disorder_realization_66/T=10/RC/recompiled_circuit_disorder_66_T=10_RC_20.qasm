OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3333617) q[0];
sx q[0];
rz(-0.97167492) q[0];
sx q[0];
rz(-1.4810286) q[0];
rz(3.0193168) q[1];
sx q[1];
rz(-3.0552157) q[1];
sx q[1];
rz(3.1187305) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6883144) q[0];
sx q[0];
rz(-0.18916057) q[0];
sx q[0];
rz(-0.90318371) q[0];
rz(0.68732287) q[2];
sx q[2];
rz(-0.63185531) q[2];
sx q[2];
rz(-2.1473715) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.80427985) q[1];
sx q[1];
rz(-0.37706456) q[1];
sx q[1];
rz(1.9878597) q[1];
x q[2];
rz(-2.1894737) q[3];
sx q[3];
rz(-1.3631571) q[3];
sx q[3];
rz(-1.4116956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.4188529) q[2];
sx q[2];
rz(-2.3143694) q[2];
sx q[2];
rz(2.375405) q[2];
rz(0.17151672) q[3];
sx q[3];
rz(-2.4092716) q[3];
sx q[3];
rz(0.58656251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8121174) q[0];
sx q[0];
rz(-0.43918878) q[0];
sx q[0];
rz(0.086659327) q[0];
rz(-0.86241972) q[1];
sx q[1];
rz(-2.5984867) q[1];
sx q[1];
rz(-0.08509732) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.800195) q[0];
sx q[0];
rz(-1.5872103) q[0];
sx q[0];
rz(2.5192758) q[0];
rz(-2.8367963) q[2];
sx q[2];
rz(-0.99064231) q[2];
sx q[2];
rz(-0.09588974) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3955235) q[1];
sx q[1];
rz(-1.4847401) q[1];
sx q[1];
rz(-1.6834016) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8220903) q[3];
sx q[3];
rz(-2.0204244) q[3];
sx q[3];
rz(-0.055469661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8791085) q[2];
sx q[2];
rz(-1.5510677) q[2];
sx q[2];
rz(1.6764199) q[2];
rz(1.1761459) q[3];
sx q[3];
rz(-0.90548235) q[3];
sx q[3];
rz(-0.24648497) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48433205) q[0];
sx q[0];
rz(-3.0069139) q[0];
sx q[0];
rz(0.9056257) q[0];
rz(-0.33272818) q[1];
sx q[1];
rz(-2.2867124) q[1];
sx q[1];
rz(1.7571626) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9860486) q[0];
sx q[0];
rz(-1.4882003) q[0];
sx q[0];
rz(0.47604576) q[0];
rz(-pi) q[1];
x q[1];
rz(2.076782) q[2];
sx q[2];
rz(-1.9188606) q[2];
sx q[2];
rz(0.51007523) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3789931) q[1];
sx q[1];
rz(-1.5484719) q[1];
sx q[1];
rz(0.75548817) q[1];
rz(2.1372041) q[3];
sx q[3];
rz(-2.0319788) q[3];
sx q[3];
rz(-0.084463488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8973792) q[2];
sx q[2];
rz(-2.801557) q[2];
sx q[2];
rz(-1.8133694) q[2];
rz(1.3978847) q[3];
sx q[3];
rz(-0.54026794) q[3];
sx q[3];
rz(3.0298997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0760647) q[0];
sx q[0];
rz(-0.17834839) q[0];
sx q[0];
rz(-0.67101014) q[0];
rz(1.3440075) q[1];
sx q[1];
rz(-1.9799415) q[1];
sx q[1];
rz(-2.9342594) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5642693) q[0];
sx q[0];
rz(-0.042357001) q[0];
sx q[0];
rz(1.7839412) q[0];
x q[1];
rz(-0.030427293) q[2];
sx q[2];
rz(-1.3969587) q[2];
sx q[2];
rz(1.1241084) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3934717) q[1];
sx q[1];
rz(-1.3337743) q[1];
sx q[1];
rz(-0.15740983) q[1];
x q[2];
rz(-2.6299423) q[3];
sx q[3];
rz(-1.5833862) q[3];
sx q[3];
rz(-2.1714641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.31071445) q[2];
sx q[2];
rz(-2.0534616) q[2];
sx q[2];
rz(-0.80835289) q[2];
rz(2.3338142) q[3];
sx q[3];
rz(-0.72179663) q[3];
sx q[3];
rz(2.9479153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59072524) q[0];
sx q[0];
rz(-2.4243675) q[0];
sx q[0];
rz(2.2461058) q[0];
rz(2.8764309) q[1];
sx q[1];
rz(-2.3064955) q[1];
sx q[1];
rz(-0.53363824) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69913188) q[0];
sx q[0];
rz(-1.8457544) q[0];
sx q[0];
rz(-1.2752227) q[0];
rz(-pi) q[1];
rz(-0.76782121) q[2];
sx q[2];
rz(-1.7349093) q[2];
sx q[2];
rz(0.62190157) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0698318) q[1];
sx q[1];
rz(-1.3511718) q[1];
sx q[1];
rz(1.7204309) q[1];
rz(-pi) q[2];
rz(1.8495464) q[3];
sx q[3];
rz(-2.4933443) q[3];
sx q[3];
rz(1.3384782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9020033) q[2];
sx q[2];
rz(-1.5782372) q[2];
sx q[2];
rz(0.56838244) q[2];
rz(-1.9249632) q[3];
sx q[3];
rz(-2.670848) q[3];
sx q[3];
rz(-0.34354982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6333273) q[0];
sx q[0];
rz(-0.87843043) q[0];
sx q[0];
rz(1.9653962) q[0];
rz(1.5559224) q[1];
sx q[1];
rz(-2.1891179) q[1];
sx q[1];
rz(2.1369381) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83197901) q[0];
sx q[0];
rz(-1.9404611) q[0];
sx q[0];
rz(-1.2020822) q[0];
rz(-0.3968233) q[2];
sx q[2];
rz(-1.6453711) q[2];
sx q[2];
rz(2.4900988) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1130044) q[1];
sx q[1];
rz(-0.60739809) q[1];
sx q[1];
rz(0.44087704) q[1];
rz(-pi) q[2];
x q[2];
rz(0.011990487) q[3];
sx q[3];
rz(-1.1569835) q[3];
sx q[3];
rz(-0.12366611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4040318) q[2];
sx q[2];
rz(-2.1514251) q[2];
sx q[2];
rz(-0.44815865) q[2];
rz(2.8435977) q[3];
sx q[3];
rz(-0.69152504) q[3];
sx q[3];
rz(0.34860778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
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
rz(0.77666831) q[0];
sx q[0];
rz(-1.8402599) q[0];
sx q[0];
rz(-0.65814322) q[0];
rz(-1.8572042) q[1];
sx q[1];
rz(-2.6874266) q[1];
sx q[1];
rz(0.67214322) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2015398) q[0];
sx q[0];
rz(-1.7367898) q[0];
sx q[0];
rz(0.0075017651) q[0];
rz(-pi) q[1];
rz(2.3357046) q[2];
sx q[2];
rz(-0.58008367) q[2];
sx q[2];
rz(0.16575955) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0587412) q[1];
sx q[1];
rz(-0.62742678) q[1];
sx q[1];
rz(0.65545603) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4268731) q[3];
sx q[3];
rz(-2.2669499) q[3];
sx q[3];
rz(1.6203984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5333574) q[2];
sx q[2];
rz(-0.63295263) q[2];
sx q[2];
rz(0.39880025) q[2];
rz(0.82459015) q[3];
sx q[3];
rz(-1.7493068) q[3];
sx q[3];
rz(0.59857541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.44052112) q[0];
sx q[0];
rz(-2.7109787) q[0];
sx q[0];
rz(-2.9833802) q[0];
rz(2.054706) q[1];
sx q[1];
rz(-0.67651665) q[1];
sx q[1];
rz(2.3349082) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0772484) q[0];
sx q[0];
rz(-1.4532538) q[0];
sx q[0];
rz(2.006152) q[0];
x q[1];
rz(-0.20571795) q[2];
sx q[2];
rz(-1.4652025) q[2];
sx q[2];
rz(-0.92668698) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.73831564) q[1];
sx q[1];
rz(-0.74790955) q[1];
sx q[1];
rz(-0.43055375) q[1];
rz(-1.4679457) q[3];
sx q[3];
rz(-0.30234435) q[3];
sx q[3];
rz(-0.24149382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.57758254) q[2];
sx q[2];
rz(-1.492604) q[2];
sx q[2];
rz(-0.99964833) q[2];
rz(3.0380761) q[3];
sx q[3];
rz(-0.13705702) q[3];
sx q[3];
rz(1.1243533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12061159) q[0];
sx q[0];
rz(-2.3636901) q[0];
sx q[0];
rz(2.4840684) q[0];
rz(0.28655562) q[1];
sx q[1];
rz(-0.92266881) q[1];
sx q[1];
rz(-0.34067571) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17250241) q[0];
sx q[0];
rz(-1.2524676) q[0];
sx q[0];
rz(0.83842917) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1358143) q[2];
sx q[2];
rz(-1.3823595) q[2];
sx q[2];
rz(-2.4609158) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.68742467) q[1];
sx q[1];
rz(-1.3892281) q[1];
sx q[1];
rz(2.1329692) q[1];
x q[2];
rz(1.4021923) q[3];
sx q[3];
rz(-1.7662893) q[3];
sx q[3];
rz(-0.80381264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3866773) q[2];
sx q[2];
rz(-1.200054) q[2];
sx q[2];
rz(2.6436451) q[2];
rz(-2.8399816) q[3];
sx q[3];
rz(-2.7409654) q[3];
sx q[3];
rz(0.51914674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4122445) q[0];
sx q[0];
rz(-0.059818581) q[0];
sx q[0];
rz(0.27858946) q[0];
rz(-0.57922286) q[1];
sx q[1];
rz(-0.93943739) q[1];
sx q[1];
rz(-3.0648807) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0522703) q[0];
sx q[0];
rz(-2.8025586) q[0];
sx q[0];
rz(-1.9498755) q[0];
x q[1];
rz(1.3086583) q[2];
sx q[2];
rz(-2.5524676) q[2];
sx q[2];
rz(-2.4208456) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0127276) q[1];
sx q[1];
rz(-2.7755133) q[1];
sx q[1];
rz(-3.0831343) q[1];
rz(-2.9647787) q[3];
sx q[3];
rz(-1.6578976) q[3];
sx q[3];
rz(2.2967695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77438337) q[2];
sx q[2];
rz(-2.4343906) q[2];
sx q[2];
rz(-2.6860766) q[2];
rz(-0.43141836) q[3];
sx q[3];
rz(-3.016267) q[3];
sx q[3];
rz(3.07807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1800304) q[0];
sx q[0];
rz(-1.4008235) q[0];
sx q[0];
rz(0.93749198) q[0];
rz(-0.58615276) q[1];
sx q[1];
rz(-1.6393607) q[1];
sx q[1];
rz(1.6929109) q[1];
rz(0.2858328) q[2];
sx q[2];
rz(-1.1163859) q[2];
sx q[2];
rz(-3.1039539) q[2];
rz(-0.69910819) q[3];
sx q[3];
rz(-0.53634488) q[3];
sx q[3];
rz(2.8041822) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
