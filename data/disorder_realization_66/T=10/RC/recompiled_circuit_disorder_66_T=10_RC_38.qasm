OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.80823094) q[0];
sx q[0];
rz(4.1132676) q[0];
sx q[0];
rz(10.905807) q[0];
rz(-0.12227585) q[1];
sx q[1];
rz(-0.08637698) q[1];
sx q[1];
rz(0.02286214) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45863736) q[0];
sx q[0];
rz(-1.6874755) q[0];
sx q[0];
rz(-1.7200243) q[0];
rz(2.6267251) q[2];
sx q[2];
rz(-1.1866743) q[2];
sx q[2];
rz(-1.9799973) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.248977) q[1];
sx q[1];
rz(-1.9141344) q[1];
sx q[1];
rz(0.15906072) q[1];
rz(2.1894737) q[3];
sx q[3];
rz(-1.7784356) q[3];
sx q[3];
rz(-1.4116956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7227398) q[2];
sx q[2];
rz(-2.3143694) q[2];
sx q[2];
rz(2.375405) q[2];
rz(0.17151672) q[3];
sx q[3];
rz(-2.4092716) q[3];
sx q[3];
rz(-2.5550301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(2.3294753) q[0];
sx q[0];
rz(-2.7024039) q[0];
sx q[0];
rz(0.086659327) q[0];
rz(-0.86241972) q[1];
sx q[1];
rz(-0.54310596) q[1];
sx q[1];
rz(-3.0564953) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25227308) q[0];
sx q[0];
rz(-0.62250455) q[0];
sx q[0];
rz(3.1134393) q[0];
x q[1];
rz(-0.30479635) q[2];
sx q[2];
rz(-2.1509503) q[2];
sx q[2];
rz(-0.09588974) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9566006) q[1];
sx q[1];
rz(-1.6829832) q[1];
sx q[1];
rz(-0.086602028) q[1];
x q[2];
rz(2.0410791) q[3];
sx q[3];
rz(-1.2840052) q[3];
sx q[3];
rz(1.6581397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.26248419) q[2];
sx q[2];
rz(-1.5905249) q[2];
sx q[2];
rz(-1.6764199) q[2];
rz(1.1761459) q[3];
sx q[3];
rz(-2.2361103) q[3];
sx q[3];
rz(-2.8951077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48433205) q[0];
sx q[0];
rz(-3.0069139) q[0];
sx q[0];
rz(-2.235967) q[0];
rz(0.33272818) q[1];
sx q[1];
rz(-0.85488027) q[1];
sx q[1];
rz(-1.3844301) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6838283) q[0];
sx q[0];
rz(-2.0450852) q[0];
sx q[0];
rz(1.6636687) q[0];
rz(-2.2134174) q[2];
sx q[2];
rz(-2.536142) q[2];
sx q[2];
rz(-1.5290934) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3123734) q[1];
sx q[1];
rz(-0.81554283) q[1];
sx q[1];
rz(-1.6014598) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3178187) q[3];
sx q[3];
rz(-0.71411055) q[3];
sx q[3];
rz(-2.2658474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.24421346) q[2];
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
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0760647) q[0];
sx q[0];
rz(-0.17834839) q[0];
sx q[0];
rz(-2.4705825) q[0];
rz(-1.7975851) q[1];
sx q[1];
rz(-1.1616511) q[1];
sx q[1];
rz(-0.20733325) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5642693) q[0];
sx q[0];
rz(-0.042357001) q[0];
sx q[0];
rz(1.7839412) q[0];
rz(-pi) q[1];
rz(-1.3992594) q[2];
sx q[2];
rz(-0.17645391) q[2];
sx q[2];
rz(1.2982969) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.281665) q[1];
sx q[1];
rz(-1.4178226) q[1];
sx q[1];
rz(-1.330919) q[1];
rz(-1.585235) q[3];
sx q[3];
rz(-1.0591905) q[3];
sx q[3];
rz(-2.5479941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8308782) q[2];
sx q[2];
rz(-2.0534616) q[2];
sx q[2];
rz(-0.80835289) q[2];
rz(0.80777848) q[3];
sx q[3];
rz(-0.72179663) q[3];
sx q[3];
rz(-2.9479153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59072524) q[0];
sx q[0];
rz(-0.71722513) q[0];
sx q[0];
rz(-0.89548683) q[0];
rz(0.26516178) q[1];
sx q[1];
rz(-2.3064955) q[1];
sx q[1];
rz(0.53363824) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1874439) q[0];
sx q[0];
rz(-1.8549518) q[0];
sx q[0];
rz(-2.8548293) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9075378) q[2];
sx q[2];
rz(-0.78163994) q[2];
sx q[2];
rz(-0.78125886) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6753908) q[1];
sx q[1];
rz(-1.4247822) q[1];
sx q[1];
rz(2.9195665) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8495464) q[3];
sx q[3];
rz(-2.4933443) q[3];
sx q[3];
rz(-1.3384782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9020033) q[2];
sx q[2];
rz(-1.5633554) q[2];
sx q[2];
rz(-0.56838244) q[2];
rz(-1.2166294) q[3];
sx q[3];
rz(-2.670848) q[3];
sx q[3];
rz(0.34354982) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50826532) q[0];
sx q[0];
rz(-0.87843043) q[0];
sx q[0];
rz(-1.1761965) q[0];
rz(1.5559224) q[1];
sx q[1];
rz(-2.1891179) q[1];
sx q[1];
rz(2.1369381) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4906824) q[0];
sx q[0];
rz(-0.51603979) q[0];
sx q[0];
rz(-2.3923621) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7447694) q[2];
sx q[2];
rz(-1.4962215) q[2];
sx q[2];
rz(0.65149388) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1130044) q[1];
sx q[1];
rz(-0.60739809) q[1];
sx q[1];
rz(0.44087704) q[1];
x q[2];
rz(0.011990487) q[3];
sx q[3];
rz(-1.1569835) q[3];
sx q[3];
rz(-0.12366611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4040318) q[2];
sx q[2];
rz(-0.99016756) q[2];
sx q[2];
rz(2.693434) q[2];
rz(-0.29799497) q[3];
sx q[3];
rz(-2.4500676) q[3];
sx q[3];
rz(-0.34860778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3649243) q[0];
sx q[0];
rz(-1.3013327) q[0];
sx q[0];
rz(-0.65814322) q[0];
rz(1.8572042) q[1];
sx q[1];
rz(-2.6874266) q[1];
sx q[1];
rz(2.4694494) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9400529) q[0];
sx q[0];
rz(-1.4048028) q[0];
sx q[0];
rz(0.0075017651) q[0];
rz(-0.42598287) q[2];
sx q[2];
rz(-1.1642712) q[2];
sx q[2];
rz(-0.68824088) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8183648) q[1];
sx q[1];
rz(-1.0867026) q[1];
sx q[1];
rz(1.9869884) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7147195) q[3];
sx q[3];
rz(-0.87464273) q[3];
sx q[3];
rz(1.6203984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5333574) q[2];
sx q[2];
rz(-2.50864) q[2];
sx q[2];
rz(-0.39880025) q[2];
rz(0.82459015) q[3];
sx q[3];
rz(-1.7493068) q[3];
sx q[3];
rz(-2.5430172) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7010715) q[0];
sx q[0];
rz(-2.7109787) q[0];
sx q[0];
rz(-0.15821247) q[0];
rz(-2.054706) q[1];
sx q[1];
rz(-0.67651665) q[1];
sx q[1];
rz(0.80668443) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8950302) q[0];
sx q[0];
rz(-0.44996214) q[0];
sx q[0];
rz(1.8438086) q[0];
x q[1];
rz(0.20571795) q[2];
sx q[2];
rz(-1.4652025) q[2];
sx q[2];
rz(0.92668698) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1572691) q[1];
sx q[1];
rz(-1.2829797) q[1];
sx q[1];
rz(2.4411574) q[1];
x q[2];
rz(-1.673647) q[3];
sx q[3];
rz(-0.30234435) q[3];
sx q[3];
rz(0.24149382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5640101) q[2];
sx q[2];
rz(-1.6489886) q[2];
sx q[2];
rz(-0.99964833) q[2];
rz(-3.0380761) q[3];
sx q[3];
rz(-3.0045356) q[3];
sx q[3];
rz(-2.0172393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0209811) q[0];
sx q[0];
rz(-0.7779026) q[0];
sx q[0];
rz(2.4840684) q[0];
rz(-2.855037) q[1];
sx q[1];
rz(-2.2189238) q[1];
sx q[1];
rz(-2.8009169) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1239615) q[0];
sx q[0];
rz(-2.258856) q[0];
sx q[0];
rz(-2.7244363) q[0];
rz(-2.1358143) q[2];
sx q[2];
rz(-1.3823595) q[2];
sx q[2];
rz(2.4609158) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9790736) q[1];
sx q[1];
rz(-2.5538429) q[1];
sx q[1];
rz(1.9025004) q[1];
x q[2];
rz(-2.4386028) q[3];
sx q[3];
rz(-0.25745108) q[3];
sx q[3];
rz(-1.6182181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.75491536) q[2];
sx q[2];
rz(-1.9415386) q[2];
sx q[2];
rz(-2.6436451) q[2];
rz(-0.30161101) q[3];
sx q[3];
rz(-2.7409654) q[3];
sx q[3];
rz(-0.51914674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4122445) q[0];
sx q[0];
rz(-0.059818581) q[0];
sx q[0];
rz(2.8630032) q[0];
rz(2.5623698) q[1];
sx q[1];
rz(-0.93943739) q[1];
sx q[1];
rz(-3.0648807) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3007606) q[0];
sx q[0];
rz(-1.4474086) q[0];
sx q[0];
rz(1.2542017) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8329343) q[2];
sx q[2];
rz(-0.58912504) q[2];
sx q[2];
rz(2.4208456) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0753239) q[1];
sx q[1];
rz(-1.2053718) q[1];
sx q[1];
rz(1.5484023) q[1];
rz(1.65927) q[3];
sx q[3];
rz(-1.3946597) q[3];
sx q[3];
rz(-2.4000771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3672093) q[2];
sx q[2];
rz(-0.70720208) q[2];
sx q[2];
rz(0.45551604) q[2];
rz(-0.43141836) q[3];
sx q[3];
rz(-3.016267) q[3];
sx q[3];
rz(3.07807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1800304) q[0];
sx q[0];
rz(-1.7407692) q[0];
sx q[0];
rz(-2.2041007) q[0];
rz(-2.5554399) q[1];
sx q[1];
rz(-1.502232) q[1];
sx q[1];
rz(-1.4486817) q[1];
rz(1.099844) q[2];
sx q[2];
rz(-1.3146613) q[2];
sx q[2];
rz(1.4801499) q[2];
rz(1.9361817) q[3];
sx q[3];
rz(-1.1689417) q[3];
sx q[3];
rz(0.4369215) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];