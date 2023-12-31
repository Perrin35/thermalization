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
rz(1.6605641) q[0];
rz(-0.12227585) q[1];
sx q[1];
rz(-0.08637698) q[1];
sx q[1];
rz(-3.1187305) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6829553) q[0];
sx q[0];
rz(-1.6874755) q[0];
sx q[0];
rz(-1.4215683) q[0];
x q[1];
rz(0.51486751) q[2];
sx q[2];
rz(-1.9549184) q[2];
sx q[2];
rz(1.1615953) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.80427985) q[1];
sx q[1];
rz(-0.37706456) q[1];
sx q[1];
rz(1.1537329) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25306563) q[3];
sx q[3];
rz(-2.1742637) q[3];
sx q[3];
rz(2.8367708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.4188529) q[2];
sx q[2];
rz(-0.8272233) q[2];
sx q[2];
rz(0.76618761) q[2];
rz(2.9700759) q[3];
sx q[3];
rz(-2.4092716) q[3];
sx q[3];
rz(-0.58656251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8121174) q[0];
sx q[0];
rz(-2.7024039) q[0];
sx q[0];
rz(-0.086659327) q[0];
rz(2.2791729) q[1];
sx q[1];
rz(-2.5984867) q[1];
sx q[1];
rz(-0.08509732) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9239685) q[0];
sx q[0];
rz(-2.1930165) q[0];
sx q[0];
rz(1.5505962) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1413967) q[2];
sx q[2];
rz(-2.4944802) q[2];
sx q[2];
rz(0.61691689) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9566006) q[1];
sx q[1];
rz(-1.4586095) q[1];
sx q[1];
rz(-3.0549906) q[1];
rz(-pi) q[2];
rz(2.1477633) q[3];
sx q[3];
rz(-2.5964063) q[3];
sx q[3];
rz(2.546437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8791085) q[2];
sx q[2];
rz(-1.5905249) q[2];
sx q[2];
rz(1.6764199) q[2];
rz(-1.1761459) q[3];
sx q[3];
rz(-0.90548235) q[3];
sx q[3];
rz(-2.8951077) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6572606) q[0];
sx q[0];
rz(-3.0069139) q[0];
sx q[0];
rz(2.235967) q[0];
rz(0.33272818) q[1];
sx q[1];
rz(-0.85488027) q[1];
sx q[1];
rz(1.7571626) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6838283) q[0];
sx q[0];
rz(-2.0450852) q[0];
sx q[0];
rz(-1.6636687) q[0];
rz(-pi) q[1];
rz(-2.7483814) q[2];
sx q[2];
rz(-1.097743) q[2];
sx q[2];
rz(-0.87393239) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.76259957) q[1];
sx q[1];
rz(-1.5931207) q[1];
sx q[1];
rz(-2.3861045) q[1];
rz(2.3178187) q[3];
sx q[3];
rz(-0.71411055) q[3];
sx q[3];
rz(2.2658474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8973792) q[2];
sx q[2];
rz(-2.801557) q[2];
sx q[2];
rz(-1.3282233) q[2];
rz(1.7437079) q[3];
sx q[3];
rz(-2.6013247) q[3];
sx q[3];
rz(-0.11169294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0760647) q[0];
sx q[0];
rz(-0.17834839) q[0];
sx q[0];
rz(0.67101014) q[0];
rz(-1.7975851) q[1];
sx q[1];
rz(-1.1616511) q[1];
sx q[1];
rz(-0.20733325) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78051356) q[0];
sx q[0];
rz(-1.5797537) q[0];
sx q[0];
rz(-1.6121959) q[0];
rz(-pi) q[1];
rz(-1.7447128) q[2];
sx q[2];
rz(-1.6007649) q[2];
sx q[2];
rz(2.7001691) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.748121) q[1];
sx q[1];
rz(-1.3337743) q[1];
sx q[1];
rz(-0.15740983) q[1];
x q[2];
rz(0.025709318) q[3];
sx q[3];
rz(-2.6298012) q[3];
sx q[3];
rz(2.5185086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.5508674) q[0];
sx q[0];
rz(-0.71722513) q[0];
sx q[0];
rz(0.89548683) q[0];
rz(-2.8764309) q[1];
sx q[1];
rz(-2.3064955) q[1];
sx q[1];
rz(-2.6079544) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4424608) q[0];
sx q[0];
rz(-1.8457544) q[0];
sx q[0];
rz(-1.2752227) q[0];
x q[1];
rz(1.7970423) q[2];
sx q[2];
rz(-2.3257253) q[2];
sx q[2];
rz(2.0362542) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.071760885) q[1];
sx q[1];
rz(-1.7904209) q[1];
sx q[1];
rz(1.4211618) q[1];
rz(-pi) q[2];
x q[2];
rz(0.20547262) q[3];
sx q[3];
rz(-2.1901202) q[3];
sx q[3];
rz(-1.6831786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2395893) q[2];
sx q[2];
rz(-1.5782372) q[2];
sx q[2];
rz(-2.5732102) q[2];
rz(-1.2166294) q[3];
sx q[3];
rz(-0.47074461) q[3];
sx q[3];
rz(-0.34354982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6333273) q[0];
sx q[0];
rz(-0.87843043) q[0];
sx q[0];
rz(-1.9653962) q[0];
rz(-1.5856702) q[1];
sx q[1];
rz(-0.95247477) q[1];
sx q[1];
rz(-2.1369381) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5414808) q[0];
sx q[0];
rz(-1.2280557) q[0];
sx q[0];
rz(-2.747885) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7447694) q[2];
sx q[2];
rz(-1.6453711) q[2];
sx q[2];
rz(0.65149388) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.028588258) q[1];
sx q[1];
rz(-0.60739809) q[1];
sx q[1];
rz(-0.44087704) q[1];
rz(-1.9846356) q[3];
sx q[3];
rz(-1.5817747) q[3];
sx q[3];
rz(-1.4423086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4040318) q[2];
sx q[2];
rz(-2.1514251) q[2];
sx q[2];
rz(-2.693434) q[2];
rz(0.29799497) q[3];
sx q[3];
rz(-0.69152504) q[3];
sx q[3];
rz(-0.34860778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.77666831) q[0];
sx q[0];
rz(-1.3013327) q[0];
sx q[0];
rz(-0.65814322) q[0];
rz(-1.8572042) q[1];
sx q[1];
rz(-0.45416608) q[1];
sx q[1];
rz(-0.67214322) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9400529) q[0];
sx q[0];
rz(-1.4048028) q[0];
sx q[0];
rz(0.0075017651) q[0];
x q[1];
rz(-1.1291814) q[2];
sx q[2];
rz(-1.9600944) q[2];
sx q[2];
rz(2.0814975) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0969095) q[1];
sx q[1];
rz(-1.9367366) q[1];
sx q[1];
rz(0.52176042) q[1];
rz(-pi) q[2];
rz(-0.16996202) q[3];
sx q[3];
rz(-0.70843491) q[3];
sx q[3];
rz(-1.8426614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5333574) q[2];
sx q[2];
rz(-2.50864) q[2];
sx q[2];
rz(0.39880025) q[2];
rz(-0.82459015) q[3];
sx q[3];
rz(-1.7493068) q[3];
sx q[3];
rz(2.5430172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7010715) q[0];
sx q[0];
rz(-0.43061391) q[0];
sx q[0];
rz(-0.15821247) q[0];
rz(-1.0868866) q[1];
sx q[1];
rz(-2.465076) q[1];
sx q[1];
rz(-2.3349082) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2465625) q[0];
sx q[0];
rz(-0.44996214) q[0];
sx q[0];
rz(-1.297784) q[0];
x q[1];
rz(2.9358747) q[2];
sx q[2];
rz(-1.6763902) q[2];
sx q[2];
rz(0.92668698) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1572691) q[1];
sx q[1];
rz(-1.8586129) q[1];
sx q[1];
rz(-2.4411574) q[1];
x q[2];
rz(1.4679457) q[3];
sx q[3];
rz(-2.8392483) q[3];
sx q[3];
rz(-0.24149382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.57758254) q[2];
sx q[2];
rz(-1.6489886) q[2];
sx q[2];
rz(0.99964833) q[2];
rz(-0.10351652) q[3];
sx q[3];
rz(-3.0045356) q[3];
sx q[3];
rz(2.0172393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0209811) q[0];
sx q[0];
rz(-2.3636901) q[0];
sx q[0];
rz(-0.65752423) q[0];
rz(-2.855037) q[1];
sx q[1];
rz(-2.2189238) q[1];
sx q[1];
rz(0.34067571) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9690902) q[0];
sx q[0];
rz(-1.889125) q[0];
sx q[0];
rz(-2.3031635) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9195243) q[2];
sx q[2];
rz(-2.1246315) q[2];
sx q[2];
rz(1.008322) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.162519) q[1];
sx q[1];
rz(-2.5538429) q[1];
sx q[1];
rz(-1.9025004) q[1];
x q[2];
rz(-1.7394003) q[3];
sx q[3];
rz(-1.3753034) q[3];
sx q[3];
rz(0.80381264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.75491536) q[2];
sx q[2];
rz(-1.9415386) q[2];
sx q[2];
rz(0.49794751) q[2];
rz(2.8399816) q[3];
sx q[3];
rz(-0.40062723) q[3];
sx q[3];
rz(0.51914674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4122445) q[0];
sx q[0];
rz(-0.059818581) q[0];
sx q[0];
rz(-2.8630032) q[0];
rz(-0.57922286) q[1];
sx q[1];
rz(-0.93943739) q[1];
sx q[1];
rz(-3.0648807) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0893223) q[0];
sx q[0];
rz(-2.8025586) q[0];
sx q[0];
rz(-1.9498755) q[0];
rz(-pi) q[1];
rz(2.1439728) q[2];
sx q[2];
rz(-1.4263037) q[2];
sx q[2];
rz(-1.0695374) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0662688) q[1];
sx q[1];
rz(-1.2053718) q[1];
sx q[1];
rz(1.5484023) q[1];
rz(-pi) q[2];
rz(-2.9647787) q[3];
sx q[3];
rz(-1.483695) q[3];
sx q[3];
rz(-2.2967695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.77438337) q[2];
sx q[2];
rz(-2.4343906) q[2];
sx q[2];
rz(-0.45551604) q[2];
rz(0.43141836) q[3];
sx q[3];
rz(-3.016267) q[3];
sx q[3];
rz(-3.07807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9615622) q[0];
sx q[0];
rz(-1.4008235) q[0];
sx q[0];
rz(0.93749198) q[0];
rz(-0.58615276) q[1];
sx q[1];
rz(-1.6393607) q[1];
sx q[1];
rz(1.6929109) q[1];
rz(2.0942681) q[2];
sx q[2];
rz(-2.6101255) q[2];
sx q[2];
rz(2.5892467) q[2];
rz(2.4424845) q[3];
sx q[3];
rz(-0.53634488) q[3];
sx q[3];
rz(2.8041822) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
