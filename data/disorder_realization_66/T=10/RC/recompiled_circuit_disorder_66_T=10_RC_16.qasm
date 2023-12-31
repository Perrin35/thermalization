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
rz(-0.12227585) q[1];
sx q[1];
rz(-0.08637698) q[1];
sx q[1];
rz(-3.1187305) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1296596) q[0];
sx q[0];
rz(-1.7190022) q[0];
sx q[0];
rz(3.0236142) q[0];
rz(-pi) q[1];
rz(0.68732287) q[2];
sx q[2];
rz(-2.5097373) q[2];
sx q[2];
rz(-0.99422115) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7658246) q[1];
sx q[1];
rz(-1.7205015) q[1];
sx q[1];
rz(1.2234115) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25306563) q[3];
sx q[3];
rz(-2.1742637) q[3];
sx q[3];
rz(-0.30482182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.4188529) q[2];
sx q[2];
rz(-0.8272233) q[2];
sx q[2];
rz(-2.375405) q[2];
rz(-0.17151672) q[3];
sx q[3];
rz(-2.4092716) q[3];
sx q[3];
rz(-0.58656251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3294753) q[0];
sx q[0];
rz(-2.7024039) q[0];
sx q[0];
rz(3.0549333) q[0];
rz(0.86241972) q[1];
sx q[1];
rz(-2.5984867) q[1];
sx q[1];
rz(0.08509732) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8893196) q[0];
sx q[0];
rz(-0.62250455) q[0];
sx q[0];
rz(0.028153367) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1413967) q[2];
sx q[2];
rz(-0.64711249) q[2];
sx q[2];
rz(0.61691689) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3955235) q[1];
sx q[1];
rz(-1.4847401) q[1];
sx q[1];
rz(1.6834016) q[1];
rz(-pi) q[2];
rz(0.99382932) q[3];
sx q[3];
rz(-2.5964063) q[3];
sx q[3];
rz(-2.546437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8791085) q[2];
sx q[2];
rz(-1.5905249) q[2];
sx q[2];
rz(1.4651728) q[2];
rz(-1.9654467) q[3];
sx q[3];
rz(-0.90548235) q[3];
sx q[3];
rz(-0.24648497) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48433205) q[0];
sx q[0];
rz(-0.13467877) q[0];
sx q[0];
rz(0.9056257) q[0];
rz(-2.8088645) q[1];
sx q[1];
rz(-0.85488027) q[1];
sx q[1];
rz(-1.3844301) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15554409) q[0];
sx q[0];
rz(-1.4882003) q[0];
sx q[0];
rz(-2.6655469) q[0];
x q[1];
rz(-1.0648107) q[2];
sx q[2];
rz(-1.2227321) q[2];
sx q[2];
rz(2.6315174) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3570909) q[1];
sx q[1];
rz(-0.75575268) q[1];
sx q[1];
rz(-0.032553629) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1372041) q[3];
sx q[3];
rz(-1.1096138) q[3];
sx q[3];
rz(-3.0571292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8973792) q[2];
sx q[2];
rz(-0.34003568) q[2];
sx q[2];
rz(-1.8133694) q[2];
rz(-1.3978847) q[3];
sx q[3];
rz(-0.54026794) q[3];
sx q[3];
rz(0.11169294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.06552799) q[0];
sx q[0];
rz(-0.17834839) q[0];
sx q[0];
rz(0.67101014) q[0];
rz(1.3440075) q[1];
sx q[1];
rz(-1.1616511) q[1];
sx q[1];
rz(2.9342594) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3509388) q[0];
sx q[0];
rz(-1.5293984) q[0];
sx q[0];
rz(3.1326276) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1111654) q[2];
sx q[2];
rz(-1.7446339) q[2];
sx q[2];
rz(-2.0174842) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.281665) q[1];
sx q[1];
rz(-1.7237701) q[1];
sx q[1];
rz(1.8106736) q[1];
rz(-pi) q[2];
rz(0.025709318) q[3];
sx q[3];
rz(-2.6298012) q[3];
sx q[3];
rz(2.5185086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.31071445) q[2];
sx q[2];
rz(-2.0534616) q[2];
sx q[2];
rz(-2.3332398) q[2];
rz(2.3338142) q[3];
sx q[3];
rz(-2.419796) q[3];
sx q[3];
rz(-2.9479153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5508674) q[0];
sx q[0];
rz(-0.71722513) q[0];
sx q[0];
rz(-0.89548683) q[0];
rz(0.26516178) q[1];
sx q[1];
rz(-2.3064955) q[1];
sx q[1];
rz(0.53363824) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95414872) q[0];
sx q[0];
rz(-1.2866409) q[0];
sx q[0];
rz(-0.28676333) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23405481) q[2];
sx q[2];
rz(-2.3599527) q[2];
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
rz(-0.22202613) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2002033) q[3];
sx q[3];
rz(-1.7377059) q[3];
sx q[3];
rz(0.0080136673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2395893) q[2];
sx q[2];
rz(-1.5633554) q[2];
sx q[2];
rz(2.5732102) q[2];
rz(-1.2166294) q[3];
sx q[3];
rz(-2.670848) q[3];
sx q[3];
rz(0.34354982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6333273) q[0];
sx q[0];
rz(-2.2631622) q[0];
sx q[0];
rz(-1.1761965) q[0];
rz(1.5856702) q[1];
sx q[1];
rz(-0.95247477) q[1];
sx q[1];
rz(2.1369381) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4906824) q[0];
sx q[0];
rz(-2.6255529) q[0];
sx q[0];
rz(-0.74923058) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.651628) q[2];
sx q[2];
rz(-1.1751375) q[2];
sx q[2];
rz(0.8880907) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1130044) q[1];
sx q[1];
rz(-2.5341946) q[1];
sx q[1];
rz(-2.7007156) q[1];
rz(-pi) q[2];
rz(1.5980914) q[3];
sx q[3];
rz(-2.7276162) q[3];
sx q[3];
rz(2.9881145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7375609) q[2];
sx q[2];
rz(-0.99016756) q[2];
sx q[2];
rz(-0.44815865) q[2];
rz(-0.29799497) q[3];
sx q[3];
rz(-2.4500676) q[3];
sx q[3];
rz(-0.34860778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77666831) q[0];
sx q[0];
rz(-1.3013327) q[0];
sx q[0];
rz(0.65814322) q[0];
rz(1.8572042) q[1];
sx q[1];
rz(-0.45416608) q[1];
sx q[1];
rz(0.67214322) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9400529) q[0];
sx q[0];
rz(-1.7367898) q[0];
sx q[0];
rz(-0.0075017651) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7156098) q[2];
sx q[2];
rz(-1.9773215) q[2];
sx q[2];
rz(0.68824088) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3232279) q[1];
sx q[1];
rz(-2.0548901) q[1];
sx q[1];
rz(1.9869884) q[1];
rz(-pi) q[2];
rz(1.7147195) q[3];
sx q[3];
rz(-2.2669499) q[3];
sx q[3];
rz(1.6203984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5333574) q[2];
sx q[2];
rz(-0.63295263) q[2];
sx q[2];
rz(-2.7427924) q[2];
rz(0.82459015) q[3];
sx q[3];
rz(-1.7493068) q[3];
sx q[3];
rz(-2.5430172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7010715) q[0];
sx q[0];
rz(-0.43061391) q[0];
sx q[0];
rz(2.9833802) q[0];
rz(2.054706) q[1];
sx q[1];
rz(-0.67651665) q[1];
sx q[1];
rz(-0.80668443) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2465625) q[0];
sx q[0];
rz(-0.44996214) q[0];
sx q[0];
rz(-1.8438086) q[0];
rz(-pi) q[1];
rz(-1.6786472) q[2];
sx q[2];
rz(-1.3662405) q[2];
sx q[2];
rz(2.4754935) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.17864171) q[1];
sx q[1];
rz(-2.2370403) q[1];
sx q[1];
rz(-1.9402177) q[1];
x q[2];
rz(-1.4679457) q[3];
sx q[3];
rz(-0.30234435) q[3];
sx q[3];
rz(2.9000988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5640101) q[2];
sx q[2];
rz(-1.6489886) q[2];
sx q[2];
rz(-2.1419443) q[2];
rz(3.0380761) q[3];
sx q[3];
rz(-3.0045356) q[3];
sx q[3];
rz(2.0172393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0209811) q[0];
sx q[0];
rz(-0.7779026) q[0];
sx q[0];
rz(2.4840684) q[0];
rz(0.28655562) q[1];
sx q[1];
rz(-0.92266881) q[1];
sx q[1];
rz(2.8009169) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4083456) q[0];
sx q[0];
rz(-0.78662965) q[0];
sx q[0];
rz(1.1128845) q[0];
rz(1.2286439) q[2];
sx q[2];
rz(-0.59235448) q[2];
sx q[2];
rz(-2.5387788) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3714911) q[1];
sx q[1];
rz(-2.1226468) q[1];
sx q[1];
rz(-2.9279207) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.70298985) q[3];
sx q[3];
rz(-2.8841416) q[3];
sx q[3];
rz(-1.6182181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3866773) q[2];
sx q[2];
rz(-1.200054) q[2];
sx q[2];
rz(2.6436451) q[2];
rz(-0.30161101) q[3];
sx q[3];
rz(-0.40062723) q[3];
sx q[3];
rz(-2.6224459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
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
rz(-2.5623698) q[1];
sx q[1];
rz(-0.93943739) q[1];
sx q[1];
rz(-0.07671193) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4519276) q[0];
sx q[0];
rz(-1.2566914) q[0];
sx q[0];
rz(3.0118224) q[0];
rz(1.8329343) q[2];
sx q[2];
rz(-2.5524676) q[2];
sx q[2];
rz(2.4208456) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.50347603) q[1];
sx q[1];
rz(-1.5917115) q[1];
sx q[1];
rz(-0.36550826) q[1];
rz(-pi) q[2];
rz(0.46080188) q[3];
sx q[3];
rz(-0.19690234) q[3];
sx q[3];
rz(-0.27289665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3672093) q[2];
sx q[2];
rz(-0.70720208) q[2];
sx q[2];
rz(-0.45551604) q[2];
rz(-2.7101743) q[3];
sx q[3];
rz(-3.016267) q[3];
sx q[3];
rz(-3.07807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(1.1800304) q[0];
sx q[0];
rz(-1.4008235) q[0];
sx q[0];
rz(0.93749198) q[0];
rz(-2.5554399) q[1];
sx q[1];
rz(-1.502232) q[1];
sx q[1];
rz(-1.4486817) q[1];
rz(2.0942681) q[2];
sx q[2];
rz(-2.6101255) q[2];
sx q[2];
rz(2.5892467) q[2];
rz(-2.7145731) q[3];
sx q[3];
rz(-1.9058766) q[3];
sx q[3];
rz(1.8591892) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
