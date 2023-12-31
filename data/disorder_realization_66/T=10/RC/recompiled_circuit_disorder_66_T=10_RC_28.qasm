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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.011933) q[0];
sx q[0];
rz(-1.7190022) q[0];
sx q[0];
rz(-0.11797842) q[0];
x q[1];
rz(-2.6267251) q[2];
sx q[2];
rz(-1.9549184) q[2];
sx q[2];
rz(-1.9799973) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.37576807) q[1];
sx q[1];
rz(-1.7205015) q[1];
sx q[1];
rz(-1.2234115) q[1];
rz(2.1894737) q[3];
sx q[3];
rz(-1.3631571) q[3];
sx q[3];
rz(-1.729897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.4188529) q[2];
sx q[2];
rz(-2.3143694) q[2];
sx q[2];
rz(-2.375405) q[2];
rz(-0.17151672) q[3];
sx q[3];
rz(-2.4092716) q[3];
sx q[3];
rz(2.5550301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8121174) q[0];
sx q[0];
rz(-2.7024039) q[0];
sx q[0];
rz(3.0549333) q[0];
rz(-2.2791729) q[1];
sx q[1];
rz(-2.5984867) q[1];
sx q[1];
rz(0.08509732) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8893196) q[0];
sx q[0];
rz(-0.62250455) q[0];
sx q[0];
rz(-3.1134393) q[0];
rz(-pi) q[1];
rz(0.9688103) q[2];
sx q[2];
rz(-1.3170871) q[2];
sx q[2];
rz(-1.837455) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9566006) q[1];
sx q[1];
rz(-1.6829832) q[1];
sx q[1];
rz(0.086602028) q[1];
rz(-0.99382932) q[3];
sx q[3];
rz(-0.54518632) q[3];
sx q[3];
rz(-2.546437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.26248419) q[2];
sx q[2];
rz(-1.5510677) q[2];
sx q[2];
rz(-1.6764199) q[2];
rz(-1.9654467) q[3];
sx q[3];
rz(-0.90548235) q[3];
sx q[3];
rz(2.8951077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48433205) q[0];
sx q[0];
rz(-0.13467877) q[0];
sx q[0];
rz(2.235967) q[0];
rz(-2.8088645) q[1];
sx q[1];
rz(-0.85488027) q[1];
sx q[1];
rz(1.7571626) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4577643) q[0];
sx q[0];
rz(-2.0450852) q[0];
sx q[0];
rz(1.477924) q[0];
rz(2.7483814) q[2];
sx q[2];
rz(-1.097743) q[2];
sx q[2];
rz(-2.2676603) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.78450173) q[1];
sx q[1];
rz(-2.38584) q[1];
sx q[1];
rz(-0.032553629) q[1];
rz(-1.0043886) q[3];
sx q[3];
rz(-1.1096138) q[3];
sx q[3];
rz(-3.0571292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.24421346) q[2];
sx q[2];
rz(-2.801557) q[2];
sx q[2];
rz(1.3282233) q[2];
rz(1.7437079) q[3];
sx q[3];
rz(-0.54026794) q[3];
sx q[3];
rz(-3.0298997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0760647) q[0];
sx q[0];
rz(-0.17834839) q[0];
sx q[0];
rz(2.4705825) q[0];
rz(1.3440075) q[1];
sx q[1];
rz(-1.9799415) q[1];
sx q[1];
rz(0.20733325) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3610791) q[0];
sx q[0];
rz(-1.5618389) q[0];
sx q[0];
rz(1.6121959) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.030427293) q[2];
sx q[2];
rz(-1.7446339) q[2];
sx q[2];
rz(-1.1241084) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.85992766) q[1];
sx q[1];
rz(-1.7237701) q[1];
sx q[1];
rz(1.8106736) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.585235) q[3];
sx q[3];
rz(-2.0824021) q[3];
sx q[3];
rz(2.5479941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.31071445) q[2];
sx q[2];
rz(-1.0881311) q[2];
sx q[2];
rz(2.3332398) q[2];
rz(-0.80777848) q[3];
sx q[3];
rz(-0.72179663) q[3];
sx q[3];
rz(-0.19367735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5508674) q[0];
sx q[0];
rz(-2.4243675) q[0];
sx q[0];
rz(0.89548683) q[0];
rz(2.8764309) q[1];
sx q[1];
rz(-0.83509713) q[1];
sx q[1];
rz(0.53363824) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4424608) q[0];
sx q[0];
rz(-1.2958382) q[0];
sx q[0];
rz(1.2752227) q[0];
rz(-pi) q[1];
rz(-0.76782121) q[2];
sx q[2];
rz(-1.4066833) q[2];
sx q[2];
rz(-0.62190157) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.071760885) q[1];
sx q[1];
rz(-1.7904209) q[1];
sx q[1];
rz(-1.7204309) q[1];
rz(-2.93612) q[3];
sx q[3];
rz(-0.95147248) q[3];
sx q[3];
rz(1.6831786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2395893) q[2];
sx q[2];
rz(-1.5633554) q[2];
sx q[2];
rz(-2.5732102) q[2];
rz(-1.9249632) q[3];
sx q[3];
rz(-0.47074461) q[3];
sx q[3];
rz(0.34354982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-2.6333273) q[0];
sx q[0];
rz(-2.2631622) q[0];
sx q[0];
rz(1.1761965) q[0];
rz(1.5856702) q[1];
sx q[1];
rz(-2.1891179) q[1];
sx q[1];
rz(-2.1369381) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5414808) q[0];
sx q[0];
rz(-1.2280557) q[0];
sx q[0];
rz(-2.747885) q[0];
x q[1];
rz(1.4899646) q[2];
sx q[2];
rz(-1.1751375) q[2];
sx q[2];
rz(-2.253502) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.55014729) q[1];
sx q[1];
rz(-1.0284371) q[1];
sx q[1];
rz(1.28246) q[1];
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
rz(-2.1514251) q[2];
sx q[2];
rz(0.44815865) q[2];
rz(-0.29799497) q[3];
sx q[3];
rz(-0.69152504) q[3];
sx q[3];
rz(0.34860778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77666831) q[0];
sx q[0];
rz(-1.8402599) q[0];
sx q[0];
rz(0.65814322) q[0];
rz(1.8572042) q[1];
sx q[1];
rz(-0.45416608) q[1];
sx q[1];
rz(0.67214322) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9400529) q[0];
sx q[0];
rz(-1.4048028) q[0];
sx q[0];
rz(0.0075017651) q[0];
rz(-pi) q[1];
x q[1];
rz(0.42598287) q[2];
sx q[2];
rz(-1.9773215) q[2];
sx q[2];
rz(2.4533518) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0587412) q[1];
sx q[1];
rz(-2.5141659) q[1];
sx q[1];
rz(-2.4861366) q[1];
rz(2.4403205) q[3];
sx q[3];
rz(-1.4605195) q[3];
sx q[3];
rz(-0.14227223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.60823524) q[2];
sx q[2];
rz(-2.50864) q[2];
sx q[2];
rz(2.7427924) q[2];
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
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7010715) q[0];
sx q[0];
rz(-0.43061391) q[0];
sx q[0];
rz(2.9833802) q[0];
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
rz(-0.2465625) q[0];
sx q[0];
rz(-2.6916305) q[0];
sx q[0];
rz(-1.8438086) q[0];
rz(-pi) q[1];
rz(0.47862349) q[2];
sx q[2];
rz(-2.9106986) q[2];
sx q[2];
rz(0.17639562) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.17864171) q[1];
sx q[1];
rz(-0.9045524) q[1];
sx q[1];
rz(1.9402177) q[1];
x q[2];
rz(-3.1095803) q[3];
sx q[3];
rz(-1.8714928) q[3];
sx q[3];
rz(-0.34919448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.57758254) q[2];
sx q[2];
rz(-1.492604) q[2];
sx q[2];
rz(-2.1419443) q[2];
rz(-0.10351652) q[3];
sx q[3];
rz(-0.13705702) q[3];
sx q[3];
rz(1.1243533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12061159) q[0];
sx q[0];
rz(-2.3636901) q[0];
sx q[0];
rz(-2.4840684) q[0];
rz(0.28655562) q[1];
sx q[1];
rz(-2.2189238) q[1];
sx q[1];
rz(-2.8009169) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9690902) q[0];
sx q[0];
rz(-1.2524676) q[0];
sx q[0];
rz(-2.3031635) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1358143) q[2];
sx q[2];
rz(-1.7592332) q[2];
sx q[2];
rz(2.4609158) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.162519) q[1];
sx q[1];
rz(-0.58774978) q[1];
sx q[1];
rz(-1.2390922) q[1];
rz(-pi) q[2];
rz(-1.4021923) q[3];
sx q[3];
rz(-1.3753034) q[3];
sx q[3];
rz(-0.80381264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3866773) q[2];
sx q[2];
rz(-1.9415386) q[2];
sx q[2];
rz(-2.6436451) q[2];
rz(0.30161101) q[3];
sx q[3];
rz(-0.40062723) q[3];
sx q[3];
rz(2.6224459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(0.07671193) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3007606) q[0];
sx q[0];
rz(-1.4474086) q[0];
sx q[0];
rz(-1.2542017) q[0];
x q[1];
rz(0.17148359) q[2];
sx q[2];
rz(-2.137261) q[2];
sx q[2];
rz(-2.7330074) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1288651) q[1];
sx q[1];
rz(-0.36607933) q[1];
sx q[1];
rz(3.0831343) q[1];
rz(-pi) q[2];
rz(-0.17681392) q[3];
sx q[3];
rz(-1.483695) q[3];
sx q[3];
rz(-0.84482312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3672093) q[2];
sx q[2];
rz(-0.70720208) q[2];
sx q[2];
rz(-0.45551604) q[2];
rz(2.7101743) q[3];
sx q[3];
rz(-0.12532561) q[3];
sx q[3];
rz(-3.07807) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1800304) q[0];
sx q[0];
rz(-1.7407692) q[0];
sx q[0];
rz(-2.2041007) q[0];
rz(0.58615276) q[1];
sx q[1];
rz(-1.502232) q[1];
sx q[1];
rz(-1.4486817) q[1];
rz(-2.8557599) q[2];
sx q[2];
rz(-1.1163859) q[2];
sx q[2];
rz(-3.1039539) q[2];
rz(-1.9361817) q[3];
sx q[3];
rz(-1.972651) q[3];
sx q[3];
rz(-2.7046711) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
