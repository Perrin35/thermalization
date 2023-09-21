OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.83710837) q[0];
sx q[0];
rz(-1.4533071) q[0];
sx q[0];
rz(0.31153554) q[0];
rz(2.7040634) q[1];
sx q[1];
rz(4.9649927) q[1];
sx q[1];
rz(8.8658219) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6214949) q[0];
sx q[0];
rz(-1.9858452) q[0];
sx q[0];
rz(0.15226224) q[0];
rz(-2.5842701) q[2];
sx q[2];
rz(-1.5814591) q[2];
sx q[2];
rz(-2.9023841) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.641687) q[1];
sx q[1];
rz(-0.99398621) q[1];
sx q[1];
rz(0.66024248) q[1];
x q[2];
rz(1.8564838) q[3];
sx q[3];
rz(-1.1477071) q[3];
sx q[3];
rz(-0.74564122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3922334) q[2];
sx q[2];
rz(-1.8584676) q[2];
sx q[2];
rz(0.63670811) q[2];
rz(0.84896815) q[3];
sx q[3];
rz(-2.520112) q[3];
sx q[3];
rz(2.9076715) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37671509) q[0];
sx q[0];
rz(-0.24704084) q[0];
sx q[0];
rz(2.9887181) q[0];
rz(-0.75694594) q[1];
sx q[1];
rz(-1.5544954) q[1];
sx q[1];
rz(2.1551932) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2043641) q[0];
sx q[0];
rz(-1.5110656) q[0];
sx q[0];
rz(1.6110957) q[0];
rz(-pi) q[1];
rz(1.0436922) q[2];
sx q[2];
rz(-0.88000789) q[2];
sx q[2];
rz(-1.3259128) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1101148) q[1];
sx q[1];
rz(-1.8077393) q[1];
sx q[1];
rz(0.38775878) q[1];
x q[2];
rz(-1.9513449) q[3];
sx q[3];
rz(-0.35223397) q[3];
sx q[3];
rz(0.16570839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5793005) q[2];
sx q[2];
rz(-1.219517) q[2];
sx q[2];
rz(0.78312773) q[2];
rz(0.018571818) q[3];
sx q[3];
rz(-1.637807) q[3];
sx q[3];
rz(2.7338681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0531533) q[0];
sx q[0];
rz(-2.8494371) q[0];
sx q[0];
rz(-2.1799178) q[0];
rz(-2.7812474) q[1];
sx q[1];
rz(-2.0397489) q[1];
sx q[1];
rz(0.12869421) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1291217) q[0];
sx q[0];
rz(-0.095329849) q[0];
sx q[0];
rz(2.096813) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1305069) q[2];
sx q[2];
rz(-2.2957544) q[2];
sx q[2];
rz(0.17979187) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8752746) q[1];
sx q[1];
rz(-0.53578636) q[1];
sx q[1];
rz(2.5712625) q[1];
x q[2];
rz(-3.0069628) q[3];
sx q[3];
rz(-2.3343146) q[3];
sx q[3];
rz(-1.5546297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.06015691) q[2];
sx q[2];
rz(-1.9443941) q[2];
sx q[2];
rz(-1.241768) q[2];
rz(0.5870108) q[3];
sx q[3];
rz(-0.95654064) q[3];
sx q[3];
rz(-2.164042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
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
rz(2.509165) q[0];
sx q[0];
rz(-0.88212633) q[0];
sx q[0];
rz(-1.084491) q[0];
rz(-1.4831316) q[1];
sx q[1];
rz(-2.5741534) q[1];
sx q[1];
rz(3.049057) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1922069) q[0];
sx q[0];
rz(-1.8398251) q[0];
sx q[0];
rz(1.2775248) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7961568) q[2];
sx q[2];
rz(-1.1159117) q[2];
sx q[2];
rz(0.60857426) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.98574084) q[1];
sx q[1];
rz(-2.7936613) q[1];
sx q[1];
rz(2.9702529) q[1];
rz(-pi) q[2];
x q[2];
rz(0.16102287) q[3];
sx q[3];
rz(-0.92217126) q[3];
sx q[3];
rz(0.30740689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.83355054) q[2];
sx q[2];
rz(-1.7001067) q[2];
sx q[2];
rz(-2.8095424) q[2];
rz(-1.0559233) q[3];
sx q[3];
rz(-0.27763405) q[3];
sx q[3];
rz(-2.5312996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8191391) q[0];
sx q[0];
rz(-1.8503014) q[0];
sx q[0];
rz(0.21155587) q[0];
rz(-1.8353204) q[1];
sx q[1];
rz(-1.2439367) q[1];
sx q[1];
rz(-0.64770118) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.081213148) q[0];
sx q[0];
rz(-2.9465284) q[0];
sx q[0];
rz(-2.5293406) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7366016) q[2];
sx q[2];
rz(-0.37618033) q[2];
sx q[2];
rz(2.1567547) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.065121) q[1];
sx q[1];
rz(-1.3741125) q[1];
sx q[1];
rz(1.0573439) q[1];
rz(-pi) q[2];
rz(1.5858438) q[3];
sx q[3];
rz(-1.5516557) q[3];
sx q[3];
rz(-1.0552989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5806879) q[2];
sx q[2];
rz(-0.40955341) q[2];
sx q[2];
rz(-2.4482751) q[2];
rz(2.4723315) q[3];
sx q[3];
rz(-1.4368493) q[3];
sx q[3];
rz(-0.0049237331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82419056) q[0];
sx q[0];
rz(-0.22739246) q[0];
sx q[0];
rz(1.2325226) q[0];
rz(2.0690074) q[1];
sx q[1];
rz(-1.0718081) q[1];
sx q[1];
rz(-0.17428621) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6264544) q[0];
sx q[0];
rz(-2.2720552) q[0];
sx q[0];
rz(-2.7324972) q[0];
rz(0.66138791) q[2];
sx q[2];
rz(-2.1629913) q[2];
sx q[2];
rz(-2.4370898) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4982521) q[1];
sx q[1];
rz(-0.87264112) q[1];
sx q[1];
rz(-0.50486418) q[1];
rz(0.11630451) q[3];
sx q[3];
rz(-1.8810086) q[3];
sx q[3];
rz(0.54799622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3198513) q[2];
sx q[2];
rz(-2.3345626) q[2];
sx q[2];
rz(0.19763395) q[2];
rz(0.28891426) q[3];
sx q[3];
rz(-2.2646326) q[3];
sx q[3];
rz(-1.7355841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.063868) q[0];
sx q[0];
rz(-0.48982319) q[0];
sx q[0];
rz(-2.9329964) q[0];
rz(0.96616191) q[1];
sx q[1];
rz(-1.9816793) q[1];
sx q[1];
rz(-1.5055515) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2710072) q[0];
sx q[0];
rz(-0.85575543) q[0];
sx q[0];
rz(2.4164651) q[0];
rz(-2.7111972) q[2];
sx q[2];
rz(-1.6039404) q[2];
sx q[2];
rz(-0.37552777) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.437285) q[1];
sx q[1];
rz(-2.0211126) q[1];
sx q[1];
rz(2.1561161) q[1];
rz(0.286245) q[3];
sx q[3];
rz(-1.4713333) q[3];
sx q[3];
rz(-0.092982987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2839526) q[2];
sx q[2];
rz(-0.74308926) q[2];
sx q[2];
rz(0.097578438) q[2];
rz(-1.7476667) q[3];
sx q[3];
rz(-1.7912309) q[3];
sx q[3];
rz(-0.71684366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39032787) q[0];
sx q[0];
rz(-1.8126235) q[0];
sx q[0];
rz(2.6275997) q[0];
rz(0.12318525) q[1];
sx q[1];
rz(-0.24736483) q[1];
sx q[1];
rz(0.93200144) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4829464) q[0];
sx q[0];
rz(-2.6771149) q[0];
sx q[0];
rz(-1.2645725) q[0];
x q[1];
rz(-2.2064662) q[2];
sx q[2];
rz(-1.7889651) q[2];
sx q[2];
rz(-0.93014923) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.92703687) q[1];
sx q[1];
rz(-2.3408457) q[1];
sx q[1];
rz(-2.3826249) q[1];
rz(-pi) q[2];
rz(1.7979513) q[3];
sx q[3];
rz(-0.72941581) q[3];
sx q[3];
rz(-1.9279355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1121858) q[2];
sx q[2];
rz(-1.1108578) q[2];
sx q[2];
rz(-0.4294447) q[2];
rz(-1.9321692) q[3];
sx q[3];
rz(-2.7691787) q[3];
sx q[3];
rz(0.69303524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3747303) q[0];
sx q[0];
rz(-1.6748036) q[0];
sx q[0];
rz(1.3388348) q[0];
rz(0.70612899) q[1];
sx q[1];
rz(-1.2495722) q[1];
sx q[1];
rz(-0.95058092) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4336193) q[0];
sx q[0];
rz(-0.74445671) q[0];
sx q[0];
rz(2.3696193) q[0];
rz(-pi) q[1];
rz(-1.4016897) q[2];
sx q[2];
rz(-0.42113129) q[2];
sx q[2];
rz(0.087547628) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.49950019) q[1];
sx q[1];
rz(-1.7576808) q[1];
sx q[1];
rz(1.0287813) q[1];
rz(-1.9825963) q[3];
sx q[3];
rz(-1.8048865) q[3];
sx q[3];
rz(1.7115418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6616228) q[2];
sx q[2];
rz(-1.491549) q[2];
sx q[2];
rz(0.48428145) q[2];
rz(2.2144923) q[3];
sx q[3];
rz(-1.3092594) q[3];
sx q[3];
rz(-0.62121975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1442239) q[0];
sx q[0];
rz(-0.08865083) q[0];
sx q[0];
rz(-2.9123059) q[0];
rz(2.7067822) q[1];
sx q[1];
rz(-1.2229342) q[1];
sx q[1];
rz(0.71892175) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2842429) q[0];
sx q[0];
rz(-1.6650668) q[0];
sx q[0];
rz(2.1988792) q[0];
x q[1];
rz(-1.747379) q[2];
sx q[2];
rz(-0.59497661) q[2];
sx q[2];
rz(2.8240311) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5303516) q[1];
sx q[1];
rz(-1.2019005) q[1];
sx q[1];
rz(1.7289274) q[1];
rz(-pi) q[2];
rz(2.6370254) q[3];
sx q[3];
rz(-2.1310398) q[3];
sx q[3];
rz(-1.954078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7583313) q[2];
sx q[2];
rz(-1.9394082) q[2];
sx q[2];
rz(-2.3948005) q[2];
rz(-0.87219277) q[3];
sx q[3];
rz(-0.82834297) q[3];
sx q[3];
rz(-0.94223589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9185716) q[0];
sx q[0];
rz(-1.3970319) q[0];
sx q[0];
rz(-1.3402517) q[0];
rz(-0.37721286) q[1];
sx q[1];
rz(-1.4706392) q[1];
sx q[1];
rz(-2.6249862) q[1];
rz(-1.2582851) q[2];
sx q[2];
rz(-1.0733114) q[2];
sx q[2];
rz(0.14807362) q[2];
rz(1.4217581) q[3];
sx q[3];
rz(-0.85902135) q[3];
sx q[3];
rz(-1.600941) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
