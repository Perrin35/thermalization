OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.902154) q[0];
sx q[0];
rz(-1.3286989) q[0];
sx q[0];
rz(-2.9297096) q[0];
rz(1.7243241) q[1];
sx q[1];
rz(-0.53242004) q[1];
sx q[1];
rz(-0.37791696) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6110395) q[0];
sx q[0];
rz(-1.5612649) q[0];
sx q[0];
rz(3.1411489) q[0];
rz(-pi) q[1];
rz(1.9340408) q[2];
sx q[2];
rz(-0.9245199) q[2];
sx q[2];
rz(0.56528795) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4204882) q[1];
sx q[1];
rz(-1.2961565) q[1];
sx q[1];
rz(2.4824449) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3232949) q[3];
sx q[3];
rz(-0.94486559) q[3];
sx q[3];
rz(1.4656386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8523031) q[2];
sx q[2];
rz(-1.0785582) q[2];
sx q[2];
rz(0.079785384) q[2];
rz(-0.96528178) q[3];
sx q[3];
rz(-1.4502757) q[3];
sx q[3];
rz(-2.4108346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(2.6642283) q[0];
sx q[0];
rz(-1.0034765) q[0];
sx q[0];
rz(0.088951237) q[0];
rz(-1.3720007) q[1];
sx q[1];
rz(-1.4274495) q[1];
sx q[1];
rz(1.1044097) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2671794) q[0];
sx q[0];
rz(-2.4738418) q[0];
sx q[0];
rz(-2.7655168) q[0];
x q[1];
rz(-0.49141617) q[2];
sx q[2];
rz(-2.2960536) q[2];
sx q[2];
rz(2.8710136) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8364077) q[1];
sx q[1];
rz(-2.6713604) q[1];
sx q[1];
rz(1.9051308) q[1];
x q[2];
rz(-2.6022909) q[3];
sx q[3];
rz(-1.4852583) q[3];
sx q[3];
rz(0.77682367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8770807) q[2];
sx q[2];
rz(-1.1254213) q[2];
sx q[2];
rz(-1.6748927) q[2];
rz(0.029021164) q[3];
sx q[3];
rz(-1.0163739) q[3];
sx q[3];
rz(2.4513054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90917176) q[0];
sx q[0];
rz(-3.0747774) q[0];
sx q[0];
rz(-2.8636041) q[0];
rz(-1.4311283) q[1];
sx q[1];
rz(-0.93228308) q[1];
sx q[1];
rz(0.40036449) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58996449) q[0];
sx q[0];
rz(-1.1271521) q[0];
sx q[0];
rz(-0.69728627) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8258926) q[2];
sx q[2];
rz(-0.9005024) q[2];
sx q[2];
rz(2.3885661) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6975721) q[1];
sx q[1];
rz(-2.6409147) q[1];
sx q[1];
rz(-0.33659192) q[1];
rz(-0.72680803) q[3];
sx q[3];
rz(-0.44379674) q[3];
sx q[3];
rz(-2.3151195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4572738) q[2];
sx q[2];
rz(-1.6088586) q[2];
sx q[2];
rz(-2.686783) q[2];
rz(-1.8999892) q[3];
sx q[3];
rz(-0.95507115) q[3];
sx q[3];
rz(0.72030592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.95551816) q[0];
sx q[0];
rz(-1.3294514) q[0];
sx q[0];
rz(-3.1410826) q[0];
rz(-2.5406802) q[1];
sx q[1];
rz(-2.3020703) q[1];
sx q[1];
rz(-2.9972163) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6116981) q[0];
sx q[0];
rz(-1.7774044) q[0];
sx q[0];
rz(-2.1312461) q[0];
x q[1];
rz(2.8881489) q[2];
sx q[2];
rz(-1.7161233) q[2];
sx q[2];
rz(-3.1325454) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3444654) q[1];
sx q[1];
rz(-2.4542311) q[1];
sx q[1];
rz(-2.7784154) q[1];
x q[2];
rz(1.6779283) q[3];
sx q[3];
rz(-1.0262353) q[3];
sx q[3];
rz(-0.0060826172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1986177) q[2];
sx q[2];
rz(-1.696442) q[2];
sx q[2];
rz(2.1790738) q[2];
rz(-1.5396384) q[3];
sx q[3];
rz(-1.4055777) q[3];
sx q[3];
rz(0.34077728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-1.9050423) q[0];
sx q[0];
rz(-1.6860697) q[0];
sx q[0];
rz(1.1849674) q[0];
rz(-0.22625893) q[1];
sx q[1];
rz(-0.87892756) q[1];
sx q[1];
rz(2.102898) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5796367) q[0];
sx q[0];
rz(-2.547894) q[0];
sx q[0];
rz(-2.9953792) q[0];
rz(-pi) q[1];
rz(2.2009497) q[2];
sx q[2];
rz(-2.2046208) q[2];
sx q[2];
rz(-1.9535106) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.86335582) q[1];
sx q[1];
rz(-1.0627295) q[1];
sx q[1];
rz(-0.79973508) q[1];
rz(-pi) q[2];
x q[2];
rz(3.076055) q[3];
sx q[3];
rz(-2.0100694) q[3];
sx q[3];
rz(-0.36079839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.431939) q[2];
sx q[2];
rz(-2.4434872) q[2];
sx q[2];
rz(0.31614885) q[2];
rz(1.6759253) q[3];
sx q[3];
rz(-1.1301872) q[3];
sx q[3];
rz(-1.3795615) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50452152) q[0];
sx q[0];
rz(-0.9032473) q[0];
sx q[0];
rz(-2.0380518) q[0];
rz(0.48031131) q[1];
sx q[1];
rz(-2.4791398) q[1];
sx q[1];
rz(-0.59741098) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53086583) q[0];
sx q[0];
rz(-1.7590176) q[0];
sx q[0];
rz(-1.782531) q[0];
rz(1.811932) q[2];
sx q[2];
rz(-2.3015907) q[2];
sx q[2];
rz(0.75474778) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.5525517) q[1];
sx q[1];
rz(-1.3197761) q[1];
sx q[1];
rz(1.4000721) q[1];
rz(-1.8036929) q[3];
sx q[3];
rz(-0.38023708) q[3];
sx q[3];
rz(1.5986795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9395113) q[2];
sx q[2];
rz(-1.5152405) q[2];
sx q[2];
rz(1.0763947) q[2];
rz(1.9988029) q[3];
sx q[3];
rz(-2.3484774) q[3];
sx q[3];
rz(-0.49016652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48653212) q[0];
sx q[0];
rz(-1.158411) q[0];
sx q[0];
rz(0.038473815) q[0];
rz(-3.0768652) q[1];
sx q[1];
rz(-1.3836626) q[1];
sx q[1];
rz(-0.23385349) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.192894) q[0];
sx q[0];
rz(-1.3854376) q[0];
sx q[0];
rz(-0.22484397) q[0];
rz(-1.9373158) q[2];
sx q[2];
rz(-1.2784174) q[2];
sx q[2];
rz(0.93380837) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.95215248) q[1];
sx q[1];
rz(-2.2169211) q[1];
sx q[1];
rz(1.4586071) q[1];
x q[2];
rz(-3.0117118) q[3];
sx q[3];
rz(-1.5135153) q[3];
sx q[3];
rz(0.86806017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6528299) q[2];
sx q[2];
rz(-1.5583928) q[2];
sx q[2];
rz(2.5433507) q[2];
rz(-0.11387842) q[3];
sx q[3];
rz(-1.3860393) q[3];
sx q[3];
rz(-2.2940476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4050196) q[0];
sx q[0];
rz(-0.84252715) q[0];
sx q[0];
rz(0.2463499) q[0];
rz(1.8440638) q[1];
sx q[1];
rz(-1.1187226) q[1];
sx q[1];
rz(0.49682239) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4637359) q[0];
sx q[0];
rz(-1.1564768) q[0];
sx q[0];
rz(-0.68092771) q[0];
rz(2.4364528) q[2];
sx q[2];
rz(-1.6821096) q[2];
sx q[2];
rz(1.706858) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0973952) q[1];
sx q[1];
rz(-2.5977511) q[1];
sx q[1];
rz(-2.0378276) q[1];
rz(-pi) q[2];
rz(-1.0786177) q[3];
sx q[3];
rz(-1.7075305) q[3];
sx q[3];
rz(-2.0355952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7184427) q[2];
sx q[2];
rz(-2.6079874) q[2];
sx q[2];
rz(1.144484) q[2];
rz(1.1540958) q[3];
sx q[3];
rz(-1.4242947) q[3];
sx q[3];
rz(0.19788876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6941187) q[0];
sx q[0];
rz(-2.9070774) q[0];
sx q[0];
rz(0.06614729) q[0];
rz(-1.9937531) q[1];
sx q[1];
rz(-1.7639953) q[1];
sx q[1];
rz(0.60417169) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8291289) q[0];
sx q[0];
rz(-1.7760217) q[0];
sx q[0];
rz(0.28136307) q[0];
rz(-pi) q[1];
rz(0.52671098) q[2];
sx q[2];
rz(-1.4350495) q[2];
sx q[2];
rz(1.2413687) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3499602) q[1];
sx q[1];
rz(-0.59795982) q[1];
sx q[1];
rz(-2.7954742) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.82138942) q[3];
sx q[3];
rz(-0.8474955) q[3];
sx q[3];
rz(-2.5625474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8784647) q[2];
sx q[2];
rz(-2.2915514) q[2];
sx q[2];
rz(-2.7086332) q[2];
rz(-1.912502) q[3];
sx q[3];
rz(-1.2058328) q[3];
sx q[3];
rz(-1.3273201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1935254) q[0];
sx q[0];
rz(-1.075241) q[0];
sx q[0];
rz(-1.2731592) q[0];
rz(-2.676447) q[1];
sx q[1];
rz(-1.3659313) q[1];
sx q[1];
rz(0.21496162) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33782712) q[0];
sx q[0];
rz(-1.8233607) q[0];
sx q[0];
rz(2.4830677) q[0];
rz(-2.4768171) q[2];
sx q[2];
rz(-0.6419581) q[2];
sx q[2];
rz(0.73981111) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7219639) q[1];
sx q[1];
rz(-0.89744324) q[1];
sx q[1];
rz(0.90760214) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.56193476) q[3];
sx q[3];
rz(-1.7896381) q[3];
sx q[3];
rz(-2.9666025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2994069) q[2];
sx q[2];
rz(-0.81373787) q[2];
sx q[2];
rz(-2.1827533) q[2];
rz(1.9793319) q[3];
sx q[3];
rz(-1.4872888) q[3];
sx q[3];
rz(-1.4452665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6017629) q[0];
sx q[0];
rz(-1.5400664) q[0];
sx q[0];
rz(-1.6590317) q[0];
rz(2.4330347) q[1];
sx q[1];
rz(-0.19150145) q[1];
sx q[1];
rz(2.3932744) q[1];
rz(0.49025771) q[2];
sx q[2];
rz(-2.5210862) q[2];
sx q[2];
rz(1.6747337) q[2];
rz(-1.3832573) q[3];
sx q[3];
rz(-2.2839727) q[3];
sx q[3];
rz(-1.4061389) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
