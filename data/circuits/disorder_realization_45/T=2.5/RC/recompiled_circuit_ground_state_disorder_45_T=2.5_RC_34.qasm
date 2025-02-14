OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.328182) q[0];
sx q[0];
rz(-1.0398749) q[0];
sx q[0];
rz(1.9591969) q[0];
rz(0.27529588) q[1];
sx q[1];
rz(-2.1405914) q[1];
sx q[1];
rz(-1.4471952) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4997542) q[0];
sx q[0];
rz(-1.9302881) q[0];
sx q[0];
rz(-0.50292869) q[0];
x q[1];
rz(1.916352) q[2];
sx q[2];
rz(-2.3961903) q[2];
sx q[2];
rz(-0.8236664) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1684678) q[1];
sx q[1];
rz(-2.458312) q[1];
sx q[1];
rz(-0.82078187) q[1];
x q[2];
rz(0.98684825) q[3];
sx q[3];
rz(-0.60923558) q[3];
sx q[3];
rz(1.9690994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.64624661) q[2];
sx q[2];
rz(-0.71141156) q[2];
sx q[2];
rz(-1.6515674) q[2];
rz(-1.9735769) q[3];
sx q[3];
rz(-1.370859) q[3];
sx q[3];
rz(-1.542154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9769576) q[0];
sx q[0];
rz(-0.63666207) q[0];
sx q[0];
rz(2.013999) q[0];
rz(-0.98059869) q[1];
sx q[1];
rz(-2.6041928) q[1];
sx q[1];
rz(1.604039) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.011833103) q[0];
sx q[0];
rz(-0.3922222) q[0];
sx q[0];
rz(-1.5935728) q[0];
rz(-pi) q[1];
rz(1.4784548) q[2];
sx q[2];
rz(-1.5870816) q[2];
sx q[2];
rz(-0.0013601842) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.455273) q[1];
sx q[1];
rz(-1.4497633) q[1];
sx q[1];
rz(-2.1081806) q[1];
x q[2];
rz(0.62894024) q[3];
sx q[3];
rz(-2.6835614) q[3];
sx q[3];
rz(-1.1093793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8424707) q[2];
sx q[2];
rz(-0.8490347) q[2];
sx q[2];
rz(0.15057286) q[2];
rz(1.2493089) q[3];
sx q[3];
rz(-0.95730296) q[3];
sx q[3];
rz(1.5875491) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4318749) q[0];
sx q[0];
rz(-1.7133602) q[0];
sx q[0];
rz(3.0336483) q[0];
rz(-0.81186324) q[1];
sx q[1];
rz(-2.4968562) q[1];
sx q[1];
rz(2.937584) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58751719) q[0];
sx q[0];
rz(-0.63236672) q[0];
sx q[0];
rz(-1.5484346) q[0];
rz(-2.7008867) q[2];
sx q[2];
rz(-1.2262029) q[2];
sx q[2];
rz(-1.454892) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.42878681) q[1];
sx q[1];
rz(-1.6920691) q[1];
sx q[1];
rz(1.0324816) q[1];
rz(-2.7864441) q[3];
sx q[3];
rz(-0.85718507) q[3];
sx q[3];
rz(-2.5220152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.418499) q[2];
sx q[2];
rz(-0.30947026) q[2];
sx q[2];
rz(1.2014368) q[2];
rz(-2.5712683) q[3];
sx q[3];
rz(-2.0968962) q[3];
sx q[3];
rz(0.23049155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.044060556) q[0];
sx q[0];
rz(-2.1756797) q[0];
sx q[0];
rz(-3.0934546) q[0];
rz(2.8037181) q[1];
sx q[1];
rz(-0.66403762) q[1];
sx q[1];
rz(-2.7545676) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52464763) q[0];
sx q[0];
rz(-0.94587406) q[0];
sx q[0];
rz(-2.8760103) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0104644) q[2];
sx q[2];
rz(-0.8789458) q[2];
sx q[2];
rz(0.95702167) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.54285262) q[1];
sx q[1];
rz(-2.4495995) q[1];
sx q[1];
rz(2.180468) q[1];
rz(-pi) q[2];
rz(-1.1293255) q[3];
sx q[3];
rz(-2.8366907) q[3];
sx q[3];
rz(-0.68561844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7016697) q[2];
sx q[2];
rz(-2.5950626) q[2];
sx q[2];
rz(-2.2129391) q[2];
rz(-2.6141911) q[3];
sx q[3];
rz(-1.4832067) q[3];
sx q[3];
rz(1.6209013) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1072024) q[0];
sx q[0];
rz(-0.90029383) q[0];
sx q[0];
rz(2.9541556) q[0];
rz(-1.772359) q[1];
sx q[1];
rz(-2.1010459) q[1];
sx q[1];
rz(1.4184515) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.040203) q[0];
sx q[0];
rz(-1.529415) q[0];
sx q[0];
rz(2.1939254) q[0];
x q[1];
rz(-2.4985053) q[2];
sx q[2];
rz(-1.4355112) q[2];
sx q[2];
rz(-2.7500918) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7781648) q[1];
sx q[1];
rz(-1.3194934) q[1];
sx q[1];
rz(-2.9775724) q[1];
rz(-pi) q[2];
rz(-1.6571331) q[3];
sx q[3];
rz(-0.94924474) q[3];
sx q[3];
rz(0.5542258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9886542) q[2];
sx q[2];
rz(-2.8230437) q[2];
sx q[2];
rz(-0.62206507) q[2];
rz(0.41811824) q[3];
sx q[3];
rz(-1.6303635) q[3];
sx q[3];
rz(0.45463425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.37586) q[0];
sx q[0];
rz(-1.3272165) q[0];
sx q[0];
rz(2.3151929) q[0];
rz(1.7835669) q[1];
sx q[1];
rz(-0.91550255) q[1];
sx q[1];
rz(3.0625536) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2717857) q[0];
sx q[0];
rz(-3.0328864) q[0];
sx q[0];
rz(2.4034017) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0108931) q[2];
sx q[2];
rz(-0.93564763) q[2];
sx q[2];
rz(3.1237313) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.88785401) q[1];
sx q[1];
rz(-1.8716545) q[1];
sx q[1];
rz(-2.1734089) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7594245) q[3];
sx q[3];
rz(-0.49083884) q[3];
sx q[3];
rz(-3.0828305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0283811) q[2];
sx q[2];
rz(-1.5225284) q[2];
sx q[2];
rz(1.9436504) q[2];
rz(2.1146663) q[3];
sx q[3];
rz(-0.17470655) q[3];
sx q[3];
rz(-3.0026156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2748264) q[0];
sx q[0];
rz(-2.3737895) q[0];
sx q[0];
rz(-2.7796032) q[0];
rz(-0.29701862) q[1];
sx q[1];
rz(-1.2535008) q[1];
sx q[1];
rz(2.6716935) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7720761) q[0];
sx q[0];
rz(-2.2945991) q[0];
sx q[0];
rz(0.93627255) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5748575) q[2];
sx q[2];
rz(-0.60535678) q[2];
sx q[2];
rz(2.6108205) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1220144) q[1];
sx q[1];
rz(-1.2570001) q[1];
sx q[1];
rz(3.1336844) q[1];
rz(-pi) q[2];
rz(0.04977241) q[3];
sx q[3];
rz(-2.1920466) q[3];
sx q[3];
rz(1.5240108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.517259) q[2];
sx q[2];
rz(-2.2225311) q[2];
sx q[2];
rz(-2.2993235) q[2];
rz(1.9159348) q[3];
sx q[3];
rz(-1.813348) q[3];
sx q[3];
rz(-1.8153502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4351742) q[0];
sx q[0];
rz(-1.59732) q[0];
sx q[0];
rz(-0.73961863) q[0];
rz(2.836272) q[1];
sx q[1];
rz(-1.0605597) q[1];
sx q[1];
rz(1.7230497) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9288257) q[0];
sx q[0];
rz(-2.1706842) q[0];
sx q[0];
rz(-1.3296574) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9786359) q[2];
sx q[2];
rz(-0.96471404) q[2];
sx q[2];
rz(1.1997171) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.67603079) q[1];
sx q[1];
rz(-1.1793765) q[1];
sx q[1];
rz(-1.2483601) q[1];
rz(-pi) q[2];
rz(1.0064686) q[3];
sx q[3];
rz(-0.25403401) q[3];
sx q[3];
rz(2.0182146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.50489396) q[2];
sx q[2];
rz(-0.86504522) q[2];
sx q[2];
rz(-2.0828341) q[2];
rz(-2.5508358) q[3];
sx q[3];
rz(-0.96610132) q[3];
sx q[3];
rz(3.0920715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.010794086) q[0];
sx q[0];
rz(-2.8191691) q[0];
sx q[0];
rz(-1.4149433) q[0];
rz(2.2825799) q[1];
sx q[1];
rz(-2.450727) q[1];
sx q[1];
rz(1.0038928) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5952565) q[0];
sx q[0];
rz(-1.7352859) q[0];
sx q[0];
rz(-0.25877742) q[0];
rz(2.48433) q[2];
sx q[2];
rz(-1.3018152) q[2];
sx q[2];
rz(1.2439234) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.8050766) q[1];
sx q[1];
rz(-2.4081701) q[1];
sx q[1];
rz(-3.0086711) q[1];
x q[2];
rz(-0.74476242) q[3];
sx q[3];
rz(-1.6339374) q[3];
sx q[3];
rz(-2.5936332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.88177195) q[2];
sx q[2];
rz(-1.5544145) q[2];
sx q[2];
rz(0.47015321) q[2];
rz(-1.638089) q[3];
sx q[3];
rz(-0.94109002) q[3];
sx q[3];
rz(-0.21465429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44899061) q[0];
sx q[0];
rz(-0.34994352) q[0];
sx q[0];
rz(2.129648) q[0];
rz(2.8590406) q[1];
sx q[1];
rz(-0.91767445) q[1];
sx q[1];
rz(0.73175398) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5846281) q[0];
sx q[0];
rz(-0.04806732) q[0];
sx q[0];
rz(0.2063423) q[0];
rz(-2.9596161) q[2];
sx q[2];
rz(-1.4018347) q[2];
sx q[2];
rz(0.5398869) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0132799) q[1];
sx q[1];
rz(-1.9802112) q[1];
sx q[1];
rz(-2.4854543) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8153207) q[3];
sx q[3];
rz(-2.7774307) q[3];
sx q[3];
rz(-0.18362602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4628576) q[2];
sx q[2];
rz(-1.5588458) q[2];
sx q[2];
rz(2.8850214) q[2];
rz(-0.76777846) q[3];
sx q[3];
rz(-2.0753746) q[3];
sx q[3];
rz(1.0919051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0302122) q[0];
sx q[0];
rz(-0.82652265) q[0];
sx q[0];
rz(-0.65435456) q[0];
rz(2.4334473) q[1];
sx q[1];
rz(-1.000052) q[1];
sx q[1];
rz(2.4194385) q[1];
rz(-2.3329493) q[2];
sx q[2];
rz(-1.8075347) q[2];
sx q[2];
rz(-1.6703477) q[2];
rz(-3.1307316) q[3];
sx q[3];
rz(-1.167611) q[3];
sx q[3];
rz(2.1394185) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
