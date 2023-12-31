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
rz(4.8298782) q[0];
sx q[0];
rz(9.7363135) q[0];
rz(2.7040634) q[1];
sx q[1];
rz(-1.3181926) q[1];
sx q[1];
rz(-0.55895609) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1564724) q[0];
sx q[0];
rz(-0.44056842) q[0];
sx q[0];
rz(-1.23929) q[0];
rz(-pi) q[1];
rz(-0.020157651) q[2];
sx q[2];
rz(-0.55741376) q[2];
sx q[2];
rz(1.7928979) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6701339) q[1];
sx q[1];
rz(-2.11073) q[1];
sx q[1];
rz(-0.88175168) q[1];
rz(-pi) q[2];
rz(0.55922237) q[3];
sx q[3];
rz(-0.50563522) q[3];
sx q[3];
rz(-1.7749744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7493593) q[2];
sx q[2];
rz(-1.8584676) q[2];
sx q[2];
rz(0.63670811) q[2];
rz(2.2926245) q[3];
sx q[3];
rz(-2.520112) q[3];
sx q[3];
rz(0.23392114) q[3];
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
rz(pi/2) q[2];
sx q[2];
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
rz(-1.5870973) q[1];
sx q[1];
rz(0.98639948) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9372285) q[0];
sx q[0];
rz(-1.5110656) q[0];
sx q[0];
rz(1.5304969) q[0];
rz(-0.7631626) q[2];
sx q[2];
rz(-1.9689416) q[2];
sx q[2];
rz(0.60000186) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.03147787) q[1];
sx q[1];
rz(-1.8077393) q[1];
sx q[1];
rz(-0.38775878) q[1];
x q[2];
rz(1.2419224) q[3];
sx q[3];
rz(-1.4423014) q[3];
sx q[3];
rz(-1.0458898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5793005) q[2];
sx q[2];
rz(-1.9220756) q[2];
sx q[2];
rz(2.3584649) q[2];
rz(0.018571818) q[3];
sx q[3];
rz(-1.637807) q[3];
sx q[3];
rz(2.7338681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0884393) q[0];
sx q[0];
rz(-0.29215559) q[0];
sx q[0];
rz(-2.1799178) q[0];
rz(-0.36034521) q[1];
sx q[1];
rz(-2.0397489) q[1];
sx q[1];
rz(3.0128984) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51552396) q[0];
sx q[0];
rz(-1.4883853) q[0];
sx q[0];
rz(-3.0936196) q[0];
x q[1];
rz(2.3340036) q[2];
sx q[2];
rz(-1.1620887) q[2];
sx q[2];
rz(-1.7847716) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.80850959) q[1];
sx q[1];
rz(-1.2915478) q[1];
sx q[1];
rz(-2.6781494) q[1];
rz(-pi) q[2];
rz(-1.7101173) q[3];
sx q[3];
rz(-0.77292597) q[3];
sx q[3];
rz(-1.7804002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0814357) q[2];
sx q[2];
rz(-1.1971985) q[2];
sx q[2];
rz(-1.8998247) q[2];
rz(2.5545819) q[3];
sx q[3];
rz(-0.95654064) q[3];
sx q[3];
rz(-0.97755066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.509165) q[0];
sx q[0];
rz(-2.2594663) q[0];
sx q[0];
rz(-1.084491) q[0];
rz(-1.658461) q[1];
sx q[1];
rz(-0.56743923) q[1];
sx q[1];
rz(3.049057) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10034427) q[0];
sx q[0];
rz(-0.39533246) q[0];
sx q[0];
rz(2.3325217) q[0];
rz(0.96524694) q[2];
sx q[2];
rz(-2.5778228) q[2];
sx q[2];
rz(0.077686003) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.98574084) q[1];
sx q[1];
rz(-0.34793138) q[1];
sx q[1];
rz(-2.9702529) q[1];
x q[2];
rz(-1.362364) q[3];
sx q[3];
rz(-2.4760893) q[3];
sx q[3];
rz(-3.0968551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.83355054) q[2];
sx q[2];
rz(-1.7001067) q[2];
sx q[2];
rz(-0.33205024) q[2];
rz(2.0856693) q[3];
sx q[3];
rz(-2.8639586) q[3];
sx q[3];
rz(2.5312996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8191391) q[0];
sx q[0];
rz(-1.2912913) q[0];
sx q[0];
rz(-2.9300368) q[0];
rz(-1.3062723) q[1];
sx q[1];
rz(-1.2439367) q[1];
sx q[1];
rz(0.64770118) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0487329) q[0];
sx q[0];
rz(-1.6824241) q[0];
sx q[0];
rz(2.9812921) q[0];
rz(1.7251882) q[2];
sx q[2];
rz(-1.2263745) q[2];
sx q[2];
rz(-1.7248578) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9808637) q[1];
sx q[1];
rz(-0.54667066) q[1];
sx q[1];
rz(1.9561808) q[1];
rz(0.019142814) q[3];
sx q[3];
rz(-1.585841) q[3];
sx q[3];
rz(-2.6263833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.56090474) q[2];
sx q[2];
rz(-2.7320392) q[2];
sx q[2];
rz(-2.4482751) q[2];
rz(-0.66926113) q[3];
sx q[3];
rz(-1.7047434) q[3];
sx q[3];
rz(-3.1366689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82419056) q[0];
sx q[0];
rz(-2.9142002) q[0];
sx q[0];
rz(-1.90907) q[0];
rz(-2.0690074) q[1];
sx q[1];
rz(-1.0718081) q[1];
sx q[1];
rz(0.17428621) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5151383) q[0];
sx q[0];
rz(-2.2720552) q[0];
sx q[0];
rz(-0.40909543) q[0];
x q[1];
rz(-0.66138791) q[2];
sx q[2];
rz(-0.97860133) q[2];
sx q[2];
rz(-2.4370898) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4982521) q[1];
sx q[1];
rz(-2.2689515) q[1];
sx q[1];
rz(-2.6367285) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[1];
rz(-0.8217414) q[2];
sx q[2];
rz(-0.80703002) q[2];
sx q[2];
rz(0.19763395) q[2];
rz(0.28891426) q[3];
sx q[3];
rz(-0.87696004) q[3];
sx q[3];
rz(1.7355841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.063868) q[0];
sx q[0];
rz(-0.48982319) q[0];
sx q[0];
rz(2.9329964) q[0];
rz(2.1754307) q[1];
sx q[1];
rz(-1.1599133) q[1];
sx q[1];
rz(1.6360412) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2710072) q[0];
sx q[0];
rz(-0.85575543) q[0];
sx q[0];
rz(-2.4164651) q[0];
x q[1];
rz(-1.534329) q[2];
sx q[2];
rz(-1.140653) q[2];
sx q[2];
rz(1.210481) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.437285) q[1];
sx q[1];
rz(-2.0211126) q[1];
sx q[1];
rz(-2.1561161) q[1];
rz(-1.6744485) q[3];
sx q[3];
rz(-1.855587) q[3];
sx q[3];
rz(-1.4485952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.85764) q[2];
sx q[2];
rz(-0.74308926) q[2];
sx q[2];
rz(-0.097578438) q[2];
rz(1.7476667) q[3];
sx q[3];
rz(-1.7912309) q[3];
sx q[3];
rz(-2.424749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39032787) q[0];
sx q[0];
rz(-1.3289691) q[0];
sx q[0];
rz(-2.6275997) q[0];
rz(3.0184074) q[1];
sx q[1];
rz(-2.8942278) q[1];
sx q[1];
rz(0.93200144) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65864627) q[0];
sx q[0];
rz(-2.6771149) q[0];
sx q[0];
rz(1.2645725) q[0];
x q[1];
rz(-0.93512647) q[2];
sx q[2];
rz(-1.7889651) q[2];
sx q[2];
rz(0.93014923) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8646647) q[1];
sx q[1];
rz(-1.0229467) q[1];
sx q[1];
rz(-2.1879556) q[1];
rz(-pi) q[2];
rz(1.3436414) q[3];
sx q[3];
rz(-2.4121768) q[3];
sx q[3];
rz(-1.9279355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1121858) q[2];
sx q[2];
rz(-2.0307348) q[2];
sx q[2];
rz(-0.4294447) q[2];
rz(-1.9321692) q[3];
sx q[3];
rz(-0.37241396) q[3];
sx q[3];
rz(2.4485574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(0.76686239) q[0];
sx q[0];
rz(-1.6748036) q[0];
sx q[0];
rz(1.8027579) q[0];
rz(2.4354637) q[1];
sx q[1];
rz(-1.2495722) q[1];
sx q[1];
rz(-2.1910117) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7079733) q[0];
sx q[0];
rz(-2.3971359) q[0];
sx q[0];
rz(0.77197335) q[0];
rz(1.1549994) q[2];
sx q[2];
rz(-1.5019413) q[2];
sx q[2];
rz(-1.5037675) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1817158) q[1];
sx q[1];
rz(-1.0392337) q[1];
sx q[1];
rz(-0.2172443) q[1];
x q[2];
rz(1.0334942) q[3];
sx q[3];
rz(-2.6712382) q[3];
sx q[3];
rz(0.62894097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6616228) q[2];
sx q[2];
rz(-1.6500436) q[2];
sx q[2];
rz(-0.48428145) q[2];
rz(2.2144923) q[3];
sx q[3];
rz(-1.3092594) q[3];
sx q[3];
rz(2.5203729) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1442239) q[0];
sx q[0];
rz(-0.08865083) q[0];
sx q[0];
rz(-0.22928672) q[0];
rz(2.7067822) q[1];
sx q[1];
rz(-1.9186585) q[1];
sx q[1];
rz(2.4226709) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2842429) q[0];
sx q[0];
rz(-1.6650668) q[0];
sx q[0];
rz(2.1988792) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3942137) q[2];
sx q[2];
rz(-0.59497661) q[2];
sx q[2];
rz(2.8240311) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.016991888) q[1];
sx q[1];
rz(-1.4233839) q[1];
sx q[1];
rz(0.37313811) q[1];
x q[2];
rz(2.2273916) q[3];
sx q[3];
rz(-0.73540348) q[3];
sx q[3];
rz(-0.38287336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7583313) q[2];
sx q[2];
rz(-1.2021844) q[2];
sx q[2];
rz(-0.74679217) q[2];
rz(0.87219277) q[3];
sx q[3];
rz(-0.82834297) q[3];
sx q[3];
rz(0.94223589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.223021) q[0];
sx q[0];
rz(-1.7445607) q[0];
sx q[0];
rz(1.8013409) q[0];
rz(2.7643798) q[1];
sx q[1];
rz(-1.4706392) q[1];
sx q[1];
rz(-2.6249862) q[1];
rz(0.5151667) q[2];
sx q[2];
rz(-2.5611521) q[2];
sx q[2];
rz(-0.44708154) q[2];
rz(-0.17046455) q[3];
sx q[3];
rz(-2.4170605) q[3];
sx q[3];
rz(-1.3749882) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
