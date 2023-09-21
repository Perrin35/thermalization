OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3044843) q[0];
sx q[0];
rz(-1.6882856) q[0];
sx q[0];
rz(-0.31153554) q[0];
rz(-0.43752924) q[1];
sx q[1];
rz(-1.8234) q[1];
sx q[1];
rz(0.55895609) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5200978) q[0];
sx q[0];
rz(-1.1557475) q[0];
sx q[0];
rz(2.9893304) q[0];
x q[1];
rz(0.55732255) q[2];
sx q[2];
rz(-1.5814591) q[2];
sx q[2];
rz(0.23920857) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.47145876) q[1];
sx q[1];
rz(-1.0308627) q[1];
sx q[1];
rz(0.88175168) q[1];
rz(-pi) q[2];
rz(0.43879126) q[3];
sx q[3];
rz(-1.3109129) q[3];
sx q[3];
rz(2.4364542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7493593) q[2];
sx q[2];
rz(-1.8584676) q[2];
sx q[2];
rz(0.63670811) q[2];
rz(-0.84896815) q[3];
sx q[3];
rz(-0.62148062) q[3];
sx q[3];
rz(2.9076715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37671509) q[0];
sx q[0];
rz(-0.24704084) q[0];
sx q[0];
rz(0.15287457) q[0];
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
rz(0.61030412) q[0];
sx q[0];
rz(-3.0695519) q[0];
sx q[0];
rz(-2.5487367) q[0];
rz(-pi) q[1];
rz(-2.37843) q[2];
sx q[2];
rz(-1.9689416) q[2];
sx q[2];
rz(2.5415908) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0609329) q[1];
sx q[1];
rz(-2.6903209) q[1];
sx q[1];
rz(-2.5732451) q[1];
x q[2];
rz(-1.1902477) q[3];
sx q[3];
rz(-0.35223397) q[3];
sx q[3];
rz(-0.16570839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5793005) q[2];
sx q[2];
rz(-1.219517) q[2];
sx q[2];
rz(0.78312773) q[2];
rz(3.1230208) q[3];
sx q[3];
rz(-1.5037856) q[3];
sx q[3];
rz(2.7338681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0884393) q[0];
sx q[0];
rz(-2.8494371) q[0];
sx q[0];
rz(2.1799178) q[0];
rz(0.36034521) q[1];
sx q[1];
rz(-1.1018437) q[1];
sx q[1];
rz(3.0128984) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0823682) q[0];
sx q[0];
rz(-1.5229862) q[0];
sx q[0];
rz(1.6533018) q[0];
rz(-2.1305069) q[2];
sx q[2];
rz(-0.84583827) q[2];
sx q[2];
rz(-0.17979187) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8752746) q[1];
sx q[1];
rz(-0.53578636) q[1];
sx q[1];
rz(0.57033013) q[1];
rz(-0.13462984) q[3];
sx q[3];
rz(-2.3343146) q[3];
sx q[3];
rz(-1.5869629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.06015691) q[2];
sx q[2];
rz(-1.9443941) q[2];
sx q[2];
rz(-1.8998247) q[2];
rz(2.5545819) q[3];
sx q[3];
rz(-2.185052) q[3];
sx q[3];
rz(-2.164042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.509165) q[0];
sx q[0];
rz(-0.88212633) q[0];
sx q[0];
rz(2.0571016) q[0];
rz(1.658461) q[1];
sx q[1];
rz(-2.5741534) q[1];
sx q[1];
rz(3.049057) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70150347) q[0];
sx q[0];
rz(-1.8532231) q[0];
sx q[0];
rz(-2.861172) q[0];
rz(1.091435) q[2];
sx q[2];
rz(-1.879868) q[2];
sx q[2];
rz(-2.022559) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7177928) q[1];
sx q[1];
rz(-1.6289627) q[1];
sx q[1];
rz(2.7983626) q[1];
rz(2.2257005) q[3];
sx q[3];
rz(-1.4426784) q[3];
sx q[3];
rz(1.3611984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.83355054) q[2];
sx q[2];
rz(-1.4414859) q[2];
sx q[2];
rz(-2.8095424) q[2];
rz(-2.0856693) q[3];
sx q[3];
rz(-2.8639586) q[3];
sx q[3];
rz(-2.5312996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8191391) q[0];
sx q[0];
rz(-1.2912913) q[0];
sx q[0];
rz(-2.9300368) q[0];
rz(-1.8353204) q[1];
sx q[1];
rz(-1.2439367) q[1];
sx q[1];
rz(2.4938915) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0487329) q[0];
sx q[0];
rz(-1.4591685) q[0];
sx q[0];
rz(2.9812921) q[0];
x q[1];
rz(2.7933502) q[2];
sx q[2];
rz(-1.7160545) q[2];
sx q[2];
rz(0.20656221) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9808637) q[1];
sx q[1];
rz(-2.594922) q[1];
sx q[1];
rz(1.9561808) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5858438) q[3];
sx q[3];
rz(-1.5516557) q[3];
sx q[3];
rz(-1.0552989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.56090474) q[2];
sx q[2];
rz(-2.7320392) q[2];
sx q[2];
rz(0.69331759) q[2];
rz(-0.66926113) q[3];
sx q[3];
rz(-1.4368493) q[3];
sx q[3];
rz(3.1366689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3174021) q[0];
sx q[0];
rz(-0.22739246) q[0];
sx q[0];
rz(1.90907) q[0];
rz(2.0690074) q[1];
sx q[1];
rz(-1.0718081) q[1];
sx q[1];
rz(2.9673064) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.106819) q[0];
sx q[0];
rz(-0.79402059) q[0];
sx q[0];
rz(1.130571) q[0];
x q[1];
rz(-2.3107489) q[2];
sx q[2];
rz(-2.2846662) q[2];
sx q[2];
rz(-1.4884399) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7277158) q[1];
sx q[1];
rz(-1.9503647) q[1];
sx q[1];
rz(-0.80645251) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2586081) q[3];
sx q[3];
rz(-1.6815261) q[3];
sx q[3];
rz(-2.1544416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.8217414) q[2];
sx q[2];
rz(-2.3345626) q[2];
sx q[2];
rz(0.19763395) q[2];
rz(-0.28891426) q[3];
sx q[3];
rz(-0.87696004) q[3];
sx q[3];
rz(1.4060085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0777247) q[0];
sx q[0];
rz(-0.48982319) q[0];
sx q[0];
rz(-0.20859627) q[0];
rz(2.1754307) q[1];
sx q[1];
rz(-1.1599133) q[1];
sx q[1];
rz(1.6360412) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.063232139) q[0];
sx q[0];
rz(-2.1713543) q[0];
sx q[0];
rz(2.2230704) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6072636) q[2];
sx q[2];
rz(-2.0009396) q[2];
sx q[2];
rz(-1.9311116) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.98547658) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4671441) q[3];
sx q[3];
rz(-1.2860057) q[3];
sx q[3];
rz(-1.6929975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2839526) q[2];
sx q[2];
rz(-0.74308926) q[2];
sx q[2];
rz(3.0440142) q[2];
rz(-1.7476667) q[3];
sx q[3];
rz(-1.3503617) q[3];
sx q[3];
rz(0.71684366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39032787) q[0];
sx q[0];
rz(-1.8126235) q[0];
sx q[0];
rz(-0.51399291) q[0];
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
rz(-0.46447771) q[0];
sx q[0];
rz(-1.8770201) q[0];
rz(-0.93512647) q[2];
sx q[2];
rz(-1.7889651) q[2];
sx q[2];
rz(0.93014923) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.276928) q[1];
sx q[1];
rz(-1.0229467) q[1];
sx q[1];
rz(-0.95363708) q[1];
x q[2];
rz(0.85429116) q[3];
sx q[3];
rz(-1.4201418) q[3];
sx q[3];
rz(-2.6137969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0294068) q[2];
sx q[2];
rz(-2.0307348) q[2];
sx q[2];
rz(2.712148) q[2];
rz(-1.2094234) q[3];
sx q[3];
rz(-2.7691787) q[3];
sx q[3];
rz(2.4485574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3747303) q[0];
sx q[0];
rz(-1.466789) q[0];
sx q[0];
rz(-1.8027579) q[0];
rz(2.4354637) q[1];
sx q[1];
rz(-1.2495722) q[1];
sx q[1];
rz(0.95058092) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4336193) q[0];
sx q[0];
rz(-2.3971359) q[0];
sx q[0];
rz(-0.77197335) q[0];
rz(-3.0663475) q[2];
sx q[2];
rz(-1.1560455) q[2];
sx q[2];
rz(-3.0441949) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6420925) q[1];
sx q[1];
rz(-1.3839118) q[1];
sx q[1];
rz(-1.0287813) q[1];
rz(-1.9825963) q[3];
sx q[3];
rz(-1.8048865) q[3];
sx q[3];
rz(-1.4300508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6616228) q[2];
sx q[2];
rz(-1.491549) q[2];
sx q[2];
rz(-0.48428145) q[2];
rz(-0.92710036) q[3];
sx q[3];
rz(-1.8323332) q[3];
sx q[3];
rz(0.62121975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.1442239) q[0];
sx q[0];
rz(-0.08865083) q[0];
sx q[0];
rz(2.9123059) q[0];
rz(0.43481049) q[1];
sx q[1];
rz(-1.9186585) q[1];
sx q[1];
rz(0.71892175) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35480354) q[0];
sx q[0];
rz(-2.1956586) q[0];
sx q[0];
rz(-3.025269) q[0];
x q[1];
rz(0.11833338) q[2];
sx q[2];
rz(-2.1552857) q[2];
sx q[2];
rz(-3.0362533) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9467266) q[1];
sx q[1];
rz(-0.39992878) q[1];
sx q[1];
rz(-2.7547794) q[1];
x q[2];
rz(-0.94902456) q[3];
sx q[3];
rz(-1.9927295) q[3];
sx q[3];
rz(2.4728647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7583313) q[2];
sx q[2];
rz(-1.9394082) q[2];
sx q[2];
rz(-2.3948005) q[2];
rz(2.2693999) q[3];
sx q[3];
rz(-2.3132497) q[3];
sx q[3];
rz(0.94223589) q[3];
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
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.223021) q[0];
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
rz(-0.71729284) q[3];
sx q[3];
rz(-1.6834696) q[3];
sx q[3];
rz(-3.0739741) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];