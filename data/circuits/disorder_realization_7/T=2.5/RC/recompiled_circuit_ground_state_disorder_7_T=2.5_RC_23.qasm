OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.0291542) q[0];
sx q[0];
rz(-1.6310863) q[0];
sx q[0];
rz(-0.36614585) q[0];
rz(-1.0547628) q[1];
sx q[1];
rz(-1.2843479) q[1];
sx q[1];
rz(-2.4020014) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8407362) q[0];
sx q[0];
rz(-1.2175517) q[0];
sx q[0];
rz(-2.059444) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.97402699) q[2];
sx q[2];
rz(-1.0218909) q[2];
sx q[2];
rz(3.10534) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.89902564) q[1];
sx q[1];
rz(-2.0787313) q[1];
sx q[1];
rz(-2.6655156) q[1];
x q[2];
rz(1.8424019) q[3];
sx q[3];
rz(-1.8823307) q[3];
sx q[3];
rz(2.2304362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4313878) q[2];
sx q[2];
rz(-1.7040161) q[2];
sx q[2];
rz(-1.0817184) q[2];
rz(-1.9234575) q[3];
sx q[3];
rz(-1.5617322) q[3];
sx q[3];
rz(1.6703828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0064938) q[0];
sx q[0];
rz(-2.3068937) q[0];
sx q[0];
rz(0.83719069) q[0];
rz(2.2039425) q[1];
sx q[1];
rz(-1.7559914) q[1];
sx q[1];
rz(-3.0167276) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3628085) q[0];
sx q[0];
rz(-1.1511401) q[0];
sx q[0];
rz(-2.5843689) q[0];
rz(-pi) q[1];
rz(2.5636531) q[2];
sx q[2];
rz(-1.8252531) q[2];
sx q[2];
rz(2.6277666) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8229577) q[1];
sx q[1];
rz(-2.1698916) q[1];
sx q[1];
rz(-1.4692719) q[1];
rz(0.11294194) q[3];
sx q[3];
rz(-2.8340351) q[3];
sx q[3];
rz(1.2483734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.79895926) q[2];
sx q[2];
rz(-2.0755167) q[2];
sx q[2];
rz(-0.40668818) q[2];
rz(-2.0090571) q[3];
sx q[3];
rz(-1.6228024) q[3];
sx q[3];
rz(2.2741545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36658528) q[0];
sx q[0];
rz(-2.5560515) q[0];
sx q[0];
rz(-2.0612234) q[0];
rz(-0.22251546) q[1];
sx q[1];
rz(-0.76970005) q[1];
sx q[1];
rz(-0.57058191) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8711609) q[0];
sx q[0];
rz(-0.50073114) q[0];
sx q[0];
rz(-0.88520925) q[0];
rz(-pi) q[1];
rz(-0.81265575) q[2];
sx q[2];
rz(-0.77746292) q[2];
sx q[2];
rz(-1.846738) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2192206) q[1];
sx q[1];
rz(-0.76792323) q[1];
sx q[1];
rz(1.6025387) q[1];
x q[2];
rz(-1.7526083) q[3];
sx q[3];
rz(-2.1382647) q[3];
sx q[3];
rz(2.770707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3065765) q[2];
sx q[2];
rz(-1.4811652) q[2];
sx q[2];
rz(1.2349077) q[2];
rz(2.0645449) q[3];
sx q[3];
rz(-1.5826694) q[3];
sx q[3];
rz(-2.7243015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1398337) q[0];
sx q[0];
rz(-1.7708906) q[0];
sx q[0];
rz(0.68029252) q[0];
rz(-1.744005) q[1];
sx q[1];
rz(-2.0060507) q[1];
sx q[1];
rz(-2.9063693) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096073791) q[0];
sx q[0];
rz(-2.1102724) q[0];
sx q[0];
rz(1.9941313) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.806947) q[2];
sx q[2];
rz(-0.88289936) q[2];
sx q[2];
rz(-2.2918224) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4724303) q[1];
sx q[1];
rz(-1.002527) q[1];
sx q[1];
rz(-1.8316339) q[1];
rz(1.6793988) q[3];
sx q[3];
rz(-1.3050526) q[3];
sx q[3];
rz(-2.7419529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.816421) q[2];
sx q[2];
rz(-2.0656526) q[2];
sx q[2];
rz(1.761033) q[2];
rz(-2.374968) q[3];
sx q[3];
rz(-2.4216757) q[3];
sx q[3];
rz(-2.5015674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0945053) q[0];
sx q[0];
rz(-2.6386059) q[0];
sx q[0];
rz(-1.267953) q[0];
rz(-0.9371593) q[1];
sx q[1];
rz(-0.93473923) q[1];
sx q[1];
rz(0.25757214) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1510096) q[0];
sx q[0];
rz(-2.5611612) q[0];
sx q[0];
rz(1.849017) q[0];
rz(-pi) q[1];
rz(2.0054162) q[2];
sx q[2];
rz(-2.4939257) q[2];
sx q[2];
rz(-0.90463624) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3993381) q[1];
sx q[1];
rz(-0.4961001) q[1];
sx q[1];
rz(-0.45544405) q[1];
rz(-2.9048728) q[3];
sx q[3];
rz(-2.0437996) q[3];
sx q[3];
rz(1.6110171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3246954) q[2];
sx q[2];
rz(-2.2163053) q[2];
sx q[2];
rz(2.2980105) q[2];
rz(-0.90513611) q[3];
sx q[3];
rz(-0.89591187) q[3];
sx q[3];
rz(0.26346537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0596685) q[0];
sx q[0];
rz(-0.34316871) q[0];
sx q[0];
rz(-1.4362417) q[0];
rz(1.2417271) q[1];
sx q[1];
rz(-1.714434) q[1];
sx q[1];
rz(2.5232975) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4390608) q[0];
sx q[0];
rz(-1.0510604) q[0];
sx q[0];
rz(0.63430877) q[0];
x q[1];
rz(0.16538362) q[2];
sx q[2];
rz(-2.2464387) q[2];
sx q[2];
rz(-2.7242416) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.30696318) q[1];
sx q[1];
rz(-1.8719881) q[1];
sx q[1];
rz(-0.42940112) q[1];
rz(-pi) q[2];
rz(0.49694903) q[3];
sx q[3];
rz(-2.5922311) q[3];
sx q[3];
rz(-0.13298154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5228086) q[2];
sx q[2];
rz(-1.8960543) q[2];
sx q[2];
rz(1.9469384) q[2];
rz(-1.6510125) q[3];
sx q[3];
rz(-1.4158764) q[3];
sx q[3];
rz(1.850261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2456583) q[0];
sx q[0];
rz(-0.36933649) q[0];
sx q[0];
rz(2.9781407) q[0];
rz(-2.5864736) q[1];
sx q[1];
rz(-1.4733543) q[1];
sx q[1];
rz(-2.3177573) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30310589) q[0];
sx q[0];
rz(-0.34922276) q[0];
sx q[0];
rz(0.43189148) q[0];
rz(-pi) q[1];
rz(-1.8386318) q[2];
sx q[2];
rz(-0.80129708) q[2];
sx q[2];
rz(1.6472226) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.40385622) q[1];
sx q[1];
rz(-1.1178167) q[1];
sx q[1];
rz(1.7688874) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7300445) q[3];
sx q[3];
rz(-1.2661066) q[3];
sx q[3];
rz(-2.2652486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0865563) q[2];
sx q[2];
rz(-1.8359416) q[2];
sx q[2];
rz(-1.6769064) q[2];
rz(0.054051789) q[3];
sx q[3];
rz(-2.2650104) q[3];
sx q[3];
rz(-2.7189253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5257877) q[0];
sx q[0];
rz(-0.38830385) q[0];
sx q[0];
rz(-1.7098606) q[0];
rz(1.8852385) q[1];
sx q[1];
rz(-1.4579371) q[1];
sx q[1];
rz(1.2725405) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4542101) q[0];
sx q[0];
rz(-2.2159528) q[0];
sx q[0];
rz(-1.7560455) q[0];
rz(-pi) q[1];
rz(0.48154837) q[2];
sx q[2];
rz(-0.56101834) q[2];
sx q[2];
rz(2.53538) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1450729) q[1];
sx q[1];
rz(-2.6710018) q[1];
sx q[1];
rz(-2.3595292) q[1];
x q[2];
rz(-2.3846517) q[3];
sx q[3];
rz(-1.310756) q[3];
sx q[3];
rz(3.0631331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.58497477) q[2];
sx q[2];
rz(-1.7375526) q[2];
sx q[2];
rz(0.72648826) q[2];
rz(-3.0226184) q[3];
sx q[3];
rz(-1.3930895) q[3];
sx q[3];
rz(-2.6208641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93100905) q[0];
sx q[0];
rz(-1.9887661) q[0];
sx q[0];
rz(-0.3535122) q[0];
rz(1.968169) q[1];
sx q[1];
rz(-1.2856154) q[1];
sx q[1];
rz(-1.2849503) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.374732) q[0];
sx q[0];
rz(-2.6832125) q[0];
sx q[0];
rz(-1.2552257) q[0];
x q[1];
rz(-1.1750269) q[2];
sx q[2];
rz(-1.9916196) q[2];
sx q[2];
rz(1.2683587) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8612449) q[1];
sx q[1];
rz(-1.2767562) q[1];
sx q[1];
rz(3.1190013) q[1];
rz(-pi) q[2];
rz(-1.7400763) q[3];
sx q[3];
rz(-0.96099412) q[3];
sx q[3];
rz(-0.15925254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6080007) q[2];
sx q[2];
rz(-2.3054391) q[2];
sx q[2];
rz(0.81544915) q[2];
rz(0.6472553) q[3];
sx q[3];
rz(-1.7630968) q[3];
sx q[3];
rz(-0.7816202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4746998) q[0];
sx q[0];
rz(-0.26645461) q[0];
sx q[0];
rz(1.4315963) q[0];
rz(-0.41150269) q[1];
sx q[1];
rz(-1.1485547) q[1];
sx q[1];
rz(3.0130951) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14376613) q[0];
sx q[0];
rz(-2.3989005) q[0];
sx q[0];
rz(-0.48074845) q[0];
rz(-2.2811878) q[2];
sx q[2];
rz(-1.6771072) q[2];
sx q[2];
rz(-0.27731976) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4849699) q[1];
sx q[1];
rz(-0.8994714) q[1];
sx q[1];
rz(0.89793082) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.64537489) q[3];
sx q[3];
rz(-0.62964343) q[3];
sx q[3];
rz(-0.30984344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0633885) q[2];
sx q[2];
rz(-1.3271164) q[2];
sx q[2];
rz(2.3627031) q[2];
rz(-1.8552467) q[3];
sx q[3];
rz(-1.0948007) q[3];
sx q[3];
rz(1.8761427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.417199) q[0];
sx q[0];
rz(-1.9211641) q[0];
sx q[0];
rz(-1.7838508) q[0];
rz(-2.4667274) q[1];
sx q[1];
rz(-1.4874896) q[1];
sx q[1];
rz(1.0543324) q[1];
rz(-0.24302287) q[2];
sx q[2];
rz(-2.7293548) q[2];
sx q[2];
rz(-1.7133452) q[2];
rz(-2.3159977) q[3];
sx q[3];
rz(-1.5271689) q[3];
sx q[3];
rz(0.76693514) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
