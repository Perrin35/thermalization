OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.11243842) q[0];
sx q[0];
rz(-1.5105063) q[0];
sx q[0];
rz(0.36614585) q[0];
rz(-1.0547628) q[1];
sx q[1];
rz(-1.2843479) q[1];
sx q[1];
rz(-2.4020014) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6897637) q[0];
sx q[0];
rz(-2.0268931) q[0];
sx q[0];
rz(-2.7460237) q[0];
rz(-0.74313626) q[2];
sx q[2];
rz(-2.3541267) q[2];
sx q[2];
rz(-2.1894388) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.89902564) q[1];
sx q[1];
rz(-1.0628614) q[1];
sx q[1];
rz(2.6655156) q[1];
rz(-pi) q[2];
rz(0.32259703) q[3];
sx q[3];
rz(-1.3125714) q[3];
sx q[3];
rz(0.74479529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.71020484) q[2];
sx q[2];
rz(-1.7040161) q[2];
sx q[2];
rz(-1.0817184) q[2];
rz(1.2181351) q[3];
sx q[3];
rz(-1.5617322) q[3];
sx q[3];
rz(1.6703828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1350988) q[0];
sx q[0];
rz(-0.83469892) q[0];
sx q[0];
rz(2.304402) q[0];
rz(-2.2039425) q[1];
sx q[1];
rz(-1.3856013) q[1];
sx q[1];
rz(-3.0167276) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21289028) q[0];
sx q[0];
rz(-0.68395185) q[0];
sx q[0];
rz(-0.70080832) q[0];
x q[1];
rz(-2.5636531) q[2];
sx q[2];
rz(-1.3163396) q[2];
sx q[2];
rz(-0.51382609) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.13990046) q[1];
sx q[1];
rz(-2.5349974) q[1];
sx q[1];
rz(0.14735518) q[1];
rz(-pi) q[2];
rz(3.0286507) q[3];
sx q[3];
rz(-0.30755755) q[3];
sx q[3];
rz(1.2483734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3426334) q[2];
sx q[2];
rz(-1.066076) q[2];
sx q[2];
rz(0.40668818) q[2];
rz(1.1325356) q[3];
sx q[3];
rz(-1.6228024) q[3];
sx q[3];
rz(2.2741545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
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
rz(2.7750074) q[0];
sx q[0];
rz(-0.58554119) q[0];
sx q[0];
rz(-1.0803692) q[0];
rz(-2.9190772) q[1];
sx q[1];
rz(-0.76970005) q[1];
sx q[1];
rz(-2.5710107) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8711609) q[0];
sx q[0];
rz(-2.6408615) q[0];
sx q[0];
rz(-2.2563834) q[0];
x q[1];
rz(-2.1913085) q[2];
sx q[2];
rz(-2.0740905) q[2];
sx q[2];
rz(2.2719943) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1751082) q[1];
sx q[1];
rz(-2.3382332) q[1];
sx q[1];
rz(-3.1109555) q[1];
x q[2];
rz(0.57502745) q[3];
sx q[3];
rz(-1.7238657) q[3];
sx q[3];
rz(1.1014155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3065765) q[2];
sx q[2];
rz(-1.6604275) q[2];
sx q[2];
rz(1.906685) q[2];
rz(-1.0770477) q[3];
sx q[3];
rz(-1.5826694) q[3];
sx q[3];
rz(-2.7243015) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0017589105) q[0];
sx q[0];
rz(-1.3707021) q[0];
sx q[0];
rz(-2.4613001) q[0];
rz(1.744005) q[1];
sx q[1];
rz(-1.1355419) q[1];
sx q[1];
rz(0.23522338) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4394191) q[0];
sx q[0];
rz(-1.2105976) q[0];
sx q[0];
rz(-2.5605306) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1905922) q[2];
sx q[2];
rz(-2.3887164) q[2];
sx q[2];
rz(2.7928758) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1295812) q[1];
sx q[1];
rz(-2.5223603) q[1];
sx q[1];
rz(0.38384743) q[1];
rz(-2.7626183) q[3];
sx q[3];
rz(-2.8550006) q[3];
sx q[3];
rz(0.0061356286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.816421) q[2];
sx q[2];
rz(-2.0656526) q[2];
sx q[2];
rz(-1.761033) q[2];
rz(0.76662463) q[3];
sx q[3];
rz(-0.719917) q[3];
sx q[3];
rz(-0.64002526) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0945053) q[0];
sx q[0];
rz(-2.6386059) q[0];
sx q[0];
rz(1.267953) q[0];
rz(0.9371593) q[1];
sx q[1];
rz(-2.2068534) q[1];
sx q[1];
rz(0.25757214) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1510096) q[0];
sx q[0];
rz(-2.5611612) q[0];
sx q[0];
rz(1.2925757) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1361764) q[2];
sx q[2];
rz(-2.4939257) q[2];
sx q[2];
rz(0.90463624) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.89119205) q[1];
sx q[1];
rz(-2.0124984) q[1];
sx q[1];
rz(1.3370727) q[1];
x q[2];
rz(-0.23671984) q[3];
sx q[3];
rz(-1.0977931) q[3];
sx q[3];
rz(1.6110171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.81689721) q[2];
sx q[2];
rz(-2.2163053) q[2];
sx q[2];
rz(2.2980105) q[2];
rz(2.2364565) q[3];
sx q[3];
rz(-2.2456808) q[3];
sx q[3];
rz(-0.26346537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0819241) q[0];
sx q[0];
rz(-2.7984239) q[0];
sx q[0];
rz(-1.4362417) q[0];
rz(-1.2417271) q[1];
sx q[1];
rz(-1.714434) q[1];
sx q[1];
rz(0.61829511) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4390608) q[0];
sx q[0];
rz(-1.0510604) q[0];
sx q[0];
rz(-0.63430877) q[0];
rz(-2.253153) q[2];
sx q[2];
rz(-1.4419781) q[2];
sx q[2];
rz(1.0494378) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.30696318) q[1];
sx q[1];
rz(-1.8719881) q[1];
sx q[1];
rz(-2.7121915) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6446436) q[3];
sx q[3];
rz(-2.5922311) q[3];
sx q[3];
rz(3.0086111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.61878407) q[2];
sx q[2];
rz(-1.8960543) q[2];
sx q[2];
rz(-1.1946542) q[2];
rz(1.4905802) q[3];
sx q[3];
rz(-1.7257163) q[3];
sx q[3];
rz(-1.850261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8959344) q[0];
sx q[0];
rz(-2.7722562) q[0];
sx q[0];
rz(-2.9781407) q[0];
rz(2.5864736) q[1];
sx q[1];
rz(-1.4733543) q[1];
sx q[1];
rz(2.3177573) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15292955) q[0];
sx q[0];
rz(-1.8867765) q[0];
sx q[0];
rz(1.7220604) q[0];
x q[1];
rz(2.3539435) q[2];
sx q[2];
rz(-1.3795492) q[2];
sx q[2];
rz(-0.2650964) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.1154309) q[1];
sx q[1];
rz(-2.6499601) q[1];
sx q[1];
rz(2.7573654) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.6662382) q[3];
sx q[3];
rz(-2.6347646) q[3];
sx q[3];
rz(-0.092286197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0550363) q[2];
sx q[2];
rz(-1.8359416) q[2];
sx q[2];
rz(1.6769064) q[2];
rz(-3.0875409) q[3];
sx q[3];
rz(-2.2650104) q[3];
sx q[3];
rz(-2.7189253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5257877) q[0];
sx q[0];
rz(-0.38830385) q[0];
sx q[0];
rz(-1.7098606) q[0];
rz(1.2563541) q[1];
sx q[1];
rz(-1.4579371) q[1];
sx q[1];
rz(1.8690522) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1521069) q[0];
sx q[0];
rz(-2.4740334) q[0];
sx q[0];
rz(2.9015673) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6334596) q[2];
sx q[2];
rz(-1.8197803) q[2];
sx q[2];
rz(0.5480042) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1570253) q[1];
sx q[1];
rz(-1.8982982) q[1];
sx q[1];
rz(1.2265601) q[1];
rz(-pi) q[2];
rz(-1.2199336) q[3];
sx q[3];
rz(-0.84515709) q[3];
sx q[3];
rz(1.8875287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.58497477) q[2];
sx q[2];
rz(-1.7375526) q[2];
sx q[2];
rz(0.72648826) q[2];
rz(3.0226184) q[3];
sx q[3];
rz(-1.7485031) q[3];
sx q[3];
rz(0.52072853) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2105836) q[0];
sx q[0];
rz(-1.9887661) q[0];
sx q[0];
rz(-0.3535122) q[0];
rz(1.968169) q[1];
sx q[1];
rz(-1.2856154) q[1];
sx q[1];
rz(1.8566424) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48088257) q[0];
sx q[0];
rz(-1.7085643) q[0];
sx q[0];
rz(2.0093927) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6899723) q[2];
sx q[2];
rz(-1.9303782) q[2];
sx q[2];
rz(3.0082085) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8576926) q[1];
sx q[1];
rz(-1.5924179) q[1];
sx q[1];
rz(1.2766854) q[1];
rz(-pi) q[2];
rz(-0.23663123) q[3];
sx q[3];
rz(-0.62997216) q[3];
sx q[3];
rz(-0.13076846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6080007) q[2];
sx q[2];
rz(-0.83615357) q[2];
sx q[2];
rz(-2.3261435) q[2];
rz(2.4943374) q[3];
sx q[3];
rz(-1.7630968) q[3];
sx q[3];
rz(-2.3599724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4746998) q[0];
sx q[0];
rz(-0.26645461) q[0];
sx q[0];
rz(1.7099963) q[0];
rz(2.73009) q[1];
sx q[1];
rz(-1.1485547) q[1];
sx q[1];
rz(-0.12849753) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9978265) q[0];
sx q[0];
rz(-2.3989005) q[0];
sx q[0];
rz(0.48074845) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7329966) q[2];
sx q[2];
rz(-0.71692978) q[2];
sx q[2];
rz(-1.7253815) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.576927) q[1];
sx q[1];
rz(-0.91178545) q[1];
sx q[1];
rz(2.4763649) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.64537489) q[3];
sx q[3];
rz(-2.5119492) q[3];
sx q[3];
rz(0.30984344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0633885) q[2];
sx q[2];
rz(-1.8144763) q[2];
sx q[2];
rz(2.3627031) q[2];
rz(1.8552467) q[3];
sx q[3];
rz(-2.0467919) q[3];
sx q[3];
rz(-1.26545) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
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
rz(-1.4659526) q[2];
sx q[2];
rz(-1.9702199) q[2];
sx q[2];
rz(1.1639845) q[2];
rz(1.5065083) q[3];
sx q[3];
rz(-0.74623204) q[3];
sx q[3];
rz(-0.85109477) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
