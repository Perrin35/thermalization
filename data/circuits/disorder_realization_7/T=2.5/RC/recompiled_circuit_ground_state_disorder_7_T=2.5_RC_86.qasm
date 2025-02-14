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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30692431) q[0];
sx q[0];
rz(-0.59446228) q[0];
sx q[0];
rz(0.90499808) q[0];
rz(0.6366836) q[2];
sx q[2];
rz(-1.0708059) q[2];
sx q[2];
rz(1.9477109) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.084433638) q[1];
sx q[1];
rz(-2.4600303) q[1];
sx q[1];
rz(-2.2595899) q[1];
rz(-pi) q[2];
rz(0.69460709) q[3];
sx q[3];
rz(-0.41037729) q[3];
sx q[3];
rz(-0.17363901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4313878) q[2];
sx q[2];
rz(-1.7040161) q[2];
sx q[2];
rz(-1.0817184) q[2];
rz(-1.9234575) q[3];
sx q[3];
rz(-1.5798605) q[3];
sx q[3];
rz(1.4712099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0064938) q[0];
sx q[0];
rz(-2.3068937) q[0];
sx q[0];
rz(-2.304402) q[0];
rz(2.2039425) q[1];
sx q[1];
rz(-1.3856013) q[1];
sx q[1];
rz(3.0167276) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0406348) q[0];
sx q[0];
rz(-2.0748108) q[0];
sx q[0];
rz(-2.0547778) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5636531) q[2];
sx q[2];
rz(-1.3163396) q[2];
sx q[2];
rz(-0.51382609) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.309545) q[1];
sx q[1];
rz(-1.654594) q[1];
sx q[1];
rz(2.540091) q[1];
x q[2];
rz(1.6065793) q[3];
sx q[3];
rz(-1.2652618) q[3];
sx q[3];
rz(1.3668253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.79895926) q[2];
sx q[2];
rz(-2.0755167) q[2];
sx q[2];
rz(-2.7349045) q[2];
rz(-2.0090571) q[3];
sx q[3];
rz(-1.5187902) q[3];
sx q[3];
rz(0.86743814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7750074) q[0];
sx q[0];
rz(-2.5560515) q[0];
sx q[0];
rz(-1.0803692) q[0];
rz(0.22251546) q[1];
sx q[1];
rz(-0.76970005) q[1];
sx q[1];
rz(-2.5710107) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0208698) q[0];
sx q[0];
rz(-1.1900702) q[0];
sx q[0];
rz(-0.33353593) q[0];
x q[1];
rz(0.81265575) q[2];
sx q[2];
rz(-0.77746292) q[2];
sx q[2];
rz(-1.2948546) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5160061) q[1];
sx q[1];
rz(-1.592844) q[1];
sx q[1];
rz(-0.80312487) q[1];
rz(-pi) q[2];
rz(-0.27640859) q[3];
sx q[3];
rz(-0.59282527) q[3];
sx q[3];
rz(-0.70044493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.83501619) q[2];
sx q[2];
rz(-1.6604275) q[2];
sx q[2];
rz(-1.2349077) q[2];
rz(-1.0770477) q[3];
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
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0017589105) q[0];
sx q[0];
rz(-1.3707021) q[0];
sx q[0];
rz(2.4613001) q[0];
rz(1.744005) q[1];
sx q[1];
rz(-2.0060507) q[1];
sx q[1];
rz(-0.23522338) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0455189) q[0];
sx q[0];
rz(-2.1102724) q[0];
sx q[0];
rz(-1.1474613) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1905922) q[2];
sx q[2];
rz(-0.75287627) q[2];
sx q[2];
rz(2.7928758) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1295812) q[1];
sx q[1];
rz(-2.5223603) q[1];
sx q[1];
rz(-0.38384743) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4621939) q[3];
sx q[3];
rz(-1.83654) q[3];
sx q[3];
rz(0.39963978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.32517165) q[2];
sx q[2];
rz(-1.07594) q[2];
sx q[2];
rz(-1.761033) q[2];
rz(-2.374968) q[3];
sx q[3];
rz(-2.4216757) q[3];
sx q[3];
rz(-2.5015674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0470873) q[0];
sx q[0];
rz(-0.50298679) q[0];
sx q[0];
rz(1.8736396) q[0];
rz(0.9371593) q[1];
sx q[1];
rz(-2.2068534) q[1];
sx q[1];
rz(0.25757214) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66143888) q[0];
sx q[0];
rz(-2.1262125) q[0];
sx q[0];
rz(-0.17819779) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1721971) q[2];
sx q[2];
rz(-1.8276518) q[2];
sx q[2];
rz(1.0207301) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.563411) q[1];
sx q[1];
rz(-1.7817307) q[1];
sx q[1];
rz(2.6891449) q[1];
x q[2];
rz(-2.0553642) q[3];
sx q[3];
rz(-1.3604829) q[3];
sx q[3];
rz(2.9919101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3246954) q[2];
sx q[2];
rz(-0.92528737) q[2];
sx q[2];
rz(0.84358215) q[2];
rz(-2.2364565) q[3];
sx q[3];
rz(-2.2456808) q[3];
sx q[3];
rz(-2.8781273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0819241) q[0];
sx q[0];
rz(-0.34316871) q[0];
sx q[0];
rz(1.7053509) q[0];
rz(-1.8998655) q[1];
sx q[1];
rz(-1.4271586) q[1];
sx q[1];
rz(-2.5232975) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.623659) q[0];
sx q[0];
rz(-2.111064) q[0];
sx q[0];
rz(0.95312693) q[0];
rz(-pi) q[1];
rz(-2.253153) q[2];
sx q[2];
rz(-1.6996146) q[2];
sx q[2];
rz(-1.0494378) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1288207) q[1];
sx q[1];
rz(-1.1619131) q[1];
sx q[1];
rz(-1.9000221) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6478751) q[3];
sx q[3];
rz(-1.3192216) q[3];
sx q[3];
rz(1.2705402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.61878407) q[2];
sx q[2];
rz(-1.2455384) q[2];
sx q[2];
rz(-1.9469384) q[2];
rz(1.6510125) q[3];
sx q[3];
rz(-1.4158764) q[3];
sx q[3];
rz(-1.850261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(1.2456583) q[0];
sx q[0];
rz(-0.36933649) q[0];
sx q[0];
rz(0.16345197) q[0];
rz(2.5864736) q[1];
sx q[1];
rz(-1.4733543) q[1];
sx q[1];
rz(-0.82383531) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4651983) q[0];
sx q[0];
rz(-1.7145183) q[0];
sx q[0];
rz(2.8222047) q[0];
rz(-pi) q[1];
rz(0.26668875) q[2];
sx q[2];
rz(-0.80563918) q[2];
sx q[2];
rz(-1.1186816) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.026161748) q[1];
sx q[1];
rz(-0.49163252) q[1];
sx q[1];
rz(-2.7573654) q[1];
rz(-pi) q[2];
rz(1.9013405) q[3];
sx q[3];
rz(-1.9623266) q[3];
sx q[3];
rz(2.3169404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0865563) q[2];
sx q[2];
rz(-1.3056511) q[2];
sx q[2];
rz(-1.6769064) q[2];
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
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61580491) q[0];
sx q[0];
rz(-0.38830385) q[0];
sx q[0];
rz(-1.4317321) q[0];
rz(-1.2563541) q[1];
sx q[1];
rz(-1.6836555) q[1];
sx q[1];
rz(1.8690522) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4542101) q[0];
sx q[0];
rz(-2.2159528) q[0];
sx q[0];
rz(-1.7560455) q[0];
x q[1];
rz(-0.48154837) q[2];
sx q[2];
rz(-0.56101834) q[2];
sx q[2];
rz(-2.53538) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9965197) q[1];
sx q[1];
rz(-2.6710018) q[1];
sx q[1];
rz(2.3595292) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.75694098) q[3];
sx q[3];
rz(-1.8308367) q[3];
sx q[3];
rz(3.0631331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5566179) q[2];
sx q[2];
rz(-1.7375526) q[2];
sx q[2];
rz(-0.72648826) q[2];
rz(0.11897421) q[3];
sx q[3];
rz(-1.3930895) q[3];
sx q[3];
rz(0.52072853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2105836) q[0];
sx q[0];
rz(-1.1528265) q[0];
sx q[0];
rz(-0.3535122) q[0];
rz(-1.968169) q[1];
sx q[1];
rz(-1.8559772) q[1];
sx q[1];
rz(-1.2849503) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6607101) q[0];
sx q[0];
rz(-1.7085643) q[0];
sx q[0];
rz(-1.1321999) q[0];
rz(-pi) q[1];
rz(0.71106203) q[2];
sx q[2];
rz(-0.56945062) q[2];
sx q[2];
rz(-1.0765778) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8576926) q[1];
sx q[1];
rz(-1.5924179) q[1];
sx q[1];
rz(-1.2766854) q[1];
rz(-pi) q[2];
x q[2];
rz(0.23663123) q[3];
sx q[3];
rz(-0.62997216) q[3];
sx q[3];
rz(0.13076846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6080007) q[2];
sx q[2];
rz(-0.83615357) q[2];
sx q[2];
rz(-2.3261435) q[2];
rz(-2.4943374) q[3];
sx q[3];
rz(-1.7630968) q[3];
sx q[3];
rz(-0.7816202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-1.6668929) q[0];
sx q[0];
rz(-2.875138) q[0];
sx q[0];
rz(1.7099963) q[0];
rz(-0.41150269) q[1];
sx q[1];
rz(-1.1485547) q[1];
sx q[1];
rz(-0.12849753) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.381739) q[0];
sx q[0];
rz(-2.2138192) q[0];
sx q[0];
rz(-1.9722776) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.86040489) q[2];
sx q[2];
rz(-1.4644854) q[2];
sx q[2];
rz(2.8642729) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5646657) q[1];
sx q[1];
rz(-0.91178545) q[1];
sx q[1];
rz(2.4763649) q[1];
rz(-pi) q[2];
rz(-2.4962178) q[3];
sx q[3];
rz(-2.5119492) q[3];
sx q[3];
rz(-0.30984344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0633885) q[2];
sx q[2];
rz(-1.3271164) q[2];
sx q[2];
rz(-2.3627031) q[2];
rz(1.286346) q[3];
sx q[3];
rz(-2.0467919) q[3];
sx q[3];
rz(-1.8761427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7243937) q[0];
sx q[0];
rz(-1.9211641) q[0];
sx q[0];
rz(-1.7838508) q[0];
rz(-2.4667274) q[1];
sx q[1];
rz(-1.4874896) q[1];
sx q[1];
rz(1.0543324) q[1];
rz(1.67564) q[2];
sx q[2];
rz(-1.9702199) q[2];
sx q[2];
rz(1.1639845) q[2];
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
