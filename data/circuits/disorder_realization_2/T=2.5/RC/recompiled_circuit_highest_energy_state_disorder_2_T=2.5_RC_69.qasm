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
rz(2.9444115) q[0];
sx q[0];
rz(-1.0721167) q[0];
sx q[0];
rz(1.5443235) q[0];
rz(0.3577258) q[1];
sx q[1];
rz(4.975714) q[1];
sx q[1];
rz(9.8006021) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25426769) q[0];
sx q[0];
rz(-1.4987117) q[0];
sx q[0];
rz(2.3743126) q[0];
x q[1];
rz(2.69015) q[2];
sx q[2];
rz(-2.5677748) q[2];
sx q[2];
rz(0.30893477) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5615422) q[1];
sx q[1];
rz(-1.0197749) q[1];
sx q[1];
rz(-0.8755133) q[1];
rz(-pi) q[2];
x q[2];
rz(2.014767) q[3];
sx q[3];
rz(-2.4587963) q[3];
sx q[3];
rz(-1.7449774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4481005) q[2];
sx q[2];
rz(-2.7585612) q[2];
sx q[2];
rz(2.7310272) q[2];
rz(0.52452123) q[3];
sx q[3];
rz(-0.21206847) q[3];
sx q[3];
rz(-1.9922403) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9583424) q[0];
sx q[0];
rz(-2.9999833) q[0];
sx q[0];
rz(-0.70567065) q[0];
rz(-1.8535829) q[1];
sx q[1];
rz(-2.5649004) q[1];
sx q[1];
rz(-2.0075331) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41578226) q[0];
sx q[0];
rz(-2.232612) q[0];
sx q[0];
rz(0.51915581) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6666849) q[2];
sx q[2];
rz(-2.7596052) q[2];
sx q[2];
rz(-1.2820006) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4422851) q[1];
sx q[1];
rz(-2.2622427) q[1];
sx q[1];
rz(2.0705018) q[1];
rz(-pi) q[2];
x q[2];
rz(2.864704) q[3];
sx q[3];
rz(-2.5586506) q[3];
sx q[3];
rz(2.049169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.80295339) q[2];
sx q[2];
rz(-0.45265517) q[2];
sx q[2];
rz(-0.53042859) q[2];
rz(-0.28702921) q[3];
sx q[3];
rz(-2.65286) q[3];
sx q[3];
rz(0.73634017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1257783) q[0];
sx q[0];
rz(-1.1017355) q[0];
sx q[0];
rz(1.0544448) q[0];
rz(-1.4796481) q[1];
sx q[1];
rz(-2.2190861) q[1];
sx q[1];
rz(-2.0747917) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9196292) q[0];
sx q[0];
rz(-2.9502075) q[0];
sx q[0];
rz(1.3758977) q[0];
rz(-pi) q[1];
rz(-2.4652723) q[2];
sx q[2];
rz(-1.498885) q[2];
sx q[2];
rz(-1.7184217) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.40854) q[1];
sx q[1];
rz(-2.618737) q[1];
sx q[1];
rz(1.1001292) q[1];
x q[2];
rz(-2.2852632) q[3];
sx q[3];
rz(-1.7995872) q[3];
sx q[3];
rz(-1.8280054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.58831424) q[2];
sx q[2];
rz(-1.30043) q[2];
sx q[2];
rz(2.039457) q[2];
rz(-0.45840248) q[3];
sx q[3];
rz(-1.6130684) q[3];
sx q[3];
rz(2.5100191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.691953) q[0];
sx q[0];
rz(-0.75943685) q[0];
sx q[0];
rz(2.6547883) q[0];
rz(1.5454769) q[1];
sx q[1];
rz(-0.38914746) q[1];
sx q[1];
rz(-2.5130689) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2038949) q[0];
sx q[0];
rz(-1.9638229) q[0];
sx q[0];
rz(-2.9980396) q[0];
rz(-pi) q[1];
rz(-2.1439379) q[2];
sx q[2];
rz(-1.6948683) q[2];
sx q[2];
rz(2.5981552) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0027517) q[1];
sx q[1];
rz(-1.3831487) q[1];
sx q[1];
rz(-1.9676588) q[1];
x q[2];
rz(-2.1194589) q[3];
sx q[3];
rz(-1.5096153) q[3];
sx q[3];
rz(-3.1312628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.57807606) q[2];
sx q[2];
rz(-3.0264644) q[2];
sx q[2];
rz(0.76187491) q[2];
rz(0.21841194) q[3];
sx q[3];
rz(-2.8438447) q[3];
sx q[3];
rz(1.0485605) q[3];
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
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8258284) q[0];
sx q[0];
rz(-2.789412) q[0];
sx q[0];
rz(-2.5370362) q[0];
rz(1.293659) q[1];
sx q[1];
rz(-0.81984055) q[1];
sx q[1];
rz(-0.34436071) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47681706) q[0];
sx q[0];
rz(-1.4882384) q[0];
sx q[0];
rz(-0.020391012) q[0];
rz(-pi) q[1];
rz(-1.9042964) q[2];
sx q[2];
rz(-2.1715559) q[2];
sx q[2];
rz(-1.5432026) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8230668) q[1];
sx q[1];
rz(-2.0289377) q[1];
sx q[1];
rz(-1.0847102) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8959534) q[3];
sx q[3];
rz(-1.1955441) q[3];
sx q[3];
rz(0.2764052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1206426) q[2];
sx q[2];
rz(-2.5072493) q[2];
sx q[2];
rz(0.20993385) q[2];
rz(2.0338992) q[3];
sx q[3];
rz(-1.6898797) q[3];
sx q[3];
rz(-2.8797188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8629859) q[0];
sx q[0];
rz(-2.5513273) q[0];
sx q[0];
rz(-0.069409542) q[0];
rz(0.033992652) q[1];
sx q[1];
rz(-2.4885204) q[1];
sx q[1];
rz(-0.45190826) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0586813) q[0];
sx q[0];
rz(-1.4103019) q[0];
sx q[0];
rz(-2.7468264) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4310125) q[2];
sx q[2];
rz(-2.3167774) q[2];
sx q[2];
rz(1.1912861) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.76594162) q[1];
sx q[1];
rz(-0.022863511) q[1];
sx q[1];
rz(-1.2506865) q[1];
rz(2.8699401) q[3];
sx q[3];
rz(-1.8166887) q[3];
sx q[3];
rz(-0.065448381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7273093) q[2];
sx q[2];
rz(-2.3392129) q[2];
sx q[2];
rz(-1.1451954) q[2];
rz(2.1833873) q[3];
sx q[3];
rz(-1.9327791) q[3];
sx q[3];
rz(2.8632979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3713643) q[0];
sx q[0];
rz(-2.6537277) q[0];
sx q[0];
rz(-2.2801953) q[0];
rz(-0.96249145) q[1];
sx q[1];
rz(-2.1884514) q[1];
sx q[1];
rz(0.15693396) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79377896) q[0];
sx q[0];
rz(-0.89317489) q[0];
sx q[0];
rz(-0.47394808) q[0];
x q[1];
rz(-1.9279582) q[2];
sx q[2];
rz(-1.4712804) q[2];
sx q[2];
rz(-0.9377756) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7879856) q[1];
sx q[1];
rz(-0.79670409) q[1];
sx q[1];
rz(0.70938797) q[1];
rz(-pi) q[2];
rz(-1.3488999) q[3];
sx q[3];
rz(-0.21824868) q[3];
sx q[3];
rz(1.9980007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.12338766) q[2];
sx q[2];
rz(-2.3473479) q[2];
sx q[2];
rz(1.0495074) q[2];
rz(2.6617995) q[3];
sx q[3];
rz(-0.19194651) q[3];
sx q[3];
rz(2.979579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70742339) q[0];
sx q[0];
rz(-2.8987085) q[0];
sx q[0];
rz(-0.49651399) q[0];
rz(-1.6560582) q[1];
sx q[1];
rz(-0.86692202) q[1];
sx q[1];
rz(-0.022580126) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.966205) q[0];
sx q[0];
rz(-1.7838571) q[0];
sx q[0];
rz(2.9609462) q[0];
x q[1];
rz(1.1316001) q[2];
sx q[2];
rz(-2.9832057) q[2];
sx q[2];
rz(-1.8902376) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.39358006) q[1];
sx q[1];
rz(-2.339254) q[1];
sx q[1];
rz(0.2549766) q[1];
rz(-pi) q[2];
x q[2];
rz(0.19666861) q[3];
sx q[3];
rz(-1.0447352) q[3];
sx q[3];
rz(-0.23516178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5217343) q[2];
sx q[2];
rz(-1.1405742) q[2];
sx q[2];
rz(-3.1263331) q[2];
rz(-0.36733019) q[3];
sx q[3];
rz(-2.6945249) q[3];
sx q[3];
rz(2.7801311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58886445) q[0];
sx q[0];
rz(-0.22607729) q[0];
sx q[0];
rz(3.0638301) q[0];
rz(0.11847682) q[1];
sx q[1];
rz(-0.34093726) q[1];
sx q[1];
rz(-0.65611398) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8828147) q[0];
sx q[0];
rz(-1.6973295) q[0];
sx q[0];
rz(1.2650498) q[0];
x q[1];
rz(-2.9703548) q[2];
sx q[2];
rz(-1.3067553) q[2];
sx q[2];
rz(1.1080826) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1826577) q[1];
sx q[1];
rz(-1.1706585) q[1];
sx q[1];
rz(-2.0939287) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5893297) q[3];
sx q[3];
rz(-1.1202462) q[3];
sx q[3];
rz(-1.0737863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.52350998) q[2];
sx q[2];
rz(-1.221035) q[2];
sx q[2];
rz(-2.6005884) q[2];
rz(-2.6909761) q[3];
sx q[3];
rz(-1.485598) q[3];
sx q[3];
rz(-0.97644067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4377874) q[0];
sx q[0];
rz(-1.3492275) q[0];
sx q[0];
rz(-0.10677862) q[0];
rz(-1.5497426) q[1];
sx q[1];
rz(-0.50250643) q[1];
sx q[1];
rz(-1.7639311) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14895591) q[0];
sx q[0];
rz(-2.0201448) q[0];
sx q[0];
rz(-0.0489671) q[0];
x q[1];
rz(2.8286727) q[2];
sx q[2];
rz(-1.7831137) q[2];
sx q[2];
rz(0.91845271) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.351016) q[1];
sx q[1];
rz(-1.3315397) q[1];
sx q[1];
rz(-1.3161216) q[1];
rz(-pi) q[2];
rz(-2.6076508) q[3];
sx q[3];
rz(-2.0755872) q[3];
sx q[3];
rz(-1.8478888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0489073) q[2];
sx q[2];
rz(-1.4212298) q[2];
sx q[2];
rz(0.29472026) q[2];
rz(-0.60485351) q[3];
sx q[3];
rz(-1.9113144) q[3];
sx q[3];
rz(0.10943432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7034364) q[0];
sx q[0];
rz(-1.3180757) q[0];
sx q[0];
rz(-0.99676589) q[0];
rz(-0.017771487) q[1];
sx q[1];
rz(-1.2635063) q[1];
sx q[1];
rz(1.8020188) q[1];
rz(-0.80452917) q[2];
sx q[2];
rz(-1.1297516) q[2];
sx q[2];
rz(-0.17937112) q[2];
rz(-0.65405398) q[3];
sx q[3];
rz(-0.97915836) q[3];
sx q[3];
rz(-0.98042515) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
