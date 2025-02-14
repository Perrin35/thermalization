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
rz(0.8402549) q[0];
sx q[0];
rz(4.1794887) q[0];
sx q[0];
rz(10.269796) q[0];
rz(-2.4294699) q[1];
sx q[1];
rz(-1.0016088) q[1];
sx q[1];
rz(-1.4955624) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5069731) q[0];
sx q[0];
rz(-1.9911626) q[0];
sx q[0];
rz(-2.2078321) q[0];
rz(0.35829131) q[2];
sx q[2];
rz(-0.40696496) q[2];
sx q[2];
rz(-0.55487061) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.597376) q[1];
sx q[1];
rz(-2.4565182) q[1];
sx q[1];
rz(2.3822576) q[1];
rz(-pi) q[2];
rz(-1.846662) q[3];
sx q[3];
rz(-1.1893236) q[3];
sx q[3];
rz(-1.671553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.40319765) q[2];
sx q[2];
rz(-1.1824111) q[2];
sx q[2];
rz(-0.5564059) q[2];
rz(2.6465936) q[3];
sx q[3];
rz(-2.7904816) q[3];
sx q[3];
rz(2.2569412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86968386) q[0];
sx q[0];
rz(-0.10232919) q[0];
sx q[0];
rz(0.62865692) q[0];
rz(-2.2643845) q[1];
sx q[1];
rz(-0.49867201) q[1];
sx q[1];
rz(-0.25310755) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49324711) q[0];
sx q[0];
rz(-1.5551994) q[0];
sx q[0];
rz(-0.007997917) q[0];
rz(-pi) q[1];
rz(-1.8688525) q[2];
sx q[2];
rz(-1.9949081) q[2];
sx q[2];
rz(-1.9587245) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4483466) q[1];
sx q[1];
rz(-2.9102059) q[1];
sx q[1];
rz(-0.53129249) q[1];
rz(-pi) q[2];
rz(-2.3023002) q[3];
sx q[3];
rz(-1.7613693) q[3];
sx q[3];
rz(-2.1803602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6796598) q[2];
sx q[2];
rz(-1.2075281) q[2];
sx q[2];
rz(2.0750462) q[2];
rz(-1.8343532) q[3];
sx q[3];
rz(-2.6756838) q[3];
sx q[3];
rz(0.90827847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9268554) q[0];
sx q[0];
rz(-2.1934788) q[0];
sx q[0];
rz(2.1562449) q[0];
rz(-0.39380479) q[1];
sx q[1];
rz(-2.1486798) q[1];
sx q[1];
rz(-0.52678144) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.55888) q[0];
sx q[0];
rz(-1.0436907) q[0];
sx q[0];
rz(-3.0452791) q[0];
x q[1];
rz(0.57931945) q[2];
sx q[2];
rz(-1.123482) q[2];
sx q[2];
rz(-1.2944702) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.56934988) q[1];
sx q[1];
rz(-1.4016289) q[1];
sx q[1];
rz(1.9985241) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6963077) q[3];
sx q[3];
rz(-2.1213343) q[3];
sx q[3];
rz(-2.8360644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7963316) q[2];
sx q[2];
rz(-2.3806206) q[2];
sx q[2];
rz(2.995028) q[2];
rz(-2.2740299) q[3];
sx q[3];
rz(-2.2680794) q[3];
sx q[3];
rz(-2.383702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4124311) q[0];
sx q[0];
rz(-0.49598345) q[0];
sx q[0];
rz(-0.55369401) q[0];
rz(-1.6665392) q[1];
sx q[1];
rz(-2.8429884) q[1];
sx q[1];
rz(2.8972304) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.037652) q[0];
sx q[0];
rz(-2.6459011) q[0];
sx q[0];
rz(-2.3858059) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3073439) q[2];
sx q[2];
rz(-1.3400199) q[2];
sx q[2];
rz(1.007387) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7556954) q[1];
sx q[1];
rz(-1.2127969) q[1];
sx q[1];
rz(-0.44803195) q[1];
rz(-pi) q[2];
rz(-1.5958173) q[3];
sx q[3];
rz(-0.85880781) q[3];
sx q[3];
rz(-0.98464363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0814521) q[2];
sx q[2];
rz(-1.1271366) q[2];
sx q[2];
rz(-0.78090182) q[2];
rz(2.7447356) q[3];
sx q[3];
rz(-3.1271827) q[3];
sx q[3];
rz(2.0675596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9020554) q[0];
sx q[0];
rz(-2.9750415) q[0];
sx q[0];
rz(-2.8848414) q[0];
rz(-0.17770879) q[1];
sx q[1];
rz(-2.5959028) q[1];
sx q[1];
rz(2.1977052) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2985757) q[0];
sx q[0];
rz(-1.1305048) q[0];
sx q[0];
rz(0.40233516) q[0];
rz(0.98942049) q[2];
sx q[2];
rz(-2.7983449) q[2];
sx q[2];
rz(-0.95232576) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0066932) q[1];
sx q[1];
rz(-1.0311145) q[1];
sx q[1];
rz(-1.5758908) q[1];
rz(-pi) q[2];
rz(-2.1274756) q[3];
sx q[3];
rz(-1.1667487) q[3];
sx q[3];
rz(-1.7168644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.93861598) q[2];
sx q[2];
rz(-2.1111033) q[2];
sx q[2];
rz(0.20982783) q[2];
rz(-0.046791568) q[3];
sx q[3];
rz(-1.6822466) q[3];
sx q[3];
rz(-0.34745026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0824579) q[0];
sx q[0];
rz(-0.81675285) q[0];
sx q[0];
rz(-1.2030075) q[0];
rz(1.5164392) q[1];
sx q[1];
rz(-0.37767437) q[1];
sx q[1];
rz(-0.032940544) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66844475) q[0];
sx q[0];
rz(-2.0347036) q[0];
sx q[0];
rz(2.8345246) q[0];
rz(-pi) q[1];
x q[1];
rz(0.98349468) q[2];
sx q[2];
rz(-2.0222531) q[2];
sx q[2];
rz(1.2304103) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6194544) q[1];
sx q[1];
rz(-0.65548766) q[1];
sx q[1];
rz(0.54490276) q[1];
x q[2];
rz(2.2516656) q[3];
sx q[3];
rz(-1.0491199) q[3];
sx q[3];
rz(1.8311178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5861627) q[2];
sx q[2];
rz(-1.0032434) q[2];
sx q[2];
rz(-2.6891151) q[2];
rz(0.14719506) q[3];
sx q[3];
rz(-0.23284027) q[3];
sx q[3];
rz(-1.2991306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.074987) q[0];
sx q[0];
rz(-2.7050278) q[0];
sx q[0];
rz(2.7845352) q[0];
rz(2.1287411) q[1];
sx q[1];
rz(-1.169299) q[1];
sx q[1];
rz(-3.1049407) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79082876) q[0];
sx q[0];
rz(-2.0022656) q[0];
sx q[0];
rz(-0.58833265) q[0];
rz(-pi) q[1];
rz(-0.22256644) q[2];
sx q[2];
rz(-0.43912008) q[2];
sx q[2];
rz(-2.665208) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4067024) q[1];
sx q[1];
rz(-1.6631366) q[1];
sx q[1];
rz(-0.86812302) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5662848) q[3];
sx q[3];
rz(-1.6016212) q[3];
sx q[3];
rz(0.84336057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6645633) q[2];
sx q[2];
rz(-2.1881115) q[2];
sx q[2];
rz(-3.0518517) q[2];
rz(1.2612032) q[3];
sx q[3];
rz(-1.829105) q[3];
sx q[3];
rz(-1.7173654) q[3];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66202128) q[0];
sx q[0];
rz(-0.57858545) q[0];
sx q[0];
rz(-1.660996) q[0];
rz(1.6914233) q[1];
sx q[1];
rz(-2.4822576) q[1];
sx q[1];
rz(-2.3816542) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63045365) q[0];
sx q[0];
rz(-1.5811111) q[0];
sx q[0];
rz(-1.6869839) q[0];
x q[1];
rz(-1.6641812) q[2];
sx q[2];
rz(-2.1969446) q[2];
sx q[2];
rz(-1.6531675) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2719136) q[1];
sx q[1];
rz(-1.1506423) q[1];
sx q[1];
rz(0.43964444) q[1];
rz(0.44260232) q[3];
sx q[3];
rz(-0.55130225) q[3];
sx q[3];
rz(-2.9568903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.30457589) q[2];
sx q[2];
rz(-1.872007) q[2];
sx q[2];
rz(3.0820091) q[2];
rz(-2.7130821) q[3];
sx q[3];
rz(-0.84463745) q[3];
sx q[3];
rz(-0.60215157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.2483599) q[0];
sx q[0];
rz(-2.8587274) q[0];
sx q[0];
rz(-2.7981113) q[0];
rz(2.9454625) q[1];
sx q[1];
rz(-0.88534147) q[1];
sx q[1];
rz(2.3416065) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6340187) q[0];
sx q[0];
rz(-2.2306666) q[0];
sx q[0];
rz(0.52652653) q[0];
x q[1];
rz(-0.56234151) q[2];
sx q[2];
rz(-1.7844756) q[2];
sx q[2];
rz(-2.4116229) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1175122) q[1];
sx q[1];
rz(-1.615106) q[1];
sx q[1];
rz(2.5818392) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2535747) q[3];
sx q[3];
rz(-1.6245442) q[3];
sx q[3];
rz(0.88784224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7949152) q[2];
sx q[2];
rz(-2.3778043) q[2];
sx q[2];
rz(-2.9128892) q[2];
rz(-1.7250569) q[3];
sx q[3];
rz(-1.4226457) q[3];
sx q[3];
rz(-0.67030877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59920853) q[0];
sx q[0];
rz(-0.51775652) q[0];
sx q[0];
rz(1.1826578) q[0];
rz(0.54284894) q[1];
sx q[1];
rz(-3.019637) q[1];
sx q[1];
rz(-0.45546946) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1802954) q[0];
sx q[0];
rz(-2.0019128) q[0];
sx q[0];
rz(-2.6102553) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9331839) q[2];
sx q[2];
rz(-2.0129856) q[2];
sx q[2];
rz(0.65763523) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.093468637) q[1];
sx q[1];
rz(-2.6542722) q[1];
sx q[1];
rz(0.29429277) q[1];
x q[2];
rz(-0.72584589) q[3];
sx q[3];
rz(-0.56768394) q[3];
sx q[3];
rz(2.3633702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.21696422) q[2];
sx q[2];
rz(-2.2483726) q[2];
sx q[2];
rz(2.7157937) q[2];
rz(-0.18481542) q[3];
sx q[3];
rz(-2.5033689) q[3];
sx q[3];
rz(-0.027801175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7040779) q[0];
sx q[0];
rz(-1.6117493) q[0];
sx q[0];
rz(3.1060863) q[0];
rz(2.6620445) q[1];
sx q[1];
rz(-1.5621114) q[1];
sx q[1];
rz(1.5893804) q[1];
rz(1.7935971) q[2];
sx q[2];
rz(-0.98735129) q[2];
sx q[2];
rz(-2.9410887) q[2];
rz(0.95994031) q[3];
sx q[3];
rz(-1.452594) q[3];
sx q[3];
rz(0.95127524) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
