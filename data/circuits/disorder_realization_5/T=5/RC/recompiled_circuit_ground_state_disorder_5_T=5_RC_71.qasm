OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.65052819) q[0];
sx q[0];
rz(-1.0125546) q[0];
sx q[0];
rz(-2.2192686) q[0];
rz(2.2946279) q[1];
sx q[1];
rz(-1.4743409) q[1];
sx q[1];
rz(2.9234731) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6309752) q[0];
sx q[0];
rz(-1.4790863) q[0];
sx q[0];
rz(-0.46268483) q[0];
rz(2.7635771) q[2];
sx q[2];
rz(-1.5326008) q[2];
sx q[2];
rz(0.60249828) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.74789819) q[1];
sx q[1];
rz(-1.5285349) q[1];
sx q[1];
rz(1.7903922) q[1];
rz(2.0006337) q[3];
sx q[3];
rz(-1.4996734) q[3];
sx q[3];
rz(-1.5971368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.073079022) q[2];
sx q[2];
rz(-1.928669) q[2];
sx q[2];
rz(-0.53654137) q[2];
rz(-2.1327298) q[3];
sx q[3];
rz(-0.77102414) q[3];
sx q[3];
rz(-2.3276276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63758481) q[0];
sx q[0];
rz(-1.8300087) q[0];
sx q[0];
rz(-3.1245533) q[0];
rz(-2.0783966) q[1];
sx q[1];
rz(-2.3737962) q[1];
sx q[1];
rz(0.95471901) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3957152) q[0];
sx q[0];
rz(-1.7999188) q[0];
sx q[0];
rz(-1.627501) q[0];
rz(-pi) q[1];
rz(-0.164523) q[2];
sx q[2];
rz(-1.2753107) q[2];
sx q[2];
rz(2.4128259) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3054637) q[1];
sx q[1];
rz(-1.3333798) q[1];
sx q[1];
rz(1.8029455) q[1];
x q[2];
rz(1.9185136) q[3];
sx q[3];
rz(-0.056611185) q[3];
sx q[3];
rz(2.268102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9118328) q[2];
sx q[2];
rz(-2.2288897) q[2];
sx q[2];
rz(1.1141874) q[2];
rz(-1.1344502) q[3];
sx q[3];
rz(-2.9454234) q[3];
sx q[3];
rz(-0.1964868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.590362) q[0];
sx q[0];
rz(-2.1901665) q[0];
sx q[0];
rz(2.9587342) q[0];
rz(3.0602449) q[1];
sx q[1];
rz(-0.75129879) q[1];
sx q[1];
rz(-2.523211) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.334916) q[0];
sx q[0];
rz(-2.3551919) q[0];
sx q[0];
rz(-1.4674835) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76128086) q[2];
sx q[2];
rz(-1.4321186) q[2];
sx q[2];
rz(-1.5468462) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.74677361) q[1];
sx q[1];
rz(-1.2319984) q[1];
sx q[1];
rz(-0.14751409) q[1];
x q[2];
rz(-0.86520536) q[3];
sx q[3];
rz(-2.3045447) q[3];
sx q[3];
rz(3.1045543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0028093) q[2];
sx q[2];
rz(-1.323779) q[2];
sx q[2];
rz(-0.36661026) q[2];
rz(-1.2159411) q[3];
sx q[3];
rz(-1.1953019) q[3];
sx q[3];
rz(-2.478157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65130305) q[0];
sx q[0];
rz(-1.2325352) q[0];
sx q[0];
rz(-0.7269727) q[0];
rz(-2.2333249) q[1];
sx q[1];
rz(-2.5636702) q[1];
sx q[1];
rz(1.04331) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20651992) q[0];
sx q[0];
rz(-1.5653725) q[0];
sx q[0];
rz(-0.66117735) q[0];
rz(1.1953148) q[2];
sx q[2];
rz(-1.6282778) q[2];
sx q[2];
rz(-0.45141803) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.336652) q[1];
sx q[1];
rz(-2.0139317) q[1];
sx q[1];
rz(0.27733485) q[1];
rz(-pi) q[2];
rz(0.38049145) q[3];
sx q[3];
rz(-2.0999319) q[3];
sx q[3];
rz(-2.7545415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.72994453) q[2];
sx q[2];
rz(-0.32117716) q[2];
sx q[2];
rz(-1.3058861) q[2];
rz(-0.11792396) q[3];
sx q[3];
rz(-0.4685466) q[3];
sx q[3];
rz(-2.0809295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1082728) q[0];
sx q[0];
rz(-0.98618996) q[0];
sx q[0];
rz(1.6075851) q[0];
rz(-2.8458505) q[1];
sx q[1];
rz(-1.5402126) q[1];
sx q[1];
rz(1.6827513) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5285501) q[0];
sx q[0];
rz(-2.3146221) q[0];
sx q[0];
rz(-3.0463329) q[0];
rz(-pi) q[1];
rz(2.2590911) q[2];
sx q[2];
rz(-0.46445981) q[2];
sx q[2];
rz(1.8732173) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4319358) q[1];
sx q[1];
rz(-1.767358) q[1];
sx q[1];
rz(-0.86039575) q[1];
x q[2];
rz(0.98990654) q[3];
sx q[3];
rz(-1.3767929) q[3];
sx q[3];
rz(-1.7542183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7532588) q[2];
sx q[2];
rz(-0.39187852) q[2];
sx q[2];
rz(-2.4957116) q[2];
rz(-1.9258457) q[3];
sx q[3];
rz(-1.4589717) q[3];
sx q[3];
rz(1.3418044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.011768613) q[0];
sx q[0];
rz(-2.3079066) q[0];
sx q[0];
rz(2.5339793) q[0];
rz(-1.1786849) q[1];
sx q[1];
rz(-1.2858398) q[1];
sx q[1];
rz(-0.36373055) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5967412) q[0];
sx q[0];
rz(-2.5646696) q[0];
sx q[0];
rz(2.0175009) q[0];
x q[1];
rz(-1.4371488) q[2];
sx q[2];
rz(-2.8060348) q[2];
sx q[2];
rz(-0.37119532) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1235001) q[1];
sx q[1];
rz(-2.315633) q[1];
sx q[1];
rz(0.1634181) q[1];
rz(-pi) q[2];
rz(3.1167459) q[3];
sx q[3];
rz(-0.97462208) q[3];
sx q[3];
rz(2.5667584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.58934775) q[2];
sx q[2];
rz(-0.10422464) q[2];
sx q[2];
rz(1.9488526) q[2];
rz(-0.73181152) q[3];
sx q[3];
rz(-1.7164427) q[3];
sx q[3];
rz(1.653479) q[3];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8530387) q[0];
sx q[0];
rz(-2.9712501) q[0];
sx q[0];
rz(-1.6188251) q[0];
rz(1.6743926) q[1];
sx q[1];
rz(-1.8312981) q[1];
sx q[1];
rz(-2.4640962) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8988494) q[0];
sx q[0];
rz(-1.6209319) q[0];
sx q[0];
rz(-2.2671347) q[0];
rz(2.9076495) q[2];
sx q[2];
rz(-1.3877227) q[2];
sx q[2];
rz(-1.2884566) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.17854161) q[1];
sx q[1];
rz(-1.1483014) q[1];
sx q[1];
rz(0.4252379) q[1];
x q[2];
rz(0.73245184) q[3];
sx q[3];
rz(-1.6897413) q[3];
sx q[3];
rz(-0.67057395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7447394) q[2];
sx q[2];
rz(-0.91529673) q[2];
sx q[2];
rz(-1.8358561) q[2];
rz(0.33852494) q[3];
sx q[3];
rz(-2.0207696) q[3];
sx q[3];
rz(-2.4566076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5386388) q[0];
sx q[0];
rz(-2.9260577) q[0];
sx q[0];
rz(-2.9576874) q[0];
rz(-1.454608) q[1];
sx q[1];
rz(-2.589476) q[1];
sx q[1];
rz(-2.6457381) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82168285) q[0];
sx q[0];
rz(-1.3047072) q[0];
sx q[0];
rz(-1.0644819) q[0];
rz(-0.74842986) q[2];
sx q[2];
rz(-2.1572621) q[2];
sx q[2];
rz(0.85916729) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.42577463) q[1];
sx q[1];
rz(-1.1905021) q[1];
sx q[1];
rz(-2.7575995) q[1];
rz(-pi) q[2];
rz(-1.3794231) q[3];
sx q[3];
rz(-2.5809605) q[3];
sx q[3];
rz(-2.3440685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0869861) q[2];
sx q[2];
rz(-0.84422529) q[2];
sx q[2];
rz(1.1620713) q[2];
rz(1.453513) q[3];
sx q[3];
rz(-0.96265692) q[3];
sx q[3];
rz(-1.8614205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.42089713) q[0];
sx q[0];
rz(-1.1524042) q[0];
sx q[0];
rz(0.56588093) q[0];
rz(-0.3903009) q[1];
sx q[1];
rz(-1.5719599) q[1];
sx q[1];
rz(1.3077259) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6193251) q[0];
sx q[0];
rz(-0.59460708) q[0];
sx q[0];
rz(-0.99647605) q[0];
rz(-pi) q[1];
rz(2.9290767) q[2];
sx q[2];
rz(-1.2005998) q[2];
sx q[2];
rz(-3.1345362) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2016328) q[1];
sx q[1];
rz(-1.8324513) q[1];
sx q[1];
rz(-2.5152339) q[1];
x q[2];
rz(1.4212178) q[3];
sx q[3];
rz(-1.6892551) q[3];
sx q[3];
rz(-1.7427904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.86964837) q[2];
sx q[2];
rz(-0.41676909) q[2];
sx q[2];
rz(-2.1040253) q[2];
rz(-1.3197673) q[3];
sx q[3];
rz(-2.3128553) q[3];
sx q[3];
rz(-2.493609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2607525) q[0];
sx q[0];
rz(-0.97761959) q[0];
sx q[0];
rz(0.64224893) q[0];
rz(-1.3308659) q[1];
sx q[1];
rz(-0.86527491) q[1];
sx q[1];
rz(-2.9790402) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35643043) q[0];
sx q[0];
rz(-2.2318062) q[0];
sx q[0];
rz(2.6399773) q[0];
rz(2.7621848) q[2];
sx q[2];
rz(-1.2343293) q[2];
sx q[2];
rz(-3.0010146) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0916718) q[1];
sx q[1];
rz(-1.0921879) q[1];
sx q[1];
rz(-2.9438627) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.58329669) q[3];
sx q[3];
rz(-1.5044893) q[3];
sx q[3];
rz(-0.60140677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6315397) q[2];
sx q[2];
rz(-1.8659464) q[2];
sx q[2];
rz(-1.0408545) q[2];
rz(-0.96796525) q[3];
sx q[3];
rz(-2.0716045) q[3];
sx q[3];
rz(2.4805243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80617245) q[0];
sx q[0];
rz(-1.5652884) q[0];
sx q[0];
rz(-0.096927222) q[0];
rz(2.9087635) q[1];
sx q[1];
rz(-2.1131344) q[1];
sx q[1];
rz(0.0075385787) q[1];
rz(-0.54388028) q[2];
sx q[2];
rz(-2.3162622) q[2];
sx q[2];
rz(0.62427229) q[2];
rz(-2.6951058) q[3];
sx q[3];
rz(-2.3904908) q[3];
sx q[3];
rz(-0.6294546) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
