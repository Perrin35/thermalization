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
rz(-2.8060198) q[0];
sx q[0];
rz(-1.6874474) q[0];
sx q[0];
rz(2.1079221) q[0];
rz(1.4312862) q[1];
sx q[1];
rz(-0.55748504) q[1];
sx q[1];
rz(2.1318336) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.566322) q[0];
sx q[0];
rz(-2.913544) q[0];
sx q[0];
rz(-2.1581677) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6160585) q[2];
sx q[2];
rz(-3.0590364) q[2];
sx q[2];
rz(-0.98708594) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3897455) q[1];
sx q[1];
rz(-1.455014) q[1];
sx q[1];
rz(0.16645653) q[1];
rz(-pi) q[2];
rz(0.030850284) q[3];
sx q[3];
rz(-2.9992963) q[3];
sx q[3];
rz(-0.22103413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3681763) q[2];
sx q[2];
rz(-2.2645617) q[2];
sx q[2];
rz(2.3167493) q[2];
rz(-0.33251897) q[3];
sx q[3];
rz(-0.85086346) q[3];
sx q[3];
rz(-0.48377812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2291718) q[0];
sx q[0];
rz(-3.0569172) q[0];
sx q[0];
rz(-2.605865) q[0];
rz(2.3748705) q[1];
sx q[1];
rz(-1.7886536) q[1];
sx q[1];
rz(-1.3923233) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.509393) q[0];
sx q[0];
rz(-2.0344816) q[0];
sx q[0];
rz(-2.7576588) q[0];
rz(3.1145623) q[2];
sx q[2];
rz(-0.34124869) q[2];
sx q[2];
rz(-0.42021593) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.079683144) q[1];
sx q[1];
rz(-2.3259729) q[1];
sx q[1];
rz(-0.029600005) q[1];
x q[2];
rz(1.4658607) q[3];
sx q[3];
rz(-1.6989048) q[3];
sx q[3];
rz(-0.48612938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2367737) q[2];
sx q[2];
rz(-0.83319131) q[2];
sx q[2];
rz(1.0986249) q[2];
rz(2.3084579) q[3];
sx q[3];
rz(-1.5435217) q[3];
sx q[3];
rz(-2.0875077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58040923) q[0];
sx q[0];
rz(-1.1495178) q[0];
sx q[0];
rz(-0.68823254) q[0];
rz(-1.2511823) q[1];
sx q[1];
rz(-2.4871608) q[1];
sx q[1];
rz(1.2122663) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0949757) q[0];
sx q[0];
rz(-0.88030926) q[0];
sx q[0];
rz(0.63027905) q[0];
rz(2.278944) q[2];
sx q[2];
rz(-0.35637636) q[2];
sx q[2];
rz(-1.4087806) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9759298) q[1];
sx q[1];
rz(-1.5651795) q[1];
sx q[1];
rz(2.4826239) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5098677) q[3];
sx q[3];
rz(-0.079515545) q[3];
sx q[3];
rz(-0.36859504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4091829) q[2];
sx q[2];
rz(-0.99486351) q[2];
sx q[2];
rz(0.44962064) q[2];
rz(-2.2879587) q[3];
sx q[3];
rz(-2.4946404) q[3];
sx q[3];
rz(2.4322815) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7483826) q[0];
sx q[0];
rz(-2.125183) q[0];
sx q[0];
rz(0.40312314) q[0];
rz(2.368811) q[1];
sx q[1];
rz(-1.2330331) q[1];
sx q[1];
rz(-2.3162139) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51096237) q[0];
sx q[0];
rz(-1.2494506) q[0];
sx q[0];
rz(0.83842959) q[0];
x q[1];
rz(1.5947444) q[2];
sx q[2];
rz(-1.5976816) q[2];
sx q[2];
rz(1.250911) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1267438) q[1];
sx q[1];
rz(-0.9909317) q[1];
sx q[1];
rz(0.090437263) q[1];
rz(-pi) q[2];
x q[2];
rz(0.89572923) q[3];
sx q[3];
rz(-0.78215996) q[3];
sx q[3];
rz(1.6767927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1021425) q[2];
sx q[2];
rz(-1.8723698) q[2];
sx q[2];
rz(2.6329182) q[2];
rz(-1.8968808) q[3];
sx q[3];
rz(-2.0679943) q[3];
sx q[3];
rz(1.3332453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4077069) q[0];
sx q[0];
rz(-1.5835967) q[0];
sx q[0];
rz(-0.35633126) q[0];
rz(-1.7286667) q[1];
sx q[1];
rz(-1.3099542) q[1];
sx q[1];
rz(1.8467356) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0644449) q[0];
sx q[0];
rz(-2.003621) q[0];
sx q[0];
rz(-1.3291784) q[0];
rz(-0.6614954) q[2];
sx q[2];
rz(-1.7057888) q[2];
sx q[2];
rz(2.1141426) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.19902557) q[1];
sx q[1];
rz(-0.46442214) q[1];
sx q[1];
rz(-0.53408043) q[1];
rz(-pi) q[2];
rz(-2.1160841) q[3];
sx q[3];
rz(-1.9615615) q[3];
sx q[3];
rz(0.39749591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.20092189) q[2];
sx q[2];
rz(-1.385043) q[2];
sx q[2];
rz(-2.5015639) q[2];
rz(2.707543) q[3];
sx q[3];
rz(-0.95165747) q[3];
sx q[3];
rz(-2.0210338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5459179) q[0];
sx q[0];
rz(-0.43788236) q[0];
sx q[0];
rz(-0.75089279) q[0];
rz(-0.76260507) q[1];
sx q[1];
rz(-1.7060988) q[1];
sx q[1];
rz(2.6079752) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2180839) q[0];
sx q[0];
rz(-1.2752011) q[0];
sx q[0];
rz(1.2696928) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.80383382) q[2];
sx q[2];
rz(-1.8920004) q[2];
sx q[2];
rz(2.8558232) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.68249629) q[1];
sx q[1];
rz(-0.67855103) q[1];
sx q[1];
rz(1.7595923) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5848006) q[3];
sx q[3];
rz(-2.0813312) q[3];
sx q[3];
rz(0.42312121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.077347191) q[2];
sx q[2];
rz(-0.40234819) q[2];
sx q[2];
rz(-0.8858436) q[2];
rz(1.3028076) q[3];
sx q[3];
rz(-1.9582483) q[3];
sx q[3];
rz(-0.47811374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9293905) q[0];
sx q[0];
rz(-0.92340702) q[0];
sx q[0];
rz(-0.59342629) q[0];
rz(1.628283) q[1];
sx q[1];
rz(-1.1791041) q[1];
sx q[1];
rz(2.5679307) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3447269) q[0];
sx q[0];
rz(-2.0484951) q[0];
sx q[0];
rz(-1.9928689) q[0];
rz(-2.7717765) q[2];
sx q[2];
rz(-2.5313583) q[2];
sx q[2];
rz(3.1265697) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6384068) q[1];
sx q[1];
rz(-1.3246623) q[1];
sx q[1];
rz(0.41001525) q[1];
x q[2];
rz(0.078646831) q[3];
sx q[3];
rz(-1.5203195) q[3];
sx q[3];
rz(0.61928643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.71184984) q[2];
sx q[2];
rz(-1.2028799) q[2];
sx q[2];
rz(1.4329866) q[2];
rz(1.9728707) q[3];
sx q[3];
rz(-2.70372) q[3];
sx q[3];
rz(1.6247113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8716005) q[0];
sx q[0];
rz(-0.90684909) q[0];
sx q[0];
rz(0.20371833) q[0];
rz(-1.8261955) q[1];
sx q[1];
rz(-1.306465) q[1];
sx q[1];
rz(1.3320097) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34623554) q[0];
sx q[0];
rz(-2.1846844) q[0];
sx q[0];
rz(2.3625663) q[0];
rz(-pi) q[1];
x q[1];
rz(0.90317995) q[2];
sx q[2];
rz(-1.1051264) q[2];
sx q[2];
rz(-0.47448928) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1696265) q[1];
sx q[1];
rz(-2.1887767) q[1];
sx q[1];
rz(1.0849668) q[1];
rz(-2.2409954) q[3];
sx q[3];
rz(-2.934747) q[3];
sx q[3];
rz(-1.2634009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.77067644) q[2];
sx q[2];
rz(-3.0597661) q[2];
sx q[2];
rz(1.3814242) q[2];
rz(2.5904739) q[3];
sx q[3];
rz(-1.8128606) q[3];
sx q[3];
rz(-0.74530017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4296752) q[0];
sx q[0];
rz(-1.4497919) q[0];
sx q[0];
rz(1.2592738) q[0];
rz(-1.7970386) q[1];
sx q[1];
rz(-0.63251907) q[1];
sx q[1];
rz(-1.22619) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0622985) q[0];
sx q[0];
rz(-0.31656893) q[0];
sx q[0];
rz(1.1226673) q[0];
rz(2.2974469) q[2];
sx q[2];
rz(-0.54960362) q[2];
sx q[2];
rz(-1.5507825) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3955955) q[1];
sx q[1];
rz(-1.5052553) q[1];
sx q[1];
rz(-0.74742527) q[1];
rz(2.1192083) q[3];
sx q[3];
rz(-2.5520241) q[3];
sx q[3];
rz(1.1353253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5598477) q[2];
sx q[2];
rz(-0.86220828) q[2];
sx q[2];
rz(1.6017412) q[2];
rz(1.7848484) q[3];
sx q[3];
rz(-2.609085) q[3];
sx q[3];
rz(1.9073568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(0.27720472) q[0];
sx q[0];
rz(-1.5437523) q[0];
sx q[0];
rz(-0.10449115) q[0];
rz(-1.3970207) q[1];
sx q[1];
rz(-1.7897768) q[1];
sx q[1];
rz(-1.6207961) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.049690954) q[0];
sx q[0];
rz(-0.77020634) q[0];
sx q[0];
rz(2.8390604) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7303329) q[2];
sx q[2];
rz(-1.9070093) q[2];
sx q[2];
rz(-1.6415063) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0238259) q[1];
sx q[1];
rz(-1.2833047) q[1];
sx q[1];
rz(1.6790798) q[1];
rz(0.057985882) q[3];
sx q[3];
rz(-2.1515905) q[3];
sx q[3];
rz(-0.58422663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.70667679) q[2];
sx q[2];
rz(-0.84182635) q[2];
sx q[2];
rz(-1.288877) q[2];
rz(1.0051109) q[3];
sx q[3];
rz(-1.0594599) q[3];
sx q[3];
rz(3.03481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2224251) q[0];
sx q[0];
rz(-1.5821624) q[0];
sx q[0];
rz(-0.17056175) q[0];
rz(-0.87986058) q[1];
sx q[1];
rz(-0.51645551) q[1];
sx q[1];
rz(-2.5194306) q[1];
rz(-0.11779412) q[2];
sx q[2];
rz(-2.6391023) q[2];
sx q[2];
rz(-0.53822623) q[2];
rz(-0.11446791) q[3];
sx q[3];
rz(-2.538402) q[3];
sx q[3];
rz(1.7414471) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
