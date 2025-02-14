OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.652777) q[0];
sx q[0];
rz(-2.1762159) q[0];
sx q[0];
rz(0.91685796) q[0];
rz(1.4404453) q[1];
sx q[1];
rz(-1.0280161) q[1];
sx q[1];
rz(-0.69579679) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.097336107) q[0];
sx q[0];
rz(-1.5166266) q[0];
sx q[0];
rz(-1.9977117) q[0];
rz(-pi) q[1];
rz(1.6184427) q[2];
sx q[2];
rz(-1.5334763) q[2];
sx q[2];
rz(2.9452528) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.6219139) q[1];
sx q[1];
rz(-1.4002454) q[1];
sx q[1];
rz(-2.7615553) q[1];
x q[2];
rz(2.0654109) q[3];
sx q[3];
rz(-1.2418223) q[3];
sx q[3];
rz(1.2124202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9838788) q[2];
sx q[2];
rz(-2.3215051) q[2];
sx q[2];
rz(-2.2411818) q[2];
rz(0.79312098) q[3];
sx q[3];
rz(-1.4916865) q[3];
sx q[3];
rz(2.4429863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
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
rz(1.6124436) q[0];
sx q[0];
rz(-1.1684893) q[0];
sx q[0];
rz(0.68552619) q[0];
rz(-0.43288747) q[1];
sx q[1];
rz(-1.8750178) q[1];
sx q[1];
rz(-1.3139075) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0493916) q[0];
sx q[0];
rz(-0.77038465) q[0];
sx q[0];
rz(-2.3197473) q[0];
x q[1];
rz(-2.3150748) q[2];
sx q[2];
rz(-1.2213301) q[2];
sx q[2];
rz(-0.068880388) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5470299) q[1];
sx q[1];
rz(-2.4395151) q[1];
sx q[1];
rz(-0.80499332) q[1];
rz(-pi) q[2];
x q[2];
rz(0.55696662) q[3];
sx q[3];
rz(-1.3560221) q[3];
sx q[3];
rz(1.0207411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.96585298) q[2];
sx q[2];
rz(-2.7312835) q[2];
sx q[2];
rz(0.24743323) q[2];
rz(0.11582173) q[3];
sx q[3];
rz(-1.426906) q[3];
sx q[3];
rz(-0.57870948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1568569) q[0];
sx q[0];
rz(-0.37608376) q[0];
sx q[0];
rz(2.0767186) q[0];
rz(0.6330601) q[1];
sx q[1];
rz(-1.9828911) q[1];
sx q[1];
rz(3.0050468) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36710237) q[0];
sx q[0];
rz(-2.3068412) q[0];
sx q[0];
rz(-1.6036877) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.94046392) q[2];
sx q[2];
rz(-0.55359888) q[2];
sx q[2];
rz(0.74650899) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1376129) q[1];
sx q[1];
rz(-2.2110327) q[1];
sx q[1];
rz(2.1709987) q[1];
rz(-pi) q[2];
rz(1.1401922) q[3];
sx q[3];
rz(-1.1893442) q[3];
sx q[3];
rz(2.7611707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.69918767) q[2];
sx q[2];
rz(-1.2523315) q[2];
sx q[2];
rz(-1.5415972) q[2];
rz(-2.048118) q[3];
sx q[3];
rz(-1.6396061) q[3];
sx q[3];
rz(-0.34539616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71444756) q[0];
sx q[0];
rz(-0.26145014) q[0];
sx q[0];
rz(-0.57681042) q[0];
rz(2.414074) q[1];
sx q[1];
rz(-1.0914165) q[1];
sx q[1];
rz(-2.2732546) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49844102) q[0];
sx q[0];
rz(-0.77385573) q[0];
sx q[0];
rz(-0.58437785) q[0];
rz(-0.84980201) q[2];
sx q[2];
rz(-2.4219208) q[2];
sx q[2];
rz(-1.6320736) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.93993856) q[1];
sx q[1];
rz(-2.2141406) q[1];
sx q[1];
rz(-2.0045947) q[1];
rz(-pi) q[2];
rz(-2.1959325) q[3];
sx q[3];
rz(-1.994761) q[3];
sx q[3];
rz(-1.5353919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.19520983) q[2];
sx q[2];
rz(-1.752744) q[2];
sx q[2];
rz(-2.4460804) q[2];
rz(2.9269311) q[3];
sx q[3];
rz(-1.6404459) q[3];
sx q[3];
rz(0.8133088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.6014086) q[0];
sx q[0];
rz(-0.61487991) q[0];
sx q[0];
rz(2.0414798) q[0];
rz(-2.9310138) q[1];
sx q[1];
rz(-2.5081878) q[1];
sx q[1];
rz(1.6483773) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8755668) q[0];
sx q[0];
rz(-1.1023836) q[0];
sx q[0];
rz(-1.0288971) q[0];
rz(1.1334992) q[2];
sx q[2];
rz(-2.9455212) q[2];
sx q[2];
rz(-0.17889775) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4752057) q[1];
sx q[1];
rz(-2.0568486) q[1];
sx q[1];
rz(2.9654825) q[1];
rz(2.8272618) q[3];
sx q[3];
rz(-2.4498457) q[3];
sx q[3];
rz(-1.7199797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7275927) q[2];
sx q[2];
rz(-1.7121345) q[2];
sx q[2];
rz(-2.2445938) q[2];
rz(-1.2796848) q[3];
sx q[3];
rz(-2.1309659) q[3];
sx q[3];
rz(3.0943833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9996662) q[0];
sx q[0];
rz(-2.5157295) q[0];
sx q[0];
rz(-2.8114787) q[0];
rz(-2.5008028) q[1];
sx q[1];
rz(-1.375744) q[1];
sx q[1];
rz(0.97553387) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3307083) q[0];
sx q[0];
rz(-2.160794) q[0];
sx q[0];
rz(-1.6390653) q[0];
rz(2.8728666) q[2];
sx q[2];
rz(-1.2626422) q[2];
sx q[2];
rz(-1.2040963) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.16224081) q[1];
sx q[1];
rz(-1.2901514) q[1];
sx q[1];
rz(2.6885217) q[1];
x q[2];
rz(1.9930494) q[3];
sx q[3];
rz(-1.7429757) q[3];
sx q[3];
rz(-2.9957221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8581802) q[2];
sx q[2];
rz(-2.8874669) q[2];
sx q[2];
rz(-2.251909) q[2];
rz(-2.5455918) q[3];
sx q[3];
rz(-1.069205) q[3];
sx q[3];
rz(0.10665756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.13570531) q[0];
sx q[0];
rz(-2.0054673) q[0];
sx q[0];
rz(-2.0086052) q[0];
rz(-0.63944447) q[1];
sx q[1];
rz(-2.0638128) q[1];
sx q[1];
rz(-2.1741672) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8113049) q[0];
sx q[0];
rz(-2.0795224) q[0];
sx q[0];
rz(2.9697493) q[0];
rz(-pi) q[1];
x q[1];
rz(0.37992875) q[2];
sx q[2];
rz(-1.8779781) q[2];
sx q[2];
rz(1.9028705) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6270476) q[1];
sx q[1];
rz(-0.25330287) q[1];
sx q[1];
rz(-0.47635081) q[1];
rz(0.65770517) q[3];
sx q[3];
rz(-1.3240361) q[3];
sx q[3];
rz(-1.2031499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7449259) q[2];
sx q[2];
rz(-1.0090088) q[2];
sx q[2];
rz(2.1972556) q[2];
rz(-0.62263292) q[3];
sx q[3];
rz(-2.1532652) q[3];
sx q[3];
rz(0.78701377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(3.0866163) q[0];
sx q[0];
rz(-2.3353307) q[0];
sx q[0];
rz(-3.1373613) q[0];
rz(1.4684234) q[1];
sx q[1];
rz(-0.76328841) q[1];
sx q[1];
rz(2.9671293) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2989104) q[0];
sx q[0];
rz(-1.7571661) q[0];
sx q[0];
rz(3.0475265) q[0];
rz(-pi) q[1];
rz(-2.9315591) q[2];
sx q[2];
rz(-0.45725098) q[2];
sx q[2];
rz(-2.5041472) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0956438) q[1];
sx q[1];
rz(-2.2679331) q[1];
sx q[1];
rz(-3.0554153) q[1];
rz(-pi) q[2];
rz(3.0765303) q[3];
sx q[3];
rz(-1.8064677) q[3];
sx q[3];
rz(-0.5103569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8909495) q[2];
sx q[2];
rz(-0.093891056) q[2];
sx q[2];
rz(0.54035652) q[2];
rz(-0.017102608) q[3];
sx q[3];
rz(-0.98186866) q[3];
sx q[3];
rz(2.4798149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0906202) q[0];
sx q[0];
rz(-2.3391188) q[0];
sx q[0];
rz(-2.354994) q[0];
rz(-2.0358918) q[1];
sx q[1];
rz(-1.4812171) q[1];
sx q[1];
rz(-2.9203663) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8241137) q[0];
sx q[0];
rz(-2.2750018) q[0];
sx q[0];
rz(-2.7735577) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7443666) q[2];
sx q[2];
rz(-1.6643057) q[2];
sx q[2];
rz(1.6261404) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.11253416) q[1];
sx q[1];
rz(-2.0529944) q[1];
sx q[1];
rz(-1.8052989) q[1];
rz(1.33696) q[3];
sx q[3];
rz(-2.2474996) q[3];
sx q[3];
rz(1.2016344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9959539) q[2];
sx q[2];
rz(-0.24253878) q[2];
sx q[2];
rz(1.185574) q[2];
rz(-1.1061741) q[3];
sx q[3];
rz(-1.0217051) q[3];
sx q[3];
rz(2.1574028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7800061) q[0];
sx q[0];
rz(-1.3573283) q[0];
sx q[0];
rz(-2.4318045) q[0];
rz(0.89372006) q[1];
sx q[1];
rz(-0.62785405) q[1];
sx q[1];
rz(2.6745904) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4295964) q[0];
sx q[0];
rz(-1.8372559) q[0];
sx q[0];
rz(-0.14273739) q[0];
rz(-pi) q[1];
rz(-2.1960718) q[2];
sx q[2];
rz(-1.3511124) q[2];
sx q[2];
rz(-1.7284808) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.71605395) q[1];
sx q[1];
rz(-1.8093093) q[1];
sx q[1];
rz(-3.1032326) q[1];
rz(0.8080838) q[3];
sx q[3];
rz(-2.2833725) q[3];
sx q[3];
rz(2.9014362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3422157) q[2];
sx q[2];
rz(-0.90505427) q[2];
sx q[2];
rz(0.13599914) q[2];
rz(-2.811725) q[3];
sx q[3];
rz(-1.0672652) q[3];
sx q[3];
rz(2.8444667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(1.8316523) q[0];
sx q[0];
rz(-2.0257873) q[0];
sx q[0];
rz(-2.5153487) q[0];
rz(-0.88824875) q[1];
sx q[1];
rz(-1.8969957) q[1];
sx q[1];
rz(-3.0096164) q[1];
rz(1.3462979) q[2];
sx q[2];
rz(-1.1830876) q[2];
sx q[2];
rz(3.0814688) q[2];
rz(-2.4129151) q[3];
sx q[3];
rz(-1.8986438) q[3];
sx q[3];
rz(1.8713634) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
