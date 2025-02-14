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
rz(4.1069694) q[0];
sx q[0];
rz(10.341636) q[0];
rz(1.4404453) q[1];
sx q[1];
rz(-1.0280161) q[1];
sx q[1];
rz(2.4457959) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3549933) q[0];
sx q[0];
rz(-0.43012857) q[0];
sx q[0];
rz(-1.7010078) q[0];
x q[1];
rz(-0.037362413) q[2];
sx q[2];
rz(-1.6184095) q[2];
sx q[2];
rz(1.7689153) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.88118824) q[1];
sx q[1];
rz(-1.1965479) q[1];
sx q[1];
rz(1.3874235) q[1];
rz(2.0654109) q[3];
sx q[3];
rz(-1.8997703) q[3];
sx q[3];
rz(-1.2124202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.15771389) q[2];
sx q[2];
rz(-2.3215051) q[2];
sx q[2];
rz(2.2411818) q[2];
rz(-0.79312098) q[3];
sx q[3];
rz(-1.6499062) q[3];
sx q[3];
rz(2.4429863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6124436) q[0];
sx q[0];
rz(-1.9731033) q[0];
sx q[0];
rz(0.68552619) q[0];
rz(0.43288747) q[1];
sx q[1];
rz(-1.2665749) q[1];
sx q[1];
rz(-1.3139075) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.092201) q[0];
sx q[0];
rz(-0.77038465) q[0];
sx q[0];
rz(-2.3197473) q[0];
x q[1];
rz(2.6816112) q[2];
sx q[2];
rz(-0.88085266) q[2];
sx q[2];
rz(-1.8073818) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8364468) q[1];
sx q[1];
rz(-2.0550108) q[1];
sx q[1];
rz(0.53026235) q[1];
rz(-pi) q[2];
rz(0.55696662) q[3];
sx q[3];
rz(-1.3560221) q[3];
sx q[3];
rz(1.0207411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.96585298) q[2];
sx q[2];
rz(-0.41030914) q[2];
sx q[2];
rz(-0.24743323) q[2];
rz(0.11582173) q[3];
sx q[3];
rz(-1.7146866) q[3];
sx q[3];
rz(-2.5628832) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1568569) q[0];
sx q[0];
rz(-0.37608376) q[0];
sx q[0];
rz(-2.0767186) q[0];
rz(-0.6330601) q[1];
sx q[1];
rz(-1.9828911) q[1];
sx q[1];
rz(0.13654581) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7744903) q[0];
sx q[0];
rz(-2.3068412) q[0];
sx q[0];
rz(1.537905) q[0];
rz(-pi) q[1];
rz(-0.94046392) q[2];
sx q[2];
rz(-0.55359888) q[2];
sx q[2];
rz(-2.3950837) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1376129) q[1];
sx q[1];
rz(-0.93055994) q[1];
sx q[1];
rz(0.97059393) q[1];
rz(2.7259215) q[3];
sx q[3];
rz(-1.9686254) q[3];
sx q[3];
rz(-1.7818539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.69918767) q[2];
sx q[2];
rz(-1.8892611) q[2];
sx q[2];
rz(1.5415972) q[2];
rz(-2.048118) q[3];
sx q[3];
rz(-1.6396061) q[3];
sx q[3];
rz(2.7961965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71444756) q[0];
sx q[0];
rz(-0.26145014) q[0];
sx q[0];
rz(-2.5647822) q[0];
rz(-2.414074) q[1];
sx q[1];
rz(-1.0914165) q[1];
sx q[1];
rz(2.2732546) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5142347) q[0];
sx q[0];
rz(-1.1749724) q[0];
sx q[0];
rz(2.4577599) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.98855726) q[2];
sx q[2];
rz(-2.0209656) q[2];
sx q[2];
rz(0.52272138) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5440923) q[1];
sx q[1];
rz(-2.3832633) q[1];
sx q[1];
rz(0.51095283) q[1];
x q[2];
rz(-0.9138388) q[3];
sx q[3];
rz(-2.4025177) q[3];
sx q[3];
rz(2.5881059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.19520983) q[2];
sx q[2];
rz(-1.3888487) q[2];
sx q[2];
rz(2.4460804) q[2];
rz(-2.9269311) q[3];
sx q[3];
rz(-1.6404459) q[3];
sx q[3];
rz(2.3282839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5401841) q[0];
sx q[0];
rz(-2.5267127) q[0];
sx q[0];
rz(-1.1001128) q[0];
rz(-0.21057883) q[1];
sx q[1];
rz(-2.5081878) q[1];
sx q[1];
rz(1.4932154) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8029047) q[0];
sx q[0];
rz(-0.7006104) q[0];
sx q[0];
rz(0.79498298) q[0];
rz(-pi) q[1];
rz(1.7488241) q[2];
sx q[2];
rz(-1.653394) q[2];
sx q[2];
rz(-0.96197739) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3201218) q[1];
sx q[1];
rz(-1.7263328) q[1];
sx q[1];
rz(-1.0782773) q[1];
rz(-pi) q[2];
rz(-0.3143309) q[3];
sx q[3];
rz(-2.4498457) q[3];
sx q[3];
rz(1.421613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4139999) q[2];
sx q[2];
rz(-1.7121345) q[2];
sx q[2];
rz(-0.89699888) q[2];
rz(1.8619079) q[3];
sx q[3];
rz(-2.1309659) q[3];
sx q[3];
rz(3.0943833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.9996662) q[0];
sx q[0];
rz(-2.5157295) q[0];
sx q[0];
rz(-2.8114787) q[0];
rz(-0.64078981) q[1];
sx q[1];
rz(-1.7658486) q[1];
sx q[1];
rz(-2.1660588) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3307083) q[0];
sx q[0];
rz(-2.160794) q[0];
sx q[0];
rz(1.5025274) q[0];
x q[1];
rz(-2.8728666) q[2];
sx q[2];
rz(-1.8789504) q[2];
sx q[2];
rz(1.9374963) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2157877) q[1];
sx q[1];
rz(-0.5277718) q[1];
sx q[1];
rz(0.58234071) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9532072) q[3];
sx q[3];
rz(-1.1551765) q[3];
sx q[3];
rz(1.3481026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8581802) q[2];
sx q[2];
rz(-2.8874669) q[2];
sx q[2];
rz(2.251909) q[2];
rz(2.5455918) q[3];
sx q[3];
rz(-1.069205) q[3];
sx q[3];
rz(3.0349351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0058873) q[0];
sx q[0];
rz(-1.1361253) q[0];
sx q[0];
rz(-2.0086052) q[0];
rz(0.63944447) q[1];
sx q[1];
rz(-1.0777799) q[1];
sx q[1];
rz(-2.1741672) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9854162) q[0];
sx q[0];
rz(-1.7207017) q[0];
sx q[0];
rz(-1.0557337) q[0];
x q[1];
rz(2.4339811) q[2];
sx q[2];
rz(-0.48383265) q[2];
sx q[2];
rz(2.8255723) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6220807) q[1];
sx q[1];
rz(-1.685962) q[1];
sx q[1];
rz(0.22611109) q[1];
x q[2];
rz(-0.65770517) q[3];
sx q[3];
rz(-1.8175565) q[3];
sx q[3];
rz(1.9384428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7449259) q[2];
sx q[2];
rz(-1.0090088) q[2];
sx q[2];
rz(0.94433707) q[2];
rz(-0.62263292) q[3];
sx q[3];
rz(-0.98832744) q[3];
sx q[3];
rz(-0.78701377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054976376) q[0];
sx q[0];
rz(-2.3353307) q[0];
sx q[0];
rz(-0.0042313519) q[0];
rz(1.4684234) q[1];
sx q[1];
rz(-2.3783042) q[1];
sx q[1];
rz(-2.9671293) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2989104) q[0];
sx q[0];
rz(-1.7571661) q[0];
sx q[0];
rz(3.0475265) q[0];
x q[1];
rz(2.6930843) q[2];
sx q[2];
rz(-1.6629728) q[2];
sx q[2];
rz(-1.1223457) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6721539) q[1];
sx q[1];
rz(-1.6368333) q[1];
sx q[1];
rz(-2.2697637) q[1];
rz(-pi) q[2];
rz(-1.3346436) q[3];
sx q[3];
rz(-1.6340578) q[3];
sx q[3];
rz(1.0452273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8909495) q[2];
sx q[2];
rz(-3.0477016) q[2];
sx q[2];
rz(-0.54035652) q[2];
rz(-0.017102608) q[3];
sx q[3];
rz(-0.98186866) q[3];
sx q[3];
rz(2.4798149) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0509725) q[0];
sx q[0];
rz(-2.3391188) q[0];
sx q[0];
rz(0.78659868) q[0];
rz(1.1057009) q[1];
sx q[1];
rz(-1.6603755) q[1];
sx q[1];
rz(2.9203663) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3174789) q[0];
sx q[0];
rz(-0.86659089) q[0];
sx q[0];
rz(0.368035) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7443666) q[2];
sx q[2];
rz(-1.6643057) q[2];
sx q[2];
rz(-1.5154523) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7936642) q[1];
sx q[1];
rz(-1.778144) q[1];
sx q[1];
rz(2.6479032) q[1];
x q[2];
rz(-2.4513793) q[3];
sx q[3];
rz(-1.3891474) q[3];
sx q[3];
rz(0.51723328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1456387) q[2];
sx q[2];
rz(-2.8990539) q[2];
sx q[2];
rz(1.9560187) q[2];
rz(-2.0354185) q[3];
sx q[3];
rz(-1.0217051) q[3];
sx q[3];
rz(-2.1574028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36158654) q[0];
sx q[0];
rz(-1.3573283) q[0];
sx q[0];
rz(0.7097882) q[0];
rz(-2.2478726) q[1];
sx q[1];
rz(-0.62785405) q[1];
sx q[1];
rz(2.6745904) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2115973) q[0];
sx q[0];
rz(-0.30147895) q[0];
sx q[0];
rz(-2.0512352) q[0];
rz(-2.1960718) q[2];
sx q[2];
rz(-1.3511124) q[2];
sx q[2];
rz(1.4131119) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4255387) q[1];
sx q[1];
rz(-1.3322833) q[1];
sx q[1];
rz(-0.038360049) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.8740467) q[3];
sx q[3];
rz(-1.0206887) q[3];
sx q[3];
rz(-2.3693905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3422157) q[2];
sx q[2];
rz(-2.2365384) q[2];
sx q[2];
rz(-3.0055935) q[2];
rz(0.32986766) q[3];
sx q[3];
rz(-2.0743275) q[3];
sx q[3];
rz(-2.8444667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8316523) q[0];
sx q[0];
rz(-1.1158054) q[0];
sx q[0];
rz(0.62624395) q[0];
rz(-2.2533439) q[1];
sx q[1];
rz(-1.2445969) q[1];
sx q[1];
rz(0.1319763) q[1];
rz(-1.7952948) q[2];
sx q[2];
rz(-1.1830876) q[2];
sx q[2];
rz(3.0814688) q[2];
rz(0.72867758) q[3];
sx q[3];
rz(-1.8986438) q[3];
sx q[3];
rz(1.8713634) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
