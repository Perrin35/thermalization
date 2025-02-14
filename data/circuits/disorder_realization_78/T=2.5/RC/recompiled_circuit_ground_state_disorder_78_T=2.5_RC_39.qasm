OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.4888157) q[0];
sx q[0];
rz(-0.96537679) q[0];
sx q[0];
rz(-0.91685796) q[0];
rz(1.4404453) q[1];
sx q[1];
rz(2.1135766) q[1];
sx q[1];
rz(10.120575) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7865993) q[0];
sx q[0];
rz(-0.43012857) q[0];
sx q[0];
rz(1.4405849) q[0];
rz(-pi) q[1];
rz(-1.6184427) q[2];
sx q[2];
rz(-1.6081164) q[2];
sx q[2];
rz(-0.19633987) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7909319) q[1];
sx q[1];
rz(-0.41484851) q[1];
sx q[1];
rz(-0.43465876) q[1];
rz(-pi) q[2];
rz(-2.1942646) q[3];
sx q[3];
rz(-0.58637324) q[3];
sx q[3];
rz(2.9602667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9838788) q[2];
sx q[2];
rz(-0.82008755) q[2];
sx q[2];
rz(-0.90041089) q[2];
rz(-0.79312098) q[3];
sx q[3];
rz(-1.6499062) q[3];
sx q[3];
rz(2.4429863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6124436) q[0];
sx q[0];
rz(-1.9731033) q[0];
sx q[0];
rz(-2.4560665) q[0];
rz(0.43288747) q[1];
sx q[1];
rz(-1.2665749) q[1];
sx q[1];
rz(-1.3139075) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.066931574) q[0];
sx q[0];
rz(-2.064813) q[0];
sx q[0];
rz(-0.95290174) q[0];
rz(2.3150748) q[2];
sx q[2];
rz(-1.9202625) q[2];
sx q[2];
rz(3.0727123) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.59456277) q[1];
sx q[1];
rz(-0.70207754) q[1];
sx q[1];
rz(-0.80499332) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.584626) q[3];
sx q[3];
rz(-1.7855705) q[3];
sx q[3];
rz(2.1208515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.96585298) q[2];
sx q[2];
rz(-2.7312835) q[2];
sx q[2];
rz(2.8941594) q[2];
rz(3.0257709) q[3];
sx q[3];
rz(-1.426906) q[3];
sx q[3];
rz(-2.5628832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1568569) q[0];
sx q[0];
rz(-2.7655089) q[0];
sx q[0];
rz(-1.0648741) q[0];
rz(-0.6330601) q[1];
sx q[1];
rz(-1.1587016) q[1];
sx q[1];
rz(3.0050468) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7744903) q[0];
sx q[0];
rz(-0.83475141) q[0];
sx q[0];
rz(1.6036877) q[0];
x q[1];
rz(0.94046392) q[2];
sx q[2];
rz(-2.5879938) q[2];
sx q[2];
rz(0.74650899) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0039797) q[1];
sx q[1];
rz(-0.93055994) q[1];
sx q[1];
rz(-2.1709987) q[1];
rz(-pi) q[2];
rz(-2.3362558) q[3];
sx q[3];
rz(-2.5743765) q[3];
sx q[3];
rz(-0.50931168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.442405) q[2];
sx q[2];
rz(-1.8892611) q[2];
sx q[2];
rz(-1.5415972) q[2];
rz(-1.0934746) q[3];
sx q[3];
rz(-1.5019865) q[3];
sx q[3];
rz(-0.34539616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4271451) q[0];
sx q[0];
rz(-2.8801425) q[0];
sx q[0];
rz(-2.5647822) q[0];
rz(0.72751865) q[1];
sx q[1];
rz(-1.0914165) q[1];
sx q[1];
rz(2.2732546) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5142347) q[0];
sx q[0];
rz(-1.9666202) q[0];
sx q[0];
rz(-2.4577599) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2917906) q[2];
sx q[2];
rz(-0.71967185) q[2];
sx q[2];
rz(-1.6320736) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.93993856) q[1];
sx q[1];
rz(-0.92745204) q[1];
sx q[1];
rz(2.0045947) q[1];
rz(-pi) q[2];
rz(-2.2277539) q[3];
sx q[3];
rz(-2.4025177) q[3];
sx q[3];
rz(-2.5881059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9463828) q[2];
sx q[2];
rz(-1.3888487) q[2];
sx q[2];
rz(-0.69551224) q[2];
rz(-0.21466151) q[3];
sx q[3];
rz(-1.6404459) q[3];
sx q[3];
rz(0.8133088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5401841) q[0];
sx q[0];
rz(-0.61487991) q[0];
sx q[0];
rz(2.0414798) q[0];
rz(-0.21057883) q[1];
sx q[1];
rz(-0.63340488) q[1];
sx q[1];
rz(1.6483773) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2660259) q[0];
sx q[0];
rz(-2.039209) q[0];
sx q[0];
rz(-2.1126956) q[0];
rz(-pi) q[1];
rz(1.7488241) q[2];
sx q[2];
rz(-1.4881987) q[2];
sx q[2];
rz(0.96197739) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.82147089) q[1];
sx q[1];
rz(-1.7263328) q[1];
sx q[1];
rz(2.0633153) q[1];
rz(-0.3143309) q[3];
sx q[3];
rz(-2.4498457) q[3];
sx q[3];
rz(1.421613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7275927) q[2];
sx q[2];
rz(-1.7121345) q[2];
sx q[2];
rz(0.89699888) q[2];
rz(-1.2796848) q[3];
sx q[3];
rz(-1.0106267) q[3];
sx q[3];
rz(-3.0943833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9996662) q[0];
sx q[0];
rz(-2.5157295) q[0];
sx q[0];
rz(-0.33011398) q[0];
rz(-0.64078981) q[1];
sx q[1];
rz(-1.375744) q[1];
sx q[1];
rz(2.1660588) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2084239) q[0];
sx q[0];
rz(-0.59346775) q[0];
sx q[0];
rz(3.0400601) q[0];
x q[1];
rz(1.8896721) q[2];
sx q[2];
rz(-1.826573) q[2];
sx q[2];
rz(-0.45003154) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2157877) q[1];
sx q[1];
rz(-0.5277718) q[1];
sx q[1];
rz(0.58234071) q[1];
rz(-pi) q[2];
rz(1.9930494) q[3];
sx q[3];
rz(-1.3986169) q[3];
sx q[3];
rz(-0.14587054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8581802) q[2];
sx q[2];
rz(-2.8874669) q[2];
sx q[2];
rz(2.251909) q[2];
rz(2.5455918) q[3];
sx q[3];
rz(-2.0723876) q[3];
sx q[3];
rz(0.10665756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13570531) q[0];
sx q[0];
rz(-1.1361253) q[0];
sx q[0];
rz(-2.0086052) q[0];
rz(-2.5021482) q[1];
sx q[1];
rz(-1.0777799) q[1];
sx q[1];
rz(-2.1741672) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4690034) q[0];
sx q[0];
rz(-2.6070507) q[0];
sx q[0];
rz(1.8683167) q[0];
x q[1];
rz(1.899951) q[2];
sx q[2];
rz(-1.2094922) q[2];
sx q[2];
rz(-2.6893534) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0248558) q[1];
sx q[1];
rz(-1.7953838) q[1];
sx q[1];
rz(-1.6889424) q[1];
x q[2];
rz(0.65770517) q[3];
sx q[3];
rz(-1.3240361) q[3];
sx q[3];
rz(-1.2031499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7449259) q[2];
sx q[2];
rz(-2.1325839) q[2];
sx q[2];
rz(2.1972556) q[2];
rz(2.5189597) q[3];
sx q[3];
rz(-0.98832744) q[3];
sx q[3];
rz(2.3545789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.054976376) q[0];
sx q[0];
rz(-2.3353307) q[0];
sx q[0];
rz(0.0042313519) q[0];
rz(-1.6731693) q[1];
sx q[1];
rz(-2.3783042) q[1];
sx q[1];
rz(-2.9671293) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3136351) q[0];
sx q[0];
rz(-2.9330755) q[0];
sx q[0];
rz(-2.0329518) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6930843) q[2];
sx q[2];
rz(-1.4786198) q[2];
sx q[2];
rz(-2.0192469) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0956438) q[1];
sx q[1];
rz(-2.2679331) q[1];
sx q[1];
rz(3.0554153) q[1];
rz(-pi) q[2];
rz(-1.8069491) q[3];
sx q[3];
rz(-1.5075348) q[3];
sx q[3];
rz(-2.0963653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8909495) q[2];
sx q[2];
rz(-3.0477016) q[2];
sx q[2];
rz(-2.6012361) q[2];
rz(3.12449) q[3];
sx q[3];
rz(-0.98186866) q[3];
sx q[3];
rz(2.4798149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0906202) q[0];
sx q[0];
rz(-0.80247387) q[0];
sx q[0];
rz(2.354994) q[0];
rz(-2.0358918) q[1];
sx q[1];
rz(-1.4812171) q[1];
sx q[1];
rz(0.22122637) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.008667058) q[0];
sx q[0];
rz(-1.2930388) q[0];
sx q[0];
rz(-2.3093669) q[0];
rz(2.0682797) q[2];
sx q[2];
rz(-0.19693298) q[2];
sx q[2];
rz(0.54468583) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7936642) q[1];
sx q[1];
rz(-1.3634487) q[1];
sx q[1];
rz(-2.6479032) q[1];
rz(0.28085162) q[3];
sx q[3];
rz(-0.70990585) q[3];
sx q[3];
rz(-0.8381697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1456387) q[2];
sx q[2];
rz(-0.24253878) q[2];
sx q[2];
rz(1.9560187) q[2];
rz(1.1061741) q[3];
sx q[3];
rz(-1.0217051) q[3];
sx q[3];
rz(-2.1574028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36158654) q[0];
sx q[0];
rz(-1.7842643) q[0];
sx q[0];
rz(0.7097882) q[0];
rz(-2.2478726) q[1];
sx q[1];
rz(-2.5137386) q[1];
sx q[1];
rz(0.46700221) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2115973) q[0];
sx q[0];
rz(-0.30147895) q[0];
sx q[0];
rz(-1.0903574) q[0];
x q[1];
rz(-2.1960718) q[2];
sx q[2];
rz(-1.3511124) q[2];
sx q[2];
rz(1.4131119) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.87709204) q[1];
sx q[1];
rz(-0.24152006) q[1];
sx q[1];
rz(-1.4143553) q[1];
x q[2];
rz(0.8740467) q[3];
sx q[3];
rz(-2.120904) q[3];
sx q[3];
rz(0.77220213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.799377) q[2];
sx q[2];
rz(-2.2365384) q[2];
sx q[2];
rz(-3.0055935) q[2];
rz(0.32986766) q[3];
sx q[3];
rz(-1.0672652) q[3];
sx q[3];
rz(2.8444667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3099404) q[0];
sx q[0];
rz(-1.1158054) q[0];
sx q[0];
rz(0.62624395) q[0];
rz(2.2533439) q[1];
sx q[1];
rz(-1.8969957) q[1];
sx q[1];
rz(-3.0096164) q[1];
rz(0.49909231) q[2];
sx q[2];
rz(-0.44514984) q[2];
sx q[2];
rz(0.48322074) q[2];
rz(-0.47223623) q[3];
sx q[3];
rz(-0.78651169) q[3];
sx q[3];
rz(3.095918) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
