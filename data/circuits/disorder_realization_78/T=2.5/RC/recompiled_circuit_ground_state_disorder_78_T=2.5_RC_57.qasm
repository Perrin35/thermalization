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
rz(-1.7011473) q[1];
sx q[1];
rz(-2.1135766) q[1];
sx q[1];
rz(-2.4457959) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0442565) q[0];
sx q[0];
rz(-1.5166266) q[0];
sx q[0];
rz(-1.143881) q[0];
rz(-pi) q[1];
rz(0.90593018) q[2];
sx q[2];
rz(-0.060513721) q[2];
sx q[2];
rz(1.1031594) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2604044) q[1];
sx q[1];
rz(-1.9450448) q[1];
sx q[1];
rz(-1.3874235) q[1];
rz(-pi) q[2];
rz(2.7715922) q[3];
sx q[3];
rz(-1.1048855) q[3];
sx q[3];
rz(2.6107058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.15771389) q[2];
sx q[2];
rz(-2.3215051) q[2];
sx q[2];
rz(0.90041089) q[2];
rz(-2.3484717) q[3];
sx q[3];
rz(-1.6499062) q[3];
sx q[3];
rz(0.69860631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5291491) q[0];
sx q[0];
rz(-1.1684893) q[0];
sx q[0];
rz(-0.68552619) q[0];
rz(0.43288747) q[1];
sx q[1];
rz(-1.2665749) q[1];
sx q[1];
rz(-1.3139075) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.066931574) q[0];
sx q[0];
rz(-2.064813) q[0];
sx q[0];
rz(-2.1886909) q[0];
rz(-0.45998145) q[2];
sx q[2];
rz(-0.88085266) q[2];
sx q[2];
rz(-1.8073818) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8364468) q[1];
sx q[1];
rz(-1.0865819) q[1];
sx q[1];
rz(0.53026235) q[1];
x q[2];
rz(-2.584626) q[3];
sx q[3];
rz(-1.7855705) q[3];
sx q[3];
rz(2.1208515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.96585298) q[2];
sx q[2];
rz(-2.7312835) q[2];
sx q[2];
rz(0.24743323) q[2];
rz(-0.11582173) q[3];
sx q[3];
rz(-1.426906) q[3];
sx q[3];
rz(-2.5628832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1568569) q[0];
sx q[0];
rz(-0.37608376) q[0];
sx q[0];
rz(1.0648741) q[0];
rz(0.6330601) q[1];
sx q[1];
rz(-1.9828911) q[1];
sx q[1];
rz(-0.13654581) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9158123) q[0];
sx q[0];
rz(-1.5951711) q[0];
sx q[0];
rz(-0.73631411) q[0];
x q[1];
rz(2.7922379) q[2];
sx q[2];
rz(-1.1321448) q[2];
sx q[2];
rz(-3.1040526) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8571843) q[1];
sx q[1];
rz(-2.2940002) q[1];
sx q[1];
rz(2.492849) q[1];
rz(-pi) q[2];
rz(-1.1401922) q[3];
sx q[3];
rz(-1.1893442) q[3];
sx q[3];
rz(-2.7611707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.442405) q[2];
sx q[2];
rz(-1.2523315) q[2];
sx q[2];
rz(-1.5999954) q[2];
rz(-1.0934746) q[3];
sx q[3];
rz(-1.5019865) q[3];
sx q[3];
rz(-0.34539616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4271451) q[0];
sx q[0];
rz(-2.8801425) q[0];
sx q[0];
rz(-0.57681042) q[0];
rz(2.414074) q[1];
sx q[1];
rz(-1.0914165) q[1];
sx q[1];
rz(-2.2732546) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.627358) q[0];
sx q[0];
rz(-1.1749724) q[0];
sx q[0];
rz(-0.68383278) q[0];
x q[1];
rz(-0.98855726) q[2];
sx q[2];
rz(-2.0209656) q[2];
sx q[2];
rz(-2.6188713) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5440923) q[1];
sx q[1];
rz(-2.3832633) q[1];
sx q[1];
rz(2.6306398) q[1];
x q[2];
rz(2.2277539) q[3];
sx q[3];
rz(-0.73907494) q[3];
sx q[3];
rz(0.55348671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.19520983) q[2];
sx q[2];
rz(-1.752744) q[2];
sx q[2];
rz(-0.69551224) q[2];
rz(-0.21466151) q[3];
sx q[3];
rz(-1.5011468) q[3];
sx q[3];
rz(-0.8133088) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6014086) q[0];
sx q[0];
rz(-2.5267127) q[0];
sx q[0];
rz(-2.0414798) q[0];
rz(0.21057883) q[1];
sx q[1];
rz(-2.5081878) q[1];
sx q[1];
rz(1.6483773) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8755668) q[0];
sx q[0];
rz(-2.039209) q[0];
sx q[0];
rz(-1.0288971) q[0];
rz(3.0576747) q[2];
sx q[2];
rz(-1.7482107) q[2];
sx q[2];
rz(0.62366297) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3201218) q[1];
sx q[1];
rz(-1.7263328) q[1];
sx q[1];
rz(1.0782773) q[1];
x q[2];
rz(0.3143309) q[3];
sx q[3];
rz(-2.4498457) q[3];
sx q[3];
rz(-1.421613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7275927) q[2];
sx q[2];
rz(-1.4294581) q[2];
sx q[2];
rz(-2.2445938) q[2];
rz(-1.8619079) q[3];
sx q[3];
rz(-2.1309659) q[3];
sx q[3];
rz(0.047209386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1419264) q[0];
sx q[0];
rz(-0.62586313) q[0];
sx q[0];
rz(-2.8114787) q[0];
rz(-2.5008028) q[1];
sx q[1];
rz(-1.375744) q[1];
sx q[1];
rz(0.97553387) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9331688) q[0];
sx q[0];
rz(-0.59346775) q[0];
sx q[0];
rz(-0.10153254) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8896721) q[2];
sx q[2];
rz(-1.3150196) q[2];
sx q[2];
rz(-2.6915611) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.16224081) q[1];
sx q[1];
rz(-1.8514412) q[1];
sx q[1];
rz(2.6885217) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.18838547) q[3];
sx q[3];
rz(-1.9864161) q[3];
sx q[3];
rz(-1.3481026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2834125) q[2];
sx q[2];
rz(-0.25412574) q[2];
sx q[2];
rz(-0.88968366) q[2];
rz(-2.5455918) q[3];
sx q[3];
rz(-1.069205) q[3];
sx q[3];
rz(0.10665756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13570531) q[0];
sx q[0];
rz(-2.0054673) q[0];
sx q[0];
rz(1.1329875) q[0];
rz(-0.63944447) q[1];
sx q[1];
rz(-2.0638128) q[1];
sx q[1];
rz(0.9674255) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9854162) q[0];
sx q[0];
rz(-1.7207017) q[0];
sx q[0];
rz(2.085859) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7616639) q[2];
sx q[2];
rz(-1.2636145) q[2];
sx q[2];
rz(1.9028705) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1167369) q[1];
sx q[1];
rz(-1.7953838) q[1];
sx q[1];
rz(-1.6889424) q[1];
rz(2.7507325) q[3];
sx q[3];
rz(-0.69598143) q[3];
sx q[3];
rz(-0.67392193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3966668) q[2];
sx q[2];
rz(-2.1325839) q[2];
sx q[2];
rz(2.1972556) q[2];
rz(0.62263292) q[3];
sx q[3];
rz(-0.98832744) q[3];
sx q[3];
rz(0.78701377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0866163) q[0];
sx q[0];
rz(-0.80626196) q[0];
sx q[0];
rz(-3.1373613) q[0];
rz(1.4684234) q[1];
sx q[1];
rz(-0.76328841) q[1];
sx q[1];
rz(-0.17446336) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28936548) q[0];
sx q[0];
rz(-1.4783637) q[0];
sx q[0];
rz(-1.7579745) q[0];
x q[1];
rz(-2.9315591) q[2];
sx q[2];
rz(-2.6843417) q[2];
sx q[2];
rz(-0.63744545) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.17970615) q[1];
sx q[1];
rz(-0.70155662) q[1];
sx q[1];
rz(1.6732193) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3063844) q[3];
sx q[3];
rz(-0.2443265) q[3];
sx q[3];
rz(0.78247386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2506432) q[2];
sx q[2];
rz(-0.093891056) q[2];
sx q[2];
rz(-0.54035652) q[2];
rz(0.017102608) q[3];
sx q[3];
rz(-0.98186866) q[3];
sx q[3];
rz(-2.4798149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0509725) q[0];
sx q[0];
rz(-2.3391188) q[0];
sx q[0];
rz(-2.354994) q[0];
rz(-2.0358918) q[1];
sx q[1];
rz(-1.6603755) q[1];
sx q[1];
rz(2.9203663) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.008667058) q[0];
sx q[0];
rz(-1.2930388) q[0];
sx q[0];
rz(0.8322258) q[0];
x q[1];
rz(-0.094927364) q[2];
sx q[2];
rz(-1.397992) q[2];
sx q[2];
rz(-3.1026187) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0290585) q[1];
sx q[1];
rz(-1.0885982) q[1];
sx q[1];
rz(-1.8052989) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8046326) q[3];
sx q[3];
rz(-2.2474996) q[3];
sx q[3];
rz(-1.2016344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9959539) q[2];
sx q[2];
rz(-0.24253878) q[2];
sx q[2];
rz(1.185574) q[2];
rz(-2.0354185) q[3];
sx q[3];
rz(-2.1198876) q[3];
sx q[3];
rz(2.1574028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7800061) q[0];
sx q[0];
rz(-1.7842643) q[0];
sx q[0];
rz(-2.4318045) q[0];
rz(0.89372006) q[1];
sx q[1];
rz(-2.5137386) q[1];
sx q[1];
rz(0.46700221) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9625682) q[0];
sx q[0];
rz(-1.4331289) q[0];
sx q[0];
rz(1.839864) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2063575) q[2];
sx q[2];
rz(-0.65783892) q[2];
sx q[2];
rz(0.45085622) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.87709204) q[1];
sx q[1];
rz(-2.9000726) q[1];
sx q[1];
rz(1.7272374) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8080838) q[3];
sx q[3];
rz(-2.2833725) q[3];
sx q[3];
rz(2.9014362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.799377) q[2];
sx q[2];
rz(-0.90505427) q[2];
sx q[2];
rz(0.13599914) q[2];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
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
rz(-2.2533439) q[1];
sx q[1];
rz(-1.2445969) q[1];
sx q[1];
rz(0.1319763) q[1];
rz(1.7952948) q[2];
sx q[2];
rz(-1.9585051) q[2];
sx q[2];
rz(-0.060123882) q[2];
rz(2.6693564) q[3];
sx q[3];
rz(-0.78651169) q[3];
sx q[3];
rz(3.095918) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
