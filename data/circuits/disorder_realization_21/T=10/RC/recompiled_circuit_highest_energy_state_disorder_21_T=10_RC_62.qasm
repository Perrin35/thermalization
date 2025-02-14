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
rz(-2.8357808) q[0];
sx q[0];
rz(-0.6337136) q[0];
sx q[0];
rz(-0.60854882) q[0];
rz(2.4701056) q[1];
sx q[1];
rz(-0.66317135) q[1];
sx q[1];
rz(1.6821678) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43435758) q[0];
sx q[0];
rz(-0.57091129) q[0];
sx q[0];
rz(2.0152604) q[0];
rz(-pi) q[1];
rz(2.0107277) q[2];
sx q[2];
rz(-0.49060433) q[2];
sx q[2];
rz(3.0705796) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0406831) q[1];
sx q[1];
rz(-0.95129943) q[1];
sx q[1];
rz(0.49609523) q[1];
x q[2];
rz(-0.75992775) q[3];
sx q[3];
rz(-1.6960521) q[3];
sx q[3];
rz(-1.6881936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3964316) q[2];
sx q[2];
rz(-2.4691984) q[2];
sx q[2];
rz(-2.3154837) q[2];
rz(-0.36469001) q[3];
sx q[3];
rz(-1.5267905) q[3];
sx q[3];
rz(2.8472692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18225886) q[0];
sx q[0];
rz(-1.2428281) q[0];
sx q[0];
rz(-2.6958595) q[0];
rz(0.590473) q[1];
sx q[1];
rz(-1.0346552) q[1];
sx q[1];
rz(-2.7139434) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70156258) q[0];
sx q[0];
rz(-2.0272581) q[0];
sx q[0];
rz(-0.052340074) q[0];
x q[1];
rz(-0.56088367) q[2];
sx q[2];
rz(-2.2429401) q[2];
sx q[2];
rz(2.667556) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.94484896) q[1];
sx q[1];
rz(-1.3102753) q[1];
sx q[1];
rz(2.1900106) q[1];
x q[2];
rz(-3.1366051) q[3];
sx q[3];
rz(-0.79646275) q[3];
sx q[3];
rz(-2.4465268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6815971) q[2];
sx q[2];
rz(-0.57969105) q[2];
sx q[2];
rz(2.106529) q[2];
rz(1.5486108) q[3];
sx q[3];
rz(-0.17954738) q[3];
sx q[3];
rz(0.88919324) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0276133) q[0];
sx q[0];
rz(-3.1359105) q[0];
sx q[0];
rz(0.11153829) q[0];
rz(1.5882209) q[1];
sx q[1];
rz(-2.7370743) q[1];
sx q[1];
rz(2.9389971) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83765471) q[0];
sx q[0];
rz(-1.9662153) q[0];
sx q[0];
rz(0.14531345) q[0];
rz(-pi) q[1];
rz(-1.8677223) q[2];
sx q[2];
rz(-1.8477173) q[2];
sx q[2];
rz(0.94147791) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2059171) q[1];
sx q[1];
rz(-0.96921235) q[1];
sx q[1];
rz(-2.4755073) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5532365) q[3];
sx q[3];
rz(-0.6783456) q[3];
sx q[3];
rz(2.7916186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.69278875) q[2];
sx q[2];
rz(-2.0710129) q[2];
sx q[2];
rz(-0.98133522) q[2];
rz(2.8547309) q[3];
sx q[3];
rz(-2.562264) q[3];
sx q[3];
rz(0.19744344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2364748) q[0];
sx q[0];
rz(-0.050967) q[0];
sx q[0];
rz(-2.4531051) q[0];
rz(-0.14432898) q[1];
sx q[1];
rz(-1.9537484) q[1];
sx q[1];
rz(-2.3932638) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28689529) q[0];
sx q[0];
rz(-2.4208491) q[0];
sx q[0];
rz(-3.0740644) q[0];
x q[1];
rz(-2.2950254) q[2];
sx q[2];
rz(-2.7477816) q[2];
sx q[2];
rz(-1.9005601) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6627745) q[1];
sx q[1];
rz(-1.2726901) q[1];
sx q[1];
rz(1.7653905) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2714628) q[3];
sx q[3];
rz(-1.3324454) q[3];
sx q[3];
rz(2.7173661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.051108483) q[2];
sx q[2];
rz(-1.313831) q[2];
sx q[2];
rz(1.2468106) q[2];
rz(-3.0221853) q[3];
sx q[3];
rz(-1.9236919) q[3];
sx q[3];
rz(-2.8363805) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7684105) q[0];
sx q[0];
rz(-2.7550582) q[0];
sx q[0];
rz(-0.61491948) q[0];
rz(-0.83456314) q[1];
sx q[1];
rz(-1.7465697) q[1];
sx q[1];
rz(-0.83647299) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8744321) q[0];
sx q[0];
rz(-1.5800486) q[0];
sx q[0];
rz(-1.6051588) q[0];
rz(-pi) q[1];
rz(-0.98149379) q[2];
sx q[2];
rz(-0.98926614) q[2];
sx q[2];
rz(-2.6133693) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.78042049) q[1];
sx q[1];
rz(-2.5675312) q[1];
sx q[1];
rz(-1.7758382) q[1];
rz(-pi) q[2];
rz(2.3990666) q[3];
sx q[3];
rz(-2.1603185) q[3];
sx q[3];
rz(1.3878915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5323083) q[2];
sx q[2];
rz(-2.534635) q[2];
sx q[2];
rz(0.92030805) q[2];
rz(-2.3406384) q[3];
sx q[3];
rz(-0.52217537) q[3];
sx q[3];
rz(0.33867684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9698708) q[0];
sx q[0];
rz(-2.0750177) q[0];
sx q[0];
rz(-0.26611662) q[0];
rz(1.7099821) q[1];
sx q[1];
rz(-2.0224729) q[1];
sx q[1];
rz(2.5439579) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9775598) q[0];
sx q[0];
rz(-2.4983239) q[0];
sx q[0];
rz(-0.32001647) q[0];
x q[1];
rz(2.9605537) q[2];
sx q[2];
rz(-1.5597876) q[2];
sx q[2];
rz(0.77048466) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1823765) q[1];
sx q[1];
rz(-0.87494779) q[1];
sx q[1];
rz(1.6435502) q[1];
rz(-1.8915755) q[3];
sx q[3];
rz(-1.4238384) q[3];
sx q[3];
rz(-3.1023417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1282244) q[2];
sx q[2];
rz(-0.2672264) q[2];
sx q[2];
rz(-0.5989778) q[2];
rz(0.61601764) q[3];
sx q[3];
rz(-0.44121656) q[3];
sx q[3];
rz(-2.0986572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1004341) q[0];
sx q[0];
rz(-0.91223311) q[0];
sx q[0];
rz(2.8225733) q[0];
rz(3.0306385) q[1];
sx q[1];
rz(-1.7496505) q[1];
sx q[1];
rz(-2.9999733) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9931147) q[0];
sx q[0];
rz(-1.6671868) q[0];
sx q[0];
rz(0.68401159) q[0];
rz(-0.0044195375) q[2];
sx q[2];
rz(-1.0194155) q[2];
sx q[2];
rz(-1.0221635) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1846022) q[1];
sx q[1];
rz(-0.58035589) q[1];
sx q[1];
rz(2.7932211) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0504248) q[3];
sx q[3];
rz(-1.7014352) q[3];
sx q[3];
rz(-1.7294693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5854599) q[2];
sx q[2];
rz(-2.1161049) q[2];
sx q[2];
rz(-3.0312313) q[2];
rz(-0.50802556) q[3];
sx q[3];
rz(-3.09943) q[3];
sx q[3];
rz(1.9917816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0069649) q[0];
sx q[0];
rz(-2.4899794) q[0];
sx q[0];
rz(1.9223258) q[0];
rz(0.02267516) q[1];
sx q[1];
rz(-0.89212787) q[1];
sx q[1];
rz(2.5882744) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1071725) q[0];
sx q[0];
rz(-1.6593455) q[0];
sx q[0];
rz(-2.2293985) q[0];
x q[1];
rz(-2.7560541) q[2];
sx q[2];
rz(-1.0304759) q[2];
sx q[2];
rz(1.413942) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0670315) q[1];
sx q[1];
rz(-2.0843049) q[1];
sx q[1];
rz(1.1719955) q[1];
x q[2];
rz(-1.8692006) q[3];
sx q[3];
rz(-1.8074028) q[3];
sx q[3];
rz(0.92577705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.40552178) q[2];
sx q[2];
rz(-0.61554474) q[2];
sx q[2];
rz(1.7919398) q[2];
rz(0.060529709) q[3];
sx q[3];
rz(-0.90977257) q[3];
sx q[3];
rz(2.5133666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8796006) q[0];
sx q[0];
rz(-2.1475726) q[0];
sx q[0];
rz(-3.1349728) q[0];
rz(2.5026542) q[1];
sx q[1];
rz(-2.9220351) q[1];
sx q[1];
rz(-0.42983291) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.00365) q[0];
sx q[0];
rz(-1.9406422) q[0];
sx q[0];
rz(-1.957264) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0833932) q[2];
sx q[2];
rz(-0.48658961) q[2];
sx q[2];
rz(-1.8712107) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2158969) q[1];
sx q[1];
rz(-1.5491423) q[1];
sx q[1];
rz(-2.3587061) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5019526) q[3];
sx q[3];
rz(-0.9864583) q[3];
sx q[3];
rz(0.39019728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2868353) q[2];
sx q[2];
rz(-1.9660051) q[2];
sx q[2];
rz(0.01469928) q[2];
rz(-2.0059351) q[3];
sx q[3];
rz(-2.0292) q[3];
sx q[3];
rz(1.8358102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.056674615) q[0];
sx q[0];
rz(-0.25923964) q[0];
sx q[0];
rz(-1.8490476) q[0];
rz(-0.1560642) q[1];
sx q[1];
rz(-0.81764692) q[1];
sx q[1];
rz(0.83052105) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53180689) q[0];
sx q[0];
rz(-1.6702819) q[0];
sx q[0];
rz(-2.2983944) q[0];
x q[1];
rz(1.2654598) q[2];
sx q[2];
rz(-0.80080253) q[2];
sx q[2];
rz(2.3126471) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9869186) q[1];
sx q[1];
rz(-0.23122825) q[1];
sx q[1];
rz(1.926941) q[1];
rz(-1.3545348) q[3];
sx q[3];
rz(-1.2118846) q[3];
sx q[3];
rz(-2.2674136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.10892756) q[2];
sx q[2];
rz(-0.69643164) q[2];
sx q[2];
rz(-0.44309524) q[2];
rz(2.9923934) q[3];
sx q[3];
rz(-0.34330338) q[3];
sx q[3];
rz(-2.8086015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14015848) q[0];
sx q[0];
rz(-2.5509111) q[0];
sx q[0];
rz(2.1668707) q[0];
rz(-2.1433266) q[1];
sx q[1];
rz(-1.9010192) q[1];
sx q[1];
rz(-0.98767282) q[1];
rz(-1.3692296) q[2];
sx q[2];
rz(-1.907363) q[2];
sx q[2];
rz(3.0628352) q[2];
rz(2.0988437) q[3];
sx q[3];
rz(-0.73560148) q[3];
sx q[3];
rz(-0.92675496) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
