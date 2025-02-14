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
rz(-1.0280161) q[1];
sx q[1];
rz(2.4457959) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4980846) q[0];
sx q[0];
rz(-1.9970446) q[0];
sx q[0];
rz(0.059498992) q[0];
rz(-pi) q[1];
rz(0.037362413) q[2];
sx q[2];
rz(-1.6184095) q[2];
sx q[2];
rz(-1.7689153) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.6219139) q[1];
sx q[1];
rz(-1.7413472) q[1];
sx q[1];
rz(2.7615553) q[1];
rz(-pi) q[2];
rz(0.94732802) q[3];
sx q[3];
rz(-0.58637324) q[3];
sx q[3];
rz(2.9602667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9838788) q[2];
sx q[2];
rz(-2.3215051) q[2];
sx q[2];
rz(0.90041089) q[2];
rz(2.3484717) q[3];
sx q[3];
rz(-1.6499062) q[3];
sx q[3];
rz(-0.69860631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6124436) q[0];
sx q[0];
rz(-1.1684893) q[0];
sx q[0];
rz(-0.68552619) q[0];
rz(2.7087052) q[1];
sx q[1];
rz(-1.2665749) q[1];
sx q[1];
rz(-1.8276851) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9627808) q[0];
sx q[0];
rz(-1.0355486) q[0];
sx q[0];
rz(0.58387941) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.45998145) q[2];
sx q[2];
rz(-0.88085266) q[2];
sx q[2];
rz(-1.8073818) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5320783) q[1];
sx q[1];
rz(-1.1066965) q[1];
sx q[1];
rz(1.0232693) q[1];
rz(0.55696662) q[3];
sx q[3];
rz(-1.3560221) q[3];
sx q[3];
rz(-2.1208515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1757397) q[2];
sx q[2];
rz(-2.7312835) q[2];
sx q[2];
rz(-2.8941594) q[2];
rz(-3.0257709) q[3];
sx q[3];
rz(-1.7146866) q[3];
sx q[3];
rz(-2.5628832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.9847357) q[0];
sx q[0];
rz(-0.37608376) q[0];
sx q[0];
rz(-1.0648741) q[0];
rz(0.6330601) q[1];
sx q[1];
rz(-1.9828911) q[1];
sx q[1];
rz(-0.13654581) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8234607) q[0];
sx q[0];
rz(-0.7366418) q[0];
sx q[0];
rz(3.1053057) q[0];
rz(-pi) q[1];
rz(2.2011287) q[2];
sx q[2];
rz(-2.5879938) q[2];
sx q[2];
rz(-0.74650899) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.28440839) q[1];
sx q[1];
rz(-0.84759241) q[1];
sx q[1];
rz(0.64874362) q[1];
rz(-1.1401922) q[3];
sx q[3];
rz(-1.1893442) q[3];
sx q[3];
rz(-2.7611707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.442405) q[2];
sx q[2];
rz(-1.2523315) q[2];
sx q[2];
rz(-1.5415972) q[2];
rz(-1.0934746) q[3];
sx q[3];
rz(-1.6396061) q[3];
sx q[3];
rz(0.34539616) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71444756) q[0];
sx q[0];
rz(-2.8801425) q[0];
sx q[0];
rz(-2.5647822) q[0];
rz(0.72751865) q[1];
sx q[1];
rz(-1.0914165) q[1];
sx q[1];
rz(2.2732546) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8936794) q[0];
sx q[0];
rz(-2.1931139) q[0];
sx q[0];
rz(1.076368) q[0];
rz(-2.1530354) q[2];
sx q[2];
rz(-2.0209656) q[2];
sx q[2];
rz(-0.52272138) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7817745) q[1];
sx q[1];
rz(-1.9137772) q[1];
sx q[1];
rz(-2.4510259) q[1];
rz(-pi) q[2];
rz(0.94566019) q[3];
sx q[3];
rz(-1.994761) q[3];
sx q[3];
rz(-1.5353919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.19520983) q[2];
sx q[2];
rz(-1.752744) q[2];
sx q[2];
rz(-0.69551224) q[2];
rz(2.9269311) q[3];
sx q[3];
rz(-1.5011468) q[3];
sx q[3];
rz(2.3282839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.5401841) q[0];
sx q[0];
rz(-2.5267127) q[0];
sx q[0];
rz(1.1001128) q[0];
rz(2.9310138) q[1];
sx q[1];
rz(-2.5081878) q[1];
sx q[1];
rz(1.4932154) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8029047) q[0];
sx q[0];
rz(-0.7006104) q[0];
sx q[0];
rz(2.3466097) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0080935) q[2];
sx q[2];
rz(-0.19607142) q[2];
sx q[2];
rz(-2.9626949) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.666387) q[1];
sx q[1];
rz(-2.0568486) q[1];
sx q[1];
rz(-2.9654825) q[1];
rz(-1.8214956) q[3];
sx q[3];
rz(-2.2226102) q[3];
sx q[3];
rz(-1.0221611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7275927) q[2];
sx q[2];
rz(-1.4294581) q[2];
sx q[2];
rz(2.2445938) q[2];
rz(-1.8619079) q[3];
sx q[3];
rz(-2.1309659) q[3];
sx q[3];
rz(-3.0943833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
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
rz(0.33011398) q[0];
rz(-0.64078981) q[1];
sx q[1];
rz(-1.7658486) q[1];
sx q[1];
rz(-2.1660588) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2084239) q[0];
sx q[0];
rz(-2.5481249) q[0];
sx q[0];
rz(0.10153254) q[0];
rz(2.2660146) q[2];
sx q[2];
rz(-2.7355614) q[2];
sx q[2];
rz(0.46689597) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9258049) q[1];
sx q[1];
rz(-0.5277718) q[1];
sx q[1];
rz(0.58234071) q[1];
rz(-1.9930494) q[3];
sx q[3];
rz(-1.7429757) q[3];
sx q[3];
rz(2.9957221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8581802) q[2];
sx q[2];
rz(-0.25412574) q[2];
sx q[2];
rz(0.88968366) q[2];
rz(2.5455918) q[3];
sx q[3];
rz(-1.069205) q[3];
sx q[3];
rz(-0.10665756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(3.0058873) q[0];
sx q[0];
rz(-1.1361253) q[0];
sx q[0];
rz(-1.1329875) q[0];
rz(0.63944447) q[1];
sx q[1];
rz(-1.0777799) q[1];
sx q[1];
rz(0.9674255) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4690034) q[0];
sx q[0];
rz(-0.53454195) q[0];
sx q[0];
rz(1.8683167) q[0];
rz(2.4339811) q[2];
sx q[2];
rz(-0.48383265) q[2];
sx q[2];
rz(2.8255723) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6270476) q[1];
sx q[1];
rz(-2.8882898) q[1];
sx q[1];
rz(2.6652418) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.65770517) q[3];
sx q[3];
rz(-1.3240361) q[3];
sx q[3];
rz(1.2031499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7449259) q[2];
sx q[2];
rz(-2.1325839) q[2];
sx q[2];
rz(0.94433707) q[2];
rz(0.62263292) q[3];
sx q[3];
rz(-0.98832744) q[3];
sx q[3];
rz(-2.3545789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.054976376) q[0];
sx q[0];
rz(-0.80626196) q[0];
sx q[0];
rz(-0.0042313519) q[0];
rz(1.6731693) q[1];
sx q[1];
rz(-0.76328841) q[1];
sx q[1];
rz(0.17446336) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28936548) q[0];
sx q[0];
rz(-1.4783637) q[0];
sx q[0];
rz(1.3836181) q[0];
x q[1];
rz(0.44850839) q[2];
sx q[2];
rz(-1.6629728) q[2];
sx q[2];
rz(1.1223457) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.045948898) q[1];
sx q[1];
rz(-2.2679331) q[1];
sx q[1];
rz(-0.086177372) q[1];
rz(-pi) q[2];
x q[2];
rz(0.065062398) q[3];
sx q[3];
rz(-1.8064677) q[3];
sx q[3];
rz(0.5103569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8909495) q[2];
sx q[2];
rz(-0.093891056) q[2];
sx q[2];
rz(-0.54035652) q[2];
rz(0.017102608) q[3];
sx q[3];
rz(-2.159724) q[3];
sx q[3];
rz(2.4798149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0906202) q[0];
sx q[0];
rz(-0.80247387) q[0];
sx q[0];
rz(-0.78659868) q[0];
rz(-1.1057009) q[1];
sx q[1];
rz(-1.4812171) q[1];
sx q[1];
rz(2.9203663) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3174789) q[0];
sx q[0];
rz(-2.2750018) q[0];
sx q[0];
rz(0.368035) q[0];
x q[1];
rz(1.7443666) q[2];
sx q[2];
rz(-1.4772869) q[2];
sx q[2];
rz(1.6261404) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7936642) q[1];
sx q[1];
rz(-1.778144) q[1];
sx q[1];
rz(2.6479032) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.28085162) q[3];
sx q[3];
rz(-2.4316868) q[3];
sx q[3];
rz(-0.8381697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9959539) q[2];
sx q[2];
rz(-2.8990539) q[2];
sx q[2];
rz(1.9560187) q[2];
rz(-2.0354185) q[3];
sx q[3];
rz(-2.1198876) q[3];
sx q[3];
rz(2.1574028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7800061) q[0];
sx q[0];
rz(-1.3573283) q[0];
sx q[0];
rz(-0.7097882) q[0];
rz(-2.2478726) q[1];
sx q[1];
rz(-0.62785405) q[1];
sx q[1];
rz(2.6745904) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7119962) q[0];
sx q[0];
rz(-1.8372559) q[0];
sx q[0];
rz(0.14273739) q[0];
rz(2.1960718) q[2];
sx q[2];
rz(-1.3511124) q[2];
sx q[2];
rz(-1.4131119) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2959173) q[1];
sx q[1];
rz(-1.6080699) q[1];
sx q[1];
rz(1.8094783) q[1];
x q[2];
rz(-0.8740467) q[3];
sx q[3];
rz(-2.120904) q[3];
sx q[3];
rz(2.3693905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.799377) q[2];
sx q[2];
rz(-2.2365384) q[2];
sx q[2];
rz(0.13599914) q[2];
rz(2.811725) q[3];
sx q[3];
rz(-1.0672652) q[3];
sx q[3];
rz(-2.8444667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
x q[2];
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
rz(2.7449079) q[2];
sx q[2];
rz(-1.7783782) q[2];
sx q[2];
rz(1.5967899) q[2];
rz(0.47223623) q[3];
sx q[3];
rz(-2.355081) q[3];
sx q[3];
rz(-0.045674617) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
