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
rz(0.723847) q[0];
sx q[0];
rz(-1.2011733) q[0];
sx q[0];
rz(0.49402753) q[0];
rz(1.5552893) q[1];
sx q[1];
rz(2.2145693) q[1];
sx q[1];
rz(8.5742843) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14334003) q[0];
sx q[0];
rz(-1.7368642) q[0];
sx q[0];
rz(-2.2448426) q[0];
rz(-0.65324983) q[2];
sx q[2];
rz(-1.946665) q[2];
sx q[2];
rz(2.1819161) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9697992) q[1];
sx q[1];
rz(-1.3487019) q[1];
sx q[1];
rz(-2.8104258) q[1];
rz(0.36601074) q[3];
sx q[3];
rz(-1.6386392) q[3];
sx q[3];
rz(-2.3930166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6583307) q[2];
sx q[2];
rz(-0.61507812) q[2];
sx q[2];
rz(-0.94493803) q[2];
rz(-3.1350709) q[3];
sx q[3];
rz(-2.3801453) q[3];
sx q[3];
rz(0.25750461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0354804) q[0];
sx q[0];
rz(-2.1529614) q[0];
sx q[0];
rz(2.535787) q[0];
rz(-2.2094191) q[1];
sx q[1];
rz(-1.7180387) q[1];
sx q[1];
rz(0.1056284) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52077928) q[0];
sx q[0];
rz(-0.70292369) q[0];
sx q[0];
rz(2.9018794) q[0];
rz(-pi) q[1];
rz(-2.1828853) q[2];
sx q[2];
rz(-2.5664133) q[2];
sx q[2];
rz(1.5633595) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2534598) q[1];
sx q[1];
rz(-1.2769298) q[1];
sx q[1];
rz(-0.47305686) q[1];
rz(-pi) q[2];
rz(0.83189957) q[3];
sx q[3];
rz(-2.6100592) q[3];
sx q[3];
rz(-3.1011776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.85670829) q[2];
sx q[2];
rz(-1.7785037) q[2];
sx q[2];
rz(1.523783) q[2];
rz(-0.81426042) q[3];
sx q[3];
rz(-1.7819449) q[3];
sx q[3];
rz(0.8849357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.043269) q[0];
sx q[0];
rz(-0.79541484) q[0];
sx q[0];
rz(-1.5765618) q[0];
rz(-0.99524975) q[1];
sx q[1];
rz(-2.1627656) q[1];
sx q[1];
rz(2.3547122) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3580619) q[0];
sx q[0];
rz(-1.770505) q[0];
sx q[0];
rz(0.13808967) q[0];
rz(-0.34721513) q[2];
sx q[2];
rz(-1.6499398) q[2];
sx q[2];
rz(1.7007507) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7461024) q[1];
sx q[1];
rz(-1.2416035) q[1];
sx q[1];
rz(-1.9026865) q[1];
x q[2];
rz(-1.5110817) q[3];
sx q[3];
rz(-1.6303326) q[3];
sx q[3];
rz(1.8465896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.83871049) q[2];
sx q[2];
rz(-1.4883214) q[2];
sx q[2];
rz(-0.6380471) q[2];
rz(-2.6436515) q[3];
sx q[3];
rz(-2.0697856) q[3];
sx q[3];
rz(1.0897442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8968673) q[0];
sx q[0];
rz(-1.7843972) q[0];
sx q[0];
rz(-3.0778399) q[0];
rz(1.6925192) q[1];
sx q[1];
rz(-1.8359102) q[1];
sx q[1];
rz(1.3099028) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1242715) q[0];
sx q[0];
rz(-1.4910264) q[0];
sx q[0];
rz(1.5333453) q[0];
rz(2.9430423) q[2];
sx q[2];
rz(-1.3891925) q[2];
sx q[2];
rz(1.1491733) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4938072) q[1];
sx q[1];
rz(-1.6946548) q[1];
sx q[1];
rz(2.3954443) q[1];
rz(0.88366644) q[3];
sx q[3];
rz(-1.851463) q[3];
sx q[3];
rz(1.3786045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0706851) q[2];
sx q[2];
rz(-2.0347774) q[2];
sx q[2];
rz(3.0359388) q[2];
rz(-1.9581155) q[3];
sx q[3];
rz(-2.2453997) q[3];
sx q[3];
rz(-2.4699874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79528177) q[0];
sx q[0];
rz(-2.2310937) q[0];
sx q[0];
rz(-1.7972535) q[0];
rz(2.4686939) q[1];
sx q[1];
rz(-1.3105323) q[1];
sx q[1];
rz(0.013462822) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5199842) q[0];
sx q[0];
rz(-2.2604687) q[0];
sx q[0];
rz(-3.1075928) q[0];
x q[1];
rz(1.156594) q[2];
sx q[2];
rz(-1.557781) q[2];
sx q[2];
rz(-0.77337298) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.33121651) q[1];
sx q[1];
rz(-1.5840285) q[1];
sx q[1];
rz(0.41535901) q[1];
x q[2];
rz(2.0125403) q[3];
sx q[3];
rz(-0.75296445) q[3];
sx q[3];
rz(-0.66679614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5992735) q[2];
sx q[2];
rz(-0.66212526) q[2];
sx q[2];
rz(-1.8947961) q[2];
rz(3.0689734) q[3];
sx q[3];
rz(-2.9473372) q[3];
sx q[3];
rz(0.3092002) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4319864) q[0];
sx q[0];
rz(-1.1259587) q[0];
sx q[0];
rz(3.0288938) q[0];
rz(-0.20507774) q[1];
sx q[1];
rz(-1.3321184) q[1];
sx q[1];
rz(2.233706) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68112198) q[0];
sx q[0];
rz(-2.2325071) q[0];
sx q[0];
rz(0.76028334) q[0];
x q[1];
rz(-2.5656469) q[2];
sx q[2];
rz(-0.69789808) q[2];
sx q[2];
rz(-2.9550748) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.19872936) q[1];
sx q[1];
rz(-0.94063771) q[1];
sx q[1];
rz(-1.6840001) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.6648074) q[3];
sx q[3];
rz(-1.4012573) q[3];
sx q[3];
rz(-2.4162718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1418566) q[2];
sx q[2];
rz(-0.33679589) q[2];
sx q[2];
rz(-3.0736308) q[2];
rz(-1.9939907) q[3];
sx q[3];
rz(-1.2679029) q[3];
sx q[3];
rz(-1.8806774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1934018) q[0];
sx q[0];
rz(-1.3113439) q[0];
sx q[0];
rz(2.6780658) q[0];
rz(2.9947128) q[1];
sx q[1];
rz(-1.905922) q[1];
sx q[1];
rz(-0.99259496) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50336058) q[0];
sx q[0];
rz(-1.322317) q[0];
sx q[0];
rz(-0.72297024) q[0];
x q[1];
rz(-1.3456773) q[2];
sx q[2];
rz(-0.29307191) q[2];
sx q[2];
rz(2.4973283) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.22659644) q[1];
sx q[1];
rz(-0.8454171) q[1];
sx q[1];
rz(-2.5195049) q[1];
rz(-pi) q[2];
rz(2.0247884) q[3];
sx q[3];
rz(-0.89056784) q[3];
sx q[3];
rz(2.6770563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.85840449) q[2];
sx q[2];
rz(-1.7121199) q[2];
sx q[2];
rz(2.6962213) q[2];
rz(1.3238268) q[3];
sx q[3];
rz(-1.0704853) q[3];
sx q[3];
rz(1.0858735) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.104798) q[0];
sx q[0];
rz(-0.0093655149) q[0];
sx q[0];
rz(2.985756) q[0];
rz(-1.2771295) q[1];
sx q[1];
rz(-1.3469478) q[1];
sx q[1];
rz(-0.015930463) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6496507) q[0];
sx q[0];
rz(-1.8109522) q[0];
sx q[0];
rz(-2.35046) q[0];
x q[1];
rz(-0.35959776) q[2];
sx q[2];
rz(-2.3162875) q[2];
sx q[2];
rz(1.8458837) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5585352) q[1];
sx q[1];
rz(-1.0638405) q[1];
sx q[1];
rz(0.81113775) q[1];
x q[2];
rz(2.5032477) q[3];
sx q[3];
rz(-1.9192291) q[3];
sx q[3];
rz(1.750647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5126123) q[2];
sx q[2];
rz(-3.0292558) q[2];
sx q[2];
rz(-3.0158499) q[2];
rz(-2.1851152) q[3];
sx q[3];
rz(-1.4230909) q[3];
sx q[3];
rz(2.8023348) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9823031) q[0];
sx q[0];
rz(-0.24022261) q[0];
sx q[0];
rz(-2.6851728) q[0];
rz(-1.5688815) q[1];
sx q[1];
rz(-1.0036889) q[1];
sx q[1];
rz(2.454954) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4149041) q[0];
sx q[0];
rz(-0.81302887) q[0];
sx q[0];
rz(0.70720478) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.94273156) q[2];
sx q[2];
rz(-0.21272993) q[2];
sx q[2];
rz(0.63791927) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2575779) q[1];
sx q[1];
rz(-2.9513689) q[1];
sx q[1];
rz(0.68489267) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7620662) q[3];
sx q[3];
rz(-2.9538493) q[3];
sx q[3];
rz(1.2953341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.75866428) q[2];
sx q[2];
rz(-1.1229346) q[2];
sx q[2];
rz(-2.0056966) q[2];
rz(0.9797594) q[3];
sx q[3];
rz(-1.8085248) q[3];
sx q[3];
rz(-2.4433344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6887688) q[0];
sx q[0];
rz(-1.8166421) q[0];
sx q[0];
rz(2.685637) q[0];
rz(-0.042757209) q[1];
sx q[1];
rz(-1.9330934) q[1];
sx q[1];
rz(-2.4483689) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7133378) q[0];
sx q[0];
rz(-1.1866259) q[0];
sx q[0];
rz(-3.0349656) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.40463368) q[2];
sx q[2];
rz(-1.9453959) q[2];
sx q[2];
rz(-0.43231479) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9849423) q[1];
sx q[1];
rz(-1.4727815) q[1];
sx q[1];
rz(-3.1344862) q[1];
rz(-1.8878804) q[3];
sx q[3];
rz(-1.1389995) q[3];
sx q[3];
rz(-2.3136239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.91528714) q[2];
sx q[2];
rz(-1.8733571) q[2];
sx q[2];
rz(2.477296) q[2];
rz(-1.9281467) q[3];
sx q[3];
rz(-2.4197141) q[3];
sx q[3];
rz(-2.2693995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6912457) q[0];
sx q[0];
rz(-2.3007614) q[0];
sx q[0];
rz(1.7578516) q[0];
rz(1.6830403) q[1];
sx q[1];
rz(-2.8589307) q[1];
sx q[1];
rz(0.9160441) q[1];
rz(3.0677879) q[2];
sx q[2];
rz(-2.7701785) q[2];
sx q[2];
rz(-2.2899173) q[2];
rz(1.0091598) q[3];
sx q[3];
rz(-2.1957993) q[3];
sx q[3];
rz(-1.5002863) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
