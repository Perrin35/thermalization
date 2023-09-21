OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6991601) q[0];
sx q[0];
rz(-1.7572829) q[0];
sx q[0];
rz(1.260489) q[0];
rz(2.1029544) q[1];
sx q[1];
rz(-1.3488052) q[1];
sx q[1];
rz(-2.2178712) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1169352) q[0];
sx q[0];
rz(-0.57514656) q[0];
sx q[0];
rz(-0.92098178) q[0];
rz(1.107723) q[2];
sx q[2];
rz(-2.827364) q[2];
sx q[2];
rz(-2.3067834) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.87286283) q[1];
sx q[1];
rz(-1.9951207) q[1];
sx q[1];
rz(0.55117589) q[1];
rz(-pi) q[2];
rz(-1.9507017) q[3];
sx q[3];
rz(-2.0348747) q[3];
sx q[3];
rz(2.3604148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.91360056) q[2];
sx q[2];
rz(-1.2486518) q[2];
sx q[2];
rz(2.9795734) q[2];
rz(0.93531936) q[3];
sx q[3];
rz(-2.155442) q[3];
sx q[3];
rz(2.4285765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0733923) q[0];
sx q[0];
rz(-0.22664264) q[0];
sx q[0];
rz(-1.9447928) q[0];
rz(0.67990047) q[1];
sx q[1];
rz(-0.49566832) q[1];
sx q[1];
rz(-1.686036) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.992422) q[0];
sx q[0];
rz(-2.143321) q[0];
sx q[0];
rz(-0.19897977) q[0];
rz(-0.20632867) q[2];
sx q[2];
rz(-0.38197877) q[2];
sx q[2];
rz(-1.0155592) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5160408) q[1];
sx q[1];
rz(-1.5843154) q[1];
sx q[1];
rz(-1.0479755) q[1];
x q[2];
rz(3.0807642) q[3];
sx q[3];
rz(-0.41234499) q[3];
sx q[3];
rz(2.6044248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.42852795) q[2];
sx q[2];
rz(-1.45168) q[2];
sx q[2];
rz(-1.3519752) q[2];
rz(-0.18243608) q[3];
sx q[3];
rz(-0.97674102) q[3];
sx q[3];
rz(0.3119719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75333726) q[0];
sx q[0];
rz(-0.68080807) q[0];
sx q[0];
rz(0.80048168) q[0];
rz(0.02877409) q[1];
sx q[1];
rz(-1.0556227) q[1];
sx q[1];
rz(1.9690537) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7045672) q[0];
sx q[0];
rz(-0.66172681) q[0];
sx q[0];
rz(2.6252803) q[0];
rz(2.3861888) q[2];
sx q[2];
rz(-1.8131776) q[2];
sx q[2];
rz(2.541045) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2492003) q[1];
sx q[1];
rz(-1.3274267) q[1];
sx q[1];
rz(-2.1057486) q[1];
x q[2];
rz(1.3451194) q[3];
sx q[3];
rz(-2.884622) q[3];
sx q[3];
rz(-0.32303177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0671493) q[2];
sx q[2];
rz(-1.477244) q[2];
sx q[2];
rz(2.2303936) q[2];
rz(-0.95101142) q[3];
sx q[3];
rz(-2.337303) q[3];
sx q[3];
rz(-0.89200154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(0.38055414) q[0];
sx q[0];
rz(-3.0111713) q[0];
sx q[0];
rz(-0.12810853) q[0];
rz(-0.076106636) q[1];
sx q[1];
rz(-1.9271306) q[1];
sx q[1];
rz(0.52350837) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2142221) q[0];
sx q[0];
rz(-0.71338755) q[0];
sx q[0];
rz(-0.58332304) q[0];
rz(0.33004327) q[2];
sx q[2];
rz(-1.4903307) q[2];
sx q[2];
rz(-0.45804322) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.88672968) q[1];
sx q[1];
rz(-2.2543395) q[1];
sx q[1];
rz(0.32943326) q[1];
rz(1.3685437) q[3];
sx q[3];
rz(-2.5407255) q[3];
sx q[3];
rz(0.7522538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6161502) q[2];
sx q[2];
rz(-1.5443065) q[2];
sx q[2];
rz(0.564044) q[2];
rz(-2.8530252) q[3];
sx q[3];
rz(-0.42268649) q[3];
sx q[3];
rz(2.585876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48150912) q[0];
sx q[0];
rz(-2.4531589) q[0];
sx q[0];
rz(1.4915285) q[0];
rz(0.87961698) q[1];
sx q[1];
rz(-1.8938226) q[1];
sx q[1];
rz(0.99194828) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46734738) q[0];
sx q[0];
rz(-1.464198) q[0];
sx q[0];
rz(0.17515134) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2679342) q[2];
sx q[2];
rz(-1.5622683) q[2];
sx q[2];
rz(-2.0451343) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3757513) q[1];
sx q[1];
rz(-1.6110794) q[1];
sx q[1];
rz(-2.8326616) q[1];
x q[2];
rz(-1.7080073) q[3];
sx q[3];
rz(-2.3793594) q[3];
sx q[3];
rz(-1.8828132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0118959) q[2];
sx q[2];
rz(-0.36964881) q[2];
sx q[2];
rz(2.8707855) q[2];
rz(-0.21823847) q[3];
sx q[3];
rz(-1.821358) q[3];
sx q[3];
rz(2.9158084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72702423) q[0];
sx q[0];
rz(-2.4232061) q[0];
sx q[0];
rz(-1.7927992) q[0];
rz(-0.38189608) q[1];
sx q[1];
rz(-0.31612879) q[1];
sx q[1];
rz(1.7165002) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8558559) q[0];
sx q[0];
rz(-2.5897313) q[0];
sx q[0];
rz(-2.7049271) q[0];
x q[1];
rz(1.8259949) q[2];
sx q[2];
rz(-1.4824502) q[2];
sx q[2];
rz(-1.3510072) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7879494) q[1];
sx q[1];
rz(-1.1218346) q[1];
sx q[1];
rz(0.073972703) q[1];
x q[2];
rz(-0.53374966) q[3];
sx q[3];
rz(-1.0373877) q[3];
sx q[3];
rz(0.48983869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0075334) q[2];
sx q[2];
rz(-0.39742658) q[2];
sx q[2];
rz(2.5777204) q[2];
rz(2.9610736) q[3];
sx q[3];
rz(-1.518395) q[3];
sx q[3];
rz(-2.738651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5086223) q[0];
sx q[0];
rz(-0.16212012) q[0];
sx q[0];
rz(-0.41931835) q[0];
rz(-1.58889) q[1];
sx q[1];
rz(-1.2607375) q[1];
sx q[1];
rz(-2.3197876) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99188995) q[0];
sx q[0];
rz(-2.2737962) q[0];
sx q[0];
rz(0.71233149) q[0];
x q[1];
rz(-1.0242545) q[2];
sx q[2];
rz(-2.3356236) q[2];
sx q[2];
rz(0.96217996) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.070548363) q[1];
sx q[1];
rz(-0.76862915) q[1];
sx q[1];
rz(2.8119836) q[1];
rz(2.5542198) q[3];
sx q[3];
rz(-1.8837187) q[3];
sx q[3];
rz(-1.8250993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8043148) q[2];
sx q[2];
rz(-2.3874805) q[2];
sx q[2];
rz(-0.24469963) q[2];
rz(0.129536) q[3];
sx q[3];
rz(-1.9774388) q[3];
sx q[3];
rz(-1.6285508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41641763) q[0];
sx q[0];
rz(-0.019151909) q[0];
sx q[0];
rz(0.82292557) q[0];
rz(-0.30934632) q[1];
sx q[1];
rz(-1.3920709) q[1];
sx q[1];
rz(-1.3051422) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6179498) q[0];
sx q[0];
rz(-0.68730132) q[0];
sx q[0];
rz(-0.53074093) q[0];
x q[1];
rz(2.7140076) q[2];
sx q[2];
rz(-2.3901849) q[2];
sx q[2];
rz(0.69185711) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0481938) q[1];
sx q[1];
rz(-1.9462002) q[1];
sx q[1];
rz(-2.4673389) q[1];
x q[2];
rz(1.5042138) q[3];
sx q[3];
rz(-0.3399907) q[3];
sx q[3];
rz(0.41223994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1398853) q[2];
sx q[2];
rz(-1.7557764) q[2];
sx q[2];
rz(-1.6513599) q[2];
rz(-1.0772609) q[3];
sx q[3];
rz(-2.1765985) q[3];
sx q[3];
rz(0.13154496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7548783) q[0];
sx q[0];
rz(-1.2972378) q[0];
sx q[0];
rz(-2.819678) q[0];
rz(-1.5362668) q[1];
sx q[1];
rz(-1.221311) q[1];
sx q[1];
rz(2.4386491) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1212595) q[0];
sx q[0];
rz(-0.54176211) q[0];
sx q[0];
rz(-0.74777491) q[0];
rz(-pi) q[1];
rz(0.039673294) q[2];
sx q[2];
rz(-2.0797605) q[2];
sx q[2];
rz(-0.85658011) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.798588) q[1];
sx q[1];
rz(-1.4061905) q[1];
sx q[1];
rz(-2.1144457) q[1];
x q[2];
rz(1.1144936) q[3];
sx q[3];
rz(-0.66702402) q[3];
sx q[3];
rz(-1.4240571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.90074173) q[2];
sx q[2];
rz(-2.8618331) q[2];
sx q[2];
rz(-1.3396324) q[2];
rz(2.83589) q[3];
sx q[3];
rz(-1.327508) q[3];
sx q[3];
rz(1.8113177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9777578) q[0];
sx q[0];
rz(-0.75755388) q[0];
sx q[0];
rz(1.9158069) q[0];
rz(2.2380791) q[1];
sx q[1];
rz(-0.61360306) q[1];
sx q[1];
rz(2.6729565) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9194473) q[0];
sx q[0];
rz(-1.7788017) q[0];
sx q[0];
rz(-0.39214765) q[0];
rz(-pi) q[1];
rz(2.2655728) q[2];
sx q[2];
rz(-0.48601905) q[2];
sx q[2];
rz(1.8234058) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.84506449) q[1];
sx q[1];
rz(-1.5025286) q[1];
sx q[1];
rz(-1.8370085) q[1];
rz(-pi) q[2];
rz(-2.7077984) q[3];
sx q[3];
rz(-1.9926096) q[3];
sx q[3];
rz(-1.0292366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.48352155) q[2];
sx q[2];
rz(-1.8489685) q[2];
sx q[2];
rz(1.1432077) q[2];
rz(-3.0269567) q[3];
sx q[3];
rz(-0.95364037) q[3];
sx q[3];
rz(1.5293998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6464012) q[0];
sx q[0];
rz(-1.1947182) q[0];
sx q[0];
rz(2.4583046) q[0];
rz(2.519683) q[1];
sx q[1];
rz(-1.6786631) q[1];
sx q[1];
rz(2.8181029) q[1];
rz(-2.3399578) q[2];
sx q[2];
rz(-0.87052204) q[2];
sx q[2];
rz(-1.082765) q[2];
rz(-2.2223496) q[3];
sx q[3];
rz(-1.6332492) q[3];
sx q[3];
rz(1.9132683) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
