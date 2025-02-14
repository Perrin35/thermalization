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
rz(2.5330438) q[0];
rz(-0.67148709) q[1];
sx q[1];
rz(3.804764) q[1];
sx q[1];
rz(10.884203) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7072351) q[0];
sx q[0];
rz(-2.5706814) q[0];
sx q[0];
rz(2.0152604) q[0];
rz(-2.9179108) q[2];
sx q[2];
rz(-2.0111901) q[2];
sx q[2];
rz(-0.5612095) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.83460036) q[1];
sx q[1];
rz(-1.9687593) q[1];
sx q[1];
rz(0.88944056) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3816649) q[3];
sx q[3];
rz(-1.6960521) q[3];
sx q[3];
rz(1.6881936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7451611) q[2];
sx q[2];
rz(-2.4691984) q[2];
sx q[2];
rz(-0.82610899) q[2];
rz(0.36469001) q[3];
sx q[3];
rz(-1.5267905) q[3];
sx q[3];
rz(-2.8472692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18225886) q[0];
sx q[0];
rz(-1.2428281) q[0];
sx q[0];
rz(-2.6958595) q[0];
rz(2.5511197) q[1];
sx q[1];
rz(-2.1069374) q[1];
sx q[1];
rz(-2.7139434) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4400301) q[0];
sx q[0];
rz(-1.1143346) q[0];
sx q[0];
rz(-3.0892526) q[0];
x q[1];
rz(-0.81645963) q[2];
sx q[2];
rz(-2.0000946) q[2];
sx q[2];
rz(-2.4176837) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.97274894) q[1];
sx q[1];
rz(-2.4764937) q[1];
sx q[1];
rz(1.14023) q[1];
rz(-1.5656972) q[3];
sx q[3];
rz(-0.77434629) q[3];
sx q[3];
rz(2.4536595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.4599956) q[2];
sx q[2];
rz(-2.5619016) q[2];
sx q[2];
rz(1.0350636) q[2];
rz(1.5929818) q[3];
sx q[3];
rz(-0.17954738) q[3];
sx q[3];
rz(2.2523994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0276133) q[0];
sx q[0];
rz(-3.1359105) q[0];
sx q[0];
rz(0.11153829) q[0];
rz(1.5533718) q[1];
sx q[1];
rz(-0.40451834) q[1];
sx q[1];
rz(-0.20259556) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.667019) q[0];
sx q[0];
rz(-2.7216421) q[0];
sx q[0];
rz(1.2368747) q[0];
rz(-pi) q[1];
x q[1];
rz(0.28891383) q[2];
sx q[2];
rz(-1.8560858) q[2];
sx q[2];
rz(2.4288175) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8823653) q[1];
sx q[1];
rz(-2.276032) q[1];
sx q[1];
rz(-2.3036876) q[1];
x q[2];
rz(-1.9913909) q[3];
sx q[3];
rz(-1.0216102) q[3];
sx q[3];
rz(-0.35850108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4488039) q[2];
sx q[2];
rz(-2.0710129) q[2];
sx q[2];
rz(0.98133522) q[2];
rz(-2.8547309) q[3];
sx q[3];
rz(-2.562264) q[3];
sx q[3];
rz(-0.19744344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2364748) q[0];
sx q[0];
rz(-0.050967) q[0];
sx q[0];
rz(0.68848759) q[0];
rz(2.9972637) q[1];
sx q[1];
rz(-1.1878443) q[1];
sx q[1];
rz(2.3932638) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2331325) q[0];
sx q[0];
rz(-1.5262506) q[0];
sx q[0];
rz(0.71961232) q[0];
rz(-pi) q[1];
x q[1];
rz(0.84656723) q[2];
sx q[2];
rz(-0.39381105) q[2];
sx q[2];
rz(1.9005601) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.47881815) q[1];
sx q[1];
rz(-1.8689026) q[1];
sx q[1];
rz(-1.7653905) q[1];
rz(-pi) q[2];
rz(2.2714628) q[3];
sx q[3];
rz(-1.8091473) q[3];
sx q[3];
rz(-2.7173661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.051108483) q[2];
sx q[2];
rz(-1.8277617) q[2];
sx q[2];
rz(-1.8947821) q[2];
rz(-0.11940739) q[3];
sx q[3];
rz(-1.9236919) q[3];
sx q[3];
rz(2.8363805) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3731821) q[0];
sx q[0];
rz(-0.38653448) q[0];
sx q[0];
rz(2.5266732) q[0];
rz(0.83456314) q[1];
sx q[1];
rz(-1.7465697) q[1];
sx q[1];
rz(0.83647299) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.040720721) q[0];
sx q[0];
rz(-0.03558579) q[0];
sx q[0];
rz(-1.8338704) q[0];
rz(-pi) q[1];
x q[1];
rz(0.98149379) q[2];
sx q[2];
rz(-0.98926614) q[2];
sx q[2];
rz(2.6133693) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.96326085) q[1];
sx q[1];
rz(-1.6815917) q[1];
sx q[1];
rz(2.1352467) q[1];
rz(-2.3079268) q[3];
sx q[3];
rz(-0.97396894) q[3];
sx q[3];
rz(-0.65469974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5323083) q[2];
sx q[2];
rz(-2.534635) q[2];
sx q[2];
rz(-0.92030805) q[2];
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
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9698708) q[0];
sx q[0];
rz(-1.0665749) q[0];
sx q[0];
rz(2.875476) q[0];
rz(1.4316106) q[1];
sx q[1];
rz(-2.0224729) q[1];
sx q[1];
rz(-2.5439579) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.584883) q[0];
sx q[0];
rz(-2.1765255) q[0];
sx q[0];
rz(-1.8023876) q[0];
rz(-pi) q[1];
x q[1];
rz(1.581988) q[2];
sx q[2];
rz(-1.7518242) q[2];
sx q[2];
rz(-2.343296) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7998593) q[1];
sx q[1];
rz(-1.5149773) q[1];
sx q[1];
rz(-0.69715211) q[1];
rz(1.8915755) q[3];
sx q[3];
rz(-1.7177542) q[3];
sx q[3];
rz(-3.1023417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0133682) q[2];
sx q[2];
rz(-0.2672264) q[2];
sx q[2];
rz(-2.5426148) q[2];
rz(-2.525575) q[3];
sx q[3];
rz(-0.44121656) q[3];
sx q[3];
rz(1.0429355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.041158572) q[0];
sx q[0];
rz(-2.2293595) q[0];
sx q[0];
rz(0.31901932) q[0];
rz(0.11095412) q[1];
sx q[1];
rz(-1.7496505) q[1];
sx q[1];
rz(2.9999733) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7975811) q[0];
sx q[0];
rz(-2.251029) q[0];
sx q[0];
rz(-1.6949095) q[0];
x q[1];
rz(-3.1371731) q[2];
sx q[2];
rz(-2.1221771) q[2];
sx q[2];
rz(2.1194292) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1846022) q[1];
sx q[1];
rz(-0.58035589) q[1];
sx q[1];
rz(2.7932211) q[1];
rz(-pi) q[2];
rz(-1.8291437) q[3];
sx q[3];
rz(-2.606539) q[3];
sx q[3];
rz(0.38219562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5854599) q[2];
sx q[2];
rz(-2.1161049) q[2];
sx q[2];
rz(-3.0312313) q[2];
rz(0.50802556) q[3];
sx q[3];
rz(-0.042162687) q[3];
sx q[3];
rz(-1.149811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0069649) q[0];
sx q[0];
rz(-0.65161324) q[0];
sx q[0];
rz(1.2192669) q[0];
rz(3.1189175) q[1];
sx q[1];
rz(-0.89212787) q[1];
sx q[1];
rz(-2.5882744) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39529453) q[0];
sx q[0];
rz(-2.2263701) q[0];
sx q[0];
rz(-3.0297999) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0108286) q[2];
sx q[2];
rz(-0.65245318) q[2];
sx q[2];
rz(2.3956232) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.07456116) q[1];
sx q[1];
rz(-2.0843049) q[1];
sx q[1];
rz(1.1719955) q[1];
x q[2];
rz(0.24711547) q[3];
sx q[3];
rz(-1.860642) q[3];
sx q[3];
rz(-0.71699504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.40552178) q[2];
sx q[2];
rz(-2.5260479) q[2];
sx q[2];
rz(-1.3496529) q[2];
rz(0.060529709) q[3];
sx q[3];
rz(-0.90977257) q[3];
sx q[3];
rz(2.5133666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26199207) q[0];
sx q[0];
rz(-0.99402004) q[0];
sx q[0];
rz(-0.0066198786) q[0];
rz(-2.5026542) q[1];
sx q[1];
rz(-0.21955755) q[1];
sx q[1];
rz(2.7117597) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2868067) q[0];
sx q[0];
rz(-1.9298975) q[0];
sx q[0];
rz(-0.39639985) q[0];
rz(-0.25385413) q[2];
sx q[2];
rz(-1.9905258) q[2];
sx q[2];
rz(-1.3042637) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4749516) q[1];
sx q[1];
rz(-2.3534497) q[1];
sx q[1];
rz(1.5402543) q[1];
rz(-pi) q[2];
rz(1.5019526) q[3];
sx q[3];
rz(-0.9864583) q[3];
sx q[3];
rz(2.7513954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2868353) q[2];
sx q[2];
rz(-1.9660051) q[2];
sx q[2];
rz(0.01469928) q[2];
rz(2.0059351) q[3];
sx q[3];
rz(-1.1123927) q[3];
sx q[3];
rz(1.8358102) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.084918) q[0];
sx q[0];
rz(-2.882353) q[0];
sx q[0];
rz(-1.2925451) q[0];
rz(2.9855285) q[1];
sx q[1];
rz(-0.81764692) q[1];
sx q[1];
rz(0.83052105) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9915402) q[0];
sx q[0];
rz(-0.73313289) q[0];
sx q[0];
rz(-1.7197648) q[0];
rz(-1.2654598) q[2];
sx q[2];
rz(-2.3407901) q[2];
sx q[2];
rz(2.3126471) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3780508) q[1];
sx q[1];
rz(-1.4908067) q[1];
sx q[1];
rz(1.3536118) q[1];
x q[2];
rz(-2.7748608) q[3];
sx q[3];
rz(-1.368513) q[3];
sx q[3];
rz(-2.3679581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.10892756) q[2];
sx q[2];
rz(-0.69643164) q[2];
sx q[2];
rz(2.6984974) q[2];
rz(0.14919925) q[3];
sx q[3];
rz(-2.7982893) q[3];
sx q[3];
rz(-2.8086015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14015848) q[0];
sx q[0];
rz(-2.5509111) q[0];
sx q[0];
rz(2.1668707) q[0];
rz(2.1433266) q[1];
sx q[1];
rz(-1.2405735) q[1];
sx q[1];
rz(2.1539198) q[1];
rz(0.51973612) q[2];
sx q[2];
rz(-0.39032857) q[2];
sx q[2];
rz(-2.6662568) q[2];
rz(-2.7137518) q[3];
sx q[3];
rz(-2.1890752) q[3];
sx q[3];
rz(-1.5933468) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
