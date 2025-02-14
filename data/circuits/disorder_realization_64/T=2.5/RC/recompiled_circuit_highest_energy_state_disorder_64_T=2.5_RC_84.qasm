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
rz(3.1177899) q[0];
sx q[0];
rz(-1.1060214) q[0];
sx q[0];
rz(2.3646781) q[0];
rz(1.39224) q[1];
sx q[1];
rz(-1.3148146) q[1];
sx q[1];
rz(2.1652752) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1084749) q[0];
sx q[0];
rz(-0.43370789) q[0];
sx q[0];
rz(1.8932305) q[0];
rz(-pi) q[1];
x q[1];
rz(0.10586057) q[2];
sx q[2];
rz(-2.5753394) q[2];
sx q[2];
rz(2.16118) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9982517) q[1];
sx q[1];
rz(-2.6282881) q[1];
sx q[1];
rz(0.99350382) q[1];
x q[2];
rz(1.5209274) q[3];
sx q[3];
rz(-2.5679776) q[3];
sx q[3];
rz(0.22465868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.50600791) q[2];
sx q[2];
rz(-3.0878461) q[2];
sx q[2];
rz(0.40679833) q[2];
rz(-0.16945101) q[3];
sx q[3];
rz(-0.5295161) q[3];
sx q[3];
rz(-1.0725526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2605543) q[0];
sx q[0];
rz(-0.23242234) q[0];
sx q[0];
rz(-3.1378003) q[0];
rz(-0.077839851) q[1];
sx q[1];
rz(-2.4816315) q[1];
sx q[1];
rz(-0.30581623) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90214848) q[0];
sx q[0];
rz(-1.3864707) q[0];
sx q[0];
rz(2.091616) q[0];
x q[1];
rz(-2.3897116) q[2];
sx q[2];
rz(-1.6691748) q[2];
sx q[2];
rz(1.4954612) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9287323) q[1];
sx q[1];
rz(-1.9668192) q[1];
sx q[1];
rz(-2.3235882) q[1];
rz(1.5068077) q[3];
sx q[3];
rz(-1.631569) q[3];
sx q[3];
rz(-1.9072541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.99016142) q[2];
sx q[2];
rz(-1.7502681) q[2];
sx q[2];
rz(-0.64275098) q[2];
rz(0.07240545) q[3];
sx q[3];
rz(-2.0824771) q[3];
sx q[3];
rz(-1.36093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0533503) q[0];
sx q[0];
rz(-2.998816) q[0];
sx q[0];
rz(-2.9852168) q[0];
rz(0.0414255) q[1];
sx q[1];
rz(-2.5138469) q[1];
sx q[1];
rz(-1.5511537) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0254733) q[0];
sx q[0];
rz(-1.0295233) q[0];
sx q[0];
rz(-1.3295637) q[0];
rz(0.088038283) q[2];
sx q[2];
rz(-2.4836342) q[2];
sx q[2];
rz(-2.7245299) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.98387466) q[1];
sx q[1];
rz(-1.0441171) q[1];
sx q[1];
rz(1.4406564) q[1];
rz(2.2699192) q[3];
sx q[3];
rz(-2.1101885) q[3];
sx q[3];
rz(-1.9339069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7330043) q[2];
sx q[2];
rz(-0.53885794) q[2];
sx q[2];
rz(3.1206701) q[2];
rz(2.9544592) q[3];
sx q[3];
rz(-2.9360866) q[3];
sx q[3];
rz(3.0278897) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9631831) q[0];
sx q[0];
rz(-0.33247501) q[0];
sx q[0];
rz(2.7030429) q[0];
rz(1.6167538) q[1];
sx q[1];
rz(-2.8068145) q[1];
sx q[1];
rz(2.8964892) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6436076) q[0];
sx q[0];
rz(-1.6542572) q[0];
sx q[0];
rz(-2.2913234) q[0];
rz(-pi) q[1];
rz(0.39056492) q[2];
sx q[2];
rz(-2.7465944) q[2];
sx q[2];
rz(-0.67827889) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.41293834) q[1];
sx q[1];
rz(-1.0333583) q[1];
sx q[1];
rz(2.9793903) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2538337) q[3];
sx q[3];
rz(-2.2989103) q[3];
sx q[3];
rz(-1.1945981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.601292) q[2];
sx q[2];
rz(-0.43593323) q[2];
sx q[2];
rz(0.33176804) q[2];
rz(2.6541384) q[3];
sx q[3];
rz(-2.09477) q[3];
sx q[3];
rz(0.99307466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29193923) q[0];
sx q[0];
rz(-1.6878457) q[0];
sx q[0];
rz(-2.3680903) q[0];
rz(-1.1812814) q[1];
sx q[1];
rz(-2.9995194) q[1];
sx q[1];
rz(1.7519417) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6365125) q[0];
sx q[0];
rz(-2.5861222) q[0];
sx q[0];
rz(-0.73969264) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2821065) q[2];
sx q[2];
rz(-1.0946678) q[2];
sx q[2];
rz(1.7993594) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.70848318) q[1];
sx q[1];
rz(-0.481284) q[1];
sx q[1];
rz(-2.9631056) q[1];
rz(-0.23791194) q[3];
sx q[3];
rz(-1.3945268) q[3];
sx q[3];
rz(-2.545994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2191849) q[2];
sx q[2];
rz(-1.9318523) q[2];
sx q[2];
rz(-2.6694471) q[2];
rz(1.2989429) q[3];
sx q[3];
rz(-1.2932212) q[3];
sx q[3];
rz(-0.81056547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2209114) q[0];
sx q[0];
rz(-0.26316106) q[0];
sx q[0];
rz(-0.26350185) q[0];
rz(1.1031411) q[1];
sx q[1];
rz(-1.3145072) q[1];
sx q[1];
rz(-2.7679494) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5140681) q[0];
sx q[0];
rz(-1.6434945) q[0];
sx q[0];
rz(0.060427314) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0194026) q[2];
sx q[2];
rz(-2.2940679) q[2];
sx q[2];
rz(3.0360589) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3195575) q[1];
sx q[1];
rz(-0.6614092) q[1];
sx q[1];
rz(-2.103785) q[1];
x q[2];
rz(-0.5881891) q[3];
sx q[3];
rz(-2.0547227) q[3];
sx q[3];
rz(-1.5671135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0282447) q[2];
sx q[2];
rz(-2.9714606) q[2];
sx q[2];
rz(-0.55383468) q[2];
rz(1.3977741) q[3];
sx q[3];
rz(-0.60269409) q[3];
sx q[3];
rz(-0.3127313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58436191) q[0];
sx q[0];
rz(-2.1350242) q[0];
sx q[0];
rz(-1.2297909) q[0];
rz(0.23904414) q[1];
sx q[1];
rz(-1.6300647) q[1];
sx q[1];
rz(-0.30034932) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2588033) q[0];
sx q[0];
rz(-0.20016328) q[0];
sx q[0];
rz(0.61997719) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0763758) q[2];
sx q[2];
rz(-2.2110614) q[2];
sx q[2];
rz(-0.51039052) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3546037) q[1];
sx q[1];
rz(-1.5704078) q[1];
sx q[1];
rz(1.5860737) q[1];
rz(0.024552931) q[3];
sx q[3];
rz(-2.2835287) q[3];
sx q[3];
rz(0.36125444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0099237) q[2];
sx q[2];
rz(-1.6841623) q[2];
sx q[2];
rz(-2.8966676) q[2];
rz(-2.6217672) q[3];
sx q[3];
rz(-0.86331415) q[3];
sx q[3];
rz(-2.4533217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35172611) q[0];
sx q[0];
rz(-2.7930197) q[0];
sx q[0];
rz(-1.9785731) q[0];
rz(3.0746958) q[1];
sx q[1];
rz(-1.6480548) q[1];
sx q[1];
rz(-2.1287207) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6352309) q[0];
sx q[0];
rz(-1.073045) q[0];
sx q[0];
rz(-2.228581) q[0];
rz(1.0032907) q[2];
sx q[2];
rz(-1.7219647) q[2];
sx q[2];
rz(-2.0625819) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7808395) q[1];
sx q[1];
rz(-1.8855125) q[1];
sx q[1];
rz(0.34784045) q[1];
rz(2.4550405) q[3];
sx q[3];
rz(-1.9007389) q[3];
sx q[3];
rz(-0.8984962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0248727) q[2];
sx q[2];
rz(-0.97815424) q[2];
sx q[2];
rz(-2.8187974) q[2];
rz(0.60574496) q[3];
sx q[3];
rz(-2.3492458) q[3];
sx q[3];
rz(2.7927223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.087273) q[0];
sx q[0];
rz(-0.1467341) q[0];
sx q[0];
rz(0.012454575) q[0];
rz(-2.3948578) q[1];
sx q[1];
rz(-0.92339271) q[1];
sx q[1];
rz(-2.8616203) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6387647) q[0];
sx q[0];
rz(-1.6831213) q[0];
sx q[0];
rz(3.0905484) q[0];
rz(-1.1367646) q[2];
sx q[2];
rz(-1.869259) q[2];
sx q[2];
rz(-3.0500183) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2586883) q[1];
sx q[1];
rz(-1.2637485) q[1];
sx q[1];
rz(-2.3767002) q[1];
x q[2];
rz(-2.1567508) q[3];
sx q[3];
rz(-2.7142314) q[3];
sx q[3];
rz(-2.7387184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.81194699) q[2];
sx q[2];
rz(-2.3880366) q[2];
sx q[2];
rz(2.9296056) q[2];
rz(0.82344615) q[3];
sx q[3];
rz(-1.4444838) q[3];
sx q[3];
rz(2.8823891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9987746) q[0];
sx q[0];
rz(-0.057567216) q[0];
sx q[0];
rz(2.4488191) q[0];
rz(-2.5686) q[1];
sx q[1];
rz(-1.3447821) q[1];
sx q[1];
rz(-0.43100345) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0045872) q[0];
sx q[0];
rz(-2.1218897) q[0];
sx q[0];
rz(-1.5020834) q[0];
x q[1];
rz(-0.60009267) q[2];
sx q[2];
rz(-1.1962435) q[2];
sx q[2];
rz(-0.77175826) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9622456) q[1];
sx q[1];
rz(-1.822346) q[1];
sx q[1];
rz(0.72466447) q[1];
rz(-pi) q[2];
rz(1.0209084) q[3];
sx q[3];
rz(-1.3088994) q[3];
sx q[3];
rz(-0.49347116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7196322) q[2];
sx q[2];
rz(-2.8675291) q[2];
sx q[2];
rz(-2.5813622) q[2];
rz(-2.6719921) q[3];
sx q[3];
rz(-2.738939) q[3];
sx q[3];
rz(0.71389055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8568759) q[0];
sx q[0];
rz(-1.7354043) q[0];
sx q[0];
rz(2.0866557) q[0];
rz(0.84125413) q[1];
sx q[1];
rz(-1.1019191) q[1];
sx q[1];
rz(-0.046774653) q[1];
rz(-0.31568676) q[2];
sx q[2];
rz(-1.9236947) q[2];
sx q[2];
rz(-2.535939) q[2];
rz(3.1145949) q[3];
sx q[3];
rz(-0.91201966) q[3];
sx q[3];
rz(-2.1569679) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
