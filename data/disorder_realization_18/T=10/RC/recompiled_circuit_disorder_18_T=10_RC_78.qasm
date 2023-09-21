OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.82233959) q[0];
sx q[0];
rz(-0.95681325) q[0];
sx q[0];
rz(1.4332888) q[0];
rz(-0.23694555) q[1];
sx q[1];
rz(4.3720923) q[1];
sx q[1];
rz(11.421539) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24554907) q[0];
sx q[0];
rz(-2.890812) q[0];
sx q[0];
rz(-1.301469) q[0];
x q[1];
rz(1.990591) q[2];
sx q[2];
rz(-1.8325873) q[2];
sx q[2];
rz(-2.7671438) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1022542) q[1];
sx q[1];
rz(-1.5896817) q[1];
sx q[1];
rz(2.0233872) q[1];
rz(-pi) q[2];
rz(-0.54361312) q[3];
sx q[3];
rz(-1.8279148) q[3];
sx q[3];
rz(-1.4527814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.11674374) q[2];
sx q[2];
rz(-1.8086834) q[2];
sx q[2];
rz(1.2343181) q[2];
rz(1.0553137) q[3];
sx q[3];
rz(-2.0827259) q[3];
sx q[3];
rz(1.9887259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18594436) q[0];
sx q[0];
rz(-1.9444822) q[0];
sx q[0];
rz(-0.82988513) q[0];
rz(0.25289598) q[1];
sx q[1];
rz(-2.3650832) q[1];
sx q[1];
rz(-2.2944962) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9053099) q[0];
sx q[0];
rz(-0.61203996) q[0];
sx q[0];
rz(0.73487868) q[0];
rz(-pi) q[1];
rz(0.77913021) q[2];
sx q[2];
rz(-1.7841633) q[2];
sx q[2];
rz(1.9386292) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9287195) q[1];
sx q[1];
rz(-1.6391014) q[1];
sx q[1];
rz(-2.5961155) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2291525) q[3];
sx q[3];
rz(-2.6077301) q[3];
sx q[3];
rz(-0.57953366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0236686) q[2];
sx q[2];
rz(-2.7682722) q[2];
sx q[2];
rz(2.4222597) q[2];
rz(-1.8524648) q[3];
sx q[3];
rz(-2.9057861) q[3];
sx q[3];
rz(-1.6524564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8711202) q[0];
sx q[0];
rz(-1.2501165) q[0];
sx q[0];
rz(-0.57139325) q[0];
rz(0.96673036) q[1];
sx q[1];
rz(-1.8833908) q[1];
sx q[1];
rz(-2.6142696) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8910599) q[0];
sx q[0];
rz(-2.8322304) q[0];
sx q[0];
rz(0.022457794) q[0];
rz(-0.40076077) q[2];
sx q[2];
rz(-0.71360613) q[2];
sx q[2];
rz(1.00373) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.61422435) q[1];
sx q[1];
rz(-2.3389611) q[1];
sx q[1];
rz(-2.7781092) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.73266352) q[3];
sx q[3];
rz(-1.3244197) q[3];
sx q[3];
rz(0.25783595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8335235) q[2];
sx q[2];
rz(-1.4493194) q[2];
sx q[2];
rz(0.72675881) q[2];
rz(-2.1440992) q[3];
sx q[3];
rz(-2.5163429) q[3];
sx q[3];
rz(1.6231026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24546394) q[0];
sx q[0];
rz(-1.6587057) q[0];
sx q[0];
rz(-2.1380651) q[0];
rz(-3.1009122) q[1];
sx q[1];
rz(-1.129312) q[1];
sx q[1];
rz(-2.3153268) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6986615) q[0];
sx q[0];
rz(-0.94416617) q[0];
sx q[0];
rz(-2.4311964) q[0];
rz(-pi) q[1];
x q[1];
rz(0.70010186) q[2];
sx q[2];
rz(-1.6119909) q[2];
sx q[2];
rz(-2.5259279) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.367825) q[1];
sx q[1];
rz(-2.7959679) q[1];
sx q[1];
rz(-1.7961545) q[1];
rz(2.5097333) q[3];
sx q[3];
rz(-1.251289) q[3];
sx q[3];
rz(0.39998049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.56759175) q[2];
sx q[2];
rz(-1.8969798) q[2];
sx q[2];
rz(0.91439247) q[2];
rz(3.0363723) q[3];
sx q[3];
rz(-1.6243694) q[3];
sx q[3];
rz(2.5523394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2545664) q[0];
sx q[0];
rz(-1.5019324) q[0];
sx q[0];
rz(-2.0078833) q[0];
rz(3.1359613) q[1];
sx q[1];
rz(-1.7040323) q[1];
sx q[1];
rz(-0.61757913) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3591946) q[0];
sx q[0];
rz(-3.087128) q[0];
sx q[0];
rz(-0.99476238) q[0];
x q[1];
rz(-0.56474833) q[2];
sx q[2];
rz(-2.4167633) q[2];
sx q[2];
rz(0.85541475) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6981814) q[1];
sx q[1];
rz(-2.4907506) q[1];
sx q[1];
rz(-1.0990259) q[1];
rz(-pi) q[2];
rz(-1.5836309) q[3];
sx q[3];
rz(-1.7463297) q[3];
sx q[3];
rz(-2.3042391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.79409838) q[2];
sx q[2];
rz(-1.2604875) q[2];
sx q[2];
rz(2.9079672) q[2];
rz(0.99308333) q[3];
sx q[3];
rz(-0.94998327) q[3];
sx q[3];
rz(-2.5036507) q[3];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96930209) q[0];
sx q[0];
rz(-3.0811221) q[0];
sx q[0];
rz(-0.0897952) q[0];
rz(0.88109294) q[1];
sx q[1];
rz(-2.6125364) q[1];
sx q[1];
rz(-0.18009137) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0952547) q[0];
sx q[0];
rz(-1.006608) q[0];
sx q[0];
rz(-0.5395704) q[0];
rz(-0.18896582) q[2];
sx q[2];
rz(-2.1660888) q[2];
sx q[2];
rz(1.8408066) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.60702) q[1];
sx q[1];
rz(-1.8587451) q[1];
sx q[1];
rz(1.401591) q[1];
rz(-pi) q[2];
rz(-2.0635707) q[3];
sx q[3];
rz(-1.5544976) q[3];
sx q[3];
rz(-1.8340045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.2580516) q[2];
sx q[2];
rz(-1.5641314) q[2];
sx q[2];
rz(1.9656666) q[2];
rz(-3.0684493) q[3];
sx q[3];
rz(-2.4172343) q[3];
sx q[3];
rz(-2.6426962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95773762) q[0];
sx q[0];
rz(-2.4860005) q[0];
sx q[0];
rz(1.7234329) q[0];
rz(-1.1649959) q[1];
sx q[1];
rz(-2.5823451) q[1];
sx q[1];
rz(-3.1013536) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2049853) q[0];
sx q[0];
rz(-1.4892704) q[0];
sx q[0];
rz(1.3315014) q[0];
rz(1.4637429) q[2];
sx q[2];
rz(-0.40823001) q[2];
sx q[2];
rz(1.5232616) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.28133) q[1];
sx q[1];
rz(-1.7107692) q[1];
sx q[1];
rz(-2.1920188) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8242565) q[3];
sx q[3];
rz(-1.8620957) q[3];
sx q[3];
rz(2.3288162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.13850257) q[2];
sx q[2];
rz(-0.78531992) q[2];
sx q[2];
rz(2.9677532) q[2];
rz(-1.396817) q[3];
sx q[3];
rz(-1.4649748) q[3];
sx q[3];
rz(-2.2824536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.023507) q[0];
sx q[0];
rz(-0.21723391) q[0];
sx q[0];
rz(-1.404495) q[0];
rz(-1.9346168) q[1];
sx q[1];
rz(-1.6563621) q[1];
sx q[1];
rz(1.6361902) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0818427) q[0];
sx q[0];
rz(-1.2568226) q[0];
sx q[0];
rz(0.44072515) q[0];
rz(0.90754857) q[2];
sx q[2];
rz(-2.2758256) q[2];
sx q[2];
rz(0.384534) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.38547388) q[1];
sx q[1];
rz(-1.5382775) q[1];
sx q[1];
rz(2.4454152) q[1];
rz(-0.96945854) q[3];
sx q[3];
rz(-1.2974032) q[3];
sx q[3];
rz(2.6858342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.76190844) q[2];
sx q[2];
rz(-0.48019335) q[2];
sx q[2];
rz(-2.9910679) q[2];
rz(1.5395509) q[3];
sx q[3];
rz(-1.9463041) q[3];
sx q[3];
rz(-2.5185744) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6983011) q[0];
sx q[0];
rz(-1.1346096) q[0];
sx q[0];
rz(0.64569965) q[0];
rz(0.49939108) q[1];
sx q[1];
rz(-1.4285409) q[1];
sx q[1];
rz(-1.8766778) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7495959) q[0];
sx q[0];
rz(-0.8526593) q[0];
sx q[0];
rz(2.0977661) q[0];
x q[1];
rz(1.132498) q[2];
sx q[2];
rz(-2.0813137) q[2];
sx q[2];
rz(-2.8709708) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4134819) q[1];
sx q[1];
rz(-2.9915447) q[1];
sx q[1];
rz(-2.8996182) q[1];
rz(-pi) q[2];
rz(-2.7154269) q[3];
sx q[3];
rz(-2.6821972) q[3];
sx q[3];
rz(-3.1118432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.634793) q[2];
sx q[2];
rz(-2.1719077) q[2];
sx q[2];
rz(-0.52465049) q[2];
rz(-2.5850885) q[3];
sx q[3];
rz(-1.6289214) q[3];
sx q[3];
rz(2.7867253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.5234579) q[0];
sx q[0];
rz(-0.83704656) q[0];
sx q[0];
rz(2.8961704) q[0];
rz(-0.55631176) q[1];
sx q[1];
rz(-1.6106771) q[1];
sx q[1];
rz(1.6773178) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3905555) q[0];
sx q[0];
rz(-1.7883736) q[0];
sx q[0];
rz(1.6839954) q[0];
x q[1];
rz(-2.1572838) q[2];
sx q[2];
rz(-1.3830292) q[2];
sx q[2];
rz(-0.078660065) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6357395) q[1];
sx q[1];
rz(-0.74843279) q[1];
sx q[1];
rz(1.1670508) q[1];
rz(3.0562923) q[3];
sx q[3];
rz(-1.395426) q[3];
sx q[3];
rz(0.264197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8993373) q[2];
sx q[2];
rz(-1.238404) q[2];
sx q[2];
rz(2.396092) q[2];
rz(1.7177104) q[3];
sx q[3];
rz(-2.8427377) q[3];
sx q[3];
rz(2.0132813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-1.2765008) q[0];
sx q[0];
rz(-1.4958953) q[0];
sx q[0];
rz(-2.4687742) q[0];
rz(-1.2561692) q[1];
sx q[1];
rz(-2.3354463) q[1];
sx q[1];
rz(-1.068402) q[1];
rz(-1.0829265) q[2];
sx q[2];
rz(-0.91960533) q[2];
sx q[2];
rz(-1.4056924) q[2];
rz(-0.54995723) q[3];
sx q[3];
rz(-0.64823845) q[3];
sx q[3];
rz(1.3241495) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];