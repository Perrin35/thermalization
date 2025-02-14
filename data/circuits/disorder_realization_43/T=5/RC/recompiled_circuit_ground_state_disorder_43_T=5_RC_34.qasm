OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2744098) q[0];
sx q[0];
rz(2.6743968) q[0];
sx q[0];
rz(13.462486) q[0];
rz(2.3236302) q[1];
sx q[1];
rz(-1.5939413) q[1];
sx q[1];
rz(-2.154248) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86292279) q[0];
sx q[0];
rz(-2.3723944) q[0];
sx q[0];
rz(2.132036) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.82168401) q[2];
sx q[2];
rz(-1.1237306) q[2];
sx q[2];
rz(0.38082235) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3707177) q[1];
sx q[1];
rz(-0.32568078) q[1];
sx q[1];
rz(2.5925255) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4127068) q[3];
sx q[3];
rz(-1.9327379) q[3];
sx q[3];
rz(-2.2666422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.92153111) q[2];
sx q[2];
rz(-2.107373) q[2];
sx q[2];
rz(-0.49911487) q[2];
rz(0.81579298) q[3];
sx q[3];
rz(-0.59320265) q[3];
sx q[3];
rz(-0.35713404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(2.4895184) q[0];
sx q[0];
rz(-0.3638142) q[0];
sx q[0];
rz(2.9079085) q[0];
rz(3.0417327) q[1];
sx q[1];
rz(-1.654511) q[1];
sx q[1];
rz(-1.608009) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8692538) q[0];
sx q[0];
rz(-1.4492479) q[0];
sx q[0];
rz(3.0442793) q[0];
rz(2.7614715) q[2];
sx q[2];
rz(-0.46762662) q[2];
sx q[2];
rz(-1.1884226) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2054384) q[1];
sx q[1];
rz(-2.1742651) q[1];
sx q[1];
rz(-0.11118576) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30939524) q[3];
sx q[3];
rz(-1.1077048) q[3];
sx q[3];
rz(-2.9675296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0443772) q[2];
sx q[2];
rz(-0.88572398) q[2];
sx q[2];
rz(3.0866747) q[2];
rz(0.6461668) q[3];
sx q[3];
rz(-1.5271527) q[3];
sx q[3];
rz(0.072877876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1845448) q[0];
sx q[0];
rz(-0.2453198) q[0];
sx q[0];
rz(-0.74485892) q[0];
rz(-1.8648196) q[1];
sx q[1];
rz(-1.0294015) q[1];
sx q[1];
rz(2.2778817) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7521356) q[0];
sx q[0];
rz(-0.4438209) q[0];
sx q[0];
rz(-3.0054166) q[0];
x q[1];
rz(-2.6361739) q[2];
sx q[2];
rz(-0.8552455) q[2];
sx q[2];
rz(0.16929132) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.59434592) q[1];
sx q[1];
rz(-2.1859096) q[1];
sx q[1];
rz(2.1662461) q[1];
rz(-pi) q[2];
rz(-2.4808472) q[3];
sx q[3];
rz(-1.9781986) q[3];
sx q[3];
rz(-1.7254207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.471571) q[2];
sx q[2];
rz(-1.6100641) q[2];
sx q[2];
rz(-0.395533) q[2];
rz(-1.0223355) q[3];
sx q[3];
rz(-2.6774355) q[3];
sx q[3];
rz(-2.6497604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1763497) q[0];
sx q[0];
rz(-1.1818161) q[0];
sx q[0];
rz(-3.0856207) q[0];
rz(-2.5296027) q[1];
sx q[1];
rz(-1.5490218) q[1];
sx q[1];
rz(2.2956119) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9214981) q[0];
sx q[0];
rz(-1.8290231) q[0];
sx q[0];
rz(-1.0357712) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1899784) q[2];
sx q[2];
rz(-1.8581094) q[2];
sx q[2];
rz(-3.0653022) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8710259) q[1];
sx q[1];
rz(-1.881384) q[1];
sx q[1];
rz(1.3695716) q[1];
rz(-pi) q[2];
rz(-2.2262283) q[3];
sx q[3];
rz(-0.5753606) q[3];
sx q[3];
rz(1.8917454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0116288) q[2];
sx q[2];
rz(-1.446898) q[2];
sx q[2];
rz(2.4212627) q[2];
rz(-0.35308853) q[3];
sx q[3];
rz(-1.947764) q[3];
sx q[3];
rz(-0.24875719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7406152) q[0];
sx q[0];
rz(-1.688513) q[0];
sx q[0];
rz(0.98689669) q[0];
rz(-1.997021) q[1];
sx q[1];
rz(-2.3026376) q[1];
sx q[1];
rz(-0.95485895) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83897774) q[0];
sx q[0];
rz(-0.90388008) q[0];
sx q[0];
rz(0.033292183) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0097505) q[2];
sx q[2];
rz(-1.2905143) q[2];
sx q[2];
rz(0.55299475) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.25453) q[1];
sx q[1];
rz(-1.945244) q[1];
sx q[1];
rz(2.6545054) q[1];
rz(-0.14820672) q[3];
sx q[3];
rz(-1.3199521) q[3];
sx q[3];
rz(-2.776752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6285051) q[2];
sx q[2];
rz(-2.2272765) q[2];
sx q[2];
rz(1.0465735) q[2];
rz(0.70932499) q[3];
sx q[3];
rz(-1.1449292) q[3];
sx q[3];
rz(0.63372248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76383048) q[0];
sx q[0];
rz(-1.4077633) q[0];
sx q[0];
rz(-0.27107006) q[0];
rz(3.0625878) q[1];
sx q[1];
rz(-2.7075691) q[1];
sx q[1];
rz(1.7105506) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9535429) q[0];
sx q[0];
rz(-2.8859647) q[0];
sx q[0];
rz(0.60582692) q[0];
rz(-pi) q[1];
rz(2.1151401) q[2];
sx q[2];
rz(-2.0664795) q[2];
sx q[2];
rz(0.069725603) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8003098) q[1];
sx q[1];
rz(-0.73004111) q[1];
sx q[1];
rz(-2.7679539) q[1];
x q[2];
rz(1.9892788) q[3];
sx q[3];
rz(-1.7177204) q[3];
sx q[3];
rz(2.4148586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.048910353) q[2];
sx q[2];
rz(-1.9581257) q[2];
sx q[2];
rz(0.023822039) q[2];
rz(-1.9134936) q[3];
sx q[3];
rz(-2.7619599) q[3];
sx q[3];
rz(-0.14321271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.989885) q[0];
sx q[0];
rz(-0.30549529) q[0];
sx q[0];
rz(-0.09224961) q[0];
rz(-2.8406738) q[1];
sx q[1];
rz(-1.9124799) q[1];
sx q[1];
rz(1.0801962) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55744074) q[0];
sx q[0];
rz(-1.6655465) q[0];
sx q[0];
rz(-1.9173724) q[0];
rz(-pi) q[1];
x q[1];
rz(1.148573) q[2];
sx q[2];
rz(-1.9034837) q[2];
sx q[2];
rz(-2.0300558) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2094592) q[1];
sx q[1];
rz(-1.5725296) q[1];
sx q[1];
rz(-0.14793795) q[1];
rz(-pi) q[2];
x q[2];
rz(0.69484922) q[3];
sx q[3];
rz(-1.3245305) q[3];
sx q[3];
rz(1.6978557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.13176189) q[2];
sx q[2];
rz(-2.7907382) q[2];
sx q[2];
rz(-2.5133613) q[2];
rz(2.175711) q[3];
sx q[3];
rz(-1.5647669) q[3];
sx q[3];
rz(-2.8204744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0308762) q[0];
sx q[0];
rz(-2.4346011) q[0];
sx q[0];
rz(1.4068756) q[0];
rz(-3.0137317) q[1];
sx q[1];
rz(-1.7117932) q[1];
sx q[1];
rz(1.1258639) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4425621) q[0];
sx q[0];
rz(-1.3872212) q[0];
sx q[0];
rz(2.4126614) q[0];
rz(-0.01861219) q[2];
sx q[2];
rz(-2.0224704) q[2];
sx q[2];
rz(-1.9960595) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5347157) q[1];
sx q[1];
rz(-1.770853) q[1];
sx q[1];
rz(1.5061395) q[1];
rz(-0.58784318) q[3];
sx q[3];
rz(-0.70826642) q[3];
sx q[3];
rz(1.4204587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.13059482) q[2];
sx q[2];
rz(-1.1329634) q[2];
sx q[2];
rz(3.0478743) q[2];
rz(-1.4334009) q[3];
sx q[3];
rz(-0.283537) q[3];
sx q[3];
rz(1.7482429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4505287) q[0];
sx q[0];
rz(-1.0858303) q[0];
sx q[0];
rz(-2.1375256) q[0];
rz(2.0380691) q[1];
sx q[1];
rz(-0.48201489) q[1];
sx q[1];
rz(1.4080661) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1188977) q[0];
sx q[0];
rz(-1.0556882) q[0];
sx q[0];
rz(-1.5777508) q[0];
rz(-2.8292921) q[2];
sx q[2];
rz(-1.4648629) q[2];
sx q[2];
rz(-0.16367975) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1939331) q[1];
sx q[1];
rz(-2.0037738) q[1];
sx q[1];
rz(0.19755944) q[1];
x q[2];
rz(0.51897612) q[3];
sx q[3];
rz(-1.5277315) q[3];
sx q[3];
rz(1.8881292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.39190009) q[2];
sx q[2];
rz(-1.6204648) q[2];
sx q[2];
rz(2.5938972) q[2];
rz(2.8294166) q[3];
sx q[3];
rz(-1.7269937) q[3];
sx q[3];
rz(-1.4148022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.80426973) q[0];
sx q[0];
rz(-1.6707358) q[0];
sx q[0];
rz(-0.7789337) q[0];
rz(-0.3745105) q[1];
sx q[1];
rz(-1.51314) q[1];
sx q[1];
rz(-0.83190727) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13369416) q[0];
sx q[0];
rz(-1.3092894) q[0];
sx q[0];
rz(-0.71870872) q[0];
x q[1];
rz(2.4291762) q[2];
sx q[2];
rz(-2.2210026) q[2];
sx q[2];
rz(0.67931108) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.1268046) q[1];
sx q[1];
rz(-2.2378057) q[1];
sx q[1];
rz(-2.2875525) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5396773) q[3];
sx q[3];
rz(-2.3672769) q[3];
sx q[3];
rz(-1.7208791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1067918) q[2];
sx q[2];
rz(-2.782244) q[2];
sx q[2];
rz(0.0084776004) q[2];
rz(0.40555412) q[3];
sx q[3];
rz(-1.4529934) q[3];
sx q[3];
rz(-1.6976374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1494898) q[0];
sx q[0];
rz(-1.533951) q[0];
sx q[0];
rz(0.78716192) q[0];
rz(-2.0472732) q[1];
sx q[1];
rz(-2.7631187) q[1];
sx q[1];
rz(-0.43597058) q[1];
rz(0.12577429) q[2];
sx q[2];
rz(-1.3604506) q[2];
sx q[2];
rz(2.4764555) q[2];
rz(0.094219128) q[3];
sx q[3];
rz(-1.5658656) q[3];
sx q[3];
rz(1.6938536) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
