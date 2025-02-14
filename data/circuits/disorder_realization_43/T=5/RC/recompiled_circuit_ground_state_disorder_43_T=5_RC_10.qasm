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
rz(-0.46719587) q[0];
sx q[0];
rz(-0.89611563) q[0];
rz(2.3236302) q[1];
sx q[1];
rz(-1.5939413) q[1];
sx q[1];
rz(0.98734468) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2786699) q[0];
sx q[0];
rz(-0.76919829) q[0];
sx q[0];
rz(-1.0095566) q[0];
x q[1];
rz(2.1842264) q[2];
sx q[2];
rz(-2.2920458) q[2];
sx q[2];
rz(0.75480513) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8670013) q[1];
sx q[1];
rz(-1.7385635) q[1];
sx q[1];
rz(-0.28047362) q[1];
rz(-1.7288858) q[3];
sx q[3];
rz(-1.9327379) q[3];
sx q[3];
rz(-0.87495041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2200615) q[2];
sx q[2];
rz(-1.0342197) q[2];
sx q[2];
rz(2.6424778) q[2];
rz(2.3257997) q[3];
sx q[3];
rz(-2.54839) q[3];
sx q[3];
rz(-0.35713404) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65207425) q[0];
sx q[0];
rz(-2.7777785) q[0];
sx q[0];
rz(-0.2336842) q[0];
rz(3.0417327) q[1];
sx q[1];
rz(-1.4870817) q[1];
sx q[1];
rz(1.608009) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2723389) q[0];
sx q[0];
rz(-1.4492479) q[0];
sx q[0];
rz(0.097313332) q[0];
rz(-pi) q[1];
rz(2.7031021) q[2];
sx q[2];
rz(-1.4027565) q[2];
sx q[2];
rz(-2.4166256) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.69791693) q[1];
sx q[1];
rz(-1.47931) q[1];
sx q[1];
rz(-0.96443022) q[1];
rz(-pi) q[2];
rz(-2.053612) q[3];
sx q[3];
rz(-1.8466966) q[3];
sx q[3];
rz(1.538572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0443772) q[2];
sx q[2];
rz(-0.88572398) q[2];
sx q[2];
rz(-0.054917939) q[2];
rz(-0.6461668) q[3];
sx q[3];
rz(-1.5271527) q[3];
sx q[3];
rz(-0.072877876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1845448) q[0];
sx q[0];
rz(-0.2453198) q[0];
sx q[0];
rz(-2.3967337) q[0];
rz(1.8648196) q[1];
sx q[1];
rz(-2.1121912) q[1];
sx q[1];
rz(-0.863711) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.389457) q[0];
sx q[0];
rz(-0.4438209) q[0];
sx q[0];
rz(-0.13617604) q[0];
rz(-pi) q[1];
rz(-0.50541877) q[2];
sx q[2];
rz(-2.2863472) q[2];
sx q[2];
rz(0.16929132) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3491213) q[1];
sx q[1];
rz(-2.0466261) q[1];
sx q[1];
rz(0.70648593) q[1];
x q[2];
rz(-0.61285783) q[3];
sx q[3];
rz(-0.75988673) q[3];
sx q[3];
rz(2.5155544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6700217) q[2];
sx q[2];
rz(-1.5315285) q[2];
sx q[2];
rz(2.7460597) q[2];
rz(2.1192571) q[3];
sx q[3];
rz(-0.46415713) q[3];
sx q[3];
rz(2.6497604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-0.96524298) q[0];
sx q[0];
rz(-1.9597766) q[0];
sx q[0];
rz(-0.055971948) q[0];
rz(-2.5296027) q[1];
sx q[1];
rz(-1.5925708) q[1];
sx q[1];
rz(0.8459808) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2200945) q[0];
sx q[0];
rz(-1.8290231) q[0];
sx q[0];
rz(2.1058215) q[0];
rz(-pi) q[1];
rz(-1.0998639) q[2];
sx q[2];
rz(-2.467017) q[2];
sx q[2];
rz(-1.872962) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8236134) q[1];
sx q[1];
rz(-0.3682963) q[1];
sx q[1];
rz(-0.55693926) q[1];
x q[2];
rz(-2.0457129) q[3];
sx q[3];
rz(-1.9088512) q[3];
sx q[3];
rz(2.2477704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1299639) q[2];
sx q[2];
rz(-1.6946946) q[2];
sx q[2];
rz(0.72032991) q[2];
rz(-2.7885041) q[3];
sx q[3];
rz(-1.947764) q[3];
sx q[3];
rz(0.24875719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40097749) q[0];
sx q[0];
rz(-1.4530797) q[0];
sx q[0];
rz(-0.98689669) q[0];
rz(1.1445716) q[1];
sx q[1];
rz(-0.83895504) q[1];
sx q[1];
rz(0.95485895) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78518822) q[0];
sx q[0];
rz(-2.4739728) q[0];
sx q[0];
rz(1.6130527) q[0];
x q[1];
rz(-2.1318421) q[2];
sx q[2];
rz(-1.8510783) q[2];
sx q[2];
rz(2.5885979) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.92063078) q[1];
sx q[1];
rz(-0.60505962) q[1];
sx q[1];
rz(0.69843881) q[1];
x q[2];
rz(-0.14820672) q[3];
sx q[3];
rz(-1.3199521) q[3];
sx q[3];
rz(-2.776752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6285051) q[2];
sx q[2];
rz(-2.2272765) q[2];
sx q[2];
rz(-2.0950192) q[2];
rz(-2.4322677) q[3];
sx q[3];
rz(-1.9966634) q[3];
sx q[3];
rz(-0.63372248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3777622) q[0];
sx q[0];
rz(-1.4077633) q[0];
sx q[0];
rz(-0.27107006) q[0];
rz(-0.079004869) q[1];
sx q[1];
rz(-0.43402356) q[1];
sx q[1];
rz(-1.7105506) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9535429) q[0];
sx q[0];
rz(-0.255628) q[0];
sx q[0];
rz(0.60582692) q[0];
rz(-pi) q[1];
rz(-1.0264525) q[2];
sx q[2];
rz(-2.0664795) q[2];
sx q[2];
rz(0.069725603) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9985314) q[1];
sx q[1];
rz(-0.9009217) q[1];
sx q[1];
rz(1.8865442) q[1];
x q[2];
rz(1.2215516) q[3];
sx q[3];
rz(-2.6995097) q[3];
sx q[3];
rz(-0.52607049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0926823) q[2];
sx q[2];
rz(-1.9581257) q[2];
sx q[2];
rz(-0.023822039) q[2];
rz(1.9134936) q[3];
sx q[3];
rz(-2.7619599) q[3];
sx q[3];
rz(0.14321271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.989885) q[0];
sx q[0];
rz(-2.8360974) q[0];
sx q[0];
rz(0.09224961) q[0];
rz(2.8406738) q[1];
sx q[1];
rz(-1.9124799) q[1];
sx q[1];
rz(-1.0801962) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55744074) q[0];
sx q[0];
rz(-1.6655465) q[0];
sx q[0];
rz(-1.9173724) q[0];
rz(-pi) q[1];
rz(-1.9930196) q[2];
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
rz(-1.9321334) q[1];
sx q[1];
rz(-1.5725296) q[1];
sx q[1];
rz(2.9936547) q[1];
rz(-pi) q[2];
rz(1.8870513) q[3];
sx q[3];
rz(-0.90083921) q[3];
sx q[3];
rz(0.073410598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.13176189) q[2];
sx q[2];
rz(-0.35085446) q[2];
sx q[2];
rz(-0.62823137) q[2];
rz(0.96588165) q[3];
sx q[3];
rz(-1.5647669) q[3];
sx q[3];
rz(-0.32111827) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0308762) q[0];
sx q[0];
rz(-2.4346011) q[0];
sx q[0];
rz(1.4068756) q[0];
rz(0.12786099) q[1];
sx q[1];
rz(-1.7117932) q[1];
sx q[1];
rz(-2.0157287) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28982535) q[0];
sx q[0];
rz(-0.85677108) q[0];
sx q[0];
rz(1.8147574) q[0];
rz(-3.1229805) q[2];
sx q[2];
rz(-2.0224704) q[2];
sx q[2];
rz(-1.1455331) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5347157) q[1];
sx q[1];
rz(-1.770853) q[1];
sx q[1];
rz(1.6354531) q[1];
rz(-1.12735) q[3];
sx q[3];
rz(-2.1428041) q[3];
sx q[3];
rz(1.0009022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0109978) q[2];
sx q[2];
rz(-1.1329634) q[2];
sx q[2];
rz(3.0478743) q[2];
rz(-1.7081918) q[3];
sx q[3];
rz(-2.8580557) q[3];
sx q[3];
rz(-1.3933498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.691064) q[0];
sx q[0];
rz(-2.0557623) q[0];
sx q[0];
rz(2.1375256) q[0];
rz(-2.0380691) q[1];
sx q[1];
rz(-2.6595778) q[1];
sx q[1];
rz(-1.7335266) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0085786) q[0];
sx q[0];
rz(-2.6264418) q[0];
sx q[0];
rz(0.012284474) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8292921) q[2];
sx q[2];
rz(-1.6767297) q[2];
sx q[2];
rz(2.9779129) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4346659) q[1];
sx q[1];
rz(-1.7499171) q[1];
sx q[1];
rz(2.0113025) q[1];
x q[2];
rz(0.086661913) q[3];
sx q[3];
rz(-0.52059697) q[3];
sx q[3];
rz(0.39256061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7496926) q[2];
sx q[2];
rz(-1.6204648) q[2];
sx q[2];
rz(2.5938972) q[2];
rz(-0.31217602) q[3];
sx q[3];
rz(-1.4145989) q[3];
sx q[3];
rz(1.4148022) q[3];
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
rz(2.3373229) q[0];
sx q[0];
rz(-1.6707358) q[0];
sx q[0];
rz(-0.7789337) q[0];
rz(-0.3745105) q[1];
sx q[1];
rz(-1.51314) q[1];
sx q[1];
rz(-0.83190727) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4171281) q[0];
sx q[0];
rz(-2.3848371) q[0];
sx q[0];
rz(0.3860851) q[0];
x q[1];
rz(2.2807715) q[2];
sx q[2];
rz(-2.2172625) q[2];
sx q[2];
rz(2.8614028) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.191997) q[1];
sx q[1];
rz(-2.1131385) q[1];
sx q[1];
rz(0.80710141) q[1];
rz(-pi) q[2];
x q[2];
rz(2.34487) q[3];
sx q[3];
rz(-1.5925538) q[3];
sx q[3];
rz(-0.12783229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1067918) q[2];
sx q[2];
rz(-2.782244) q[2];
sx q[2];
rz(-3.1331151) q[2];
rz(-2.7360385) q[3];
sx q[3];
rz(-1.6885992) q[3];
sx q[3];
rz(1.6976374) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99210284) q[0];
sx q[0];
rz(-1.6076417) q[0];
sx q[0];
rz(-2.3544307) q[0];
rz(-2.0472732) q[1];
sx q[1];
rz(-2.7631187) q[1];
sx q[1];
rz(-0.43597058) q[1];
rz(-0.12577429) q[2];
sx q[2];
rz(-1.7811421) q[2];
sx q[2];
rz(-0.66513715) q[2];
rz(3.0892298) q[3];
sx q[3];
rz(-3.047245) q[3];
sx q[3];
rz(0.17518763) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
