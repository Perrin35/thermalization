OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.39188448) q[0];
sx q[0];
rz(-0.19667974) q[0];
sx q[0];
rz(-1.952202) q[0];
rz(0.2285129) q[1];
sx q[1];
rz(-0.84140468) q[1];
sx q[1];
rz(-2.7639311) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084378622) q[0];
sx q[0];
rz(-2.9998261) q[0];
sx q[0];
rz(-1.0300107) q[0];
x q[1];
rz(1.471465) q[2];
sx q[2];
rz(-2.8566395) q[2];
sx q[2];
rz(3.1379267) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9610112) q[1];
sx q[1];
rz(-0.86997021) q[1];
sx q[1];
rz(0.58971528) q[1];
x q[2];
rz(0.41892003) q[3];
sx q[3];
rz(-1.6994611) q[3];
sx q[3];
rz(-0.045623771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1315786) q[2];
sx q[2];
rz(-2.6476314) q[2];
sx q[2];
rz(-2.4689891) q[2];
rz(0.16942313) q[3];
sx q[3];
rz(-0.38893458) q[3];
sx q[3];
rz(-1.8030362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.910903) q[0];
sx q[0];
rz(-0.90086532) q[0];
sx q[0];
rz(0.22856523) q[0];
rz(2.9810492) q[1];
sx q[1];
rz(-1.7030145) q[1];
sx q[1];
rz(-0.28796089) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9790968) q[0];
sx q[0];
rz(-1.393035) q[0];
sx q[0];
rz(-1.3773247) q[0];
rz(-pi) q[1];
rz(-0.9777074) q[2];
sx q[2];
rz(-1.2427254) q[2];
sx q[2];
rz(-1.6239945) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2634695) q[1];
sx q[1];
rz(-0.97349226) q[1];
sx q[1];
rz(1.0695446) q[1];
rz(0.1045707) q[3];
sx q[3];
rz(-2.7894756) q[3];
sx q[3];
rz(0.62192384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5471389) q[2];
sx q[2];
rz(-1.889166) q[2];
sx q[2];
rz(-1.4734369) q[2];
rz(-2.2041221) q[3];
sx q[3];
rz(-2.6963186) q[3];
sx q[3];
rz(-2.766585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32524747) q[0];
sx q[0];
rz(-2.6419817) q[0];
sx q[0];
rz(2.8741799) q[0];
rz(-1.422241) q[1];
sx q[1];
rz(-1.1317252) q[1];
sx q[1];
rz(2.1898988) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5121582) q[0];
sx q[0];
rz(-0.21244563) q[0];
sx q[0];
rz(-2.0149219) q[0];
rz(-1.0551664) q[2];
sx q[2];
rz(-1.3736758) q[2];
sx q[2];
rz(-0.65161639) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.051085) q[1];
sx q[1];
rz(-1.8790434) q[1];
sx q[1];
rz(-0.80866637) q[1];
x q[2];
rz(-0.88055196) q[3];
sx q[3];
rz(-0.75062245) q[3];
sx q[3];
rz(-0.39117884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0114228) q[2];
sx q[2];
rz(-1.495785) q[2];
sx q[2];
rz(-0.34417957) q[2];
rz(0.88642818) q[3];
sx q[3];
rz(-2.8635946) q[3];
sx q[3];
rz(-0.56604958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13720559) q[0];
sx q[0];
rz(-1.6442278) q[0];
sx q[0];
rz(-1.244506) q[0];
rz(-3.0124774) q[1];
sx q[1];
rz(-1.3543509) q[1];
sx q[1];
rz(-0.37277645) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.059797212) q[0];
sx q[0];
rz(-2.3764388) q[0];
sx q[0];
rz(2.9761936) q[0];
rz(1.1431085) q[2];
sx q[2];
rz(-2.3017985) q[2];
sx q[2];
rz(-2.1172303) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8744295) q[1];
sx q[1];
rz(-1.0871372) q[1];
sx q[1];
rz(-2.8664385) q[1];
x q[2];
rz(1.0944486) q[3];
sx q[3];
rz(-1.2762478) q[3];
sx q[3];
rz(-2.0457207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.63561511) q[2];
sx q[2];
rz(-1.4900692) q[2];
sx q[2];
rz(-2.2367031) q[2];
rz(2.7010226) q[3];
sx q[3];
rz(-2.6456656) q[3];
sx q[3];
rz(2.1499965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47914094) q[0];
sx q[0];
rz(-0.15563706) q[0];
sx q[0];
rz(-1.8537846) q[0];
rz(1.6632535) q[1];
sx q[1];
rz(-1.9961424) q[1];
sx q[1];
rz(-0.53422654) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9984765) q[0];
sx q[0];
rz(-1.8792218) q[0];
sx q[0];
rz(0.31750676) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0787557) q[2];
sx q[2];
rz(-1.2675708) q[2];
sx q[2];
rz(0.505503) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.91671645) q[1];
sx q[1];
rz(-1.4955048) q[1];
sx q[1];
rz(0.063898357) q[1];
rz(-pi) q[2];
rz(0.71877919) q[3];
sx q[3];
rz(-1.6703509) q[3];
sx q[3];
rz(1.0574785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5010895) q[2];
sx q[2];
rz(-1.9260294) q[2];
sx q[2];
rz(-0.14349288) q[2];
rz(1.7701373) q[3];
sx q[3];
rz(-1.2797132) q[3];
sx q[3];
rz(-2.9218856) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7025529) q[0];
sx q[0];
rz(-2.2655903) q[0];
sx q[0];
rz(-0.23705661) q[0];
rz(-1.2409695) q[1];
sx q[1];
rz(-0.82890141) q[1];
sx q[1];
rz(-1.3751078) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9000589) q[0];
sx q[0];
rz(-2.7068479) q[0];
sx q[0];
rz(0.66381201) q[0];
rz(-pi) q[1];
rz(3.0639406) q[2];
sx q[2];
rz(-2.3911871) q[2];
sx q[2];
rz(0.99496182) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3875229) q[1];
sx q[1];
rz(-1.9807528) q[1];
sx q[1];
rz(1.5441262) q[1];
rz(-pi) q[2];
rz(2.6508413) q[3];
sx q[3];
rz(-1.9487582) q[3];
sx q[3];
rz(2.2839387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6087626) q[2];
sx q[2];
rz(-0.71981788) q[2];
sx q[2];
rz(-0.14870816) q[2];
rz(-3.1249629) q[3];
sx q[3];
rz(-2.7745268) q[3];
sx q[3];
rz(-0.19255157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49611133) q[0];
sx q[0];
rz(-0.2922903) q[0];
sx q[0];
rz(-0.87316978) q[0];
rz(2.219615) q[1];
sx q[1];
rz(-2.0261804) q[1];
sx q[1];
rz(-0.08392863) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8221995) q[0];
sx q[0];
rz(-2.5376352) q[0];
sx q[0];
rz(-1.1711867) q[0];
x q[1];
rz(0.35328816) q[2];
sx q[2];
rz(-1.366426) q[2];
sx q[2];
rz(-2.3452961) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2578656) q[1];
sx q[1];
rz(-2.1818433) q[1];
sx q[1];
rz(2.3962767) q[1];
rz(3.0725669) q[3];
sx q[3];
rz(-2.2996174) q[3];
sx q[3];
rz(-0.1544827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.38368791) q[2];
sx q[2];
rz(-0.44054511) q[2];
sx q[2];
rz(2.596358) q[2];
rz(0.43045726) q[3];
sx q[3];
rz(-2.0038219) q[3];
sx q[3];
rz(-2.8366413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51171821) q[0];
sx q[0];
rz(-0.015462333) q[0];
sx q[0];
rz(1.9301201) q[0];
rz(-0.97310549) q[1];
sx q[1];
rz(-2.548023) q[1];
sx q[1];
rz(2.1957695) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66447542) q[0];
sx q[0];
rz(-0.23453377) q[0];
sx q[0];
rz(-2.7995336) q[0];
rz(0.81994762) q[2];
sx q[2];
rz(-0.87265271) q[2];
sx q[2];
rz(-1.9422216) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9981873) q[1];
sx q[1];
rz(-2.7618976) q[1];
sx q[1];
rz(0.077296301) q[1];
x q[2];
rz(2.4754727) q[3];
sx q[3];
rz(-1.5770116) q[3];
sx q[3];
rz(-0.70762779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2716081) q[2];
sx q[2];
rz(-1.4166069) q[2];
sx q[2];
rz(-2.8611709) q[2];
rz(-2.5366606) q[3];
sx q[3];
rz(-1.0792024) q[3];
sx q[3];
rz(-2.159507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9873001) q[0];
sx q[0];
rz(-0.91006088) q[0];
sx q[0];
rz(0.61532414) q[0];
rz(-2.212021) q[1];
sx q[1];
rz(-2.2832182) q[1];
sx q[1];
rz(-2.5659134) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3424123) q[0];
sx q[0];
rz(-0.77501446) q[0];
sx q[0];
rz(-2.4393625) q[0];
rz(-pi) q[1];
rz(-1.0506389) q[2];
sx q[2];
rz(-0.22503223) q[2];
sx q[2];
rz(2.2573543) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5513735) q[1];
sx q[1];
rz(-3.0312523) q[1];
sx q[1];
rz(1.8054086) q[1];
rz(-pi) q[2];
rz(-1.8768164) q[3];
sx q[3];
rz(-0.56193202) q[3];
sx q[3];
rz(-2.7276873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.82970396) q[2];
sx q[2];
rz(-0.76278937) q[2];
sx q[2];
rz(-3.1196307) q[2];
rz(-2.9711376) q[3];
sx q[3];
rz(-2.1178092) q[3];
sx q[3];
rz(-2.8767265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4158674) q[0];
sx q[0];
rz(-2.2150345) q[0];
sx q[0];
rz(2.5073994) q[0];
rz(-3.1126853) q[1];
sx q[1];
rz(-2.3529265) q[1];
sx q[1];
rz(-0.96910563) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4934851) q[0];
sx q[0];
rz(-3.0457975) q[0];
sx q[0];
rz(-2.2936054) q[0];
rz(-pi) q[1];
rz(1*pi/15) q[2];
sx q[2];
rz(-0.75762118) q[2];
sx q[2];
rz(0.19291887) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9173968) q[1];
sx q[1];
rz(-2.2317413) q[1];
sx q[1];
rz(-1.7978653) q[1];
rz(0.44090791) q[3];
sx q[3];
rz(-1.3135859) q[3];
sx q[3];
rz(-0.5514901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4518296) q[2];
sx q[2];
rz(-1.4844866) q[2];
sx q[2];
rz(2.7815212) q[2];
rz(-1.6137971) q[3];
sx q[3];
rz(-0.66643047) q[3];
sx q[3];
rz(2.578919) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9344899) q[0];
sx q[0];
rz(-1.5705382) q[0];
sx q[0];
rz(-1.6194153) q[0];
rz(-0.044152505) q[1];
sx q[1];
rz(-1.4587198) q[1];
sx q[1];
rz(-1.1062467) q[1];
rz(-0.33640484) q[2];
sx q[2];
rz(-1.207926) q[2];
sx q[2];
rz(1.5628857) q[2];
rz(2.5526657) q[3];
sx q[3];
rz(-2.6116461) q[3];
sx q[3];
rz(-2.6245821) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
