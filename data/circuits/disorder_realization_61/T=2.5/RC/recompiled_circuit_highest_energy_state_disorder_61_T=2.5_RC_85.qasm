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
rz(-2.6391368) q[0];
sx q[0];
rz(-1.1075736) q[0];
sx q[0];
rz(1.7662319) q[0];
rz(2.8227168) q[1];
sx q[1];
rz(-1.5006637) q[1];
sx q[1];
rz(-0.73135102) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99217447) q[0];
sx q[0];
rz(-0.96620027) q[0];
sx q[0];
rz(0.3913619) q[0];
x q[1];
rz(-0.0053623036) q[2];
sx q[2];
rz(-0.66197936) q[2];
sx q[2];
rz(2.299451) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9125376) q[1];
sx q[1];
rz(-0.78733912) q[1];
sx q[1];
rz(-0.064219193) q[1];
rz(2.6177789) q[3];
sx q[3];
rz(-0.18693811) q[3];
sx q[3];
rz(0.65526774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2607164) q[2];
sx q[2];
rz(-1.9274351) q[2];
sx q[2];
rz(-0.94240776) q[2];
rz(-2.2830394) q[3];
sx q[3];
rz(-0.61045727) q[3];
sx q[3];
rz(0.71949351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.079916) q[0];
sx q[0];
rz(-0.55773568) q[0];
sx q[0];
rz(0.058636531) q[0];
rz(-0.05642852) q[1];
sx q[1];
rz(-1.3899048) q[1];
sx q[1];
rz(2.0243534) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53848828) q[0];
sx q[0];
rz(-1.9637917) q[0];
sx q[0];
rz(0.80727838) q[0];
rz(-pi) q[1];
rz(-1.2710167) q[2];
sx q[2];
rz(-1.4477056) q[2];
sx q[2];
rz(-1.276615) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8694508) q[1];
sx q[1];
rz(-1.488416) q[1];
sx q[1];
rz(-1.3314962) q[1];
rz(-1.6916709) q[3];
sx q[3];
rz(-2.0434922) q[3];
sx q[3];
rz(2.8877088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3650018) q[2];
sx q[2];
rz(-2.8114909) q[2];
sx q[2];
rz(-0.4064202) q[2];
rz(0.31600076) q[3];
sx q[3];
rz(-1.6812811) q[3];
sx q[3];
rz(-2.3587904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2773748) q[0];
sx q[0];
rz(-0.45806956) q[0];
sx q[0];
rz(1.3360485) q[0];
rz(-0.27851963) q[1];
sx q[1];
rz(-1.3986992) q[1];
sx q[1];
rz(-2.1727402) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8271885) q[0];
sx q[0];
rz(-0.97100706) q[0];
sx q[0];
rz(-1.1161854) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1325705) q[2];
sx q[2];
rz(-2.2780905) q[2];
sx q[2];
rz(2.2826441) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4668256) q[1];
sx q[1];
rz(-1.5286501) q[1];
sx q[1];
rz(-1.274363) q[1];
rz(-pi) q[2];
rz(-1.3125456) q[3];
sx q[3];
rz(-0.86229339) q[3];
sx q[3];
rz(1.4227305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.99732533) q[2];
sx q[2];
rz(-2.6582025) q[2];
sx q[2];
rz(0.17511314) q[2];
rz(-1.8363606) q[3];
sx q[3];
rz(-2.1988018) q[3];
sx q[3];
rz(1.2828627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99821943) q[0];
sx q[0];
rz(-1.136919) q[0];
sx q[0];
rz(-2.112222) q[0];
rz(-2.3221305) q[1];
sx q[1];
rz(-1.9792604) q[1];
sx q[1];
rz(2.3447461) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0734755) q[0];
sx q[0];
rz(-1.1321018) q[0];
sx q[0];
rz(1.7832157) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0060002319) q[2];
sx q[2];
rz(-1.0428535) q[2];
sx q[2];
rz(-2.1354298) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0174344) q[1];
sx q[1];
rz(-0.34338152) q[1];
sx q[1];
rz(2.2209973) q[1];
rz(-2.6096973) q[3];
sx q[3];
rz(-0.73557702) q[3];
sx q[3];
rz(0.07594219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0113819) q[2];
sx q[2];
rz(-3.104976) q[2];
sx q[2];
rz(-1.2288564) q[2];
rz(-0.43904385) q[3];
sx q[3];
rz(-1.565758) q[3];
sx q[3];
rz(1.9108093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48290408) q[0];
sx q[0];
rz(-1.6084325) q[0];
sx q[0];
rz(0.80832344) q[0];
rz(1.3574379) q[1];
sx q[1];
rz(-2.3517377) q[1];
sx q[1];
rz(-2.435991) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7019136) q[0];
sx q[0];
rz(-1.2703623) q[0];
sx q[0];
rz(-1.6070532) q[0];
rz(-pi) q[1];
rz(-0.23106261) q[2];
sx q[2];
rz(-2.2602644) q[2];
sx q[2];
rz(-1.5901515) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8568154) q[1];
sx q[1];
rz(-1.5734623) q[1];
sx q[1];
rz(-1.8239914) q[1];
rz(-pi) q[2];
rz(-2.5217149) q[3];
sx q[3];
rz(-2.2408463) q[3];
sx q[3];
rz(-1.0244964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.48050532) q[2];
sx q[2];
rz(-1.8578119) q[2];
sx q[2];
rz(0.60574469) q[2];
rz(-3.0211966) q[3];
sx q[3];
rz(-1.8431289) q[3];
sx q[3];
rz(-0.81407434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26009387) q[0];
sx q[0];
rz(-2.3277178) q[0];
sx q[0];
rz(0.0023728097) q[0];
rz(-0.35898769) q[1];
sx q[1];
rz(-1.2009883) q[1];
sx q[1];
rz(-0.40828362) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15015442) q[0];
sx q[0];
rz(-1.5350193) q[0];
sx q[0];
rz(1.4442025) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9469366) q[2];
sx q[2];
rz(-2.500211) q[2];
sx q[2];
rz(2.1958786) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3107016) q[1];
sx q[1];
rz(-2.9017555) q[1];
sx q[1];
rz(2.8562921) q[1];
x q[2];
rz(2.5465132) q[3];
sx q[3];
rz(-2.236955) q[3];
sx q[3];
rz(0.72197589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6458873) q[2];
sx q[2];
rz(-2.3453823) q[2];
sx q[2];
rz(0.1296981) q[2];
rz(0.4194704) q[3];
sx q[3];
rz(-1.2129815) q[3];
sx q[3];
rz(0.81370846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7523338) q[0];
sx q[0];
rz(-0.77521721) q[0];
sx q[0];
rz(0.68354052) q[0];
rz(0.93841249) q[1];
sx q[1];
rz(-1.3886195) q[1];
sx q[1];
rz(1.9155115) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3049604) q[0];
sx q[0];
rz(-1.9373405) q[0];
sx q[0];
rz(-0.62303294) q[0];
x q[1];
rz(-1.7581045) q[2];
sx q[2];
rz(-2.1626327) q[2];
sx q[2];
rz(0.96697742) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9270806) q[1];
sx q[1];
rz(-2.0190372) q[1];
sx q[1];
rz(0.31772504) q[1];
rz(2.5874373) q[3];
sx q[3];
rz(-0.99923493) q[3];
sx q[3];
rz(0.038427834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.73416) q[2];
sx q[2];
rz(-1.4340883) q[2];
sx q[2];
rz(-2.4343991) q[2];
rz(1.6678984) q[3];
sx q[3];
rz(-2.4652822) q[3];
sx q[3];
rz(-1.7128806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68607512) q[0];
sx q[0];
rz(-2.7954743) q[0];
sx q[0];
rz(0.92380512) q[0];
rz(-1.0651945) q[1];
sx q[1];
rz(-1.9520452) q[1];
sx q[1];
rz(-2.0069897) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3649639) q[0];
sx q[0];
rz(-2.7641649) q[0];
sx q[0];
rz(1.6799742) q[0];
rz(0.61750268) q[2];
sx q[2];
rz(-0.98101014) q[2];
sx q[2];
rz(1.2794354) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2345786) q[1];
sx q[1];
rz(-1.5173787) q[1];
sx q[1];
rz(2.1018886) q[1];
x q[2];
rz(1.1884965) q[3];
sx q[3];
rz(-0.99671066) q[3];
sx q[3];
rz(-2.7758383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.137546) q[2];
sx q[2];
rz(-1.4038439) q[2];
sx q[2];
rz(-0.64296067) q[2];
rz(-0.17063394) q[3];
sx q[3];
rz(-0.65244397) q[3];
sx q[3];
rz(-2.046106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3951025) q[0];
sx q[0];
rz(-3.009142) q[0];
sx q[0];
rz(0.80628959) q[0];
rz(0.0099446615) q[1];
sx q[1];
rz(-2.5237623) q[1];
sx q[1];
rz(-0.26201216) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.952103) q[0];
sx q[0];
rz(-1.4052009) q[0];
sx q[0];
rz(-2.1013592) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9295983) q[2];
sx q[2];
rz(-2.4851481) q[2];
sx q[2];
rz(-0.5613297) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1478473) q[1];
sx q[1];
rz(-0.48300749) q[1];
sx q[1];
rz(2.0266533) q[1];
rz(-pi) q[2];
rz(-0.32829653) q[3];
sx q[3];
rz(-2.0496164) q[3];
sx q[3];
rz(-1.6678068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.749873) q[2];
sx q[2];
rz(-2.5842857) q[2];
sx q[2];
rz(2.9987175) q[2];
rz(-2.2376132) q[3];
sx q[3];
rz(-1.6429792) q[3];
sx q[3];
rz(2.2601295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7664465) q[0];
sx q[0];
rz(-1.2297933) q[0];
sx q[0];
rz(-0.65650702) q[0];
rz(-0.87908602) q[1];
sx q[1];
rz(-1.1708941) q[1];
sx q[1];
rz(0.62969977) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64163369) q[0];
sx q[0];
rz(-2.2370506) q[0];
sx q[0];
rz(1.0532265) q[0];
rz(-2.1037322) q[2];
sx q[2];
rz(-1.3383972) q[2];
sx q[2];
rz(2.971422) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.15805298) q[1];
sx q[1];
rz(-1/(6*pi)) q[1];
sx q[1];
rz(-1.0178465) q[1];
rz(-1.042243) q[3];
sx q[3];
rz(-1.8404598) q[3];
sx q[3];
rz(0.93664133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5188344) q[2];
sx q[2];
rz(-0.82509416) q[2];
sx q[2];
rz(0.68010124) q[2];
rz(-2.425219) q[3];
sx q[3];
rz(-3.0383737) q[3];
sx q[3];
rz(-1.1187925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5685365) q[0];
sx q[0];
rz(-1.2328883) q[0];
sx q[0];
rz(0.20986025) q[0];
rz(-2.9987891) q[1];
sx q[1];
rz(-1.0719943) q[1];
sx q[1];
rz(0.73678585) q[1];
rz(1.8129195) q[2];
sx q[2];
rz(-1.5151146) q[2];
sx q[2];
rz(1.7641868) q[2];
rz(-1.5132001) q[3];
sx q[3];
rz(-2.5313898) q[3];
sx q[3];
rz(1.678086) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
