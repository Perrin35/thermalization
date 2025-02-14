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
rz(0.90647107) q[0];
sx q[0];
rz(4.8405092) q[0];
sx q[0];
rz(9.7066896) q[0];
rz(0.52892041) q[1];
sx q[1];
rz(4.79098) q[1];
sx q[1];
rz(11.00287) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0391386) q[0];
sx q[0];
rz(-2.730577) q[0];
sx q[0];
rz(-2.0356112) q[0];
rz(-2.8667502) q[2];
sx q[2];
rz(-1.5759528) q[2];
sx q[2];
rz(-2.8066563) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1100816) q[1];
sx q[1];
rz(-0.64267297) q[1];
sx q[1];
rz(2.9307632) q[1];
rz(0.36564499) q[3];
sx q[3];
rz(-2.4198101) q[3];
sx q[3];
rz(-0.66058285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6866744) q[2];
sx q[2];
rz(-0.29759559) q[2];
sx q[2];
rz(-0.61304027) q[2];
rz(0.4736627) q[3];
sx q[3];
rz(-1.9434171) q[3];
sx q[3];
rz(-1.7090428) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7365725) q[0];
sx q[0];
rz(-2.9614145) q[0];
sx q[0];
rz(0.75102425) q[0];
rz(-0.48149064) q[1];
sx q[1];
rz(-1.0844237) q[1];
sx q[1];
rz(-0.96985936) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4879726) q[0];
sx q[0];
rz(-2.3832085) q[0];
sx q[0];
rz(2.5874596) q[0];
rz(-pi) q[1];
rz(-2.9257751) q[2];
sx q[2];
rz(-1.6105798) q[2];
sx q[2];
rz(-1.7173187) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9958933) q[1];
sx q[1];
rz(-2.4384192) q[1];
sx q[1];
rz(1.165676) q[1];
rz(-1.1284351) q[3];
sx q[3];
rz(-2.1214888) q[3];
sx q[3];
rz(2.9912419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2227309) q[2];
sx q[2];
rz(-1.4758045) q[2];
sx q[2];
rz(-0.53885031) q[2];
rz(3.1239964) q[3];
sx q[3];
rz(-0.16418695) q[3];
sx q[3];
rz(2.2331451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4410412) q[0];
sx q[0];
rz(-2.7563162) q[0];
sx q[0];
rz(0.061263099) q[0];
rz(1.6368658) q[1];
sx q[1];
rz(-0.69925362) q[1];
sx q[1];
rz(-1.1963074) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.969097) q[0];
sx q[0];
rz(-0.75706702) q[0];
sx q[0];
rz(-1.9226546) q[0];
x q[1];
rz(-2.8878651) q[2];
sx q[2];
rz(-0.27058935) q[2];
sx q[2];
rz(1.9288127) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8316782) q[1];
sx q[1];
rz(-1.0161012) q[1];
sx q[1];
rz(-0.38387232) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.34516224) q[3];
sx q[3];
rz(-2.2335839) q[3];
sx q[3];
rz(-1.8927285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4521744) q[2];
sx q[2];
rz(-1.3070561) q[2];
sx q[2];
rz(0.43251953) q[2];
rz(-2.050926) q[3];
sx q[3];
rz(-0.64041036) q[3];
sx q[3];
rz(2.045491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5525621) q[0];
sx q[0];
rz(-0.93695372) q[0];
sx q[0];
rz(0.20137782) q[0];
rz(2.6248113) q[1];
sx q[1];
rz(-2.7613381) q[1];
sx q[1];
rz(0.94863272) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1051335) q[0];
sx q[0];
rz(-2.3984342) q[0];
sx q[0];
rz(-1.0418329) q[0];
rz(-pi) q[1];
rz(-2.2791512) q[2];
sx q[2];
rz(-1.7956327) q[2];
sx q[2];
rz(-2.1612957) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3624787) q[1];
sx q[1];
rz(-2.6370905) q[1];
sx q[1];
rz(0.027536784) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9734083) q[3];
sx q[3];
rz(-0.52709377) q[3];
sx q[3];
rz(0.2179365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8676694) q[2];
sx q[2];
rz(-2.7322768) q[2];
sx q[2];
rz(2.5841827) q[2];
rz(3.0607767) q[3];
sx q[3];
rz(-1.9010474) q[3];
sx q[3];
rz(1.2062937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4163365) q[0];
sx q[0];
rz(-1.4755604) q[0];
sx q[0];
rz(-1.0303372) q[0];
rz(0.15580767) q[1];
sx q[1];
rz(-1.067433) q[1];
sx q[1];
rz(0.20419289) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1847398) q[0];
sx q[0];
rz(-2.889688) q[0];
sx q[0];
rz(0.27916081) q[0];
x q[1];
rz(0.21303265) q[2];
sx q[2];
rz(-1.4430337) q[2];
sx q[2];
rz(-0.13612533) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0818644) q[1];
sx q[1];
rz(-1.7769741) q[1];
sx q[1];
rz(-0.49525201) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3685826) q[3];
sx q[3];
rz(-1.7436308) q[3];
sx q[3];
rz(2.5484249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.92945176) q[2];
sx q[2];
rz(-1.4931623) q[2];
sx q[2];
rz(0.41352752) q[2];
rz(-0.84189576) q[3];
sx q[3];
rz(-1.3827518) q[3];
sx q[3];
rz(-0.97257417) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0464762) q[0];
sx q[0];
rz(-2.0501417) q[0];
sx q[0];
rz(3.1360151) q[0];
rz(-1.7976409) q[1];
sx q[1];
rz(-2.9426212) q[1];
sx q[1];
rz(2.4406348) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0049135) q[0];
sx q[0];
rz(-2.4304996) q[0];
sx q[0];
rz(-2.2836766) q[0];
rz(-pi) q[1];
rz(-1.0089325) q[2];
sx q[2];
rz(-0.72266662) q[2];
sx q[2];
rz(0.36107963) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0030539) q[1];
sx q[1];
rz(-1.8459311) q[1];
sx q[1];
rz(1.5281926) q[1];
x q[2];
rz(-1.0526377) q[3];
sx q[3];
rz(-2.0817882) q[3];
sx q[3];
rz(0.68447733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.31412101) q[2];
sx q[2];
rz(-0.33385971) q[2];
sx q[2];
rz(-1.9742879) q[2];
rz(-2.2589034) q[3];
sx q[3];
rz(-2.3736931) q[3];
sx q[3];
rz(0.32514969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2348787) q[0];
sx q[0];
rz(-1.119708) q[0];
sx q[0];
rz(-0.085302189) q[0];
rz(-0.75434297) q[1];
sx q[1];
rz(-2.6383548) q[1];
sx q[1];
rz(1.0275966) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9011544) q[0];
sx q[0];
rz(-1.5605956) q[0];
sx q[0];
rz(-0.21252327) q[0];
x q[1];
rz(1.5761887) q[2];
sx q[2];
rz(-0.9741592) q[2];
sx q[2];
rz(0.83534681) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.69313777) q[1];
sx q[1];
rz(-1.4488646) q[1];
sx q[1];
rz(-1.2454525) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.87966921) q[3];
sx q[3];
rz(-2.4604049) q[3];
sx q[3];
rz(-2.3613514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9895642) q[2];
sx q[2];
rz(-2.1542408) q[2];
sx q[2];
rz(2.7222471) q[2];
rz(-2.5069405) q[3];
sx q[3];
rz(-2.3366163) q[3];
sx q[3];
rz(1.1254719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(0.80970508) q[0];
sx q[0];
rz(-2.2655847) q[0];
sx q[0];
rz(-2.4494655) q[0];
rz(-0.57811111) q[1];
sx q[1];
rz(-2.7222241) q[1];
sx q[1];
rz(-1.7587761) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7958) q[0];
sx q[0];
rz(-1.6383071) q[0];
sx q[0];
rz(0.019492143) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2302551) q[2];
sx q[2];
rz(-2.132086) q[2];
sx q[2];
rz(2.5709465) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.67685157) q[1];
sx q[1];
rz(-0.60982043) q[1];
sx q[1];
rz(-1.9863542) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.47494294) q[3];
sx q[3];
rz(-0.90427665) q[3];
sx q[3];
rz(2.8996244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77070037) q[2];
sx q[2];
rz(-2.6110677) q[2];
sx q[2];
rz(-1.4250866) q[2];
rz(-2.6254081) q[3];
sx q[3];
rz(-0.97847146) q[3];
sx q[3];
rz(-1.9019351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4952963) q[0];
sx q[0];
rz(-1.7072059) q[0];
sx q[0];
rz(-2.7700951) q[0];
rz(-2.3067572) q[1];
sx q[1];
rz(-0.8131665) q[1];
sx q[1];
rz(3.1239948) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.427582) q[0];
sx q[0];
rz(-2.0907864) q[0];
sx q[0];
rz(0.31975694) q[0];
x q[1];
rz(1.8210635) q[2];
sx q[2];
rz(-2.2907718) q[2];
sx q[2];
rz(-1.7257476) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7323276) q[1];
sx q[1];
rz(-1.0068934) q[1];
sx q[1];
rz(1.1281518) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5334355) q[3];
sx q[3];
rz(-1.3038692) q[3];
sx q[3];
rz(-1.1692695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1060433) q[2];
sx q[2];
rz(-2.273166) q[2];
sx q[2];
rz(-3.0412728) q[2];
rz(-0.22942461) q[3];
sx q[3];
rz(-2.094163) q[3];
sx q[3];
rz(-2.6917698) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68596524) q[0];
sx q[0];
rz(-2.8534511) q[0];
sx q[0];
rz(0.57383865) q[0];
rz(-2.0876743) q[1];
sx q[1];
rz(-2.0951447) q[1];
sx q[1];
rz(0.67646772) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94714245) q[0];
sx q[0];
rz(-1.3279755) q[0];
sx q[0];
rz(0.28698289) q[0];
rz(-pi) q[1];
rz(-2.8847938) q[2];
sx q[2];
rz(-1.1780103) q[2];
sx q[2];
rz(-1.8671672) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.45470787) q[1];
sx q[1];
rz(-2.096972) q[1];
sx q[1];
rz(-1.4684907) q[1];
x q[2];
rz(-0.19837899) q[3];
sx q[3];
rz(-1.0470445) q[3];
sx q[3];
rz(-2.3959121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.68710589) q[2];
sx q[2];
rz(-2.3837619) q[2];
sx q[2];
rz(-1.6472316) q[2];
rz(2.2729661) q[3];
sx q[3];
rz(-2.4354911) q[3];
sx q[3];
rz(2.6715265) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0710707) q[0];
sx q[0];
rz(-1.1916397) q[0];
sx q[0];
rz(1.2217039) q[0];
rz(1.8078049) q[1];
sx q[1];
rz(-1.8232657) q[1];
sx q[1];
rz(2.5008536) q[1];
rz(-3.0467544) q[2];
sx q[2];
rz(-0.76764501) q[2];
sx q[2];
rz(-2.1697247) q[2];
rz(1.815769) q[3];
sx q[3];
rz(-0.180937) q[3];
sx q[3];
rz(-2.141249) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
