OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2090981) q[0];
sx q[0];
rz(-3.0734835) q[0];
sx q[0];
rz(-0.77686247) q[0];
rz(-4.4486899) q[1];
sx q[1];
rz(-2.1721462) q[1];
sx q[1];
rz(8.1028508) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7711219) q[0];
sx q[0];
rz(-1.4323455) q[0];
sx q[0];
rz(2.2884503) q[0];
rz(-pi) q[1];
rz(2.6451557) q[2];
sx q[2];
rz(-1.264071) q[2];
sx q[2];
rz(-2.6076743) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4522534) q[1];
sx q[1];
rz(-1.6000556) q[1];
sx q[1];
rz(-1.9797146) q[1];
rz(1.8613937) q[3];
sx q[3];
rz(-1.7589671) q[3];
sx q[3];
rz(-0.0083323697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0778568) q[2];
sx q[2];
rz(-1.7813762) q[2];
sx q[2];
rz(-2.784101) q[2];
rz(-2.5743971) q[3];
sx q[3];
rz(-1.4128069) q[3];
sx q[3];
rz(2.7020057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1883063) q[0];
sx q[0];
rz(-2.8585241) q[0];
sx q[0];
rz(1.333492) q[0];
rz(-0.4370583) q[1];
sx q[1];
rz(-0.60085618) q[1];
sx q[1];
rz(-2.3016047) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8020426) q[0];
sx q[0];
rz(-1.7239445) q[0];
sx q[0];
rz(-2.0846372) q[0];
x q[1];
rz(2.5469668) q[2];
sx q[2];
rz(-1.42808) q[2];
sx q[2];
rz(3.0171301) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6691554) q[1];
sx q[1];
rz(-2.1320813) q[1];
sx q[1];
rz(-0.88772933) q[1];
x q[2];
rz(0.038957314) q[3];
sx q[3];
rz(-2.1202728) q[3];
sx q[3];
rz(-0.4645068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4240894) q[2];
sx q[2];
rz(-1.4428029) q[2];
sx q[2];
rz(-1.7355512) q[2];
rz(2.9794335) q[3];
sx q[3];
rz(-1.5806961) q[3];
sx q[3];
rz(-0.46419188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9185987) q[0];
sx q[0];
rz(-0.080568947) q[0];
sx q[0];
rz(2.719847) q[0];
rz(2.2987507) q[1];
sx q[1];
rz(-1.3858567) q[1];
sx q[1];
rz(-2.8657894) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11711794) q[0];
sx q[0];
rz(-1.3399276) q[0];
sx q[0];
rz(0.33322115) q[0];
x q[1];
rz(-2.6709072) q[2];
sx q[2];
rz(-1.7690725) q[2];
sx q[2];
rz(-2.7338365) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2308064) q[1];
sx q[1];
rz(-1.2503997) q[1];
sx q[1];
rz(1.6256534) q[1];
x q[2];
rz(-2.2725676) q[3];
sx q[3];
rz(-1.6459673) q[3];
sx q[3];
rz(1.4787256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3383823) q[2];
sx q[2];
rz(-0.79757491) q[2];
sx q[2];
rz(-2.4513643) q[2];
rz(0.35987443) q[3];
sx q[3];
rz(-1.7706324) q[3];
sx q[3];
rz(2.7580822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2279219) q[0];
sx q[0];
rz(-2.0991195) q[0];
sx q[0];
rz(-0.87164718) q[0];
rz(0.40920416) q[1];
sx q[1];
rz(-2.6458461) q[1];
sx q[1];
rz(-0.31235487) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.271628) q[0];
sx q[0];
rz(-1.6883435) q[0];
sx q[0];
rz(0.44674504) q[0];
rz(2.9351685) q[2];
sx q[2];
rz(-1.6684434) q[2];
sx q[2];
rz(-2.4149778) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.67955454) q[1];
sx q[1];
rz(-1.3674539) q[1];
sx q[1];
rz(0.36851369) q[1];
x q[2];
rz(-0.64849017) q[3];
sx q[3];
rz(-1.7269508) q[3];
sx q[3];
rz(0.79532571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.006134) q[2];
sx q[2];
rz(-1.4486855) q[2];
sx q[2];
rz(0.82497605) q[2];
rz(-1.4247591) q[3];
sx q[3];
rz(-0.32130876) q[3];
sx q[3];
rz(0.43535522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7216126) q[0];
sx q[0];
rz(-1.551832) q[0];
sx q[0];
rz(-2.2671674) q[0];
rz(-2.8127316) q[1];
sx q[1];
rz(-1.3855653) q[1];
sx q[1];
rz(2.0430476) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6499574) q[0];
sx q[0];
rz(-1.1960317) q[0];
sx q[0];
rz(1.6894421) q[0];
rz(-0.83580612) q[2];
sx q[2];
rz(-1.5975738) q[2];
sx q[2];
rz(1.4876897) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.85531536) q[1];
sx q[1];
rz(-2.554727) q[1];
sx q[1];
rz(-1.0081968) q[1];
rz(-pi) q[2];
rz(0.9807113) q[3];
sx q[3];
rz(-1.3385337) q[3];
sx q[3];
rz(-1.4651608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.141779) q[2];
sx q[2];
rz(-1.1835316) q[2];
sx q[2];
rz(-0.67406526) q[2];
rz(1.6541803) q[3];
sx q[3];
rz(-1.4451278) q[3];
sx q[3];
rz(2.8580247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0312626) q[0];
sx q[0];
rz(-1.0227579) q[0];
sx q[0];
rz(1.3856101) q[0];
rz(-2.022187) q[1];
sx q[1];
rz(-1.6740572) q[1];
sx q[1];
rz(-1.3853692) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7840609) q[0];
sx q[0];
rz(-0.94793944) q[0];
sx q[0];
rz(0.31185598) q[0];
x q[1];
rz(2.643061) q[2];
sx q[2];
rz(-1.6462925) q[2];
sx q[2];
rz(2.1488862) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.54810235) q[1];
sx q[1];
rz(-1.2498651) q[1];
sx q[1];
rz(-1.3671419) q[1];
rz(-pi) q[2];
rz(1.5343666) q[3];
sx q[3];
rz(-1.4163989) q[3];
sx q[3];
rz(1.2702219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0783656) q[2];
sx q[2];
rz(-2.207022) q[2];
sx q[2];
rz(-2.7654977) q[2];
rz(-1.1345351) q[3];
sx q[3];
rz(-2.7781656) q[3];
sx q[3];
rz(-2.334972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1217839) q[0];
sx q[0];
rz(-2.5406295) q[0];
sx q[0];
rz(3.0390749) q[0];
rz(-2.5566697) q[1];
sx q[1];
rz(-2.1939317) q[1];
sx q[1];
rz(-1.9580511) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69307704) q[0];
sx q[0];
rz(-1.428837) q[0];
sx q[0];
rz(-1.4407115) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41924119) q[2];
sx q[2];
rz(-0.70690522) q[2];
sx q[2];
rz(-0.32419328) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8475323) q[1];
sx q[1];
rz(-0.69662635) q[1];
sx q[1];
rz(1.1678334) q[1];
x q[2];
rz(-1.8374341) q[3];
sx q[3];
rz(-2.023306) q[3];
sx q[3];
rz(0.072991144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3465603) q[2];
sx q[2];
rz(-0.9298032) q[2];
sx q[2];
rz(0.19980508) q[2];
rz(-2.0810769) q[3];
sx q[3];
rz(-0.54148713) q[3];
sx q[3];
rz(-3.0919891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26315966) q[0];
sx q[0];
rz(-1.4819772) q[0];
sx q[0];
rz(-0.31295452) q[0];
rz(2.1921659) q[1];
sx q[1];
rz(-1.7890309) q[1];
sx q[1];
rz(1.4535646) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0003371) q[0];
sx q[0];
rz(-2.9337058) q[0];
sx q[0];
rz(-0.45335404) q[0];
rz(-pi) q[1];
rz(-1.6498927) q[2];
sx q[2];
rz(-0.95273113) q[2];
sx q[2];
rz(-1.5072418) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1968699) q[1];
sx q[1];
rz(-2.2579282) q[1];
sx q[1];
rz(1.8623167) q[1];
x q[2];
rz(1.231915) q[3];
sx q[3];
rz(-2.0032855) q[3];
sx q[3];
rz(-0.56604715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.35385418) q[2];
sx q[2];
rz(-2.165803) q[2];
sx q[2];
rz(1.3346416) q[2];
rz(-1.3245964) q[3];
sx q[3];
rz(-0.66671222) q[3];
sx q[3];
rz(-1.9273531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.036309328) q[0];
sx q[0];
rz(-2.5258625) q[0];
sx q[0];
rz(-0.69806725) q[0];
rz(-0.87567323) q[1];
sx q[1];
rz(-2.0798637) q[1];
sx q[1];
rz(2.7811513) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1669641) q[0];
sx q[0];
rz(-0.95050838) q[0];
sx q[0];
rz(1.0247158) q[0];
rz(-1.5064042) q[2];
sx q[2];
rz(-1.1887822) q[2];
sx q[2];
rz(-0.40859336) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7061758) q[1];
sx q[1];
rz(-2.1565394) q[1];
sx q[1];
rz(-0.46153586) q[1];
rz(-pi) q[2];
rz(-1.7385941) q[3];
sx q[3];
rz(-1.0281963) q[3];
sx q[3];
rz(-2.8436273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.91161072) q[2];
sx q[2];
rz(-1.4260099) q[2];
sx q[2];
rz(-0.81725517) q[2];
rz(0.83003712) q[3];
sx q[3];
rz(-2.6006112) q[3];
sx q[3];
rz(0.219492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.5068186) q[0];
sx q[0];
rz(-2.2932678) q[0];
sx q[0];
rz(0.032055227) q[0];
rz(1.8219148) q[1];
sx q[1];
rz(-1.3751043) q[1];
sx q[1];
rz(0.65418902) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7830156) q[0];
sx q[0];
rz(-0.23046432) q[0];
sx q[0];
rz(2.5925267) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0243353) q[2];
sx q[2];
rz(-2.3058866) q[2];
sx q[2];
rz(-1.0356366) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5947551) q[1];
sx q[1];
rz(-2.0528626) q[1];
sx q[1];
rz(-2.099311) q[1];
rz(-pi) q[2];
rz(-1.9496253) q[3];
sx q[3];
rz(-2.4740268) q[3];
sx q[3];
rz(0.51578705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0111982) q[2];
sx q[2];
rz(-2.0989213) q[2];
sx q[2];
rz(2.1350258) q[2];
rz(-0.69342962) q[3];
sx q[3];
rz(-1.3207685) q[3];
sx q[3];
rz(1.8600195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4257767) q[0];
sx q[0];
rz(-2.4632813) q[0];
sx q[0];
rz(3.0671469) q[0];
rz(2.2609932) q[1];
sx q[1];
rz(-1.1136628) q[1];
sx q[1];
rz(-2.640092) q[1];
rz(2.1660317) q[2];
sx q[2];
rz(-2.3481993) q[2];
sx q[2];
rz(1.3994155) q[2];
rz(-1.1893336) q[3];
sx q[3];
rz(-2.2305924) q[3];
sx q[3];
rz(2.3371405) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
