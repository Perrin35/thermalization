OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5623915) q[0];
sx q[0];
rz(-2.3810823) q[0];
sx q[0];
rz(-1.0150681) q[0];
rz(1.9782344) q[1];
sx q[1];
rz(-0.42127633) q[1];
sx q[1];
rz(1.2383229) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2531812) q[0];
sx q[0];
rz(-2.0813353) q[0];
sx q[0];
rz(-2.2967413) q[0];
rz(-0.36739393) q[2];
sx q[2];
rz(-1.1720554) q[2];
sx q[2];
rz(-0.90786394) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0933669) q[1];
sx q[1];
rz(-1.9328572) q[1];
sx q[1];
rz(1.0735379) q[1];
rz(0.26217272) q[3];
sx q[3];
rz(-0.57973639) q[3];
sx q[3];
rz(-1.4822823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6282661) q[2];
sx q[2];
rz(-1.790975) q[2];
sx q[2];
rz(-0.70293054) q[2];
rz(-2.2859196) q[3];
sx q[3];
rz(-0.84508768) q[3];
sx q[3];
rz(1.8446911) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4502451) q[0];
sx q[0];
rz(-1.0991199) q[0];
sx q[0];
rz(-2.3968089) q[0];
rz(-1.5738457) q[1];
sx q[1];
rz(-0.66190043) q[1];
sx q[1];
rz(-0.94211284) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73977913) q[0];
sx q[0];
rz(-1.7727858) q[0];
sx q[0];
rz(-2.3752579) q[0];
rz(-1.6555384) q[2];
sx q[2];
rz(-2.3872445) q[2];
sx q[2];
rz(0.95065439) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5444138) q[1];
sx q[1];
rz(-2.326818) q[1];
sx q[1];
rz(-2.9544403) q[1];
rz(-pi) q[2];
rz(1.6309687) q[3];
sx q[3];
rz(-0.98543692) q[3];
sx q[3];
rz(0.28096052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.88910237) q[2];
sx q[2];
rz(-0.95688755) q[2];
sx q[2];
rz(-1.0955742) q[2];
rz(-1.4787632) q[3];
sx q[3];
rz(-2.3507599) q[3];
sx q[3];
rz(1.3862632) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0304994) q[0];
sx q[0];
rz(-1.4249304) q[0];
sx q[0];
rz(1.4053364) q[0];
rz(1.7449215) q[1];
sx q[1];
rz(-1.3537355) q[1];
sx q[1];
rz(-2.2415846) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8843061) q[0];
sx q[0];
rz(-2.2955756) q[0];
sx q[0];
rz(-1.0777362) q[0];
rz(-pi) q[1];
rz(-0.58061231) q[2];
sx q[2];
rz(-2.0502649) q[2];
sx q[2];
rz(-1.1632077) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1088037) q[1];
sx q[1];
rz(-0.76251635) q[1];
sx q[1];
rz(-0.25188781) q[1];
x q[2];
rz(2.4056959) q[3];
sx q[3];
rz(-0.80965878) q[3];
sx q[3];
rz(-2.5084605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.46093837) q[2];
sx q[2];
rz(-2.9260981) q[2];
sx q[2];
rz(-0.94949618) q[2];
rz(2.5203868) q[3];
sx q[3];
rz(-2.216279) q[3];
sx q[3];
rz(-1.9882103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42156521) q[0];
sx q[0];
rz(-1.8864487) q[0];
sx q[0];
rz(-3.0928639) q[0];
rz(-0.85982927) q[1];
sx q[1];
rz(-2.8099334) q[1];
sx q[1];
rz(-2.8299832) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4954715) q[0];
sx q[0];
rz(-2.5866383) q[0];
sx q[0];
rz(1.036219) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0199058) q[2];
sx q[2];
rz(-0.60014987) q[2];
sx q[2];
rz(-0.7687591) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8383465) q[1];
sx q[1];
rz(-0.7184808) q[1];
sx q[1];
rz(1.9796728) q[1];
x q[2];
rz(2.1457003) q[3];
sx q[3];
rz(-2.6290335) q[3];
sx q[3];
rz(0.14581524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.17778808) q[2];
sx q[2];
rz(-2.9298941) q[2];
sx q[2];
rz(-1.6507899) q[2];
rz(-2.7866411) q[3];
sx q[3];
rz(-1.7345411) q[3];
sx q[3];
rz(1.8420334) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0896924) q[0];
sx q[0];
rz(-0.35744748) q[0];
sx q[0];
rz(-1.2667013) q[0];
rz(-3.025324) q[1];
sx q[1];
rz(-0.99028844) q[1];
sx q[1];
rz(-1.6273392) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8775692) q[0];
sx q[0];
rz(-0.91673393) q[0];
sx q[0];
rz(2.7436849) q[0];
x q[1];
rz(-3.051018) q[2];
sx q[2];
rz(-1.7901058) q[2];
sx q[2];
rz(0.62790576) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5473289) q[1];
sx q[1];
rz(-2.1143066) q[1];
sx q[1];
rz(3.0819403) q[1];
rz(-pi) q[2];
rz(-1.6125525) q[3];
sx q[3];
rz(-2.4126518) q[3];
sx q[3];
rz(1.4545914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9194455) q[2];
sx q[2];
rz(-0.47480348) q[2];
sx q[2];
rz(-1.1754645) q[2];
rz(2.2373824) q[3];
sx q[3];
rz(-1.2491106) q[3];
sx q[3];
rz(2.1632532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0745875) q[0];
sx q[0];
rz(-0.33127221) q[0];
sx q[0];
rz(2.7401155) q[0];
rz(1.1270771) q[1];
sx q[1];
rz(-1.8311484) q[1];
sx q[1];
rz(-2.3289767) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81174034) q[0];
sx q[0];
rz(-0.63416687) q[0];
sx q[0];
rz(1.5777753) q[0];
rz(-2.866686) q[2];
sx q[2];
rz(-1.3454352) q[2];
sx q[2];
rz(-1.4336841) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.26906313) q[1];
sx q[1];
rz(-1.7751964) q[1];
sx q[1];
rz(-0.60592954) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3964368) q[3];
sx q[3];
rz(-2.0260915) q[3];
sx q[3];
rz(2.8139092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.20069417) q[2];
sx q[2];
rz(-2.5914067) q[2];
sx q[2];
rz(-1.5852488) q[2];
rz(-0.19041348) q[3];
sx q[3];
rz(-0.75473458) q[3];
sx q[3];
rz(-1.93369) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82454005) q[0];
sx q[0];
rz(-0.38924488) q[0];
sx q[0];
rz(-0.12271605) q[0];
rz(-1.9484733) q[1];
sx q[1];
rz(-0.90843186) q[1];
sx q[1];
rz(0.76748031) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90119367) q[0];
sx q[0];
rz(-1.5195822) q[0];
sx q[0];
rz(0.9238433) q[0];
x q[1];
rz(-1.6558455) q[2];
sx q[2];
rz(-1.1901996) q[2];
sx q[2];
rz(-0.88262651) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7414416) q[1];
sx q[1];
rz(-1.5269244) q[1];
sx q[1];
rz(-2.2248101) q[1];
x q[2];
rz(-0.40624491) q[3];
sx q[3];
rz(-1.5749388) q[3];
sx q[3];
rz(2.3914482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3976589) q[2];
sx q[2];
rz(-0.021641061) q[2];
sx q[2];
rz(-1.5519315) q[2];
rz(1.4895561) q[3];
sx q[3];
rz(-1.4068312) q[3];
sx q[3];
rz(2.0261197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3234696) q[0];
sx q[0];
rz(-2.7796845) q[0];
sx q[0];
rz(-0.37539151) q[0];
rz(0.88343945) q[1];
sx q[1];
rz(-1.6049623) q[1];
sx q[1];
rz(0.72867957) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8189296) q[0];
sx q[0];
rz(-2.0897341) q[0];
sx q[0];
rz(0.82459088) q[0];
rz(-1.402114) q[2];
sx q[2];
rz(-1.8034435) q[2];
sx q[2];
rz(-2.6017435) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0496484) q[1];
sx q[1];
rz(-0.94251213) q[1];
sx q[1];
rz(-0.44558427) q[1];
rz(-pi) q[2];
rz(-0.94447926) q[3];
sx q[3];
rz(-0.91942838) q[3];
sx q[3];
rz(-0.29069549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4484619) q[2];
sx q[2];
rz(-1.5784266) q[2];
sx q[2];
rz(0.7473839) q[2];
rz(-1.4061617) q[3];
sx q[3];
rz(-1.9236671) q[3];
sx q[3];
rz(1.5860484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2324227) q[0];
sx q[0];
rz(-0.94038832) q[0];
sx q[0];
rz(-2.4160093) q[0];
rz(1.8046509) q[1];
sx q[1];
rz(-1.2533816) q[1];
sx q[1];
rz(0.57428378) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0030356) q[0];
sx q[0];
rz(-1.4113562) q[0];
sx q[0];
rz(-0.11388679) q[0];
rz(-pi) q[1];
rz(0.97522511) q[2];
sx q[2];
rz(-2.4831366) q[2];
sx q[2];
rz(1.083388) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.060238801) q[1];
sx q[1];
rz(-1.684184) q[1];
sx q[1];
rz(-0.84073034) q[1];
rz(-0.58135191) q[3];
sx q[3];
rz(-0.87331088) q[3];
sx q[3];
rz(-0.70840981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1672704) q[2];
sx q[2];
rz(-1.3784626) q[2];
sx q[2];
rz(0.66217011) q[2];
rz(-2.3769489) q[3];
sx q[3];
rz(-1.60631) q[3];
sx q[3];
rz(-1.8173328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26226703) q[0];
sx q[0];
rz(-0.58795324) q[0];
sx q[0];
rz(1.7437438) q[0];
rz(-2.0715879) q[1];
sx q[1];
rz(-1.1281697) q[1];
sx q[1];
rz(1.3319344) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4313244) q[0];
sx q[0];
rz(-1.6013751) q[0];
sx q[0];
rz(1.4200743) q[0];
rz(1.0584303) q[2];
sx q[2];
rz(-2.67423) q[2];
sx q[2];
rz(-0.99603727) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.20132682) q[1];
sx q[1];
rz(-1.3405217) q[1];
sx q[1];
rz(-2.9058232) q[1];
rz(-2.5601848) q[3];
sx q[3];
rz(-2.1801342) q[3];
sx q[3];
rz(-1.7096303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8269044) q[2];
sx q[2];
rz(-1.1004227) q[2];
sx q[2];
rz(-2.9617214) q[2];
rz(1.3966857) q[3];
sx q[3];
rz(-1.7161918) q[3];
sx q[3];
rz(-0.21656187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3708645) q[0];
sx q[0];
rz(-2.6612119) q[0];
sx q[0];
rz(-2.2330855) q[0];
rz(-2.6188359) q[1];
sx q[1];
rz(-2.0261384) q[1];
sx q[1];
rz(-1.1631858) q[1];
rz(1.1473473) q[2];
sx q[2];
rz(-0.97931391) q[2];
sx q[2];
rz(-0.58936832) q[2];
rz(-1.9369851) q[3];
sx q[3];
rz(-1.1259176) q[3];
sx q[3];
rz(2.068145) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
