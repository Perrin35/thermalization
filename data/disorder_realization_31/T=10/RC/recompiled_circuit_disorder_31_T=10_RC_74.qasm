OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6150317) q[0];
sx q[0];
rz(2.5685413) q[0];
sx q[0];
rz(11.723784) q[0];
rz(2.1057582) q[1];
sx q[1];
rz(8.3254568) q[1];
sx q[1];
rz(7.96666) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6095088) q[0];
sx q[0];
rz(-0.77021399) q[0];
sx q[0];
rz(0.30755933) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5277532) q[2];
sx q[2];
rz(-1.5868574) q[2];
sx q[2];
rz(-1.2889372) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.1229413) q[1];
sx q[1];
rz(-0.39995799) q[1];
sx q[1];
rz(0.33756983) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59431608) q[3];
sx q[3];
rz(-1.672097) q[3];
sx q[3];
rz(0.017410226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.52790102) q[2];
sx q[2];
rz(-1.0062904) q[2];
sx q[2];
rz(-0.17949417) q[2];
rz(-1.9159296) q[3];
sx q[3];
rz(-1.7951199) q[3];
sx q[3];
rz(-2.3195482) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74801385) q[0];
sx q[0];
rz(-2.2606235) q[0];
sx q[0];
rz(0.32546145) q[0];
rz(-1.7851967) q[1];
sx q[1];
rz(-2.0929095) q[1];
sx q[1];
rz(-1.9869841) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3852859) q[0];
sx q[0];
rz(-3.1161838) q[0];
sx q[0];
rz(-0.73861648) q[0];
x q[1];
rz(-2.1968368) q[2];
sx q[2];
rz(-1.895004) q[2];
sx q[2];
rz(-2.7446483) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.11101152) q[1];
sx q[1];
rz(-1.6554553) q[1];
sx q[1];
rz(-0.76534033) q[1];
rz(-pi) q[2];
rz(-0.18493821) q[3];
sx q[3];
rz(-1.3122845) q[3];
sx q[3];
rz(2.1188494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6894199) q[2];
sx q[2];
rz(-1.8916811) q[2];
sx q[2];
rz(-2.2581805) q[2];
rz(0.47131053) q[3];
sx q[3];
rz(-1.4383957) q[3];
sx q[3];
rz(-0.78770351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8283591) q[0];
sx q[0];
rz(-1.4947083) q[0];
sx q[0];
rz(-1.5154243) q[0];
rz(2.5405163) q[1];
sx q[1];
rz(-2.5939012) q[1];
sx q[1];
rz(-1.0916969) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1356782) q[0];
sx q[0];
rz(-2.1164829) q[0];
sx q[0];
rz(2.2360327) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41468427) q[2];
sx q[2];
rz(-0.95699246) q[2];
sx q[2];
rz(2.2874122) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.25635438) q[1];
sx q[1];
rz(-1.2043722) q[1];
sx q[1];
rz(0.79856915) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12797171) q[3];
sx q[3];
rz(-1.5068753) q[3];
sx q[3];
rz(1.2325866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8213356) q[2];
sx q[2];
rz(-0.50575033) q[2];
sx q[2];
rz(0.88095218) q[2];
rz(1.7679924) q[3];
sx q[3];
rz(-1.526984) q[3];
sx q[3];
rz(-1.0176456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3110733) q[0];
sx q[0];
rz(-1.3922465) q[0];
sx q[0];
rz(-0.4367035) q[0];
rz(2.9084335) q[1];
sx q[1];
rz(-1.2522839) q[1];
sx q[1];
rz(-2.8312347) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6757641) q[0];
sx q[0];
rz(-0.40256631) q[0];
sx q[0];
rz(-2.7990544) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6967625) q[2];
sx q[2];
rz(-2.2519886) q[2];
sx q[2];
rz(1.7484401) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0536641) q[1];
sx q[1];
rz(-1.2128608) q[1];
sx q[1];
rz(3.1217561) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.035590812) q[3];
sx q[3];
rz(-1.6944052) q[3];
sx q[3];
rz(-1.8087169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.13005304) q[2];
sx q[2];
rz(-0.71802846) q[2];
sx q[2];
rz(-2.0641573) q[2];
rz(-0.056190101) q[3];
sx q[3];
rz(-0.63779938) q[3];
sx q[3];
rz(1.5475387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.9524277) q[0];
sx q[0];
rz(-1.0452894) q[0];
sx q[0];
rz(2.8919343) q[0];
rz(-1.5646308) q[1];
sx q[1];
rz(-2.3639634) q[1];
sx q[1];
rz(-0.87019428) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10660431) q[0];
sx q[0];
rz(-0.37811324) q[0];
sx q[0];
rz(2.5614221) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.22360794) q[2];
sx q[2];
rz(-0.70906559) q[2];
sx q[2];
rz(2.2774334) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.732547) q[1];
sx q[1];
rz(-0.88102075) q[1];
sx q[1];
rz(-1.9802666) q[1];
rz(-pi) q[2];
rz(-1.0088483) q[3];
sx q[3];
rz(-1.7147439) q[3];
sx q[3];
rz(0.89509237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.27328086) q[2];
sx q[2];
rz(-1.8153278) q[2];
sx q[2];
rz(-2.4678521) q[2];
rz(-2.8379748) q[3];
sx q[3];
rz(-1.9165336) q[3];
sx q[3];
rz(1.822086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34981397) q[0];
sx q[0];
rz(-0.93739167) q[0];
sx q[0];
rz(0.2579903) q[0];
rz(0.42516431) q[1];
sx q[1];
rz(-0.95562569) q[1];
sx q[1];
rz(-1.649883) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0598037) q[0];
sx q[0];
rz(-2.7002618) q[0];
sx q[0];
rz(-2.9102737) q[0];
rz(-pi) q[1];
rz(2.4620373) q[2];
sx q[2];
rz(-1.9050042) q[2];
sx q[2];
rz(-2.213775) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7266453) q[1];
sx q[1];
rz(-0.84268314) q[1];
sx q[1];
rz(2.3689518) q[1];
rz(3.0400279) q[3];
sx q[3];
rz(-1.2079117) q[3];
sx q[3];
rz(0.18728072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1288746) q[2];
sx q[2];
rz(-2.1779163) q[2];
sx q[2];
rz(-0.091726124) q[2];
rz(-2.2979459) q[3];
sx q[3];
rz(-2.1618312) q[3];
sx q[3];
rz(-0.89404026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3180852) q[0];
sx q[0];
rz(-1.2473236) q[0];
sx q[0];
rz(0.41123018) q[0];
rz(0.86589083) q[1];
sx q[1];
rz(-2.829268) q[1];
sx q[1];
rz(-3.1076028) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54995178) q[0];
sx q[0];
rz(-2.3410428) q[0];
sx q[0];
rz(0.22649015) q[0];
rz(-pi) q[1];
rz(0.93416874) q[2];
sx q[2];
rz(-1.120943) q[2];
sx q[2];
rz(-2.5069782) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9719203) q[1];
sx q[1];
rz(-1.5202513) q[1];
sx q[1];
rz(2.2488942) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5266685) q[3];
sx q[3];
rz(-1.4498386) q[3];
sx q[3];
rz(1.3161591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5380481) q[2];
sx q[2];
rz(-0.59843439) q[2];
sx q[2];
rz(-0.87654385) q[2];
rz(-2.792568) q[3];
sx q[3];
rz(-1.9411496) q[3];
sx q[3];
rz(2.9984737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2440764) q[0];
sx q[0];
rz(-1.709047) q[0];
sx q[0];
rz(-0.39392719) q[0];
rz(-0.36755964) q[1];
sx q[1];
rz(-1.7575248) q[1];
sx q[1];
rz(1.4454909) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60364265) q[0];
sx q[0];
rz(-1.1374439) q[0];
sx q[0];
rz(-0.005714697) q[0];
rz(-pi) q[1];
rz(2.1922853) q[2];
sx q[2];
rz(-2.0645803) q[2];
sx q[2];
rz(-0.91044237) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1576924) q[1];
sx q[1];
rz(-2.1196113) q[1];
sx q[1];
rz(-1.3859315) q[1];
rz(-pi) q[2];
x q[2];
rz(0.74650713) q[3];
sx q[3];
rz(-0.81614796) q[3];
sx q[3];
rz(-1.5619123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8470856) q[2];
sx q[2];
rz(-2.2512348) q[2];
sx q[2];
rz(-0.40714804) q[2];
rz(-1.6242705) q[3];
sx q[3];
rz(-1.9842691) q[3];
sx q[3];
rz(2.8919162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80609926) q[0];
sx q[0];
rz(-2.6265916) q[0];
sx q[0];
rz(-1.8898213) q[0];
rz(0.66954008) q[1];
sx q[1];
rz(-1.1839097) q[1];
sx q[1];
rz(0.30977419) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4315902) q[0];
sx q[0];
rz(-2.6012523) q[0];
sx q[0];
rz(-0.24775981) q[0];
x q[1];
rz(1.7477112) q[2];
sx q[2];
rz(-1.5614911) q[2];
sx q[2];
rz(0.44405802) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.27948353) q[1];
sx q[1];
rz(-2.9276507) q[1];
sx q[1];
rz(1.8474384) q[1];
rz(-pi) q[2];
rz(-2.5449949) q[3];
sx q[3];
rz(-1.784424) q[3];
sx q[3];
rz(1.9006157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4391675) q[2];
sx q[2];
rz(-0.71321407) q[2];
sx q[2];
rz(-1.2072198) q[2];
rz(-1.0369982) q[3];
sx q[3];
rz(-1.2456649) q[3];
sx q[3];
rz(2.4859378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.050215125) q[0];
sx q[0];
rz(-1.3239048) q[0];
sx q[0];
rz(-1.2058831) q[0];
rz(-2.5559015) q[1];
sx q[1];
rz(-1.0810477) q[1];
sx q[1];
rz(-1.4996128) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4831055) q[0];
sx q[0];
rz(-2.8327836) q[0];
sx q[0];
rz(-1.9103861) q[0];
rz(-pi) q[1];
rz(3.0707804) q[2];
sx q[2];
rz(-2.376308) q[2];
sx q[2];
rz(0.27809696) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.26350281) q[1];
sx q[1];
rz(-2.6260758) q[1];
sx q[1];
rz(2.0104736) q[1];
x q[2];
rz(-1.5851192) q[3];
sx q[3];
rz(-2.9687772) q[3];
sx q[3];
rz(0.9848136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2575834) q[2];
sx q[2];
rz(-1.7929701) q[2];
sx q[2];
rz(2.1949027) q[2];
rz(-0.36869129) q[3];
sx q[3];
rz(-1.576141) q[3];
sx q[3];
rz(-2.6856016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35836999) q[0];
sx q[0];
rz(-1.9932278) q[0];
sx q[0];
rz(2.7182462) q[0];
rz(0.070925698) q[1];
sx q[1];
rz(-1.6880886) q[1];
sx q[1];
rz(-0.26500519) q[1];
rz(-2.0680188) q[2];
sx q[2];
rz(-1.1560658) q[2];
sx q[2];
rz(1.8552468) q[2];
rz(-0.55340135) q[3];
sx q[3];
rz(-0.80080606) q[3];
sx q[3];
rz(-1.3047119) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
