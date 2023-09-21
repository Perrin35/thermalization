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
rz(2.9449129) q[0];
sx q[0];
rz(11.37698) q[0];
rz(-2.9130798) q[1];
sx q[1];
rz(-2.300188) q[1];
sx q[1];
rz(2.7639311) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084378622) q[0];
sx q[0];
rz(-2.9998261) q[0];
sx q[0];
rz(2.111582) q[0];
rz(-3.112552) q[2];
sx q[2];
rz(-1.8543058) q[2];
sx q[2];
rz(-0.107142) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.79757534) q[1];
sx q[1];
rz(-2.0098149) q[1];
sx q[1];
rz(-0.77787351) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.30794413) q[3];
sx q[3];
rz(-2.7044798) q[3];
sx q[3];
rz(-1.2446158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1315786) q[2];
sx q[2];
rz(-2.6476314) q[2];
sx q[2];
rz(-0.67260355) q[2];
rz(-0.16942313) q[3];
sx q[3];
rz(-0.38893458) q[3];
sx q[3];
rz(-1.3385564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23068962) q[0];
sx q[0];
rz(-0.90086532) q[0];
sx q[0];
rz(-0.22856523) q[0];
rz(-2.9810492) q[1];
sx q[1];
rz(-1.7030145) q[1];
sx q[1];
rz(-2.8536318) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44293091) q[0];
sx q[0];
rz(-1.3804111) q[0];
sx q[0];
rz(-0.18106826) q[0];
rz(2.1177883) q[2];
sx q[2];
rz(-2.4734801) q[2];
sx q[2];
rz(-0.39272768) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6504753) q[1];
sx q[1];
rz(-0.75956356) q[1];
sx q[1];
rz(0.61504765) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5324627) q[3];
sx q[3];
rz(-1.9209071) q[3];
sx q[3];
rz(0.73327524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.59445375) q[2];
sx q[2];
rz(-1.889166) q[2];
sx q[2];
rz(-1.4734369) q[2];
rz(-2.2041221) q[3];
sx q[3];
rz(-2.6963186) q[3];
sx q[3];
rz(0.37500769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32524747) q[0];
sx q[0];
rz(-2.6419817) q[0];
sx q[0];
rz(-2.8741799) q[0];
rz(-1.7193517) q[1];
sx q[1];
rz(-1.1317252) q[1];
sx q[1];
rz(0.95169383) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6355977) q[0];
sx q[0];
rz(-1.4800758) q[0];
sx q[0];
rz(-1.7631625) q[0];
x q[1];
rz(-2.9159413) q[2];
sx q[2];
rz(-1.0661085) q[2];
sx q[2];
rz(2.3329608) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3141331) q[1];
sx q[1];
rz(-0.810312) q[1];
sx q[1];
rz(2.00287) q[1];
rz(-pi) q[2];
rz(-2.1943201) q[3];
sx q[3];
rz(-2.0200649) q[3];
sx q[3];
rz(-2.505213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0114228) q[2];
sx q[2];
rz(-1.6458076) q[2];
sx q[2];
rz(0.34417957) q[2];
rz(0.88642818) q[3];
sx q[3];
rz(-2.8635946) q[3];
sx q[3];
rz(2.5755431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.13720559) q[0];
sx q[0];
rz(-1.6442278) q[0];
sx q[0];
rz(1.244506) q[0];
rz(3.0124774) q[1];
sx q[1];
rz(-1.3543509) q[1];
sx q[1];
rz(0.37277645) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5107721) q[0];
sx q[0];
rz(-1.6850867) q[0];
sx q[0];
rz(2.3833016) q[0];
rz(2.7083621) q[2];
sx q[2];
rz(-0.82651143) q[2];
sx q[2];
rz(0.42529688) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7074234) q[1];
sx q[1];
rz(-1.8137099) q[1];
sx q[1];
rz(-2.0704107) q[1];
rz(-pi) q[2];
rz(2.0471441) q[3];
sx q[3];
rz(-1.2762478) q[3];
sx q[3];
rz(2.0457207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5059775) q[2];
sx q[2];
rz(-1.6515235) q[2];
sx q[2];
rz(2.2367031) q[2];
rz(0.44057009) q[3];
sx q[3];
rz(-0.4959271) q[3];
sx q[3];
rz(-0.99159616) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47914094) q[0];
sx q[0];
rz(-0.15563706) q[0];
sx q[0];
rz(1.8537846) q[0];
rz(-1.6632535) q[1];
sx q[1];
rz(-1.9961424) q[1];
sx q[1];
rz(-2.6073661) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4596817) q[0];
sx q[0];
rz(-2.7026183) q[0];
sx q[0];
rz(-2.3460593) q[0];
x q[1];
rz(-1.0787557) q[2];
sx q[2];
rz(-1.2675708) q[2];
sx q[2];
rz(-2.6360896) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.6492669) q[1];
sx q[1];
rz(-1.5070792) q[1];
sx q[1];
rz(1.6462412) q[1];
rz(2.4228135) q[3];
sx q[3];
rz(-1.4712417) q[3];
sx q[3];
rz(-2.0841141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6405032) q[2];
sx q[2];
rz(-1.2155632) q[2];
sx q[2];
rz(0.14349288) q[2];
rz(1.7701373) q[3];
sx q[3];
rz(-1.2797132) q[3];
sx q[3];
rz(-2.9218856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7025529) q[0];
sx q[0];
rz(-0.87600231) q[0];
sx q[0];
rz(-0.23705661) q[0];
rz(-1.9006231) q[1];
sx q[1];
rz(-0.82890141) q[1];
sx q[1];
rz(1.3751078) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2415337) q[0];
sx q[0];
rz(-0.43474475) q[0];
sx q[0];
rz(0.66381201) q[0];
rz(0.74890045) q[2];
sx q[2];
rz(-1.5178711) q[2];
sx q[2];
rz(2.5089094) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.82735862) q[1];
sx q[1];
rz(-1.5463366) q[1];
sx q[1];
rz(0.41008653) q[1];
rz(-pi) q[2];
rz(0.70011219) q[3];
sx q[3];
rz(-2.5317319) q[3];
sx q[3];
rz(1.8240579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6087626) q[2];
sx q[2];
rz(-0.71981788) q[2];
sx q[2];
rz(0.14870816) q[2];
rz(3.1249629) q[3];
sx q[3];
rz(-0.36706585) q[3];
sx q[3];
rz(-0.19255157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6454813) q[0];
sx q[0];
rz(-2.8493024) q[0];
sx q[0];
rz(0.87316978) q[0];
rz(-0.9219777) q[1];
sx q[1];
rz(-2.0261804) q[1];
sx q[1];
rz(-0.08392863) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31939313) q[0];
sx q[0];
rz(-2.5376352) q[0];
sx q[0];
rz(-1.970406) q[0];
rz(-1.3533808) q[2];
sx q[2];
rz(-1.916421) q[2];
sx q[2];
rz(-0.84920041) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2578656) q[1];
sx q[1];
rz(-0.95974937) q[1];
sx q[1];
rz(2.3962767) q[1];
rz(-pi) q[2];
rz(-0.069025741) q[3];
sx q[3];
rz(-0.8419753) q[3];
sx q[3];
rz(0.1544827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.38368791) q[2];
sx q[2];
rz(-0.44054511) q[2];
sx q[2];
rz(0.54523462) q[2];
rz(-0.43045726) q[3];
sx q[3];
rz(-2.0038219) q[3];
sx q[3];
rz(2.8366413) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51171821) q[0];
sx q[0];
rz(-0.015462333) q[0];
sx q[0];
rz(1.9301201) q[0];
rz(0.97310549) q[1];
sx q[1];
rz(-0.5935697) q[1];
sx q[1];
rz(2.1957695) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0154008) q[0];
sx q[0];
rz(-1.7915103) q[0];
sx q[0];
rz(-1.4908233) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.81994762) q[2];
sx q[2];
rz(-0.87265271) q[2];
sx q[2];
rz(1.9422216) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.91499) q[1];
sx q[1];
rz(-1.1922925) q[1];
sx q[1];
rz(1.6016017) q[1];
rz(0.66611992) q[3];
sx q[3];
rz(-1.5770116) q[3];
sx q[3];
rz(-2.4339649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8699845) q[2];
sx q[2];
rz(-1.7249858) q[2];
sx q[2];
rz(0.28042173) q[2];
rz(2.5366606) q[3];
sx q[3];
rz(-1.0792024) q[3];
sx q[3];
rz(-0.98208565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1542926) q[0];
sx q[0];
rz(-2.2315318) q[0];
sx q[0];
rz(-0.61532414) q[0];
rz(2.212021) q[1];
sx q[1];
rz(-0.85837448) q[1];
sx q[1];
rz(0.5756793) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0711813) q[0];
sx q[0];
rz(-1.0072664) q[0];
sx q[0];
rz(2.1348743) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0909537) q[2];
sx q[2];
rz(-2.9165604) q[2];
sx q[2];
rz(2.2573543) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.5902192) q[1];
sx q[1];
rz(-3.0312523) q[1];
sx q[1];
rz(-1.8054086) q[1];
rz(1.2647763) q[3];
sx q[3];
rz(-0.56193202) q[3];
sx q[3];
rz(-2.7276873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3118887) q[2];
sx q[2];
rz(-2.3788033) q[2];
sx q[2];
rz(-0.021961948) q[2];
rz(-2.9711376) q[3];
sx q[3];
rz(-2.1178092) q[3];
sx q[3];
rz(0.26486614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4158674) q[0];
sx q[0];
rz(-2.2150345) q[0];
sx q[0];
rz(-2.5073994) q[0];
rz(-0.028907396) q[1];
sx q[1];
rz(-2.3529265) q[1];
sx q[1];
rz(-2.172487) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0646107) q[0];
sx q[0];
rz(-1.6425898) q[0];
sx q[0];
rz(0.063477091) q[0];
rz(0.74659851) q[2];
sx q[2];
rz(-1.4274297) q[2];
sx q[2];
rz(-1.5310841) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.58418729) q[1];
sx q[1];
rz(-0.69328847) q[1];
sx q[1];
rz(0.28179817) q[1];
x q[2];
rz(2.7006847) q[3];
sx q[3];
rz(-1.3135859) q[3];
sx q[3];
rz(-2.5901026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4518296) q[2];
sx q[2];
rz(-1.4844866) q[2];
sx q[2];
rz(-0.36007145) q[2];
rz(1.6137971) q[3];
sx q[3];
rz(-2.4751622) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9344899) q[0];
sx q[0];
rz(-1.5710545) q[0];
sx q[0];
rz(1.5221773) q[0];
rz(0.044152505) q[1];
sx q[1];
rz(-1.6828729) q[1];
sx q[1];
rz(2.0353459) q[1];
rz(-2.8051878) q[2];
sx q[2];
rz(-1.9336666) q[2];
sx q[2];
rz(-1.578707) q[2];
rz(1.2561856) q[3];
sx q[3];
rz(-1.1369858) q[3];
sx q[3];
rz(2.9997957) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
