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
rz(1.8013826) q[0];
sx q[0];
rz(-0.27685452) q[0];
sx q[0];
rz(-1.0387596) q[0];
rz(-0.34173319) q[1];
sx q[1];
rz(-2.3050397) q[1];
sx q[1];
rz(0.41681448) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87931235) q[0];
sx q[0];
rz(-1.4646155) q[0];
sx q[0];
rz(-0.72632974) q[0];
rz(-pi) q[1];
rz(-2.5854163) q[2];
sx q[2];
rz(-1.7814629) q[2];
sx q[2];
rz(-0.81708252) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9824901) q[1];
sx q[1];
rz(-2.1381408) q[1];
sx q[1];
rz(2.2654387) q[1];
rz(-pi) q[2];
rz(2.9800426) q[3];
sx q[3];
rz(-1.6027502) q[3];
sx q[3];
rz(-0.7940426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1699528) q[2];
sx q[2];
rz(-1.6728741) q[2];
sx q[2];
rz(-1.2539585) q[2];
rz(-1.5818671) q[3];
sx q[3];
rz(-2.4032148) q[3];
sx q[3];
rz(2.019465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1935683) q[0];
sx q[0];
rz(-1.3688315) q[0];
sx q[0];
rz(2.3728306) q[0];
rz(1.7747152) q[1];
sx q[1];
rz(-1.8194852) q[1];
sx q[1];
rz(2.1034525) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5032673) q[0];
sx q[0];
rz(-1.5621788) q[0];
sx q[0];
rz(1.5609972) q[0];
x q[1];
rz(2.4203796) q[2];
sx q[2];
rz(-1.291847) q[2];
sx q[2];
rz(2.9245289) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8703192) q[1];
sx q[1];
rz(-2.4044577) q[1];
sx q[1];
rz(2.5802274) q[1];
rz(-pi) q[2];
rz(3.1225864) q[3];
sx q[3];
rz(-1.0822625) q[3];
sx q[3];
rz(-1.8754539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.64608964) q[2];
sx q[2];
rz(-1.8547408) q[2];
sx q[2];
rz(0.82026473) q[2];
rz(-2.8881554) q[3];
sx q[3];
rz(-0.35496747) q[3];
sx q[3];
rz(2.7804815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8609817) q[0];
sx q[0];
rz(-2.6237223) q[0];
sx q[0];
rz(-0.88731998) q[0];
rz(-0.53030983) q[1];
sx q[1];
rz(-2.2128426) q[1];
sx q[1];
rz(-1.1393772) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4654008) q[0];
sx q[0];
rz(-1.4290819) q[0];
sx q[0];
rz(-0.34872524) q[0];
rz(-pi) q[1];
rz(2.1007295) q[2];
sx q[2];
rz(-1.8609253) q[2];
sx q[2];
rz(-1.1159971) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.85998409) q[1];
sx q[1];
rz(-1.3251332) q[1];
sx q[1];
rz(-1.5857693) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5185131) q[3];
sx q[3];
rz(-1.0526592) q[3];
sx q[3];
rz(-2.5123346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1246216) q[2];
sx q[2];
rz(-1.6861702) q[2];
sx q[2];
rz(-2.7023081) q[2];
rz(-0.47075054) q[3];
sx q[3];
rz(-0.079340383) q[3];
sx q[3];
rz(1.9973756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73862326) q[0];
sx q[0];
rz(-2.2082177) q[0];
sx q[0];
rz(0.11548197) q[0];
rz(3.0237517) q[1];
sx q[1];
rz(-1.9385447) q[1];
sx q[1];
rz(1.3062564) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8354315) q[0];
sx q[0];
rz(-1.8580313) q[0];
sx q[0];
rz(-1.1348508) q[0];
rz(-pi) q[1];
rz(1.1924065) q[2];
sx q[2];
rz(-2.1956964) q[2];
sx q[2];
rz(1.6922127) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.072767898) q[1];
sx q[1];
rz(-2.1842274) q[1];
sx q[1];
rz(0.83534209) q[1];
rz(-pi) q[2];
rz(0.0751817) q[3];
sx q[3];
rz(-0.99937056) q[3];
sx q[3];
rz(-0.21440766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8370886) q[2];
sx q[2];
rz(-1.8279165) q[2];
sx q[2];
rz(-0.15466776) q[2];
rz(-1.539544) q[3];
sx q[3];
rz(-1.8013026) q[3];
sx q[3];
rz(-2.535517) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57347572) q[0];
sx q[0];
rz(-1.0433759) q[0];
sx q[0];
rz(0.83531761) q[0];
rz(-2.4256445) q[1];
sx q[1];
rz(-0.42778152) q[1];
sx q[1];
rz(-3.0221525) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41573856) q[0];
sx q[0];
rz(-1.6266168) q[0];
sx q[0];
rz(3.124696) q[0];
rz(-pi) q[1];
rz(-2.6921656) q[2];
sx q[2];
rz(-1.4327421) q[2];
sx q[2];
rz(1.4269478) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.24158289) q[1];
sx q[1];
rz(-2.0323583) q[1];
sx q[1];
rz(0.66630967) q[1];
rz(-pi) q[2];
rz(1.9914845) q[3];
sx q[3];
rz(-2.0644232) q[3];
sx q[3];
rz(-2.8378262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2401838) q[2];
sx q[2];
rz(-0.26075026) q[2];
sx q[2];
rz(3.0805947) q[2];
rz(-3.0933464) q[3];
sx q[3];
rz(-1.2333074) q[3];
sx q[3];
rz(-2.4129996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7399087) q[0];
sx q[0];
rz(-2.4329199) q[0];
sx q[0];
rz(2.0106864) q[0];
rz(-0.4785969) q[1];
sx q[1];
rz(-2.5391948) q[1];
sx q[1];
rz(-1.4071677) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3597108) q[0];
sx q[0];
rz(-0.053584307) q[0];
sx q[0];
rz(-0.028938541) q[0];
x q[1];
rz(-0.28193177) q[2];
sx q[2];
rz(-0.86840668) q[2];
sx q[2];
rz(1.3739746) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1651014) q[1];
sx q[1];
rz(-2.8860984) q[1];
sx q[1];
rz(2.0352023) q[1];
x q[2];
rz(-2.9286372) q[3];
sx q[3];
rz(-1.3939438) q[3];
sx q[3];
rz(0.18421728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.56167928) q[2];
sx q[2];
rz(-0.40739569) q[2];
sx q[2];
rz(1.9631867) q[2];
rz(-1.0493578) q[3];
sx q[3];
rz(-0.5368084) q[3];
sx q[3];
rz(1.8572846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2105763) q[0];
sx q[0];
rz(-1.9518305) q[0];
sx q[0];
rz(1.3354906) q[0];
rz(-1.3212851) q[1];
sx q[1];
rz(-2.0704465) q[1];
sx q[1];
rz(-2.8335422) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0548693) q[0];
sx q[0];
rz(-0.44741524) q[0];
sx q[0];
rz(-2.6347876) q[0];
rz(-pi) q[1];
rz(2.6263138) q[2];
sx q[2];
rz(-1.9204428) q[2];
sx q[2];
rz(-0.51008979) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.25702204) q[1];
sx q[1];
rz(-2.0557457) q[1];
sx q[1];
rz(0.91241769) q[1];
rz(-0.04754504) q[3];
sx q[3];
rz(-1.1851839) q[3];
sx q[3];
rz(-1.6507208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.41219741) q[2];
sx q[2];
rz(-2.2008379) q[2];
sx q[2];
rz(-0.13128734) q[2];
rz(0.73733759) q[3];
sx q[3];
rz(-1.8359343) q[3];
sx q[3];
rz(2.9837515) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3103631) q[0];
sx q[0];
rz(-2.0518301) q[0];
sx q[0];
rz(-0.53037733) q[0];
rz(1.7165548) q[1];
sx q[1];
rz(-2.0030463) q[1];
sx q[1];
rz(-2.4200965) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24151267) q[0];
sx q[0];
rz(-0.78966138) q[0];
sx q[0];
rz(-2.4898743) q[0];
x q[1];
rz(-2.8558735) q[2];
sx q[2];
rz(-2.6556394) q[2];
sx q[2];
rz(-3.1348117) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4911982) q[1];
sx q[1];
rz(-1.7868687) q[1];
sx q[1];
rz(2.9465527) q[1];
rz(-pi) q[2];
rz(-3.0305392) q[3];
sx q[3];
rz(-0.78563443) q[3];
sx q[3];
rz(-1.3291886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.43712744) q[2];
sx q[2];
rz(-2.7361054) q[2];
sx q[2];
rz(1.6638157) q[2];
rz(-1.1635228) q[3];
sx q[3];
rz(-1.082837) q[3];
sx q[3];
rz(1.5014974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8116233) q[0];
sx q[0];
rz(-1.4412619) q[0];
sx q[0];
rz(-0.60923088) q[0];
rz(1.5902663) q[1];
sx q[1];
rz(-0.32569277) q[1];
sx q[1];
rz(-1.4564266) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0349755) q[0];
sx q[0];
rz(-1.9272695) q[0];
sx q[0];
rz(2.8750505) q[0];
x q[1];
rz(-1.6316192) q[2];
sx q[2];
rz(-2.0552539) q[2];
sx q[2];
rz(-1.6866419) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.506914) q[1];
sx q[1];
rz(-1.7326446) q[1];
sx q[1];
rz(-2.0353725) q[1];
rz(-pi) q[2];
x q[2];
rz(0.84734579) q[3];
sx q[3];
rz(-1.0718126) q[3];
sx q[3];
rz(2.8049198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7234601) q[2];
sx q[2];
rz(-2.2432566) q[2];
sx q[2];
rz(2.2612803) q[2];
rz(2.9355925) q[3];
sx q[3];
rz(-2.7166631) q[3];
sx q[3];
rz(-0.4900842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77077615) q[0];
sx q[0];
rz(-0.66148615) q[0];
sx q[0];
rz(-3.1296375) q[0];
rz(-1.5097584) q[1];
sx q[1];
rz(-2.6963574) q[1];
sx q[1];
rz(-2.5451122) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79967989) q[0];
sx q[0];
rz(-2.0401067) q[0];
sx q[0];
rz(1.5494359) q[0];
x q[1];
rz(0.61052236) q[2];
sx q[2];
rz(-0.43402616) q[2];
sx q[2];
rz(-2.4422925) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.66377992) q[1];
sx q[1];
rz(-0.48250178) q[1];
sx q[1];
rz(1.7412118) q[1];
rz(-pi) q[2];
x q[2];
rz(0.29774547) q[3];
sx q[3];
rz(-1.4199054) q[3];
sx q[3];
rz(-0.45805629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.97800469) q[2];
sx q[2];
rz(-0.96077335) q[2];
sx q[2];
rz(-0.19700024) q[2];
rz(1.152285) q[3];
sx q[3];
rz(-0.14324337) q[3];
sx q[3];
rz(2.6939189) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2780509) q[0];
sx q[0];
rz(-1.0095689) q[0];
sx q[0];
rz(2.9472245) q[0];
rz(2.9327783) q[1];
sx q[1];
rz(-1.5369692) q[1];
sx q[1];
rz(2.1280638) q[1];
rz(-0.88125689) q[2];
sx q[2];
rz(-1.4980346) q[2];
sx q[2];
rz(-3.0537506) q[2];
rz(2.267425) q[3];
sx q[3];
rz(-2.4258191) q[3];
sx q[3];
rz(1.3620472) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
