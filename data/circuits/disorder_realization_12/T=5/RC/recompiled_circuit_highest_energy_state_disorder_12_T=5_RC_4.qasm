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
rz(-0.15859088) q[0];
sx q[0];
rz(2.6245485) q[0];
sx q[0];
rz(11.256097) q[0];
rz(-0.34215555) q[1];
sx q[1];
rz(-1.181239) q[1];
sx q[1];
rz(0.96460834) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2268596) q[0];
sx q[0];
rz(-0.64382416) q[0];
sx q[0];
rz(-0.22636803) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7959441) q[2];
sx q[2];
rz(-1.7549577) q[2];
sx q[2];
rz(0.08376567) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8339616) q[1];
sx q[1];
rz(-1.2738859) q[1];
sx q[1];
rz(2.3351257) q[1];
rz(-pi) q[2];
rz(-2.2418666) q[3];
sx q[3];
rz(-1.6613071) q[3];
sx q[3];
rz(-2.8495726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.0059119314) q[2];
sx q[2];
rz(-0.19401208) q[2];
sx q[2];
rz(1.4646336) q[2];
rz(1.3095193) q[3];
sx q[3];
rz(-2.194761) q[3];
sx q[3];
rz(2.8215698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86318535) q[0];
sx q[0];
rz(-1.139737) q[0];
sx q[0];
rz(1.3673258) q[0];
rz(1.4890081) q[1];
sx q[1];
rz(-2.3431578) q[1];
sx q[1];
rz(2.8489825) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3037138) q[0];
sx q[0];
rz(-2.7927164) q[0];
sx q[0];
rz(-0.29033355) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4762223) q[2];
sx q[2];
rz(-2.9369246) q[2];
sx q[2];
rz(0.39835793) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3019588) q[1];
sx q[1];
rz(-1.3657161) q[1];
sx q[1];
rz(-0.2409711) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8752453) q[3];
sx q[3];
rz(-1.5070591) q[3];
sx q[3];
rz(1.6132465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1613529) q[2];
sx q[2];
rz(-0.62675256) q[2];
sx q[2];
rz(-2.9577067) q[2];
rz(-0.45267496) q[3];
sx q[3];
rz(-1.5225182) q[3];
sx q[3];
rz(-0.12990738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1995131) q[0];
sx q[0];
rz(-0.71490723) q[0];
sx q[0];
rz(2.3364501) q[0];
rz(-1.0792271) q[1];
sx q[1];
rz(-0.77003038) q[1];
sx q[1];
rz(1.8260746) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41322177) q[0];
sx q[0];
rz(-0.72997366) q[0];
sx q[0];
rz(-3.0527924) q[0];
rz(-pi) q[1];
rz(2.0657077) q[2];
sx q[2];
rz(-0.80782164) q[2];
sx q[2];
rz(3.0095095) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6650369) q[1];
sx q[1];
rz(-1.8522693) q[1];
sx q[1];
rz(0.095515619) q[1];
x q[2];
rz(-0.27045111) q[3];
sx q[3];
rz(-0.31336774) q[3];
sx q[3];
rz(-2.6033786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2074073) q[2];
sx q[2];
rz(-1.4855874) q[2];
sx q[2];
rz(-2.5918813) q[2];
rz(3.1304729) q[3];
sx q[3];
rz(-2.34237) q[3];
sx q[3];
rz(1.8196677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84871197) q[0];
sx q[0];
rz(-2.1826545) q[0];
sx q[0];
rz(0.30353656) q[0];
rz(2.983298) q[1];
sx q[1];
rz(-1.0946495) q[1];
sx q[1];
rz(1.0096445) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5490259) q[0];
sx q[0];
rz(-1.1404622) q[0];
sx q[0];
rz(2.8293912) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.539297) q[2];
sx q[2];
rz(-1.3118852) q[2];
sx q[2];
rz(2.7342755) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9480966) q[1];
sx q[1];
rz(-2.0442914) q[1];
sx q[1];
rz(1.3934474) q[1];
rz(-pi) q[2];
x q[2];
rz(0.49780881) q[3];
sx q[3];
rz(-1.5890997) q[3];
sx q[3];
rz(2.8895484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.9369072) q[2];
sx q[2];
rz(-0.82131177) q[2];
sx q[2];
rz(2.8724907) q[2];
rz(1.6875632) q[3];
sx q[3];
rz(-1.5917835) q[3];
sx q[3];
rz(-1.5627741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8708385) q[0];
sx q[0];
rz(-2.0814867) q[0];
sx q[0];
rz(0.36218542) q[0];
rz(-2.4902792) q[1];
sx q[1];
rz(-1.4910699) q[1];
sx q[1];
rz(1.6168894) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0391738) q[0];
sx q[0];
rz(-2.1764767) q[0];
sx q[0];
rz(1.4132947) q[0];
x q[1];
rz(-0.94793041) q[2];
sx q[2];
rz(-1.6375223) q[2];
sx q[2];
rz(0.78363505) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.406639) q[1];
sx q[1];
rz(-1.2124774) q[1];
sx q[1];
rz(2.1291158) q[1];
rz(-pi) q[2];
rz(1.175143) q[3];
sx q[3];
rz(-0.86651245) q[3];
sx q[3];
rz(-2.8362398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0321956) q[2];
sx q[2];
rz(-2.1869662) q[2];
sx q[2];
rz(-1.7677914) q[2];
rz(0.43921709) q[3];
sx q[3];
rz(-2.4924811) q[3];
sx q[3];
rz(-0.60905987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0405025) q[0];
sx q[0];
rz(-0.73223615) q[0];
sx q[0];
rz(2.8771583) q[0];
rz(-0.6791555) q[1];
sx q[1];
rz(-0.87398386) q[1];
sx q[1];
rz(2.1566379) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.934498) q[0];
sx q[0];
rz(-0.45299981) q[0];
sx q[0];
rz(-2.3191602) q[0];
x q[1];
rz(0.28552766) q[2];
sx q[2];
rz(-2.0047024) q[2];
sx q[2];
rz(-1.3615695) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.86414908) q[1];
sx q[1];
rz(-1.207749) q[1];
sx q[1];
rz(2.6226642) q[1];
x q[2];
rz(-1.156928) q[3];
sx q[3];
rz(-2.664251) q[3];
sx q[3];
rz(-1.9628003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1017477) q[2];
sx q[2];
rz(-2.9797649) q[2];
sx q[2];
rz(2.7092773) q[2];
rz(1.013422) q[3];
sx q[3];
rz(-1.6067959) q[3];
sx q[3];
rz(2.5920946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6559615) q[0];
sx q[0];
rz(-0.39325842) q[0];
sx q[0];
rz(-1.1095169) q[0];
rz(-0.79311496) q[1];
sx q[1];
rz(-1.7879281) q[1];
sx q[1];
rz(1.5464334) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8414559) q[0];
sx q[0];
rz(-1.0086035) q[0];
sx q[0];
rz(-2.4676222) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4199615) q[2];
sx q[2];
rz(-1.411937) q[2];
sx q[2];
rz(1.5312486) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.90518752) q[1];
sx q[1];
rz(-1.3705472) q[1];
sx q[1];
rz(-0.88242857) q[1];
rz(3.0199354) q[3];
sx q[3];
rz(-2.4215048) q[3];
sx q[3];
rz(2.4984604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3194797) q[2];
sx q[2];
rz(-0.18222465) q[2];
sx q[2];
rz(-1.3694084) q[2];
rz(-0.93044126) q[3];
sx q[3];
rz(-1.7934099) q[3];
sx q[3];
rz(-1.5038917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1996138) q[0];
sx q[0];
rz(-1.1811341) q[0];
sx q[0];
rz(-2.8379295) q[0];
rz(0.97967255) q[1];
sx q[1];
rz(-2.517608) q[1];
sx q[1];
rz(0.40506515) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9054149) q[0];
sx q[0];
rz(-1.1943895) q[0];
sx q[0];
rz(1.3628528) q[0];
rz(-pi) q[1];
rz(-1.9900277) q[2];
sx q[2];
rz(-1.7287325) q[2];
sx q[2];
rz(-1.7309703) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.36282) q[1];
sx q[1];
rz(-1.3937794) q[1];
sx q[1];
rz(-2.5702012) q[1];
x q[2];
rz(-1.1478506) q[3];
sx q[3];
rz(-0.4627403) q[3];
sx q[3];
rz(3.1266318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.27641174) q[2];
sx q[2];
rz(-2.1750735) q[2];
sx q[2];
rz(2.442404) q[2];
rz(-0.80896038) q[3];
sx q[3];
rz(-1.304108) q[3];
sx q[3];
rz(0.32304025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9998099) q[0];
sx q[0];
rz(-1.9843822) q[0];
sx q[0];
rz(3.004177) q[0];
rz(0.049086463) q[1];
sx q[1];
rz(-2.0259435) q[1];
sx q[1];
rz(-0.5079937) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2823454) q[0];
sx q[0];
rz(-2.0657592) q[0];
sx q[0];
rz(-1.540394) q[0];
rz(-1.5829344) q[2];
sx q[2];
rz(-1.7430787) q[2];
sx q[2];
rz(-2.7207295) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0145871) q[1];
sx q[1];
rz(-1.0795322) q[1];
sx q[1];
rz(2.7738357) q[1];
rz(-pi) q[2];
rz(-2.2035962) q[3];
sx q[3];
rz(-0.15060234) q[3];
sx q[3];
rz(1.6737398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5892443) q[2];
sx q[2];
rz(-1.9231223) q[2];
sx q[2];
rz(-0.98769665) q[2];
rz(2.6622631) q[3];
sx q[3];
rz(-1.597581) q[3];
sx q[3];
rz(-0.89232579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4722897) q[0];
sx q[0];
rz(-2.9030114) q[0];
sx q[0];
rz(-0.14078374) q[0];
rz(0.36551481) q[1];
sx q[1];
rz(-2.3194158) q[1];
sx q[1];
rz(-2.4929094) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4869944) q[0];
sx q[0];
rz(-0.93385044) q[0];
sx q[0];
rz(2.0084698) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6874403) q[2];
sx q[2];
rz(-2.3320227) q[2];
sx q[2];
rz(-2.4621689) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.89542298) q[1];
sx q[1];
rz(-1.514558) q[1];
sx q[1];
rz(-0.59631056) q[1];
x q[2];
rz(1.8126103) q[3];
sx q[3];
rz(-1.2529272) q[3];
sx q[3];
rz(2.1235025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6792128) q[2];
sx q[2];
rz(-1.0054192) q[2];
sx q[2];
rz(-3.0688378) q[2];
rz(-0.66271979) q[3];
sx q[3];
rz(-1.8279671) q[3];
sx q[3];
rz(-2.976118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0872021) q[0];
sx q[0];
rz(-1.6086171) q[0];
sx q[0];
rz(1.4679012) q[0];
rz(1.6064593) q[1];
sx q[1];
rz(-2.2538593) q[1];
sx q[1];
rz(-1.5182553) q[1];
rz(-2.9130878) q[2];
sx q[2];
rz(-1.6863556) q[2];
sx q[2];
rz(-0.40876331) q[2];
rz(-1.0853416) q[3];
sx q[3];
rz(-2.3997636) q[3];
sx q[3];
rz(2.6482481) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
