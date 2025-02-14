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
rz(3.1004768) q[0];
sx q[0];
rz(-1.2007204) q[0];
sx q[0];
rz(1.4708941) q[0];
rz(-0.38493758) q[1];
sx q[1];
rz(2.5632783) q[1];
sx q[1];
rz(10.332528) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0931514) q[0];
sx q[0];
rz(-2.127383) q[0];
sx q[0];
rz(2.0970035) q[0];
rz(-pi) q[1];
rz(-0.29405354) q[2];
sx q[2];
rz(-1.8584677) q[2];
sx q[2];
rz(1.8881544) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5564087) q[1];
sx q[1];
rz(-2.6165462) q[1];
sx q[1];
rz(-0.37918143) q[1];
x q[2];
rz(3.0201245) q[3];
sx q[3];
rz(-1.6252015) q[3];
sx q[3];
rz(-1.3682981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3140728) q[2];
sx q[2];
rz(-1.4853518) q[2];
sx q[2];
rz(0.32076389) q[2];
rz(3.1406506) q[3];
sx q[3];
rz(-0.92206803) q[3];
sx q[3];
rz(2.6607386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28474057) q[0];
sx q[0];
rz(-2.4212403) q[0];
sx q[0];
rz(0.62069935) q[0];
rz(-1.8547828) q[1];
sx q[1];
rz(-1.4871253) q[1];
sx q[1];
rz(0.99498814) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33546977) q[0];
sx q[0];
rz(-0.54278446) q[0];
sx q[0];
rz(1.1837028) q[0];
rz(-0.95717818) q[2];
sx q[2];
rz(-0.52268302) q[2];
sx q[2];
rz(-1.3764639) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7409119) q[1];
sx q[1];
rz(-1.3865292) q[1];
sx q[1];
rz(-2.0234985) q[1];
x q[2];
rz(-0.80886474) q[3];
sx q[3];
rz(-1.1877664) q[3];
sx q[3];
rz(-1.3749141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2113125) q[2];
sx q[2];
rz(-1.2418396) q[2];
sx q[2];
rz(-2.9105183) q[2];
rz(-1.5288345) q[3];
sx q[3];
rz(-1.4569747) q[3];
sx q[3];
rz(2.5902364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.424161) q[0];
sx q[0];
rz(-1.9743974) q[0];
sx q[0];
rz(0.396808) q[0];
rz(-1.1948168) q[1];
sx q[1];
rz(-1.2620986) q[1];
sx q[1];
rz(-1.4746812) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88746027) q[0];
sx q[0];
rz(-1.4829221) q[0];
sx q[0];
rz(-0.011482419) q[0];
rz(-pi) q[1];
rz(-2.4184995) q[2];
sx q[2];
rz(-2.4567219) q[2];
sx q[2];
rz(-1.9346646) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.84473064) q[1];
sx q[1];
rz(-0.76557388) q[1];
sx q[1];
rz(-1.0699478) q[1];
rz(3.1297333) q[3];
sx q[3];
rz(-1.1766772) q[3];
sx q[3];
rz(0.32255641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.84612238) q[2];
sx q[2];
rz(-2.2961871) q[2];
sx q[2];
rz(1.4836614) q[2];
rz(-2.1814003) q[3];
sx q[3];
rz(-0.63513297) q[3];
sx q[3];
rz(2.4471762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.33261499) q[0];
sx q[0];
rz(-1.3354477) q[0];
sx q[0];
rz(-0.72840011) q[0];
rz(1.2737466) q[1];
sx q[1];
rz(-1.1796917) q[1];
sx q[1];
rz(-1.5789998) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.584986) q[0];
sx q[0];
rz(-1.5674563) q[0];
sx q[0];
rz(-1.5933611) q[0];
x q[1];
rz(2.118763) q[2];
sx q[2];
rz(-1.5336516) q[2];
sx q[2];
rz(1.9859344) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.77550239) q[1];
sx q[1];
rz(-1.3555158) q[1];
sx q[1];
rz(-1.7083113) q[1];
x q[2];
rz(0.32890851) q[3];
sx q[3];
rz(-2.2136627) q[3];
sx q[3];
rz(2.7907284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1100715) q[2];
sx q[2];
rz(-1.4036125) q[2];
sx q[2];
rz(0.20315476) q[2];
rz(-1.3567694) q[3];
sx q[3];
rz(-1.2099096) q[3];
sx q[3];
rz(-1.8206966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88306952) q[0];
sx q[0];
rz(-1.8061545) q[0];
sx q[0];
rz(-1.6254599) q[0];
rz(-2.7698703) q[1];
sx q[1];
rz(-0.55551353) q[1];
sx q[1];
rz(2.1518167) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55732841) q[0];
sx q[0];
rz(-0.66270486) q[0];
sx q[0];
rz(1.8813547) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11595102) q[2];
sx q[2];
rz(-0.77590307) q[2];
sx q[2];
rz(2.3052892) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6367148) q[1];
sx q[1];
rz(-1.9002007) q[1];
sx q[1];
rz(0.18784951) q[1];
rz(-1.2618213) q[3];
sx q[3];
rz(-1.1011964) q[3];
sx q[3];
rz(-2.0275786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.524579) q[2];
sx q[2];
rz(-2.9144574) q[2];
sx q[2];
rz(3.0830834) q[2];
rz(1.8578364) q[3];
sx q[3];
rz(-2.3144898) q[3];
sx q[3];
rz(-1.1684928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7002895) q[0];
sx q[0];
rz(-2.0952201) q[0];
sx q[0];
rz(0.48536479) q[0];
rz(-2.6742588) q[1];
sx q[1];
rz(-2.554437) q[1];
sx q[1];
rz(-2.4553518) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8507535) q[0];
sx q[0];
rz(-0.79054442) q[0];
sx q[0];
rz(0.45773069) q[0];
rz(-pi) q[1];
rz(2.9250018) q[2];
sx q[2];
rz(-0.95568839) q[2];
sx q[2];
rz(-1.5768346) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8279449) q[1];
sx q[1];
rz(-1.0961431) q[1];
sx q[1];
rz(-2.5728435) q[1];
rz(2.521311) q[3];
sx q[3];
rz(-0.93343319) q[3];
sx q[3];
rz(-0.31951527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.83963362) q[2];
sx q[2];
rz(-2.0531824) q[2];
sx q[2];
rz(-2.9097617) q[2];
rz(-2.2202282) q[3];
sx q[3];
rz(-1.5270343) q[3];
sx q[3];
rz(0.025402633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63474083) q[0];
sx q[0];
rz(-2.1379819) q[0];
sx q[0];
rz(-2.4941709) q[0];
rz(2.1000775) q[1];
sx q[1];
rz(-1.2728649) q[1];
sx q[1];
rz(-2.0288859) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0395687) q[0];
sx q[0];
rz(-1.0728076) q[0];
sx q[0];
rz(2.6302393) q[0];
rz(-pi) q[1];
rz(3.0810705) q[2];
sx q[2];
rz(-1.4281445) q[2];
sx q[2];
rz(1.5885274) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3379593) q[1];
sx q[1];
rz(-2.5542959) q[1];
sx q[1];
rz(-0.31897591) q[1];
x q[2];
rz(1.5912368) q[3];
sx q[3];
rz(-2.6764565) q[3];
sx q[3];
rz(0.49275574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.57704321) q[2];
sx q[2];
rz(-0.81114045) q[2];
sx q[2];
rz(-2.2570293) q[2];
rz(2.1029643) q[3];
sx q[3];
rz(-2.0197155) q[3];
sx q[3];
rz(1.5326726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7811979) q[0];
sx q[0];
rz(-2.9304657) q[0];
sx q[0];
rz(-1.1398844) q[0];
rz(-2.3431011) q[1];
sx q[1];
rz(-1.4356828) q[1];
sx q[1];
rz(-3.0671157) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7388106) q[0];
sx q[0];
rz(-0.96091369) q[0];
sx q[0];
rz(-2.8426811) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3788667) q[2];
sx q[2];
rz(-1.8755091) q[2];
sx q[2];
rz(-0.21108) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.18596126) q[1];
sx q[1];
rz(-0.73195405) q[1];
sx q[1];
rz(0.47853985) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2747669) q[3];
sx q[3];
rz(-0.4496963) q[3];
sx q[3];
rz(-1.7148643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8188339) q[2];
sx q[2];
rz(-0.61046159) q[2];
sx q[2];
rz(2.9739213) q[2];
rz(-0.75336114) q[3];
sx q[3];
rz(-1.2246917) q[3];
sx q[3];
rz(2.0577478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3566256) q[0];
sx q[0];
rz(-1.6870808) q[0];
sx q[0];
rz(-2.0941358) q[0];
rz(-2.6197703) q[1];
sx q[1];
rz(-1.1818886) q[1];
sx q[1];
rz(2.3375915) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8270644) q[0];
sx q[0];
rz(-2.6834163) q[0];
sx q[0];
rz(1.1660775) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.1754403) q[2];
sx q[2];
rz(-1.0217624) q[2];
sx q[2];
rz(0.4859449) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9404133) q[1];
sx q[1];
rz(-2.6087954) q[1];
sx q[1];
rz(-2.6032107) q[1];
rz(-1.9639765) q[3];
sx q[3];
rz(-2.2344974) q[3];
sx q[3];
rz(0.8908602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1890586) q[2];
sx q[2];
rz(-1.4905832) q[2];
sx q[2];
rz(-2.789759) q[2];
rz(3.0686038) q[3];
sx q[3];
rz(-1.9609541) q[3];
sx q[3];
rz(0.73968691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9214582) q[0];
sx q[0];
rz(-2.298482) q[0];
sx q[0];
rz(1.5716918) q[0];
rz(0.72228557) q[1];
sx q[1];
rz(-1.1415569) q[1];
sx q[1];
rz(1.7438181) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36436468) q[0];
sx q[0];
rz(-1.1684864) q[0];
sx q[0];
rz(1.2978294) q[0];
rz(-pi) q[1];
rz(1.3438002) q[2];
sx q[2];
rz(-0.8778866) q[2];
sx q[2];
rz(-0.60688144) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7745516) q[1];
sx q[1];
rz(-1.4907717) q[1];
sx q[1];
rz(-2.5888279) q[1];
x q[2];
rz(1.151284) q[3];
sx q[3];
rz(-1.691201) q[3];
sx q[3];
rz(-1.6199978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7950644) q[2];
sx q[2];
rz(-1.3115735) q[2];
sx q[2];
rz(-1.2947003) q[2];
rz(-2.0857701) q[3];
sx q[3];
rz(-2.081223) q[3];
sx q[3];
rz(2.136371) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1718564) q[0];
sx q[0];
rz(-2.2315401) q[0];
sx q[0];
rz(1.6730614) q[0];
rz(-0.83364529) q[1];
sx q[1];
rz(-1.2794762) q[1];
sx q[1];
rz(1.5247482) q[1];
rz(0.8236105) q[2];
sx q[2];
rz(-1.2678185) q[2];
sx q[2];
rz(2.5662255) q[2];
rz(-0.81832073) q[3];
sx q[3];
rz(-2.5723895) q[3];
sx q[3];
rz(1.7583917) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
