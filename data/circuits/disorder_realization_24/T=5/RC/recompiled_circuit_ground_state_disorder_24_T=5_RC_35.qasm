OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.4234023) q[0];
sx q[0];
rz(-0.0027522491) q[0];
sx q[0];
rz(-1.0116853) q[0];
rz(0.45971316) q[1];
sx q[1];
rz(4.4432321) q[1];
sx q[1];
rz(9.2530773) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.190024) q[0];
sx q[0];
rz(-2.5309238) q[0];
sx q[0];
rz(0.14279731) q[0];
rz(-0.19440513) q[2];
sx q[2];
rz(-1.462359) q[2];
sx q[2];
rz(-2.6284503) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5039035) q[1];
sx q[1];
rz(-1.9963963) q[1];
sx q[1];
rz(-2.5017159) q[1];
rz(-3.0875077) q[3];
sx q[3];
rz(-1.4776609) q[3];
sx q[3];
rz(-2.2786811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.95919886) q[2];
sx q[2];
rz(-0.48754075) q[2];
sx q[2];
rz(-1.7215151) q[2];
rz(-1.9567664) q[3];
sx q[3];
rz(-1.5337557) q[3];
sx q[3];
rz(1.0678631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.083953388) q[0];
sx q[0];
rz(-1.5204484) q[0];
sx q[0];
rz(0.49348304) q[0];
rz(-0.77589846) q[1];
sx q[1];
rz(-0.50965613) q[1];
sx q[1];
rz(2.6367771) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0579266) q[0];
sx q[0];
rz(-2.2424881) q[0];
sx q[0];
rz(0.85606411) q[0];
rz(-pi) q[1];
rz(-2.6146982) q[2];
sx q[2];
rz(-1.3119446) q[2];
sx q[2];
rz(-0.35824725) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9318092) q[1];
sx q[1];
rz(-1.9701013) q[1];
sx q[1];
rz(0.65815355) q[1];
rz(-pi) q[2];
x q[2];
rz(1.633058) q[3];
sx q[3];
rz(-1.8362507) q[3];
sx q[3];
rz(0.070814781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.33727553) q[2];
sx q[2];
rz(-1.8307468) q[2];
sx q[2];
rz(2.6709225) q[2];
rz(0.63052952) q[3];
sx q[3];
rz(-1.0258521) q[3];
sx q[3];
rz(-1.4001018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(3.0640963) q[0];
sx q[0];
rz(-3.016576) q[0];
sx q[0];
rz(0.65573829) q[0];
rz(-2.8302622) q[1];
sx q[1];
rz(-0.89404023) q[1];
sx q[1];
rz(-3.0701367) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27256672) q[0];
sx q[0];
rz(-1.2203958) q[0];
sx q[0];
rz(0.053215543) q[0];
x q[1];
rz(1.9391483) q[2];
sx q[2];
rz(-1.8873653) q[2];
sx q[2];
rz(-2.9376415) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1972547) q[1];
sx q[1];
rz(-1.5052915) q[1];
sx q[1];
rz(2.0599026) q[1];
x q[2];
rz(-2.7903665) q[3];
sx q[3];
rz(-2.5820508) q[3];
sx q[3];
rz(-0.21332394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.91325703) q[2];
sx q[2];
rz(-1.5568638) q[2];
sx q[2];
rz(-3.081591) q[2];
rz(1.9289198) q[3];
sx q[3];
rz(-2.4433177) q[3];
sx q[3];
rz(-0.76446271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5948828) q[0];
sx q[0];
rz(-1.8130274) q[0];
sx q[0];
rz(-2.5643964) q[0];
rz(2.3120841) q[1];
sx q[1];
rz(-1.0521768) q[1];
sx q[1];
rz(-2.3037691) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4234377) q[0];
sx q[0];
rz(-0.95036794) q[0];
sx q[0];
rz(-1.6654832) q[0];
x q[1];
rz(-0.98060645) q[2];
sx q[2];
rz(-1.5204805) q[2];
sx q[2];
rz(2.3286164) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.29463331) q[1];
sx q[1];
rz(-2.3851352) q[1];
sx q[1];
rz(-2.129617) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1076848) q[3];
sx q[3];
rz(-0.75124012) q[3];
sx q[3];
rz(0.43207016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.025658) q[2];
sx q[2];
rz(-1.6178774) q[2];
sx q[2];
rz(-1.9177829) q[2];
rz(0.4661679) q[3];
sx q[3];
rz(-0.40188447) q[3];
sx q[3];
rz(0.60552067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8849628) q[0];
sx q[0];
rz(-1.4361359) q[0];
sx q[0];
rz(-0.10109854) q[0];
rz(-0.413232) q[1];
sx q[1];
rz(-2.7006472) q[1];
sx q[1];
rz(-1.0879999) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1285217) q[0];
sx q[0];
rz(-2.5422342) q[0];
sx q[0];
rz(0.96512633) q[0];
x q[1];
rz(-1.1674983) q[2];
sx q[2];
rz(-2.3638862) q[2];
sx q[2];
rz(-2.4571927) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1697547) q[1];
sx q[1];
rz(-2.7248451) q[1];
sx q[1];
rz(-0.93081148) q[1];
rz(0.089121295) q[3];
sx q[3];
rz(-0.67135076) q[3];
sx q[3];
rz(-1.3446087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2130412) q[2];
sx q[2];
rz(-2.2553359) q[2];
sx q[2];
rz(1.5555596) q[2];
rz(-1.0726311) q[3];
sx q[3];
rz(-1.5242679) q[3];
sx q[3];
rz(0.13256375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.9685386) q[0];
sx q[0];
rz(-2.2569077) q[0];
sx q[0];
rz(1.435085) q[0];
rz(2.3834719) q[1];
sx q[1];
rz(-0.70223141) q[1];
sx q[1];
rz(-0.10163669) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91824965) q[0];
sx q[0];
rz(-2.0416284) q[0];
sx q[0];
rz(0.74145326) q[0];
rz(-2.0254214) q[2];
sx q[2];
rz(-1.897942) q[2];
sx q[2];
rz(-1.0516372) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.95127997) q[1];
sx q[1];
rz(-1.4644831) q[1];
sx q[1];
rz(-0.28369686) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2072565) q[3];
sx q[3];
rz(-1.1256256) q[3];
sx q[3];
rz(1.0220774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.3397843) q[2];
sx q[2];
rz(-2.0399751) q[2];
sx q[2];
rz(-1.3327117) q[2];
rz(-1.0821651) q[3];
sx q[3];
rz(-0.75282955) q[3];
sx q[3];
rz(1.7180721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4448755) q[0];
sx q[0];
rz(-1.2514021) q[0];
sx q[0];
rz(-2.3984997) q[0];
rz(-2.3785059) q[1];
sx q[1];
rz(-2.5021195) q[1];
sx q[1];
rz(2.2230164) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10906405) q[0];
sx q[0];
rz(-1.5957513) q[0];
sx q[0];
rz(3.1335013) q[0];
x q[1];
rz(2.1821676) q[2];
sx q[2];
rz(-0.55771962) q[2];
sx q[2];
rz(3.0447247) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.1700701) q[1];
sx q[1];
rz(-2.5405209) q[1];
sx q[1];
rz(2.5843589) q[1];
rz(2.4174446) q[3];
sx q[3];
rz(-1.6586496) q[3];
sx q[3];
rz(1.7334565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.7243728) q[2];
sx q[2];
rz(-1.1823187) q[2];
sx q[2];
rz(2.7057538) q[2];
rz(3.0271652) q[3];
sx q[3];
rz(-0.2468214) q[3];
sx q[3];
rz(-0.59823263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8082387) q[0];
sx q[0];
rz(-0.55753189) q[0];
sx q[0];
rz(2.5575141) q[0];
rz(0.43560478) q[1];
sx q[1];
rz(-1.4608811) q[1];
sx q[1];
rz(1.6216507) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10669691) q[0];
sx q[0];
rz(-0.55904065) q[0];
sx q[0];
rz(-0.99026545) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9233515) q[2];
sx q[2];
rz(-1.9710961) q[2];
sx q[2];
rz(1.2117653) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3585289) q[1];
sx q[1];
rz(-1.364768) q[1];
sx q[1];
rz(-1.5611783) q[1];
rz(-pi) q[2];
rz(-2.8826113) q[3];
sx q[3];
rz(-1.4121488) q[3];
sx q[3];
rz(2.4939031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0674151) q[2];
sx q[2];
rz(-1.7895074) q[2];
sx q[2];
rz(1.1161067) q[2];
rz(1.762278) q[3];
sx q[3];
rz(-1.1032871) q[3];
sx q[3];
rz(0.11837676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8922358) q[0];
sx q[0];
rz(-3.0638969) q[0];
sx q[0];
rz(0.78999162) q[0];
rz(3.0780011) q[1];
sx q[1];
rz(-0.60706943) q[1];
sx q[1];
rz(2.7008609) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4249464) q[0];
sx q[0];
rz(-1.2727203) q[0];
sx q[0];
rz(-1.9008725) q[0];
x q[1];
rz(-0.67040261) q[2];
sx q[2];
rz(-0.60302654) q[2];
sx q[2];
rz(0.5926026) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3778119) q[1];
sx q[1];
rz(-2.0449989) q[1];
sx q[1];
rz(-2.9567316) q[1];
rz(-pi) q[2];
rz(3.1001631) q[3];
sx q[3];
rz(-1.7600312) q[3];
sx q[3];
rz(-1.7570868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.26238394) q[2];
sx q[2];
rz(-0.79664207) q[2];
sx q[2];
rz(0.66514307) q[2];
rz(-2.3696259) q[3];
sx q[3];
rz(-0.77163458) q[3];
sx q[3];
rz(-1.6316679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7803698) q[0];
sx q[0];
rz(-2.4936115) q[0];
sx q[0];
rz(1.004647) q[0];
rz(-1.4077582) q[1];
sx q[1];
rz(-1.6561457) q[1];
sx q[1];
rz(1.2900603) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0697851) q[0];
sx q[0];
rz(-1.8172853) q[0];
sx q[0];
rz(-2.3355961) q[0];
rz(-pi) q[1];
rz(1.2783706) q[2];
sx q[2];
rz(-0.62033451) q[2];
sx q[2];
rz(2.8166556) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5651144) q[1];
sx q[1];
rz(-1.8449515) q[1];
sx q[1];
rz(1.3512011) q[1];
rz(-0.1741039) q[3];
sx q[3];
rz(-2.1645567) q[3];
sx q[3];
rz(-1.4623778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1349692) q[2];
sx q[2];
rz(-1.9776191) q[2];
sx q[2];
rz(-2.8912365) q[2];
rz(-0.97744673) q[3];
sx q[3];
rz(-1.2075295) q[3];
sx q[3];
rz(1.4704963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.1866495) q[0];
sx q[0];
rz(-1.432812) q[0];
sx q[0];
rz(-1.1695255) q[0];
rz(2.3969338) q[1];
sx q[1];
rz(-2.3458993) q[1];
sx q[1];
rz(0.69534272) q[1];
rz(-2.9982243) q[2];
sx q[2];
rz(-1.558254) q[2];
sx q[2];
rz(-2.8009453) q[2];
rz(-2.2779989) q[3];
sx q[3];
rz(-0.61932388) q[3];
sx q[3];
rz(-1.2854734) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
