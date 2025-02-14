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
rz(0.69005203) q[0];
sx q[0];
rz(-0.85059387) q[0];
sx q[0];
rz(0.2001065) q[0];
rz(-0.19835681) q[1];
sx q[1];
rz(5.7601647) q[1];
sx q[1];
rz(10.756607) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1855186) q[0];
sx q[0];
rz(-1.8022493) q[0];
sx q[0];
rz(-1.1197907) q[0];
rz(-pi) q[1];
rz(0.66352377) q[2];
sx q[2];
rz(-2.4276456) q[2];
sx q[2];
rz(0.51807846) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7847474) q[1];
sx q[1];
rz(-0.6797528) q[1];
sx q[1];
rz(3.0650782) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0588714) q[3];
sx q[3];
rz(-2.045407) q[3];
sx q[3];
rz(0.079291346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0693822) q[2];
sx q[2];
rz(-2.78077) q[2];
sx q[2];
rz(-1.8270095) q[2];
rz(2.2383111) q[3];
sx q[3];
rz(-1.7777781) q[3];
sx q[3];
rz(0.34590736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2717993) q[0];
sx q[0];
rz(-1.6558187) q[0];
sx q[0];
rz(-0.85579175) q[0];
rz(-0.77404147) q[1];
sx q[1];
rz(-1.1639405) q[1];
sx q[1];
rz(-2.3537297) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.60074) q[0];
sx q[0];
rz(-0.2517646) q[0];
sx q[0];
rz(-1.4815848) q[0];
rz(3.0237979) q[2];
sx q[2];
rz(-0.34644768) q[2];
sx q[2];
rz(-2.7730377) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.012177906) q[1];
sx q[1];
rz(-1.9099292) q[1];
sx q[1];
rz(-2.0268834) q[1];
rz(-pi) q[2];
rz(2.8773852) q[3];
sx q[3];
rz(-2.4572861) q[3];
sx q[3];
rz(0.58240376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0756125) q[2];
sx q[2];
rz(-1.2300666) q[2];
sx q[2];
rz(-2.7454929) q[2];
rz(-1.5710477) q[3];
sx q[3];
rz(-1.6377662) q[3];
sx q[3];
rz(1.3554696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.9539255) q[0];
sx q[0];
rz(-1.4952156) q[0];
sx q[0];
rz(-0.060610108) q[0];
rz(2.4568779) q[1];
sx q[1];
rz(-1.4002607) q[1];
sx q[1];
rz(-2.8025467) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7386757) q[0];
sx q[0];
rz(-2.5581723) q[0];
sx q[0];
rz(-0.03761272) q[0];
rz(2.5671447) q[2];
sx q[2];
rz(-0.22413218) q[2];
sx q[2];
rz(1.4719065) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0319034) q[1];
sx q[1];
rz(-1.2983783) q[1];
sx q[1];
rz(1.2485571) q[1];
x q[2];
rz(-3.1296785) q[3];
sx q[3];
rz(-2.3647226) q[3];
sx q[3];
rz(0.35038951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.48245779) q[2];
sx q[2];
rz(-1.8637916) q[2];
sx q[2];
rz(2.6711312) q[2];
rz(2.2792234) q[3];
sx q[3];
rz(-0.44895288) q[3];
sx q[3];
rz(1.5615777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7793133) q[0];
sx q[0];
rz(-1.1184432) q[0];
sx q[0];
rz(-0.40801868) q[0];
rz(-1.7904003) q[1];
sx q[1];
rz(-1.2974757) q[1];
sx q[1];
rz(1.321235) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6141339) q[0];
sx q[0];
rz(-1.2686994) q[0];
sx q[0];
rz(-2.6149349) q[0];
rz(3.0555326) q[2];
sx q[2];
rz(-2.2241908) q[2];
sx q[2];
rz(-1.6376405) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.8157573) q[1];
sx q[1];
rz(-1.4495295) q[1];
sx q[1];
rz(-1.497529) q[1];
rz(0.42709728) q[3];
sx q[3];
rz(-1.1833041) q[3];
sx q[3];
rz(2.1807293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9185751) q[2];
sx q[2];
rz(-1.7097946) q[2];
sx q[2];
rz(2.4187386) q[2];
rz(2.2919848) q[3];
sx q[3];
rz(-1.9690211) q[3];
sx q[3];
rz(2.1037219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9549114) q[0];
sx q[0];
rz(-0.51690042) q[0];
sx q[0];
rz(-0.55996672) q[0];
rz(2.9365183) q[1];
sx q[1];
rz(-1.2566902) q[1];
sx q[1];
rz(-2.8505039) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3223461) q[0];
sx q[0];
rz(-0.3194321) q[0];
sx q[0];
rz(-2.5752905) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5076324) q[2];
sx q[2];
rz(-0.77253714) q[2];
sx q[2];
rz(2.87839) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4704125) q[1];
sx q[1];
rz(-0.7268097) q[1];
sx q[1];
rz(-1.8720759) q[1];
rz(-1.7311312) q[3];
sx q[3];
rz(-2.3162127) q[3];
sx q[3];
rz(1.3813409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2228955) q[2];
sx q[2];
rz(-1.6632068) q[2];
sx q[2];
rz(2.843294) q[2];
rz(-1.1693303) q[3];
sx q[3];
rz(-0.9669286) q[3];
sx q[3];
rz(-1.6303308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30763141) q[0];
sx q[0];
rz(-1.6809604) q[0];
sx q[0];
rz(-2.5391915) q[0];
rz(2.7700453) q[1];
sx q[1];
rz(-1.5967775) q[1];
sx q[1];
rz(2.1098302) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69799549) q[0];
sx q[0];
rz(-1.5549608) q[0];
sx q[0];
rz(1.3965194) q[0];
rz(2.736892) q[2];
sx q[2];
rz(-1.7621303) q[2];
sx q[2];
rz(1.2197242) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6556485) q[1];
sx q[1];
rz(-1.3820501) q[1];
sx q[1];
rz(-2.6552378) q[1];
rz(-pi) q[2];
rz(0.57547456) q[3];
sx q[3];
rz(-2.1361975) q[3];
sx q[3];
rz(2.7421212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3370257) q[2];
sx q[2];
rz(-2.0402543) q[2];
sx q[2];
rz(1.5186914) q[2];
rz(2.2931781) q[3];
sx q[3];
rz(-2.266326) q[3];
sx q[3];
rz(3.1019822) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1161716) q[0];
sx q[0];
rz(-3.0905368) q[0];
sx q[0];
rz(0.99789944) q[0];
rz(-0.2746703) q[1];
sx q[1];
rz(-0.88013595) q[1];
sx q[1];
rz(1.3708699) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8557963) q[0];
sx q[0];
rz(-1.5842452) q[0];
sx q[0];
rz(-2.5105259) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9657126) q[2];
sx q[2];
rz(-0.38287258) q[2];
sx q[2];
rz(-1.9667786) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.26105598) q[1];
sx q[1];
rz(-0.66268259) q[1];
sx q[1];
rz(-0.39775325) q[1];
rz(-2.3418535) q[3];
sx q[3];
rz(-2.1465786) q[3];
sx q[3];
rz(0.48279253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2844598) q[2];
sx q[2];
rz(-2.1993115) q[2];
sx q[2];
rz(3.1119463) q[2];
rz(-0.61521411) q[3];
sx q[3];
rz(-1.3968202) q[3];
sx q[3];
rz(0.342338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0520332) q[0];
sx q[0];
rz(-2.679306) q[0];
sx q[0];
rz(-2.3585368) q[0];
rz(1.1198593) q[1];
sx q[1];
rz(-2.4744787) q[1];
sx q[1];
rz(1.7436183) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.072294839) q[0];
sx q[0];
rz(-1.1279449) q[0];
sx q[0];
rz(-2.686727) q[0];
x q[1];
rz(-0.88770788) q[2];
sx q[2];
rz(-0.36044932) q[2];
sx q[2];
rz(-0.53149022) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3191171) q[1];
sx q[1];
rz(-0.90293316) q[1];
sx q[1];
rz(0.0112337) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.21381883) q[3];
sx q[3];
rz(-2.345746) q[3];
sx q[3];
rz(0.70094943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.29584259) q[2];
sx q[2];
rz(-2.5551899) q[2];
sx q[2];
rz(-1.3207377) q[2];
rz(-2.1218421) q[3];
sx q[3];
rz(-1.4239862) q[3];
sx q[3];
rz(1.4199055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19972292) q[0];
sx q[0];
rz(-2.0482735) q[0];
sx q[0];
rz(0.95440188) q[0];
rz(1.154254) q[1];
sx q[1];
rz(-0.64607611) q[1];
sx q[1];
rz(-2.0571713) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0592151) q[0];
sx q[0];
rz(-2.431916) q[0];
sx q[0];
rz(-1.6355455) q[0];
x q[1];
rz(-2.886359) q[2];
sx q[2];
rz(-1.298578) q[2];
sx q[2];
rz(1.24025) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3855648) q[1];
sx q[1];
rz(-1.6431019) q[1];
sx q[1];
rz(-1.0977618) q[1];
rz(-pi) q[2];
rz(1.3865269) q[3];
sx q[3];
rz(-2.3109155) q[3];
sx q[3];
rz(-0.6669464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0584917) q[2];
sx q[2];
rz(-0.19187555) q[2];
sx q[2];
rz(2.7044738) q[2];
rz(1.3433749) q[3];
sx q[3];
rz(-1.790204) q[3];
sx q[3];
rz(-0.04537151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0992391) q[0];
sx q[0];
rz(-1.0004685) q[0];
sx q[0];
rz(1.4105256) q[0];
rz(-0.22077416) q[1];
sx q[1];
rz(-2.6103554) q[1];
sx q[1];
rz(-2.208362) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8338884) q[0];
sx q[0];
rz(-1.5844463) q[0];
sx q[0];
rz(0.74539124) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2104697) q[2];
sx q[2];
rz(-2.1585352) q[2];
sx q[2];
rz(2.7992942) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2360742) q[1];
sx q[1];
rz(-2.337444) q[1];
sx q[1];
rz(-2.9574636) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0296928) q[3];
sx q[3];
rz(-0.94116941) q[3];
sx q[3];
rz(0.96523413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.71020469) q[2];
sx q[2];
rz(-0.34222558) q[2];
sx q[2];
rz(2.5101275) q[2];
rz(-2.4225875) q[3];
sx q[3];
rz(-1.3471666) q[3];
sx q[3];
rz(-0.24148153) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63856335) q[0];
sx q[0];
rz(-1.6521492) q[0];
sx q[0];
rz(-1.8304992) q[0];
rz(-0.44166625) q[1];
sx q[1];
rz(-1.3413981) q[1];
sx q[1];
rz(-1.9955019) q[1];
rz(-1.491811) q[2];
sx q[2];
rz(-0.38182237) q[2];
sx q[2];
rz(0.67425722) q[2];
rz(1.4680223) q[3];
sx q[3];
rz(-1.0036052) q[3];
sx q[3];
rz(-2.6617692) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
