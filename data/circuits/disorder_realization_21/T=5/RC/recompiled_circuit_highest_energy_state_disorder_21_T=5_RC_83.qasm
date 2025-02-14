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
rz(1.1462829) q[0];
sx q[0];
rz(-3.1132071) q[0];
sx q[0];
rz(2.1838768) q[0];
rz(-2.7081642) q[1];
sx q[1];
rz(-1.6038916) q[1];
sx q[1];
rz(-2.5133207) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7564297) q[0];
sx q[0];
rz(-1.9397501) q[0];
sx q[0];
rz(-1.4403309) q[0];
x q[1];
rz(-2.069491) q[2];
sx q[2];
rz(-1.9737502) q[2];
sx q[2];
rz(-3.0080331) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.5124197) q[1];
sx q[1];
rz(-2.3827374) q[1];
sx q[1];
rz(-2.0461798) q[1];
x q[2];
rz(0.099277012) q[3];
sx q[3];
rz(-1.160272) q[3];
sx q[3];
rz(-1.6308189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6750703) q[2];
sx q[2];
rz(-2.1212981) q[2];
sx q[2];
rz(-2.7543219) q[2];
rz(-1.9668503) q[3];
sx q[3];
rz(-1.6956804) q[3];
sx q[3];
rz(2.2430879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77884498) q[0];
sx q[0];
rz(-1.7411106) q[0];
sx q[0];
rz(-1.269302) q[0];
rz(-1.6954039) q[1];
sx q[1];
rz(-1.2934928) q[1];
sx q[1];
rz(2.2206025) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35996485) q[0];
sx q[0];
rz(-2.4422997) q[0];
sx q[0];
rz(-1.101732) q[0];
x q[1];
rz(-1.3878421) q[2];
sx q[2];
rz(-2.2277846) q[2];
sx q[2];
rz(-0.25568257) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2650406) q[1];
sx q[1];
rz(-1.3807793) q[1];
sx q[1];
rz(-1.9273619) q[1];
x q[2];
rz(0.22590206) q[3];
sx q[3];
rz(-1.5776565) q[3];
sx q[3];
rz(-0.021319162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9903119) q[2];
sx q[2];
rz(-1.9419779) q[2];
sx q[2];
rz(-0.82378236) q[2];
rz(-1.524205) q[3];
sx q[3];
rz(-1.5533605) q[3];
sx q[3];
rz(-1.2113781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23564944) q[0];
sx q[0];
rz(-2.2610569) q[0];
sx q[0];
rz(-0.5916416) q[0];
rz(0.94332424) q[1];
sx q[1];
rz(-2.0843518) q[1];
sx q[1];
rz(1.0234157) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8013894) q[0];
sx q[0];
rz(-1.4638299) q[0];
sx q[0];
rz(1.68856) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7696174) q[2];
sx q[2];
rz(-1.0468654) q[2];
sx q[2];
rz(1.9934479) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8152961) q[1];
sx q[1];
rz(-0.097541428) q[1];
sx q[1];
rz(1.8952712) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.39822762) q[3];
sx q[3];
rz(-2.4169528) q[3];
sx q[3];
rz(-0.46996024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8817875) q[2];
sx q[2];
rz(-1.2299579) q[2];
sx q[2];
rz(-1.7477431) q[2];
rz(-1.4272089) q[3];
sx q[3];
rz(-0.87427679) q[3];
sx q[3];
rz(-2.8425596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3706936) q[0];
sx q[0];
rz(-2.735266) q[0];
sx q[0];
rz(2.4942177) q[0];
rz(2.8992843) q[1];
sx q[1];
rz(-1.7055885) q[1];
sx q[1];
rz(-1.4586331) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1636617) q[0];
sx q[0];
rz(-0.35172281) q[0];
sx q[0];
rz(-0.39885421) q[0];
rz(-pi) q[1];
rz(-2.3056612) q[2];
sx q[2];
rz(-2.3337337) q[2];
sx q[2];
rz(-0.69617803) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.35136552) q[1];
sx q[1];
rz(-2.0560762) q[1];
sx q[1];
rz(-0.041260551) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.69497739) q[3];
sx q[3];
rz(-2.2314592) q[3];
sx q[3];
rz(1.4042749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.709939) q[2];
sx q[2];
rz(-0.84448758) q[2];
sx q[2];
rz(-1.9169774) q[2];
rz(-2.8905408) q[3];
sx q[3];
rz(-0.71728388) q[3];
sx q[3];
rz(2.203598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5994023) q[0];
sx q[0];
rz(-1.1325862) q[0];
sx q[0];
rz(2.9436924) q[0];
rz(1.9056162) q[1];
sx q[1];
rz(-0.37207741) q[1];
sx q[1];
rz(3.1370251) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76659996) q[0];
sx q[0];
rz(-2.1439634) q[0];
sx q[0];
rz(-1.8660924) q[0];
rz(-1.0367583) q[2];
sx q[2];
rz(-1.4994086) q[2];
sx q[2];
rz(-0.56278261) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.98036042) q[1];
sx q[1];
rz(-1.2975818) q[1];
sx q[1];
rz(1.7262162) q[1];
x q[2];
rz(-1.3602348) q[3];
sx q[3];
rz(-1.9145085) q[3];
sx q[3];
rz(0.76219073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.68134754) q[2];
sx q[2];
rz(-2.6474417) q[2];
sx q[2];
rz(-0.34026185) q[2];
rz(-1.6081238) q[3];
sx q[3];
rz(-1.3932649) q[3];
sx q[3];
rz(1.6218012) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22572868) q[0];
sx q[0];
rz(-0.69750834) q[0];
sx q[0];
rz(1.1874143) q[0];
rz(-1.4215218) q[1];
sx q[1];
rz(-0.41821304) q[1];
sx q[1];
rz(0.13030599) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1811838) q[0];
sx q[0];
rz(-1.6941923) q[0];
sx q[0];
rz(0.058870319) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2537986) q[2];
sx q[2];
rz(-2.0468132) q[2];
sx q[2];
rz(0.43080518) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9355311) q[1];
sx q[1];
rz(-1.5714679) q[1];
sx q[1];
rz(2.9737581) q[1];
x q[2];
rz(1.6076419) q[3];
sx q[3];
rz(-0.26290694) q[3];
sx q[3];
rz(-2.2927566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.80798951) q[2];
sx q[2];
rz(-1.5895867) q[2];
sx q[2];
rz(1.0398593) q[2];
rz(0.17397675) q[3];
sx q[3];
rz(-2.0128553) q[3];
sx q[3];
rz(-2.4719293) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.533605) q[0];
sx q[0];
rz(-1.0796115) q[0];
sx q[0];
rz(-2.3002891) q[0];
rz(1.816642) q[1];
sx q[1];
rz(-0.94580301) q[1];
sx q[1];
rz(1.4680877) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3260088) q[0];
sx q[0];
rz(-1.2213314) q[0];
sx q[0];
rz(-0.38527003) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.011558455) q[2];
sx q[2];
rz(-2.0894139) q[2];
sx q[2];
rz(-2.8457038) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0206179) q[1];
sx q[1];
rz(-2.1281895) q[1];
sx q[1];
rz(-2.1305898) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2135336) q[3];
sx q[3];
rz(-1.284984) q[3];
sx q[3];
rz(0.88102007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2738721) q[2];
sx q[2];
rz(-2.2404859) q[2];
sx q[2];
rz(-2.4967398) q[2];
rz(-1.0817179) q[3];
sx q[3];
rz(-2.7223301) q[3];
sx q[3];
rz(-0.7209512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.69407392) q[0];
sx q[0];
rz(-2.9677291) q[0];
sx q[0];
rz(-0.41282594) q[0];
rz(1.6089571) q[1];
sx q[1];
rz(-0.91333476) q[1];
sx q[1];
rz(0.70721165) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6462881) q[0];
sx q[0];
rz(-2.2997059) q[0];
sx q[0];
rz(-1.1524617) q[0];
rz(0.17036713) q[2];
sx q[2];
rz(-2.3922386) q[2];
sx q[2];
rz(-0.10766497) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.78145786) q[1];
sx q[1];
rz(-2.4943922) q[1];
sx q[1];
rz(-0.57768808) q[1];
rz(-pi) q[2];
rz(-1.3363289) q[3];
sx q[3];
rz(-2.5904561) q[3];
sx q[3];
rz(-0.41195991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1327208) q[2];
sx q[2];
rz(-0.29942313) q[2];
sx q[2];
rz(2.6297165) q[2];
rz(-0.68073186) q[3];
sx q[3];
rz(-1.3635037) q[3];
sx q[3];
rz(2.9752922) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7452288) q[0];
sx q[0];
rz(-1.5928716) q[0];
sx q[0];
rz(-2.6863099) q[0];
rz(-2.0782754) q[1];
sx q[1];
rz(-0.89439193) q[1];
sx q[1];
rz(-3.0696226) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44350478) q[0];
sx q[0];
rz(-3.1060954) q[0];
sx q[0];
rz(0.3548236) q[0];
rz(-pi) q[1];
rz(-2.7006702) q[2];
sx q[2];
rz(-1.5498232) q[2];
sx q[2];
rz(1.6787337) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.2769952) q[1];
sx q[1];
rz(-2.2114843) q[1];
sx q[1];
rz(2.6990141) q[1];
x q[2];
rz(-2.2293363) q[3];
sx q[3];
rz(-2.7317762) q[3];
sx q[3];
rz(-1.2959105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3528184) q[2];
sx q[2];
rz(-2.7249551) q[2];
sx q[2];
rz(-2.519506) q[2];
rz(-0.82428733) q[3];
sx q[3];
rz(-1.0930748) q[3];
sx q[3];
rz(-0.85339671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.033009919) q[0];
sx q[0];
rz(-2.0350631) q[0];
sx q[0];
rz(2.1906817) q[0];
rz(1.7225601) q[1];
sx q[1];
rz(-2.6585237) q[1];
sx q[1];
rz(-2.5249265) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71910673) q[0];
sx q[0];
rz(-0.93084836) q[0];
sx q[0];
rz(1.6011222) q[0];
rz(-pi) q[1];
rz(-0.51787106) q[2];
sx q[2];
rz(-2.4642337) q[2];
sx q[2];
rz(-0.022938722) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2211799) q[1];
sx q[1];
rz(-1.3592459) q[1];
sx q[1];
rz(-1.5679564) q[1];
rz(-pi) q[2];
rz(2.4973013) q[3];
sx q[3];
rz(-0.81874311) q[3];
sx q[3];
rz(-0.14547023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0906494) q[2];
sx q[2];
rz(-1.1530777) q[2];
sx q[2];
rz(2.0743745) q[2];
rz(2.0459335) q[3];
sx q[3];
rz(-1.9039543) q[3];
sx q[3];
rz(0.9637951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.511956) q[0];
sx q[0];
rz(-2.1558599) q[0];
sx q[0];
rz(2.070367) q[0];
rz(-1.6596972) q[1];
sx q[1];
rz(-1.0922468) q[1];
sx q[1];
rz(-2.560871) q[1];
rz(-1.5459677) q[2];
sx q[2];
rz(-1.3083184) q[2];
sx q[2];
rz(2.0561894) q[2];
rz(0.4023424) q[3];
sx q[3];
rz(-0.33807031) q[3];
sx q[3];
rz(-2.8344179) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
