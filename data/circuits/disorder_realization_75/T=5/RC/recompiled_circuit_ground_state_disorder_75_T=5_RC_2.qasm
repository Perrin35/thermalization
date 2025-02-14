OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6899941) q[0];
sx q[0];
rz(-2.8328083) q[0];
sx q[0];
rz(-2.9017594) q[0];
rz(-0.56107768) q[1];
sx q[1];
rz(-2.3993888) q[1];
sx q[1];
rz(-1.5129369) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5133249) q[0];
sx q[0];
rz(-1.6986934) q[0];
sx q[0];
rz(1.2012175) q[0];
x q[1];
rz(-2.0432908) q[2];
sx q[2];
rz(-1.3368946) q[2];
sx q[2];
rz(2.5462674) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1354243) q[1];
sx q[1];
rz(-2.0574155) q[1];
sx q[1];
rz(1.7371337) q[1];
x q[2];
rz(1.2414819) q[3];
sx q[3];
rz(-2.9639177) q[3];
sx q[3];
rz(2.2840281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1538887) q[2];
sx q[2];
rz(-2.14553) q[2];
sx q[2];
rz(-1.055701) q[2];
rz(-0.40134564) q[3];
sx q[3];
rz(-1.515712) q[3];
sx q[3];
rz(1.6373985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5098679) q[0];
sx q[0];
rz(-0.26917502) q[0];
sx q[0];
rz(-1.4058231) q[0];
rz(-2.3954605) q[1];
sx q[1];
rz(-1.965062) q[1];
sx q[1];
rz(-0.064780898) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61363214) q[0];
sx q[0];
rz(-2.0191231) q[0];
sx q[0];
rz(2.0125858) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9322036) q[2];
sx q[2];
rz(-0.69062606) q[2];
sx q[2];
rz(-2.2356981) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.83103115) q[1];
sx q[1];
rz(-1.6509202) q[1];
sx q[1];
rz(-1.7145654) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4218036) q[3];
sx q[3];
rz(-2.9284366) q[3];
sx q[3];
rz(-3.0557291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.88227162) q[2];
sx q[2];
rz(-0.52053014) q[2];
sx q[2];
rz(-0.041291324) q[2];
rz(-1.1193554) q[3];
sx q[3];
rz(-1.3675523) q[3];
sx q[3];
rz(1.1411427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.816788) q[0];
sx q[0];
rz(-2.7140706) q[0];
sx q[0];
rz(-2.4422755) q[0];
rz(-0.62272561) q[1];
sx q[1];
rz(-2.3284349) q[1];
sx q[1];
rz(-0.25845382) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9759443) q[0];
sx q[0];
rz(-1.5830402) q[0];
sx q[0];
rz(-0.6338288) q[0];
x q[1];
rz(-1.7202366) q[2];
sx q[2];
rz(-2.4838313) q[2];
sx q[2];
rz(-2.2542816) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.093100637) q[1];
sx q[1];
rz(-0.44934595) q[1];
sx q[1];
rz(-0.046358422) q[1];
rz(-pi) q[2];
rz(-2.0007802) q[3];
sx q[3];
rz(-1.5376435) q[3];
sx q[3];
rz(1.5371109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.93566018) q[2];
sx q[2];
rz(-1.460133) q[2];
sx q[2];
rz(0.70651954) q[2];
rz(-0.62890729) q[3];
sx q[3];
rz(-2.3157412) q[3];
sx q[3];
rz(-1.1227054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.127447) q[0];
sx q[0];
rz(-1.4268459) q[0];
sx q[0];
rz(-0.97417796) q[0];
rz(-1.6575419) q[1];
sx q[1];
rz(-1.0933135) q[1];
sx q[1];
rz(-2.1407703) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38353048) q[0];
sx q[0];
rz(-2.5208726) q[0];
sx q[0];
rz(-1.8412983) q[0];
rz(0.8128667) q[2];
sx q[2];
rz(-2.1846601) q[2];
sx q[2];
rz(-2.6809106) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9467869) q[1];
sx q[1];
rz(-2.0344224) q[1];
sx q[1];
rz(-2.5083816) q[1];
rz(-pi) q[2];
x q[2];
rz(2.820904) q[3];
sx q[3];
rz(-1.5237892) q[3];
sx q[3];
rz(-0.21002029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.68690825) q[2];
sx q[2];
rz(-1.3409706) q[2];
sx q[2];
rz(-2.3243813) q[2];
rz(2.5456083) q[3];
sx q[3];
rz(-1.7242566) q[3];
sx q[3];
rz(-1.3982841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3572094) q[0];
sx q[0];
rz(-1.7359808) q[0];
sx q[0];
rz(-1.0719517) q[0];
rz(1.1373854) q[1];
sx q[1];
rz(-1.5147361) q[1];
sx q[1];
rz(0.17328182) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17006341) q[0];
sx q[0];
rz(-1.1779629) q[0];
sx q[0];
rz(2.6512572) q[0];
x q[1];
rz(2.2229175) q[2];
sx q[2];
rz(-0.5917509) q[2];
sx q[2];
rz(2.7482035) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3751602) q[1];
sx q[1];
rz(-0.81534895) q[1];
sx q[1];
rz(-0.098618193) q[1];
x q[2];
rz(0.16290476) q[3];
sx q[3];
rz(-1.4789733) q[3];
sx q[3];
rz(-0.11549982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5002354) q[2];
sx q[2];
rz(-0.45318979) q[2];
sx q[2];
rz(-0.87453169) q[2];
rz(2.4747961) q[3];
sx q[3];
rz(-1.4614481) q[3];
sx q[3];
rz(2.3165406) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85103971) q[0];
sx q[0];
rz(-2.9930826) q[0];
sx q[0];
rz(0.17459757) q[0];
rz(-2.5945276) q[1];
sx q[1];
rz(-2.6262296) q[1];
sx q[1];
rz(-1.863716) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39659835) q[0];
sx q[0];
rz(-1.6327259) q[0];
sx q[0];
rz(1.8577736) q[0];
rz(-pi) q[1];
rz(-0.49680423) q[2];
sx q[2];
rz(-2.3251109) q[2];
sx q[2];
rz(-1.7330488) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6656832) q[1];
sx q[1];
rz(-1.2194677) q[1];
sx q[1];
rz(0.76900469) q[1];
rz(-pi) q[2];
rz(2.5310263) q[3];
sx q[3];
rz(-0.99733099) q[3];
sx q[3];
rz(-1.2966882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3829019) q[2];
sx q[2];
rz(-0.26472696) q[2];
sx q[2];
rz(-2.7721789) q[2];
rz(-0.82768011) q[3];
sx q[3];
rz(-1.8281432) q[3];
sx q[3];
rz(0.15414342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5740042) q[0];
sx q[0];
rz(-0.91667691) q[0];
sx q[0];
rz(0.23707238) q[0];
rz(1.392662) q[1];
sx q[1];
rz(-1.9475513) q[1];
sx q[1];
rz(1.0038092) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8477551) q[0];
sx q[0];
rz(-1.7373573) q[0];
sx q[0];
rz(-2.7012347) q[0];
rz(-pi) q[1];
x q[1];
rz(0.71356524) q[2];
sx q[2];
rz(-1.2110365) q[2];
sx q[2];
rz(-0.52116115) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6663346) q[1];
sx q[1];
rz(-1.043519) q[1];
sx q[1];
rz(1.1858681) q[1];
x q[2];
rz(-0.5587033) q[3];
sx q[3];
rz(-1.1538343) q[3];
sx q[3];
rz(-0.95181634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6134593) q[2];
sx q[2];
rz(-1.5626835) q[2];
sx q[2];
rz(0.51631874) q[2];
rz(-2.3033219) q[3];
sx q[3];
rz(-1.6603989) q[3];
sx q[3];
rz(-0.18394884) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2876005) q[0];
sx q[0];
rz(-0.28165278) q[0];
sx q[0];
rz(-2.2273492) q[0];
rz(2.5573348) q[1];
sx q[1];
rz(-1.7353568) q[1];
sx q[1];
rz(0.90726888) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4205243) q[0];
sx q[0];
rz(-2.4926909) q[0];
sx q[0];
rz(1.6008928) q[0];
rz(-1.7444939) q[2];
sx q[2];
rz(-0.96250421) q[2];
sx q[2];
rz(2.4243958) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.71675352) q[1];
sx q[1];
rz(-2.3100634) q[1];
sx q[1];
rz(-0.6296954) q[1];
rz(-2.1694555) q[3];
sx q[3];
rz(-1.3479509) q[3];
sx q[3];
rz(-3.0226662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.82627901) q[2];
sx q[2];
rz(-1.7118914) q[2];
sx q[2];
rz(-0.86177525) q[2];
rz(1.2365384) q[3];
sx q[3];
rz(-3.0573513) q[3];
sx q[3];
rz(1.2110565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5009907) q[0];
sx q[0];
rz(-2.3589098) q[0];
sx q[0];
rz(2.1114517) q[0];
rz(-2.8252699) q[1];
sx q[1];
rz(-1.373469) q[1];
sx q[1];
rz(-2.8533459) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.41053) q[0];
sx q[0];
rz(-1.1708492) q[0];
sx q[0];
rz(1.7284786) q[0];
x q[1];
rz(0.64149842) q[2];
sx q[2];
rz(-2.3833123) q[2];
sx q[2];
rz(-1.7816133) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2001674) q[1];
sx q[1];
rz(-1.1516478) q[1];
sx q[1];
rz(0.9113542) q[1];
rz(1.796631) q[3];
sx q[3];
rz(-0.91167456) q[3];
sx q[3];
rz(0.82443217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2890702) q[2];
sx q[2];
rz(-1.4129637) q[2];
sx q[2];
rz(0.22656013) q[2];
rz(1.6864927) q[3];
sx q[3];
rz(-0.89613599) q[3];
sx q[3];
rz(-2.4494825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15014547) q[0];
sx q[0];
rz(-1.2757855) q[0];
sx q[0];
rz(2.5164497) q[0];
rz(-1.9236247) q[1];
sx q[1];
rz(-0.70911276) q[1];
sx q[1];
rz(-1.2120754) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21975133) q[0];
sx q[0];
rz(-2.3318687) q[0];
sx q[0];
rz(0.68897665) q[0];
x q[1];
rz(-0.67135129) q[2];
sx q[2];
rz(-1.7749447) q[2];
sx q[2];
rz(2.8189557) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6663078) q[1];
sx q[1];
rz(-1.5130677) q[1];
sx q[1];
rz(2.0121196) q[1];
x q[2];
rz(-0.76833581) q[3];
sx q[3];
rz(-2.7339122) q[3];
sx q[3];
rz(-0.74484315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.43442279) q[2];
sx q[2];
rz(-1.3753128) q[2];
sx q[2];
rz(1.0506857) q[2];
rz(0.55189842) q[3];
sx q[3];
rz(-2.4514908) q[3];
sx q[3];
rz(-1.2845854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5156749) q[0];
sx q[0];
rz(-2.3139625) q[0];
sx q[0];
rz(2.3232842) q[0];
rz(-0.7863518) q[1];
sx q[1];
rz(-0.66023371) q[1];
sx q[1];
rz(0.22088851) q[1];
rz(1.2073866) q[2];
sx q[2];
rz(-1.6061693) q[2];
sx q[2];
rz(2.1314878) q[2];
rz(-0.68610739) q[3];
sx q[3];
rz(-1.6391564) q[3];
sx q[3];
rz(2.2166722) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
