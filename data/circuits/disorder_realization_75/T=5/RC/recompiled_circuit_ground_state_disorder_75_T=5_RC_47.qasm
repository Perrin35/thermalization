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
rz(6.5919696) q[0];
sx q[0];
rz(6.0433521) q[0];
rz(2.580515) q[1];
sx q[1];
rz(-0.74220389) q[1];
sx q[1];
rz(1.5129369) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1334907) q[0];
sx q[0];
rz(-1.9372131) q[0];
sx q[0];
rz(-3.0045463) q[0];
rz(1.0983019) q[2];
sx q[2];
rz(-1.8046981) q[2];
sx q[2];
rz(-2.5462674) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4862772) q[1];
sx q[1];
rz(-1.4239171) q[1];
sx q[1];
rz(2.6492061) q[1];
x q[2];
rz(-1.2414819) q[3];
sx q[3];
rz(-0.17767492) q[3];
sx q[3];
rz(-0.8575646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.98770398) q[2];
sx q[2];
rz(-2.14553) q[2];
sx q[2];
rz(1.055701) q[2];
rz(2.740247) q[3];
sx q[3];
rz(-1.6258806) q[3];
sx q[3];
rz(1.5041941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6317247) q[0];
sx q[0];
rz(-0.26917502) q[0];
sx q[0];
rz(1.4058231) q[0];
rz(-2.3954605) q[1];
sx q[1];
rz(-1.965062) q[1];
sx q[1];
rz(3.0768118) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5279605) q[0];
sx q[0];
rz(-1.1224696) q[0];
sx q[0];
rz(-1.1290068) q[0];
rz(-pi) q[1];
x q[1];
rz(0.20938907) q[2];
sx q[2];
rz(-2.4509666) q[2];
sx q[2];
rz(-2.2356981) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3105615) q[1];
sx q[1];
rz(-1.4906724) q[1];
sx q[1];
rz(-1.4270272) q[1];
x q[2];
rz(3.1094743) q[3];
sx q[3];
rz(-1.7815551) q[3];
sx q[3];
rz(-3.075066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.88227162) q[2];
sx q[2];
rz(-2.6210625) q[2];
sx q[2];
rz(3.1003013) q[2];
rz(1.1193554) q[3];
sx q[3];
rz(-1.7740403) q[3];
sx q[3];
rz(-2.0004499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3248046) q[0];
sx q[0];
rz(-0.4275221) q[0];
sx q[0];
rz(-2.4422755) q[0];
rz(-0.62272561) q[1];
sx q[1];
rz(-0.8131578) q[1];
sx q[1];
rz(0.25845382) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1656483) q[0];
sx q[0];
rz(-1.5585525) q[0];
sx q[0];
rz(-0.6338288) q[0];
rz(-pi) q[1];
rz(1.4213561) q[2];
sx q[2];
rz(-0.65776134) q[2];
sx q[2];
rz(-0.8873111) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.041641673) q[1];
sx q[1];
rz(-1.1219684) q[1];
sx q[1];
rz(-1.5931409) q[1];
rz(3.1051226) q[3];
sx q[3];
rz(-2.0005282) q[3];
sx q[3];
rz(-0.018485395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2059325) q[2];
sx q[2];
rz(-1.460133) q[2];
sx q[2];
rz(-0.70651954) q[2];
rz(-0.62890729) q[3];
sx q[3];
rz(-2.3157412) q[3];
sx q[3];
rz(-1.1227054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.0141456) q[0];
sx q[0];
rz(-1.4268459) q[0];
sx q[0];
rz(0.97417796) q[0];
rz(1.4840508) q[1];
sx q[1];
rz(-1.0933135) q[1];
sx q[1];
rz(-2.1407703) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7580622) q[0];
sx q[0];
rz(-0.62072004) q[0];
sx q[0];
rz(-1.8412983) q[0];
rz(-pi) q[1];
rz(-0.8128667) q[2];
sx q[2];
rz(-2.1846601) q[2];
sx q[2];
rz(-0.46068207) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.19480579) q[1];
sx q[1];
rz(-2.0344224) q[1];
sx q[1];
rz(2.5083816) q[1];
rz(-pi) q[2];
rz(-0.14814143) q[3];
sx q[3];
rz(-2.8175948) q[3];
sx q[3];
rz(1.9213541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4546844) q[2];
sx q[2];
rz(-1.3409706) q[2];
sx q[2];
rz(2.3243813) q[2];
rz(0.59598437) q[3];
sx q[3];
rz(-1.7242566) q[3];
sx q[3];
rz(-1.7433085) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78438321) q[0];
sx q[0];
rz(-1.4056118) q[0];
sx q[0];
rz(-1.0719517) q[0];
rz(2.0042073) q[1];
sx q[1];
rz(-1.5147361) q[1];
sx q[1];
rz(2.9683108) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77861518) q[0];
sx q[0];
rz(-2.5234875) q[0];
sx q[0];
rz(-0.72160665) q[0];
rz(-2.7543147) q[2];
sx q[2];
rz(-2.0301295) q[2];
sx q[2];
rz(-1.1371374) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2318423) q[1];
sx q[1];
rz(-0.76059231) q[1];
sx q[1];
rz(-1.4666345) q[1];
rz(-0.51637465) q[3];
sx q[3];
rz(-0.18680113) q[3];
sx q[3];
rz(1.1774225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.6413573) q[2];
sx q[2];
rz(-2.6884029) q[2];
sx q[2];
rz(0.87453169) q[2];
rz(-2.4747961) q[3];
sx q[3];
rz(-1.6801445) q[3];
sx q[3];
rz(2.3165406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85103971) q[0];
sx q[0];
rz(-0.14851004) q[0];
sx q[0];
rz(0.17459757) q[0];
rz(0.54706508) q[1];
sx q[1];
rz(-0.51536307) q[1];
sx q[1];
rz(1.863716) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7449943) q[0];
sx q[0];
rz(-1.5088668) q[0];
sx q[0];
rz(-1.2838191) q[0];
x q[1];
rz(2.6447884) q[2];
sx q[2];
rz(-0.81648177) q[2];
sx q[2];
rz(1.7330488) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7050347) q[1];
sx q[1];
rz(-2.3113657) q[1];
sx q[1];
rz(-2.6565246) q[1];
x q[2];
rz(2.296769) q[3];
sx q[3];
rz(-0.81171821) q[3];
sx q[3];
rz(2.7561848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7586907) q[2];
sx q[2];
rz(-0.26472696) q[2];
sx q[2];
rz(-2.7721789) q[2];
rz(2.3139125) q[3];
sx q[3];
rz(-1.8281432) q[3];
sx q[3];
rz(-2.9874492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5740042) q[0];
sx q[0];
rz(-2.2249157) q[0];
sx q[0];
rz(2.9045203) q[0];
rz(1.7489307) q[1];
sx q[1];
rz(-1.9475513) q[1];
sx q[1];
rz(-1.0038092) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0802409) q[0];
sx q[0];
rz(-2.6727242) q[0];
sx q[0];
rz(0.37567015) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1091484) q[2];
sx q[2];
rz(-2.2301815) q[2];
sx q[2];
rz(-2.387799) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4752581) q[1];
sx q[1];
rz(-1.043519) q[1];
sx q[1];
rz(-1.1858681) q[1];
rz(-pi) q[2];
x q[2];
rz(0.5587033) q[3];
sx q[3];
rz(-1.1538343) q[3];
sx q[3];
rz(-2.1897763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.5281333) q[2];
sx q[2];
rz(-1.5626835) q[2];
sx q[2];
rz(-0.51631874) q[2];
rz(-0.83827072) q[3];
sx q[3];
rz(-1.6603989) q[3];
sx q[3];
rz(-2.9576438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2876005) q[0];
sx q[0];
rz(-0.28165278) q[0];
sx q[0];
rz(0.91424346) q[0];
rz(0.5842579) q[1];
sx q[1];
rz(-1.4062358) q[1];
sx q[1];
rz(0.90726888) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4205243) q[0];
sx q[0];
rz(-0.64890175) q[0];
sx q[0];
rz(1.5406999) q[0];
x q[1];
rz(1.3970988) q[2];
sx q[2];
rz(-0.96250421) q[2];
sx q[2];
rz(2.4243958) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.71675352) q[1];
sx q[1];
rz(-0.83152926) q[1];
sx q[1];
rz(2.5118973) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9531293) q[3];
sx q[3];
rz(-2.5075932) q[3];
sx q[3];
rz(1.3764189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.82627901) q[2];
sx q[2];
rz(-1.7118914) q[2];
sx q[2];
rz(0.86177525) q[2];
rz(-1.9050542) q[3];
sx q[3];
rz(-0.084241353) q[3];
sx q[3];
rz(1.9305362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5009907) q[0];
sx q[0];
rz(-0.7826829) q[0];
sx q[0];
rz(-2.1114517) q[0];
rz(0.31632272) q[1];
sx q[1];
rz(-1.7681237) q[1];
sx q[1];
rz(-0.28824678) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0228393) q[0];
sx q[0];
rz(-2.7132479) q[0];
sx q[0];
rz(0.35566766) q[0];
rz(1.0551528) q[2];
sx q[2];
rz(-2.1543157) q[2];
sx q[2];
rz(2.1595364) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.32375755) q[1];
sx q[1];
rz(-2.164748) q[1];
sx q[1];
rz(-0.51333921) q[1];
rz(-0.28137286) q[3];
sx q[3];
rz(-0.69123879) q[3];
sx q[3];
rz(0.465525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.85252243) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9914472) q[0];
sx q[0];
rz(-1.8658072) q[0];
sx q[0];
rz(0.62514296) q[0];
rz(-1.9236247) q[1];
sx q[1];
rz(-2.4324799) q[1];
sx q[1];
rz(-1.9295173) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9218413) q[0];
sx q[0];
rz(-2.3318687) q[0];
sx q[0];
rz(0.68897665) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3122968) q[2];
sx q[2];
rz(-0.91583672) q[2];
sx q[2];
rz(1.0884681) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6663078) q[1];
sx q[1];
rz(-1.5130677) q[1];
sx q[1];
rz(-2.0121196) q[1];
rz(-pi) q[2];
rz(-1.2792222) q[3];
sx q[3];
rz(-1.2816888) q[3];
sx q[3];
rz(1.5559199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7071699) q[2];
sx q[2];
rz(-1.3753128) q[2];
sx q[2];
rz(-2.0909069) q[2];
rz(-0.55189842) q[3];
sx q[3];
rz(-0.69010186) q[3];
sx q[3];
rz(-1.2845854) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62591775) q[0];
sx q[0];
rz(-0.82763012) q[0];
sx q[0];
rz(-0.81830842) q[0];
rz(2.3552409) q[1];
sx q[1];
rz(-0.66023371) q[1];
sx q[1];
rz(0.22088851) q[1];
rz(-1.4715696) q[2];
sx q[2];
rz(-0.36505112) q[2];
sx q[2];
rz(0.65341841) q[2];
rz(0.68610739) q[3];
sx q[3];
rz(-1.5024363) q[3];
sx q[3];
rz(-0.92492044) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
