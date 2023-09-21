OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4560661) q[0];
sx q[0];
rz(-0.38903061) q[0];
sx q[0];
rz(2.2580137) q[0];
rz(3.1318624) q[1];
sx q[1];
rz(-1.6844123) q[1];
sx q[1];
rz(-1.943346) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53579522) q[0];
sx q[0];
rz(-1.747316) q[0];
sx q[0];
rz(-1.4730886) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.892957) q[2];
sx q[2];
rz(-0.82167168) q[2];
sx q[2];
rz(-1.2920213) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8625496) q[1];
sx q[1];
rz(-1.057813) q[1];
sx q[1];
rz(-1.1898477) q[1];
x q[2];
rz(-1.8331068) q[3];
sx q[3];
rz(-2.9648818) q[3];
sx q[3];
rz(-1.3486598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9464232) q[2];
sx q[2];
rz(-2.158458) q[2];
sx q[2];
rz(-2.9602489) q[2];
rz(0.26120734) q[3];
sx q[3];
rz(-1.8758592) q[3];
sx q[3];
rz(2.3852824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8392035) q[0];
sx q[0];
rz(-2.8518682) q[0];
sx q[0];
rz(2.7547577) q[0];
rz(0.50239262) q[1];
sx q[1];
rz(-0.97351176) q[1];
sx q[1];
rz(-1.5418672) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4138448) q[0];
sx q[0];
rz(-1.9808597) q[0];
sx q[0];
rz(-0.64322612) q[0];
rz(-pi) q[1];
rz(1.1142271) q[2];
sx q[2];
rz(-1.7184988) q[2];
sx q[2];
rz(-0.6870803) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9215556) q[1];
sx q[1];
rz(-1.6065292) q[1];
sx q[1];
rz(-1.5807371) q[1];
rz(-2.3719671) q[3];
sx q[3];
rz(-2.5105021) q[3];
sx q[3];
rz(2.1345994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5388422) q[2];
sx q[2];
rz(-0.87783146) q[2];
sx q[2];
rz(-1.7269469) q[2];
rz(-2.291262) q[3];
sx q[3];
rz(-2.7089705) q[3];
sx q[3];
rz(-1.6833646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-0.26281115) q[0];
sx q[0];
rz(-1.5314064) q[0];
sx q[0];
rz(-2.5313654) q[0];
rz(1.8521076) q[1];
sx q[1];
rz(-0.97924966) q[1];
sx q[1];
rz(0.99197018) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2881644) q[0];
sx q[0];
rz(-0.58840226) q[0];
sx q[0];
rz(1.9073652) q[0];
rz(-pi) q[1];
rz(2.8875071) q[2];
sx q[2];
rz(-2.457329) q[2];
sx q[2];
rz(2.4917045) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.75533463) q[1];
sx q[1];
rz(-1.660368) q[1];
sx q[1];
rz(0.44134015) q[1];
rz(0.45735995) q[3];
sx q[3];
rz(-1.7914346) q[3];
sx q[3];
rz(-1.8822576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.758574) q[2];
sx q[2];
rz(-2.980361) q[2];
sx q[2];
rz(-2.8733011) q[2];
rz(0.39408436) q[3];
sx q[3];
rz(-1.2309309) q[3];
sx q[3];
rz(2.9530853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2878993) q[0];
sx q[0];
rz(-2.6278966) q[0];
sx q[0];
rz(-2.7365141) q[0];
rz(-2.4515117) q[1];
sx q[1];
rz(-1.9837374) q[1];
sx q[1];
rz(-2.4437723) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9245802) q[0];
sx q[0];
rz(-0.092886535) q[0];
sx q[0];
rz(-1.8266982) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.54450808) q[2];
sx q[2];
rz(-1.240057) q[2];
sx q[2];
rz(-2.7059908) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.3623558) q[1];
sx q[1];
rz(-2.6168565) q[1];
sx q[1];
rz(1.9441821) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9911733) q[3];
sx q[3];
rz(-1.8432518) q[3];
sx q[3];
rz(3.0416995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4219389) q[2];
sx q[2];
rz(-1.3976588) q[2];
sx q[2];
rz(-1.5092124) q[2];
rz(-0.40431067) q[3];
sx q[3];
rz(-2.4590838) q[3];
sx q[3];
rz(-1.6633165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25800911) q[0];
sx q[0];
rz(-1.3556577) q[0];
sx q[0];
rz(1.0255381) q[0];
rz(2.569596) q[1];
sx q[1];
rz(-2.0472186) q[1];
sx q[1];
rz(-2.5122723) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.660491) q[0];
sx q[0];
rz(-1.2718624) q[0];
sx q[0];
rz(-0.42994182) q[0];
x q[1];
rz(2.8867678) q[2];
sx q[2];
rz(-2.2307768) q[2];
sx q[2];
rz(-3.0072336) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3931261) q[1];
sx q[1];
rz(-1.6854291) q[1];
sx q[1];
rz(0.51490358) q[1];
rz(-pi) q[2];
rz(2.878177) q[3];
sx q[3];
rz(-0.98589555) q[3];
sx q[3];
rz(-1.5358621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.59262529) q[2];
sx q[2];
rz(-2.7719345) q[2];
sx q[2];
rz(0.34234753) q[2];
rz(-1.7957548) q[3];
sx q[3];
rz(-1.6941518) q[3];
sx q[3];
rz(2.9455744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.489007) q[0];
sx q[0];
rz(-1.2359897) q[0];
sx q[0];
rz(-0.18519369) q[0];
rz(-1.406503) q[1];
sx q[1];
rz(-1.0909189) q[1];
sx q[1];
rz(-1.7746183) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6706086) q[0];
sx q[0];
rz(-2.1489722) q[0];
sx q[0];
rz(1.3060119) q[0];
x q[1];
rz(1.9063437) q[2];
sx q[2];
rz(-2.3830072) q[2];
sx q[2];
rz(2.2461265) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5806611) q[1];
sx q[1];
rz(-1.7226189) q[1];
sx q[1];
rz(3.0488357) q[1];
rz(2.5227138) q[3];
sx q[3];
rz(-0.66580171) q[3];
sx q[3];
rz(2.562079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.24484816) q[2];
sx q[2];
rz(-2.5881793) q[2];
sx q[2];
rz(0.883376) q[2];
rz(1.7287438) q[3];
sx q[3];
rz(-2.4491375) q[3];
sx q[3];
rz(0.81378716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8900523) q[0];
sx q[0];
rz(-3.0391356) q[0];
sx q[0];
rz(-1.2782156) q[0];
rz(-0.037840769) q[1];
sx q[1];
rz(-0.81532878) q[1];
sx q[1];
rz(-1.3758804) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85526953) q[0];
sx q[0];
rz(-2.3648242) q[0];
sx q[0];
rz(1.1458678) q[0];
rz(0.059869754) q[2];
sx q[2];
rz(-2.0436358) q[2];
sx q[2];
rz(-1.5945895) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3657276) q[1];
sx q[1];
rz(-2.9850246) q[1];
sx q[1];
rz(2.4003029) q[1];
rz(-pi) q[2];
rz(-1.4010299) q[3];
sx q[3];
rz(-2.8068636) q[3];
sx q[3];
rz(1.7498121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4574796) q[2];
sx q[2];
rz(-0.71762466) q[2];
sx q[2];
rz(-2.4105371) q[2];
rz(-3.030792) q[3];
sx q[3];
rz(-1.585107) q[3];
sx q[3];
rz(0.69537648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5450127) q[0];
sx q[0];
rz(-2.2667363) q[0];
sx q[0];
rz(-0.73356432) q[0];
rz(-0.60797524) q[1];
sx q[1];
rz(-1.9476451) q[1];
sx q[1];
rz(-0.2342934) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60093555) q[0];
sx q[0];
rz(-2.0458851) q[0];
sx q[0];
rz(2.8500772) q[0];
x q[1];
rz(2.2393164) q[2];
sx q[2];
rz(-0.96792816) q[2];
sx q[2];
rz(2.5831646) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.81208166) q[1];
sx q[1];
rz(-2.3595464) q[1];
sx q[1];
rz(1.300632) q[1];
rz(2.8041744) q[3];
sx q[3];
rz(-1.0113582) q[3];
sx q[3];
rz(2.9448201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0155448) q[2];
sx q[2];
rz(-1.3092821) q[2];
sx q[2];
rz(-1.4253433) q[2];
rz(1.6783293) q[3];
sx q[3];
rz(-2.3571456) q[3];
sx q[3];
rz(-1.6459758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5995246) q[0];
sx q[0];
rz(-2.8064089) q[0];
sx q[0];
rz(1.9482127) q[0];
rz(1.880973) q[1];
sx q[1];
rz(-1.7766989) q[1];
sx q[1];
rz(2.1967922) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1856857) q[0];
sx q[0];
rz(-0.80230306) q[0];
sx q[0];
rz(0.54922744) q[0];
rz(-pi) q[1];
x q[1];
rz(0.94061942) q[2];
sx q[2];
rz(-0.47626469) q[2];
sx q[2];
rz(0.95048743) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8972842) q[1];
sx q[1];
rz(-2.8123887) q[1];
sx q[1];
rz(0.55397482) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9267843) q[3];
sx q[3];
rz(-1.3857406) q[3];
sx q[3];
rz(-1.2026279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.21645674) q[2];
sx q[2];
rz(-1.0711121) q[2];
sx q[2];
rz(-0.28820583) q[2];
rz(2.6618585) q[3];
sx q[3];
rz(-2.0917442) q[3];
sx q[3];
rz(-1.6335999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11186803) q[0];
sx q[0];
rz(-0.27619633) q[0];
sx q[0];
rz(0.9129886) q[0];
rz(2.7669725) q[1];
sx q[1];
rz(-1.4034142) q[1];
sx q[1];
rz(-0.8909117) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64542949) q[0];
sx q[0];
rz(-1.9053639) q[0];
sx q[0];
rz(-1.3076925) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.95821361) q[2];
sx q[2];
rz(-1.805086) q[2];
sx q[2];
rz(2.2189552) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.76162321) q[1];
sx q[1];
rz(-0.39561158) q[1];
sx q[1];
rz(-0.53669866) q[1];
rz(0.18817801) q[3];
sx q[3];
rz(-1.6450226) q[3];
sx q[3];
rz(-1.0159514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9051819) q[2];
sx q[2];
rz(-2.2832401) q[2];
sx q[2];
rz(-2.6399844) q[2];
rz(-1.8524528) q[3];
sx q[3];
rz(-1.6882378) q[3];
sx q[3];
rz(2.6954209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5065153) q[0];
sx q[0];
rz(-1.4415393) q[0];
sx q[0];
rz(-2.517979) q[0];
rz(-2.0093909) q[1];
sx q[1];
rz(-0.75695801) q[1];
sx q[1];
rz(-3.0523041) q[1];
rz(2.8749309) q[2];
sx q[2];
rz(-1.4229792) q[2];
sx q[2];
rz(-0.40001043) q[2];
rz(3.0922079) q[3];
sx q[3];
rz(-2.1168843) q[3];
sx q[3];
rz(-2.4612853) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];