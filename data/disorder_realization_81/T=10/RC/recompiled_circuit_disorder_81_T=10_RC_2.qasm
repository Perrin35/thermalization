OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4047591) q[0];
sx q[0];
rz(-1.7801378) q[0];
sx q[0];
rz(-1.7629495) q[0];
rz(-0.8575851) q[1];
sx q[1];
rz(-1.4839988) q[1];
sx q[1];
rz(-2.690697) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33559088) q[0];
sx q[0];
rz(-1.5561034) q[0];
sx q[0];
rz(3.0528085) q[0];
rz(-pi) q[1];
rz(-3.0300006) q[2];
sx q[2];
rz(-1.0942232) q[2];
sx q[2];
rz(3.1389719) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8494107) q[1];
sx q[1];
rz(-1.0725478) q[1];
sx q[1];
rz(-1.8241747) q[1];
rz(-2.2545635) q[3];
sx q[3];
rz(-1.4966045) q[3];
sx q[3];
rz(2.2825953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1575872) q[2];
sx q[2];
rz(-1.6820587) q[2];
sx q[2];
rz(2.297304) q[2];
rz(2.700581) q[3];
sx q[3];
rz(-0.35566548) q[3];
sx q[3];
rz(0.60602337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59250295) q[0];
sx q[0];
rz(-1.2298158) q[0];
sx q[0];
rz(-2.8785008) q[0];
rz(2.198055) q[1];
sx q[1];
rz(-0.5967921) q[1];
sx q[1];
rz(-1.1862322) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3514254) q[0];
sx q[0];
rz(-1.5895491) q[0];
sx q[0];
rz(3.0856531) q[0];
rz(-pi) q[1];
rz(1.6329174) q[2];
sx q[2];
rz(-1.7698235) q[2];
sx q[2];
rz(-1.8387427) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9895049) q[1];
sx q[1];
rz(-2.4651335) q[1];
sx q[1];
rz(-1.3701887) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.22963345) q[3];
sx q[3];
rz(-2.4335055) q[3];
sx q[3];
rz(-1.1585483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1295604) q[2];
sx q[2];
rz(-1.0027145) q[2];
sx q[2];
rz(-1.1594695) q[2];
rz(2.7705079) q[3];
sx q[3];
rz(-1.6371195) q[3];
sx q[3];
rz(2.8306567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3804669) q[0];
sx q[0];
rz(-2.0115871) q[0];
sx q[0];
rz(-0.80672112) q[0];
rz(0.21356788) q[1];
sx q[1];
rz(-0.49626207) q[1];
sx q[1];
rz(2.321373) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7263111) q[0];
sx q[0];
rz(-1.4584686) q[0];
sx q[0];
rz(0.88322722) q[0];
rz(-pi) q[1];
rz(2.9224612) q[2];
sx q[2];
rz(-2.0543155) q[2];
sx q[2];
rz(2.6205274) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1097475) q[1];
sx q[1];
rz(-2.577563) q[1];
sx q[1];
rz(2.2721223) q[1];
x q[2];
rz(0.1180325) q[3];
sx q[3];
rz(-2.0327838) q[3];
sx q[3];
rz(3.1363917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.31072581) q[2];
sx q[2];
rz(-1.5006289) q[2];
sx q[2];
rz(2.2107928) q[2];
rz(0.15549774) q[3];
sx q[3];
rz(-1.6379387) q[3];
sx q[3];
rz(-2.8500407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.528462) q[0];
sx q[0];
rz(-0.72137946) q[0];
sx q[0];
rz(-0.91127515) q[0];
rz(2.7032734) q[1];
sx q[1];
rz(-1.8194018) q[1];
sx q[1];
rz(1.320425) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94538044) q[0];
sx q[0];
rz(-1.5811265) q[0];
sx q[0];
rz(1.1611847) q[0];
rz(-pi) q[1];
rz(-1.3643866) q[2];
sx q[2];
rz(-2.0136535) q[2];
sx q[2];
rz(-0.70476156) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.98852324) q[1];
sx q[1];
rz(-1.8338025) q[1];
sx q[1];
rz(-0.70446976) q[1];
rz(2.3539691) q[3];
sx q[3];
rz(-1.1155323) q[3];
sx q[3];
rz(0.61592197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0358255) q[2];
sx q[2];
rz(-2.2129009) q[2];
sx q[2];
rz(0.34238112) q[2];
rz(0.17677447) q[3];
sx q[3];
rz(-0.43313679) q[3];
sx q[3];
rz(-2.0006196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
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
rz(0.3115561) q[0];
sx q[0];
rz(-2.4139068) q[0];
sx q[0];
rz(-0.86529055) q[0];
rz(-1.9150437) q[1];
sx q[1];
rz(-0.98926917) q[1];
sx q[1];
rz(1.3006166) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0375992) q[0];
sx q[0];
rz(-2.7866057) q[0];
sx q[0];
rz(-3.0689737) q[0];
rz(-2.956203) q[2];
sx q[2];
rz(-2.7692147) q[2];
sx q[2];
rz(2.7951954) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.865766) q[1];
sx q[1];
rz(-1.2424801) q[1];
sx q[1];
rz(2.5715716) q[1];
rz(1.424765) q[3];
sx q[3];
rz(-1.8304123) q[3];
sx q[3];
rz(2.8433593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0098003) q[2];
sx q[2];
rz(-1.1494145) q[2];
sx q[2];
rz(-0.47719964) q[2];
rz(-0.19208433) q[3];
sx q[3];
rz(-1.6936857) q[3];
sx q[3];
rz(0.93311667) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3451097) q[0];
sx q[0];
rz(-0.61426291) q[0];
sx q[0];
rz(-0.011750301) q[0];
rz(2.5911962) q[1];
sx q[1];
rz(-1.3563211) q[1];
sx q[1];
rz(1.5531497) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1588622) q[0];
sx q[0];
rz(-1.9059062) q[0];
sx q[0];
rz(-1.9431252) q[0];
x q[1];
rz(1.2049098) q[2];
sx q[2];
rz(-1.5503746) q[2];
sx q[2];
rz(2.6423955) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.99332422) q[1];
sx q[1];
rz(-0.71422186) q[1];
sx q[1];
rz(3.0118914) q[1];
x q[2];
rz(-1.6586967) q[3];
sx q[3];
rz(-2.4368736) q[3];
sx q[3];
rz(-0.3375012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.50756303) q[2];
sx q[2];
rz(-0.66036779) q[2];
sx q[2];
rz(-1.2825512) q[2];
rz(1.3698618) q[3];
sx q[3];
rz(-1.3953352) q[3];
sx q[3];
rz(-2.0231358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58105528) q[0];
sx q[0];
rz(-0.16796172) q[0];
sx q[0];
rz(0.67725956) q[0];
rz(0.15180763) q[1];
sx q[1];
rz(-1.3744524) q[1];
sx q[1];
rz(2.1645434) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86977406) q[0];
sx q[0];
rz(-0.5973814) q[0];
sx q[0];
rz(3.0074044) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5708959) q[2];
sx q[2];
rz(-1.4387555) q[2];
sx q[2];
rz(3.0299203) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9338785) q[1];
sx q[1];
rz(-1.2320476) q[1];
sx q[1];
rz(-3.0599041) q[1];
x q[2];
rz(2.3903923) q[3];
sx q[3];
rz(-0.45414543) q[3];
sx q[3];
rz(-2.5129012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3892422) q[2];
sx q[2];
rz(-0.82169473) q[2];
sx q[2];
rz(-1.0127257) q[2];
rz(1.1879454) q[3];
sx q[3];
rz(-1.0725189) q[3];
sx q[3];
rz(-0.48721203) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0174719) q[0];
sx q[0];
rz(-0.033360632) q[0];
sx q[0];
rz(2.4429328) q[0];
rz(-1.1220804) q[1];
sx q[1];
rz(-0.84609234) q[1];
sx q[1];
rz(-1.2493856) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.218924) q[0];
sx q[0];
rz(-1.3601174) q[0];
sx q[0];
rz(-1.523449) q[0];
rz(-pi) q[1];
rz(2.2970389) q[2];
sx q[2];
rz(-0.65659467) q[2];
sx q[2];
rz(-3.055228) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.318995) q[1];
sx q[1];
rz(-1.5486451) q[1];
sx q[1];
rz(-1.5652565) q[1];
rz(-1.1864248) q[3];
sx q[3];
rz(-0.66925183) q[3];
sx q[3];
rz(-0.8347019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.32968783) q[2];
sx q[2];
rz(-0.78616443) q[2];
sx q[2];
rz(-1.9630986) q[2];
rz(1.4568436) q[3];
sx q[3];
rz(-1.0624351) q[3];
sx q[3];
rz(0.38213521) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6417398) q[0];
sx q[0];
rz(-1.3795744) q[0];
sx q[0];
rz(-1.8485803) q[0];
rz(1.4216084) q[1];
sx q[1];
rz(-2.1052108) q[1];
sx q[1];
rz(-0.59757772) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4353838) q[0];
sx q[0];
rz(-1.7058813) q[0];
sx q[0];
rz(-3.1027017) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0963247) q[2];
sx q[2];
rz(-2.0647486) q[2];
sx q[2];
rz(1.8906821) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7087047) q[1];
sx q[1];
rz(-1.27379) q[1];
sx q[1];
rz(-1.2328641) q[1];
x q[2];
rz(-0.58595539) q[3];
sx q[3];
rz(-1.5065985) q[3];
sx q[3];
rz(-0.5639329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.22275816) q[2];
sx q[2];
rz(-1.4514048) q[2];
sx q[2];
rz(-1.9082327) q[2];
rz(-2.2402066) q[3];
sx q[3];
rz(-0.12005761) q[3];
sx q[3];
rz(-1.6433158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0937061) q[0];
sx q[0];
rz(-2.369635) q[0];
sx q[0];
rz(-3.1179324) q[0];
rz(-0.95611447) q[1];
sx q[1];
rz(-1.3095983) q[1];
sx q[1];
rz(2.4694209) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2490847) q[0];
sx q[0];
rz(-1.5605643) q[0];
sx q[0];
rz(0.0052878629) q[0];
rz(1.5615084) q[2];
sx q[2];
rz(-1.0592959) q[2];
sx q[2];
rz(-1.4052504) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.086239554) q[1];
sx q[1];
rz(-1.2330016) q[1];
sx q[1];
rz(0.71725459) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4168596) q[3];
sx q[3];
rz(-1.142475) q[3];
sx q[3];
rz(-1.1101013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.51222926) q[2];
sx q[2];
rz(-1.9146634) q[2];
sx q[2];
rz(0.36995861) q[2];
rz(-1.5036748) q[3];
sx q[3];
rz(-2.2556997) q[3];
sx q[3];
rz(1.9406208) q[3];
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
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5794012) q[0];
sx q[0];
rz(-0.36407064) q[0];
sx q[0];
rz(-1.9343485) q[0];
rz(2.4178986) q[1];
sx q[1];
rz(-2.1543398) q[1];
sx q[1];
rz(2.2347246) q[1];
rz(1.9359246) q[2];
sx q[2];
rz(-2.7177313) q[2];
sx q[2];
rz(2.360366) q[2];
rz(-1.9772114) q[3];
sx q[3];
rz(-1.2953399) q[3];
sx q[3];
rz(-3.0084707) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];