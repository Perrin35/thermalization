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
rz(1.6588563) q[0];
sx q[0];
rz(-0.98134494) q[0];
sx q[0];
rz(-2.0438097) q[0];
rz(2.3535347) q[1];
sx q[1];
rz(-1.0507974) q[1];
sx q[1];
rz(0.0069590574) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46411447) q[0];
sx q[0];
rz(-0.65271806) q[0];
sx q[0];
rz(-1.4783036) q[0];
rz(-pi) q[1];
rz(-1.263515) q[2];
sx q[2];
rz(-2.4715965) q[2];
sx q[2];
rz(-0.43596632) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4930058) q[1];
sx q[1];
rz(-1.8064335) q[1];
sx q[1];
rz(-0.034502397) q[1];
rz(-2.801729) q[3];
sx q[3];
rz(-1.6714665) q[3];
sx q[3];
rz(1.0235746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2385345) q[2];
sx q[2];
rz(-1.7944585) q[2];
sx q[2];
rz(-1.6713589) q[2];
rz(-3.0155731) q[3];
sx q[3];
rz(-2.6991548) q[3];
sx q[3];
rz(0.68388763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0266492) q[0];
sx q[0];
rz(-2.6925955) q[0];
sx q[0];
rz(2.5057416) q[0];
rz(0.77330971) q[1];
sx q[1];
rz(-2.4644303) q[1];
sx q[1];
rz(-1.4453452) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6757386) q[0];
sx q[0];
rz(-0.93198085) q[0];
sx q[0];
rz(0.56446989) q[0];
rz(-pi) q[1];
x q[1];
rz(1.097963) q[2];
sx q[2];
rz(-2.2961628) q[2];
sx q[2];
rz(-2.5918317) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6856001) q[1];
sx q[1];
rz(-1.4135135) q[1];
sx q[1];
rz(-2.434751) q[1];
x q[2];
rz(-0.068447114) q[3];
sx q[3];
rz(-1.4395778) q[3];
sx q[3];
rz(0.4438627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9952952) q[2];
sx q[2];
rz(-1.3714906) q[2];
sx q[2];
rz(3.0038317) q[2];
rz(2.6855101) q[3];
sx q[3];
rz(-0.59499732) q[3];
sx q[3];
rz(-0.47541398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59738961) q[0];
sx q[0];
rz(-1.5837639) q[0];
sx q[0];
rz(2.5624045) q[0];
rz(0.10313615) q[1];
sx q[1];
rz(-1.8201647) q[1];
sx q[1];
rz(2.3613222) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4057377) q[0];
sx q[0];
rz(-0.097022382) q[0];
sx q[0];
rz(2.9461622) q[0];
rz(2.9452376) q[2];
sx q[2];
rz(-0.69075023) q[2];
sx q[2];
rz(0.029411246) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.48371485) q[1];
sx q[1];
rz(-2.6838852) q[1];
sx q[1];
rz(-1.0545516) q[1];
rz(-pi) q[2];
rz(0.87522755) q[3];
sx q[3];
rz(-1.7415294) q[3];
sx q[3];
rz(3.1306992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.66037336) q[2];
sx q[2];
rz(-1.6796835) q[2];
sx q[2];
rz(3.050728) q[2];
rz(-1.4800492) q[3];
sx q[3];
rz(-1.3250947) q[3];
sx q[3];
rz(-0.81502325) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12260967) q[0];
sx q[0];
rz(-0.18773395) q[0];
sx q[0];
rz(2.6222099) q[0];
rz(-1.9620365) q[1];
sx q[1];
rz(-0.68125454) q[1];
sx q[1];
rz(0.048351668) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40019401) q[0];
sx q[0];
rz(-2.1693008) q[0];
sx q[0];
rz(2.5116176) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7555356) q[2];
sx q[2];
rz(-2.0252844) q[2];
sx q[2];
rz(2.7234361) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.4614233) q[1];
sx q[1];
rz(-1.5986414) q[1];
sx q[1];
rz(-0.48770406) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3465041) q[3];
sx q[3];
rz(-1.7582736) q[3];
sx q[3];
rz(-0.80611967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4590596) q[2];
sx q[2];
rz(-0.30632633) q[2];
sx q[2];
rz(-1.691386) q[2];
rz(2.0356483) q[3];
sx q[3];
rz(-1.148843) q[3];
sx q[3];
rz(-0.86627427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76097101) q[0];
sx q[0];
rz(-2.3411317) q[0];
sx q[0];
rz(2.3642484) q[0];
rz(-0.49579534) q[1];
sx q[1];
rz(-0.40758857) q[1];
sx q[1];
rz(-0.19283238) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7161761) q[0];
sx q[0];
rz(-1.8766073) q[0];
sx q[0];
rz(-0.56434631) q[0];
x q[1];
rz(-1.8106789) q[2];
sx q[2];
rz(-2.3238306) q[2];
sx q[2];
rz(0.32599923) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.3461756) q[1];
sx q[1];
rz(-2.8320791) q[1];
sx q[1];
rz(-1.2796106) q[1];
x q[2];
rz(-2.9787872) q[3];
sx q[3];
rz(-1.1737524) q[3];
sx q[3];
rz(2.5474078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5567646) q[2];
sx q[2];
rz(-2.3636221) q[2];
sx q[2];
rz(-2.3053816) q[2];
rz(-0.95544514) q[3];
sx q[3];
rz(-1.7183869) q[3];
sx q[3];
rz(-2.6098765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29109508) q[0];
sx q[0];
rz(-1.2245155) q[0];
sx q[0];
rz(-3.1078872) q[0];
rz(2.2233502) q[1];
sx q[1];
rz(-2.4372209) q[1];
sx q[1];
rz(-2.4423626) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79141533) q[0];
sx q[0];
rz(-0.79424131) q[0];
sx q[0];
rz(2.3167531) q[0];
rz(-pi) q[1];
rz(2.8644833) q[2];
sx q[2];
rz(-2.94706) q[2];
sx q[2];
rz(-2.3152318) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0306176) q[1];
sx q[1];
rz(-0.4781107) q[1];
sx q[1];
rz(0.39076372) q[1];
rz(-0.27856234) q[3];
sx q[3];
rz(-1.2435561) q[3];
sx q[3];
rz(0.95765169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0596727) q[2];
sx q[2];
rz(-0.6531738) q[2];
sx q[2];
rz(-0.095452249) q[2];
rz(2.0898315) q[3];
sx q[3];
rz(-1.2442518) q[3];
sx q[3];
rz(0.89422798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28602257) q[0];
sx q[0];
rz(-2.6595071) q[0];
sx q[0];
rz(1.1676189) q[0];
rz(0.35320148) q[1];
sx q[1];
rz(-0.51274061) q[1];
sx q[1];
rz(-2.8599427) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62171157) q[0];
sx q[0];
rz(-2.6586652) q[0];
sx q[0];
rz(-1.3819225) q[0];
x q[1];
rz(0.1408351) q[2];
sx q[2];
rz(-0.65208921) q[2];
sx q[2];
rz(1.2229133) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0822009) q[1];
sx q[1];
rz(-1.4363465) q[1];
sx q[1];
rz(-1.7198635) q[1];
x q[2];
rz(2.6572833) q[3];
sx q[3];
rz(-2.0748027) q[3];
sx q[3];
rz(-0.33226099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.99792751) q[2];
sx q[2];
rz(-1.9196271) q[2];
sx q[2];
rz(-1.3227051) q[2];
rz(-1.6392684) q[3];
sx q[3];
rz(-2.0570677) q[3];
sx q[3];
rz(0.098202078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99763501) q[0];
sx q[0];
rz(-1.8956381) q[0];
sx q[0];
rz(2.8676046) q[0];
rz(0.59448376) q[1];
sx q[1];
rz(-1.4130054) q[1];
sx q[1];
rz(0.53057539) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0218567) q[0];
sx q[0];
rz(-1.6101735) q[0];
sx q[0];
rz(-1.5669797) q[0];
x q[1];
rz(0.55270393) q[2];
sx q[2];
rz(-1.6356907) q[2];
sx q[2];
rz(-0.12285025) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6484687) q[1];
sx q[1];
rz(-1.1965355) q[1];
sx q[1];
rz(-1.3918274) q[1];
x q[2];
rz(-2.9782615) q[3];
sx q[3];
rz(-2.2769351) q[3];
sx q[3];
rz(-0.85757366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.24255594) q[2];
sx q[2];
rz(-1.9344923) q[2];
sx q[2];
rz(2.8846018) q[2];
rz(-1.9422003) q[3];
sx q[3];
rz(-1.2446087) q[3];
sx q[3];
rz(-1.5132343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5892107) q[0];
sx q[0];
rz(-0.95443812) q[0];
sx q[0];
rz(0.5300262) q[0];
rz(1.5026622) q[1];
sx q[1];
rz(-1.9866147) q[1];
sx q[1];
rz(-0.93856215) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28806557) q[0];
sx q[0];
rz(-2.1479283) q[0];
sx q[0];
rz(-2.1398628) q[0];
x q[1];
rz(0.1436032) q[2];
sx q[2];
rz(-1.4521027) q[2];
sx q[2];
rz(-2.375756) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5903428) q[1];
sx q[1];
rz(-1.1444725) q[1];
sx q[1];
rz(0.23386441) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1562002) q[3];
sx q[3];
rz(-1.4568351) q[3];
sx q[3];
rz(1.0572019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3878801) q[2];
sx q[2];
rz(-0.91999274) q[2];
sx q[2];
rz(-2.06125) q[2];
rz(-0.8148109) q[3];
sx q[3];
rz(-1.1956513) q[3];
sx q[3];
rz(1.4174392) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.998488) q[0];
sx q[0];
rz(-1.657722) q[0];
sx q[0];
rz(2.0527573) q[0];
rz(0.067527436) q[1];
sx q[1];
rz(-1.8636401) q[1];
sx q[1];
rz(-2.3823104) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6367594) q[0];
sx q[0];
rz(-2.8593316) q[0];
sx q[0];
rz(-2.0105848) q[0];
x q[1];
rz(2.6443321) q[2];
sx q[2];
rz(-2.8858375) q[2];
sx q[2];
rz(-1.349468) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2226505) q[1];
sx q[1];
rz(-0.83178751) q[1];
sx q[1];
rz(1.8636501) q[1];
rz(-pi) q[2];
rz(1.3410232) q[3];
sx q[3];
rz(-1.9073351) q[3];
sx q[3];
rz(2.0391109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.33409432) q[2];
sx q[2];
rz(-1.0958025) q[2];
sx q[2];
rz(-2.5661772) q[2];
rz(2.5790162) q[3];
sx q[3];
rz(-2.176599) q[3];
sx q[3];
rz(-1.7687198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.060487735) q[0];
sx q[0];
rz(-1.2860379) q[0];
sx q[0];
rz(2.2338569) q[0];
rz(2.0545215) q[1];
sx q[1];
rz(-2.0094951) q[1];
sx q[1];
rz(3.1241945) q[1];
rz(0.74093735) q[2];
sx q[2];
rz(-2.9334684) q[2];
sx q[2];
rz(2.8192782) q[2];
rz(-0.51856507) q[3];
sx q[3];
rz(-0.66558481) q[3];
sx q[3];
rz(-0.44105327) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
