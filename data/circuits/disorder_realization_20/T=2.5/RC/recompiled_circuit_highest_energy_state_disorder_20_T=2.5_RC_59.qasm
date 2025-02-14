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
rz(2.4169154) q[0];
sx q[0];
rz(-1.3114572) q[0];
sx q[0];
rz(0.80479446) q[0];
rz(0.59504396) q[1];
sx q[1];
rz(-2.9437328) q[1];
sx q[1];
rz(1.2208389) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0752995) q[0];
sx q[0];
rz(-1.6907215) q[0];
sx q[0];
rz(-1.6060702) q[0];
rz(-0.27313613) q[2];
sx q[2];
rz(-0.47800666) q[2];
sx q[2];
rz(-0.82551685) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2455622) q[1];
sx q[1];
rz(-0.82791019) q[1];
sx q[1];
rz(-1.5106535) q[1];
rz(-pi) q[2];
rz(-0.63296707) q[3];
sx q[3];
rz(-1.4598914) q[3];
sx q[3];
rz(-1.9680915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4294942) q[2];
sx q[2];
rz(-0.50769371) q[2];
sx q[2];
rz(-1.8053767) q[2];
rz(-1.3974238) q[3];
sx q[3];
rz(-1.6349399) q[3];
sx q[3];
rz(-1.5558745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76032388) q[0];
sx q[0];
rz(-1.6345432) q[0];
sx q[0];
rz(0.67833483) q[0];
rz(-2.5256269) q[1];
sx q[1];
rz(-1.0266285) q[1];
sx q[1];
rz(-2.9982627) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.029328922) q[0];
sx q[0];
rz(-1.8929795) q[0];
sx q[0];
rz(-2.9759745) q[0];
rz(-pi) q[1];
x q[1];
rz(2.289166) q[2];
sx q[2];
rz(-1.4839236) q[2];
sx q[2];
rz(1.014682) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.60494938) q[1];
sx q[1];
rz(-0.50770137) q[1];
sx q[1];
rz(-3.0439218) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.67875989) q[3];
sx q[3];
rz(-2.6142411) q[3];
sx q[3];
rz(0.55298978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4914638) q[2];
sx q[2];
rz(-0.77892196) q[2];
sx q[2];
rz(-1.9219363) q[2];
rz(0.31798142) q[3];
sx q[3];
rz(-2.3173083) q[3];
sx q[3];
rz(-0.73941755) q[3];
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
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3984482) q[0];
sx q[0];
rz(-0.92107451) q[0];
sx q[0];
rz(0.63968023) q[0];
rz(1.3321053) q[1];
sx q[1];
rz(-2.2001241) q[1];
sx q[1];
rz(2.926362) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4217241) q[0];
sx q[0];
rz(-2.6009237) q[0];
sx q[0];
rz(2.7175277) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.136342) q[2];
sx q[2];
rz(-1.286478) q[2];
sx q[2];
rz(2.9470987) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.135268) q[1];
sx q[1];
rz(-2.746114) q[1];
sx q[1];
rz(-1.5979824) q[1];
rz(-1.8120469) q[3];
sx q[3];
rz(-2.6729306) q[3];
sx q[3];
rz(1.4392343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8362063) q[2];
sx q[2];
rz(-2.6658194) q[2];
sx q[2];
rz(-0.22030182) q[2];
rz(-2.3588755) q[3];
sx q[3];
rz(-1.3681151) q[3];
sx q[3];
rz(-0.14557423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43133217) q[0];
sx q[0];
rz(-2.1556957) q[0];
sx q[0];
rz(1.3166991) q[0];
rz(2.6915164) q[1];
sx q[1];
rz(-1.4349667) q[1];
sx q[1];
rz(-3.0302474) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6346719) q[0];
sx q[0];
rz(-1.823018) q[0];
sx q[0];
rz(0.66703779) q[0];
rz(-pi) q[1];
rz(1.8315122) q[2];
sx q[2];
rz(-0.33277724) q[2];
sx q[2];
rz(0.24713255) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.01035168) q[1];
sx q[1];
rz(-1.8496173) q[1];
sx q[1];
rz(2.1479021) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5841949) q[3];
sx q[3];
rz(-1.9950997) q[3];
sx q[3];
rz(2.4941924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4075809) q[2];
sx q[2];
rz(-2.0110726) q[2];
sx q[2];
rz(3.0111175) q[2];
rz(-2.981251) q[3];
sx q[3];
rz(-0.66157833) q[3];
sx q[3];
rz(-2.9714238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57123667) q[0];
sx q[0];
rz(-2.9354876) q[0];
sx q[0];
rz(-0.11827949) q[0];
rz(2.2453399) q[1];
sx q[1];
rz(-1.4786485) q[1];
sx q[1];
rz(0.55312696) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7161067) q[0];
sx q[0];
rz(-2.5488459) q[0];
sx q[0];
rz(0.61381189) q[0];
rz(-pi) q[1];
x q[1];
rz(0.73605717) q[2];
sx q[2];
rz(-1.9070574) q[2];
sx q[2];
rz(0.93341174) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0748079) q[1];
sx q[1];
rz(-0.60943595) q[1];
sx q[1];
rz(-1.273073) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0712264) q[3];
sx q[3];
rz(-2.4132846) q[3];
sx q[3];
rz(-2.0508931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3502189) q[2];
sx q[2];
rz(-0.72177902) q[2];
sx q[2];
rz(-2.1283894) q[2];
rz(-0.31747097) q[3];
sx q[3];
rz(-1.443913) q[3];
sx q[3];
rz(2.1875994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0233699) q[0];
sx q[0];
rz(-0.88339266) q[0];
sx q[0];
rz(2.6364442) q[0];
rz(-0.39816868) q[1];
sx q[1];
rz(-1.8759556) q[1];
sx q[1];
rz(1.8123951) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1627462) q[0];
sx q[0];
rz(-1.1067451) q[0];
sx q[0];
rz(-0.18382182) q[0];
x q[1];
rz(-1.9301435) q[2];
sx q[2];
rz(-2.076175) q[2];
sx q[2];
rz(1.1352678) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8230629) q[1];
sx q[1];
rz(-1.6803871) q[1];
sx q[1];
rz(-1.0062045) q[1];
rz(-1.2990713) q[3];
sx q[3];
rz(-1.8193805) q[3];
sx q[3];
rz(0.18928537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.053981606) q[2];
sx q[2];
rz(-1.3302646) q[2];
sx q[2];
rz(-0.83016738) q[2];
rz(-1.6603445) q[3];
sx q[3];
rz(-1.0426499) q[3];
sx q[3];
rz(-1.9771077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8692577) q[0];
sx q[0];
rz(-0.91599688) q[0];
sx q[0];
rz(-1.2087615) q[0];
rz(0.2441949) q[1];
sx q[1];
rz(-1.9701651) q[1];
sx q[1];
rz(-0.29921439) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0335145) q[0];
sx q[0];
rz(-0.21768269) q[0];
sx q[0];
rz(1.6924477) q[0];
rz(-2.5317279) q[2];
sx q[2];
rz(-2.0195896) q[2];
sx q[2];
rz(-0.14022174) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.17822325) q[1];
sx q[1];
rz(-0.91919628) q[1];
sx q[1];
rz(-0.061151932) q[1];
x q[2];
rz(-2.4387909) q[3];
sx q[3];
rz(-1.2497823) q[3];
sx q[3];
rz(1.7197844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9485335) q[2];
sx q[2];
rz(-1.4024573) q[2];
sx q[2];
rz(2.194727) q[2];
rz(1.7638505) q[3];
sx q[3];
rz(-2.6008714) q[3];
sx q[3];
rz(2.4989959) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68719012) q[0];
sx q[0];
rz(-2.7070422) q[0];
sx q[0];
rz(-0.53892556) q[0];
rz(-2.9680805) q[1];
sx q[1];
rz(-2.3850846) q[1];
sx q[1];
rz(-2.7573746) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0464041) q[0];
sx q[0];
rz(-2.6676819) q[0];
sx q[0];
rz(-0.35081151) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.27663934) q[2];
sx q[2];
rz(-2.5211272) q[2];
sx q[2];
rz(2.7850604) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9335971) q[1];
sx q[1];
rz(-1.5969445) q[1];
sx q[1];
rz(1.6297865) q[1];
rz(-1.614607) q[3];
sx q[3];
rz(-0.5543983) q[3];
sx q[3];
rz(1.8917494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3197799) q[2];
sx q[2];
rz(-0.76924789) q[2];
sx q[2];
rz(-2.5452851) q[2];
rz(-2.1042306) q[3];
sx q[3];
rz(-0.77330247) q[3];
sx q[3];
rz(1.1843225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7647112) q[0];
sx q[0];
rz(-0.75279555) q[0];
sx q[0];
rz(-2.9994614) q[0];
rz(-2.0243952) q[1];
sx q[1];
rz(-1.6551207) q[1];
sx q[1];
rz(1.9140917) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27782521) q[0];
sx q[0];
rz(-1.2949048) q[0];
sx q[0];
rz(-0.49580144) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6864482) q[2];
sx q[2];
rz(-0.1384379) q[2];
sx q[2];
rz(-1.9659245) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.043316226) q[1];
sx q[1];
rz(-2.195561) q[1];
sx q[1];
rz(-3.0618145) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7015721) q[3];
sx q[3];
rz(-0.87910637) q[3];
sx q[3];
rz(-1.7753851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.050798945) q[2];
sx q[2];
rz(-2.8992081) q[2];
sx q[2];
rz(2.3432815) q[2];
rz(1.5618886) q[3];
sx q[3];
rz(-1.7589933) q[3];
sx q[3];
rz(-0.17670512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(2.0853737) q[0];
sx q[0];
rz(-0.64301411) q[0];
sx q[0];
rz(-3.1363078) q[0];
rz(-3.058744) q[1];
sx q[1];
rz(-1.1033892) q[1];
sx q[1];
rz(-2.8740035) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9339211) q[0];
sx q[0];
rz(-2.2333642) q[0];
sx q[0];
rz(1.1712058) q[0];
rz(-pi) q[1];
rz(-0.29163714) q[2];
sx q[2];
rz(-1.1979629) q[2];
sx q[2];
rz(0.68260869) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5983683) q[1];
sx q[1];
rz(-2.7562047) q[1];
sx q[1];
rz(-0.97167991) q[1];
x q[2];
rz(-2.6676099) q[3];
sx q[3];
rz(-2.0360395) q[3];
sx q[3];
rz(1.3657692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5679607) q[2];
sx q[2];
rz(-0.48144123) q[2];
sx q[2];
rz(-1.2130223) q[2];
rz(1.8968449) q[3];
sx q[3];
rz(-2.6601578) q[3];
sx q[3];
rz(1.1904233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5382814) q[0];
sx q[0];
rz(-2.17157) q[0];
sx q[0];
rz(-3.0441913) q[0];
rz(-0.010460336) q[1];
sx q[1];
rz(-0.86581007) q[1];
sx q[1];
rz(2.5444358) q[1];
rz(-1.6594748) q[2];
sx q[2];
rz(-1.0368928) q[2];
sx q[2];
rz(-1.3925585) q[2];
rz(2.1599471) q[3];
sx q[3];
rz(-1.4245778) q[3];
sx q[3];
rz(2.1733976) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
