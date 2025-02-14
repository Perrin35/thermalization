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
rz(0.53606755) q[0];
sx q[0];
rz(-1.9789088) q[0];
sx q[0];
rz(2.7106078) q[0];
rz(-2.7886136) q[1];
sx q[1];
rz(-1.2231491) q[1];
sx q[1];
rz(-0.42981848) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0864705) q[0];
sx q[0];
rz(-1.4807379) q[0];
sx q[0];
rz(-0.49205972) q[0];
rz(-pi) q[1];
rz(-2.0031646) q[2];
sx q[2];
rz(-1.585344) q[2];
sx q[2];
rz(0.61847875) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7215772) q[1];
sx q[1];
rz(-0.41999451) q[1];
sx q[1];
rz(0.53424044) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9517035) q[3];
sx q[3];
rz(-1.105068) q[3];
sx q[3];
rz(1.7349752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.063804403) q[2];
sx q[2];
rz(-0.52738515) q[2];
sx q[2];
rz(1.0112313) q[2];
rz(2.1003335) q[3];
sx q[3];
rz(-0.51308378) q[3];
sx q[3];
rz(-2.7019555) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0332727) q[0];
sx q[0];
rz(-0.99650538) q[0];
sx q[0];
rz(-0.74747768) q[0];
rz(-0.29391089) q[1];
sx q[1];
rz(-1.1312609) q[1];
sx q[1];
rz(0.25516137) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0145484) q[0];
sx q[0];
rz(-0.60227312) q[0];
sx q[0];
rz(1.8667579) q[0];
rz(-pi) q[1];
rz(1.3841278) q[2];
sx q[2];
rz(-1.871284) q[2];
sx q[2];
rz(1.0050424) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0415062) q[1];
sx q[1];
rz(-1.1250682) q[1];
sx q[1];
rz(0.84560945) q[1];
x q[2];
rz(2.1952755) q[3];
sx q[3];
rz(-1.7229652) q[3];
sx q[3];
rz(2.5077903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1110288) q[2];
sx q[2];
rz(-3.0869637) q[2];
sx q[2];
rz(0.89152208) q[2];
rz(0.24957481) q[3];
sx q[3];
rz(-0.88281837) q[3];
sx q[3];
rz(-2.7454624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4075277) q[0];
sx q[0];
rz(-2.0052818) q[0];
sx q[0];
rz(-3.0974645) q[0];
rz(1.8924425) q[1];
sx q[1];
rz(-0.42611486) q[1];
sx q[1];
rz(-0.485802) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42949781) q[0];
sx q[0];
rz(-0.22616874) q[0];
sx q[0];
rz(-3.0543213) q[0];
rz(2.2052231) q[2];
sx q[2];
rz(-0.23072019) q[2];
sx q[2];
rz(-1.3619193) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9456182) q[1];
sx q[1];
rz(-0.83688078) q[1];
sx q[1];
rz(0.7463781) q[1];
rz(-pi) q[2];
rz(-1.1385192) q[3];
sx q[3];
rz(-1.9345934) q[3];
sx q[3];
rz(1.0897374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8140063) q[2];
sx q[2];
rz(-1.0884716) q[2];
sx q[2];
rz(2.4468454) q[2];
rz(0.57573777) q[3];
sx q[3];
rz(-1.271194) q[3];
sx q[3];
rz(1.7339138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(2.1908252) q[0];
sx q[0];
rz(-2.8468843) q[0];
sx q[0];
rz(-2.8488391) q[0];
rz(-0.5689019) q[1];
sx q[1];
rz(-0.82745537) q[1];
sx q[1];
rz(-2.2946766) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6349061) q[0];
sx q[0];
rz(-1.6483848) q[0];
sx q[0];
rz(1.1617178) q[0];
rz(-1.4796094) q[2];
sx q[2];
rz(-0.6715179) q[2];
sx q[2];
rz(-2.6330122) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.49063045) q[1];
sx q[1];
rz(-1.8086901) q[1];
sx q[1];
rz(1.7693273) q[1];
x q[2];
rz(1.2388703) q[3];
sx q[3];
rz(-2.9369825) q[3];
sx q[3];
rz(1.0910891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.34951052) q[2];
sx q[2];
rz(-1.258472) q[2];
sx q[2];
rz(-0.2087896) q[2];
rz(-2.4236692) q[3];
sx q[3];
rz(-0.28306285) q[3];
sx q[3];
rz(-2.2010402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0402886) q[0];
sx q[0];
rz(-0.52241075) q[0];
sx q[0];
rz(2.8029602) q[0];
rz(0.70308095) q[1];
sx q[1];
rz(-2.4026726) q[1];
sx q[1];
rz(0.83831659) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0378483) q[0];
sx q[0];
rz(-0.5068584) q[0];
sx q[0];
rz(0.54308191) q[0];
x q[1];
rz(-2.9432326) q[2];
sx q[2];
rz(-2.3075929) q[2];
sx q[2];
rz(-1.552358) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1271124) q[1];
sx q[1];
rz(-0.70074425) q[1];
sx q[1];
rz(-1.5190009) q[1];
rz(-pi) q[2];
rz(-0.60689599) q[3];
sx q[3];
rz(-1.2147153) q[3];
sx q[3];
rz(-2.0739561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4579953) q[2];
sx q[2];
rz(-1.0883051) q[2];
sx q[2];
rz(1.0028769) q[2];
rz(-2.9966127) q[3];
sx q[3];
rz(-3.0478015) q[3];
sx q[3];
rz(-0.063044757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1231287) q[0];
sx q[0];
rz(-2.5557684) q[0];
sx q[0];
rz(0.96328324) q[0];
rz(0.18114289) q[1];
sx q[1];
rz(-0.89814848) q[1];
sx q[1];
rz(-1.2394261) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9267777) q[0];
sx q[0];
rz(-0.96450761) q[0];
sx q[0];
rz(-2.0303594) q[0];
rz(-pi) q[1];
rz(-1.2935898) q[2];
sx q[2];
rz(-1.1133901) q[2];
sx q[2];
rz(-0.99062571) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9628004) q[1];
sx q[1];
rz(-1.1988678) q[1];
sx q[1];
rz(-2.7375601) q[1];
rz(-pi) q[2];
rz(-3.0125727) q[3];
sx q[3];
rz(-1.7212528) q[3];
sx q[3];
rz(0.64808382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2565101) q[2];
sx q[2];
rz(-2.2233267) q[2];
sx q[2];
rz(0.93290848) q[2];
rz(-1.7289303) q[3];
sx q[3];
rz(-1.7191073) q[3];
sx q[3];
rz(0.2093813) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6949718) q[0];
sx q[0];
rz(-0.3681204) q[0];
sx q[0];
rz(-2.3387261) q[0];
rz(0.74202263) q[1];
sx q[1];
rz(-1.4890393) q[1];
sx q[1];
rz(-2.7313357) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3543676) q[0];
sx q[0];
rz(-2.4929308) q[0];
sx q[0];
rz(-2.5543268) q[0];
rz(-pi) q[1];
rz(-1.4531936) q[2];
sx q[2];
rz(-0.80999331) q[2];
sx q[2];
rz(-1.9339028) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3222625) q[1];
sx q[1];
rz(-2.9633489) q[1];
sx q[1];
rz(-2.001754) q[1];
x q[2];
rz(2.0461844) q[3];
sx q[3];
rz(-2.1719526) q[3];
sx q[3];
rz(2.9132995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2731169) q[2];
sx q[2];
rz(-1.5584402) q[2];
sx q[2];
rz(-2.9203019) q[2];
rz(0.41915974) q[3];
sx q[3];
rz(-0.48210382) q[3];
sx q[3];
rz(2.6656849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9047423) q[0];
sx q[0];
rz(-2.1433647) q[0];
sx q[0];
rz(-1.5297484) q[0];
rz(-1.3329685) q[1];
sx q[1];
rz(-0.95104853) q[1];
sx q[1];
rz(-1.453368) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8771031) q[0];
sx q[0];
rz(-2.1801688) q[0];
sx q[0];
rz(-0.86848702) q[0];
rz(2.58415) q[2];
sx q[2];
rz(-0.22391437) q[2];
sx q[2];
rz(-2.4049408) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.61078) q[1];
sx q[1];
rz(-1.4207834) q[1];
sx q[1];
rz(0.71041815) q[1];
x q[2];
rz(-2.0734208) q[3];
sx q[3];
rz(-0.31538559) q[3];
sx q[3];
rz(2.634545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.86614418) q[2];
sx q[2];
rz(-0.25442213) q[2];
sx q[2];
rz(3.1308543) q[2];
rz(0.3206611) q[3];
sx q[3];
rz(-1.2977707) q[3];
sx q[3];
rz(-0.58967725) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43309942) q[0];
sx q[0];
rz(-2.1864317) q[0];
sx q[0];
rz(0.47472111) q[0];
rz(0.2933329) q[1];
sx q[1];
rz(-1.1055929) q[1];
sx q[1];
rz(1.6620103) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0728481) q[0];
sx q[0];
rz(-1.5340051) q[0];
sx q[0];
rz(-1.5307013) q[0];
rz(2.4655645) q[2];
sx q[2];
rz(-1.3814298) q[2];
sx q[2];
rz(2.2631136) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8822215) q[1];
sx q[1];
rz(-2.6135635) q[1];
sx q[1];
rz(-2.6709983) q[1];
rz(-2.9607356) q[3];
sx q[3];
rz(-1.81362) q[3];
sx q[3];
rz(2.9202785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.097229615) q[2];
sx q[2];
rz(-1.8624511) q[2];
sx q[2];
rz(-0.43107671) q[2];
rz(0.19032446) q[3];
sx q[3];
rz(-0.35076916) q[3];
sx q[3];
rz(0.88461191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97656074) q[0];
sx q[0];
rz(-2.850552) q[0];
sx q[0];
rz(-2.9227559) q[0];
rz(2.1826375) q[1];
sx q[1];
rz(-2.405165) q[1];
sx q[1];
rz(1.3921907) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.099571596) q[0];
sx q[0];
rz(-2.9416658) q[0];
sx q[0];
rz(2.8801092) q[0];
x q[1];
rz(1.7151681) q[2];
sx q[2];
rz(-1.1544268) q[2];
sx q[2];
rz(-1.8954111) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4569623) q[1];
sx q[1];
rz(-0.70856386) q[1];
sx q[1];
rz(-0.4672017) q[1];
x q[2];
rz(2.1479761) q[3];
sx q[3];
rz(-2.4419067) q[3];
sx q[3];
rz(-2.4260739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2739111) q[2];
sx q[2];
rz(-0.76649222) q[2];
sx q[2];
rz(-1.319818) q[2];
rz(0.53226081) q[3];
sx q[3];
rz(-1.7936734) q[3];
sx q[3];
rz(1.5490612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.559448) q[0];
sx q[0];
rz(-1.5912709) q[0];
sx q[0];
rz(-2.7406319) q[0];
rz(-0.22790146) q[1];
sx q[1];
rz(-1.0453929) q[1];
sx q[1];
rz(1.4452404) q[1];
rz(0.41517045) q[2];
sx q[2];
rz(-1.3317458) q[2];
sx q[2];
rz(1.7785704) q[2];
rz(-1.3263973) q[3];
sx q[3];
rz(-2.5483589) q[3];
sx q[3];
rz(0.15729558) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
