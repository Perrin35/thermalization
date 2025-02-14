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
rz(-4.0487657) q[0];
sx q[0];
rz(2.2450759) q[0];
sx q[0];
rz(9.4656691) q[0];
rz(1.8812802) q[1];
sx q[1];
rz(-0.44668302) q[1];
sx q[1];
rz(-0.094951542) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68317682) q[0];
sx q[0];
rz(-1.9703034) q[0];
sx q[0];
rz(1.3241121) q[0];
rz(1.8008194) q[2];
sx q[2];
rz(-1.025295) q[2];
sx q[2];
rz(-0.70822424) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5490492) q[1];
sx q[1];
rz(-0.84225878) q[1];
sx q[1];
rz(-2.0768354) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.818472) q[3];
sx q[3];
rz(-2.5305037) q[3];
sx q[3];
rz(2.1095534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.2824668) q[2];
sx q[2];
rz(-0.60255113) q[2];
sx q[2];
rz(-1.9762543) q[2];
rz(-2.2230542) q[3];
sx q[3];
rz(-3.0375752) q[3];
sx q[3];
rz(1.9281841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4024046) q[0];
sx q[0];
rz(-0.63888752) q[0];
sx q[0];
rz(0.30493394) q[0];
rz(2.1091499) q[1];
sx q[1];
rz(-2.86125) q[1];
sx q[1];
rz(-0.84411821) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88972487) q[0];
sx q[0];
rz(-1.3322404) q[0];
sx q[0];
rz(-0.15695928) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3621037) q[2];
sx q[2];
rz(-1.6911397) q[2];
sx q[2];
rz(-0.99882767) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5059974) q[1];
sx q[1];
rz(-1.330101) q[1];
sx q[1];
rz(2.7240915) q[1];
x q[2];
rz(-3.1202597) q[3];
sx q[3];
rz(-2.5144684) q[3];
sx q[3];
rz(0.54652484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.9742763) q[2];
sx q[2];
rz(-2.4582489) q[2];
sx q[2];
rz(0.58504504) q[2];
rz(2.282418) q[3];
sx q[3];
rz(-1.1480568) q[3];
sx q[3];
rz(-2.008321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-1.0228731) q[0];
sx q[0];
rz(-2.3206503) q[0];
sx q[0];
rz(0.098544772) q[0];
rz(-0.56602829) q[1];
sx q[1];
rz(-2.2955344) q[1];
sx q[1];
rz(2.6354094) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0289259) q[0];
sx q[0];
rz(-2.5124164) q[0];
sx q[0];
rz(0.58626808) q[0];
rz(-pi) q[1];
rz(-2.3213941) q[2];
sx q[2];
rz(-1.9026105) q[2];
sx q[2];
rz(0.92727509) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.652388) q[1];
sx q[1];
rz(-1.4474157) q[1];
sx q[1];
rz(-1.7921993) q[1];
rz(2.7485756) q[3];
sx q[3];
rz(-2.1073494) q[3];
sx q[3];
rz(2.4483333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.11518662) q[2];
sx q[2];
rz(-1.2135442) q[2];
sx q[2];
rz(-0.085065993) q[2];
rz(0.75299844) q[3];
sx q[3];
rz(-0.66998657) q[3];
sx q[3];
rz(1.5299214) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.531115) q[0];
sx q[0];
rz(-1.6401289) q[0];
sx q[0];
rz(-0.53989545) q[0];
rz(-3.1119697) q[1];
sx q[1];
rz(-0.74631515) q[1];
sx q[1];
rz(1.7866887) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6209117) q[0];
sx q[0];
rz(-1.3518466) q[0];
sx q[0];
rz(2.3205873) q[0];
rz(1.2446885) q[2];
sx q[2];
rz(-2.2549013) q[2];
sx q[2];
rz(-1.5665485) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3554509) q[1];
sx q[1];
rz(-2.1340573) q[1];
sx q[1];
rz(-1.2652629) q[1];
x q[2];
rz(-0.088884066) q[3];
sx q[3];
rz(-2.349472) q[3];
sx q[3];
rz(-1.3490335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4386091) q[2];
sx q[2];
rz(-2.1712124) q[2];
sx q[2];
rz(-0.89540974) q[2];
rz(1.4488975) q[3];
sx q[3];
rz(-1.1619032) q[3];
sx q[3];
rz(2.7321775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92622906) q[0];
sx q[0];
rz(-2.138593) q[0];
sx q[0];
rz(0.0005501752) q[0];
rz(0.5254566) q[1];
sx q[1];
rz(-0.84392396) q[1];
sx q[1];
rz(0.97132436) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7346974) q[0];
sx q[0];
rz(-2.3427123) q[0];
sx q[0];
rz(-2.6968234) q[0];
rz(-pi) q[1];
rz(1.1575562) q[2];
sx q[2];
rz(-1.2451425) q[2];
sx q[2];
rz(-2.535274) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.57572039) q[1];
sx q[1];
rz(-1.0219058) q[1];
sx q[1];
rz(2.4574005) q[1];
rz(-pi) q[2];
rz(2.9266403) q[3];
sx q[3];
rz(-1.8773517) q[3];
sx q[3];
rz(-0.34270278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1725258) q[2];
sx q[2];
rz(-1.0774287) q[2];
sx q[2];
rz(1.3544918) q[2];
rz(-2.5403678) q[3];
sx q[3];
rz(-1.0433334) q[3];
sx q[3];
rz(-1.575298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-2.4300267) q[0];
sx q[0];
rz(-0.30421782) q[0];
sx q[0];
rz(2.8416204) q[0];
rz(-0.9777588) q[1];
sx q[1];
rz(-0.58039665) q[1];
sx q[1];
rz(2.207644) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47563206) q[0];
sx q[0];
rz(-0.36132672) q[0];
sx q[0];
rz(-0.40478171) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63542346) q[2];
sx q[2];
rz(-2.4197856) q[2];
sx q[2];
rz(-1.3363163) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0505414) q[1];
sx q[1];
rz(-0.81650298) q[1];
sx q[1];
rz(1.2495118) q[1];
x q[2];
rz(0.40075366) q[3];
sx q[3];
rz(-1.2448896) q[3];
sx q[3];
rz(-3.000695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0705491) q[2];
sx q[2];
rz(-1.0087548) q[2];
sx q[2];
rz(1.1318995) q[2];
rz(1.2567358) q[3];
sx q[3];
rz(-2.3102424) q[3];
sx q[3];
rz(-1.4164213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8202332) q[0];
sx q[0];
rz(-1.2437404) q[0];
sx q[0];
rz(0.61068049) q[0];
rz(-2.3848379) q[1];
sx q[1];
rz(-0.38833955) q[1];
sx q[1];
rz(1.6441708) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0038549) q[0];
sx q[0];
rz(-3.014999) q[0];
sx q[0];
rz(-0.20124443) q[0];
rz(-2.7779105) q[2];
sx q[2];
rz(-1.1712345) q[2];
sx q[2];
rz(1.3043208) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.037379384) q[1];
sx q[1];
rz(-2.2168405) q[1];
sx q[1];
rz(-1.7039596) q[1];
rz(0.7216924) q[3];
sx q[3];
rz(-1.6462197) q[3];
sx q[3];
rz(1.9957927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.18621592) q[2];
sx q[2];
rz(-2.4775041) q[2];
sx q[2];
rz(-2.5353298) q[2];
rz(2.9210505) q[3];
sx q[3];
rz(-1.8044148) q[3];
sx q[3];
rz(-0.36673275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99899387) q[0];
sx q[0];
rz(-1.6479011) q[0];
sx q[0];
rz(2.2338474) q[0];
rz(1.6084464) q[1];
sx q[1];
rz(-2.1870859) q[1];
sx q[1];
rz(-0.018891637) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7291537) q[0];
sx q[0];
rz(-1.6866418) q[0];
sx q[0];
rz(-0.47286589) q[0];
rz(-0.29952403) q[2];
sx q[2];
rz(-1.5148485) q[2];
sx q[2];
rz(-2.6331226) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0844473) q[1];
sx q[1];
rz(-1.6020163) q[1];
sx q[1];
rz(0.46226661) q[1];
rz(-pi) q[2];
rz(-0.5762866) q[3];
sx q[3];
rz(-0.91720795) q[3];
sx q[3];
rz(1.6393076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7682401) q[2];
sx q[2];
rz(-1.5529996) q[2];
sx q[2];
rz(-2.8949883) q[2];
rz(1.8433579) q[3];
sx q[3];
rz(-1.6876829) q[3];
sx q[3];
rz(-1.2687782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75803718) q[0];
sx q[0];
rz(-1.0724496) q[0];
sx q[0];
rz(2.2776336) q[0];
rz(-1.2147238) q[1];
sx q[1];
rz(-2.0236969) q[1];
sx q[1];
rz(2.5606959) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7621135) q[0];
sx q[0];
rz(-0.91766463) q[0];
sx q[0];
rz(-2.1996621) q[0];
x q[1];
rz(1.6132108) q[2];
sx q[2];
rz(-0.73684947) q[2];
sx q[2];
rz(2.9880817) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5553415) q[1];
sx q[1];
rz(-2.1650392) q[1];
sx q[1];
rz(-1.1583561) q[1];
rz(-pi) q[2];
rz(-0.77928752) q[3];
sx q[3];
rz(-1.7699935) q[3];
sx q[3];
rz(-1.3785386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8860127) q[2];
sx q[2];
rz(-2.6801127) q[2];
sx q[2];
rz(0.18730051) q[2];
rz(0.97697941) q[3];
sx q[3];
rz(-1.5005485) q[3];
sx q[3];
rz(1.8627082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.5115857) q[0];
sx q[0];
rz(-2.4749909) q[0];
sx q[0];
rz(1.2641719) q[0];
rz(2.1695747) q[1];
sx q[1];
rz(-0.7518026) q[1];
sx q[1];
rz(-1.972563) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8147217) q[0];
sx q[0];
rz(-1.7059312) q[0];
sx q[0];
rz(-1.2343455) q[0];
rz(-0.75836597) q[2];
sx q[2];
rz(-2.6347199) q[2];
sx q[2];
rz(3.1257983) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0824521) q[1];
sx q[1];
rz(-0.881221) q[1];
sx q[1];
rz(0.13549681) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6852387) q[3];
sx q[3];
rz(-1.2565359) q[3];
sx q[3];
rz(2.9357167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6675889) q[2];
sx q[2];
rz(-1.7550125) q[2];
sx q[2];
rz(-2.7645195) q[2];
rz(0.22370473) q[3];
sx q[3];
rz(-2.5535899) q[3];
sx q[3];
rz(1.2288176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(1.8857464) q[0];
sx q[0];
rz(-1.3930014) q[0];
sx q[0];
rz(-0.2572671) q[0];
rz(0.59355758) q[1];
sx q[1];
rz(-2.3496353) q[1];
sx q[1];
rz(-1.4812462) q[1];
rz(2.3667546) q[2];
sx q[2];
rz(-0.88817468) q[2];
sx q[2];
rz(-1.9356288) q[2];
rz(1.9695342) q[3];
sx q[3];
rz(-1.8137365) q[3];
sx q[3];
rz(0.89567281) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
