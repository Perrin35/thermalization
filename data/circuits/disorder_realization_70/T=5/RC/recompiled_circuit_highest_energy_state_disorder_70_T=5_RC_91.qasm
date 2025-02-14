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
rz(-2.1844644) q[0];
sx q[0];
rz(-1.6532093) q[0];
sx q[0];
rz(-1.1487577) q[0];
rz(0.59612885) q[1];
sx q[1];
rz(-0.40607536) q[1];
sx q[1];
rz(1.6869071) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3341951) q[0];
sx q[0];
rz(-0.79803033) q[0];
sx q[0];
rz(-1.0502771) q[0];
rz(-pi) q[1];
rz(-1.6802915) q[2];
sx q[2];
rz(-2.3286162) q[2];
sx q[2];
rz(-2.2535498) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0815989) q[1];
sx q[1];
rz(-2.3761056) q[1];
sx q[1];
rz(0.19123921) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1249073) q[3];
sx q[3];
rz(-0.66185274) q[3];
sx q[3];
rz(-3.007909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.181432) q[2];
sx q[2];
rz(-0.97603193) q[2];
sx q[2];
rz(-1.9350447) q[2];
rz(-1.1613965) q[3];
sx q[3];
rz(-2.5738218) q[3];
sx q[3];
rz(1.5706221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1777765) q[0];
sx q[0];
rz(-2.0160567) q[0];
sx q[0];
rz(2.1695082) q[0];
rz(-2.8731335) q[1];
sx q[1];
rz(-1.4413036) q[1];
sx q[1];
rz(0.45480248) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3773012) q[0];
sx q[0];
rz(-1.2746841) q[0];
sx q[0];
rz(0.13237093) q[0];
x q[1];
rz(1.9399547) q[2];
sx q[2];
rz(-2.7281986) q[2];
sx q[2];
rz(-1.3184551) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0510682) q[1];
sx q[1];
rz(-1.6220434) q[1];
sx q[1];
rz(0.4649051) q[1];
rz(1.0271026) q[3];
sx q[3];
rz(-2.0683401) q[3];
sx q[3];
rz(-3.080201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1099757) q[2];
sx q[2];
rz(-2.4133108) q[2];
sx q[2];
rz(-0.99962437) q[2];
rz(3.055618) q[3];
sx q[3];
rz(-1.5857668) q[3];
sx q[3];
rz(3.1302248) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7407532) q[0];
sx q[0];
rz(-1.4680306) q[0];
sx q[0];
rz(-2.9174347) q[0];
rz(-1.594918) q[1];
sx q[1];
rz(-1.5983351) q[1];
sx q[1];
rz(-1.8348144) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6431893) q[0];
sx q[0];
rz(-1.0224332) q[0];
sx q[0];
rz(1.4275309) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4073657) q[2];
sx q[2];
rz(-0.9817183) q[2];
sx q[2];
rz(-0.40225077) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3114709) q[1];
sx q[1];
rz(-1.1563627) q[1];
sx q[1];
rz(0.45763514) q[1];
rz(-pi) q[2];
rz(2.8257224) q[3];
sx q[3];
rz(-2.4474553) q[3];
sx q[3];
rz(0.71198502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.33009067) q[2];
sx q[2];
rz(-1.5218488) q[2];
sx q[2];
rz(2.1211993) q[2];
rz(1.8118106) q[3];
sx q[3];
rz(-1.1251757) q[3];
sx q[3];
rz(-1.6167195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40600768) q[0];
sx q[0];
rz(-2.0841053) q[0];
sx q[0];
rz(1.1757346) q[0];
rz(0.27578393) q[1];
sx q[1];
rz(-2.3049054) q[1];
sx q[1];
rz(2.5036459) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0887749) q[0];
sx q[0];
rz(-0.034398641) q[0];
sx q[0];
rz(-0.31487353) q[0];
rz(1.3512072) q[2];
sx q[2];
rz(-0.96115005) q[2];
sx q[2];
rz(0.7555421) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.13220683) q[1];
sx q[1];
rz(-1.2458382) q[1];
sx q[1];
rz(-2.4924949) q[1];
rz(-0.27880554) q[3];
sx q[3];
rz(-2.8099217) q[3];
sx q[3];
rz(-1.3607894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5059169) q[2];
sx q[2];
rz(-1.5779147) q[2];
sx q[2];
rz(1.6944616) q[2];
rz(2.9315089) q[3];
sx q[3];
rz(-2.688372) q[3];
sx q[3];
rz(-2.213018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2337445) q[0];
sx q[0];
rz(-2.9892428) q[0];
sx q[0];
rz(1.0688758) q[0];
rz(-1.2999889) q[1];
sx q[1];
rz(-1.4062107) q[1];
sx q[1];
rz(-1.2058421) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17888363) q[0];
sx q[0];
rz(-2.7307352) q[0];
sx q[0];
rz(-1.7715447) q[0];
rz(-3.1320747) q[2];
sx q[2];
rz(-2.2235907) q[2];
sx q[2];
rz(1.643484) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0413549) q[1];
sx q[1];
rz(-0.85548988) q[1];
sx q[1];
rz(3.0042786) q[1];
rz(-pi) q[2];
x q[2];
rz(0.48128328) q[3];
sx q[3];
rz(-1.4247155) q[3];
sx q[3];
rz(1.2969601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.892889) q[2];
sx q[2];
rz(-2.7619669) q[2];
sx q[2];
rz(-3.0730263) q[2];
rz(2.2019703) q[3];
sx q[3];
rz(-1.4087804) q[3];
sx q[3];
rz(1.9127964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3786316) q[0];
sx q[0];
rz(-1.3968422) q[0];
sx q[0];
rz(-0.5214386) q[0];
rz(1.6261082) q[1];
sx q[1];
rz(-1.7518967) q[1];
sx q[1];
rz(-2.5159786) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9245878) q[0];
sx q[0];
rz(-0.5408322) q[0];
sx q[0];
rz(1.6312253) q[0];
x q[1];
rz(-1.4418774) q[2];
sx q[2];
rz(-1.9471418) q[2];
sx q[2];
rz(0.75808816) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.292899) q[1];
sx q[1];
rz(-2.0239415) q[1];
sx q[1];
rz(0.21993665) q[1];
x q[2];
rz(-2.7196521) q[3];
sx q[3];
rz(-1.5152389) q[3];
sx q[3];
rz(-2.0618604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1336512) q[2];
sx q[2];
rz(-0.59591728) q[2];
sx q[2];
rz(0.1055183) q[2];
rz(0.96406913) q[3];
sx q[3];
rz(-1.151842) q[3];
sx q[3];
rz(-0.9849557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.069020011) q[0];
sx q[0];
rz(-0.21518406) q[0];
sx q[0];
rz(1.4129114) q[0];
rz(2.7627796) q[1];
sx q[1];
rz(-1.8870995) q[1];
sx q[1];
rz(0.55317318) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7460497) q[0];
sx q[0];
rz(-0.93727124) q[0];
sx q[0];
rz(2.6468011) q[0];
rz(3.000053) q[2];
sx q[2];
rz(-1.4813378) q[2];
sx q[2];
rz(-0.78023655) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4338816) q[1];
sx q[1];
rz(-1.5305145) q[1];
sx q[1];
rz(-2.8580185) q[1];
rz(-pi) q[2];
rz(-1.8314701) q[3];
sx q[3];
rz(-1.7672774) q[3];
sx q[3];
rz(0.5082265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6806014) q[2];
sx q[2];
rz(-2.1059771) q[2];
sx q[2];
rz(2.7294532) q[2];
rz(-2.4814217) q[3];
sx q[3];
rz(-0.5414525) q[3];
sx q[3];
rz(0.49303833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93765813) q[0];
sx q[0];
rz(-1.0451319) q[0];
sx q[0];
rz(-2.9397021) q[0];
rz(2.0837325) q[1];
sx q[1];
rz(-2.7767534) q[1];
sx q[1];
rz(-0.26022628) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79151151) q[0];
sx q[0];
rz(-0.66257325) q[0];
sx q[0];
rz(-2.4988079) q[0];
x q[1];
rz(-0.8809907) q[2];
sx q[2];
rz(-0.63369232) q[2];
sx q[2];
rz(1.6735759) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7937745) q[1];
sx q[1];
rz(-1.291569) q[1];
sx q[1];
rz(-2.7814968) q[1];
x q[2];
rz(-0.66047858) q[3];
sx q[3];
rz(-1.704543) q[3];
sx q[3];
rz(1.4123358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1883833) q[2];
sx q[2];
rz(-1.4346432) q[2];
sx q[2];
rz(-2.5592213) q[2];
rz(0.44238704) q[3];
sx q[3];
rz(-1.9059537) q[3];
sx q[3];
rz(-0.83627397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22170947) q[0];
sx q[0];
rz(-1.5279122) q[0];
sx q[0];
rz(-1.4269933) q[0];
rz(-1.3555948) q[1];
sx q[1];
rz(-1.6537063) q[1];
sx q[1];
rz(-1.7048763) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5557818) q[0];
sx q[0];
rz(-1.8904875) q[0];
sx q[0];
rz(0.021110708) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.605805) q[2];
sx q[2];
rz(-2.5965207) q[2];
sx q[2];
rz(1.4591914) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1940665) q[1];
sx q[1];
rz(-1.8119651) q[1];
sx q[1];
rz(3.126866) q[1];
rz(-pi) q[2];
rz(-2.0974227) q[3];
sx q[3];
rz(-1.4344525) q[3];
sx q[3];
rz(0.52140331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.23472486) q[2];
sx q[2];
rz(-1.8610672) q[2];
sx q[2];
rz(1.6024164) q[2];
rz(-0.60798821) q[3];
sx q[3];
rz(-1.4092849) q[3];
sx q[3];
rz(2.4517945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7749087) q[0];
sx q[0];
rz(-0.66119778) q[0];
sx q[0];
rz(2.3181424) q[0];
rz(-2.2460294) q[1];
sx q[1];
rz(-1.3500373) q[1];
sx q[1];
rz(-2.7604738) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6247647) q[0];
sx q[0];
rz(-2.4104558) q[0];
sx q[0];
rz(-0.69964377) q[0];
x q[1];
rz(-0.40171774) q[2];
sx q[2];
rz(-1.7251996) q[2];
sx q[2];
rz(0.32333514) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1824519) q[1];
sx q[1];
rz(-2.509626) q[1];
sx q[1];
rz(2.0672634) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1249843) q[3];
sx q[3];
rz(-1.018152) q[3];
sx q[3];
rz(1.2070398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8054008) q[2];
sx q[2];
rz(-0.41173428) q[2];
sx q[2];
rz(2.1545048) q[2];
rz(1.5385212) q[3];
sx q[3];
rz(-0.58731949) q[3];
sx q[3];
rz(2.9015598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22454746) q[0];
sx q[0];
rz(-1.9018835) q[0];
sx q[0];
rz(-1.6214669) q[0];
rz(-2.4317901) q[1];
sx q[1];
rz(-2.5820422) q[1];
sx q[1];
rz(-1.6834264) q[1];
rz(0.29322704) q[2];
sx q[2];
rz(-2.7623889) q[2];
sx q[2];
rz(1.8017906) q[2];
rz(-0.62623528) q[3];
sx q[3];
rz(-2.2953485) q[3];
sx q[3];
rz(2.9957537) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
