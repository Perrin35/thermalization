OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2064535) q[0];
sx q[0];
rz(-0.78092617) q[0];
sx q[0];
rz(2.934802) q[0];
rz(-2.7517154) q[1];
sx q[1];
rz(-2.0808527) q[1];
sx q[1];
rz(-2.9614255) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5313523) q[0];
sx q[0];
rz(-0.58824476) q[0];
sx q[0];
rz(2.504185) q[0];
x q[1];
rz(2.6287659) q[2];
sx q[2];
rz(-1.4564118) q[2];
sx q[2];
rz(2.0778542) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.7000007) q[1];
sx q[1];
rz(-2.4543426) q[1];
sx q[1];
rz(1.6874466) q[1];
x q[2];
rz(-0.22523017) q[3];
sx q[3];
rz(-2.4472144) q[3];
sx q[3];
rz(0.089689342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7303598) q[2];
sx q[2];
rz(-2.2683472) q[2];
sx q[2];
rz(-1.8475378) q[2];
rz(-2.7358352) q[3];
sx q[3];
rz(-1.6399222) q[3];
sx q[3];
rz(2.7348203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2186573) q[0];
sx q[0];
rz(-2.1056392) q[0];
sx q[0];
rz(-0.12582114) q[0];
rz(-0.80548349) q[1];
sx q[1];
rz(-2.3352354) q[1];
sx q[1];
rz(1.7696101) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94700891) q[0];
sx q[0];
rz(-0.8964552) q[0];
sx q[0];
rz(-1.6499004) q[0];
rz(-pi) q[1];
rz(2.9397474) q[2];
sx q[2];
rz(-2.5864374) q[2];
sx q[2];
rz(-1.3943878) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5139211) q[1];
sx q[1];
rz(-1.5253592) q[1];
sx q[1];
rz(-2.8915878) q[1];
x q[2];
rz(1.814718) q[3];
sx q[3];
rz(-1.1141277) q[3];
sx q[3];
rz(-1.4671385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.25386086) q[2];
sx q[2];
rz(-2.4527206) q[2];
sx q[2];
rz(-0.99622336) q[2];
rz(-2.0455202) q[3];
sx q[3];
rz(-1.6888065) q[3];
sx q[3];
rz(-0.85038275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(0.4662194) q[0];
sx q[0];
rz(-0.96005625) q[0];
sx q[0];
rz(-1.0590142) q[0];
rz(1.9937218) q[1];
sx q[1];
rz(-0.73906001) q[1];
sx q[1];
rz(-2.0770729) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76751906) q[0];
sx q[0];
rz(-0.78157434) q[0];
sx q[0];
rz(-2.563345) q[0];
x q[1];
rz(-1.1883005) q[2];
sx q[2];
rz(-1.4244392) q[2];
sx q[2];
rz(2.086703) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3036084) q[1];
sx q[1];
rz(-1.8808937) q[1];
sx q[1];
rz(3.1293948) q[1];
rz(0.84633175) q[3];
sx q[3];
rz(-1.6958106) q[3];
sx q[3];
rz(0.66305977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0718096) q[2];
sx q[2];
rz(-1.1867384) q[2];
sx q[2];
rz(-2.8783669) q[2];
rz(-1.1188544) q[3];
sx q[3];
rz(-1.2759195) q[3];
sx q[3];
rz(-2.3222205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6333106) q[0];
sx q[0];
rz(-3.0968956) q[0];
sx q[0];
rz(1.3942962) q[0];
rz(-2.1030203) q[1];
sx q[1];
rz(-1.690052) q[1];
sx q[1];
rz(1.0669473) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70908961) q[0];
sx q[0];
rz(-1.5320008) q[0];
sx q[0];
rz(2.1006656) q[0];
x q[1];
rz(-1.489584) q[2];
sx q[2];
rz(-1.3021384) q[2];
sx q[2];
rz(3.0693698) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9001138) q[1];
sx q[1];
rz(-2.334324) q[1];
sx q[1];
rz(-0.57358731) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.22575836) q[3];
sx q[3];
rz(-2.0421713) q[3];
sx q[3];
rz(-0.86172047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2944494) q[2];
sx q[2];
rz(-1.7913982) q[2];
sx q[2];
rz(0.28953826) q[2];
rz(-0.52044049) q[3];
sx q[3];
rz(-1.8013022) q[3];
sx q[3];
rz(-0.7590487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6510058) q[0];
sx q[0];
rz(-0.0087954272) q[0];
sx q[0];
rz(0.83475137) q[0];
rz(-1.2443776) q[1];
sx q[1];
rz(-1.8908187) q[1];
sx q[1];
rz(-1.7117737) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5127038) q[0];
sx q[0];
rz(-1.5682724) q[0];
sx q[0];
rz(-0.062688962) q[0];
rz(-pi) q[1];
x q[1];
rz(0.28508913) q[2];
sx q[2];
rz(-1.8008917) q[2];
sx q[2];
rz(2.6094112) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6068078) q[1];
sx q[1];
rz(-0.77960289) q[1];
sx q[1];
rz(-1.603655) q[1];
rz(2.963644) q[3];
sx q[3];
rz(-0.92015172) q[3];
sx q[3];
rz(-1.734317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4804068) q[2];
sx q[2];
rz(-1.0057665) q[2];
sx q[2];
rz(1.8583813) q[2];
rz(0.13218203) q[3];
sx q[3];
rz(-2.81288) q[3];
sx q[3];
rz(-2.8018518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4955687) q[0];
sx q[0];
rz(-0.060957242) q[0];
sx q[0];
rz(-0.47750372) q[0];
rz(1.6409138) q[1];
sx q[1];
rz(-1.532282) q[1];
sx q[1];
rz(-2.5710411) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7659089) q[0];
sx q[0];
rz(-1.8609957) q[0];
sx q[0];
rz(1.1955839) q[0];
rz(2.4628377) q[2];
sx q[2];
rz(-2.5828913) q[2];
sx q[2];
rz(0.15685454) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3444896) q[1];
sx q[1];
rz(-2.2625766) q[1];
sx q[1];
rz(0.89747353) q[1];
rz(-pi) q[2];
rz(-0.029383226) q[3];
sx q[3];
rz(-0.64850649) q[3];
sx q[3];
rz(2.9103968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.70696124) q[2];
sx q[2];
rz(-2.091445) q[2];
sx q[2];
rz(-0.79664191) q[2];
rz(2.8213275) q[3];
sx q[3];
rz(-1.0864778) q[3];
sx q[3];
rz(-1.2789352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0758078) q[0];
sx q[0];
rz(-2.1053574) q[0];
sx q[0];
rz(2.7835223) q[0];
rz(0.2886731) q[1];
sx q[1];
rz(-0.52983785) q[1];
sx q[1];
rz(-1.3105062) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.34328) q[0];
sx q[0];
rz(-2.5710921) q[0];
sx q[0];
rz(-0.35240726) q[0];
rz(-2.1475122) q[2];
sx q[2];
rz(-2.1834282) q[2];
sx q[2];
rz(2.6660369) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.14168921) q[1];
sx q[1];
rz(-1.1232166) q[1];
sx q[1];
rz(2.5738641) q[1];
rz(-pi) q[2];
rz(0.86705039) q[3];
sx q[3];
rz(-0.82074814) q[3];
sx q[3];
rz(-1.5501319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8075809) q[2];
sx q[2];
rz(-0.58471218) q[2];
sx q[2];
rz(-0.99747783) q[2];
rz(-0.57972646) q[3];
sx q[3];
rz(-2.4119191) q[3];
sx q[3];
rz(-1.4355481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8758133) q[0];
sx q[0];
rz(-1.6515323) q[0];
sx q[0];
rz(0.34061256) q[0];
rz(1.2365201) q[1];
sx q[1];
rz(-2.3842936) q[1];
sx q[1];
rz(-1.1901201) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1163568) q[0];
sx q[0];
rz(-1.3394757) q[0];
sx q[0];
rz(-3.0781834) q[0];
rz(-0.068433381) q[2];
sx q[2];
rz(-1.4901731) q[2];
sx q[2];
rz(-1.2863408) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.78672532) q[1];
sx q[1];
rz(-0.39199542) q[1];
sx q[1];
rz(2.3945432) q[1];
rz(-2.6885919) q[3];
sx q[3];
rz(-2.0815597) q[3];
sx q[3];
rz(-2.9651027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.74165806) q[2];
sx q[2];
rz(-0.88230336) q[2];
sx q[2];
rz(-2.7271872) q[2];
rz(1.3828145) q[3];
sx q[3];
rz(-1.7410991) q[3];
sx q[3];
rz(-2.8167021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25093108) q[0];
sx q[0];
rz(-2.119976) q[0];
sx q[0];
rz(-0.1517621) q[0];
rz(-1.4258619) q[1];
sx q[1];
rz(-1.7646004) q[1];
sx q[1];
rz(-0.94917667) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28578368) q[0];
sx q[0];
rz(-1.4306418) q[0];
sx q[0];
rz(-2.6967718) q[0];
rz(-1.3504215) q[2];
sx q[2];
rz(-2.1890963) q[2];
sx q[2];
rz(-0.097620336) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6712499) q[1];
sx q[1];
rz(-1.8103231) q[1];
sx q[1];
rz(2.4688979) q[1];
rz(-pi) q[2];
rz(-0.48524951) q[3];
sx q[3];
rz(-2.1059548) q[3];
sx q[3];
rz(2.8458965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7342547) q[2];
sx q[2];
rz(-2.7463425) q[2];
sx q[2];
rz(2.7091743) q[2];
rz(1.5405103) q[3];
sx q[3];
rz(-1.6316905) q[3];
sx q[3];
rz(0.66175118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(1.0712873) q[0];
sx q[0];
rz(-2.8685331) q[0];
sx q[0];
rz(-0.37316698) q[0];
rz(-0.35692731) q[1];
sx q[1];
rz(-2.280805) q[1];
sx q[1];
rz(-0.7235136) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7454119) q[0];
sx q[0];
rz(-1.6999082) q[0];
sx q[0];
rz(-0.73208916) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4805484) q[2];
sx q[2];
rz(-2.2482578) q[2];
sx q[2];
rz(-1.2573164) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8404624) q[1];
sx q[1];
rz(-0.81043078) q[1];
sx q[1];
rz(-2.7644972) q[1];
x q[2];
rz(-3.1240508) q[3];
sx q[3];
rz(-0.27046698) q[3];
sx q[3];
rz(-1.5568352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2202806) q[2];
sx q[2];
rz(-1.3833412) q[2];
sx q[2];
rz(-0.55345654) q[2];
rz(2.4297595) q[3];
sx q[3];
rz(-0.40829855) q[3];
sx q[3];
rz(-2.0994983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2198467) q[0];
sx q[0];
rz(-0.60315673) q[0];
sx q[0];
rz(-3.0043816) q[0];
rz(-0.28221054) q[1];
sx q[1];
rz(-2.3574775) q[1];
sx q[1];
rz(-2.2025253) q[1];
rz(2.4137905) q[2];
sx q[2];
rz(-0.62725485) q[2];
sx q[2];
rz(-2.2841452) q[2];
rz(1.8525193) q[3];
sx q[3];
rz(-1.1829794) q[3];
sx q[3];
rz(1.1869528) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
