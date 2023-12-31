OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5946755) q[0];
sx q[0];
rz(-1.0008873) q[0];
sx q[0];
rz(2.9291908) q[0];
rz(-2.4266333) q[1];
sx q[1];
rz(-0.7874878) q[1];
sx q[1];
rz(1.8600872) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0948439) q[0];
sx q[0];
rz(-1.9108512) q[0];
sx q[0];
rz(2.9283471) q[0];
rz(-pi) q[1];
rz(1.6002866) q[2];
sx q[2];
rz(-0.87848308) q[2];
sx q[2];
rz(-0.32683795) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.233853) q[1];
sx q[1];
rz(-2.6820161) q[1];
sx q[1];
rz(1.185226) q[1];
rz(-pi) q[2];
rz(-1.699502) q[3];
sx q[3];
rz(-1.5081815) q[3];
sx q[3];
rz(-2.6821729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8193801) q[2];
sx q[2];
rz(-3.0266422) q[2];
sx q[2];
rz(-1.6248576) q[2];
rz(0.24762282) q[3];
sx q[3];
rz(-1.473749) q[3];
sx q[3];
rz(2.4860399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.4441929) q[0];
sx q[0];
rz(-1.4780761) q[0];
sx q[0];
rz(-2.9852988) q[0];
rz(0.63931757) q[1];
sx q[1];
rz(-2.1289861) q[1];
sx q[1];
rz(-1.4888391) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1600906) q[0];
sx q[0];
rz(-1.0233876) q[0];
sx q[0];
rz(1.7574969) q[0];
rz(-2.3180914) q[2];
sx q[2];
rz(-2.245792) q[2];
sx q[2];
rz(-2.0958971) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3540346) q[1];
sx q[1];
rz(-0.84317849) q[1];
sx q[1];
rz(-0.80227279) q[1];
rz(-pi) q[2];
rz(-0.21807166) q[3];
sx q[3];
rz(-0.13987386) q[3];
sx q[3];
rz(-1.1034031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.8504101) q[2];
sx q[2];
rz(-3.0427142) q[2];
sx q[2];
rz(-0.29671159) q[2];
rz(-0.3324278) q[3];
sx q[3];
rz(-1.5184831) q[3];
sx q[3];
rz(-1.2338352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6333703) q[0];
sx q[0];
rz(-1.2587661) q[0];
sx q[0];
rz(0.53031522) q[0];
rz(-0.5439533) q[1];
sx q[1];
rz(-2.0201611) q[1];
sx q[1];
rz(1.189032) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9246763) q[0];
sx q[0];
rz(-2.8042256) q[0];
sx q[0];
rz(-0.49305537) q[0];
rz(-pi) q[1];
x q[1];
rz(0.24818111) q[2];
sx q[2];
rz(-1.6096874) q[2];
sx q[2];
rz(-1.789202) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5570453) q[1];
sx q[1];
rz(-2.097633) q[1];
sx q[1];
rz(0.79750632) q[1];
rz(2.6044106) q[3];
sx q[3];
rz(-2.3677285) q[3];
sx q[3];
rz(-2.8230132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.90551463) q[2];
sx q[2];
rz(-1.6685852) q[2];
sx q[2];
rz(-1.2515602) q[2];
rz(0.45423147) q[3];
sx q[3];
rz(-0.43332228) q[3];
sx q[3];
rz(2.7864578) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77804756) q[0];
sx q[0];
rz(-2.8319034) q[0];
sx q[0];
rz(1.3409412) q[0];
rz(2.5355133) q[1];
sx q[1];
rz(-0.93901912) q[1];
sx q[1];
rz(1.9569424) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9979447) q[0];
sx q[0];
rz(-2.8638683) q[0];
sx q[0];
rz(-2.4585637) q[0];
rz(0.84817727) q[2];
sx q[2];
rz(-0.64372534) q[2];
sx q[2];
rz(-1.4959469) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6726917) q[1];
sx q[1];
rz(-0.60201445) q[1];
sx q[1];
rz(-1.7778346) q[1];
rz(2.437856) q[3];
sx q[3];
rz(-0.53547137) q[3];
sx q[3];
rz(0.4243917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9081395) q[2];
sx q[2];
rz(-0.88580695) q[2];
sx q[2];
rz(0.63956368) q[2];
rz(-0.8574287) q[3];
sx q[3];
rz(-1.9968417) q[3];
sx q[3];
rz(-0.48979428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63067591) q[0];
sx q[0];
rz(-2.5998721) q[0];
sx q[0];
rz(0.52039352) q[0];
rz(0.10294542) q[1];
sx q[1];
rz(-1.0504477) q[1];
sx q[1];
rz(0.72881126) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33671185) q[0];
sx q[0];
rz(-1.6687972) q[0];
sx q[0];
rz(0.20551198) q[0];
rz(0.8466709) q[2];
sx q[2];
rz(-1.6933105) q[2];
sx q[2];
rz(-0.12144897) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5908969) q[1];
sx q[1];
rz(-1.3662852) q[1];
sx q[1];
rz(-1.900161) q[1];
rz(-0.54069424) q[3];
sx q[3];
rz(-2.1139507) q[3];
sx q[3];
rz(2.3054996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0985428) q[2];
sx q[2];
rz(-2.7816732) q[2];
sx q[2];
rz(0.70971242) q[2];
rz(-2.5085311) q[3];
sx q[3];
rz(-1.9934306) q[3];
sx q[3];
rz(3.0098797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0387886) q[0];
sx q[0];
rz(-3.0392201) q[0];
sx q[0];
rz(-2.2578755) q[0];
rz(-2.2209514) q[1];
sx q[1];
rz(-1.0620774) q[1];
sx q[1];
rz(-2.9439435) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6186002) q[0];
sx q[0];
rz(-2.1654961) q[0];
sx q[0];
rz(0.80070514) q[0];
rz(-pi) q[1];
rz(2.501776) q[2];
sx q[2];
rz(-3.0090927) q[2];
sx q[2];
rz(-0.25623955) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.27053988) q[1];
sx q[1];
rz(-2.6846243) q[1];
sx q[1];
rz(-3.0448826) q[1];
rz(-2.382155) q[3];
sx q[3];
rz(-1.0348088) q[3];
sx q[3];
rz(-2.939455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1137696) q[2];
sx q[2];
rz(-1.343507) q[2];
sx q[2];
rz(3.0659884) q[2];
rz(1.4893701) q[3];
sx q[3];
rz(-0.39947709) q[3];
sx q[3];
rz(0.82908019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7224715) q[0];
sx q[0];
rz(-2.2561181) q[0];
sx q[0];
rz(-0.042073123) q[0];
rz(-1.178859) q[1];
sx q[1];
rz(-1.983164) q[1];
sx q[1];
rz(-1.9721608) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.021335) q[0];
sx q[0];
rz(-1.406167) q[0];
sx q[0];
rz(0.067532587) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7514236) q[2];
sx q[2];
rz(-1.3069659) q[2];
sx q[2];
rz(-0.90261501) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2028423) q[1];
sx q[1];
rz(-2.4882567) q[1];
sx q[1];
rz(0.57196879) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.38798214) q[3];
sx q[3];
rz(-1.6098607) q[3];
sx q[3];
rz(-2.8475873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3183257) q[2];
sx q[2];
rz(-1.1151168) q[2];
sx q[2];
rz(-3.1271093) q[2];
rz(-1.5197586) q[3];
sx q[3];
rz(-1.2967916) q[3];
sx q[3];
rz(0.55317944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13977215) q[0];
sx q[0];
rz(-2.791239) q[0];
sx q[0];
rz(2.5792504) q[0];
rz(-1.7100122) q[1];
sx q[1];
rz(-1.7446012) q[1];
sx q[1];
rz(0.45773488) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.021488) q[0];
sx q[0];
rz(-1.095626) q[0];
sx q[0];
rz(0.094046353) q[0];
rz(2.5896795) q[2];
sx q[2];
rz(-2.1898824) q[2];
sx q[2];
rz(1.0362831) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9131635) q[1];
sx q[1];
rz(-2.404681) q[1];
sx q[1];
rz(1.5494898) q[1];
rz(-pi) q[2];
rz(0.3433414) q[3];
sx q[3];
rz(-1.0700873) q[3];
sx q[3];
rz(0.2688558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.83595014) q[2];
sx q[2];
rz(-2.2017411) q[2];
sx q[2];
rz(0.28253728) q[2];
rz(-0.95747581) q[3];
sx q[3];
rz(-1.3922393) q[3];
sx q[3];
rz(-2.6628475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95865059) q[0];
sx q[0];
rz(-0.68529737) q[0];
sx q[0];
rz(1.7392993) q[0];
rz(-0.69933403) q[1];
sx q[1];
rz(-1.3032841) q[1];
sx q[1];
rz(1.8797849) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14784797) q[0];
sx q[0];
rz(-1.9980668) q[0];
sx q[0];
rz(-1.2742313) q[0];
x q[1];
rz(-2.4883534) q[2];
sx q[2];
rz(-1.7368894) q[2];
sx q[2];
rz(1.7476029) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0314434) q[1];
sx q[1];
rz(-1.823339) q[1];
sx q[1];
rz(-0.7048216) q[1];
x q[2];
rz(2.7997167) q[3];
sx q[3];
rz(-2.5856527) q[3];
sx q[3];
rz(3.0059909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4385779) q[2];
sx q[2];
rz(-1.5038749) q[2];
sx q[2];
rz(-1.4578488) q[2];
rz(-0.66649377) q[3];
sx q[3];
rz(-1.7137182) q[3];
sx q[3];
rz(2.3898747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89501971) q[0];
sx q[0];
rz(-0.81037766) q[0];
sx q[0];
rz(-1.1599468) q[0];
rz(-3.0874522) q[1];
sx q[1];
rz(-1.4777947) q[1];
sx q[1];
rz(-1.0704401) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1150072) q[0];
sx q[0];
rz(-2.3837649) q[0];
sx q[0];
rz(1.8146145) q[0];
rz(-2.3355867) q[2];
sx q[2];
rz(-1.5984048) q[2];
sx q[2];
rz(-2.1105786) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.66026238) q[1];
sx q[1];
rz(-2.3412625) q[1];
sx q[1];
rz(-0.34923133) q[1];
rz(-2.8937267) q[3];
sx q[3];
rz(-1.7409179) q[3];
sx q[3];
rz(-0.47249139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.76876172) q[2];
sx q[2];
rz(-2.5925345) q[2];
sx q[2];
rz(-1.5268415) q[2];
rz(2.5620143) q[3];
sx q[3];
rz(-1.7947581) q[3];
sx q[3];
rz(-0.65210623) q[3];
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
rz(-2.6282745) q[0];
sx q[0];
rz(-1.6074629) q[0];
sx q[0];
rz(0.71832023) q[0];
rz(2.1075481) q[1];
sx q[1];
rz(-0.68987344) q[1];
sx q[1];
rz(-0.72763163) q[1];
rz(2.3028159) q[2];
sx q[2];
rz(-1.598874) q[2];
sx q[2];
rz(1.4882006) q[2];
rz(-0.53229971) q[3];
sx q[3];
rz(-1.1333864) q[3];
sx q[3];
rz(-0.65896853) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
