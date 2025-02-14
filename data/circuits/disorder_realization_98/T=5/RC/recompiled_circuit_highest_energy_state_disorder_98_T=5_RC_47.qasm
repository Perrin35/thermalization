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
rz(-0.67796081) q[0];
sx q[0];
rz(-0.27422658) q[0];
sx q[0];
rz(1.8507313) q[0];
rz(-0.21386799) q[1];
sx q[1];
rz(-2.8254421) q[1];
sx q[1];
rz(1.4996127) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8052657) q[0];
sx q[0];
rz(-1.299843) q[0];
sx q[0];
rz(1.3265557) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2798645) q[2];
sx q[2];
rz(-0.81255823) q[2];
sx q[2];
rz(-1.9951374) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5668402) q[1];
sx q[1];
rz(-0.87927239) q[1];
sx q[1];
rz(1.994446) q[1];
x q[2];
rz(-2.2003907) q[3];
sx q[3];
rz(-1.968002) q[3];
sx q[3];
rz(-1.8091699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9700254) q[2];
sx q[2];
rz(-1.0415404) q[2];
sx q[2];
rz(-2.2402666) q[2];
rz(-0.46402913) q[3];
sx q[3];
rz(-1.8601067) q[3];
sx q[3];
rz(3.071781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.0483911) q[0];
sx q[0];
rz(-3.1049325) q[0];
sx q[0];
rz(-1.2777591) q[0];
rz(0.094206421) q[1];
sx q[1];
rz(-0.46368805) q[1];
sx q[1];
rz(-1.6069848) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33714715) q[0];
sx q[0];
rz(-0.60071105) q[0];
sx q[0];
rz(-1.5128193) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8707451) q[2];
sx q[2];
rz(-0.19057628) q[2];
sx q[2];
rz(-1.9924763) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6834641) q[1];
sx q[1];
rz(-2.5741815) q[1];
sx q[1];
rz(-0.50000425) q[1];
rz(-pi) q[2];
rz(-0.8022763) q[3];
sx q[3];
rz(-0.55219383) q[3];
sx q[3];
rz(1.7084165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4216807) q[2];
sx q[2];
rz(-1.2016502) q[2];
sx q[2];
rz(0.16033944) q[2];
rz(-0.73873377) q[3];
sx q[3];
rz(-0.84638798) q[3];
sx q[3];
rz(1.4275449) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5147603) q[0];
sx q[0];
rz(-1.4381831) q[0];
sx q[0];
rz(0.50977388) q[0];
rz(1.431142) q[1];
sx q[1];
rz(-1.9606083) q[1];
sx q[1];
rz(-0.74657718) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4228446) q[0];
sx q[0];
rz(-1.5507409) q[0];
sx q[0];
rz(-1.0182747) q[0];
rz(-pi) q[1];
rz(-2.4785537) q[2];
sx q[2];
rz(-0.84868492) q[2];
sx q[2];
rz(-2.7977365) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.4697527) q[1];
sx q[1];
rz(-2.3418509) q[1];
sx q[1];
rz(2.0861162) q[1];
rz(2.2516817) q[3];
sx q[3];
rz(-1.8783675) q[3];
sx q[3];
rz(-2.6079082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9868682) q[2];
sx q[2];
rz(-0.26686033) q[2];
sx q[2];
rz(1.9179087) q[2];
rz(-2.0679421) q[3];
sx q[3];
rz(-1.7045538) q[3];
sx q[3];
rz(0.8980208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2735485) q[0];
sx q[0];
rz(-0.88669625) q[0];
sx q[0];
rz(-2.7622188) q[0];
rz(-2.9647478) q[1];
sx q[1];
rz(-1.6601446) q[1];
sx q[1];
rz(-0.79536974) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5569755) q[0];
sx q[0];
rz(-1.4780227) q[0];
sx q[0];
rz(1.2014821) q[0];
rz(-2.3248939) q[2];
sx q[2];
rz(-1.3743322) q[2];
sx q[2];
rz(2.8981371) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6713555) q[1];
sx q[1];
rz(-1.0507601) q[1];
sx q[1];
rz(0.29754559) q[1];
x q[2];
rz(-0.42511149) q[3];
sx q[3];
rz(-0.13775857) q[3];
sx q[3];
rz(-0.83168304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7103601) q[2];
sx q[2];
rz(-1.3453588) q[2];
sx q[2];
rz(-1.7809407) q[2];
rz(-2.1121173) q[3];
sx q[3];
rz(-1.6607213) q[3];
sx q[3];
rz(1.8195389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65115702) q[0];
sx q[0];
rz(-1.5999726) q[0];
sx q[0];
rz(-1.5486451) q[0];
rz(-2.7410638) q[1];
sx q[1];
rz(-1.6533886) q[1];
sx q[1];
rz(-1.7318447) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.627305) q[0];
sx q[0];
rz(-2.3804733) q[0];
sx q[0];
rz(-2.850015) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7821252) q[2];
sx q[2];
rz(-1.7181944) q[2];
sx q[2];
rz(-1.1254252) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.61297902) q[1];
sx q[1];
rz(-0.76957146) q[1];
sx q[1];
rz(-2.2320094) q[1];
rz(-pi) q[2];
rz(-2.1928113) q[3];
sx q[3];
rz(-1.6262486) q[3];
sx q[3];
rz(2.2353719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8283525) q[2];
sx q[2];
rz(-1.7007549) q[2];
sx q[2];
rz(0.43133119) q[2];
rz(3.0865772) q[3];
sx q[3];
rz(-0.31156817) q[3];
sx q[3];
rz(-0.55606786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7677652) q[0];
sx q[0];
rz(-2.8887833) q[0];
sx q[0];
rz(2.1164236) q[0];
rz(-0.45285666) q[1];
sx q[1];
rz(-0.57944524) q[1];
sx q[1];
rz(0.90389171) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4898281) q[0];
sx q[0];
rz(-1.2625041) q[0];
sx q[0];
rz(-0.90109189) q[0];
x q[1];
rz(2.2270062) q[2];
sx q[2];
rz(-1.0433955) q[2];
sx q[2];
rz(1.5114776) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3749668) q[1];
sx q[1];
rz(-1.300525) q[1];
sx q[1];
rz(-2.388776) q[1];
rz(-pi) q[2];
rz(2.3748906) q[3];
sx q[3];
rz(-2.9057716) q[3];
sx q[3];
rz(2.9028149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.30770939) q[2];
sx q[2];
rz(-1.4557975) q[2];
sx q[2];
rz(2.9252388) q[2];
rz(0.89573914) q[3];
sx q[3];
rz(-2.4560865) q[3];
sx q[3];
rz(-1.640813) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1322587) q[0];
sx q[0];
rz(-1.1067156) q[0];
sx q[0];
rz(-2.4427781) q[0];
rz(1.9578594) q[1];
sx q[1];
rz(-1.4212757) q[1];
sx q[1];
rz(-0.85174495) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4119371) q[0];
sx q[0];
rz(-1.0348399) q[0];
sx q[0];
rz(-3.0974749) q[0];
rz(-pi) q[1];
rz(-2.5642407) q[2];
sx q[2];
rz(-1.5492348) q[2];
sx q[2];
rz(0.93761629) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5761583) q[1];
sx q[1];
rz(-1.9322188) q[1];
sx q[1];
rz(-0.71908497) q[1];
rz(-pi) q[2];
rz(-0.61148879) q[3];
sx q[3];
rz(-0.73966714) q[3];
sx q[3];
rz(-0.31114331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.304004) q[2];
sx q[2];
rz(-1.3112661) q[2];
sx q[2];
rz(-0.50789976) q[2];
rz(1.0287644) q[3];
sx q[3];
rz(-2.2977836) q[3];
sx q[3];
rz(0.60104162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-0.69865882) q[0];
sx q[0];
rz(-1.3154987) q[0];
sx q[0];
rz(-0.0096631924) q[0];
rz(1.6161605) q[1];
sx q[1];
rz(-1.7639672) q[1];
sx q[1];
rz(-0.38472167) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54668173) q[0];
sx q[0];
rz(-1.4230886) q[0];
sx q[0];
rz(-1.8099996) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.184395) q[2];
sx q[2];
rz(-1.3829872) q[2];
sx q[2];
rz(-1.1593429) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.87054756) q[1];
sx q[1];
rz(-1.3636075) q[1];
sx q[1];
rz(-1.8533587) q[1];
rz(0.99920239) q[3];
sx q[3];
rz(-2.2561142) q[3];
sx q[3];
rz(0.1911605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.5821417) q[2];
sx q[2];
rz(-2.1820575) q[2];
sx q[2];
rz(3.1257296) q[2];
rz(-1.9150241) q[3];
sx q[3];
rz(-2.3358986) q[3];
sx q[3];
rz(-2.5198643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3847619) q[0];
sx q[0];
rz(-2.4767196) q[0];
sx q[0];
rz(2.050198) q[0];
rz(2.3233844) q[1];
sx q[1];
rz(-1.810377) q[1];
sx q[1];
rz(-0.6689201) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10169928) q[0];
sx q[0];
rz(-2.153206) q[0];
sx q[0];
rz(-0.82836133) q[0];
x q[1];
rz(2.2833334) q[2];
sx q[2];
rz(-1.2405292) q[2];
sx q[2];
rz(-2.1338322) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0099677) q[1];
sx q[1];
rz(-2.5796081) q[1];
sx q[1];
rz(-1.4373006) q[1];
x q[2];
rz(-0.33031611) q[3];
sx q[3];
rz(-2.1282176) q[3];
sx q[3];
rz(1.6176318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6341256) q[2];
sx q[2];
rz(-0.52906817) q[2];
sx q[2];
rz(0.96366209) q[2];
rz(1.4108747) q[3];
sx q[3];
rz(-1.7129292) q[3];
sx q[3];
rz(1.7760407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.016634781) q[0];
sx q[0];
rz(-2.5732915) q[0];
sx q[0];
rz(-1.4547263) q[0];
rz(-1.1085054) q[1];
sx q[1];
rz(-1.5364372) q[1];
sx q[1];
rz(-0.32807168) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5660368) q[0];
sx q[0];
rz(-1.456161) q[0];
sx q[0];
rz(3.0835955) q[0];
rz(1.6871618) q[2];
sx q[2];
rz(-1.7504217) q[2];
sx q[2];
rz(0.76919523) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.1075503) q[1];
sx q[1];
rz(-2.3880516) q[1];
sx q[1];
rz(-0.57489245) q[1];
x q[2];
rz(0.107268) q[3];
sx q[3];
rz(-1.3984496) q[3];
sx q[3];
rz(-1.8037667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.46128094) q[2];
sx q[2];
rz(-1.842247) q[2];
sx q[2];
rz(-0.4168365) q[2];
rz(-2.4428115) q[3];
sx q[3];
rz(-0.67658934) q[3];
sx q[3];
rz(-0.39066395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1998491) q[0];
sx q[0];
rz(-1.9244292) q[0];
sx q[0];
rz(1.695965) q[0];
rz(-2.4327714) q[1];
sx q[1];
rz(-0.91239057) q[1];
sx q[1];
rz(-0.053587996) q[1];
rz(-1.3016635) q[2];
sx q[2];
rz(-0.65432815) q[2];
sx q[2];
rz(2.339044) q[2];
rz(-0.56545267) q[3];
sx q[3];
rz(-0.91309091) q[3];
sx q[3];
rz(2.0500195) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
