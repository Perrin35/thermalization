OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.86201) q[0];
sx q[0];
rz(-0.53505889) q[0];
sx q[0];
rz(0.23393272) q[0];
rz(0.84199953) q[1];
sx q[1];
rz(-1.7276126) q[1];
sx q[1];
rz(1.6230621) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5748225) q[0];
sx q[0];
rz(-1.8810417) q[0];
sx q[0];
rz(2.6022018) q[0];
rz(-pi) q[1];
rz(-2.9286381) q[2];
sx q[2];
rz(-1.1446272) q[2];
sx q[2];
rz(2.1338685) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9649992) q[1];
sx q[1];
rz(-1.0142583) q[1];
sx q[1];
rz(-2.2189369) q[1];
rz(-0.47187658) q[3];
sx q[3];
rz(-0.96376824) q[3];
sx q[3];
rz(0.31479731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5412377) q[2];
sx q[2];
rz(-1.4404094) q[2];
sx q[2];
rz(-2.3360628) q[2];
rz(0.39189288) q[3];
sx q[3];
rz(-1.3561748) q[3];
sx q[3];
rz(2.384757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58981744) q[0];
sx q[0];
rz(-0.68142319) q[0];
sx q[0];
rz(-0.43847325) q[0];
rz(-1.27502) q[1];
sx q[1];
rz(-1.6908815) q[1];
sx q[1];
rz(-0.10800392) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40168821) q[0];
sx q[0];
rz(-1.5956912) q[0];
sx q[0];
rz(-1.4598926) q[0];
x q[1];
rz(-0.63814068) q[2];
sx q[2];
rz(-0.9305939) q[2];
sx q[2];
rz(-1.3253044) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7234791) q[1];
sx q[1];
rz(-1.8239841) q[1];
sx q[1];
rz(-3.0550356) q[1];
rz(-2.0612848) q[3];
sx q[3];
rz(-1.401859) q[3];
sx q[3];
rz(0.47034697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7763623) q[2];
sx q[2];
rz(-0.47351101) q[2];
sx q[2];
rz(3.0253809) q[2];
rz(2.727437) q[3];
sx q[3];
rz(-1.073758) q[3];
sx q[3];
rz(-0.80348408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6666343) q[0];
sx q[0];
rz(-1.8284429) q[0];
sx q[0];
rz(2.3057002) q[0];
rz(-1.6235141) q[1];
sx q[1];
rz(-1.6267136) q[1];
sx q[1];
rz(1.2380884) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9515563) q[0];
sx q[0];
rz(-2.073003) q[0];
sx q[0];
rz(-2.7156208) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.34332163) q[2];
sx q[2];
rz(-1.3385217) q[2];
sx q[2];
rz(0.49988036) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7459384) q[1];
sx q[1];
rz(-1.4534344) q[1];
sx q[1];
rz(3.0836041) q[1];
rz(-pi) q[2];
rz(-2.8594871) q[3];
sx q[3];
rz(-1.8477401) q[3];
sx q[3];
rz(-1.1462584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8572924) q[2];
sx q[2];
rz(-2.9967873) q[2];
sx q[2];
rz(-2.4227552) q[2];
rz(2.5607064) q[3];
sx q[3];
rz(-1.418117) q[3];
sx q[3];
rz(-0.7287997) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.523664) q[0];
sx q[0];
rz(-1.5476462) q[0];
sx q[0];
rz(-2.0148328) q[0];
rz(1.2976546) q[1];
sx q[1];
rz(-1.6512066) q[1];
sx q[1];
rz(-1.0430956) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1820871) q[0];
sx q[0];
rz(-2.8694177) q[0];
sx q[0];
rz(-1.4578117) q[0];
rz(1.1642493) q[2];
sx q[2];
rz(-0.57465982) q[2];
sx q[2];
rz(0.80917796) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9981737) q[1];
sx q[1];
rz(-2.7980248) q[1];
sx q[1];
rz(-2.0433389) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11123379) q[3];
sx q[3];
rz(-2.4832442) q[3];
sx q[3];
rz(1.8452132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.8635233) q[2];
sx q[2];
rz(-1.3845283) q[2];
sx q[2];
rz(-2.1045904) q[2];
rz(-2.1841124) q[3];
sx q[3];
rz(-2.0786395) q[3];
sx q[3];
rz(-2.3069042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1400414) q[0];
sx q[0];
rz(-0.20714864) q[0];
sx q[0];
rz(3.0539883) q[0];
rz(2.7032848) q[1];
sx q[1];
rz(-2.0594845) q[1];
sx q[1];
rz(-1.3137438) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7676108) q[0];
sx q[0];
rz(-2.0519407) q[0];
sx q[0];
rz(2.8767881) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6694267) q[2];
sx q[2];
rz(-1.6079788) q[2];
sx q[2];
rz(1.94869) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1723154) q[1];
sx q[1];
rz(-1.7427497) q[1];
sx q[1];
rz(0.090222619) q[1];
rz(-0.16745407) q[3];
sx q[3];
rz(-1.823638) q[3];
sx q[3];
rz(1.188864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6614512) q[2];
sx q[2];
rz(-1.2378614) q[2];
sx q[2];
rz(-2.5782149) q[2];
rz(-1.9208469) q[3];
sx q[3];
rz(-1.7505587) q[3];
sx q[3];
rz(-0.63108546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0193943) q[0];
sx q[0];
rz(-0.65827426) q[0];
sx q[0];
rz(3.1296545) q[0];
rz(2.7519233) q[1];
sx q[1];
rz(-2.4395112) q[1];
sx q[1];
rz(0.24519244) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.747904) q[0];
sx q[0];
rz(-0.98391279) q[0];
sx q[0];
rz(-2.6176207) q[0];
x q[1];
rz(-0.71227534) q[2];
sx q[2];
rz(-0.33604188) q[2];
sx q[2];
rz(3.0672376) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.17055146) q[1];
sx q[1];
rz(-1.028113) q[1];
sx q[1];
rz(-1.0566864) q[1];
rz(-pi) q[2];
rz(2.526029) q[3];
sx q[3];
rz(-0.48071733) q[3];
sx q[3];
rz(-1.8970888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.37990722) q[2];
sx q[2];
rz(-1.3096389) q[2];
sx q[2];
rz(-1.7725819) q[2];
rz(0.53064972) q[3];
sx q[3];
rz(-2.4574418) q[3];
sx q[3];
rz(1.0859038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.948792) q[0];
sx q[0];
rz(-1.3230319) q[0];
sx q[0];
rz(0.2970933) q[0];
rz(-1.5104712) q[1];
sx q[1];
rz(-0.62989569) q[1];
sx q[1];
rz(-2.7105601) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1265701) q[0];
sx q[0];
rz(-1.7807142) q[0];
sx q[0];
rz(-2.8464343) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6234342) q[2];
sx q[2];
rz(-1.3402914) q[2];
sx q[2];
rz(-0.33256691) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.7073133) q[1];
sx q[1];
rz(-0.61798862) q[1];
sx q[1];
rz(-2.1531578) q[1];
rz(-pi) q[2];
x q[2];
rz(0.2934612) q[3];
sx q[3];
rz(-2.0803422) q[3];
sx q[3];
rz(-0.58125416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3057574) q[2];
sx q[2];
rz(-2.367986) q[2];
sx q[2];
rz(2.8744899) q[2];
rz(-0.032020656) q[3];
sx q[3];
rz(-1.1705541) q[3];
sx q[3];
rz(0.10572461) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.06374643) q[0];
sx q[0];
rz(-1.2910605) q[0];
sx q[0];
rz(-2.5015976) q[0];
rz(2.2867639) q[1];
sx q[1];
rz(-1.4850441) q[1];
sx q[1];
rz(1.4124426) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8518191) q[0];
sx q[0];
rz(-0.69727325) q[0];
sx q[0];
rz(-0.66582219) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0630025) q[2];
sx q[2];
rz(-1.6160496) q[2];
sx q[2];
rz(2.947888) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.58929491) q[1];
sx q[1];
rz(-1.4865685) q[1];
sx q[1];
rz(-0.34119795) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9962003) q[3];
sx q[3];
rz(-1.914905) q[3];
sx q[3];
rz(-0.95208012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.60014805) q[2];
sx q[2];
rz(-1.4358127) q[2];
sx q[2];
rz(-2.7195462) q[2];
rz(-2.6680434) q[3];
sx q[3];
rz(-0.62255064) q[3];
sx q[3];
rz(0.1951018) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37050978) q[0];
sx q[0];
rz(-0.93733731) q[0];
sx q[0];
rz(-1.7899845) q[0];
rz(-2.8260258) q[1];
sx q[1];
rz(-1.4548929) q[1];
sx q[1];
rz(-1.1531166) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2446072) q[0];
sx q[0];
rz(-0.073052064) q[0];
sx q[0];
rz(0.79985072) q[0];
x q[1];
rz(3.0598203) q[2];
sx q[2];
rz(-0.53682971) q[2];
sx q[2];
rz(0.54807907) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7866384) q[1];
sx q[1];
rz(-1.6980722) q[1];
sx q[1];
rz(2.5580225) q[1];
x q[2];
rz(0.38001506) q[3];
sx q[3];
rz(-1.1855864) q[3];
sx q[3];
rz(-2.4917701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6844668) q[2];
sx q[2];
rz(-1.9292597) q[2];
sx q[2];
rz(2.7678164) q[2];
rz(-1.9155546) q[3];
sx q[3];
rz(-0.91807476) q[3];
sx q[3];
rz(1.8019684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95947295) q[0];
sx q[0];
rz(-2.4612893) q[0];
sx q[0];
rz(1.3207588) q[0];
rz(0.93718115) q[1];
sx q[1];
rz(-1.3359741) q[1];
sx q[1];
rz(2.6722867) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3166817) q[0];
sx q[0];
rz(-1.6880205) q[0];
sx q[0];
rz(3.0646695) q[0];
x q[1];
rz(1.2861757) q[2];
sx q[2];
rz(-0.60065833) q[2];
sx q[2];
rz(2.9702368) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.093897029) q[1];
sx q[1];
rz(-2.9710431) q[1];
sx q[1];
rz(-0.49255126) q[1];
rz(-pi) q[2];
rz(-1.8176043) q[3];
sx q[3];
rz(-1.863552) q[3];
sx q[3];
rz(-0.61652641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4361973) q[2];
sx q[2];
rz(-0.86279482) q[2];
sx q[2];
rz(1.9311284) q[2];
rz(-2.5734899) q[3];
sx q[3];
rz(-1.3848687) q[3];
sx q[3];
rz(2.1657522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2198915) q[0];
sx q[0];
rz(-0.85048631) q[0];
sx q[0];
rz(1.462107) q[0];
rz(-2.2484491) q[1];
sx q[1];
rz(-1.8291263) q[1];
sx q[1];
rz(-1.6222454) q[1];
rz(-0.90441119) q[2];
sx q[2];
rz(-0.36199649) q[2];
sx q[2];
rz(0.71142759) q[2];
rz(-2.1410971) q[3];
sx q[3];
rz(-1.9420997) q[3];
sx q[3];
rz(1.6122769) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
