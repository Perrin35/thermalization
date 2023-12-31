OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6991601) q[0];
sx q[0];
rz(4.5259024) q[0];
sx q[0];
rz(10.685267) q[0];
rz(-1.0386382) q[1];
sx q[1];
rz(4.4903978) q[1];
sx q[1];
rz(8.5010565) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1169352) q[0];
sx q[0];
rz(-0.57514656) q[0];
sx q[0];
rz(-2.2206109) q[0];
rz(-pi) q[1];
x q[1];
rz(0.14416868) q[2];
sx q[2];
rz(-1.2906133) q[2];
sx q[2];
rz(-0.35137128) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.87286283) q[1];
sx q[1];
rz(-1.1464719) q[1];
sx q[1];
rz(0.55117589) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6472595) q[3];
sx q[3];
rz(-1.9088073) q[3];
sx q[3];
rz(0.61275834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2279921) q[2];
sx q[2];
rz(-1.8929409) q[2];
sx q[2];
rz(0.16201924) q[2];
rz(-2.2062733) q[3];
sx q[3];
rz(-0.98615065) q[3];
sx q[3];
rz(0.71301618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0682003) q[0];
sx q[0];
rz(-2.91495) q[0];
sx q[0];
rz(1.9447928) q[0];
rz(0.67990047) q[1];
sx q[1];
rz(-2.6459243) q[1];
sx q[1];
rz(1.686036) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14917063) q[0];
sx q[0];
rz(-0.99827168) q[0];
sx q[0];
rz(0.19897977) q[0];
x q[1];
rz(1.4886841) q[2];
sx q[2];
rz(-1.1973235) q[2];
sx q[2];
rz(0.79370802) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1728954) q[1];
sx q[1];
rz(-0.52297938) q[1];
sx q[1];
rz(1.5978659) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5442113) q[3];
sx q[3];
rz(-1.1592602) q[3];
sx q[3];
rz(-0.4707903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.42852795) q[2];
sx q[2];
rz(-1.45168) q[2];
sx q[2];
rz(-1.7896174) q[2];
rz(2.9591566) q[3];
sx q[3];
rz(-0.97674102) q[3];
sx q[3];
rz(0.3119719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(-2.3882554) q[0];
sx q[0];
rz(-2.4607846) q[0];
sx q[0];
rz(-2.341111) q[0];
rz(-3.1128186) q[1];
sx q[1];
rz(-2.0859699) q[1];
sx q[1];
rz(1.172539) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5867509) q[0];
sx q[0];
rz(-1.8790073) q[0];
sx q[0];
rz(-0.59535938) q[0];
rz(-pi) q[1];
rz(2.3861888) q[2];
sx q[2];
rz(-1.3284151) q[2];
sx q[2];
rz(0.6005477) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.4634358) q[1];
sx q[1];
rz(-1.0532182) q[1];
sx q[1];
rz(0.28097681) q[1];
rz(0.058733744) q[3];
sx q[3];
rz(-1.8211094) q[3];
sx q[3];
rz(-0.55610031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0671493) q[2];
sx q[2];
rz(-1.6643486) q[2];
sx q[2];
rz(0.91119901) q[2];
rz(-0.95101142) q[3];
sx q[3];
rz(-2.337303) q[3];
sx q[3];
rz(-0.89200154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38055414) q[0];
sx q[0];
rz(-3.0111713) q[0];
sx q[0];
rz(0.12810853) q[0];
rz(-0.076106636) q[1];
sx q[1];
rz(-1.9271306) q[1];
sx q[1];
rz(-2.6180843) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9273705) q[0];
sx q[0];
rz(-2.4282051) q[0];
sx q[0];
rz(2.5582696) q[0];
rz(-pi) q[1];
rz(-1.4857616) q[2];
sx q[2];
rz(-1.8997314) q[2];
sx q[2];
rz(-2.0563682) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.254863) q[1];
sx q[1];
rz(-2.2543395) q[1];
sx q[1];
rz(-2.8121594) q[1];
x q[2];
rz(-2.162096) q[3];
sx q[3];
rz(-1.6846091) q[3];
sx q[3];
rz(2.1554961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.52544242) q[2];
sx q[2];
rz(-1.5443065) q[2];
sx q[2];
rz(-0.564044) q[2];
rz(-2.8530252) q[3];
sx q[3];
rz(-2.7189062) q[3];
sx q[3];
rz(-2.585876) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6600835) q[0];
sx q[0];
rz(-0.68843377) q[0];
sx q[0];
rz(-1.4915285) q[0];
rz(0.87961698) q[1];
sx q[1];
rz(-1.8938226) q[1];
sx q[1];
rz(0.99194828) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0569699) q[0];
sx q[0];
rz(-1.7449433) q[0];
sx q[0];
rz(1.4625545) q[0];
x q[1];
rz(-1.5993824) q[2];
sx q[2];
rz(-0.30297849) q[2];
sx q[2];
rz(0.44705331) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3494898) q[1];
sx q[1];
rz(-1.2621242) q[1];
sx q[1];
rz(1.5285138) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0117412) q[3];
sx q[3];
rz(-0.81749812) q[3];
sx q[3];
rz(-2.0714456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1296967) q[2];
sx q[2];
rz(-2.7719438) q[2];
sx q[2];
rz(-0.27080718) q[2];
rz(2.9233542) q[3];
sx q[3];
rz(-1.821358) q[3];
sx q[3];
rz(-0.22578421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72702423) q[0];
sx q[0];
rz(-0.71838656) q[0];
sx q[0];
rz(1.3487934) q[0];
rz(-0.38189608) q[1];
sx q[1];
rz(-2.8254639) q[1];
sx q[1];
rz(-1.7165002) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3544918) q[0];
sx q[0];
rz(-2.0658501) q[0];
sx q[0];
rz(-1.3160734) q[0];
x q[1];
rz(0.09128696) q[2];
sx q[2];
rz(-1.316615) q[2];
sx q[2];
rz(-0.19677481) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6188366) q[1];
sx q[1];
rz(-2.6869876) q[1];
sx q[1];
rz(1.7230117) q[1];
rz(-pi) q[2];
rz(-0.53374966) q[3];
sx q[3];
rz(-2.1042049) q[3];
sx q[3];
rz(2.651754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1340593) q[2];
sx q[2];
rz(-2.7441661) q[2];
sx q[2];
rz(2.5777204) q[2];
rz(0.18051906) q[3];
sx q[3];
rz(-1.518395) q[3];
sx q[3];
rz(2.738651) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6329704) q[0];
sx q[0];
rz(-2.9794725) q[0];
sx q[0];
rz(-2.7222743) q[0];
rz(-1.58889) q[1];
sx q[1];
rz(-1.2607375) q[1];
sx q[1];
rz(0.82180506) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0718943) q[0];
sx q[0];
rz(-2.0928203) q[0];
sx q[0];
rz(-0.72899039) q[0];
rz(-pi) q[1];
rz(1.0242545) q[2];
sx q[2];
rz(-2.3356236) q[2];
sx q[2];
rz(-0.96217996) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2591178) q[1];
sx q[1];
rz(-1.3438517) q[1];
sx q[1];
rz(2.4005753) q[1];
rz(-pi) q[2];
rz(-2.6131367) q[3];
sx q[3];
rz(-0.65675694) q[3];
sx q[3];
rz(-2.9627851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8043148) q[2];
sx q[2];
rz(-2.3874805) q[2];
sx q[2];
rz(-2.896893) q[2];
rz(-3.0120567) q[3];
sx q[3];
rz(-1.9774388) q[3];
sx q[3];
rz(1.5130419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41641763) q[0];
sx q[0];
rz(-0.019151909) q[0];
sx q[0];
rz(-2.3186671) q[0];
rz(2.8322463) q[1];
sx q[1];
rz(-1.7495218) q[1];
sx q[1];
rz(-1.8364505) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5236429) q[0];
sx q[0];
rz(-2.4542913) q[0];
sx q[0];
rz(-0.53074093) q[0];
rz(0.70456409) q[2];
sx q[2];
rz(-1.8578055) q[2];
sx q[2];
rz(1.2003843) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0481938) q[1];
sx q[1];
rz(-1.1953925) q[1];
sx q[1];
rz(-2.4673389) q[1];
rz(-pi) q[2];
rz(1.2315024) q[3];
sx q[3];
rz(-1.5929856) q[3];
sx q[3];
rz(-2.0458178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0017073) q[2];
sx q[2];
rz(-1.7557764) q[2];
sx q[2];
rz(-1.4902327) q[2];
rz(-2.0643318) q[3];
sx q[3];
rz(-2.1765985) q[3];
sx q[3];
rz(-0.13154496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3867144) q[0];
sx q[0];
rz(-1.2972378) q[0];
sx q[0];
rz(0.3219147) q[0];
rz(-1.6053258) q[1];
sx q[1];
rz(-1.9202817) q[1];
sx q[1];
rz(-0.70294356) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2625933) q[0];
sx q[0];
rz(-1.2125373) q[0];
sx q[0];
rz(-0.4155638) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0800955) q[2];
sx q[2];
rz(-1.5361538) q[2];
sx q[2];
rz(0.73355567) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.326509) q[1];
sx q[1];
rz(-2.106296) q[1];
sx q[1];
rz(0.19170796) q[1];
rz(-pi) q[2];
rz(0.33396696) q[3];
sx q[3];
rz(-2.1595862) q[3];
sx q[3];
rz(0.86563084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2408509) q[2];
sx q[2];
rz(-2.8618331) q[2];
sx q[2];
rz(1.8019603) q[2];
rz(-2.83589) q[3];
sx q[3];
rz(-1.8140847) q[3];
sx q[3];
rz(-1.3302749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.16383485) q[0];
sx q[0];
rz(-2.3840388) q[0];
sx q[0];
rz(1.2257858) q[0];
rz(-0.90351358) q[1];
sx q[1];
rz(-0.61360306) q[1];
sx q[1];
rz(-0.46863619) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2221453) q[0];
sx q[0];
rz(-1.3627909) q[0];
sx q[0];
rz(-0.39214765) q[0];
rz(-pi) q[1];
rz(-1.956316) q[2];
sx q[2];
rz(-1.8744933) q[2];
sx q[2];
rz(-0.38244837) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.48060265) q[1];
sx q[1];
rz(-2.8669679) q[1];
sx q[1];
rz(1.8250699) q[1];
rz(-pi) q[2];
rz(-0.43379421) q[3];
sx q[3];
rz(-1.9926096) q[3];
sx q[3];
rz(1.0292366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.48352155) q[2];
sx q[2];
rz(-1.8489685) q[2];
sx q[2];
rz(1.998385) q[2];
rz(-0.11463595) q[3];
sx q[3];
rz(-0.95364037) q[3];
sx q[3];
rz(-1.5293998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4951915) q[0];
sx q[0];
rz(-1.1947182) q[0];
sx q[0];
rz(2.4583046) q[0];
rz(-2.519683) q[1];
sx q[1];
rz(-1.4629296) q[1];
sx q[1];
rz(-0.32348979) q[1];
rz(-2.4516104) q[2];
sx q[2];
rz(-0.98914115) q[2];
sx q[2];
rz(3.0422899) q[2];
rz(0.078483742) q[3];
sx q[3];
rz(-2.2208636) q[3];
sx q[3];
rz(0.29490864) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
