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
rz(-2.4009268) q[0];
sx q[0];
rz(-1.8494777) q[0];
sx q[0];
rz(-2.3722755) q[0];
rz(-1.5462592) q[1];
sx q[1];
rz(-0.21472628) q[1];
sx q[1];
rz(0.62155849) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82418782) q[0];
sx q[0];
rz(-2.3185668) q[0];
sx q[0];
rz(-2.4322741) q[0];
x q[1];
rz(-1.5311422) q[2];
sx q[2];
rz(-0.41799212) q[2];
sx q[2];
rz(2.5681873) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4159704) q[1];
sx q[1];
rz(-2.2752004) q[1];
sx q[1];
rz(-2.7654057) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6042799) q[3];
sx q[3];
rz(-2.9859324) q[3];
sx q[3];
rz(0.90720219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8605211) q[2];
sx q[2];
rz(-2.8513384) q[2];
sx q[2];
rz(-0.93112913) q[2];
rz(-0.18533254) q[3];
sx q[3];
rz(-2.1087746) q[3];
sx q[3];
rz(1.7271656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0463878) q[0];
sx q[0];
rz(-2.8075908) q[0];
sx q[0];
rz(-0.97292501) q[0];
rz(-2.3880549) q[1];
sx q[1];
rz(-2.2851508) q[1];
sx q[1];
rz(-3.0389752) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0282183) q[0];
sx q[0];
rz(-1.4732142) q[0];
sx q[0];
rz(3.0563838) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5061812) q[2];
sx q[2];
rz(-2.9757068) q[2];
sx q[2];
rz(-0.36402853) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2778138) q[1];
sx q[1];
rz(-1.4881388) q[1];
sx q[1];
rz(1.7277) q[1];
rz(-pi) q[2];
rz(1.9512734) q[3];
sx q[3];
rz(-0.35934908) q[3];
sx q[3];
rz(0.41216601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7332581) q[2];
sx q[2];
rz(-1.7034986) q[2];
sx q[2];
rz(-2.5737838) q[2];
rz(-0.37219498) q[3];
sx q[3];
rz(-1.8122383) q[3];
sx q[3];
rz(1.6224434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5837625) q[0];
sx q[0];
rz(-0.74302858) q[0];
sx q[0];
rz(1.8999735) q[0];
rz(-0.81635967) q[1];
sx q[1];
rz(-0.38874778) q[1];
sx q[1];
rz(0.98966086) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1802487) q[0];
sx q[0];
rz(-1.3575866) q[0];
sx q[0];
rz(-2.7556987) q[0];
rz(-pi) q[1];
rz(-1.4098566) q[2];
sx q[2];
rz(-1.66203) q[2];
sx q[2];
rz(0.021878069) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1376393) q[1];
sx q[1];
rz(-1.7491566) q[1];
sx q[1];
rz(0.41282005) q[1];
rz(-1.3157522) q[3];
sx q[3];
rz(-2.7560184) q[3];
sx q[3];
rz(0.60690597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.376754) q[2];
sx q[2];
rz(-1.9003442) q[2];
sx q[2];
rz(-1.6468916) q[2];
rz(-2.3946848) q[3];
sx q[3];
rz(-0.61498314) q[3];
sx q[3];
rz(0.99228215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84585369) q[0];
sx q[0];
rz(-2.5339412) q[0];
sx q[0];
rz(2.1920152) q[0];
rz(1.9751366) q[1];
sx q[1];
rz(-1.8828853) q[1];
sx q[1];
rz(2.9393401) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1631781) q[0];
sx q[0];
rz(-1.5401773) q[0];
sx q[0];
rz(-3.091673) q[0];
rz(-pi) q[1];
rz(1.7385029) q[2];
sx q[2];
rz(-2.2532092) q[2];
sx q[2];
rz(-2.3273327) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4829887) q[1];
sx q[1];
rz(-2.3301396) q[1];
sx q[1];
rz(-0.8354018) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1465591) q[3];
sx q[3];
rz(-1.0594133) q[3];
sx q[3];
rz(0.56727876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4919081) q[2];
sx q[2];
rz(-1.4338355) q[2];
sx q[2];
rz(-1.8053619) q[2];
rz(2.6797471) q[3];
sx q[3];
rz(-1.0577842) q[3];
sx q[3];
rz(-2.113078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.860054) q[0];
sx q[0];
rz(-3.0335732) q[0];
sx q[0];
rz(-2.7463013) q[0];
rz(0.95834243) q[1];
sx q[1];
rz(-1.3805362) q[1];
sx q[1];
rz(0.39452943) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4636587) q[0];
sx q[0];
rz(-1.5884841) q[0];
sx q[0];
rz(1.3050446) q[0];
rz(-pi) q[1];
rz(-0.44107159) q[2];
sx q[2];
rz(-1.8243378) q[2];
sx q[2];
rz(2.2833786) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1198152) q[1];
sx q[1];
rz(-1.8938601) q[1];
sx q[1];
rz(-2.3141239) q[1];
rz(-pi) q[2];
rz(0.29629956) q[3];
sx q[3];
rz(-2.5023016) q[3];
sx q[3];
rz(3.1217255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7390274) q[2];
sx q[2];
rz(-1.5274916) q[2];
sx q[2];
rz(-1.019545) q[2];
rz(-1.839365) q[3];
sx q[3];
rz(-3.0890833) q[3];
sx q[3];
rz(-0.70560613) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5810982) q[0];
sx q[0];
rz(-2.5769233) q[0];
sx q[0];
rz(-2.1288921) q[0];
rz(-2.6794491) q[1];
sx q[1];
rz(-1.4116762) q[1];
sx q[1];
rz(-0.046646811) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29829866) q[0];
sx q[0];
rz(-3.0180535) q[0];
sx q[0];
rz(-1.4840829) q[0];
rz(2.88575) q[2];
sx q[2];
rz(-1.2835955) q[2];
sx q[2];
rz(-2.246496) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0477284) q[1];
sx q[1];
rz(-1.0055491) q[1];
sx q[1];
rz(-0.39611343) q[1];
rz(-pi) q[2];
rz(1.9479475) q[3];
sx q[3];
rz(-1.5483678) q[3];
sx q[3];
rz(2.7183661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3095653) q[2];
sx q[2];
rz(-3.053061) q[2];
sx q[2];
rz(0.95477611) q[2];
rz(2.471762) q[3];
sx q[3];
rz(-1.9570743) q[3];
sx q[3];
rz(2.7772016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1466115) q[0];
sx q[0];
rz(-2.8611188) q[0];
sx q[0];
rz(1.9914419) q[0];
rz(3.0448044) q[1];
sx q[1];
rz(-2.0014747) q[1];
sx q[1];
rz(-1.6614301) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87931765) q[0];
sx q[0];
rz(-1.9497383) q[0];
sx q[0];
rz(-2.3728601) q[0];
rz(2.509233) q[2];
sx q[2];
rz(-0.39146921) q[2];
sx q[2];
rz(2.7923358) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.11058775) q[1];
sx q[1];
rz(-1.8613986) q[1];
sx q[1];
rz(0.1433934) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2477996) q[3];
sx q[3];
rz(-2.6867001) q[3];
sx q[3];
rz(-0.98857075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.73411554) q[2];
sx q[2];
rz(-0.87782562) q[2];
sx q[2];
rz(-2.5131098) q[2];
rz(-2.8001522) q[3];
sx q[3];
rz(-2.1554558) q[3];
sx q[3];
rz(-2.3054874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8516561) q[0];
sx q[0];
rz(-2.4735232) q[0];
sx q[0];
rz(0.063902721) q[0];
rz(-0.54197657) q[1];
sx q[1];
rz(-1.130645) q[1];
sx q[1];
rz(0.37905395) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.097103216) q[0];
sx q[0];
rz(-1.3987204) q[0];
sx q[0];
rz(-1.8015566) q[0];
rz(-1.7070243) q[2];
sx q[2];
rz(-0.76644015) q[2];
sx q[2];
rz(2.5640783) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8833911) q[1];
sx q[1];
rz(-2.5978226) q[1];
sx q[1];
rz(0.030042458) q[1];
x q[2];
rz(-0.42551252) q[3];
sx q[3];
rz(-2.1149016) q[3];
sx q[3];
rz(-1.5489006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1698251) q[2];
sx q[2];
rz(-1.1173893) q[2];
sx q[2];
rz(-2.8525412) q[2];
rz(1.0158094) q[3];
sx q[3];
rz(-0.93520516) q[3];
sx q[3];
rz(-0.22739205) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0357901) q[0];
sx q[0];
rz(-0.38689026) q[0];
sx q[0];
rz(0.33045688) q[0];
rz(-0.48078787) q[1];
sx q[1];
rz(-1.5590706) q[1];
sx q[1];
rz(-0.50723433) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.071166) q[0];
sx q[0];
rz(-0.88001635) q[0];
sx q[0];
rz(0.69657495) q[0];
rz(-pi) q[1];
rz(1.9418409) q[2];
sx q[2];
rz(-2.0026653) q[2];
sx q[2];
rz(-0.88954207) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.005269231) q[1];
sx q[1];
rz(-1.988896) q[1];
sx q[1];
rz(0.24009248) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0546419) q[3];
sx q[3];
rz(-0.6020012) q[3];
sx q[3];
rz(-0.68893637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.5759739) q[2];
sx q[2];
rz(-2.2810292) q[2];
sx q[2];
rz(-0.38798517) q[2];
rz(3.0787789) q[3];
sx q[3];
rz(-0.7395491) q[3];
sx q[3];
rz(1.1886103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0906618) q[0];
sx q[0];
rz(-0.020314038) q[0];
sx q[0];
rz(0.055572979) q[0];
rz(-0.22124258) q[1];
sx q[1];
rz(-2.1317) q[1];
sx q[1];
rz(1.6411068) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7941754) q[0];
sx q[0];
rz(-2.5567434) q[0];
sx q[0];
rz(1.0561159) q[0];
rz(-pi) q[1];
rz(-0.13672853) q[2];
sx q[2];
rz(-1.5216344) q[2];
sx q[2];
rz(-1.8135742) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4909497) q[1];
sx q[1];
rz(-1.2526647) q[1];
sx q[1];
rz(-1.1699808) q[1];
rz(-0.59751038) q[3];
sx q[3];
rz(-2.198285) q[3];
sx q[3];
rz(1.6399872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1859493) q[2];
sx q[2];
rz(-2.9538437) q[2];
sx q[2];
rz(1.1090247) q[2];
rz(-3.1254613) q[3];
sx q[3];
rz(-1.707209) q[3];
sx q[3];
rz(2.6764638) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3817417) q[0];
sx q[0];
rz(-1.009059) q[0];
sx q[0];
rz(-0.41643634) q[0];
rz(-0.79553678) q[1];
sx q[1];
rz(-1.7557314) q[1];
sx q[1];
rz(3.0605127) q[1];
rz(3.0007985) q[2];
sx q[2];
rz(-1.9119605) q[2];
sx q[2];
rz(-0.20902363) q[2];
rz(2.1063188) q[3];
sx q[3];
rz(-2.2761619) q[3];
sx q[3];
rz(-0.63009562) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
