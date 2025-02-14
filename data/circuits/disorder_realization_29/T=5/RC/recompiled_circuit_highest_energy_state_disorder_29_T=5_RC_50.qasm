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
rz(-0.40773243) q[0];
sx q[0];
rz(4.3507504) q[0];
sx q[0];
rz(11.103295) q[0];
rz(-2.8830124) q[1];
sx q[1];
rz(-1.2376031) q[1];
sx q[1];
rz(-3.0078476) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8664854) q[0];
sx q[0];
rz(-1.5767528) q[0];
sx q[0];
rz(1.4136825) q[0];
x q[1];
rz(2.5630579) q[2];
sx q[2];
rz(-1.0867439) q[2];
sx q[2];
rz(-0.59147385) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.7514389) q[1];
sx q[1];
rz(-0.62733106) q[1];
sx q[1];
rz(-0.52951685) q[1];
rz(-pi) q[2];
rz(-3.1414019) q[3];
sx q[3];
rz(-1.4179363) q[3];
sx q[3];
rz(-0.78814414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.75970834) q[2];
sx q[2];
rz(-0.81177652) q[2];
sx q[2];
rz(-1.831057) q[2];
rz(-1.7740907) q[3];
sx q[3];
rz(-2.01367) q[3];
sx q[3];
rz(0.016157063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.1644208) q[0];
sx q[0];
rz(-1.1559957) q[0];
sx q[0];
rz(1.0154065) q[0];
rz(0.39335355) q[1];
sx q[1];
rz(-1.669408) q[1];
sx q[1];
rz(-0.70158395) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9727952) q[0];
sx q[0];
rz(-1.4864362) q[0];
sx q[0];
rz(2.286351) q[0];
rz(0.010741269) q[2];
sx q[2];
rz(-1.2678384) q[2];
sx q[2];
rz(0.078127351) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7483731) q[1];
sx q[1];
rz(-1.6459568) q[1];
sx q[1];
rz(1.2542226) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2982992) q[3];
sx q[3];
rz(-2.1873584) q[3];
sx q[3];
rz(-2.9697501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.027448805) q[2];
sx q[2];
rz(-1.0237209) q[2];
sx q[2];
rz(-2.0360428) q[2];
rz(-0.28087273) q[3];
sx q[3];
rz(-0.61087817) q[3];
sx q[3];
rz(-0.15225473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0356782) q[0];
sx q[0];
rz(-0.67258251) q[0];
sx q[0];
rz(2.0689082) q[0];
rz(-3.1283123) q[1];
sx q[1];
rz(-1.599879) q[1];
sx q[1];
rz(2.5650011) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7130797) q[0];
sx q[0];
rz(-1.858485) q[0];
sx q[0];
rz(1.5031394) q[0];
x q[1];
rz(2.1303612) q[2];
sx q[2];
rz(-0.15655993) q[2];
sx q[2];
rz(-2.8383534) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2943337) q[1];
sx q[1];
rz(-1.0953961) q[1];
sx q[1];
rz(-0.21902276) q[1];
x q[2];
rz(2.1284038) q[3];
sx q[3];
rz(-1.7806781) q[3];
sx q[3];
rz(2.5887973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.85792929) q[2];
sx q[2];
rz(-0.89202213) q[2];
sx q[2];
rz(2.499495) q[2];
rz(-0.56504956) q[3];
sx q[3];
rz(-0.27500209) q[3];
sx q[3];
rz(2.9716085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86244407) q[0];
sx q[0];
rz(-1.2553517) q[0];
sx q[0];
rz(-2.8992262) q[0];
rz(0.53388059) q[1];
sx q[1];
rz(-2.4568074) q[1];
sx q[1];
rz(-1.6530316) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80099312) q[0];
sx q[0];
rz(-1.6417608) q[0];
sx q[0];
rz(0.59619716) q[0];
rz(-pi) q[1];
rz(-0.42885355) q[2];
sx q[2];
rz(-1.2786713) q[2];
sx q[2];
rz(0.090306239) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.10746126) q[1];
sx q[1];
rz(-0.80781579) q[1];
sx q[1];
rz(-0.40523578) q[1];
x q[2];
rz(-0.29976254) q[3];
sx q[3];
rz(-1.5325755) q[3];
sx q[3];
rz(2.770407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7786467) q[2];
sx q[2];
rz(-1.2641509) q[2];
sx q[2];
rz(1.6710336) q[2];
rz(-0.10739022) q[3];
sx q[3];
rz(-1.8348179) q[3];
sx q[3];
rz(-1.2764527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.030468) q[0];
sx q[0];
rz(-2.680439) q[0];
sx q[0];
rz(1.9386468) q[0];
rz(0.89490926) q[1];
sx q[1];
rz(-1.1824965) q[1];
sx q[1];
rz(-2.9295909) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1231736) q[0];
sx q[0];
rz(-0.97476649) q[0];
sx q[0];
rz(-0.53431781) q[0];
rz(-1.2707527) q[2];
sx q[2];
rz(-2.0738389) q[2];
sx q[2];
rz(-1.4936269) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.743269) q[1];
sx q[1];
rz(-1.124427) q[1];
sx q[1];
rz(-0.94029398) q[1];
rz(-pi) q[2];
rz(1.5863442) q[3];
sx q[3];
rz(-0.55335303) q[3];
sx q[3];
rz(-1.4912332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8733069) q[2];
sx q[2];
rz(-2.3640859) q[2];
sx q[2];
rz(-3.0262465) q[2];
rz(-2.1558732) q[3];
sx q[3];
rz(-1.7306354) q[3];
sx q[3];
rz(-2.7058097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2301521) q[0];
sx q[0];
rz(-1.6807012) q[0];
sx q[0];
rz(0.69754115) q[0];
rz(-0.41360924) q[1];
sx q[1];
rz(-1.6802843) q[1];
sx q[1];
rz(-2.4381309) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0516104) q[0];
sx q[0];
rz(-1.713366) q[0];
sx q[0];
rz(2.8683314) q[0];
rz(-0.38833924) q[2];
sx q[2];
rz(-2.1415798) q[2];
sx q[2];
rz(-2.6663189) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5825962) q[1];
sx q[1];
rz(-1.1038053) q[1];
sx q[1];
rz(-1.9860616) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2715152) q[3];
sx q[3];
rz(-1.8095152) q[3];
sx q[3];
rz(-1.6629926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0387663) q[2];
sx q[2];
rz(-2.8732754) q[2];
sx q[2];
rz(1.212567) q[2];
rz(-2.8596527) q[3];
sx q[3];
rz(-1.599267) q[3];
sx q[3];
rz(-0.50659242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3459699) q[0];
sx q[0];
rz(-2.3891734) q[0];
sx q[0];
rz(2.2947626) q[0];
rz(-2.1961191) q[1];
sx q[1];
rz(-1.5676326) q[1];
sx q[1];
rz(0.6792773) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4436809) q[0];
sx q[0];
rz(-0.4831995) q[0];
sx q[0];
rz(-0.16541055) q[0];
rz(-pi) q[1];
x q[1];
rz(0.34946673) q[2];
sx q[2];
rz(-2.4878056) q[2];
sx q[2];
rz(0.9330627) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.32545391) q[1];
sx q[1];
rz(-1.7569764) q[1];
sx q[1];
rz(2.0912311) q[1];
rz(-pi) q[2];
rz(-0.019370989) q[3];
sx q[3];
rz(-1.8550411) q[3];
sx q[3];
rz(1.5694048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.88963738) q[2];
sx q[2];
rz(-2.1869982) q[2];
sx q[2];
rz(-1.6471222) q[2];
rz(0.55388081) q[3];
sx q[3];
rz(-1.5364105) q[3];
sx q[3];
rz(1.4166098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3225591) q[0];
sx q[0];
rz(-2.366134) q[0];
sx q[0];
rz(-2.0490647) q[0];
rz(-0.93387261) q[1];
sx q[1];
rz(-1.4051508) q[1];
sx q[1];
rz(-2.3563103) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0781495) q[0];
sx q[0];
rz(-2.5012996) q[0];
sx q[0];
rz(2.3904353) q[0];
x q[1];
rz(1.0243591) q[2];
sx q[2];
rz(-1.1866202) q[2];
sx q[2];
rz(2.667493) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.685647) q[1];
sx q[1];
rz(-0.37044493) q[1];
sx q[1];
rz(0.73729314) q[1];
rz(2.0075032) q[3];
sx q[3];
rz(-0.99815449) q[3];
sx q[3];
rz(-2.9331895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2698722) q[2];
sx q[2];
rz(-1.4464658) q[2];
sx q[2];
rz(-2.5093057) q[2];
rz(2.2271683) q[3];
sx q[3];
rz(-1.881003) q[3];
sx q[3];
rz(0.15973346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3226427) q[0];
sx q[0];
rz(-0.85700789) q[0];
sx q[0];
rz(0.44347611) q[0];
rz(-0.30287287) q[1];
sx q[1];
rz(-0.58285204) q[1];
sx q[1];
rz(1.3076967) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46638495) q[0];
sx q[0];
rz(-1.2699915) q[0];
sx q[0];
rz(-1.3928901) q[0];
rz(1.8820178) q[2];
sx q[2];
rz(-1.0838795) q[2];
sx q[2];
rz(2.4608909) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2949867) q[1];
sx q[1];
rz(-1.8364803) q[1];
sx q[1];
rz(-2.0779266) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0293209) q[3];
sx q[3];
rz(-0.53526894) q[3];
sx q[3];
rz(-0.9891555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0206454) q[2];
sx q[2];
rz(-2.1006613) q[2];
sx q[2];
rz(-0.61332235) q[2];
rz(0.81765085) q[3];
sx q[3];
rz(-2.652467) q[3];
sx q[3];
rz(2.8221655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15449512) q[0];
sx q[0];
rz(-1.7648062) q[0];
sx q[0];
rz(-2.7823271) q[0];
rz(2.477395) q[1];
sx q[1];
rz(-1.8134873) q[1];
sx q[1];
rz(0.42116234) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65972787) q[0];
sx q[0];
rz(-2.7057428) q[0];
sx q[0];
rz(-2.1207366) q[0];
rz(3.1362651) q[2];
sx q[2];
rz(-1.6432646) q[2];
sx q[2];
rz(-0.42124149) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5790184) q[1];
sx q[1];
rz(-1.1179325) q[1];
sx q[1];
rz(0.29154204) q[1];
x q[2];
rz(-2.3844658) q[3];
sx q[3];
rz(-2.3806678) q[3];
sx q[3];
rz(-1.8426551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.26587129) q[2];
sx q[2];
rz(-3.044812) q[2];
sx q[2];
rz(1.9270012) q[2];
rz(-2.752979) q[3];
sx q[3];
rz(-0.92258421) q[3];
sx q[3];
rz(1.0987561) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6471967) q[0];
sx q[0];
rz(-1.5955441) q[0];
sx q[0];
rz(1.2936976) q[0];
rz(-2.9248059) q[1];
sx q[1];
rz(-2.4516791) q[1];
sx q[1];
rz(-2.6286415) q[1];
rz(2.6224455) q[2];
sx q[2];
rz(-2.3943974) q[2];
sx q[2];
rz(1.129528) q[2];
rz(2.67542) q[3];
sx q[3];
rz(-1.331658) q[3];
sx q[3];
rz(-0.030203947) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
