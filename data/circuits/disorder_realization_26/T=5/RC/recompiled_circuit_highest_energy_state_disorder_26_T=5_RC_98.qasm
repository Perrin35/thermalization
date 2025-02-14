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
rz(0.17231365) q[0];
sx q[0];
rz(-2.9336689) q[0];
sx q[0];
rz(1.2907668) q[0];
rz(-2.9895904) q[1];
sx q[1];
rz(-0.64270371) q[1];
sx q[1];
rz(-0.42833498) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5458459) q[0];
sx q[0];
rz(-1.4662678) q[0];
sx q[0];
rz(-1.4609817) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8646127) q[2];
sx q[2];
rz(-0.91962469) q[2];
sx q[2];
rz(3.0319954) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0103156) q[1];
sx q[1];
rz(-2.4020264) q[1];
sx q[1];
rz(0.063994813) q[1];
rz(-0.21875225) q[3];
sx q[3];
rz(-1.7667337) q[3];
sx q[3];
rz(-2.4234555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.20994818) q[2];
sx q[2];
rz(-2.1817709) q[2];
sx q[2];
rz(1.2098562) q[2];
rz(0.079553902) q[3];
sx q[3];
rz(-0.39593655) q[3];
sx q[3];
rz(-2.615926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2570268) q[0];
sx q[0];
rz(-0.50978065) q[0];
sx q[0];
rz(2.7575745) q[0];
rz(2.2513023) q[1];
sx q[1];
rz(-1.9046116) q[1];
sx q[1];
rz(0.30337897) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7913032) q[0];
sx q[0];
rz(-1.3801738) q[0];
sx q[0];
rz(2.1146569) q[0];
x q[1];
rz(-2.882134) q[2];
sx q[2];
rz(-2.0739809) q[2];
sx q[2];
rz(0.34215701) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1282004) q[1];
sx q[1];
rz(-1.6329995) q[1];
sx q[1];
rz(1.661646) q[1];
rz(-pi) q[2];
rz(1.1214662) q[3];
sx q[3];
rz(-0.96689618) q[3];
sx q[3];
rz(2.3788798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7552135) q[2];
sx q[2];
rz(-1.1073802) q[2];
sx q[2];
rz(-2.4066822) q[2];
rz(-0.84878659) q[3];
sx q[3];
rz(-2.3990192) q[3];
sx q[3];
rz(-2.7809704) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81247771) q[0];
sx q[0];
rz(-0.37549967) q[0];
sx q[0];
rz(-3.062881) q[0];
rz(-0.82376897) q[1];
sx q[1];
rz(-1.0853826) q[1];
sx q[1];
rz(-1.2944006) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20939182) q[0];
sx q[0];
rz(-0.086361445) q[0];
sx q[0];
rz(-2.3461653) q[0];
x q[1];
rz(2.2357777) q[2];
sx q[2];
rz(-2.40529) q[2];
sx q[2];
rz(1.5018963) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.489276) q[1];
sx q[1];
rz(-1.4536263) q[1];
sx q[1];
rz(2.8588309) q[1];
rz(-0.23842509) q[3];
sx q[3];
rz(-2.2608071) q[3];
sx q[3];
rz(2.391614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2417629) q[2];
sx q[2];
rz(-0.37909847) q[2];
sx q[2];
rz(-2.3251593) q[2];
rz(-1.5166616) q[3];
sx q[3];
rz(-2.1818325) q[3];
sx q[3];
rz(2.4412156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65275943) q[0];
sx q[0];
rz(-2.3214898) q[0];
sx q[0];
rz(0.56831992) q[0];
rz(1.9991416) q[1];
sx q[1];
rz(-1.7894952) q[1];
sx q[1];
rz(-1.4521339) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39740085) q[0];
sx q[0];
rz(-2.4699587) q[0];
sx q[0];
rz(-0.36104843) q[0];
rz(-pi) q[1];
rz(-2.6645489) q[2];
sx q[2];
rz(-0.20186587) q[2];
sx q[2];
rz(-2.9143726) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.3209705) q[1];
sx q[1];
rz(-1.8838716) q[1];
sx q[1];
rz(-2.2105107) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9098572) q[3];
sx q[3];
rz(-1.4022458) q[3];
sx q[3];
rz(1.7730561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.51231724) q[2];
sx q[2];
rz(-2.6111111) q[2];
sx q[2];
rz(-0.076889195) q[2];
rz(-0.39938375) q[3];
sx q[3];
rz(-2.2279584) q[3];
sx q[3];
rz(0.54401773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9870616) q[0];
sx q[0];
rz(-2.7847544) q[0];
sx q[0];
rz(2.8416908) q[0];
rz(1.9918282) q[1];
sx q[1];
rz(-1.0118142) q[1];
sx q[1];
rz(0.48847517) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91575275) q[0];
sx q[0];
rz(-1.7560648) q[0];
sx q[0];
rz(-0.033323296) q[0];
rz(2.4534493) q[2];
sx q[2];
rz(-2.4240085) q[2];
sx q[2];
rz(-0.46844278) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.27389535) q[1];
sx q[1];
rz(-0.58389837) q[1];
sx q[1];
rz(-0.03429596) q[1];
rz(-pi) q[2];
rz(-1.8397055) q[3];
sx q[3];
rz(-2.1950414) q[3];
sx q[3];
rz(-1.3866803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.38146314) q[2];
sx q[2];
rz(-0.13801408) q[2];
sx q[2];
rz(1.6628954) q[2];
rz(1.1100769) q[3];
sx q[3];
rz(-0.9980945) q[3];
sx q[3];
rz(2.4983675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7396624) q[0];
sx q[0];
rz(-1.1818385) q[0];
sx q[0];
rz(-2.8231743) q[0];
rz(-0.034612522) q[1];
sx q[1];
rz(-0.56762677) q[1];
sx q[1];
rz(-2.6908223) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0885267) q[0];
sx q[0];
rz(-0.51061106) q[0];
sx q[0];
rz(-1.5998163) q[0];
x q[1];
rz(-0.87574739) q[2];
sx q[2];
rz(-1.0330457) q[2];
sx q[2];
rz(-1.328383) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.83799441) q[1];
sx q[1];
rz(-0.7571836) q[1];
sx q[1];
rz(-0.96214575) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9616021) q[3];
sx q[3];
rz(-2.3473661) q[3];
sx q[3];
rz(-3.1193747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0938809) q[2];
sx q[2];
rz(-0.1897976) q[2];
sx q[2];
rz(0.06165687) q[2];
rz(0.098585248) q[3];
sx q[3];
rz(-2.3969789) q[3];
sx q[3];
rz(-0.70257598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3848569) q[0];
sx q[0];
rz(-2.0282133) q[0];
sx q[0];
rz(0.039948832) q[0];
rz(-2.3710251) q[1];
sx q[1];
rz(-0.71208411) q[1];
sx q[1];
rz(-0.22824731) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21439274) q[0];
sx q[0];
rz(-1.6480371) q[0];
sx q[0];
rz(-1.6595675) q[0];
x q[1];
rz(-1.8769774) q[2];
sx q[2];
rz(-1.5723229) q[2];
sx q[2];
rz(-0.36097872) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.952164) q[1];
sx q[1];
rz(-1.4888114) q[1];
sx q[1];
rz(1.6240955) q[1];
rz(-pi) q[2];
rz(2.1114717) q[3];
sx q[3];
rz(-0.70899963) q[3];
sx q[3];
rz(-1.0741155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0780636) q[2];
sx q[2];
rz(-1.0588131) q[2];
sx q[2];
rz(0.25827363) q[2];
rz(2.9410948) q[3];
sx q[3];
rz(-0.79810464) q[3];
sx q[3];
rz(2.958278) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9789326) q[0];
sx q[0];
rz(-1.7698092) q[0];
sx q[0];
rz(-0.3072511) q[0];
rz(1.2112674) q[1];
sx q[1];
rz(-2.6478421) q[1];
sx q[1];
rz(0.58120751) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4330313) q[0];
sx q[0];
rz(-1.1662467) q[0];
sx q[0];
rz(-0.033971196) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9517854) q[2];
sx q[2];
rz(-2.4229089) q[2];
sx q[2];
rz(-2.725012) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8381184) q[1];
sx q[1];
rz(-1.3708769) q[1];
sx q[1];
rz(-2.3568627) q[1];
rz(-pi) q[2];
rz(0.56116207) q[3];
sx q[3];
rz(-0.61120874) q[3];
sx q[3];
rz(-1.2964013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6792949) q[2];
sx q[2];
rz(-2.0079948) q[2];
sx q[2];
rz(0.031166859) q[2];
rz(2.9160685) q[3];
sx q[3];
rz(-1.3616819) q[3];
sx q[3];
rz(-1.0303191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0533326) q[0];
sx q[0];
rz(-3.0971165) q[0];
sx q[0];
rz(0.70892507) q[0];
rz(-0.21253474) q[1];
sx q[1];
rz(-2.1268851) q[1];
sx q[1];
rz(-2.2841891) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.014985059) q[0];
sx q[0];
rz(-1.7005655) q[0];
sx q[0];
rz(3.1364822) q[0];
rz(-pi) q[1];
rz(2.6231758) q[2];
sx q[2];
rz(-2.5441493) q[2];
sx q[2];
rz(-0.13472508) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.58536968) q[1];
sx q[1];
rz(-1.1818647) q[1];
sx q[1];
rz(-0.57805581) q[1];
rz(-pi) q[2];
rz(2.4444432) q[3];
sx q[3];
rz(-1.4065885) q[3];
sx q[3];
rz(-2.177161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4400441) q[2];
sx q[2];
rz(-2.7893119) q[2];
sx q[2];
rz(0.80553833) q[2];
rz(2.7696179) q[3];
sx q[3];
rz(-1.6476846) q[3];
sx q[3];
rz(-0.39723799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.032967903) q[0];
sx q[0];
rz(-1.911835) q[0];
sx q[0];
rz(0.97880542) q[0];
rz(-0.72273123) q[1];
sx q[1];
rz(-2.0131854) q[1];
sx q[1];
rz(-0.61789787) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5586097) q[0];
sx q[0];
rz(-1.8454058) q[0];
sx q[0];
rz(2.3586169) q[0];
rz(-pi) q[1];
rz(-2.5064836) q[2];
sx q[2];
rz(-2.375787) q[2];
sx q[2];
rz(2.2857411) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7405211) q[1];
sx q[1];
rz(-2.8943672) q[1];
sx q[1];
rz(-0.81584357) q[1];
rz(2.3506451) q[3];
sx q[3];
rz(-3.0229048) q[3];
sx q[3];
rz(-3.0968248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.939398) q[2];
sx q[2];
rz(-1.3556182) q[2];
sx q[2];
rz(0.00016577684) q[2];
rz(0.58445066) q[3];
sx q[3];
rz(-2.1355459) q[3];
sx q[3];
rz(0.59529006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38681876) q[0];
sx q[0];
rz(-1.5703572) q[0];
sx q[0];
rz(1.5686709) q[0];
rz(1.8008925) q[1];
sx q[1];
rz(-1.1005713) q[1];
sx q[1];
rz(1.5250199) q[1];
rz(-0.98967057) q[2];
sx q[2];
rz(-0.72643092) q[2];
sx q[2];
rz(0.55813172) q[2];
rz(-1.9289005) q[3];
sx q[3];
rz(-1.9869589) q[3];
sx q[3];
rz(-0.12369894) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
