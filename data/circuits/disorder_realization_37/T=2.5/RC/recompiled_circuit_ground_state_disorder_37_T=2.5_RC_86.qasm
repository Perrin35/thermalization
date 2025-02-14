OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3545544) q[0];
sx q[0];
rz(0.7631425) q[0];
sx q[0];
rz(7.2416303) q[0];
rz(-2.7148442) q[1];
sx q[1];
rz(3.5659748) q[1];
sx q[1];
rz(10.397484) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1866731) q[0];
sx q[0];
rz(-1.3698319) q[0];
sx q[0];
rz(-1.475564) q[0];
rz(-pi) q[1];
rz(-0.17225687) q[2];
sx q[2];
rz(-1.5035727) q[2];
sx q[2];
rz(0.40861118) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8805566) q[1];
sx q[1];
rz(-0.86819217) q[1];
sx q[1];
rz(-1.2289341) q[1];
rz(-pi) q[2];
rz(-1.7250502) q[3];
sx q[3];
rz(-1.8498382) q[3];
sx q[3];
rz(-1.0629625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.080231754) q[2];
sx q[2];
rz(-1.1876567) q[2];
sx q[2];
rz(2.3415671) q[2];
rz(0.13036615) q[3];
sx q[3];
rz(-0.76821199) q[3];
sx q[3];
rz(-0.31726328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.80113634) q[0];
sx q[0];
rz(-1.2202593) q[0];
sx q[0];
rz(-2.1521547) q[0];
rz(-0.89775741) q[1];
sx q[1];
rz(-1.0549301) q[1];
sx q[1];
rz(-2.3602233) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5637276) q[0];
sx q[0];
rz(-1.0672309) q[0];
sx q[0];
rz(0.2204075) q[0];
rz(0.23405646) q[2];
sx q[2];
rz(-1.4849201) q[2];
sx q[2];
rz(1.1669005) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9354089) q[1];
sx q[1];
rz(-2.0061135) q[1];
sx q[1];
rz(-1.0401506) q[1];
x q[2];
rz(0.3831634) q[3];
sx q[3];
rz(-1.8472611) q[3];
sx q[3];
rz(-0.081646669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.117131) q[2];
sx q[2];
rz(-1.3776366) q[2];
sx q[2];
rz(-0.77829877) q[2];
rz(0.69027573) q[3];
sx q[3];
rz(-0.92167753) q[3];
sx q[3];
rz(-1.5811623) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86163259) q[0];
sx q[0];
rz(-2.3065688) q[0];
sx q[0];
rz(0.41197187) q[0];
rz(2.1906134) q[1];
sx q[1];
rz(-0.6528267) q[1];
sx q[1];
rz(-1.2118118) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6854343) q[0];
sx q[0];
rz(-3.0904909) q[0];
sx q[0];
rz(-0.85275485) q[0];
rz(1.3479718) q[2];
sx q[2];
rz(-1.4744024) q[2];
sx q[2];
rz(-2.4501462) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.40719051) q[1];
sx q[1];
rz(-2.4311594) q[1];
sx q[1];
rz(-3.0991952) q[1];
x q[2];
rz(0.36529593) q[3];
sx q[3];
rz(-1.4202227) q[3];
sx q[3];
rz(3.0677861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1872306) q[2];
sx q[2];
rz(-0.53202859) q[2];
sx q[2];
rz(1.5032035) q[2];
rz(3.1014118) q[3];
sx q[3];
rz(-0.92622042) q[3];
sx q[3];
rz(0.23475501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16315854) q[0];
sx q[0];
rz(-1.5143159) q[0];
sx q[0];
rz(3.0117595) q[0];
rz(1.8561329) q[1];
sx q[1];
rz(-1.9655656) q[1];
sx q[1];
rz(1.0349549) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.267535) q[0];
sx q[0];
rz(-2.3685936) q[0];
sx q[0];
rz(-2.2202047) q[0];
rz(-3.0327243) q[2];
sx q[2];
rz(-2.2316859) q[2];
sx q[2];
rz(-0.79273293) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.538547) q[1];
sx q[1];
rz(-1.6563452) q[1];
sx q[1];
rz(-0.24876068) q[1];
x q[2];
rz(1.0134936) q[3];
sx q[3];
rz(-2.3689007) q[3];
sx q[3];
rz(1.1232291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3639565) q[2];
sx q[2];
rz(-1.3866321) q[2];
sx q[2];
rz(3.0470972) q[2];
rz(2.7206521) q[3];
sx q[3];
rz(-1.4237483) q[3];
sx q[3];
rz(-1.4360992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41973758) q[0];
sx q[0];
rz(-2.5696745) q[0];
sx q[0];
rz(-0.016059248) q[0];
rz(2.2968966) q[1];
sx q[1];
rz(-0.41929308) q[1];
sx q[1];
rz(0.90763456) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1354271) q[0];
sx q[0];
rz(-2.3434397) q[0];
sx q[0];
rz(0.31240518) q[0];
x q[1];
rz(0.24036562) q[2];
sx q[2];
rz(-1.1615001) q[2];
sx q[2];
rz(1.5549007) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6542977) q[1];
sx q[1];
rz(-1.134344) q[1];
sx q[1];
rz(-0.36466332) q[1];
rz(-pi) q[2];
rz(-1.5441555) q[3];
sx q[3];
rz(-1.2268042) q[3];
sx q[3];
rz(-0.87061239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.57546651) q[2];
sx q[2];
rz(-1.0593654) q[2];
sx q[2];
rz(0.063684138) q[2];
rz(0.56695402) q[3];
sx q[3];
rz(-0.83438116) q[3];
sx q[3];
rz(-2.5440192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31396922) q[0];
sx q[0];
rz(-2.9997928) q[0];
sx q[0];
rz(-1.1578479) q[0];
rz(1.6572378) q[1];
sx q[1];
rz(-1.4446222) q[1];
sx q[1];
rz(-0.2690014) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3753178) q[0];
sx q[0];
rz(-2.9937389) q[0];
sx q[0];
rz(0.30884858) q[0];
rz(3.0458955) q[2];
sx q[2];
rz(-1.2996593) q[2];
sx q[2];
rz(-2.1094131) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.77226247) q[1];
sx q[1];
rz(-0.20552615) q[1];
sx q[1];
rz(1.0632681) q[1];
rz(-0.46697843) q[3];
sx q[3];
rz(-0.45736936) q[3];
sx q[3];
rz(0.46588313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.77702648) q[2];
sx q[2];
rz(-0.51743162) q[2];
sx q[2];
rz(0.88097921) q[2];
rz(3.0173054) q[3];
sx q[3];
rz(-1.7139939) q[3];
sx q[3];
rz(-2.4833) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2800901) q[0];
sx q[0];
rz(-2.8741591) q[0];
sx q[0];
rz(-0.68914831) q[0];
rz(2.6742477) q[1];
sx q[1];
rz(-0.57599774) q[1];
sx q[1];
rz(2.0435832) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3680237) q[0];
sx q[0];
rz(-1.7360285) q[0];
sx q[0];
rz(1.6501897) q[0];
rz(-pi) q[1];
rz(-1.271209) q[2];
sx q[2];
rz(-2.0986084) q[2];
sx q[2];
rz(2.6115548) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4773405) q[1];
sx q[1];
rz(-2.1108529) q[1];
sx q[1];
rz(-1.6083205) q[1];
x q[2];
rz(-2.6379273) q[3];
sx q[3];
rz(-2.1379469) q[3];
sx q[3];
rz(0.43225542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5757307) q[2];
sx q[2];
rz(-0.87811676) q[2];
sx q[2];
rz(2.511054) q[2];
rz(2.5343043) q[3];
sx q[3];
rz(-2.2346965) q[3];
sx q[3];
rz(-1.0926532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31834114) q[0];
sx q[0];
rz(-2.9956151) q[0];
sx q[0];
rz(-1.3463705) q[0];
rz(1.3660376) q[1];
sx q[1];
rz(-0.64907688) q[1];
sx q[1];
rz(0.62873658) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7766984) q[0];
sx q[0];
rz(-1.8275675) q[0];
sx q[0];
rz(2.5904694) q[0];
rz(-pi) q[1];
rz(0.89968483) q[2];
sx q[2];
rz(-1.1555399) q[2];
sx q[2];
rz(-2.496831) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3676074) q[1];
sx q[1];
rz(-2.1091568) q[1];
sx q[1];
rz(2.143285) q[1];
rz(-pi) q[2];
rz(0.076185779) q[3];
sx q[3];
rz(-1.3722808) q[3];
sx q[3];
rz(3.1358058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2748572) q[2];
sx q[2];
rz(-1.3254415) q[2];
sx q[2];
rz(-2.2692915) q[2];
rz(-2.9288779) q[3];
sx q[3];
rz(-1.745696) q[3];
sx q[3];
rz(-0.58342903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23671959) q[0];
sx q[0];
rz(-1.6107591) q[0];
sx q[0];
rz(-1.8778296) q[0];
rz(1.0172552) q[1];
sx q[1];
rz(-2.2433498) q[1];
sx q[1];
rz(0.57473007) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.547366) q[0];
sx q[0];
rz(-1.5061139) q[0];
sx q[0];
rz(0.5314188) q[0];
rz(-pi) q[1];
rz(-0.2391441) q[2];
sx q[2];
rz(-2.2986293) q[2];
sx q[2];
rz(0.69918888) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8338477) q[1];
sx q[1];
rz(-2.7258432) q[1];
sx q[1];
rz(-0.47885696) q[1];
rz(-0.25463661) q[3];
sx q[3];
rz(-2.0395425) q[3];
sx q[3];
rz(0.96118467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.42402521) q[2];
sx q[2];
rz(-1.8342476) q[2];
sx q[2];
rz(-0.032111017) q[2];
rz(1.1561681) q[3];
sx q[3];
rz(-0.38438946) q[3];
sx q[3];
rz(-1.9305852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1011937) q[0];
sx q[0];
rz(-2.5960584) q[0];
sx q[0];
rz(1.1241166) q[0];
rz(-0.54565564) q[1];
sx q[1];
rz(-0.35174313) q[1];
sx q[1];
rz(0.55145946) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0113946) q[0];
sx q[0];
rz(-1.4470248) q[0];
sx q[0];
rz(1.3817203) q[0];
rz(-pi) q[1];
rz(-1.8449502) q[2];
sx q[2];
rz(-1.5197496) q[2];
sx q[2];
rz(2.5921043) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3797326) q[1];
sx q[1];
rz(-1.1071536) q[1];
sx q[1];
rz(-2.597888) q[1];
rz(-pi) q[2];
rz(2.6607108) q[3];
sx q[3];
rz(-2.7989498) q[3];
sx q[3];
rz(-0.51183701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.55190279) q[2];
sx q[2];
rz(-2.6305113) q[2];
sx q[2];
rz(1.4116633) q[2];
rz(0.74665135) q[3];
sx q[3];
rz(-1.6801497) q[3];
sx q[3];
rz(2.1667229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-0.31155561) q[0];
sx q[0];
rz(-1.5811601) q[0];
sx q[0];
rz(1.5720221) q[0];
rz(-0.39881067) q[1];
sx q[1];
rz(-1.8067982) q[1];
sx q[1];
rz(-1.6511818) q[1];
rz(3.1337784) q[2];
sx q[2];
rz(-1.0823156) q[2];
sx q[2];
rz(-1.2132614) q[2];
rz(-0.79979782) q[3];
sx q[3];
rz(-1.3577438) q[3];
sx q[3];
rz(1.8997026) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
