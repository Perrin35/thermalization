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
rz(2.7338602) q[0];
sx q[0];
rz(-1.2091577) q[0];
sx q[0];
rz(-1.6785167) q[0];
rz(-2.8830124) q[1];
sx q[1];
rz(-1.2376031) q[1];
sx q[1];
rz(-3.0078476) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2966327) q[0];
sx q[0];
rz(-1.4136853) q[0];
sx q[0];
rz(-3.1355619) q[0];
rz(-pi) q[1];
rz(2.5630579) q[2];
sx q[2];
rz(-2.0548487) q[2];
sx q[2];
rz(-2.5501188) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.37686302) q[1];
sx q[1];
rz(-1.8718182) q[1];
sx q[1];
rz(-2.5824598) q[1];
rz(-pi) q[2];
rz(-1.5720346) q[3];
sx q[3];
rz(-0.15286013) q[3];
sx q[3];
rz(2.3521956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.75970834) q[2];
sx q[2];
rz(-0.81177652) q[2];
sx q[2];
rz(-1.3105357) q[2];
rz(1.367502) q[3];
sx q[3];
rz(-2.01367) q[3];
sx q[3];
rz(0.016157063) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1644208) q[0];
sx q[0];
rz(-1.1559957) q[0];
sx q[0];
rz(-2.1261862) q[0];
rz(-2.7482391) q[1];
sx q[1];
rz(-1.669408) q[1];
sx q[1];
rz(2.4400087) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4751101) q[0];
sx q[0];
rz(-0.85832867) q[0];
sx q[0];
rz(-0.11157596) q[0];
rz(1.6051454) q[2];
sx q[2];
rz(-2.8384502) q[2];
sx q[2];
rz(0.042138635) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.39321956) q[1];
sx q[1];
rz(-1.4956359) q[1];
sx q[1];
rz(1.88737) q[1];
x q[2];
rz(0.84329347) q[3];
sx q[3];
rz(-2.1873584) q[3];
sx q[3];
rz(-0.17184251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.027448805) q[2];
sx q[2];
rz(-1.0237209) q[2];
sx q[2];
rz(2.0360428) q[2];
rz(-2.8607199) q[3];
sx q[3];
rz(-2.5307145) q[3];
sx q[3];
rz(2.9893379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1059145) q[0];
sx q[0];
rz(-0.67258251) q[0];
sx q[0];
rz(-1.0726844) q[0];
rz(-3.1283123) q[1];
sx q[1];
rz(-1.5417136) q[1];
sx q[1];
rz(0.57659155) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4786561) q[0];
sx q[0];
rz(-0.29532224) q[0];
sx q[0];
rz(0.22462019) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0835952) q[2];
sx q[2];
rz(-1.4382677) q[2];
sx q[2];
rz(-0.26187632) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.84725897) q[1];
sx q[1];
rz(-2.0461966) q[1];
sx q[1];
rz(2.9225699) q[1];
rz(-pi) q[2];
rz(1.9535096) q[3];
sx q[3];
rz(-2.5497304) q[3];
sx q[3];
rz(-2.446021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.85792929) q[2];
sx q[2];
rz(-2.2495705) q[2];
sx q[2];
rz(-0.64209765) q[2];
rz(-0.56504956) q[3];
sx q[3];
rz(-2.8665906) q[3];
sx q[3];
rz(-2.9716085) q[3];
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
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2791486) q[0];
sx q[0];
rz(-1.886241) q[0];
sx q[0];
rz(-0.24236648) q[0];
rz(-2.6077121) q[1];
sx q[1];
rz(-2.4568074) q[1];
sx q[1];
rz(-1.6530316) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87391728) q[0];
sx q[0];
rz(-2.5416964) q[0];
sx q[0];
rz(3.0156662) q[0];
rz(0.62612957) q[2];
sx q[2];
rz(-0.51373791) q[2];
sx q[2];
rz(-1.0990248) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0341314) q[1];
sx q[1];
rz(-0.80781579) q[1];
sx q[1];
rz(0.40523578) q[1];
rz(-pi) q[2];
rz(-1.6107992) q[3];
sx q[3];
rz(-1.8703331) q[3];
sx q[3];
rz(-1.9301723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.36294595) q[2];
sx q[2];
rz(-1.2641509) q[2];
sx q[2];
rz(-1.470559) q[2];
rz(3.0342024) q[3];
sx q[3];
rz(-1.3067747) q[3];
sx q[3];
rz(-1.86514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.030468) q[0];
sx q[0];
rz(-0.46115369) q[0];
sx q[0];
rz(1.2029458) q[0];
rz(0.89490926) q[1];
sx q[1];
rz(-1.9590961) q[1];
sx q[1];
rz(2.9295909) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0184191) q[0];
sx q[0];
rz(-0.97476649) q[0];
sx q[0];
rz(2.6072748) q[0];
x q[1];
rz(-1.8708399) q[2];
sx q[2];
rz(-2.0738389) q[2];
sx q[2];
rz(-1.6479657) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.70658503) q[1];
sx q[1];
rz(-2.3870578) q[1];
sx q[1];
rz(2.2526787) q[1];
rz(-pi) q[2];
x q[2];
rz(3.131989) q[3];
sx q[3];
rz(-1.017518) q[3];
sx q[3];
rz(1.6320848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2682858) q[2];
sx q[2];
rz(-2.3640859) q[2];
sx q[2];
rz(-3.0262465) q[2];
rz(-2.1558732) q[3];
sx q[3];
rz(-1.4109572) q[3];
sx q[3];
rz(2.7058097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91144052) q[0];
sx q[0];
rz(-1.6807012) q[0];
sx q[0];
rz(0.69754115) q[0];
rz(-2.7279834) q[1];
sx q[1];
rz(-1.6802843) q[1];
sx q[1];
rz(2.4381309) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.441012) q[0];
sx q[0];
rz(-1.3003775) q[0];
sx q[0];
rz(-1.4228113) q[0];
x q[1];
rz(-1.037961) q[2];
sx q[2];
rz(-0.67801266) q[2];
sx q[2];
rz(0.17282669) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5589964) q[1];
sx q[1];
rz(-2.0377874) q[1];
sx q[1];
rz(-1.9860616) q[1];
x q[2];
rz(-0.88075068) q[3];
sx q[3];
rz(-2.7610169) q[3];
sx q[3];
rz(-0.56169034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0387663) q[2];
sx q[2];
rz(-2.8732754) q[2];
sx q[2];
rz(-1.9290257) q[2];
rz(2.8596527) q[3];
sx q[3];
rz(-1.599267) q[3];
sx q[3];
rz(0.50659242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7956227) q[0];
sx q[0];
rz(-0.75241929) q[0];
sx q[0];
rz(-0.84683007) q[0];
rz(2.1961191) q[1];
sx q[1];
rz(-1.5676326) q[1];
sx q[1];
rz(-0.6792773) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4436809) q[0];
sx q[0];
rz(-2.6583932) q[0];
sx q[0];
rz(-0.16541055) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7921259) q[2];
sx q[2];
rz(-0.65378705) q[2];
sx q[2];
rz(0.9330627) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5577364) q[1];
sx q[1];
rz(-2.5917555) q[1];
sx q[1];
rz(-1.932895) q[1];
rz(-pi) q[2];
x q[2];
rz(0.019370989) q[3];
sx q[3];
rz(-1.2865515) q[3];
sx q[3];
rz(-1.5721878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.88963738) q[2];
sx q[2];
rz(-0.95459443) q[2];
sx q[2];
rz(-1.6471222) q[2];
rz(-2.5877118) q[3];
sx q[3];
rz(-1.6051822) q[3];
sx q[3];
rz(-1.4166098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3225591) q[0];
sx q[0];
rz(-2.366134) q[0];
sx q[0];
rz(2.0490647) q[0];
rz(0.93387261) q[1];
sx q[1];
rz(-1.4051508) q[1];
sx q[1];
rz(2.3563103) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0781495) q[0];
sx q[0];
rz(-0.64029303) q[0];
sx q[0];
rz(0.75115738) q[0];
rz(-0.44194989) q[2];
sx q[2];
rz(-1.0681249) q[2];
sx q[2];
rz(-1.3208226) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8251961) q[1];
sx q[1];
rz(-1.2994718) q[1];
sx q[1];
rz(1.8261938) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0075032) q[3];
sx q[3];
rz(-2.1434382) q[3];
sx q[3];
rz(-2.9331895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2698722) q[2];
sx q[2];
rz(-1.6951268) q[2];
sx q[2];
rz(-0.63228697) q[2];
rz(2.2271683) q[3];
sx q[3];
rz(-1.881003) q[3];
sx q[3];
rz(-2.9818592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3226427) q[0];
sx q[0];
rz(-0.85700789) q[0];
sx q[0];
rz(-0.44347611) q[0];
rz(-0.30287287) q[1];
sx q[1];
rz(-0.58285204) q[1];
sx q[1];
rz(-1.833896) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0904059) q[0];
sx q[0];
rz(-1.4009579) q[0];
sx q[0];
rz(0.30533653) q[0];
rz(1.2595749) q[2];
sx q[2];
rz(-1.0838795) q[2];
sx q[2];
rz(-2.4608909) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9759226) q[1];
sx q[1];
rz(-0.56708401) q[1];
sx q[1];
rz(-2.0815064) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1122718) q[3];
sx q[3];
rz(-0.53526894) q[3];
sx q[3];
rz(2.1524371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.12094721) q[2];
sx q[2];
rz(-1.0409313) q[2];
sx q[2];
rz(2.5282703) q[2];
rz(-2.3239418) q[3];
sx q[3];
rz(-0.48912564) q[3];
sx q[3];
rz(0.31942719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9870975) q[0];
sx q[0];
rz(-1.7648062) q[0];
sx q[0];
rz(0.35926551) q[0];
rz(-0.66419762) q[1];
sx q[1];
rz(-1.3281053) q[1];
sx q[1];
rz(-0.42116234) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0764685) q[0];
sx q[0];
rz(-1.9389922) q[0];
sx q[0];
rz(2.9028331) q[0];
rz(-pi) q[1];
rz(0.0053275544) q[2];
sx q[2];
rz(-1.6432646) q[2];
sx q[2];
rz(0.42124149) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9778455) q[1];
sx q[1];
rz(-2.6085269) q[1];
sx q[1];
rz(-1.0372437) q[1];
rz(-0.60539104) q[3];
sx q[3];
rz(-2.0642114) q[3];
sx q[3];
rz(-0.32829729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8757214) q[2];
sx q[2];
rz(-0.096780626) q[2];
sx q[2];
rz(-1.2145915) q[2];
rz(-0.38861361) q[3];
sx q[3];
rz(-0.92258421) q[3];
sx q[3];
rz(2.0428366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.494396) q[0];
sx q[0];
rz(-1.5460486) q[0];
sx q[0];
rz(-1.8478951) q[0];
rz(-2.9248059) q[1];
sx q[1];
rz(-2.4516791) q[1];
sx q[1];
rz(-2.6286415) q[1];
rz(-0.67736652) q[2];
sx q[2];
rz(-1.9147041) q[2];
sx q[2];
rz(-0.83821581) q[2];
rz(-0.46617266) q[3];
sx q[3];
rz(-1.331658) q[3];
sx q[3];
rz(-0.030203947) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
