OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.225086) q[0];
sx q[0];
rz(-0.14296159) q[0];
sx q[0];
rz(11.077865) q[0];
rz(-2.0090964) q[1];
sx q[1];
rz(5.8512591) q[1];
sx q[1];
rz(6.9195256) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37381946) q[0];
sx q[0];
rz(-2.1507906) q[0];
sx q[0];
rz(2.679002) q[0];
rz(-pi) q[1];
x q[1];
rz(0.77677988) q[2];
sx q[2];
rz(-0.52195364) q[2];
sx q[2];
rz(-1.7372434) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.072187034) q[1];
sx q[1];
rz(-1.9652307) q[1];
sx q[1];
rz(-3.0071114) q[1];
x q[2];
rz(2.2872055) q[3];
sx q[3];
rz(-2.1463089) q[3];
sx q[3];
rz(-0.42785409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1119614) q[2];
sx q[2];
rz(-1.1250857) q[2];
sx q[2];
rz(-0.65043989) q[2];
rz(-1.8247617) q[3];
sx q[3];
rz(-1.7951169) q[3];
sx q[3];
rz(0.10183798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0541075) q[0];
sx q[0];
rz(-0.58687812) q[0];
sx q[0];
rz(0.45463872) q[0];
rz(3.1184323) q[1];
sx q[1];
rz(-1.6975479) q[1];
sx q[1];
rz(2.815411) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21837337) q[0];
sx q[0];
rz(-2.120864) q[0];
sx q[0];
rz(2.0327507) q[0];
rz(-pi) q[1];
x q[1];
rz(1.147338) q[2];
sx q[2];
rz(-1.8427765) q[2];
sx q[2];
rz(-0.47951298) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7190951) q[1];
sx q[1];
rz(-2.0764253) q[1];
sx q[1];
rz(-2.0431594) q[1];
rz(-pi) q[2];
rz(-2.2315027) q[3];
sx q[3];
rz(-1.6636563) q[3];
sx q[3];
rz(1.4610491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0252016) q[2];
sx q[2];
rz(-1.2025183) q[2];
sx q[2];
rz(0.95428673) q[2];
rz(1.8918234) q[3];
sx q[3];
rz(-0.3229177) q[3];
sx q[3];
rz(-3.1402816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(1.2713476) q[0];
sx q[0];
rz(-1.988669) q[0];
sx q[0];
rz(1.1767607) q[0];
rz(-0.36508834) q[1];
sx q[1];
rz(-1.8389713) q[1];
sx q[1];
rz(-1.6129859) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9823526) q[0];
sx q[0];
rz(-1.5499347) q[0];
sx q[0];
rz(-3.1283035) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4055785) q[2];
sx q[2];
rz(-2.6927136) q[2];
sx q[2];
rz(0.53476483) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4903725) q[1];
sx q[1];
rz(-1.4191322) q[1];
sx q[1];
rz(-0.030812736) q[1];
rz(-pi) q[2];
rz(3.0729483) q[3];
sx q[3];
rz(-1.513283) q[3];
sx q[3];
rz(0.36348331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8927346) q[2];
sx q[2];
rz(-0.55344075) q[2];
sx q[2];
rz(-1.9888606) q[2];
rz(-1.8966127) q[3];
sx q[3];
rz(-0.98074073) q[3];
sx q[3];
rz(2.8583756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1588441) q[0];
sx q[0];
rz(-11*pi/12) q[0];
sx q[0];
rz(-2.4520279) q[0];
rz(-0.025029643) q[1];
sx q[1];
rz(-1.6233147) q[1];
sx q[1];
rz(1.7207346) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8129028) q[0];
sx q[0];
rz(-2.5089009) q[0];
sx q[0];
rz(2.9660874) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1406607) q[2];
sx q[2];
rz(-1.4073512) q[2];
sx q[2];
rz(2.8317833) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.017548718) q[1];
sx q[1];
rz(-2.6376403) q[1];
sx q[1];
rz(-0.057226463) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2989276) q[3];
sx q[3];
rz(-1.2449578) q[3];
sx q[3];
rz(-3.1414023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.20874061) q[2];
sx q[2];
rz(-0.17540652) q[2];
sx q[2];
rz(-1.9746732) q[2];
rz(-2.6973727) q[3];
sx q[3];
rz(-1.2362213) q[3];
sx q[3];
rz(2.8319632) q[3];
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
rz(-1.2883478) q[0];
sx q[0];
rz(-1.7720368) q[0];
sx q[0];
rz(2.7254768) q[0];
rz(-2.6314645) q[1];
sx q[1];
rz(-2.3704539) q[1];
sx q[1];
rz(2.3663734) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32019553) q[0];
sx q[0];
rz(-2.1761083) q[0];
sx q[0];
rz(1.8138252) q[0];
rz(2.7363766) q[2];
sx q[2];
rz(-1.7381258) q[2];
sx q[2];
rz(-2.2115603) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0803804) q[1];
sx q[1];
rz(-0.16875544) q[1];
sx q[1];
rz(-0.83909281) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.52948496) q[3];
sx q[3];
rz(-2.4765402) q[3];
sx q[3];
rz(-2.3163026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.2148718) q[2];
sx q[2];
rz(-2.1472609) q[2];
sx q[2];
rz(2.1057687) q[2];
rz(0.18946798) q[3];
sx q[3];
rz(-2.6697956) q[3];
sx q[3];
rz(-0.50260472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.254461) q[0];
sx q[0];
rz(-0.04763617) q[0];
sx q[0];
rz(-2.7857842) q[0];
rz(1.4954781) q[1];
sx q[1];
rz(-0.9551841) q[1];
sx q[1];
rz(-1.1411508) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4116652) q[0];
sx q[0];
rz(-2.228274) q[0];
sx q[0];
rz(-0.2667) q[0];
rz(3.1325794) q[2];
sx q[2];
rz(-2.1823332) q[2];
sx q[2];
rz(2.1076941) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.49887744) q[1];
sx q[1];
rz(-1.6756454) q[1];
sx q[1];
rz(0.14123209) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1294501) q[3];
sx q[3];
rz(-2.0046765) q[3];
sx q[3];
rz(-1.3557243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5100539) q[2];
sx q[2];
rz(-2.4271991) q[2];
sx q[2];
rz(-1.4368524) q[2];
rz(-2.0884183) q[3];
sx q[3];
rz(-1.5318233) q[3];
sx q[3];
rz(0.50301445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.961504) q[0];
sx q[0];
rz(-2.8453974) q[0];
sx q[0];
rz(-3.0190813) q[0];
rz(-0.65613121) q[1];
sx q[1];
rz(-1.8729112) q[1];
sx q[1];
rz(1.8001385) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65402664) q[0];
sx q[0];
rz(-1.5254505) q[0];
sx q[0];
rz(-1.6066172) q[0];
rz(-pi) q[1];
rz(-2.7846863) q[2];
sx q[2];
rz(-2.8404337) q[2];
sx q[2];
rz(-0.44277175) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.70133658) q[1];
sx q[1];
rz(-2.1783268) q[1];
sx q[1];
rz(-2.8154147) q[1];
rz(-pi) q[2];
rz(-1.6844325) q[3];
sx q[3];
rz(-1.1294147) q[3];
sx q[3];
rz(1.4809004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5682257) q[2];
sx q[2];
rz(-1.2218916) q[2];
sx q[2];
rz(1.6778256) q[2];
rz(0.73417869) q[3];
sx q[3];
rz(-0.28407431) q[3];
sx q[3];
rz(1.8370139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
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
rz(-1.7952591) q[0];
sx q[0];
rz(-1.7743552) q[0];
sx q[0];
rz(-2.6829868) q[0];
rz(-2.1268225) q[1];
sx q[1];
rz(-1.1943123) q[1];
sx q[1];
rz(1.0626622) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.053772702) q[0];
sx q[0];
rz(-1.4045949) q[0];
sx q[0];
rz(-1.6292162) q[0];
rz(2.3331679) q[2];
sx q[2];
rz(-1.171441) q[2];
sx q[2];
rz(-2.6724919) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6148551) q[1];
sx q[1];
rz(-1.7018082) q[1];
sx q[1];
rz(3.0796771) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0902781) q[3];
sx q[3];
rz(-1.2881345) q[3];
sx q[3];
rz(1.1804898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8899272) q[2];
sx q[2];
rz(-1.7532316) q[2];
sx q[2];
rz(1.3512705) q[2];
rz(-1.9225559) q[3];
sx q[3];
rz(-0.42870298) q[3];
sx q[3];
rz(0.8333227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0625921) q[0];
sx q[0];
rz(-1.7575678) q[0];
sx q[0];
rz(-2.2341527) q[0];
rz(0.22706789) q[1];
sx q[1];
rz(-1.0457958) q[1];
sx q[1];
rz(2.2423832) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.718841) q[0];
sx q[0];
rz(-2.1678574) q[0];
sx q[0];
rz(0.62946749) q[0];
rz(-pi) q[1];
rz(-1.4883409) q[2];
sx q[2];
rz(-1.6489779) q[2];
sx q[2];
rz(-2.9417335) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2503629) q[1];
sx q[1];
rz(-2.0768407) q[1];
sx q[1];
rz(2.9292468) q[1];
rz(3.090631) q[3];
sx q[3];
rz(-0.82333857) q[3];
sx q[3];
rz(-2.2281102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1850618) q[2];
sx q[2];
rz(-1.1659634) q[2];
sx q[2];
rz(0.60066191) q[2];
rz(-1.0968084) q[3];
sx q[3];
rz(-2.5513702) q[3];
sx q[3];
rz(-0.87305951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7507062) q[0];
sx q[0];
rz(-2.0426671) q[0];
sx q[0];
rz(-2.1571958) q[0];
rz(2.6495433) q[1];
sx q[1];
rz(-1.7285408) q[1];
sx q[1];
rz(1.9459928) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0988783) q[0];
sx q[0];
rz(-2.1568568) q[0];
sx q[0];
rz(1.1436966) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2480505) q[2];
sx q[2];
rz(-1.1718501) q[2];
sx q[2];
rz(0.91541399) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.25567935) q[1];
sx q[1];
rz(-1.1885841) q[1];
sx q[1];
rz(-0.087319386) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.672555) q[3];
sx q[3];
rz(-0.82995633) q[3];
sx q[3];
rz(-1.1560464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.055858) q[2];
sx q[2];
rz(-1.5208289) q[2];
sx q[2];
rz(0.98304191) q[2];
rz(-2.8899657) q[3];
sx q[3];
rz(-0.42767891) q[3];
sx q[3];
rz(2.1380641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6428103) q[0];
sx q[0];
rz(-2.1693873) q[0];
sx q[0];
rz(0.65709773) q[0];
rz(-1.2596754) q[1];
sx q[1];
rz(-1.7301662) q[1];
sx q[1];
rz(0.87283254) q[1];
rz(1.8173366) q[2];
sx q[2];
rz(-2.4132072) q[2];
sx q[2];
rz(-2.5843399) q[2];
rz(0.76446492) q[3];
sx q[3];
rz(-0.44776147) q[3];
sx q[3];
rz(0.54816435) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
