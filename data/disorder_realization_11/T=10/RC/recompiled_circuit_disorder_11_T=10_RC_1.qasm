OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.55460632) q[0];
sx q[0];
rz(2.245683) q[0];
sx q[0];
rz(10.829344) q[0];
rz(0.70146927) q[1];
sx q[1];
rz(-0.62086064) q[1];
sx q[1];
rz(-1.6552742) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1264869) q[0];
sx q[0];
rz(-1.4529714) q[0];
sx q[0];
rz(-2.9563222) q[0];
rz(-pi) q[1];
rz(1.4990184) q[2];
sx q[2];
rz(-1.992618) q[2];
sx q[2];
rz(2.3274802) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7662498) q[1];
sx q[1];
rz(-1.9613593) q[1];
sx q[1];
rz(1.4669963) q[1];
rz(0.08333929) q[3];
sx q[3];
rz(-3.0623751) q[3];
sx q[3];
rz(3.0230429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0554589) q[2];
sx q[2];
rz(-1.3749342) q[2];
sx q[2];
rz(-1.7479744) q[2];
rz(-2.3404775) q[3];
sx q[3];
rz(-1.5812185) q[3];
sx q[3];
rz(-1.0629517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3515117) q[0];
sx q[0];
rz(-1.4485437) q[0];
sx q[0];
rz(0.077985667) q[0];
rz(2.4474735) q[1];
sx q[1];
rz(-1.1342987) q[1];
sx q[1];
rz(-2.7819113) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7168424) q[0];
sx q[0];
rz(-1.0263838) q[0];
sx q[0];
rz(-2.398688) q[0];
rz(-pi) q[1];
rz(0.23736575) q[2];
sx q[2];
rz(-0.73420213) q[2];
sx q[2];
rz(-1.3993625) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.8453796) q[1];
sx q[1];
rz(-1.737533) q[1];
sx q[1];
rz(-1.0389894) q[1];
x q[2];
rz(-1.5248604) q[3];
sx q[3];
rz(-0.71456281) q[3];
sx q[3];
rz(-0.050690953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.89547196) q[2];
sx q[2];
rz(-1.3161696) q[2];
sx q[2];
rz(-0.94334156) q[2];
rz(-1.881276) q[3];
sx q[3];
rz(-0.80788079) q[3];
sx q[3];
rz(-2.4356306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.869732) q[0];
sx q[0];
rz(-1.6534709) q[0];
sx q[0];
rz(-2.673972) q[0];
rz(-0.46472654) q[1];
sx q[1];
rz(-2.6578891) q[1];
sx q[1];
rz(0.31276774) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8395961) q[0];
sx q[0];
rz(-1.490864) q[0];
sx q[0];
rz(1.3586112) q[0];
rz(-1.9415641) q[2];
sx q[2];
rz(-0.76075483) q[2];
sx q[2];
rz(0.2955546) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.86381972) q[1];
sx q[1];
rz(-1.4295271) q[1];
sx q[1];
rz(-1.3978676) q[1];
rz(-pi) q[2];
rz(0.4069445) q[3];
sx q[3];
rz(-0.35396468) q[3];
sx q[3];
rz(-2.0719761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.23190817) q[2];
sx q[2];
rz(-1.568305) q[2];
sx q[2];
rz(-2.8364733) q[2];
rz(-1.8049847) q[3];
sx q[3];
rz(-2.2620585) q[3];
sx q[3];
rz(2.6763776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49622789) q[0];
sx q[0];
rz(-0.35571686) q[0];
sx q[0];
rz(3.1135476) q[0];
rz(1.6758945) q[1];
sx q[1];
rz(-0.60246712) q[1];
sx q[1];
rz(-2.4688597) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4417277) q[0];
sx q[0];
rz(-1.2727951) q[0];
sx q[0];
rz(0.55158789) q[0];
rz(-pi) q[1];
rz(-1.8276617) q[2];
sx q[2];
rz(-2.1055524) q[2];
sx q[2];
rz(1.9321439) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9949918) q[1];
sx q[1];
rz(-0.25895893) q[1];
sx q[1];
rz(-1.109757) q[1];
rz(-pi) q[2];
rz(2.1898515) q[3];
sx q[3];
rz(-1.4674867) q[3];
sx q[3];
rz(2.2593474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4611886) q[2];
sx q[2];
rz(-0.79128069) q[2];
sx q[2];
rz(2.380774) q[2];
rz(-0.86756724) q[3];
sx q[3];
rz(-1.3191728) q[3];
sx q[3];
rz(-0.76604617) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5701533) q[0];
sx q[0];
rz(-2.0338991) q[0];
sx q[0];
rz(-2.0779579) q[0];
rz(1.6617552) q[1];
sx q[1];
rz(-0.96264833) q[1];
sx q[1];
rz(-0.038287727) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3064561) q[0];
sx q[0];
rz(-1.6922173) q[0];
sx q[0];
rz(-1.8655213) q[0];
rz(-2.5262881) q[2];
sx q[2];
rz(-1.1593137) q[2];
sx q[2];
rz(0.40508168) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8748862) q[1];
sx q[1];
rz(-1.453207) q[1];
sx q[1];
rz(1.3959195) q[1];
x q[2];
rz(-1.1039626) q[3];
sx q[3];
rz(-1.8133014) q[3];
sx q[3];
rz(1.8777101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9866207) q[2];
sx q[2];
rz(-1.203323) q[2];
sx q[2];
rz(1.0920452) q[2];
rz(1.6070222) q[3];
sx q[3];
rz(-0.97073308) q[3];
sx q[3];
rz(0.78305125) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9606278) q[0];
sx q[0];
rz(-1.6406849) q[0];
sx q[0];
rz(1.3076179) q[0];
rz(-3.0788105) q[1];
sx q[1];
rz(-1.1761913) q[1];
sx q[1];
rz(2.9506156) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.068114) q[0];
sx q[0];
rz(-1.7811791) q[0];
sx q[0];
rz(0.14260261) q[0];
rz(-pi) q[1];
rz(-0.88626185) q[2];
sx q[2];
rz(-1.2949416) q[2];
sx q[2];
rz(0.72826284) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.29589614) q[1];
sx q[1];
rz(-1.0277883) q[1];
sx q[1];
rz(3.1177109) q[1];
rz(-pi) q[2];
x q[2];
rz(1.005638) q[3];
sx q[3];
rz(-1.2747545) q[3];
sx q[3];
rz(-0.63488301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6993616) q[2];
sx q[2];
rz(-1.8403534) q[2];
sx q[2];
rz(-2.7065275) q[2];
rz(-2.2502031) q[3];
sx q[3];
rz(-1.9458709) q[3];
sx q[3];
rz(1.5589327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1973535) q[0];
sx q[0];
rz(-2.9159912) q[0];
sx q[0];
rz(0.42022589) q[0];
rz(2.6603783) q[1];
sx q[1];
rz(-2.2599506) q[1];
sx q[1];
rz(-2.1113254) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39474836) q[0];
sx q[0];
rz(-1.0943824) q[0];
sx q[0];
rz(-1.1498515) q[0];
rz(-pi) q[1];
rz(-1.2375564) q[2];
sx q[2];
rz(-0.2313279) q[2];
sx q[2];
rz(1.3598134) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.87228862) q[1];
sx q[1];
rz(-1.270712) q[1];
sx q[1];
rz(-1.3809204) q[1];
rz(2.35595) q[3];
sx q[3];
rz(-2.8757186) q[3];
sx q[3];
rz(3.0674988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.39945012) q[2];
sx q[2];
rz(-2.7059677) q[2];
sx q[2];
rz(2.8263261) q[2];
rz(1.0366084) q[3];
sx q[3];
rz(-1.8866115) q[3];
sx q[3];
rz(2.172327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2675562) q[0];
sx q[0];
rz(-2.3836305) q[0];
sx q[0];
rz(-1.1918921) q[0];
rz(-2.2299956) q[1];
sx q[1];
rz(-2.4688768) q[1];
sx q[1];
rz(-1.0620767) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05957219) q[0];
sx q[0];
rz(-1.0351512) q[0];
sx q[0];
rz(2.9938712) q[0];
x q[1];
rz(1.1658792) q[2];
sx q[2];
rz(-1.9654044) q[2];
sx q[2];
rz(1.504606) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.33151367) q[1];
sx q[1];
rz(-2.1760097) q[1];
sx q[1];
rz(0.64085754) q[1];
x q[2];
rz(1.6154556) q[3];
sx q[3];
rz(-0.66590819) q[3];
sx q[3];
rz(0.33223104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.186211) q[2];
sx q[2];
rz(-2.8583128) q[2];
sx q[2];
rz(-0.051076802) q[2];
rz(-2.4222899) q[3];
sx q[3];
rz(-1.6254814) q[3];
sx q[3];
rz(-2.0843845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99701571) q[0];
sx q[0];
rz(-0.41467312) q[0];
sx q[0];
rz(2.4966519) q[0];
rz(1.1357409) q[1];
sx q[1];
rz(-0.93808162) q[1];
sx q[1];
rz(-0.84916806) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6309755) q[0];
sx q[0];
rz(-1.7946587) q[0];
sx q[0];
rz(-1.1918681) q[0];
x q[1];
rz(1.7230942) q[2];
sx q[2];
rz(-1.2010152) q[2];
sx q[2];
rz(2.5068482) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4832284) q[1];
sx q[1];
rz(-1.798097) q[1];
sx q[1];
rz(1.8140273) q[1];
rz(2.455515) q[3];
sx q[3];
rz(-2.8842673) q[3];
sx q[3];
rz(-2.0285332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.34716216) q[2];
sx q[2];
rz(-0.55134761) q[2];
sx q[2];
rz(0.17865044) q[2];
rz(2.3000439) q[3];
sx q[3];
rz(-1.0431362) q[3];
sx q[3];
rz(2.6508022) q[3];
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
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0586044) q[0];
sx q[0];
rz(-1.937979) q[0];
sx q[0];
rz(0.9151181) q[0];
rz(-1.2471584) q[1];
sx q[1];
rz(-1.1727389) q[1];
sx q[1];
rz(0.21249214) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1157163) q[0];
sx q[0];
rz(-1.3256782) q[0];
sx q[0];
rz(-0.91761597) q[0];
rz(0.041065644) q[2];
sx q[2];
rz(-1.2588871) q[2];
sx q[2];
rz(1.8598156) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5455268) q[1];
sx q[1];
rz(-0.83546987) q[1];
sx q[1];
rz(0.86931418) q[1];
rz(-pi) q[2];
x q[2];
rz(2.859697) q[3];
sx q[3];
rz(-0.75337871) q[3];
sx q[3];
rz(2.5582112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1797669) q[2];
sx q[2];
rz(-3.0111854) q[2];
sx q[2];
rz(1.4546222) q[2];
rz(-1.1949332) q[3];
sx q[3];
rz(-1.3062198) q[3];
sx q[3];
rz(2.6132244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7385948) q[0];
sx q[0];
rz(-1.7398555) q[0];
sx q[0];
rz(1.6237988) q[0];
rz(-0.075642792) q[1];
sx q[1];
rz(-1.6041926) q[1];
sx q[1];
rz(1.4354482) q[1];
rz(-1.5008925) q[2];
sx q[2];
rz(-1.5528266) q[2];
sx q[2];
rz(2.3334353) q[2];
rz(0.71804071) q[3];
sx q[3];
rz(-2.4951039) q[3];
sx q[3];
rz(-1.390355) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
