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
rz(3.8430619) q[1];
sx q[1];
rz(3.7624533) q[1];
sx q[1];
rz(7.9384595) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1264869) q[0];
sx q[0];
rz(-1.6886212) q[0];
sx q[0];
rz(0.18527041) q[0];
x q[1];
rz(-0.42278554) q[2];
sx q[2];
rz(-1.6362731) q[2];
sx q[2];
rz(-0.72725429) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.15581407) q[1];
sx q[1];
rz(-1.6667546) q[1];
sx q[1];
rz(2.7491261) q[1];
rz(-pi) q[2];
rz(1.5774044) q[3];
sx q[3];
rz(-1.6497383) q[3];
sx q[3];
rz(-3.1066432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0554589) q[2];
sx q[2];
rz(-1.7666585) q[2];
sx q[2];
rz(-1.3936183) q[2];
rz(2.3404775) q[3];
sx q[3];
rz(-1.5812185) q[3];
sx q[3];
rz(1.0629517) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
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
rz(-3.063607) q[0];
rz(0.69411913) q[1];
sx q[1];
rz(-1.1342987) q[1];
sx q[1];
rz(2.7819113) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7168424) q[0];
sx q[0];
rz(-1.0263838) q[0];
sx q[0];
rz(-0.74290468) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9042269) q[2];
sx q[2];
rz(-0.73420213) q[2];
sx q[2];
rz(-1.7422301) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.8453796) q[1];
sx q[1];
rz(-1.4040596) q[1];
sx q[1];
rz(2.1026033) q[1];
rz(-pi) q[2];
rz(-0.039814063) q[3];
sx q[3];
rz(-2.2844444) q[3];
sx q[3];
rz(-3.0301222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2461207) q[2];
sx q[2];
rz(-1.8254231) q[2];
sx q[2];
rz(2.1982511) q[2];
rz(1.881276) q[3];
sx q[3];
rz(-0.80788079) q[3];
sx q[3];
rz(-0.70596203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.869732) q[0];
sx q[0];
rz(-1.4881217) q[0];
sx q[0];
rz(0.46762064) q[0];
rz(-2.6768661) q[1];
sx q[1];
rz(-2.6578891) q[1];
sx q[1];
rz(-0.31276774) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8899925) q[0];
sx q[0];
rz(-1.359299) q[0];
sx q[0];
rz(3.0598346) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8094693) q[2];
sx q[2];
rz(-2.2687074) q[2];
sx q[2];
rz(2.3534564) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.028542472) q[1];
sx q[1];
rz(-0.22284914) q[1];
sx q[1];
rz(-2.261496) q[1];
rz(-pi) q[2];
rz(2.7346482) q[3];
sx q[3];
rz(-0.35396468) q[3];
sx q[3];
rz(2.0719761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.23190817) q[2];
sx q[2];
rz(-1.568305) q[2];
sx q[2];
rz(-2.8364733) q[2];
rz(-1.336608) q[3];
sx q[3];
rz(-0.87953416) q[3];
sx q[3];
rz(2.6763776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6453648) q[0];
sx q[0];
rz(-0.35571686) q[0];
sx q[0];
rz(-3.1135476) q[0];
rz(-1.4656981) q[1];
sx q[1];
rz(-0.60246712) q[1];
sx q[1];
rz(0.67273295) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4417277) q[0];
sx q[0];
rz(-1.2727951) q[0];
sx q[0];
rz(2.5900048) q[0];
rz(-pi) q[1];
rz(-1.3139309) q[2];
sx q[2];
rz(-2.1055524) q[2];
sx q[2];
rz(-1.9321439) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2697301) q[1];
sx q[1];
rz(-1.6849663) q[1];
sx q[1];
rz(-1.3378548) q[1];
rz(-pi) q[2];
rz(-1.7476051) q[3];
sx q[3];
rz(-2.5150931) q[3];
sx q[3];
rz(0.83229438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4611886) q[2];
sx q[2];
rz(-0.79128069) q[2];
sx q[2];
rz(0.76081863) q[2];
rz(2.2740254) q[3];
sx q[3];
rz(-1.3191728) q[3];
sx q[3];
rz(-0.76604617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5714394) q[0];
sx q[0];
rz(-1.1076936) q[0];
sx q[0];
rz(-1.0636348) q[0];
rz(-1.6617552) q[1];
sx q[1];
rz(-0.96264833) q[1];
sx q[1];
rz(0.038287727) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1153039) q[0];
sx q[0];
rz(-2.8235108) q[0];
sx q[0];
rz(-1.9684857) q[0];
rz(-2.5262881) q[2];
sx q[2];
rz(-1.9822789) q[2];
sx q[2];
rz(-0.40508168) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8582279) q[1];
sx q[1];
rz(-1.7444532) q[1];
sx q[1];
rz(-0.11939343) q[1];
rz(1.0682265) q[3];
sx q[3];
rz(-2.6196819) q[3];
sx q[3];
rz(0.13773242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9866207) q[2];
sx q[2];
rz(-1.9382696) q[2];
sx q[2];
rz(-2.0495474) q[2];
rz(1.5345705) q[3];
sx q[3];
rz(-2.1708596) q[3];
sx q[3];
rz(-2.3585414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9606278) q[0];
sx q[0];
rz(-1.5009078) q[0];
sx q[0];
rz(-1.8339748) q[0];
rz(3.0788105) q[1];
sx q[1];
rz(-1.9654013) q[1];
sx q[1];
rz(-0.19097701) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67577772) q[0];
sx q[0];
rz(-2.8880279) q[0];
sx q[0];
rz(-2.158014) q[0];
rz(-1.9917166) q[2];
sx q[2];
rz(-2.4119666) q[2];
sx q[2];
rz(-0.52044496) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8456965) q[1];
sx q[1];
rz(-1.0277883) q[1];
sx q[1];
rz(-0.023881749) q[1];
rz(-1.005638) q[3];
sx q[3];
rz(-1.2747545) q[3];
sx q[3];
rz(-2.5067096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.442231) q[2];
sx q[2];
rz(-1.3012393) q[2];
sx q[2];
rz(2.7065275) q[2];
rz(2.2502031) q[3];
sx q[3];
rz(-1.1957217) q[3];
sx q[3];
rz(1.5589327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1680465) q[0];
sx q[0];
rz(-1.1991812) q[0];
sx q[0];
rz(2.6269873) q[0];
x q[1];
rz(1.2375564) q[2];
sx q[2];
rz(-2.9102647) q[2];
sx q[2];
rz(-1.7817792) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8458) q[1];
sx q[1];
rz(-0.35357057) q[1];
sx q[1];
rz(0.54770623) q[1];
x q[2];
rz(0.78564268) q[3];
sx q[3];
rz(-0.26587405) q[3];
sx q[3];
rz(3.0674988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7421425) q[2];
sx q[2];
rz(-0.43562499) q[2];
sx q[2];
rz(2.8263261) q[2];
rz(1.0366084) q[3];
sx q[3];
rz(-1.8866115) q[3];
sx q[3];
rz(-0.96926564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8740365) q[0];
sx q[0];
rz(-0.75796217) q[0];
sx q[0];
rz(-1.1918921) q[0];
rz(-0.9115971) q[1];
sx q[1];
rz(-0.67271581) q[1];
sx q[1];
rz(-1.0620767) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22411352) q[0];
sx q[0];
rz(-0.55372059) q[0];
sx q[0];
rz(-1.327716) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9757134) q[2];
sx q[2];
rz(-1.9654044) q[2];
sx q[2];
rz(-1.504606) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2506822) q[1];
sx q[1];
rz(-0.85103304) q[1];
sx q[1];
rz(0.85810424) q[1];
rz(-1.5261371) q[3];
sx q[3];
rz(-2.4756845) q[3];
sx q[3];
rz(2.8093616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.186211) q[2];
sx q[2];
rz(-0.28327981) q[2];
sx q[2];
rz(-0.051076802) q[2];
rz(-0.71930277) q[3];
sx q[3];
rz(-1.5161113) q[3];
sx q[3];
rz(-2.0843845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99701571) q[0];
sx q[0];
rz(-0.41467312) q[0];
sx q[0];
rz(-0.64494079) q[0];
rz(-2.0058517) q[1];
sx q[1];
rz(-2.203511) q[1];
sx q[1];
rz(-2.2924246) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1695771) q[0];
sx q[0];
rz(-1.2017842) q[0];
sx q[0];
rz(-0.24032648) q[0];
x q[1];
rz(-1.7230942) q[2];
sx q[2];
rz(-1.2010152) q[2];
sx q[2];
rz(0.6347444) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2850212) q[1];
sx q[1];
rz(-1.3339431) q[1];
sx q[1];
rz(0.23394886) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9407223) q[3];
sx q[3];
rz(-1.7327274) q[3];
sx q[3];
rz(-2.014132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.34716216) q[2];
sx q[2];
rz(-0.55134761) q[2];
sx q[2];
rz(2.9629422) q[2];
rz(-2.3000439) q[3];
sx q[3];
rz(-2.0984564) q[3];
sx q[3];
rz(-0.49079045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.08298824) q[0];
sx q[0];
rz(-1.937979) q[0];
sx q[0];
rz(0.9151181) q[0];
rz(1.8944342) q[1];
sx q[1];
rz(-1.1727389) q[1];
sx q[1];
rz(-2.9291005) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2895848) q[0];
sx q[0];
rz(-0.6913018) q[0];
sx q[0];
rz(1.9612802) q[0];
x q[1];
rz(1.882952) q[2];
sx q[2];
rz(-1.6098795) q[2];
sx q[2];
rz(-2.8651819) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4957628) q[1];
sx q[1];
rz(-0.96853515) q[1];
sx q[1];
rz(-0.61969238) q[1];
rz(-pi) q[2];
rz(-0.28189567) q[3];
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
rz(2.1797669) q[2];
sx q[2];
rz(-3.0111854) q[2];
sx q[2];
rz(-1.4546222) q[2];
rz(-1.1949332) q[3];
sx q[3];
rz(-1.8353728) q[3];
sx q[3];
rz(0.52836829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-1.4029978) q[0];
sx q[0];
rz(-1.7398555) q[0];
sx q[0];
rz(1.6237988) q[0];
rz(-0.075642792) q[1];
sx q[1];
rz(-1.6041926) q[1];
sx q[1];
rz(1.4354482) q[1];
rz(0.018013714) q[2];
sx q[2];
rz(-1.6406888) q[2];
sx q[2];
rz(-2.3802118) q[2];
rz(-2.0316488) q[3];
sx q[3];
rz(-2.0416595) q[3];
sx q[3];
rz(2.5817081) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
