OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5869863) q[0];
sx q[0];
rz(-2.245683) q[0];
sx q[0];
rz(-1.7370261) q[0];
rz(-2.4401234) q[1];
sx q[1];
rz(-2.520732) q[1];
sx q[1];
rz(-1.4863185) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0044999997) q[0];
sx q[0];
rz(-0.21919964) q[0];
sx q[0];
rz(-0.57114925) q[0];
x q[1];
rz(-1.6425743) q[2];
sx q[2];
rz(-1.1489747) q[2];
sx q[2];
rz(0.81411241) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6424375) q[1];
sx q[1];
rz(-2.7381574) q[1];
sx q[1];
rz(2.8950476) q[1];
rz(-pi) q[2];
rz(0.078943723) q[3];
sx q[3];
rz(-1.5642089) q[3];
sx q[3];
rz(-1.5353257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0861337) q[2];
sx q[2];
rz(-1.3749342) q[2];
sx q[2];
rz(-1.7479744) q[2];
rz(-0.80111516) q[3];
sx q[3];
rz(-1.5812185) q[3];
sx q[3];
rz(1.0629517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(-3.063607) q[0];
rz(-0.69411913) q[1];
sx q[1];
rz(-1.1342987) q[1];
sx q[1];
rz(-2.7819113) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29794824) q[0];
sx q[0];
rz(-2.1878562) q[0];
sx q[0];
rz(-0.88275568) q[0];
rz(-pi) q[1];
x q[1];
rz(1.361679) q[2];
sx q[2];
rz(-0.86161999) q[2];
sx q[2];
rz(-1.7143955) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.45047346) q[1];
sx q[1];
rz(-0.55492655) q[1];
sx q[1];
rz(1.8912485) q[1];
rz(3.1017786) q[3];
sx q[3];
rz(-2.2844444) q[3];
sx q[3];
rz(0.1114705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2461207) q[2];
sx q[2];
rz(-1.8254231) q[2];
sx q[2];
rz(0.94334156) q[2];
rz(-1.2603166) q[3];
sx q[3];
rz(-2.3337119) q[3];
sx q[3];
rz(0.70596203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27186069) q[0];
sx q[0];
rz(-1.6534709) q[0];
sx q[0];
rz(-2.673972) q[0];
rz(2.6768661) q[1];
sx q[1];
rz(-2.6578891) q[1];
sx q[1];
rz(-2.8288249) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62374672) q[0];
sx q[0];
rz(-2.9150634) q[0];
sx q[0];
rz(1.9342599) q[0];
rz(-pi) q[1];
x q[1];
rz(0.84509387) q[2];
sx q[2];
rz(-1.8232864) q[2];
sx q[2];
rz(1.0007728) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2777729) q[1];
sx q[1];
rz(-1.7120655) q[1];
sx q[1];
rz(-1.743725) q[1];
rz(-pi) q[2];
rz(1.7160277) q[3];
sx q[3];
rz(-1.2468474) q[3];
sx q[3];
rz(2.5026929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.23190817) q[2];
sx q[2];
rz(-1.5732876) q[2];
sx q[2];
rz(2.8364733) q[2];
rz(1.8049847) q[3];
sx q[3];
rz(-2.2620585) q[3];
sx q[3];
rz(0.46521503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6453648) q[0];
sx q[0];
rz(-2.7858758) q[0];
sx q[0];
rz(-0.028045068) q[0];
rz(-1.4656981) q[1];
sx q[1];
rz(-2.5391255) q[1];
sx q[1];
rz(-0.67273295) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4493895) q[0];
sx q[0];
rz(-1.0461079) q[0];
sx q[0];
rz(-1.2246818) q[0];
rz(-pi) q[1];
rz(-2.5920934) q[2];
sx q[2];
rz(-1.7911583) q[2];
sx q[2];
rz(-2.9133177) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2697301) q[1];
sx q[1];
rz(-1.4566263) q[1];
sx q[1];
rz(-1.8037379) q[1];
x q[2];
rz(3.014971) q[3];
sx q[3];
rz(-2.1860578) q[3];
sx q[3];
rz(-0.61520731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.68040401) q[2];
sx q[2];
rz(-0.79128069) q[2];
sx q[2];
rz(2.380774) q[2];
rz(0.86756724) q[3];
sx q[3];
rz(-1.3191728) q[3];
sx q[3];
rz(0.76604617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5701533) q[0];
sx q[0];
rz(-1.1076936) q[0];
sx q[0];
rz(-1.0636348) q[0];
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
rz(0.83513658) q[0];
sx q[0];
rz(-1.4493754) q[0];
sx q[0];
rz(1.2760713) q[0];
x q[1];
rz(-0.61530453) q[2];
sx q[2];
rz(-1.1593137) q[2];
sx q[2];
rz(2.736511) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8582279) q[1];
sx q[1];
rz(-1.3971395) q[1];
sx q[1];
rz(0.11939343) q[1];
rz(-2.0376301) q[3];
sx q[3];
rz(-1.8133014) q[3];
sx q[3];
rz(1.2638826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9866207) q[2];
sx q[2];
rz(-1.203323) q[2];
sx q[2];
rz(-1.0920452) q[2];
rz(-1.6070222) q[3];
sx q[3];
rz(-0.97073308) q[3];
sx q[3];
rz(-0.78305125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9606278) q[0];
sx q[0];
rz(-1.6406849) q[0];
sx q[0];
rz(-1.8339748) q[0];
rz(0.062782137) q[1];
sx q[1];
rz(-1.9654013) q[1];
sx q[1];
rz(0.19097701) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.068114) q[0];
sx q[0];
rz(-1.7811791) q[0];
sx q[0];
rz(2.99899) q[0];
x q[1];
rz(-2.7912748) q[2];
sx q[2];
rz(-2.2248474) q[2];
sx q[2];
rz(2.0803114) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.2497017) q[1];
sx q[1];
rz(-0.54348031) q[1];
sx q[1];
rz(1.6103423) q[1];
rz(1.005638) q[3];
sx q[3];
rz(-1.2747545) q[3];
sx q[3];
rz(2.5067096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.442231) q[2];
sx q[2];
rz(-1.8403534) q[2];
sx q[2];
rz(2.7065275) q[2];
rz(0.89138952) q[3];
sx q[3];
rz(-1.9458709) q[3];
sx q[3];
rz(1.5589327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-2.1973535) q[0];
sx q[0];
rz(-2.9159912) q[0];
sx q[0];
rz(-2.7213668) q[0];
rz(0.48121437) q[1];
sx q[1];
rz(-0.88164202) q[1];
sx q[1];
rz(-2.1113254) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9734605) q[0];
sx q[0];
rz(-2.5167743) q[0];
sx q[0];
rz(-2.4718651) q[0];
x q[1];
rz(-1.7898125) q[2];
sx q[2];
rz(-1.6458626) q[2];
sx q[2];
rz(2.605627) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.75525857) q[1];
sx q[1];
rz(-1.3895021) q[1];
sx q[1];
rz(-2.8363486) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9514063) q[3];
sx q[3];
rz(-1.38387) q[3];
sx q[3];
rz(-0.72894064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7421425) q[2];
sx q[2];
rz(-0.43562499) q[2];
sx q[2];
rz(-2.8263261) q[2];
rz(2.1049843) q[3];
sx q[3];
rz(-1.8866115) q[3];
sx q[3];
rz(0.96926564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2675562) q[0];
sx q[0];
rz(-0.75796217) q[0];
sx q[0];
rz(1.9497005) q[0];
rz(0.9115971) q[1];
sx q[1];
rz(-2.4688768) q[1];
sx q[1];
rz(2.079516) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22411352) q[0];
sx q[0];
rz(-2.5878721) q[0];
sx q[0];
rz(1.327716) q[0];
rz(-pi) q[1];
rz(0.75762962) q[2];
sx q[2];
rz(-0.55765753) q[2];
sx q[2];
rz(2.4766395) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.83795611) q[1];
sx q[1];
rz(-1.0567697) q[1];
sx q[1];
rz(2.282826) q[1];
rz(0.90537269) q[3];
sx q[3];
rz(-1.5983799) q[3];
sx q[3];
rz(-1.2736922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9553817) q[2];
sx q[2];
rz(-0.28327981) q[2];
sx q[2];
rz(0.051076802) q[2];
rz(-2.4222899) q[3];
sx q[3];
rz(-1.5161113) q[3];
sx q[3];
rz(-1.0572082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1445769) q[0];
sx q[0];
rz(-0.41467312) q[0];
sx q[0];
rz(2.4966519) q[0];
rz(-1.1357409) q[1];
sx q[1];
rz(-0.93808162) q[1];
sx q[1];
rz(0.84916806) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5687953) q[0];
sx q[0];
rz(-2.704247) q[0];
sx q[0];
rz(1.019078) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7230942) q[2];
sx q[2];
rz(-1.2010152) q[2];
sx q[2];
rz(0.6347444) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.65836425) q[1];
sx q[1];
rz(-1.3434957) q[1];
sx q[1];
rz(1.3275654) q[1];
rz(-pi) q[2];
rz(0.68607761) q[3];
sx q[3];
rz(-0.25732532) q[3];
sx q[3];
rz(1.1130594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.34716216) q[2];
sx q[2];
rz(-0.55134761) q[2];
sx q[2];
rz(2.9629422) q[2];
rz(2.3000439) q[3];
sx q[3];
rz(-2.0984564) q[3];
sx q[3];
rz(0.49079045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.08298824) q[0];
sx q[0];
rz(-1.937979) q[0];
sx q[0];
rz(-2.2264746) q[0];
rz(-1.8944342) q[1];
sx q[1];
rz(-1.1727389) q[1];
sx q[1];
rz(2.9291005) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.025876306) q[0];
sx q[0];
rz(-1.8159144) q[0];
sx q[0];
rz(2.2239767) q[0];
x q[1];
rz(1.2586406) q[2];
sx q[2];
rz(-1.6098795) q[2];
sx q[2];
rz(-0.27641073) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5455268) q[1];
sx q[1];
rz(-0.83546987) q[1];
sx q[1];
rz(-0.86931418) q[1];
rz(1.8260164) q[3];
sx q[3];
rz(-2.2877684) q[3];
sx q[3];
rz(-2.9361801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1797669) q[2];
sx q[2];
rz(-0.13040725) q[2];
sx q[2];
rz(-1.6869705) q[2];
rz(1.1949332) q[3];
sx q[3];
rz(-1.3062198) q[3];
sx q[3];
rz(-2.6132244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7385948) q[0];
sx q[0];
rz(-1.7398555) q[0];
sx q[0];
rz(1.6237988) q[0];
rz(3.0659499) q[1];
sx q[1];
rz(-1.6041926) q[1];
sx q[1];
rz(1.4354482) q[1];
rz(-1.8226345) q[2];
sx q[2];
rz(-3.0694198) q[2];
sx q[2];
rz(1.0138489) q[2];
rz(-2.4235519) q[3];
sx q[3];
rz(-2.4951039) q[3];
sx q[3];
rz(-1.390355) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
