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
rz(-2.2505724) q[0];
sx q[0];
rz(-1.855259) q[0];
sx q[0];
rz(2.2556055) q[0];
rz(-2.6180144) q[1];
sx q[1];
rz(-0.55966592) q[1];
sx q[1];
rz(1.3203415) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.490075) q[0];
sx q[0];
rz(-1.0371672) q[0];
sx q[0];
rz(-2.9367111) q[0];
rz(0.77729887) q[2];
sx q[2];
rz(-1.8172998) q[2];
sx q[2];
rz(0.2709301) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.34397555) q[1];
sx q[1];
rz(-1.6317751) q[1];
sx q[1];
rz(1.0073376) q[1];
rz(-0.76372778) q[3];
sx q[3];
rz(-0.9850174) q[3];
sx q[3];
rz(1.2645185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3628799) q[2];
sx q[2];
rz(-1.820463) q[2];
sx q[2];
rz(-2.2134181) q[2];
rz(1.5556395) q[3];
sx q[3];
rz(-1.0471683) q[3];
sx q[3];
rz(-1.9716523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.25598079) q[0];
sx q[0];
rz(-2.8035127) q[0];
sx q[0];
rz(1.0611435) q[0];
rz(1.9288918) q[1];
sx q[1];
rz(-2.3586528) q[1];
sx q[1];
rz(-1.3139668) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6550094) q[0];
sx q[0];
rz(-0.81196852) q[0];
sx q[0];
rz(-0.49488189) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.92483123) q[2];
sx q[2];
rz(-1.926117) q[2];
sx q[2];
rz(-2.3424847) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.81345859) q[1];
sx q[1];
rz(-0.5603921) q[1];
sx q[1];
rz(-1.8916393) q[1];
x q[2];
rz(-2.6077929) q[3];
sx q[3];
rz(-2.2225923) q[3];
sx q[3];
rz(2.8894412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5706011) q[2];
sx q[2];
rz(-0.97390318) q[2];
sx q[2];
rz(0.904733) q[2];
rz(-1.5500801) q[3];
sx q[3];
rz(-1.2856057) q[3];
sx q[3];
rz(1.2929644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0631113) q[0];
sx q[0];
rz(-0.08992973) q[0];
sx q[0];
rz(-2.2233546) q[0];
rz(-0.69084424) q[1];
sx q[1];
rz(-1.3965239) q[1];
sx q[1];
rz(-2.3450559) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9555003) q[0];
sx q[0];
rz(-1.4791795) q[0];
sx q[0];
rz(1.3618401) q[0];
rz(1.630322) q[2];
sx q[2];
rz(-1.5280485) q[2];
sx q[2];
rz(1.8190365) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.84948375) q[1];
sx q[1];
rz(-0.70569084) q[1];
sx q[1];
rz(1.73805) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1311319) q[3];
sx q[3];
rz(-2.4395203) q[3];
sx q[3];
rz(-1.0855826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.010926509) q[2];
sx q[2];
rz(-0.9504168) q[2];
sx q[2];
rz(1.8939023) q[2];
rz(0.55825663) q[3];
sx q[3];
rz(-1.6853251) q[3];
sx q[3];
rz(-0.021350967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70304865) q[0];
sx q[0];
rz(-0.049228638) q[0];
sx q[0];
rz(0.7274279) q[0];
rz(-2.9797331) q[1];
sx q[1];
rz(-1.9700123) q[1];
sx q[1];
rz(1.1579827) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3129565) q[0];
sx q[0];
rz(-1.4738171) q[0];
sx q[0];
rz(3.0752379) q[0];
x q[1];
rz(0.00096614758) q[2];
sx q[2];
rz(-1.7100934) q[2];
sx q[2];
rz(1.3749706) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.47745142) q[1];
sx q[1];
rz(-0.22813561) q[1];
sx q[1];
rz(-0.99598187) q[1];
rz(-pi) q[2];
rz(-1.4876266) q[3];
sx q[3];
rz(-2.0953878) q[3];
sx q[3];
rz(0.34357854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.27266476) q[2];
sx q[2];
rz(-1.2020035) q[2];
sx q[2];
rz(1.7499917) q[2];
rz(2.9113655) q[3];
sx q[3];
rz(-1.1743436) q[3];
sx q[3];
rz(0.19190425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4153083) q[0];
sx q[0];
rz(-0.10157651) q[0];
sx q[0];
rz(-2.0368077) q[0];
rz(-0.11100189) q[1];
sx q[1];
rz(-2.0260725) q[1];
sx q[1];
rz(-1.312779) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7953781) q[0];
sx q[0];
rz(-1.7722478) q[0];
sx q[0];
rz(-0.18152118) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4823303) q[2];
sx q[2];
rz(-2.1328204) q[2];
sx q[2];
rz(2.1254181) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5178419) q[1];
sx q[1];
rz(-1.1718796) q[1];
sx q[1];
rz(-2.9241882) q[1];
x q[2];
rz(-0.41683414) q[3];
sx q[3];
rz(-2.2456944) q[3];
sx q[3];
rz(1.0639497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.63867265) q[2];
sx q[2];
rz(-2.0827677) q[2];
sx q[2];
rz(-2.1144833) q[2];
rz(0.98661679) q[3];
sx q[3];
rz(-0.28237453) q[3];
sx q[3];
rz(1.7178887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33243772) q[0];
sx q[0];
rz(-1.2323392) q[0];
sx q[0];
rz(-0.80068457) q[0];
rz(2.370131) q[1];
sx q[1];
rz(-2.0636676) q[1];
sx q[1];
rz(-1.7652184) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3478394) q[0];
sx q[0];
rz(-2.3476566) q[0];
sx q[0];
rz(-0.57200045) q[0];
x q[1];
rz(-0.034997392) q[2];
sx q[2];
rz(-1.1049602) q[2];
sx q[2];
rz(-0.83110561) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.47471646) q[1];
sx q[1];
rz(-2.2762269) q[1];
sx q[1];
rz(1.5882467) q[1];
x q[2];
rz(0.8326859) q[3];
sx q[3];
rz(-0.89431864) q[3];
sx q[3];
rz(0.075084075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.8419522) q[2];
sx q[2];
rz(-1.9349293) q[2];
sx q[2];
rz(-3.1237777) q[2];
rz(-0.31339112) q[3];
sx q[3];
rz(-1.0217977) q[3];
sx q[3];
rz(2.4221086) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3571091) q[0];
sx q[0];
rz(-1.5672368) q[0];
sx q[0];
rz(0.16726476) q[0];
rz(-0.74370614) q[1];
sx q[1];
rz(-0.88631648) q[1];
sx q[1];
rz(3.0053265) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24700227) q[0];
sx q[0];
rz(-2.407529) q[0];
sx q[0];
rz(1.9904925) q[0];
rz(-pi) q[1];
rz(-2.9685814) q[2];
sx q[2];
rz(-1.8699416) q[2];
sx q[2];
rz(1.1227705) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.6667585) q[1];
sx q[1];
rz(-1.5038361) q[1];
sx q[1];
rz(2.7826742) q[1];
x q[2];
rz(-0.15382556) q[3];
sx q[3];
rz(-0.3496146) q[3];
sx q[3];
rz(-1.459801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3027489) q[2];
sx q[2];
rz(-0.15041298) q[2];
sx q[2];
rz(2.0602843) q[2];
rz(0.14351621) q[3];
sx q[3];
rz(-1.2986526) q[3];
sx q[3];
rz(1.6941841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3926587) q[0];
sx q[0];
rz(-1.8086139) q[0];
sx q[0];
rz(-0.64312154) q[0];
rz(0.92203036) q[1];
sx q[1];
rz(-2.8755201) q[1];
sx q[1];
rz(1.5111074) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47925835) q[0];
sx q[0];
rz(-1.3915724) q[0];
sx q[0];
rz(-1.916612) q[0];
rz(-pi) q[1];
rz(1.8505356) q[2];
sx q[2];
rz(-1.4609173) q[2];
sx q[2];
rz(2.2834275) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0135611) q[1];
sx q[1];
rz(-1.3741387) q[1];
sx q[1];
rz(0.0081045919) q[1];
rz(-1.891252) q[3];
sx q[3];
rz(-1.081574) q[3];
sx q[3];
rz(1.7699522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.94613218) q[2];
sx q[2];
rz(-1.5798502) q[2];
sx q[2];
rz(1.8411609) q[2];
rz(-0.91222936) q[3];
sx q[3];
rz(-1.2765086) q[3];
sx q[3];
rz(-0.31360489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4881956) q[0];
sx q[0];
rz(-2.6618239) q[0];
sx q[0];
rz(-2.524014) q[0];
rz(-0.081534475) q[1];
sx q[1];
rz(-0.50059861) q[1];
sx q[1];
rz(-0.68731442) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7813635) q[0];
sx q[0];
rz(-1.6183252) q[0];
sx q[0];
rz(-1.8076872) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7527037) q[2];
sx q[2];
rz(-2.8538879) q[2];
sx q[2];
rz(3.1294745) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8246714) q[1];
sx q[1];
rz(-1.8536398) q[1];
sx q[1];
rz(0.43277312) q[1];
x q[2];
rz(-0.048316794) q[3];
sx q[3];
rz(-2.1135277) q[3];
sx q[3];
rz(-0.84191546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.11757892) q[2];
sx q[2];
rz(-1.964183) q[2];
sx q[2];
rz(3.0702316) q[2];
rz(1.9060382) q[3];
sx q[3];
rz(-0.33696285) q[3];
sx q[3];
rz(-0.56628847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5736893) q[0];
sx q[0];
rz(-2.8901143) q[0];
sx q[0];
rz(-0.13791826) q[0];
rz(-2.5714286) q[1];
sx q[1];
rz(-2.0591996) q[1];
sx q[1];
rz(-0.015448419) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.041067657) q[0];
sx q[0];
rz(-0.78162949) q[0];
sx q[0];
rz(-1.1055205) q[0];
rz(-pi) q[1];
rz(-0.76064588) q[2];
sx q[2];
rz(-1.2181251) q[2];
sx q[2];
rz(-1.2404397) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8739104) q[1];
sx q[1];
rz(-1.7355669) q[1];
sx q[1];
rz(-0.49652892) q[1];
x q[2];
rz(-1.3872765) q[3];
sx q[3];
rz(-1.6257111) q[3];
sx q[3];
rz(-2.4650989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.40921679) q[2];
sx q[2];
rz(-0.78745431) q[2];
sx q[2];
rz(2.709205) q[2];
rz(0.89808291) q[3];
sx q[3];
rz(-0.87528527) q[3];
sx q[3];
rz(1.557365) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6569923) q[0];
sx q[0];
rz(-2.0130172) q[0];
sx q[0];
rz(0.76242557) q[0];
rz(2.5905329) q[1];
sx q[1];
rz(-1.3403799) q[1];
sx q[1];
rz(-0.47245477) q[1];
rz(2.77825) q[2];
sx q[2];
rz(-1.6132334) q[2];
sx q[2];
rz(1.2237534) q[2];
rz(-0.47598015) q[3];
sx q[3];
rz(-1.0110216) q[3];
sx q[3];
rz(2.4612751) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
