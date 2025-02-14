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
rz(0.28111464) q[0];
sx q[0];
rz(1.7908362) q[0];
sx q[0];
rz(10.375441) q[0];
rz(-2.9188393) q[1];
sx q[1];
rz(-0.25382257) q[1];
sx q[1];
rz(0.66787994) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53846516) q[0];
sx q[0];
rz(-1.5422945) q[0];
sx q[0];
rz(2.7817581) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4052497) q[2];
sx q[2];
rz(-1.5864463) q[2];
sx q[2];
rz(2.1652255) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4797629) q[1];
sx q[1];
rz(-2.9905028) q[1];
sx q[1];
rz(-2.2041049) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1228074) q[3];
sx q[3];
rz(-2.9308056) q[3];
sx q[3];
rz(2.2700406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5369947) q[2];
sx q[2];
rz(-0.94702417) q[2];
sx q[2];
rz(-2.9941518) q[2];
rz(1.7668746) q[3];
sx q[3];
rz(-1.3876029) q[3];
sx q[3];
rz(-1.714777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.073567063) q[0];
sx q[0];
rz(-1.0140714) q[0];
sx q[0];
rz(-2.9378939) q[0];
rz(3.0286466) q[1];
sx q[1];
rz(-2.4599383) q[1];
sx q[1];
rz(3.122094) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9984922) q[0];
sx q[0];
rz(-1.84853) q[0];
sx q[0];
rz(-0.84622519) q[0];
x q[1];
rz(-1.9505368) q[2];
sx q[2];
rz(-1.4342208) q[2];
sx q[2];
rz(1.4771763) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9128117) q[1];
sx q[1];
rz(-1.8789018) q[1];
sx q[1];
rz(0.56477115) q[1];
x q[2];
rz(-2.2381435) q[3];
sx q[3];
rz(-2.3424582) q[3];
sx q[3];
rz(2.7064221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.514275) q[2];
sx q[2];
rz(-1.5655727) q[2];
sx q[2];
rz(-2.9493098) q[2];
rz(-0.2187885) q[3];
sx q[3];
rz(-1.8407121) q[3];
sx q[3];
rz(-1.7910262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4718276) q[0];
sx q[0];
rz(-0.28194031) q[0];
sx q[0];
rz(-0.48926735) q[0];
rz(1.4504112) q[1];
sx q[1];
rz(-1.1510808) q[1];
sx q[1];
rz(2.6077008) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0301186) q[0];
sx q[0];
rz(-1.5970802) q[0];
sx q[0];
rz(2.8199955) q[0];
rz(-pi) q[1];
rz(0.30812906) q[2];
sx q[2];
rz(-2.9957383) q[2];
sx q[2];
rz(0.062383555) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3870176) q[1];
sx q[1];
rz(-1.1711517) q[1];
sx q[1];
rz(-0.15367457) q[1];
x q[2];
rz(1.3891641) q[3];
sx q[3];
rz(-2.2881621) q[3];
sx q[3];
rz(-2.7434012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.78407946) q[2];
sx q[2];
rz(-1.0369022) q[2];
sx q[2];
rz(-2.4927523) q[2];
rz(2.8283548) q[3];
sx q[3];
rz(-0.99010885) q[3];
sx q[3];
rz(-1.8296957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.677815) q[0];
sx q[0];
rz(-2.2721993) q[0];
sx q[0];
rz(-0.75611269) q[0];
rz(-2.6847878) q[1];
sx q[1];
rz(-1.1242547) q[1];
sx q[1];
rz(-2.4809428) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27460262) q[0];
sx q[0];
rz(-1.9714234) q[0];
sx q[0];
rz(-1.7652049) q[0];
rz(2.9161952) q[2];
sx q[2];
rz(-1.2092012) q[2];
sx q[2];
rz(-2.0394675) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7380437) q[1];
sx q[1];
rz(-1.3894086) q[1];
sx q[1];
rz(-0.84457835) q[1];
rz(-0.33922555) q[3];
sx q[3];
rz(-1.1407099) q[3];
sx q[3];
rz(-0.11232377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.75055355) q[2];
sx q[2];
rz(-1.4918574) q[2];
sx q[2];
rz(0.98461119) q[2];
rz(-1.5084958) q[3];
sx q[3];
rz(-1.0909785) q[3];
sx q[3];
rz(0.35198894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19531798) q[0];
sx q[0];
rz(-2.1066505) q[0];
sx q[0];
rz(-2.4786733) q[0];
rz(1.7341057) q[1];
sx q[1];
rz(-2.2588142) q[1];
sx q[1];
rz(2.1972806) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86032077) q[0];
sx q[0];
rz(-2.4594735) q[0];
sx q[0];
rz(-0.25231326) q[0];
x q[1];
rz(2.495079) q[2];
sx q[2];
rz(-1.0295838) q[2];
sx q[2];
rz(2.3346287) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.27054003) q[1];
sx q[1];
rz(-2.4118198) q[1];
sx q[1];
rz(-0.60128984) q[1];
x q[2];
rz(1.8098894) q[3];
sx q[3];
rz(-1.0464365) q[3];
sx q[3];
rz(-0.6839377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.044109) q[2];
sx q[2];
rz(-2.5554843) q[2];
sx q[2];
rz(0.07494542) q[2];
rz(-0.99503851) q[3];
sx q[3];
rz(-1.1040265) q[3];
sx q[3];
rz(1.9869355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.5465882) q[0];
sx q[0];
rz(-1.1052479) q[0];
sx q[0];
rz(0.30584359) q[0];
rz(-2.2587237) q[1];
sx q[1];
rz(-2.5038033) q[1];
sx q[1];
rz(2.2960704) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5954593) q[0];
sx q[0];
rz(-2.2499878) q[0];
sx q[0];
rz(0.29588411) q[0];
rz(-pi) q[1];
rz(3.1404308) q[2];
sx q[2];
rz(-1.2865304) q[2];
sx q[2];
rz(2.6068316) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9041992) q[1];
sx q[1];
rz(-0.9353726) q[1];
sx q[1];
rz(0.91604427) q[1];
rz(-pi) q[2];
rz(-1.7843397) q[3];
sx q[3];
rz(-0.72493267) q[3];
sx q[3];
rz(1.2675347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.66158295) q[2];
sx q[2];
rz(-0.84605399) q[2];
sx q[2];
rz(-2.0704849) q[2];
rz(0.24460159) q[3];
sx q[3];
rz(-2.1961803) q[3];
sx q[3];
rz(1.5379813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59413183) q[0];
sx q[0];
rz(-0.56532156) q[0];
sx q[0];
rz(-2.2070337) q[0];
rz(2.4681828) q[1];
sx q[1];
rz(-2.5118561) q[1];
sx q[1];
rz(-3.0455132) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.061961147) q[0];
sx q[0];
rz(-1.5966822) q[0];
sx q[0];
rz(-1.4560115) q[0];
x q[1];
rz(0.78676486) q[2];
sx q[2];
rz(-0.25100916) q[2];
sx q[2];
rz(-0.9326084) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2404732) q[1];
sx q[1];
rz(-0.34209004) q[1];
sx q[1];
rz(-0.53804496) q[1];
rz(-pi) q[2];
rz(0.9137093) q[3];
sx q[3];
rz(-1.4944701) q[3];
sx q[3];
rz(0.17599353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.38976321) q[2];
sx q[2];
rz(-1.7182173) q[2];
sx q[2];
rz(-2.5982889) q[2];
rz(3.1144888) q[3];
sx q[3];
rz(-0.47309858) q[3];
sx q[3];
rz(-1.0075587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(-2.2946224) q[0];
sx q[0];
rz(-2.4610418) q[0];
sx q[0];
rz(-0.74817014) q[0];
rz(-1.6698042) q[1];
sx q[1];
rz(-1.4133778) q[1];
sx q[1];
rz(-0.69563037) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75891209) q[0];
sx q[0];
rz(-1.4551388) q[0];
sx q[0];
rz(-1.8291436) q[0];
rz(-pi) q[1];
rz(-1.3514481) q[2];
sx q[2];
rz(-0.92479127) q[2];
sx q[2];
rz(-1.9771119) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7456949) q[1];
sx q[1];
rz(-2.0217748) q[1];
sx q[1];
rz(-1.9459501) q[1];
rz(-pi) q[2];
x q[2];
rz(0.78277709) q[3];
sx q[3];
rz(-2.9874688) q[3];
sx q[3];
rz(-1.085404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2601605) q[2];
sx q[2];
rz(-2.3023534) q[2];
sx q[2];
rz(0.053675573) q[2];
rz(-1.3850877) q[3];
sx q[3];
rz(-1.799492) q[3];
sx q[3];
rz(-0.16689859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78228918) q[0];
sx q[0];
rz(-1.8817236) q[0];
sx q[0];
rz(2.2480929) q[0];
rz(2.9529849) q[1];
sx q[1];
rz(-0.74917561) q[1];
sx q[1];
rz(-2.9395054) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8477488) q[0];
sx q[0];
rz(-1.5408775) q[0];
sx q[0];
rz(2.317753) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6568606) q[2];
sx q[2];
rz(-0.52495365) q[2];
sx q[2];
rz(2.2579272) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.79574078) q[1];
sx q[1];
rz(-0.58064532) q[1];
sx q[1];
rz(-1.1827512) q[1];
rz(-pi) q[2];
rz(2.7703448) q[3];
sx q[3];
rz(-1.8168257) q[3];
sx q[3];
rz(0.061316874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9933652) q[2];
sx q[2];
rz(-1.8348285) q[2];
sx q[2];
rz(0.6692878) q[2];
rz(0.78684849) q[3];
sx q[3];
rz(-1.9653178) q[3];
sx q[3];
rz(-2.440786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65449077) q[0];
sx q[0];
rz(-0.50906068) q[0];
sx q[0];
rz(1.2855726) q[0];
rz(-0.14699832) q[1];
sx q[1];
rz(-1.9963943) q[1];
sx q[1];
rz(0.95050341) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6089382) q[0];
sx q[0];
rz(-1.0883347) q[0];
sx q[0];
rz(-2.5368693) q[0];
x q[1];
rz(0.062531907) q[2];
sx q[2];
rz(-1.034063) q[2];
sx q[2];
rz(-1.5944531) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9626395) q[1];
sx q[1];
rz(-0.8118642) q[1];
sx q[1];
rz(0.78412263) q[1];
rz(-pi) q[2];
rz(-1.4531636) q[3];
sx q[3];
rz(-1.7087987) q[3];
sx q[3];
rz(-1.3641629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6805083) q[2];
sx q[2];
rz(-0.66404873) q[2];
sx q[2];
rz(-0.39592478) q[2];
rz(2.8094273) q[3];
sx q[3];
rz(-1.9931404) q[3];
sx q[3];
rz(-2.240326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.111515) q[0];
sx q[0];
rz(-0.60965309) q[0];
sx q[0];
rz(-1.6994221) q[0];
rz(-0.036399966) q[1];
sx q[1];
rz(-1.2926688) q[1];
sx q[1];
rz(-1.7784437) q[1];
rz(2.8424765) q[2];
sx q[2];
rz(-0.80807297) q[2];
sx q[2];
rz(-0.29812625) q[2];
rz(-2.4717109) q[3];
sx q[3];
rz(-1.5276147) q[3];
sx q[3];
rz(2.0338175) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
