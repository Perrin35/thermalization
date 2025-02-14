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
rz(-0.32831353) q[0];
sx q[0];
rz(-1.1348731) q[0];
sx q[0];
rz(-2.5757134) q[0];
rz(-2.8973051) q[1];
sx q[1];
rz(-1.7841508) q[1];
sx q[1];
rz(2.607333) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.256828) q[0];
sx q[0];
rz(-0.74640154) q[0];
sx q[0];
rz(-1.9749667) q[0];
rz(-pi) q[1];
x q[1];
rz(0.59492971) q[2];
sx q[2];
rz(-1.8244566) q[2];
sx q[2];
rz(2.4640535) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6837343) q[1];
sx q[1];
rz(-2.8401655) q[1];
sx q[1];
rz(3.0617759) q[1];
rz(-pi) q[2];
x q[2];
rz(2.684428) q[3];
sx q[3];
rz(-1.346753) q[3];
sx q[3];
rz(-2.8347297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.013022097) q[2];
sx q[2];
rz(-2.8910922) q[2];
sx q[2];
rz(2.8924083) q[2];
rz(1.3439882) q[3];
sx q[3];
rz(-1.598571) q[3];
sx q[3];
rz(-0.070076076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4942112) q[0];
sx q[0];
rz(-1.9152315) q[0];
sx q[0];
rz(0.98980728) q[0];
rz(-0.79065943) q[1];
sx q[1];
rz(-0.76690563) q[1];
sx q[1];
rz(-2.8299423) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30325952) q[0];
sx q[0];
rz(-2.0796392) q[0];
sx q[0];
rz(-0.80755393) q[0];
rz(1.0490408) q[2];
sx q[2];
rz(-0.21460303) q[2];
sx q[2];
rz(1.8282481) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6271882) q[1];
sx q[1];
rz(-1.3986821) q[1];
sx q[1];
rz(-1.3283511) q[1];
rz(-pi) q[2];
rz(-1.1053322) q[3];
sx q[3];
rz(-0.66879818) q[3];
sx q[3];
rz(-1.5580524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1404169) q[2];
sx q[2];
rz(-1.1602465) q[2];
sx q[2];
rz(-2.1902093) q[2];
rz(2.9133255) q[3];
sx q[3];
rz(-2.164866) q[3];
sx q[3];
rz(-1.2736646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9839639) q[0];
sx q[0];
rz(-1.552859) q[0];
sx q[0];
rz(0.11446318) q[0];
rz(1.0022256) q[1];
sx q[1];
rz(-0.4267692) q[1];
sx q[1];
rz(2.9797629) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67437088) q[0];
sx q[0];
rz(-0.97839744) q[0];
sx q[0];
rz(-0.47426275) q[0];
rz(-pi) q[1];
rz(-0.21185565) q[2];
sx q[2];
rz(-2.7447034) q[2];
sx q[2];
rz(-0.8241764) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.964965) q[1];
sx q[1];
rz(-1.3313659) q[1];
sx q[1];
rz(2.5188732) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.80064818) q[3];
sx q[3];
rz(-1.8654239) q[3];
sx q[3];
rz(1.1931452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3942922) q[2];
sx q[2];
rz(-2.6960399) q[2];
sx q[2];
rz(3.0008345) q[2];
rz(-0.55177871) q[3];
sx q[3];
rz(-1.2740302) q[3];
sx q[3];
rz(-0.75132918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20525876) q[0];
sx q[0];
rz(-2.8755499) q[0];
sx q[0];
rz(-2.4667013) q[0];
rz(-0.15448013) q[1];
sx q[1];
rz(-1.7325502) q[1];
sx q[1];
rz(2.6836269) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5602011) q[0];
sx q[0];
rz(-2.0580252) q[0];
sx q[0];
rz(1.1283406) q[0];
x q[1];
rz(-2.1764917) q[2];
sx q[2];
rz(-2.6332246) q[2];
sx q[2];
rz(-1.3077715) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.47487101) q[1];
sx q[1];
rz(-2.7242766) q[1];
sx q[1];
rz(-0.46837193) q[1];
x q[2];
rz(-0.19019698) q[3];
sx q[3];
rz(-2.6256997) q[3];
sx q[3];
rz(-2.7864393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.13350479) q[2];
sx q[2];
rz(-2.4613481) q[2];
sx q[2];
rz(2.8221455) q[2];
rz(-1.8114629) q[3];
sx q[3];
rz(-1.0165241) q[3];
sx q[3];
rz(0.5602079) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1595681) q[0];
sx q[0];
rz(-2.0052795) q[0];
sx q[0];
rz(-1.9200448) q[0];
rz(-0.23280652) q[1];
sx q[1];
rz(-2.5722952) q[1];
sx q[1];
rz(0.98308841) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.074551) q[0];
sx q[0];
rz(-2.3877151) q[0];
sx q[0];
rz(-0.75410944) q[0];
rz(-pi) q[1];
rz(-0.95011656) q[2];
sx q[2];
rz(-1.663637) q[2];
sx q[2];
rz(-0.29618057) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.59683387) q[1];
sx q[1];
rz(-1.7167517) q[1];
sx q[1];
rz(0.18827855) q[1];
x q[2];
rz(1.1484723) q[3];
sx q[3];
rz(-1.7230801) q[3];
sx q[3];
rz(-1.9773169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.38272196) q[2];
sx q[2];
rz(-1.4843586) q[2];
sx q[2];
rz(0.37528428) q[2];
rz(-2.0659857) q[3];
sx q[3];
rz(-2.7552102) q[3];
sx q[3];
rz(2.6066499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7406834) q[0];
sx q[0];
rz(-1.704819) q[0];
sx q[0];
rz(1.4373454) q[0];
rz(3.0628693) q[1];
sx q[1];
rz(-1.3767786) q[1];
sx q[1];
rz(-2.5934503) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21273147) q[0];
sx q[0];
rz(-2.2240046) q[0];
sx q[0];
rz(-1.8381339) q[0];
rz(-pi) q[1];
rz(2.5896543) q[2];
sx q[2];
rz(-1.1480398) q[2];
sx q[2];
rz(-1.960159) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2386475) q[1];
sx q[1];
rz(-0.79552197) q[1];
sx q[1];
rz(2.5832157) q[1];
x q[2];
rz(0.94977149) q[3];
sx q[3];
rz(-0.91663893) q[3];
sx q[3];
rz(-1.1266687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6846201) q[2];
sx q[2];
rz(-0.62070864) q[2];
sx q[2];
rz(-0.6753298) q[2];
rz(-1.5843377) q[3];
sx q[3];
rz(-1.3074338) q[3];
sx q[3];
rz(2.5244782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-2.5976582) q[0];
sx q[0];
rz(-1.958853) q[0];
sx q[0];
rz(-2.6455998) q[0];
rz(1.687382) q[1];
sx q[1];
rz(-1.0554689) q[1];
sx q[1];
rz(-1.1167663) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1710715) q[0];
sx q[0];
rz(-1.9200033) q[0];
sx q[0];
rz(1.1961351) q[0];
rz(-pi) q[1];
x q[1];
rz(0.26979763) q[2];
sx q[2];
rz(-2.0922263) q[2];
sx q[2];
rz(-0.21873378) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.36193725) q[1];
sx q[1];
rz(-1.2314267) q[1];
sx q[1];
rz(1.9091868) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12567606) q[3];
sx q[3];
rz(-2.6824441) q[3];
sx q[3];
rz(2.7530376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.83069673) q[2];
sx q[2];
rz(-2.1572957) q[2];
sx q[2];
rz(2.7867219) q[2];
rz(0.76534671) q[3];
sx q[3];
rz(-0.86330515) q[3];
sx q[3];
rz(-0.089546831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4527721) q[0];
sx q[0];
rz(-1.4790164) q[0];
sx q[0];
rz(0.76559693) q[0];
rz(2.2701524) q[1];
sx q[1];
rz(-2.3759418) q[1];
sx q[1];
rz(-1.3227468) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.046503456) q[0];
sx q[0];
rz(-0.26581598) q[0];
sx q[0];
rz(2.9429716) q[0];
rz(0.23643146) q[2];
sx q[2];
rz(-2.2364103) q[2];
sx q[2];
rz(-2.1588391) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.31963667) q[1];
sx q[1];
rz(-1.0624891) q[1];
sx q[1];
rz(1.5896569) q[1];
rz(-pi) q[2];
rz(0.83006184) q[3];
sx q[3];
rz(-0.30127159) q[3];
sx q[3];
rz(0.64727441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4837997) q[2];
sx q[2];
rz(-1.0427661) q[2];
sx q[2];
rz(-0.22171177) q[2];
rz(2.0486369) q[3];
sx q[3];
rz(-1.4128128) q[3];
sx q[3];
rz(2.4634585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.073864989) q[0];
sx q[0];
rz(-1.8949969) q[0];
sx q[0];
rz(-2.8644323) q[0];
rz(0.74835888) q[1];
sx q[1];
rz(-0.44735083) q[1];
sx q[1];
rz(1.998273) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3353334) q[0];
sx q[0];
rz(-1.8758095) q[0];
sx q[0];
rz(-0.2103896) q[0];
x q[1];
rz(2.9306445) q[2];
sx q[2];
rz(-2.2648513) q[2];
sx q[2];
rz(1.8927434) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.34830293) q[1];
sx q[1];
rz(-1.7277246) q[1];
sx q[1];
rz(1.8918397) q[1];
rz(-pi) q[2];
rz(0.36181776) q[3];
sx q[3];
rz(-1.4103977) q[3];
sx q[3];
rz(-3.0830864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2271759) q[2];
sx q[2];
rz(-1.02355) q[2];
sx q[2];
rz(0.048967036) q[2];
rz(2.07552) q[3];
sx q[3];
rz(-0.62838307) q[3];
sx q[3];
rz(-0.15570417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0186036) q[0];
sx q[0];
rz(-1.5276696) q[0];
sx q[0];
rz(0.019512026) q[0];
rz(-2.3628269) q[1];
sx q[1];
rz(-0.6929144) q[1];
sx q[1];
rz(-1.5787517) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39017856) q[0];
sx q[0];
rz(-1.4407684) q[0];
sx q[0];
rz(0.089827248) q[0];
rz(-pi) q[1];
rz(1.3124002) q[2];
sx q[2];
rz(-2.5157197) q[2];
sx q[2];
rz(-1.3385119) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9924888) q[1];
sx q[1];
rz(-1.1380151) q[1];
sx q[1];
rz(0.91723587) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3046523) q[3];
sx q[3];
rz(-1.4599953) q[3];
sx q[3];
rz(1.5594123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.19266263) q[2];
sx q[2];
rz(-0.86343416) q[2];
sx q[2];
rz(-2.9760402) q[2];
rz(-2.5340269) q[3];
sx q[3];
rz(-1.8568042) q[3];
sx q[3];
rz(2.6732388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7623357) q[0];
sx q[0];
rz(-0.48185928) q[0];
sx q[0];
rz(3.0961105) q[0];
rz(1.634585) q[1];
sx q[1];
rz(-1.6804809) q[1];
sx q[1];
rz(1.5483821) q[1];
rz(1.7138557) q[2];
sx q[2];
rz(-0.80052994) q[2];
sx q[2];
rz(1.1442727) q[2];
rz(2.8589917) q[3];
sx q[3];
rz(-0.62487124) q[3];
sx q[3];
rz(1.665653) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
