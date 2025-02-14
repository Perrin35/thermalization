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
rz(0.87585706) q[0];
sx q[0];
rz(2.1729204) q[0];
sx q[0];
rz(7.2090413) q[0];
rz(0.61697382) q[1];
sx q[1];
rz(-0.66520912) q[1];
sx q[1];
rz(-1.8170504) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5228341) q[0];
sx q[0];
rz(-1.6346187) q[0];
sx q[0];
rz(-0.60019779) q[0];
rz(-pi) q[1];
rz(1.6683031) q[2];
sx q[2];
rz(-2.153399) q[2];
sx q[2];
rz(-2.1511457) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1086667) q[1];
sx q[1];
rz(-0.98331988) q[1];
sx q[1];
rz(-1.0125748) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.76530568) q[3];
sx q[3];
rz(-2.0767143) q[3];
sx q[3];
rz(1.7047061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2609743) q[2];
sx q[2];
rz(-1.7558492) q[2];
sx q[2];
rz(0.79208881) q[2];
rz(2.2657307) q[3];
sx q[3];
rz(-3.0328817) q[3];
sx q[3];
rz(-0.66332269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51139128) q[0];
sx q[0];
rz(-1.8237317) q[0];
sx q[0];
rz(0.9414916) q[0];
rz(-0.9990274) q[1];
sx q[1];
rz(-0.92728725) q[1];
sx q[1];
rz(0.082854465) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0339399) q[0];
sx q[0];
rz(-1.6710738) q[0];
sx q[0];
rz(-2.4058008) q[0];
rz(-pi) q[1];
rz(-2.8993334) q[2];
sx q[2];
rz(-1.2238811) q[2];
sx q[2];
rz(-2.0478947) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8585457) q[1];
sx q[1];
rz(-1.0632005) q[1];
sx q[1];
rz(-2.0920442) q[1];
x q[2];
rz(0.64506809) q[3];
sx q[3];
rz(-1.0655902) q[3];
sx q[3];
rz(-0.92407214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2443709) q[2];
sx q[2];
rz(-1.3542078) q[2];
sx q[2];
rz(1.8841891) q[2];
rz(0.8463549) q[3];
sx q[3];
rz(-2.9823163) q[3];
sx q[3];
rz(0.67811051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9427247) q[0];
sx q[0];
rz(-1.4582448) q[0];
sx q[0];
rz(0.22920907) q[0];
rz(-0.21408679) q[1];
sx q[1];
rz(-2.2702718) q[1];
sx q[1];
rz(2.7600938) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7575398) q[0];
sx q[0];
rz(-2.4017576) q[0];
sx q[0];
rz(0.53348855) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7842641) q[2];
sx q[2];
rz(-0.74356438) q[2];
sx q[2];
rz(-1.5539757) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2974071) q[1];
sx q[1];
rz(-1.5995889) q[1];
sx q[1];
rz(1.5178568) q[1];
rz(-0.34708628) q[3];
sx q[3];
rz(-1.5659837) q[3];
sx q[3];
rz(0.50932415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7233589) q[2];
sx q[2];
rz(-2.8853719) q[2];
sx q[2];
rz(-2.5733433) q[2];
rz(1.4189643) q[3];
sx q[3];
rz(-1.7122372) q[3];
sx q[3];
rz(-1.0770146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-3.1133465) q[0];
sx q[0];
rz(-1.1320817) q[0];
sx q[0];
rz(2.8908253) q[0];
rz(-3.057042) q[1];
sx q[1];
rz(-1.0944347) q[1];
sx q[1];
rz(-1.9409404) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47287175) q[0];
sx q[0];
rz(-0.78470147) q[0];
sx q[0];
rz(1.0502104) q[0];
rz(1.0798825) q[2];
sx q[2];
rz(-2.6497095) q[2];
sx q[2];
rz(-1.5747923) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6509962) q[1];
sx q[1];
rz(-0.8621434) q[1];
sx q[1];
rz(-1.7017822) q[1];
rz(-pi) q[2];
rz(2.8765466) q[3];
sx q[3];
rz(-2.0816021) q[3];
sx q[3];
rz(0.068451133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.67479977) q[2];
sx q[2];
rz(-1.1226706) q[2];
sx q[2];
rz(-1.2910845) q[2];
rz(0.42168266) q[3];
sx q[3];
rz(-2.6899874) q[3];
sx q[3];
rz(-1.565056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52294937) q[0];
sx q[0];
rz(-2.4186501) q[0];
sx q[0];
rz(1.500754) q[0];
rz(-0.067642637) q[1];
sx q[1];
rz(-1.5875971) q[1];
sx q[1];
rz(2.0134451) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90445176) q[0];
sx q[0];
rz(-1.6078976) q[0];
sx q[0];
rz(-0.18795975) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7893547) q[2];
sx q[2];
rz(-2.0697429) q[2];
sx q[2];
rz(1.3832472) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6349244) q[1];
sx q[1];
rz(-1.4574058) q[1];
sx q[1];
rz(3.0874814) q[1];
rz(-pi) q[2];
rz(-0.23388548) q[3];
sx q[3];
rz(-1.8429759) q[3];
sx q[3];
rz(-0.50258499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1077914) q[2];
sx q[2];
rz(-2.0267603) q[2];
sx q[2];
rz(2.4647253) q[2];
rz(-0.90773165) q[3];
sx q[3];
rz(-1.2264484) q[3];
sx q[3];
rz(-2.7270253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9229729) q[0];
sx q[0];
rz(-0.75160471) q[0];
sx q[0];
rz(2.3722017) q[0];
rz(-1.9006624) q[1];
sx q[1];
rz(-2.2878094) q[1];
sx q[1];
rz(-0.21533899) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0701987) q[0];
sx q[0];
rz(-2.5339273) q[0];
sx q[0];
rz(1.9032701) q[0];
rz(0.017069503) q[2];
sx q[2];
rz(-1.2400377) q[2];
sx q[2];
rz(-2.0278553) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.78241888) q[1];
sx q[1];
rz(-2.4451588) q[1];
sx q[1];
rz(2.7553808) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.69778473) q[3];
sx q[3];
rz(-1.3078193) q[3];
sx q[3];
rz(1.3828204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.26663366) q[2];
sx q[2];
rz(-1.5077488) q[2];
sx q[2];
rz(-0.61666644) q[2];
rz(2.941361) q[3];
sx q[3];
rz(-2.4112406) q[3];
sx q[3];
rz(-2.0088137) q[3];
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
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.013833372) q[0];
sx q[0];
rz(-2.5795689) q[0];
sx q[0];
rz(3.1264937) q[0];
rz(-3.0191782) q[1];
sx q[1];
rz(-1.2950803) q[1];
sx q[1];
rz(-1.0282358) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6262344) q[0];
sx q[0];
rz(-1.5152103) q[0];
sx q[0];
rz(0.65599156) q[0];
x q[1];
rz(1.2790658) q[2];
sx q[2];
rz(-2.5030067) q[2];
sx q[2];
rz(1.8420458) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4738481) q[1];
sx q[1];
rz(-1.5989002) q[1];
sx q[1];
rz(2.7098118) q[1];
rz(-2.4373152) q[3];
sx q[3];
rz(-2.7368494) q[3];
sx q[3];
rz(-1.2621439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5369109) q[2];
sx q[2];
rz(-1.9313507) q[2];
sx q[2];
rz(-0.49986419) q[2];
rz(-2.8258421) q[3];
sx q[3];
rz(-0.67086589) q[3];
sx q[3];
rz(-1.0189112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.71378088) q[0];
sx q[0];
rz(-1.918387) q[0];
sx q[0];
rz(2.1977303) q[0];
rz(1.7367412) q[1];
sx q[1];
rz(-1.8753139) q[1];
sx q[1];
rz(-0.92165438) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6674468) q[0];
sx q[0];
rz(-1.6255127) q[0];
sx q[0];
rz(-2.823141) q[0];
x q[1];
rz(-0.87041847) q[2];
sx q[2];
rz(-0.38258115) q[2];
sx q[2];
rz(-2.062881) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6347305) q[1];
sx q[1];
rz(-0.94351879) q[1];
sx q[1];
rz(-0.92232134) q[1];
x q[2];
rz(-1.5783903) q[3];
sx q[3];
rz(-0.51285997) q[3];
sx q[3];
rz(0.22554413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1450119) q[2];
sx q[2];
rz(-1.8700446) q[2];
sx q[2];
rz(-1.5550295) q[2];
rz(-1.0541213) q[3];
sx q[3];
rz(-1.328822) q[3];
sx q[3];
rz(1.740295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9911875) q[0];
sx q[0];
rz(-1.0973955) q[0];
sx q[0];
rz(-2.4990668) q[0];
rz(1.4379028) q[1];
sx q[1];
rz(-2.3846886) q[1];
sx q[1];
rz(-2.9978851) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7724975) q[0];
sx q[0];
rz(-1.0155627) q[0];
sx q[0];
rz(1.9291718) q[0];
x q[1];
rz(-2.2238067) q[2];
sx q[2];
rz(-2.3812903) q[2];
sx q[2];
rz(-2.0298634) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9160812) q[1];
sx q[1];
rz(-2.4298432) q[1];
sx q[1];
rz(2.7307778) q[1];
rz(-2.4293184) q[3];
sx q[3];
rz(-2.5872719) q[3];
sx q[3];
rz(-1.8795151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5370499) q[2];
sx q[2];
rz(-1.1979251) q[2];
sx q[2];
rz(-2.878888) q[2];
rz(0.25775868) q[3];
sx q[3];
rz(-2.7596605) q[3];
sx q[3];
rz(-2.2885382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8464552) q[0];
sx q[0];
rz(-1.6520123) q[0];
sx q[0];
rz(1.2204131) q[0];
rz(-2.659761) q[1];
sx q[1];
rz(-2.316663) q[1];
sx q[1];
rz(0.89734546) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37384826) q[0];
sx q[0];
rz(-1.6929099) q[0];
sx q[0];
rz(2.6227345) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4608896) q[2];
sx q[2];
rz(-1.8898003) q[2];
sx q[2];
rz(1.4665335) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7792977) q[1];
sx q[1];
rz(-0.26302281) q[1];
sx q[1];
rz(-1.7630601) q[1];
rz(-pi) q[2];
rz(0.10372377) q[3];
sx q[3];
rz(-2.1770766) q[3];
sx q[3];
rz(2.1192354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.64409488) q[2];
sx q[2];
rz(-2.2150025) q[2];
sx q[2];
rz(3.0268055) q[2];
rz(-2.3980906) q[3];
sx q[3];
rz(-1.8758352) q[3];
sx q[3];
rz(-0.44873294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7262309) q[0];
sx q[0];
rz(-0.76114934) q[0];
sx q[0];
rz(2.1834955) q[0];
rz(1.6356946) q[1];
sx q[1];
rz(-1.1264569) q[1];
sx q[1];
rz(1.1503848) q[1];
rz(-2.220357) q[2];
sx q[2];
rz(-0.9365281) q[2];
sx q[2];
rz(-2.3021883) q[2];
rz(-2.791849) q[3];
sx q[3];
rz(-2.6833224) q[3];
sx q[3];
rz(0.53862017) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
