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
rz(1.4242564) q[0];
sx q[0];
rz(-1.2196701) q[0];
sx q[0];
rz(2.5665459) q[0];
rz(-2.9991034) q[1];
sx q[1];
rz(-1.2187076) q[1];
sx q[1];
rz(-2.277318) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40022093) q[0];
sx q[0];
rz(-1.3793886) q[0];
sx q[0];
rz(1.7975397) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2021212) q[2];
sx q[2];
rz(-1.4465904) q[2];
sx q[2];
rz(-0.4865464) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.083347224) q[1];
sx q[1];
rz(-1.5728523) q[1];
sx q[1];
rz(0.53002341) q[1];
x q[2];
rz(2.5745886) q[3];
sx q[3];
rz(-1.1409014) q[3];
sx q[3];
rz(0.56741949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.076866604) q[2];
sx q[2];
rz(-1.6419502) q[2];
sx q[2];
rz(-1.0092674) q[2];
rz(-2.3276954) q[3];
sx q[3];
rz(-0.32865694) q[3];
sx q[3];
rz(0.3391268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-3.0547884) q[0];
sx q[0];
rz(-1.7330994) q[0];
sx q[0];
rz(0.24222294) q[0];
rz(1.9107266) q[1];
sx q[1];
rz(-2.5098398) q[1];
sx q[1];
rz(0.79808527) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4674007) q[0];
sx q[0];
rz(-0.62819203) q[0];
sx q[0];
rz(0.19285658) q[0];
rz(-pi) q[1];
rz(1.961444) q[2];
sx q[2];
rz(-2.7728348) q[2];
sx q[2];
rz(0.40290305) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4646513) q[1];
sx q[1];
rz(-1.2124227) q[1];
sx q[1];
rz(2.6531069) q[1];
x q[2];
rz(1.4105878) q[3];
sx q[3];
rz(-1.2168365) q[3];
sx q[3];
rz(-2.8067971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4855087) q[2];
sx q[2];
rz(-1.069243) q[2];
sx q[2];
rz(-0.83267027) q[2];
rz(-0.63140702) q[3];
sx q[3];
rz(-2.4000945) q[3];
sx q[3];
rz(0.00028636534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1239419) q[0];
sx q[0];
rz(-2.232382) q[0];
sx q[0];
rz(-0.18774524) q[0];
rz(0.84367696) q[1];
sx q[1];
rz(-1.2742821) q[1];
sx q[1];
rz(-2.5086596) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62719856) q[0];
sx q[0];
rz(-2.7717675) q[0];
sx q[0];
rz(-1.373795) q[0];
rz(-pi) q[1];
rz(-2.454921) q[2];
sx q[2];
rz(-1.8455659) q[2];
sx q[2];
rz(-0.49897721) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5714076) q[1];
sx q[1];
rz(-1.270789) q[1];
sx q[1];
rz(-2.6064998) q[1];
x q[2];
rz(-0.12944451) q[3];
sx q[3];
rz(-1.2264381) q[3];
sx q[3];
rz(-2.5710921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8741499) q[2];
sx q[2];
rz(-0.94620693) q[2];
sx q[2];
rz(1.8702742) q[2];
rz(-1.6455796) q[3];
sx q[3];
rz(-1.4964024) q[3];
sx q[3];
rz(0.97150826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2520168) q[0];
sx q[0];
rz(-2.3887964) q[0];
sx q[0];
rz(0.65943199) q[0];
rz(-1.8239498) q[1];
sx q[1];
rz(-1.8023856) q[1];
sx q[1];
rz(-0.79016322) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1135318) q[0];
sx q[0];
rz(-1.4469997) q[0];
sx q[0];
rz(-2.8040299) q[0];
rz(-pi) q[1];
rz(0.019823337) q[2];
sx q[2];
rz(-1.6035282) q[2];
sx q[2];
rz(-2.6898877) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.31393889) q[1];
sx q[1];
rz(-1.0943593) q[1];
sx q[1];
rz(1.6446132) q[1];
rz(-0.42415027) q[3];
sx q[3];
rz(-2.5451535) q[3];
sx q[3];
rz(-3.0354478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0822175) q[2];
sx q[2];
rz(-1.6148115) q[2];
sx q[2];
rz(-0.29602948) q[2];
rz(2.5791903) q[3];
sx q[3];
rz(-0.86807576) q[3];
sx q[3];
rz(0.11399046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.247308) q[0];
sx q[0];
rz(-1.9534651) q[0];
sx q[0];
rz(-0.2071912) q[0];
rz(-2.1098792) q[1];
sx q[1];
rz(-1.9887911) q[1];
sx q[1];
rz(1.4854887) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88865796) q[0];
sx q[0];
rz(-2.2307848) q[0];
sx q[0];
rz(-1.0696017) q[0];
x q[1];
rz(1.0096512) q[2];
sx q[2];
rz(-1.818294) q[2];
sx q[2];
rz(-2.0405318) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7149026) q[1];
sx q[1];
rz(-1.443092) q[1];
sx q[1];
rz(-0.35121484) q[1];
rz(-pi) q[2];
rz(1.6909825) q[3];
sx q[3];
rz(-0.73420364) q[3];
sx q[3];
rz(-0.42163532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1032054) q[2];
sx q[2];
rz(-0.81192553) q[2];
sx q[2];
rz(2.0337598) q[2];
rz(-2.2191018) q[3];
sx q[3];
rz(-1.4100217) q[3];
sx q[3];
rz(0.24423519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6237727) q[0];
sx q[0];
rz(-0.51384846) q[0];
sx q[0];
rz(-2.9129831) q[0];
rz(-0.27433431) q[1];
sx q[1];
rz(-2.3286596) q[1];
sx q[1];
rz(2.6913604) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73295702) q[0];
sx q[0];
rz(-1.3191603) q[0];
sx q[0];
rz(0.77045124) q[0];
rz(-pi) q[1];
rz(-1.3415496) q[2];
sx q[2];
rz(-1.5783797) q[2];
sx q[2];
rz(2.6784865) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.47856146) q[1];
sx q[1];
rz(-1.3197462) q[1];
sx q[1];
rz(0.51086225) q[1];
rz(-pi) q[2];
rz(2.8725876) q[3];
sx q[3];
rz(-1.9719567) q[3];
sx q[3];
rz(1.3022193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7776514) q[2];
sx q[2];
rz(-2.5886017) q[2];
sx q[2];
rz(-3.1138163) q[2];
rz(-2.3751496) q[3];
sx q[3];
rz(-1.981363) q[3];
sx q[3];
rz(0.69507039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.9195093) q[0];
sx q[0];
rz(-1.5064025) q[0];
sx q[0];
rz(-0.62244225) q[0];
rz(0.70603236) q[1];
sx q[1];
rz(-0.89203867) q[1];
sx q[1];
rz(-2.5546254) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76510274) q[0];
sx q[0];
rz(-0.73552948) q[0];
sx q[0];
rz(-0.88659783) q[0];
rz(-pi) q[1];
rz(-0.48611792) q[2];
sx q[2];
rz(-1.3335506) q[2];
sx q[2];
rz(-0.051818661) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7365954) q[1];
sx q[1];
rz(-2.6301503) q[1];
sx q[1];
rz(-2.9963125) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2828279) q[3];
sx q[3];
rz(-1.5629349) q[3];
sx q[3];
rz(-0.9031537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.065757699) q[2];
sx q[2];
rz(-1.4328052) q[2];
sx q[2];
rz(-0.59448057) q[2];
rz(-2.8822656) q[3];
sx q[3];
rz(-1.2649053) q[3];
sx q[3];
rz(1.2172786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.038430564) q[0];
sx q[0];
rz(-0.90534627) q[0];
sx q[0];
rz(2.5627947) q[0];
rz(-1.3102866) q[1];
sx q[1];
rz(-1.2625887) q[1];
sx q[1];
rz(2.328918) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90234251) q[0];
sx q[0];
rz(-2.6391811) q[0];
sx q[0];
rz(-0.30186304) q[0];
rz(-2.8825106) q[2];
sx q[2];
rz(-2.1171085) q[2];
sx q[2];
rz(1.7958884) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.26749565) q[1];
sx q[1];
rz(-2.2363642) q[1];
sx q[1];
rz(2.0366621) q[1];
rz(-pi) q[2];
rz(-2.1401494) q[3];
sx q[3];
rz(-1.1318996) q[3];
sx q[3];
rz(2.2341408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3506713) q[2];
sx q[2];
rz(-1.2988043) q[2];
sx q[2];
rz(2.217963) q[2];
rz(-1.0203429) q[3];
sx q[3];
rz(-2.2444921) q[3];
sx q[3];
rz(-3.0237107) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0608805) q[0];
sx q[0];
rz(-1.4360282) q[0];
sx q[0];
rz(1.4870148) q[0];
rz(-3.0912073) q[1];
sx q[1];
rz(-2.3245508) q[1];
sx q[1];
rz(2.494716) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4934702) q[0];
sx q[0];
rz(-2.1645344) q[0];
sx q[0];
rz(-2.6208682) q[0];
x q[1];
rz(-1.0912618) q[2];
sx q[2];
rz(-2.1796436) q[2];
sx q[2];
rz(0.11907585) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3072435) q[1];
sx q[1];
rz(-1.3615047) q[1];
sx q[1];
rz(2.7268975) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.15022962) q[3];
sx q[3];
rz(-1.7818799) q[3];
sx q[3];
rz(-0.37841132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6369624) q[2];
sx q[2];
rz(-1.5550104) q[2];
sx q[2];
rz(-2.385425) q[2];
rz(1.6715096) q[3];
sx q[3];
rz(-0.78592891) q[3];
sx q[3];
rz(2.3362931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1038372) q[0];
sx q[0];
rz(-1.2498195) q[0];
sx q[0];
rz(-2.0334429) q[0];
rz(0.26421079) q[1];
sx q[1];
rz(-1.4690396) q[1];
sx q[1];
rz(-2.5392551) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7450934) q[0];
sx q[0];
rz(-1.709793) q[0];
sx q[0];
rz(-1.2953912) q[0];
x q[1];
rz(-1.7246288) q[2];
sx q[2];
rz(-1.7454426) q[2];
sx q[2];
rz(-1.5690116) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4385637) q[1];
sx q[1];
rz(-1.5224341) q[1];
sx q[1];
rz(2.2579184) q[1];
rz(-pi) q[2];
rz(0.29903166) q[3];
sx q[3];
rz(-2.4686738) q[3];
sx q[3];
rz(2.3367019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5084874) q[2];
sx q[2];
rz(-0.54163951) q[2];
sx q[2];
rz(-2.2701021) q[2];
rz(1.2095215) q[3];
sx q[3];
rz(-2.8612374) q[3];
sx q[3];
rz(2.5476294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4393944) q[0];
sx q[0];
rz(-1.3855423) q[0];
sx q[0];
rz(2.6227797) q[0];
rz(-2.5595472) q[1];
sx q[1];
rz(-2.0379635) q[1];
sx q[1];
rz(1.4720974) q[1];
rz(-1.1956566) q[2];
sx q[2];
rz(-2.2715501) q[2];
sx q[2];
rz(3.0987433) q[2];
rz(-0.20584917) q[3];
sx q[3];
rz(-2.4424853) q[3];
sx q[3];
rz(1.427099) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
