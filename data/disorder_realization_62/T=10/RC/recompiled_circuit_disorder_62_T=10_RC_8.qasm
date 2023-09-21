OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2405038) q[0];
sx q[0];
rz(-2.9641889) q[0];
sx q[0];
rz(2.0071964) q[0];
rz(-1.9534684) q[1];
sx q[1];
rz(-1.0367353) q[1];
sx q[1];
rz(-2.477975) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7952607) q[0];
sx q[0];
rz(-1.5471317) q[0];
sx q[0];
rz(2.6803826) q[0];
rz(-pi) q[1];
rz(0.77779777) q[2];
sx q[2];
rz(-1.3133089) q[2];
sx q[2];
rz(0.68629005) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.44890468) q[1];
sx q[1];
rz(-1.2346039) q[1];
sx q[1];
rz(0.2775788) q[1];
rz(0.10856467) q[3];
sx q[3];
rz(-1.3285471) q[3];
sx q[3];
rz(0.066075174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.87876451) q[2];
sx q[2];
rz(-0.44439134) q[2];
sx q[2];
rz(3.0905241) q[2];
rz(0.55705327) q[3];
sx q[3];
rz(-2.3414108) q[3];
sx q[3];
rz(-1.5548271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58650815) q[0];
sx q[0];
rz(-2.3095135) q[0];
sx q[0];
rz(-2.5449975) q[0];
rz(0.82582981) q[1];
sx q[1];
rz(-1.700371) q[1];
sx q[1];
rz(-1.2260431) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7694089) q[0];
sx q[0];
rz(-2.6499977) q[0];
sx q[0];
rz(-0.43222897) q[0];
x q[1];
rz(-0.82310279) q[2];
sx q[2];
rz(-0.96379333) q[2];
sx q[2];
rz(-0.58128231) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5920168) q[1];
sx q[1];
rz(-1.5703652) q[1];
sx q[1];
rz(1.7906584) q[1];
rz(-pi) q[2];
rz(-1.5561043) q[3];
sx q[3];
rz(-2.5689295) q[3];
sx q[3];
rz(0.0088012561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3423959) q[2];
sx q[2];
rz(-1.9689955) q[2];
sx q[2];
rz(-2.8105695) q[2];
rz(2.3349169) q[3];
sx q[3];
rz(-1.2253864) q[3];
sx q[3];
rz(-1.4276918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9817292) q[0];
sx q[0];
rz(-1.9112497) q[0];
sx q[0];
rz(1.249041) q[0];
rz(-3.0535835) q[1];
sx q[1];
rz(-1.1227337) q[1];
sx q[1];
rz(2.1121315) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7824161) q[0];
sx q[0];
rz(-2.4826907) q[0];
sx q[0];
rz(1.2078148) q[0];
rz(-pi) q[1];
rz(2.4918409) q[2];
sx q[2];
rz(-1.3403112) q[2];
sx q[2];
rz(0.31313716) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.50476915) q[1];
sx q[1];
rz(-1.4154589) q[1];
sx q[1];
rz(-1.4298646) q[1];
rz(-1.6920407) q[3];
sx q[3];
rz(-2.3070934) q[3];
sx q[3];
rz(2.8672225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.007894667) q[2];
sx q[2];
rz(-1.4174613) q[2];
sx q[2];
rz(-0.50764817) q[2];
rz(1.3890022) q[3];
sx q[3];
rz(-2.8116083) q[3];
sx q[3];
rz(1.1631789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4841109) q[0];
sx q[0];
rz(-2.8634475) q[0];
sx q[0];
rz(-1.5959651) q[0];
rz(1.0428628) q[1];
sx q[1];
rz(-1.9680126) q[1];
sx q[1];
rz(1.5159336) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8221178) q[0];
sx q[0];
rz(-1.428077) q[0];
sx q[0];
rz(-0.6832173) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0117202) q[2];
sx q[2];
rz(-1.498073) q[2];
sx q[2];
rz(-0.59414547) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.792946) q[1];
sx q[1];
rz(-0.9325087) q[1];
sx q[1];
rz(-1.9116886) q[1];
rz(-0.32228542) q[3];
sx q[3];
rz(-0.61338131) q[3];
sx q[3];
rz(1.8005288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.76379124) q[2];
sx q[2];
rz(-2.2968473) q[2];
sx q[2];
rz(-1.2949004) q[2];
rz(0.14136782) q[3];
sx q[3];
rz(-2.5989792) q[3];
sx q[3];
rz(-2.1550089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7871053) q[0];
sx q[0];
rz(-1.961345) q[0];
sx q[0];
rz(-1.6217344) q[0];
rz(2.5095818) q[1];
sx q[1];
rz(-2.4193587) q[1];
sx q[1];
rz(-2.246726) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1631854) q[0];
sx q[0];
rz(-0.95046959) q[0];
sx q[0];
rz(-1.0582256) q[0];
rz(-pi) q[1];
rz(-1.0084444) q[2];
sx q[2];
rz(-2.214606) q[2];
sx q[2];
rz(-1.2959727) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.289031) q[1];
sx q[1];
rz(-0.82419306) q[1];
sx q[1];
rz(2.186071) q[1];
rz(-pi) q[2];
rz(-1.5785061) q[3];
sx q[3];
rz(-1.8572516) q[3];
sx q[3];
rz(0.25659284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.52508369) q[2];
sx q[2];
rz(-0.36642763) q[2];
sx q[2];
rz(1.9990702) q[2];
rz(2.3948495) q[3];
sx q[3];
rz(-1.2544422) q[3];
sx q[3];
rz(2.0660627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.557945) q[0];
sx q[0];
rz(-0.32145158) q[0];
sx q[0];
rz(-1.3775795) q[0];
rz(0.47239834) q[1];
sx q[1];
rz(-2.6230085) q[1];
sx q[1];
rz(0.46498743) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87100055) q[0];
sx q[0];
rz(-2.5809079) q[0];
sx q[0];
rz(0.90765783) q[0];
x q[1];
rz(2.9372413) q[2];
sx q[2];
rz(-0.34826476) q[2];
sx q[2];
rz(1.5262926) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7558414) q[1];
sx q[1];
rz(-2.5045536) q[1];
sx q[1];
rz(-1.4904651) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3607849) q[3];
sx q[3];
rz(-0.87267733) q[3];
sx q[3];
rz(0.39892808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.650699) q[2];
sx q[2];
rz(-1.2630254) q[2];
sx q[2];
rz(-0.4450376) q[2];
rz(-0.93368357) q[3];
sx q[3];
rz(-1.4327587) q[3];
sx q[3];
rz(-2.8745108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12748195) q[0];
sx q[0];
rz(-1.575379) q[0];
sx q[0];
rz(1.4200462) q[0];
rz(0.02380112) q[1];
sx q[1];
rz(-0.61444608) q[1];
sx q[1];
rz(-0.15596095) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0359081) q[0];
sx q[0];
rz(-1.3996291) q[0];
sx q[0];
rz(-2.8842584) q[0];
x q[1];
rz(1.9975125) q[2];
sx q[2];
rz(-2.0442171) q[2];
sx q[2];
rz(-2.0063426) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1289039) q[1];
sx q[1];
rz(-1.745599) q[1];
sx q[1];
rz(0.054794475) q[1];
rz(-pi) q[2];
rz(0.74219269) q[3];
sx q[3];
rz(-2.4626197) q[3];
sx q[3];
rz(1.2364482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.069313958) q[2];
sx q[2];
rz(-0.57702714) q[2];
sx q[2];
rz(-2.0689266) q[2];
rz(-2.8159451) q[3];
sx q[3];
rz(-1.9653392) q[3];
sx q[3];
rz(1.6252888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9496562) q[0];
sx q[0];
rz(-2.7392116) q[0];
sx q[0];
rz(0.33777133) q[0];
rz(1.0900963) q[1];
sx q[1];
rz(-2.1665159) q[1];
sx q[1];
rz(2.8930194) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8817431) q[0];
sx q[0];
rz(-1.8200841) q[0];
sx q[0];
rz(2.6134406) q[0];
rz(-pi) q[1];
x q[1];
rz(0.48034251) q[2];
sx q[2];
rz(-1.1104274) q[2];
sx q[2];
rz(-0.72887052) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3358811) q[1];
sx q[1];
rz(-2.1599249) q[1];
sx q[1];
rz(-0.59405234) q[1];
x q[2];
rz(1.9931273) q[3];
sx q[3];
rz(-0.18581192) q[3];
sx q[3];
rz(-1.8728135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8481855) q[2];
sx q[2];
rz(-0.53331393) q[2];
sx q[2];
rz(3.0548813) q[2];
rz(2.6596206) q[3];
sx q[3];
rz(-0.50589365) q[3];
sx q[3];
rz(-1.988525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6329353) q[0];
sx q[0];
rz(-1.0077227) q[0];
sx q[0];
rz(0.28276643) q[0];
rz(2.4400318) q[1];
sx q[1];
rz(-0.82156721) q[1];
sx q[1];
rz(1.3185906) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1222262) q[0];
sx q[0];
rz(-1.1302233) q[0];
sx q[0];
rz(-3.0526572) q[0];
rz(-pi) q[1];
rz(2.4852072) q[2];
sx q[2];
rz(-2.0771386) q[2];
sx q[2];
rz(2.2040747) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.56837592) q[1];
sx q[1];
rz(-1.425256) q[1];
sx q[1];
rz(-2.2578866) q[1];
rz(0.75026476) q[3];
sx q[3];
rz(-1.2816458) q[3];
sx q[3];
rz(-2.9671448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41137722) q[2];
sx q[2];
rz(-0.29835478) q[2];
sx q[2];
rz(0.39917699) q[2];
rz(-2.2579851) q[3];
sx q[3];
rz(-1.4727605) q[3];
sx q[3];
rz(-1.9201027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76319641) q[0];
sx q[0];
rz(-0.77532399) q[0];
sx q[0];
rz(-1.5426853) q[0];
rz(1.0653161) q[1];
sx q[1];
rz(-2.1677446) q[1];
sx q[1];
rz(-1.261196) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3189145) q[0];
sx q[0];
rz(-0.42352522) q[0];
sx q[0];
rz(2.8828997) q[0];
rz(-pi) q[1];
x q[1];
rz(0.59826675) q[2];
sx q[2];
rz(-1.5524128) q[2];
sx q[2];
rz(2.604904) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0295769) q[1];
sx q[1];
rz(-0.50260168) q[1];
sx q[1];
rz(0.13346787) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4114755) q[3];
sx q[3];
rz(-2.2485002) q[3];
sx q[3];
rz(2.2417828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.5579055) q[2];
sx q[2];
rz(-0.62290278) q[2];
sx q[2];
rz(0.56979257) q[2];
rz(1.2184881) q[3];
sx q[3];
rz(-2.4945538) q[3];
sx q[3];
rz(-0.56263721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8284843) q[0];
sx q[0];
rz(-0.88043558) q[0];
sx q[0];
rz(1.2783929) q[0];
rz(2.5333511) q[1];
sx q[1];
rz(-0.47641644) q[1];
sx q[1];
rz(0.48412916) q[1];
rz(-3.1153395) q[2];
sx q[2];
rz(-2.142341) q[2];
sx q[2];
rz(-0.061948902) q[2];
rz(-1.7096056) q[3];
sx q[3];
rz(-2.1574253) q[3];
sx q[3];
rz(0.16711259) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];