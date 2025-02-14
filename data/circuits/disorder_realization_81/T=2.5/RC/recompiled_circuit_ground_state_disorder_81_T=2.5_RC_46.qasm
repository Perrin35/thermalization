OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9658907) q[0];
sx q[0];
rz(-1.8718636) q[0];
sx q[0];
rz(-1.9667392) q[0];
rz(2.0570316) q[1];
sx q[1];
rz(4.6133981) q[1];
sx q[1];
rz(2.5785799) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63516894) q[0];
sx q[0];
rz(-0.20316589) q[0];
sx q[0];
rz(-2.8750393) q[0];
rz(-2.6649551) q[2];
sx q[2];
rz(-0.99299201) q[2];
sx q[2];
rz(-3.1057242) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4952132) q[1];
sx q[1];
rz(-2.7732012) q[1];
sx q[1];
rz(0.45763335) q[1];
rz(2.3894541) q[3];
sx q[3];
rz(-2.6400551) q[3];
sx q[3];
rz(-0.88830321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.53434831) q[2];
sx q[2];
rz(-1.4899985) q[2];
sx q[2];
rz(-1.472817) q[2];
rz(-1.2502753) q[3];
sx q[3];
rz(-1.6142802) q[3];
sx q[3];
rz(-0.97243029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.2431353) q[0];
sx q[0];
rz(-0.055483015) q[0];
sx q[0];
rz(-2.6836416) q[0];
rz(-0.10174879) q[1];
sx q[1];
rz(-1.190217) q[1];
sx q[1];
rz(2.8153458) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9407584) q[0];
sx q[0];
rz(-2.5837008) q[0];
sx q[0];
rz(0.87403654) q[0];
rz(2.5102551) q[2];
sx q[2];
rz(-2.7777618) q[2];
sx q[2];
rz(-1.199786) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.40405289) q[1];
sx q[1];
rz(-2.0758325) q[1];
sx q[1];
rz(-2.2801599) q[1];
x q[2];
rz(2.7986854) q[3];
sx q[3];
rz(-0.77292934) q[3];
sx q[3];
rz(-0.29706222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.0094256224) q[2];
sx q[2];
rz(-0.99401179) q[2];
sx q[2];
rz(-1.8780635) q[2];
rz(-3.1215014) q[3];
sx q[3];
rz(-0.76954904) q[3];
sx q[3];
rz(1.9115537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1175213) q[0];
sx q[0];
rz(-0.063952359) q[0];
sx q[0];
rz(-2.5653895) q[0];
rz(-1.4441215) q[1];
sx q[1];
rz(-1.2497808) q[1];
sx q[1];
rz(1.261927) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0281953) q[0];
sx q[0];
rz(-1.602766) q[0];
sx q[0];
rz(1.8799923) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6143965) q[2];
sx q[2];
rz(-1.9771075) q[2];
sx q[2];
rz(0.57002178) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0351959) q[1];
sx q[1];
rz(-0.96265618) q[1];
sx q[1];
rz(1.2713096) q[1];
x q[2];
rz(-0.71407302) q[3];
sx q[3];
rz(-1.1311817) q[3];
sx q[3];
rz(3.1187268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.55804092) q[2];
sx q[2];
rz(-2.5463107) q[2];
sx q[2];
rz(-2.9884647) q[2];
rz(0.88614744) q[3];
sx q[3];
rz(-1.359442) q[3];
sx q[3];
rz(-1.9396293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8837638) q[0];
sx q[0];
rz(-1.9572636) q[0];
sx q[0];
rz(0.15876874) q[0];
rz(1.0348882) q[1];
sx q[1];
rz(-0.95304573) q[1];
sx q[1];
rz(0.75611702) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52126982) q[0];
sx q[0];
rz(-1.6286932) q[0];
sx q[0];
rz(-1.5204563) q[0];
x q[1];
rz(-1.6768084) q[2];
sx q[2];
rz(-0.25463018) q[2];
sx q[2];
rz(2.1531596) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.084001) q[1];
sx q[1];
rz(-2.822031) q[1];
sx q[1];
rz(-2.4366385) q[1];
rz(-pi) q[2];
rz(2.4231195) q[3];
sx q[3];
rz(-0.30820307) q[3];
sx q[3];
rz(-2.001345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5990344) q[2];
sx q[2];
rz(-0.90347806) q[2];
sx q[2];
rz(2.2912045) q[2];
rz(1.6709857) q[3];
sx q[3];
rz(-1.8149523) q[3];
sx q[3];
rz(-2.3175122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65664148) q[0];
sx q[0];
rz(-1.4018207) q[0];
sx q[0];
rz(2.9630419) q[0];
rz(-1.2483596) q[1];
sx q[1];
rz(-1.6105885) q[1];
sx q[1];
rz(0.53370968) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7631012) q[0];
sx q[0];
rz(-2.7437177) q[0];
sx q[0];
rz(0.48872013) q[0];
rz(-2.1458964) q[2];
sx q[2];
rz(-0.786869) q[2];
sx q[2];
rz(-0.58955075) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.935993) q[1];
sx q[1];
rz(-2.471711) q[1];
sx q[1];
rz(1.655446) q[1];
rz(-1.7469095) q[3];
sx q[3];
rz(-0.51100547) q[3];
sx q[3];
rz(-2.3677625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8278213) q[2];
sx q[2];
rz(-1.2510108) q[2];
sx q[2];
rz(-1.0531462) q[2];
rz(0.18038067) q[3];
sx q[3];
rz(-0.80612055) q[3];
sx q[3];
rz(2.7789796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.5617705) q[0];
sx q[0];
rz(-2.1403911) q[0];
sx q[0];
rz(1.2286105) q[0];
rz(-2.6749532) q[1];
sx q[1];
rz(-1.2184315) q[1];
sx q[1];
rz(0.12900464) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59959847) q[0];
sx q[0];
rz(-2.9739673) q[0];
sx q[0];
rz(0.92131281) q[0];
x q[1];
rz(-2.3006347) q[2];
sx q[2];
rz(-1.691868) q[2];
sx q[2];
rz(0.88545495) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1946757) q[1];
sx q[1];
rz(-1.8652152) q[1];
sx q[1];
rz(-0.63398449) q[1];
rz(3.0016162) q[3];
sx q[3];
rz(-1.8365321) q[3];
sx q[3];
rz(1.1418726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1351607) q[2];
sx q[2];
rz(-0.60144037) q[2];
sx q[2];
rz(-0.50430164) q[2];
rz(-2.048061) q[3];
sx q[3];
rz(-1.3239599) q[3];
sx q[3];
rz(-0.95660153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6495551) q[0];
sx q[0];
rz(-1.2832337) q[0];
sx q[0];
rz(1.5851703) q[0];
rz(0.57836142) q[1];
sx q[1];
rz(-1.882694) q[1];
sx q[1];
rz(-0.10231054) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1083173) q[0];
sx q[0];
rz(-2.1731097) q[0];
sx q[0];
rz(0.48745103) q[0];
rz(-pi) q[1];
rz(1.8004216) q[2];
sx q[2];
rz(-1.183325) q[2];
sx q[2];
rz(-1.4619399) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.31033406) q[1];
sx q[1];
rz(-1.4142904) q[1];
sx q[1];
rz(0.308728) q[1];
rz(-pi) q[2];
rz(0.7500521) q[3];
sx q[3];
rz(-2.1382209) q[3];
sx q[3];
rz(-0.22818434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9461225) q[2];
sx q[2];
rz(-1.2637694) q[2];
sx q[2];
rz(0.61402399) q[2];
rz(-2.2139003) q[3];
sx q[3];
rz(-1.5427019) q[3];
sx q[3];
rz(-2.4641002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5880661) q[0];
sx q[0];
rz(-2.8480242) q[0];
sx q[0];
rz(1.199383) q[0];
rz(1.1792432) q[1];
sx q[1];
rz(-1.0284547) q[1];
sx q[1];
rz(-0.95338043) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74443734) q[0];
sx q[0];
rz(-1.0375451) q[0];
sx q[0];
rz(0.6102194) q[0];
rz(-pi) q[1];
rz(-2.1087476) q[2];
sx q[2];
rz(-1.6643833) q[2];
sx q[2];
rz(-2.6757698) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9645365) q[1];
sx q[1];
rz(-2.3971618) q[1];
sx q[1];
rz(1.4076648) q[1];
rz(-pi) q[2];
rz(-2.9913919) q[3];
sx q[3];
rz(-0.48637342) q[3];
sx q[3];
rz(1.9155765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.392268) q[2];
sx q[2];
rz(-0.12568036) q[2];
sx q[2];
rz(-0.18088642) q[2];
rz(-1.1130029) q[3];
sx q[3];
rz(-1.0180611) q[3];
sx q[3];
rz(0.38016144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0522633) q[0];
sx q[0];
rz(-0.90173975) q[0];
sx q[0];
rz(-2.3774636) q[0];
rz(-1.0230505) q[1];
sx q[1];
rz(-1.5307531) q[1];
sx q[1];
rz(-1.6463564) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4039072) q[0];
sx q[0];
rz(-1.4171184) q[0];
sx q[0];
rz(-0.30074461) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33319117) q[2];
sx q[2];
rz(-0.99634561) q[2];
sx q[2];
rz(1.7111749) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.047840441) q[1];
sx q[1];
rz(-1.6421659) q[1];
sx q[1];
rz(2.755295) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9633308) q[3];
sx q[3];
rz(-1.9826673) q[3];
sx q[3];
rz(0.28366006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9493147) q[2];
sx q[2];
rz(-1.01769) q[2];
sx q[2];
rz(-0.44753543) q[2];
rz(0.13628515) q[3];
sx q[3];
rz(-2.648573) q[3];
sx q[3];
rz(-0.59806699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78551453) q[0];
sx q[0];
rz(-2.1165753) q[0];
sx q[0];
rz(0.4726952) q[0];
rz(-0.39974943) q[1];
sx q[1];
rz(-0.70416299) q[1];
sx q[1];
rz(-0.8185111) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21441701) q[0];
sx q[0];
rz(-2.7821988) q[0];
sx q[0];
rz(-2.6136616) q[0];
rz(-pi) q[1];
rz(2.632896) q[2];
sx q[2];
rz(-2.3983068) q[2];
sx q[2];
rz(-2.7270779) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5732729) q[1];
sx q[1];
rz(-1.5740054) q[1];
sx q[1];
rz(-1.8934956) q[1];
x q[2];
rz(0.34719877) q[3];
sx q[3];
rz(-1.7337017) q[3];
sx q[3];
rz(-2.3293975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5588536) q[2];
sx q[2];
rz(-2.0595198) q[2];
sx q[2];
rz(1.2782798) q[2];
rz(-1.9418955) q[3];
sx q[3];
rz(-1.4408828) q[3];
sx q[3];
rz(2.8619158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3661135) q[0];
sx q[0];
rz(-1.4237325) q[0];
sx q[0];
rz(1.7970418) q[0];
rz(1.7300425) q[1];
sx q[1];
rz(-0.81137864) q[1];
sx q[1];
rz(2.484533) q[1];
rz(-0.1409762) q[2];
sx q[2];
rz(-1.5715512) q[2];
sx q[2];
rz(-2.0748009) q[2];
rz(0.90669244) q[3];
sx q[3];
rz(-1.2994874) q[3];
sx q[3];
rz(-0.52728925) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
