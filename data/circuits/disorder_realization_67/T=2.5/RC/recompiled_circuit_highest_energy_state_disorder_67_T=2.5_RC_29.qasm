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
rz(1.2713852) q[0];
sx q[0];
rz(3.1280018) q[0];
sx q[0];
rz(9.4511436) q[0];
rz(0.68323505) q[1];
sx q[1];
rz(-1.3807715) q[1];
sx q[1];
rz(0.071579054) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3535093) q[0];
sx q[0];
rz(-0.026876315) q[0];
sx q[0];
rz(-2.7719017) q[0];
x q[1];
rz(-2.6144892) q[2];
sx q[2];
rz(-0.44473916) q[2];
sx q[2];
rz(1.9183733) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2944379) q[1];
sx q[1];
rz(-1.5757474) q[1];
sx q[1];
rz(1.5520384) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.31622131) q[3];
sx q[3];
rz(-1.8762778) q[3];
sx q[3];
rz(-2.9169634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8611531) q[2];
sx q[2];
rz(-0.70715487) q[2];
sx q[2];
rz(2.6748778) q[2];
rz(0.47433445) q[3];
sx q[3];
rz(-3.1204087) q[3];
sx q[3];
rz(3.1202313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83653432) q[0];
sx q[0];
rz(-2.6439522) q[0];
sx q[0];
rz(-3.1217788) q[0];
rz(1.548798) q[1];
sx q[1];
rz(-2.9217547) q[1];
sx q[1];
rz(-1.6682495) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90813808) q[0];
sx q[0];
rz(-2.0950982) q[0];
sx q[0];
rz(2.3154696) q[0];
x q[1];
rz(-0.85313932) q[2];
sx q[2];
rz(-0.38333508) q[2];
sx q[2];
rz(-0.16259532) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.32343601) q[1];
sx q[1];
rz(-1.4940573) q[1];
sx q[1];
rz(-1.7652926) q[1];
x q[2];
rz(2.3109644) q[3];
sx q[3];
rz(-0.98139709) q[3];
sx q[3];
rz(0.43260655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8389429) q[2];
sx q[2];
rz(-0.5984211) q[2];
sx q[2];
rz(-1.2897276) q[2];
rz(-1.2368115) q[3];
sx q[3];
rz(-0.33331063) q[3];
sx q[3];
rz(0.70241565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1692093) q[0];
sx q[0];
rz(-1.1595668) q[0];
sx q[0];
rz(1.5709391) q[0];
rz(-1.4717357) q[1];
sx q[1];
rz(-1.5262628) q[1];
sx q[1];
rz(0.45438802) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7189302) q[0];
sx q[0];
rz(-2.5276428) q[0];
sx q[0];
rz(-1.802868) q[0];
rz(-pi) q[1];
rz(1.7538007) q[2];
sx q[2];
rz(-0.14438528) q[2];
sx q[2];
rz(2.3600277) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2889298) q[1];
sx q[1];
rz(-1.689012) q[1];
sx q[1];
rz(0.022799076) q[1];
rz(-pi) q[2];
rz(-0.79409365) q[3];
sx q[3];
rz(-2.3183868) q[3];
sx q[3];
rz(-0.87827194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4001974) q[2];
sx q[2];
rz(-0.016308451) q[2];
sx q[2];
rz(-2.7941217) q[2];
rz(2.6853284) q[3];
sx q[3];
rz(-3.1268692) q[3];
sx q[3];
rz(-0.99304503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4149813) q[0];
sx q[0];
rz(-1.9009637) q[0];
sx q[0];
rz(1.4062784) q[0];
rz(-2.6982488) q[1];
sx q[1];
rz(-2.1163546) q[1];
sx q[1];
rz(-1.5700856) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8271874) q[0];
sx q[0];
rz(-1.6256204) q[0];
sx q[0];
rz(-0.64029764) q[0];
x q[1];
rz(1.6353959) q[2];
sx q[2];
rz(-1.6396785) q[2];
sx q[2];
rz(1.2749102) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0519331) q[1];
sx q[1];
rz(-1.5610236) q[1];
sx q[1];
rz(-1.2925413) q[1];
rz(-pi) q[2];
rz(-1.0138233) q[3];
sx q[3];
rz(-0.49684769) q[3];
sx q[3];
rz(-2.5574977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7683679) q[2];
sx q[2];
rz(-2.7478605) q[2];
sx q[2];
rz(3.0623398) q[2];
rz(-1.9521889) q[3];
sx q[3];
rz(-1.7704084) q[3];
sx q[3];
rz(-1.5090548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4382512) q[0];
sx q[0];
rz(-2.6282613) q[0];
sx q[0];
rz(0.80605036) q[0];
rz(2.3015859) q[1];
sx q[1];
rz(-3.1286616) q[1];
sx q[1];
rz(2.3654225) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38473338) q[0];
sx q[0];
rz(-0.98664588) q[0];
sx q[0];
rz(1.7279051) q[0];
rz(-pi) q[1];
rz(-0.010439053) q[2];
sx q[2];
rz(-1.5722646) q[2];
sx q[2];
rz(-1.8515585) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.2649721) q[1];
sx q[1];
rz(-1.435346) q[1];
sx q[1];
rz(-1.5587864) q[1];
x q[2];
rz(-0.40078409) q[3];
sx q[3];
rz(-2.8607607) q[3];
sx q[3];
rz(-3.0147482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.025803056) q[2];
sx q[2];
rz(-1.5610521) q[2];
sx q[2];
rz(2.4157794) q[2];
rz(0.18802655) q[3];
sx q[3];
rz(-3.0811716) q[3];
sx q[3];
rz(2.4316725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1793154) q[0];
sx q[0];
rz(-0.55914068) q[0];
sx q[0];
rz(2.6002) q[0];
rz(-2.9560282) q[1];
sx q[1];
rz(-1.5488397) q[1];
sx q[1];
rz(3.013179) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5207883) q[0];
sx q[0];
rz(-2.6628011) q[0];
sx q[0];
rz(0.6995246) q[0];
rz(-1.6096079) q[2];
sx q[2];
rz(-3.0201742) q[2];
sx q[2];
rz(-1.5281637) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5137607) q[1];
sx q[1];
rz(-2.7134088) q[1];
sx q[1];
rz(-2.9437888) q[1];
rz(1.7592247) q[3];
sx q[3];
rz(-1.5811833) q[3];
sx q[3];
rz(-0.56558181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7554756) q[2];
sx q[2];
rz(-3.0839034) q[2];
sx q[2];
rz(0.84121394) q[2];
rz(-2.9403213) q[3];
sx q[3];
rz(-1.6095251) q[3];
sx q[3];
rz(-2.8819528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5614618) q[0];
sx q[0];
rz(-0.78829563) q[0];
sx q[0];
rz(-1.5808251) q[0];
rz(2.4687817) q[1];
sx q[1];
rz(-1.7377868) q[1];
sx q[1];
rz(3.1127081) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2530841) q[0];
sx q[0];
rz(-0.38034236) q[0];
sx q[0];
rz(2.4230291) q[0];
rz(-pi) q[1];
rz(2.7682226) q[2];
sx q[2];
rz(-2.6334116) q[2];
sx q[2];
rz(2.1477107) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7812271) q[1];
sx q[1];
rz(-2.659428) q[1];
sx q[1];
rz(-2.1025679) q[1];
x q[2];
rz(-2.6050909) q[3];
sx q[3];
rz(-2.287545) q[3];
sx q[3];
rz(-2.0889747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9081356) q[2];
sx q[2];
rz(-0.59671777) q[2];
sx q[2];
rz(2.2133568) q[2];
rz(0.2975896) q[3];
sx q[3];
rz(-2.987515) q[3];
sx q[3];
rz(-0.84460622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.069227844) q[0];
sx q[0];
rz(-0.2437676) q[0];
sx q[0];
rz(0.099076554) q[0];
rz(-2.15436) q[1];
sx q[1];
rz(-1.8376553) q[1];
sx q[1];
rz(-2.6027021) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7659788) q[0];
sx q[0];
rz(-1.5642316) q[0];
sx q[0];
rz(-1.6536941) q[0];
rz(-pi) q[1];
x q[1];
rz(0.024709021) q[2];
sx q[2];
rz(-1.6201783) q[2];
sx q[2];
rz(-1.0339689) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.423279) q[1];
sx q[1];
rz(-0.65650249) q[1];
sx q[1];
rz(-0.61415192) q[1];
rz(-pi) q[2];
rz(-1.5032926) q[3];
sx q[3];
rz(-1.5933523) q[3];
sx q[3];
rz(-1.2872788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.146356) q[2];
sx q[2];
rz(-3.1260999) q[2];
sx q[2];
rz(2.8323979) q[2];
rz(0.63979465) q[3];
sx q[3];
rz(-3.1412509) q[3];
sx q[3];
rz(3.0777847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(0.30534202) q[0];
sx q[0];
rz(-0.59665614) q[0];
sx q[0];
rz(-3.0869361) q[0];
rz(1.1326185) q[1];
sx q[1];
rz(-1.1575969) q[1];
sx q[1];
rz(-1.2861015) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15988201) q[0];
sx q[0];
rz(-1.3927841) q[0];
sx q[0];
rz(-3.0973439) q[0];
rz(-pi) q[1];
rz(-1.6129095) q[2];
sx q[2];
rz(-1.6115808) q[2];
sx q[2];
rz(-1.5342174) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1680345) q[1];
sx q[1];
rz(-1.4806627) q[1];
sx q[1];
rz(-0.22613495) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.70657) q[3];
sx q[3];
rz(-2.6505436) q[3];
sx q[3];
rz(-0.74573475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.22141156) q[2];
sx q[2];
rz(-0.56500089) q[2];
sx q[2];
rz(-2.0378225) q[2];
rz(-1.53299) q[3];
sx q[3];
rz(-0.04403232) q[3];
sx q[3];
rz(-2.485763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.00103818) q[0];
sx q[0];
rz(-0.1796722) q[0];
sx q[0];
rz(3.1371064) q[0];
rz(-1.5525612) q[1];
sx q[1];
rz(-1.6922502) q[1];
sx q[1];
rz(-3.0844614) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4085081) q[0];
sx q[0];
rz(-1.7501236) q[0];
sx q[0];
rz(-3.0441557) q[0];
rz(1.7148027) q[2];
sx q[2];
rz(-1.3504354) q[2];
sx q[2];
rz(-1.8437459) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2559214) q[1];
sx q[1];
rz(-2.2388865) q[1];
sx q[1];
rz(0.1541962) q[1];
x q[2];
rz(2.9895498) q[3];
sx q[3];
rz(-2.3052944) q[3];
sx q[3];
rz(1.1610069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.39157465) q[2];
sx q[2];
rz(-0.027316814) q[2];
sx q[2];
rz(-0.86891437) q[2];
rz(-0.9817552) q[3];
sx q[3];
rz(-3.1120286) q[3];
sx q[3];
rz(0.50299197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-0.045573087) q[0];
sx q[0];
rz(-1.6572784) q[0];
sx q[0];
rz(-1.4838765) q[0];
rz(0.45687301) q[1];
sx q[1];
rz(-2.9869106) q[1];
sx q[1];
rz(3.0966495) q[1];
rz(2.9484684) q[2];
sx q[2];
rz(-2.3681691) q[2];
sx q[2];
rz(-2.9409627) q[2];
rz(2.4784935) q[3];
sx q[3];
rz(-1.6840006) q[3];
sx q[3];
rz(-1.4406289) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
