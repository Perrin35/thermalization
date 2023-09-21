OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.59453073) q[0];
sx q[0];
rz(-1.1214331) q[0];
sx q[0];
rz(-2.9601331) q[0];
rz(2.060086) q[1];
sx q[1];
rz(5.6097538) q[1];
sx q[1];
rz(14.654832) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6183137) q[0];
sx q[0];
rz(-1.7821454) q[0];
sx q[0];
rz(-1.1002512) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.62280957) q[2];
sx q[2];
rz(-1.8138759) q[2];
sx q[2];
rz(0.18261766) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.92703544) q[1];
sx q[1];
rz(-1.7662449) q[1];
sx q[1];
rz(-2.3656225) q[1];
rz(-pi) q[2];
rz(2.0041204) q[3];
sx q[3];
rz(-0.37131272) q[3];
sx q[3];
rz(0.565688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9782605) q[2];
sx q[2];
rz(-2.0330727) q[2];
sx q[2];
rz(-1.367761) q[2];
rz(2.1286428) q[3];
sx q[3];
rz(-2.2949341) q[3];
sx q[3];
rz(3.0701385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0435836) q[0];
sx q[0];
rz(-1.4494891) q[0];
sx q[0];
rz(0.59511551) q[0];
rz(-1.0455421) q[1];
sx q[1];
rz(-1.4141934) q[1];
sx q[1];
rz(1.6275303) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1925416) q[0];
sx q[0];
rz(-2.4727614) q[0];
sx q[0];
rz(0.87829725) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3081231) q[2];
sx q[2];
rz(-1.5916628) q[2];
sx q[2];
rz(-1.2534864) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0552057) q[1];
sx q[1];
rz(-1.3761763) q[1];
sx q[1];
rz(-2.9662532) q[1];
rz(-pi) q[2];
rz(-1.4684832) q[3];
sx q[3];
rz(-1.521109) q[3];
sx q[3];
rz(-1.4440086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.65511584) q[2];
sx q[2];
rz(-2.1003508) q[2];
sx q[2];
rz(-2.3454323) q[2];
rz(2.1697309) q[3];
sx q[3];
rz(-2.4367417) q[3];
sx q[3];
rz(0.17175737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3765091) q[0];
sx q[0];
rz(-0.77873814) q[0];
sx q[0];
rz(0.064095108) q[0];
rz(0.31072101) q[1];
sx q[1];
rz(-1.4704082) q[1];
sx q[1];
rz(1.4583189) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38760936) q[0];
sx q[0];
rz(-1.4235272) q[0];
sx q[0];
rz(0.047588451) q[0];
rz(-0.96785737) q[2];
sx q[2];
rz(-1.0087011) q[2];
sx q[2];
rz(-2.2130655) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.872936) q[1];
sx q[1];
rz(-1.5169414) q[1];
sx q[1];
rz(0.40952803) q[1];
rz(-0.052555368) q[3];
sx q[3];
rz(-0.31523809) q[3];
sx q[3];
rz(1.0759575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.98214275) q[2];
sx q[2];
rz(-0.84656707) q[2];
sx q[2];
rz(-0.67908755) q[2];
rz(-1.7689765) q[3];
sx q[3];
rz(-1.8434098) q[3];
sx q[3];
rz(2.1931271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2963592) q[0];
sx q[0];
rz(-0.56931749) q[0];
sx q[0];
rz(1.1244208) q[0];
rz(1.2202948) q[1];
sx q[1];
rz(-1.9492457) q[1];
sx q[1];
rz(1.6569998) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35614466) q[0];
sx q[0];
rz(-1.6011642) q[0];
sx q[0];
rz(-1.5201475) q[0];
rz(-pi) q[1];
rz(-1.053327) q[2];
sx q[2];
rz(-1.4033068) q[2];
sx q[2];
rz(-2.6094764) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.209219) q[1];
sx q[1];
rz(-2.0062431) q[1];
sx q[1];
rz(-2.2030764) q[1];
rz(-pi) q[2];
rz(-2.7709511) q[3];
sx q[3];
rz(-2.3834043) q[3];
sx q[3];
rz(-1.3850118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1241887) q[2];
sx q[2];
rz(-1.8576531) q[2];
sx q[2];
rz(-1.0162639) q[2];
rz(-1.7381564) q[3];
sx q[3];
rz(-1.4973463) q[3];
sx q[3];
rz(1.0884292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6338585) q[0];
sx q[0];
rz(-1.8554747) q[0];
sx q[0];
rz(-2.741709) q[0];
rz(1.1625066) q[1];
sx q[1];
rz(-1.8116654) q[1];
sx q[1];
rz(-2.9679325) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0807053) q[0];
sx q[0];
rz(-1.6433435) q[0];
sx q[0];
rz(1.5990431) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7816254) q[2];
sx q[2];
rz(-1.9569009) q[2];
sx q[2];
rz(-2.2112276) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.225519) q[1];
sx q[1];
rz(-1.5991296) q[1];
sx q[1];
rz(-1.4654935) q[1];
x q[2];
rz(-0.55032702) q[3];
sx q[3];
rz(-0.384207) q[3];
sx q[3];
rz(-1.7526383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6614723) q[2];
sx q[2];
rz(-1.9813333) q[2];
sx q[2];
rz(0.76812569) q[2];
rz(-2.2875732) q[3];
sx q[3];
rz(-1.7212399) q[3];
sx q[3];
rz(0.97222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0356692) q[0];
sx q[0];
rz(-0.24374715) q[0];
sx q[0];
rz(-1.4703898) q[0];
rz(-0.51180965) q[1];
sx q[1];
rz(-2.6302331) q[1];
sx q[1];
rz(-1.1434198) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15018806) q[0];
sx q[0];
rz(-2.9692869) q[0];
sx q[0];
rz(0.50921391) q[0];
rz(-pi) q[1];
rz(-0.60322275) q[2];
sx q[2];
rz(-1.2957186) q[2];
sx q[2];
rz(-0.075866931) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1089576) q[1];
sx q[1];
rz(-2.7320478) q[1];
sx q[1];
rz(-2.0202548) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0295168) q[3];
sx q[3];
rz(-1.8347077) q[3];
sx q[3];
rz(1.2332066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.64289552) q[2];
sx q[2];
rz(-0.84474793) q[2];
sx q[2];
rz(1.5303622) q[2];
rz(1.6879843) q[3];
sx q[3];
rz(-1.0691103) q[3];
sx q[3];
rz(-0.24916515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11944184) q[0];
sx q[0];
rz(-2.3972153) q[0];
sx q[0];
rz(1.7013593) q[0];
rz(-0.72921905) q[1];
sx q[1];
rz(-1.9897285) q[1];
sx q[1];
rz(1.1332606) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92111174) q[0];
sx q[0];
rz(-0.32520121) q[0];
sx q[0];
rz(-0.13334206) q[0];
x q[1];
rz(1.5845756) q[2];
sx q[2];
rz(-2.2829977) q[2];
sx q[2];
rz(2.6028002) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.95722317) q[1];
sx q[1];
rz(-1.2979227) q[1];
sx q[1];
rz(2.2437614) q[1];
x q[2];
rz(-0.53448581) q[3];
sx q[3];
rz(-2.839698) q[3];
sx q[3];
rz(-1.2045977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7363654) q[2];
sx q[2];
rz(-2.0183125) q[2];
sx q[2];
rz(-0.48842946) q[2];
rz(1.3119665) q[3];
sx q[3];
rz(-3.0977111) q[3];
sx q[3];
rz(0.64129889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064780386) q[0];
sx q[0];
rz(-1.8490054) q[0];
sx q[0];
rz(-0.07117614) q[0];
rz(-3.1094303) q[1];
sx q[1];
rz(-1.3379438) q[1];
sx q[1];
rz(-1.9326928) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6541518) q[0];
sx q[0];
rz(-1.7983266) q[0];
sx q[0];
rz(-2.3339416) q[0];
rz(-2.2279943) q[2];
sx q[2];
rz(-0.29619869) q[2];
sx q[2];
rz(1.2072472) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2541131) q[1];
sx q[1];
rz(-1.5548318) q[1];
sx q[1];
rz(1.70114) q[1];
rz(-pi) q[2];
rz(-0.62187059) q[3];
sx q[3];
rz(-1.8727881) q[3];
sx q[3];
rz(2.5627476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.20748392) q[2];
sx q[2];
rz(-2.948163) q[2];
sx q[2];
rz(1.5709546) q[2];
rz(0.87336826) q[3];
sx q[3];
rz(-1.4040754) q[3];
sx q[3];
rz(0.35203716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1677925) q[0];
sx q[0];
rz(-1.5252824) q[0];
sx q[0];
rz(-2.8299676) q[0];
rz(0.82178003) q[1];
sx q[1];
rz(-0.59097925) q[1];
sx q[1];
rz(1.6315546) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25699297) q[0];
sx q[0];
rz(-1.8879226) q[0];
sx q[0];
rz(-1.3604128) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2276149) q[2];
sx q[2];
rz(-1.3859663) q[2];
sx q[2];
rz(-2.1798346) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.77126399) q[1];
sx q[1];
rz(-0.9653829) q[1];
sx q[1];
rz(1.8172812) q[1];
rz(-0.23267965) q[3];
sx q[3];
rz(-1.7494697) q[3];
sx q[3];
rz(-1.0089547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8682378) q[2];
sx q[2];
rz(-0.53871012) q[2];
sx q[2];
rz(-0.81595016) q[2];
rz(0.50968918) q[3];
sx q[3];
rz(-0.78521252) q[3];
sx q[3];
rz(-1.2379439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8907392) q[0];
sx q[0];
rz(-1.7128523) q[0];
sx q[0];
rz(-2.9113286) q[0];
rz(-2.5157805) q[1];
sx q[1];
rz(-0.9451378) q[1];
sx q[1];
rz(-2.4831916) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2021159) q[0];
sx q[0];
rz(-2.3411223) q[0];
sx q[0];
rz(-2.4329484) q[0];
rz(1.2344822) q[2];
sx q[2];
rz(-2.2173777) q[2];
sx q[2];
rz(-0.99415776) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.876993) q[1];
sx q[1];
rz(-2.0473285) q[1];
sx q[1];
rz(2.9546253) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.55024054) q[3];
sx q[3];
rz(-2.5329258) q[3];
sx q[3];
rz(-0.88702162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4663503) q[2];
sx q[2];
rz(-0.8224951) q[2];
sx q[2];
rz(0.33774439) q[2];
rz(2.1181469) q[3];
sx q[3];
rz(-1.7815536) q[3];
sx q[3];
rz(-1.2158998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(0.52453775) q[0];
sx q[0];
rz(-0.90159566) q[0];
sx q[0];
rz(3.0774075) q[0];
rz(1.9032003) q[1];
sx q[1];
rz(-1.4307784) q[1];
sx q[1];
rz(1.4684114) q[1];
rz(-2.9624883) q[2];
sx q[2];
rz(-0.95218198) q[2];
sx q[2];
rz(1.7330963) q[2];
rz(-0.099048793) q[3];
sx q[3];
rz(-2.5669813) q[3];
sx q[3];
rz(3.0909227) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
