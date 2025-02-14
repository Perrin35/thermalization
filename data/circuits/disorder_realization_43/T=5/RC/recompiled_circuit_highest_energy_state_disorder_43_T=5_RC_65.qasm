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
rz(2.9146258) q[0];
sx q[0];
rz(-0.6114971) q[0];
sx q[0];
rz(0.55846941) q[0];
rz(0.59977579) q[1];
sx q[1];
rz(-1.7668626) q[1];
sx q[1];
rz(-3.0629646) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0723018) q[0];
sx q[0];
rz(-2.3218394) q[0];
sx q[0];
rz(2.3425774) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15462523) q[2];
sx q[2];
rz(-2.4569656) q[2];
sx q[2];
rz(-1.3473617) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.78863555) q[1];
sx q[1];
rz(-2.87559) q[1];
sx q[1];
rz(0.31600829) q[1];
rz(0.73689245) q[3];
sx q[3];
rz(-2.6813563) q[3];
sx q[3];
rz(-0.62084197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.89995304) q[2];
sx q[2];
rz(-0.048154801) q[2];
sx q[2];
rz(-2.107035) q[2];
rz(-3.0311846) q[3];
sx q[3];
rz(-1.4805099) q[3];
sx q[3];
rz(0.30606562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0050874) q[0];
sx q[0];
rz(-1.445048) q[0];
sx q[0];
rz(-1.1862296) q[0];
rz(-1.2677445) q[1];
sx q[1];
rz(-0.45919752) q[1];
sx q[1];
rz(1.7523821) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1339552) q[0];
sx q[0];
rz(-1.8764155) q[0];
sx q[0];
rz(1.0108749) q[0];
x q[1];
rz(-2.9302338) q[2];
sx q[2];
rz(-0.49079259) q[2];
sx q[2];
rz(-2.9826814) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4965044) q[1];
sx q[1];
rz(-1.2821272) q[1];
sx q[1];
rz(-1.225807) q[1];
rz(-pi) q[2];
rz(0.94654437) q[3];
sx q[3];
rz(-1.4930176) q[3];
sx q[3];
rz(2.3619224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9691951) q[2];
sx q[2];
rz(-1.1149422) q[2];
sx q[2];
rz(-1.7051075) q[2];
rz(-1.2596333) q[3];
sx q[3];
rz(-0.45537046) q[3];
sx q[3];
rz(-0.026738515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36419511) q[0];
sx q[0];
rz(-1.3854249) q[0];
sx q[0];
rz(-1.2170323) q[0];
rz(0.2150391) q[1];
sx q[1];
rz(-1.8937078) q[1];
sx q[1];
rz(-1.7020114) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.091451784) q[0];
sx q[0];
rz(-0.9147075) q[0];
sx q[0];
rz(1.4138282) q[0];
x q[1];
rz(-1.0595752) q[2];
sx q[2];
rz(-1.6523835) q[2];
sx q[2];
rz(-2.4961584) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2405458) q[1];
sx q[1];
rz(-2.8786297) q[1];
sx q[1];
rz(-2.0790624) q[1];
rz(-pi) q[2];
rz(1.8441299) q[3];
sx q[3];
rz(-1.36051) q[3];
sx q[3];
rz(3.1137636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3415459) q[2];
sx q[2];
rz(-1.505082) q[2];
sx q[2];
rz(2.3254507) q[2];
rz(-0.0084361313) q[3];
sx q[3];
rz(-2.4237027) q[3];
sx q[3];
rz(-1.5661904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7208045) q[0];
sx q[0];
rz(-1.1710465) q[0];
sx q[0];
rz(-2.8724443) q[0];
rz(-1.7587657) q[1];
sx q[1];
rz(-2.5803284) q[1];
sx q[1];
rz(2.6699578) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7155701) q[0];
sx q[0];
rz(-1.4310808) q[0];
sx q[0];
rz(0.84265291) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.083375562) q[2];
sx q[2];
rz(-2.7365587) q[2];
sx q[2];
rz(2.0138559) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.91995222) q[1];
sx q[1];
rz(-0.62949179) q[1];
sx q[1];
rz(2.5821177) q[1];
rz(-pi) q[2];
x q[2];
rz(0.90507953) q[3];
sx q[3];
rz(-1.8536012) q[3];
sx q[3];
rz(1.0156516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.65583324) q[2];
sx q[2];
rz(-1.7344319) q[2];
sx q[2];
rz(1.887623) q[2];
rz(-2.8433825) q[3];
sx q[3];
rz(-1.2757653) q[3];
sx q[3];
rz(1.5378753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1042079) q[0];
sx q[0];
rz(-2.2315114) q[0];
sx q[0];
rz(-0.75697672) q[0];
rz(-1.0589927) q[1];
sx q[1];
rz(-1.4451566) q[1];
sx q[1];
rz(-0.34245488) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9152074) q[0];
sx q[0];
rz(-1.7614375) q[0];
sx q[0];
rz(1.1676428) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2724897) q[2];
sx q[2];
rz(-1.3192186) q[2];
sx q[2];
rz(0.053843018) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7360644) q[1];
sx q[1];
rz(-1.9384512) q[1];
sx q[1];
rz(-0.50761948) q[1];
rz(-pi) q[2];
rz(-2.5352372) q[3];
sx q[3];
rz(-2.6487051) q[3];
sx q[3];
rz(-0.78745251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4464438) q[2];
sx q[2];
rz(-1.6764287) q[2];
sx q[2];
rz(-2.5214419) q[2];
rz(0.54689637) q[3];
sx q[3];
rz(-0.20709012) q[3];
sx q[3];
rz(2.0641522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3243489) q[0];
sx q[0];
rz(-0.20340782) q[0];
sx q[0];
rz(1.9203583) q[0];
rz(-1.1034032) q[1];
sx q[1];
rz(-1.1627448) q[1];
sx q[1];
rz(0.81072909) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2317144) q[0];
sx q[0];
rz(-2.8389769) q[0];
sx q[0];
rz(0.55993863) q[0];
rz(-pi) q[1];
rz(-0.9378433) q[2];
sx q[2];
rz(-1.019291) q[2];
sx q[2];
rz(2.1126975) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.4524881) q[1];
sx q[1];
rz(-2.9924722) q[1];
sx q[1];
rz(1.945687) q[1];
rz(-pi) q[2];
x q[2];
rz(0.91355904) q[3];
sx q[3];
rz(-2.676042) q[3];
sx q[3];
rz(0.54009932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.80060426) q[2];
sx q[2];
rz(-1.1977414) q[2];
sx q[2];
rz(1.1654589) q[2];
rz(-0.10945877) q[3];
sx q[3];
rz(-1.0215003) q[3];
sx q[3];
rz(-0.060733184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3133746) q[0];
sx q[0];
rz(-0.010638588) q[0];
sx q[0];
rz(-2.4995372) q[0];
rz(2.8984046) q[1];
sx q[1];
rz(-2.7311192) q[1];
sx q[1];
rz(2.5004255) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.015804) q[0];
sx q[0];
rz(-2.1282853) q[0];
sx q[0];
rz(2.1899738) q[0];
rz(0.38491048) q[2];
sx q[2];
rz(-2.353047) q[2];
sx q[2];
rz(-0.53764082) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7542759) q[1];
sx q[1];
rz(-0.85016221) q[1];
sx q[1];
rz(2.6971223) q[1];
rz(-1.0585467) q[3];
sx q[3];
rz(-2.6918758) q[3];
sx q[3];
rz(-1.095677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0082561) q[2];
sx q[2];
rz(-1.9828321) q[2];
sx q[2];
rz(0.94375098) q[2];
rz(-2.4933695) q[3];
sx q[3];
rz(-2.4002878) q[3];
sx q[3];
rz(-2.4751723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.7159395) q[0];
sx q[0];
rz(-2.3305927) q[0];
sx q[0];
rz(-0.42914036) q[0];
rz(-0.40402135) q[1];
sx q[1];
rz(-2.7262913) q[1];
sx q[1];
rz(-2.2228352) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12660927) q[0];
sx q[0];
rz(-1.0246687) q[0];
sx q[0];
rz(-1.2864497) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1303004) q[2];
sx q[2];
rz(-2.8545032) q[2];
sx q[2];
rz(-2.8237791) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.67099735) q[1];
sx q[1];
rz(-1.2609224) q[1];
sx q[1];
rz(2.9122374) q[1];
x q[2];
rz(-1.2567344) q[3];
sx q[3];
rz(-1.6410284) q[3];
sx q[3];
rz(2.8633871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1227485) q[2];
sx q[2];
rz(-0.65905535) q[2];
sx q[2];
rz(1.7913691) q[2];
rz(-2.8090779) q[3];
sx q[3];
rz(-2.2612488) q[3];
sx q[3];
rz(1.7727218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3933082) q[0];
sx q[0];
rz(-0.64993334) q[0];
sx q[0];
rz(0.45676029) q[0];
rz(1.7204174) q[1];
sx q[1];
rz(-1.2799809) q[1];
sx q[1];
rz(3.0866887) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53518045) q[0];
sx q[0];
rz(-0.73343432) q[0];
sx q[0];
rz(1.5277521) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0012399) q[2];
sx q[2];
rz(-1.9483742) q[2];
sx q[2];
rz(2.0098639) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.6206606) q[1];
sx q[1];
rz(-1.0865313) q[1];
sx q[1];
rz(-2.8961532) q[1];
rz(2.6087481) q[3];
sx q[3];
rz(-0.60361629) q[3];
sx q[3];
rz(-0.6288022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.88825893) q[2];
sx q[2];
rz(-2.4312225) q[2];
sx q[2];
rz(-0.19676512) q[2];
rz(-1.0737859) q[3];
sx q[3];
rz(-1.1353227) q[3];
sx q[3];
rz(0.73608583) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1573023) q[0];
sx q[0];
rz(-0.8534011) q[0];
sx q[0];
rz(-2.2421457) q[0];
rz(0.31117123) q[1];
sx q[1];
rz(-1.2093465) q[1];
sx q[1];
rz(-1.7412294) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3547826) q[0];
sx q[0];
rz(-1.7799095) q[0];
sx q[0];
rz(1.0037006) q[0];
x q[1];
rz(0.24263361) q[2];
sx q[2];
rz(-1.0972692) q[2];
sx q[2];
rz(0.99019105) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0687814) q[1];
sx q[1];
rz(-2.1635219) q[1];
sx q[1];
rz(-1.0981506) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1651044) q[3];
sx q[3];
rz(-2.1955238) q[3];
sx q[3];
rz(2.4710666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0805936) q[2];
sx q[2];
rz(-0.62005764) q[2];
sx q[2];
rz(1.1617804) q[2];
rz(2.0628085) q[3];
sx q[3];
rz(-0.99083841) q[3];
sx q[3];
rz(-2.259528) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5921191) q[0];
sx q[0];
rz(-1.619512) q[0];
sx q[0];
rz(-1.6721538) q[0];
rz(3.0081765) q[1];
sx q[1];
rz(-1.3795556) q[1];
sx q[1];
rz(2.1605927) q[1];
rz(0.0010112671) q[2];
sx q[2];
rz(-1.5173923) q[2];
sx q[2];
rz(-1.5569729) q[2];
rz(1.4992989) q[3];
sx q[3];
rz(-1.0701978) q[3];
sx q[3];
rz(1.3842907) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
