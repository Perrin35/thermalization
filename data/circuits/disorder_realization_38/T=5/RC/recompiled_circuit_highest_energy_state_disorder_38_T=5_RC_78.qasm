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
rz(0.44904798) q[0];
sx q[0];
rz(-2.0023876) q[0];
sx q[0];
rz(0.44266242) q[0];
rz(-2.9709385) q[1];
sx q[1];
rz(-0.79167384) q[1];
sx q[1];
rz(2.4125374) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9937482) q[0];
sx q[0];
rz(-1.9707075) q[0];
sx q[0];
rz(-2.4104959) q[0];
rz(-pi) q[1];
x q[1];
rz(0.38972028) q[2];
sx q[2];
rz(-1.5699631) q[2];
sx q[2];
rz(-2.4454947) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2393136) q[1];
sx q[1];
rz(-2.2444176) q[1];
sx q[1];
rz(0.31142733) q[1];
x q[2];
rz(-2.1074082) q[3];
sx q[3];
rz(-0.70176673) q[3];
sx q[3];
rz(0.90493363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2810104) q[2];
sx q[2];
rz(-3.0594825) q[2];
sx q[2];
rz(2.3561467) q[2];
rz(0.9663409) q[3];
sx q[3];
rz(-1.0680501) q[3];
sx q[3];
rz(2.5016224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6099089) q[0];
sx q[0];
rz(-0.56782472) q[0];
sx q[0];
rz(0.63572836) q[0];
rz(-0.6334148) q[1];
sx q[1];
rz(-0.76786357) q[1];
sx q[1];
rz(-2.8107218) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4269955) q[0];
sx q[0];
rz(-1.2408371) q[0];
sx q[0];
rz(0.64100463) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0378626) q[2];
sx q[2];
rz(-1.5946009) q[2];
sx q[2];
rz(-1.3184774) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.1304266) q[1];
sx q[1];
rz(-0.067421801) q[1];
sx q[1];
rz(2.0934589) q[1];
rz(-pi) q[2];
rz(-3.0569836) q[3];
sx q[3];
rz(-1.5831909) q[3];
sx q[3];
rz(2.6680714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6364381) q[2];
sx q[2];
rz(-0.81710368) q[2];
sx q[2];
rz(-0.46112296) q[2];
rz(-2.4165706) q[3];
sx q[3];
rz(-0.45261639) q[3];
sx q[3];
rz(-0.66864526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4918936) q[0];
sx q[0];
rz(-1.3839586) q[0];
sx q[0];
rz(2.6032676) q[0];
rz(0.0093983924) q[1];
sx q[1];
rz(-2.5098269) q[1];
sx q[1];
rz(-1.9422772) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.050735) q[0];
sx q[0];
rz(-1.6155302) q[0];
sx q[0];
rz(1.0629774) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4583425) q[2];
sx q[2];
rz(-1.6804192) q[2];
sx q[2];
rz(2.8124968) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6965683) q[1];
sx q[1];
rz(-2.4557132) q[1];
sx q[1];
rz(-2.1522151) q[1];
rz(-0.25377804) q[3];
sx q[3];
rz(-2.2814085) q[3];
sx q[3];
rz(-2.7108278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.37375307) q[2];
sx q[2];
rz(-1.491863) q[2];
sx q[2];
rz(0.72244942) q[2];
rz(-0.59822285) q[3];
sx q[3];
rz(-3.1066419) q[3];
sx q[3];
rz(-1.2178864) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55906051) q[0];
sx q[0];
rz(-1.8836972) q[0];
sx q[0];
rz(-0.019056888) q[0];
rz(1.6825153) q[1];
sx q[1];
rz(-2.9190813) q[1];
sx q[1];
rz(2.6313475) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.048240926) q[0];
sx q[0];
rz(-0.66931319) q[0];
sx q[0];
rz(2.0055683) q[0];
rz(2.7925582) q[2];
sx q[2];
rz(-2.3885997) q[2];
sx q[2];
rz(0.53607625) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6905247) q[1];
sx q[1];
rz(-2.1803586) q[1];
sx q[1];
rz(1.068348) q[1];
rz(2.727577) q[3];
sx q[3];
rz(-1.7018205) q[3];
sx q[3];
rz(1.4690635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.14429188) q[2];
sx q[2];
rz(-1.68953) q[2];
sx q[2];
rz(2.9193381) q[2];
rz(0.079515919) q[3];
sx q[3];
rz(-2.5955718) q[3];
sx q[3];
rz(0.63681805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6154489) q[0];
sx q[0];
rz(-2.0455102) q[0];
sx q[0];
rz(1.8002321) q[0];
rz(-2.6306131) q[1];
sx q[1];
rz(-1.363441) q[1];
sx q[1];
rz(2.708639) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.110033) q[0];
sx q[0];
rz(-1.3136787) q[0];
sx q[0];
rz(-1.5118309) q[0];
x q[1];
rz(1.1673981) q[2];
sx q[2];
rz(-1.7678542) q[2];
sx q[2];
rz(-2.6464406) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4652327) q[1];
sx q[1];
rz(-1.3238878) q[1];
sx q[1];
rz(1.2311155) q[1];
rz(-pi) q[2];
rz(0.28750395) q[3];
sx q[3];
rz(-2.3477481) q[3];
sx q[3];
rz(0.35121976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0114835) q[2];
sx q[2];
rz(-2.2138962) q[2];
sx q[2];
rz(-2.1437272) q[2];
rz(-0.70895553) q[3];
sx q[3];
rz(-0.38781375) q[3];
sx q[3];
rz(0.97878218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.3157432) q[0];
sx q[0];
rz(-2.5542927) q[0];
sx q[0];
rz(2.24776) q[0];
rz(1.029344) q[1];
sx q[1];
rz(-2.4004332) q[1];
sx q[1];
rz(-3.0267874) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.608716) q[0];
sx q[0];
rz(-0.32793697) q[0];
sx q[0];
rz(0.84285835) q[0];
rz(-0.64115142) q[2];
sx q[2];
rz(-1.8617407) q[2];
sx q[2];
rz(0.095130446) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5789168) q[1];
sx q[1];
rz(-1.8496184) q[1];
sx q[1];
rz(-1.2827966) q[1];
rz(-pi) q[2];
rz(-1.5070262) q[3];
sx q[3];
rz(-1.8511417) q[3];
sx q[3];
rz(-0.75322039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.43742314) q[2];
sx q[2];
rz(-0.88714522) q[2];
sx q[2];
rz(2.9278921) q[2];
rz(2.5050488) q[3];
sx q[3];
rz(-0.53037363) q[3];
sx q[3];
rz(-2.4056733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32695025) q[0];
sx q[0];
rz(-2.3982168) q[0];
sx q[0];
rz(0.11006926) q[0];
rz(3.0649109) q[1];
sx q[1];
rz(-0.86014599) q[1];
sx q[1];
rz(2.0032047) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0698893) q[0];
sx q[0];
rz(-1.919863) q[0];
sx q[0];
rz(1.3368692) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5745671) q[2];
sx q[2];
rz(-2.3585883) q[2];
sx q[2];
rz(1.4418999) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.79087154) q[1];
sx q[1];
rz(-2.6199503) q[1];
sx q[1];
rz(1.9463368) q[1];
rz(2.8951169) q[3];
sx q[3];
rz(-1.3630056) q[3];
sx q[3];
rz(2.7026619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.25973213) q[2];
sx q[2];
rz(-2.0596108) q[2];
sx q[2];
rz(-2.760375) q[2];
rz(3.0242053) q[3];
sx q[3];
rz(-2.6365247) q[3];
sx q[3];
rz(0.086732619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7748902) q[0];
sx q[0];
rz(-0.77228868) q[0];
sx q[0];
rz(-0.49852398) q[0];
rz(2.3122834) q[1];
sx q[1];
rz(-2.3259951) q[1];
sx q[1];
rz(3.0438429) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9691188) q[0];
sx q[0];
rz(-1.1398672) q[0];
sx q[0];
rz(-0.33719595) q[0];
x q[1];
rz(1.421861) q[2];
sx q[2];
rz(-1.907299) q[2];
sx q[2];
rz(0.83068427) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.48673442) q[1];
sx q[1];
rz(-0.86198258) q[1];
sx q[1];
rz(-2.4177119) q[1];
x q[2];
rz(-0.0093383647) q[3];
sx q[3];
rz(-1.2492325) q[3];
sx q[3];
rz(-3.0967979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8349614) q[2];
sx q[2];
rz(-1.1502879) q[2];
sx q[2];
rz(2.5508896) q[2];
rz(-0.090204209) q[3];
sx q[3];
rz(-0.43961757) q[3];
sx q[3];
rz(2.2265767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10839323) q[0];
sx q[0];
rz(-2.2985701) q[0];
sx q[0];
rz(2.3868308) q[0];
rz(-0.33590487) q[1];
sx q[1];
rz(-2.5867808) q[1];
sx q[1];
rz(-1.0364484) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1514318) q[0];
sx q[0];
rz(-0.88967547) q[0];
sx q[0];
rz(1.7220884) q[0];
rz(-pi) q[1];
rz(3.1024643) q[2];
sx q[2];
rz(-2.6600398) q[2];
sx q[2];
rz(1.0796384) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4200585) q[1];
sx q[1];
rz(-1.0544525) q[1];
sx q[1];
rz(1.615585) q[1];
rz(0.01364661) q[3];
sx q[3];
rz(-1.8004775) q[3];
sx q[3];
rz(-2.1170634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.094940946) q[2];
sx q[2];
rz(-1.4163821) q[2];
sx q[2];
rz(2.6851192) q[2];
rz(-0.60232317) q[3];
sx q[3];
rz(-2.9538302) q[3];
sx q[3];
rz(2.203234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8896821) q[0];
sx q[0];
rz(-1.2924117) q[0];
sx q[0];
rz(-0.46345261) q[0];
rz(0.015333029) q[1];
sx q[1];
rz(-0.066611193) q[1];
sx q[1];
rz(2.8212246) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8149016) q[0];
sx q[0];
rz(-1.6477655) q[0];
sx q[0];
rz(-0.14369731) q[0];
rz(-pi) q[1];
rz(1.8202844) q[2];
sx q[2];
rz(-1.87566) q[2];
sx q[2];
rz(-2.9314205) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7586472) q[1];
sx q[1];
rz(-1.0898814) q[1];
sx q[1];
rz(2.2702899) q[1];
rz(-pi) q[2];
rz(-2.4418751) q[3];
sx q[3];
rz(-2.7072002) q[3];
sx q[3];
rz(-0.35067973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.028712332) q[2];
sx q[2];
rz(-0.72673231) q[2];
sx q[2];
rz(-1.0090562) q[2];
rz(2.3446337) q[3];
sx q[3];
rz(-1.8877441) q[3];
sx q[3];
rz(-0.33826452) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0924031) q[0];
sx q[0];
rz(-1.5143464) q[0];
sx q[0];
rz(-1.2595246) q[0];
rz(3.0567723) q[1];
sx q[1];
rz(-1.4823352) q[1];
sx q[1];
rz(-1.7099554) q[1];
rz(-3.0778733) q[2];
sx q[2];
rz(-1.3294217) q[2];
sx q[2];
rz(0.64634993) q[2];
rz(-2.6081035) q[3];
sx q[3];
rz(-2.7937952) q[3];
sx q[3];
rz(-1.4813012) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
