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
rz(-2.6925447) q[0];
sx q[0];
rz(-1.1392051) q[0];
sx q[0];
rz(2.6989302) q[0];
rz(-2.9709385) q[1];
sx q[1];
rz(-0.79167384) q[1];
sx q[1];
rz(-0.72905529) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1478445) q[0];
sx q[0];
rz(-1.1708852) q[0];
sx q[0];
rz(0.73109674) q[0];
x q[1];
rz(-0.38972028) q[2];
sx q[2];
rz(-1.5699631) q[2];
sx q[2];
rz(-0.69609797) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2393136) q[1];
sx q[1];
rz(-2.2444176) q[1];
sx q[1];
rz(2.8301653) q[1];
rz(-pi) q[2];
rz(2.1990875) q[3];
sx q[3];
rz(-1.2344591) q[3];
sx q[3];
rz(-2.9021341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.86058229) q[2];
sx q[2];
rz(-3.0594825) q[2];
sx q[2];
rz(0.78544593) q[2];
rz(-0.9663409) q[3];
sx q[3];
rz(-1.0680501) q[3];
sx q[3];
rz(-2.5016224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.6099089) q[0];
sx q[0];
rz(-0.56782472) q[0];
sx q[0];
rz(-0.63572836) q[0];
rz(2.5081778) q[1];
sx q[1];
rz(-2.3737291) q[1];
sx q[1];
rz(-0.33087081) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71459717) q[0];
sx q[0];
rz(-1.2408371) q[0];
sx q[0];
rz(0.64100463) q[0];
rz(0.026658807) q[2];
sx q[2];
rz(-2.0377198) q[2];
sx q[2];
rz(0.24031642) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4875191) q[1];
sx q[1];
rz(-1.6292058) q[1];
sx q[1];
rz(-3.1078981) q[1];
x q[2];
rz(0.14563609) q[3];
sx q[3];
rz(-0.085509954) q[3];
sx q[3];
rz(-0.95216361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6364381) q[2];
sx q[2];
rz(-2.324489) q[2];
sx q[2];
rz(-2.6804697) q[2];
rz(-0.72502208) q[3];
sx q[3];
rz(-2.6889763) q[3];
sx q[3];
rz(2.4729474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.649699) q[0];
sx q[0];
rz(-1.3839586) q[0];
sx q[0];
rz(0.53832501) q[0];
rz(-3.1321943) q[1];
sx q[1];
rz(-2.5098269) q[1];
sx q[1];
rz(1.1993154) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5814686) q[0];
sx q[0];
rz(-0.5096137) q[0];
sx q[0];
rz(1.4790003) q[0];
rz(-1.7117768) q[2];
sx q[2];
rz(-2.2491697) q[2];
sx q[2];
rz(1.3305261) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6965683) q[1];
sx q[1];
rz(-0.68587947) q[1];
sx q[1];
rz(0.98937757) q[1];
rz(-0.84405293) q[3];
sx q[3];
rz(-1.3793325) q[3];
sx q[3];
rz(-1.8339616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.37375307) q[2];
sx q[2];
rz(-1.491863) q[2];
sx q[2];
rz(2.4191432) q[2];
rz(-2.5433698) q[3];
sx q[3];
rz(-0.0349508) q[3];
sx q[3];
rz(-1.2178864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5825321) q[0];
sx q[0];
rz(-1.2578955) q[0];
sx q[0];
rz(3.1225358) q[0];
rz(1.6825153) q[1];
sx q[1];
rz(-0.2225114) q[1];
sx q[1];
rz(0.51024514) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8718349) q[0];
sx q[0];
rz(-1.8352011) q[0];
sx q[0];
rz(-0.94839903) q[0];
rz(-pi) q[1];
x q[1];
rz(0.34903445) q[2];
sx q[2];
rz(-0.75299298) q[2];
sx q[2];
rz(0.53607625) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8277792) q[1];
sx q[1];
rz(-0.76906068) q[1];
sx q[1];
rz(2.5378347) q[1];
rz(-pi) q[2];
rz(0.41401569) q[3];
sx q[3];
rz(-1.7018205) q[3];
sx q[3];
rz(-1.4690635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9973008) q[2];
sx q[2];
rz(-1.68953) q[2];
sx q[2];
rz(0.22225456) q[2];
rz(-0.079515919) q[3];
sx q[3];
rz(-2.5955718) q[3];
sx q[3];
rz(2.5047746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52614373) q[0];
sx q[0];
rz(-1.0960824) q[0];
sx q[0];
rz(1.8002321) q[0];
rz(0.51097956) q[1];
sx q[1];
rz(-1.7781517) q[1];
sx q[1];
rz(0.43295369) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.110033) q[0];
sx q[0];
rz(-1.3136787) q[0];
sx q[0];
rz(1.5118309) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1673981) q[2];
sx q[2];
rz(-1.3737384) q[2];
sx q[2];
rz(-2.6464406) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9498747) q[1];
sx q[1];
rz(-1.8997802) q[1];
sx q[1];
rz(-2.8803746) q[1];
rz(-2.8540887) q[3];
sx q[3];
rz(-2.3477481) q[3];
sx q[3];
rz(-2.7903729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.13010919) q[2];
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
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3157432) q[0];
sx q[0];
rz(-2.5542927) q[0];
sx q[0];
rz(2.24776) q[0];
rz(2.1122487) q[1];
sx q[1];
rz(-2.4004332) q[1];
sx q[1];
rz(3.0267874) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.288088) q[0];
sx q[0];
rz(-1.813632) q[0];
sx q[0];
rz(0.22260861) q[0];
x q[1];
rz(0.64115142) q[2];
sx q[2];
rz(-1.8617407) q[2];
sx q[2];
rz(3.0464622) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0521212) q[1];
sx q[1];
rz(-1.8473746) q[1];
sx q[1];
rz(0.29015038) q[1];
x q[2];
rz(-1.6345665) q[3];
sx q[3];
rz(-1.8511417) q[3];
sx q[3];
rz(0.75322039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43742314) q[2];
sx q[2];
rz(-2.2544474) q[2];
sx q[2];
rz(2.9278921) q[2];
rz(-2.5050488) q[3];
sx q[3];
rz(-0.53037363) q[3];
sx q[3];
rz(2.4056733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32695025) q[0];
sx q[0];
rz(-0.7433759) q[0];
sx q[0];
rz(-3.0315234) q[0];
rz(3.0649109) q[1];
sx q[1];
rz(-0.86014599) q[1];
sx q[1];
rz(2.0032047) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46318836) q[0];
sx q[0];
rz(-2.7240755) q[0];
sx q[0];
rz(-2.5745086) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5670256) q[2];
sx q[2];
rz(-2.3585883) q[2];
sx q[2];
rz(-1.6996927) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3507211) q[1];
sx q[1];
rz(-0.52164237) q[1];
sx q[1];
rz(1.1952559) q[1];
x q[2];
rz(-1.3567233) q[3];
sx q[3];
rz(-1.8118637) q[3];
sx q[3];
rz(1.9578707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8818605) q[2];
sx q[2];
rz(-2.0596108) q[2];
sx q[2];
rz(0.38121769) q[2];
rz(-3.0242053) q[3];
sx q[3];
rz(-2.6365247) q[3];
sx q[3];
rz(-0.086732619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7748902) q[0];
sx q[0];
rz(-2.369304) q[0];
sx q[0];
rz(0.49852398) q[0];
rz(-0.82930928) q[1];
sx q[1];
rz(-0.81559759) q[1];
sx q[1];
rz(0.097749762) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17247385) q[0];
sx q[0];
rz(-2.0017255) q[0];
sx q[0];
rz(0.33719595) q[0];
rz(-pi) q[1];
rz(0.4011863) q[2];
sx q[2];
rz(-2.7747535) q[2];
sx q[2];
rz(0.40415472) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.5619288) q[1];
sx q[1];
rz(-2.0975998) q[1];
sx q[1];
rz(2.4234523) q[1];
rz(-pi) q[2];
rz(-0.0093383647) q[3];
sx q[3];
rz(-1.2492325) q[3];
sx q[3];
rz(0.044794785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8349614) q[2];
sx q[2];
rz(-1.9913048) q[2];
sx q[2];
rz(-0.59070307) q[2];
rz(3.0513884) q[3];
sx q[3];
rz(-0.43961757) q[3];
sx q[3];
rz(2.2265767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0331994) q[0];
sx q[0];
rz(-2.2985701) q[0];
sx q[0];
rz(-0.75476187) q[0];
rz(-2.8056878) q[1];
sx q[1];
rz(-0.55481189) q[1];
sx q[1];
rz(2.1051443) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2277136) q[0];
sx q[0];
rz(-2.4464985) q[0];
sx q[0];
rz(0.18385017) q[0];
rz(0.039128379) q[2];
sx q[2];
rz(-2.6600398) q[2];
sx q[2];
rz(2.0619543) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7215342) q[1];
sx q[1];
rz(-1.0544525) q[1];
sx q[1];
rz(-1.615585) q[1];
x q[2];
rz(1.3410946) q[3];
sx q[3];
rz(-1.5575081) q[3];
sx q[3];
rz(0.54315996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0466517) q[2];
sx q[2];
rz(-1.4163821) q[2];
sx q[2];
rz(-2.6851192) q[2];
rz(0.60232317) q[3];
sx q[3];
rz(-0.18776247) q[3];
sx q[3];
rz(-0.93835866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8896821) q[0];
sx q[0];
rz(-1.849181) q[0];
sx q[0];
rz(-0.46345261) q[0];
rz(-0.015333029) q[1];
sx q[1];
rz(-3.0749815) q[1];
sx q[1];
rz(2.8212246) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8972337) q[0];
sx q[0];
rz(-2.9787052) q[0];
sx q[0];
rz(-0.49400671) q[0];
x q[1];
rz(0.3139852) q[2];
sx q[2];
rz(-1.3330402) q[2];
sx q[2];
rz(-1.2842922) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.31481701) q[1];
sx q[1];
rz(-0.82524521) q[1];
sx q[1];
rz(-2.25186) q[1];
x q[2];
rz(-0.69971754) q[3];
sx q[3];
rz(-0.43439242) q[3];
sx q[3];
rz(2.7909129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1128803) q[2];
sx q[2];
rz(-2.4148603) q[2];
sx q[2];
rz(-1.0090562) q[2];
rz(2.3446337) q[3];
sx q[3];
rz(-1.2538486) q[3];
sx q[3];
rz(-2.8033281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0491895) q[0];
sx q[0];
rz(-1.5143464) q[0];
sx q[0];
rz(-1.2595246) q[0];
rz(-0.084820329) q[1];
sx q[1];
rz(-1.4823352) q[1];
sx q[1];
rz(-1.7099554) q[1];
rz(-1.3289497) q[2];
sx q[2];
rz(-1.6326661) q[2];
sx q[2];
rz(2.2018955) q[2];
rz(0.30257135) q[3];
sx q[3];
rz(-1.745001) q[3];
sx q[3];
rz(0.59636084) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
