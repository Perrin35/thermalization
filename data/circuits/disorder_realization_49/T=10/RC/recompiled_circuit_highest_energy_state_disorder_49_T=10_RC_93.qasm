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
rz(-2.7558514) q[0];
sx q[0];
rz(-2.1585611) q[0];
sx q[0];
rz(-2.5842343) q[0];
rz(-1.9269257) q[1];
sx q[1];
rz(-2.2872556) q[1];
sx q[1];
rz(0.94205034) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5621126) q[0];
sx q[0];
rz(-1.3619553) q[0];
sx q[0];
rz(-0.12128854) q[0];
x q[1];
rz(0.73150191) q[2];
sx q[2];
rz(-2.3580551) q[2];
sx q[2];
rz(0.48043007) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.57191932) q[1];
sx q[1];
rz(-1.633857) q[1];
sx q[1];
rz(2.75534) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1325087) q[3];
sx q[3];
rz(-1.0367437) q[3];
sx q[3];
rz(0.81919794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4152834) q[2];
sx q[2];
rz(-1.7674663) q[2];
sx q[2];
rz(-1.3709925) q[2];
rz(0.81909424) q[3];
sx q[3];
rz(-2.2947125) q[3];
sx q[3];
rz(-0.11025652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46721989) q[0];
sx q[0];
rz(-1.9085566) q[0];
sx q[0];
rz(-2.9357173) q[0];
rz(-1.3174093) q[1];
sx q[1];
rz(-1.2818047) q[1];
sx q[1];
rz(0.46897108) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7531573) q[0];
sx q[0];
rz(-0.65555619) q[0];
sx q[0];
rz(-0.038319328) q[0];
rz(1.1577206) q[2];
sx q[2];
rz(-2.6402355) q[2];
sx q[2];
rz(0.92228466) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6862985) q[1];
sx q[1];
rz(-1.0166234) q[1];
sx q[1];
rz(-0.55880736) q[1];
rz(2.3875565) q[3];
sx q[3];
rz(-2.0181351) q[3];
sx q[3];
rz(-3.0424117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5466902) q[2];
sx q[2];
rz(-2.3089843) q[2];
sx q[2];
rz(-0.62758315) q[2];
rz(-1.0707431) q[3];
sx q[3];
rz(-2.6379733) q[3];
sx q[3];
rz(2.812775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4048432) q[0];
sx q[0];
rz(-0.81031814) q[0];
sx q[0];
rz(-2.3531083) q[0];
rz(2.5704747) q[1];
sx q[1];
rz(-1.8126789) q[1];
sx q[1];
rz(-1.4459389) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2803422) q[0];
sx q[0];
rz(-0.26901562) q[0];
sx q[0];
rz(1.6205257) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6094016) q[2];
sx q[2];
rz(-1.0330832) q[2];
sx q[2];
rz(-0.78372389) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7608628) q[1];
sx q[1];
rz(-0.31425409) q[1];
sx q[1];
rz(-0.32829378) q[1];
rz(2.6258008) q[3];
sx q[3];
rz(-0.63871562) q[3];
sx q[3];
rz(-1.1209008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0995522) q[2];
sx q[2];
rz(-2.6991762) q[2];
sx q[2];
rz(0.37910795) q[2];
rz(-1.3742617) q[3];
sx q[3];
rz(-0.97135025) q[3];
sx q[3];
rz(-0.9837392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9336201) q[0];
sx q[0];
rz(-2.3868028) q[0];
sx q[0];
rz(2.2291613) q[0];
rz(1.3937996) q[1];
sx q[1];
rz(-0.50519609) q[1];
sx q[1];
rz(2.8392653) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5057004) q[0];
sx q[0];
rz(-1.6306584) q[0];
sx q[0];
rz(-3.0144431) q[0];
rz(-pi) q[1];
rz(-2.7458397) q[2];
sx q[2];
rz(-0.89598334) q[2];
sx q[2];
rz(-0.84854919) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.15790882) q[1];
sx q[1];
rz(-1.6140466) q[1];
sx q[1];
rz(-1.3614348) q[1];
rz(-pi) q[2];
rz(0.793806) q[3];
sx q[3];
rz(-0.86177877) q[3];
sx q[3];
rz(-0.10464917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6572774) q[2];
sx q[2];
rz(-2.4399098) q[2];
sx q[2];
rz(2.5122128) q[2];
rz(1.7221919) q[3];
sx q[3];
rz(-1.8604167) q[3];
sx q[3];
rz(-0.70845503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.092954271) q[0];
sx q[0];
rz(-1.0372256) q[0];
sx q[0];
rz(-1.8222437) q[0];
rz(-2.0535779) q[1];
sx q[1];
rz(-1.8698147) q[1];
sx q[1];
rz(1.4647269) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3528093) q[0];
sx q[0];
rz(-1.6029198) q[0];
sx q[0];
rz(-0.87091586) q[0];
rz(-0.044243926) q[2];
sx q[2];
rz(-2.508456) q[2];
sx q[2];
rz(0.69685941) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.32680997) q[1];
sx q[1];
rz(-1.9091354) q[1];
sx q[1];
rz(1.7717351) q[1];
x q[2];
rz(-1.7978908) q[3];
sx q[3];
rz(-0.97714564) q[3];
sx q[3];
rz(-2.32292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.69270837) q[2];
sx q[2];
rz(-3.0227737) q[2];
sx q[2];
rz(2.9046655) q[2];
rz(-2.1610625) q[3];
sx q[3];
rz(-1.7122995) q[3];
sx q[3];
rz(0.94394365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9834845) q[0];
sx q[0];
rz(-0.10949245) q[0];
sx q[0];
rz(0.14516251) q[0];
rz(-1.521184) q[1];
sx q[1];
rz(-2.3416134) q[1];
sx q[1];
rz(-0.0072341166) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63165346) q[0];
sx q[0];
rz(-1.8064033) q[0];
sx q[0];
rz(0.8610691) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5356338) q[2];
sx q[2];
rz(-2.3176386) q[2];
sx q[2];
rz(0.29203019) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0238432) q[1];
sx q[1];
rz(-2.8850318) q[1];
sx q[1];
rz(1.9909977) q[1];
rz(0.16611734) q[3];
sx q[3];
rz(-1.5990077) q[3];
sx q[3];
rz(-1.3835554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2468557) q[2];
sx q[2];
rz(-1.1892908) q[2];
sx q[2];
rz(2.5640633) q[2];
rz(1.9891116) q[3];
sx q[3];
rz(-2.5566176) q[3];
sx q[3];
rz(-1.729689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57303992) q[0];
sx q[0];
rz(-2.9459406) q[0];
sx q[0];
rz(1.0618807) q[0];
rz(-1.7537687) q[1];
sx q[1];
rz(-1.2304708) q[1];
sx q[1];
rz(1.498163) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6681435) q[0];
sx q[0];
rz(-1.7034344) q[0];
sx q[0];
rz(2.7154865) q[0];
x q[1];
rz(1.241356) q[2];
sx q[2];
rz(-0.95166517) q[2];
sx q[2];
rz(-0.090290221) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0727783) q[1];
sx q[1];
rz(-2.1238951) q[1];
sx q[1];
rz(-0.26550737) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2467674) q[3];
sx q[3];
rz(-1.1041485) q[3];
sx q[3];
rz(2.4396536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3987223) q[2];
sx q[2];
rz(-3.0169432) q[2];
sx q[2];
rz(2.2275662) q[2];
rz(-2.3877609) q[3];
sx q[3];
rz(-1.2316615) q[3];
sx q[3];
rz(0.45647538) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1330426) q[0];
sx q[0];
rz(-2.8164016) q[0];
sx q[0];
rz(-3.131026) q[0];
rz(-2.1376624) q[1];
sx q[1];
rz(-1.7251451) q[1];
sx q[1];
rz(-3.1026057) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9462601) q[0];
sx q[0];
rz(-1.8882013) q[0];
sx q[0];
rz(-0.18407777) q[0];
rz(-pi) q[1];
rz(-3.1124074) q[2];
sx q[2];
rz(-2.8246636) q[2];
sx q[2];
rz(-2.2259797) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.4036364) q[1];
sx q[1];
rz(-1.1925329) q[1];
sx q[1];
rz(2.4191816) q[1];
x q[2];
rz(-1.4812864) q[3];
sx q[3];
rz(-2.2319372) q[3];
sx q[3];
rz(1.4941708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5616592) q[2];
sx q[2];
rz(-1.9114405) q[2];
sx q[2];
rz(0.099543355) q[2];
rz(-1.2933412) q[3];
sx q[3];
rz(-1.2852531) q[3];
sx q[3];
rz(-0.37555638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.471591) q[0];
sx q[0];
rz(-1.3674068) q[0];
sx q[0];
rz(-1.8294096) q[0];
rz(-0.52681628) q[1];
sx q[1];
rz(-2.3878038) q[1];
sx q[1];
rz(0.83676338) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0055702607) q[0];
sx q[0];
rz(-1.3107398) q[0];
sx q[0];
rz(-1.2664766) q[0];
x q[1];
rz(-2.2808321) q[2];
sx q[2];
rz(-0.75131159) q[2];
sx q[2];
rz(-1.4850791) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.93406127) q[1];
sx q[1];
rz(-0.69489266) q[1];
sx q[1];
rz(-2.7259105) q[1];
x q[2];
rz(-2.2403735) q[3];
sx q[3];
rz(-2.5680084) q[3];
sx q[3];
rz(0.18665126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7709363) q[2];
sx q[2];
rz(-0.61345658) q[2];
sx q[2];
rz(1.9343617) q[2];
rz(1.3197445) q[3];
sx q[3];
rz(-1.5650257) q[3];
sx q[3];
rz(-2.3453662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.094548263) q[0];
sx q[0];
rz(-2.6667509) q[0];
sx q[0];
rz(2.4249401) q[0];
rz(-1.5776177) q[1];
sx q[1];
rz(-1.8378704) q[1];
sx q[1];
rz(-0.61734739) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3302688) q[0];
sx q[0];
rz(-1.7416546) q[0];
sx q[0];
rz(-2.2854684) q[0];
rz(-pi) q[1];
rz(-0.1495628) q[2];
sx q[2];
rz(-2.0345275) q[2];
sx q[2];
rz(0.92926393) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6177298) q[1];
sx q[1];
rz(-0.81392127) q[1];
sx q[1];
rz(0.62266962) q[1];
rz(-pi) q[2];
rz(2.8598911) q[3];
sx q[3];
rz(-0.94958491) q[3];
sx q[3];
rz(-3.1007183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8120332) q[2];
sx q[2];
rz(-2.7541408) q[2];
sx q[2];
rz(-2.6623902) q[2];
rz(-1.6241578) q[3];
sx q[3];
rz(-1.7370217) q[3];
sx q[3];
rz(-1.6421389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0064405) q[0];
sx q[0];
rz(-1.3545481) q[0];
sx q[0];
rz(-0.56726278) q[0];
rz(2.3334423) q[1];
sx q[1];
rz(-1.6193401) q[1];
sx q[1];
rz(-1.5338939) q[1];
rz(1.7980747) q[2];
sx q[2];
rz(-0.53126104) q[2];
sx q[2];
rz(1.9781611) q[2];
rz(2.3324689) q[3];
sx q[3];
rz(-1.7715142) q[3];
sx q[3];
rz(0.51668744) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
