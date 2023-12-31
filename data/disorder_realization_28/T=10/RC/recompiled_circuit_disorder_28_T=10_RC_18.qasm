OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.93513918) q[0];
sx q[0];
rz(3.9225188) q[0];
sx q[0];
rz(9.6315686) q[0];
rz(-2.7517154) q[1];
sx q[1];
rz(-2.0808527) q[1];
sx q[1];
rz(0.18016711) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6102403) q[0];
sx q[0];
rz(-0.58824476) q[0];
sx q[0];
rz(-2.504185) q[0];
rz(-2.6287659) q[2];
sx q[2];
rz(-1.4564118) q[2];
sx q[2];
rz(1.0637384) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2911421) q[1];
sx q[1];
rz(-0.88911118) q[1];
sx q[1];
rz(-3.046361) q[1];
rz(-pi) q[2];
x q[2];
rz(0.68182919) q[3];
sx q[3];
rz(-1.7141984) q[3];
sx q[3];
rz(1.6553866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.41123286) q[2];
sx q[2];
rz(-0.87324548) q[2];
sx q[2];
rz(1.8475378) q[2];
rz(2.7358352) q[3];
sx q[3];
rz(-1.6399222) q[3];
sx q[3];
rz(-2.7348203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2186573) q[0];
sx q[0];
rz(-2.1056392) q[0];
sx q[0];
rz(-3.0157715) q[0];
rz(-2.3361092) q[1];
sx q[1];
rz(-0.80635726) q[1];
sx q[1];
rz(1.7696101) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1945837) q[0];
sx q[0];
rz(-0.8964552) q[0];
sx q[0];
rz(1.6499004) q[0];
rz(-2.9397474) q[2];
sx q[2];
rz(-2.5864374) q[2];
sx q[2];
rz(1.3943878) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6276715) q[1];
sx q[1];
rz(-1.5253592) q[1];
sx q[1];
rz(-0.25000484) q[1];
rz(-pi) q[2];
rz(-1.3268746) q[3];
sx q[3];
rz(-1.1141277) q[3];
sx q[3];
rz(1.6744542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8877318) q[2];
sx q[2];
rz(-0.68887201) q[2];
sx q[2];
rz(-2.1453693) q[2];
rz(-2.0455202) q[3];
sx q[3];
rz(-1.4527861) q[3];
sx q[3];
rz(0.85038275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6753733) q[0];
sx q[0];
rz(-2.1815364) q[0];
sx q[0];
rz(2.0825785) q[0];
rz(-1.1478708) q[1];
sx q[1];
rz(-0.73906001) q[1];
sx q[1];
rz(-2.0770729) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36944593) q[0];
sx q[0];
rz(-1.9659974) q[0];
sx q[0];
rz(2.4482083) q[0];
rz(-pi) q[1];
rz(-1.9469444) q[2];
sx q[2];
rz(-0.40824879) q[2];
sx q[2];
rz(-0.86366913) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.7980305) q[1];
sx q[1];
rz(-0.31032944) q[1];
sx q[1];
rz(1.5327492) q[1];
x q[2];
rz(-2.9753261) q[3];
sx q[3];
rz(-0.85321745) q[3];
sx q[3];
rz(-2.1239514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.069783) q[2];
sx q[2];
rz(-1.9548543) q[2];
sx q[2];
rz(0.26322571) q[2];
rz(1.1188544) q[3];
sx q[3];
rz(-1.8656732) q[3];
sx q[3];
rz(0.81937218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5082821) q[0];
sx q[0];
rz(-0.044697035) q[0];
sx q[0];
rz(-1.7472965) q[0];
rz(-1.0385723) q[1];
sx q[1];
rz(-1.690052) q[1];
sx q[1];
rz(-1.0669473) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79558668) q[0];
sx q[0];
rz(-0.53115244) q[0];
sx q[0];
rz(-1.6474433) q[0];
rz(0.26950403) q[2];
sx q[2];
rz(-1.4925033) q[2];
sx q[2];
rz(1.6646202) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3921515) q[1];
sx q[1];
rz(-1.9736119) q[1];
sx q[1];
rz(0.72026003) q[1];
x q[2];
rz(-1.0889441) q[3];
sx q[3];
rz(-1.7715766) q[3];
sx q[3];
rz(0.60515412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.84714326) q[2];
sx q[2];
rz(-1.7913982) q[2];
sx q[2];
rz(2.8520544) q[2];
rz(2.6211522) q[3];
sx q[3];
rz(-1.8013022) q[3];
sx q[3];
rz(2.382544) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6510058) q[0];
sx q[0];
rz(-0.0087954272) q[0];
sx q[0];
rz(0.83475137) q[0];
rz(-1.2443776) q[1];
sx q[1];
rz(-1.8908187) q[1];
sx q[1];
rz(-1.7117737) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.017905047) q[0];
sx q[0];
rz(-3.078853) q[0];
sx q[0];
rz(-3.1013261) q[0];
x q[1];
rz(0.28508913) q[2];
sx q[2];
rz(-1.8008917) q[2];
sx q[2];
rz(2.6094112) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1289542) q[1];
sx q[1];
rz(-1.5476989) q[1];
sx q[1];
rz(0.79146339) q[1];
rz(-pi) q[2];
rz(2.2291282) q[3];
sx q[3];
rz(-1.7121127) q[3];
sx q[3];
rz(0.055012881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4804068) q[2];
sx q[2];
rz(-2.1358261) q[2];
sx q[2];
rz(1.8583813) q[2];
rz(-0.13218203) q[3];
sx q[3];
rz(-0.32871267) q[3];
sx q[3];
rz(-2.8018518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4955687) q[0];
sx q[0];
rz(-0.060957242) q[0];
sx q[0];
rz(-2.6640889) q[0];
rz(-1.5006789) q[1];
sx q[1];
rz(-1.532282) q[1];
sx q[1];
rz(0.57055155) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8342499) q[0];
sx q[0];
rz(-1.2120005) q[0];
sx q[0];
rz(-2.8310199) q[0];
rz(-pi) q[1];
rz(-2.4628377) q[2];
sx q[2];
rz(-0.55870134) q[2];
sx q[2];
rz(0.15685454) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.693336) q[1];
sx q[1];
rz(-2.2168471) q[1];
sx q[1];
rz(0.64530428) q[1];
x q[2];
rz(-2.493294) q[3];
sx q[3];
rz(-1.5530506) q[3];
sx q[3];
rz(-1.363021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4346314) q[2];
sx q[2];
rz(-1.0501477) q[2];
sx q[2];
rz(-2.3449507) q[2];
rz(-2.8213275) q[3];
sx q[3];
rz(-2.0551149) q[3];
sx q[3];
rz(-1.2789352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0657848) q[0];
sx q[0];
rz(-1.0362352) q[0];
sx q[0];
rz(-0.35807034) q[0];
rz(-0.2886731) q[1];
sx q[1];
rz(-2.6117548) q[1];
sx q[1];
rz(-1.3105062) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.072648777) q[0];
sx q[0];
rz(-1.3832958) q[0];
sx q[0];
rz(-2.5995273) q[0];
rz(-pi) q[1];
rz(-2.4438379) q[2];
sx q[2];
rz(-2.0332094) q[2];
sx q[2];
rz(0.73730872) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4431745) q[1];
sx q[1];
rz(-2.0768754) q[1];
sx q[1];
rz(1.0531823) q[1];
rz(2.2745423) q[3];
sx q[3];
rz(-0.82074814) q[3];
sx q[3];
rz(1.5501319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.33401176) q[2];
sx q[2];
rz(-2.5568805) q[2];
sx q[2];
rz(2.1441148) q[2];
rz(0.57972646) q[3];
sx q[3];
rz(-0.72967356) q[3];
sx q[3];
rz(-1.4355481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8758133) q[0];
sx q[0];
rz(-1.6515323) q[0];
sx q[0];
rz(-2.8009801) q[0];
rz(-1.9050725) q[1];
sx q[1];
rz(-2.3842936) q[1];
sx q[1];
rz(-1.1901201) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0252359) q[0];
sx q[0];
rz(-1.3394757) q[0];
sx q[0];
rz(-3.0781834) q[0];
rz(-1.4899848) q[2];
sx q[2];
rz(-1.5025856) q[2];
sx q[2];
rz(-2.8516172) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5684143) q[1];
sx q[1];
rz(-1.8548994) q[1];
sx q[1];
rz(1.2969639) q[1];
rz(-0.90772273) q[3];
sx q[3];
rz(-0.66909664) q[3];
sx q[3];
rz(-2.1818386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3999346) q[2];
sx q[2];
rz(-2.2592893) q[2];
sx q[2];
rz(2.7271872) q[2];
rz(-1.7587781) q[3];
sx q[3];
rz(-1.4004935) q[3];
sx q[3];
rz(-0.32489052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8906616) q[0];
sx q[0];
rz(-2.119976) q[0];
sx q[0];
rz(2.9898306) q[0];
rz(1.7157308) q[1];
sx q[1];
rz(-1.7646004) q[1];
sx q[1];
rz(2.192416) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9230726) q[0];
sx q[0];
rz(-2.0109482) q[0];
sx q[0];
rz(-1.4157622) q[0];
x q[1];
rz(1.7911712) q[2];
sx q[2];
rz(-2.1890963) q[2];
sx q[2];
rz(3.0439723) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9526082) q[1];
sx q[1];
rz(-2.4338255) q[1];
sx q[1];
rz(-2.7680552) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6563431) q[3];
sx q[3];
rz(-2.1059548) q[3];
sx q[3];
rz(-2.8458965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7342547) q[2];
sx q[2];
rz(-0.39525017) q[2];
sx q[2];
rz(0.43241832) q[2];
rz(-1.5405103) q[3];
sx q[3];
rz(-1.5099022) q[3];
sx q[3];
rz(0.66175118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0712873) q[0];
sx q[0];
rz(-2.8685331) q[0];
sx q[0];
rz(0.37316698) q[0];
rz(2.7846653) q[1];
sx q[1];
rz(-0.86078763) q[1];
sx q[1];
rz(0.7235136) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7454119) q[0];
sx q[0];
rz(-1.4416845) q[0];
sx q[0];
rz(-2.4095035) q[0];
x q[1];
rz(-2.4621387) q[2];
sx q[2];
rz(-1.6410769) q[2];
sx q[2];
rz(-0.37014222) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3011303) q[1];
sx q[1];
rz(-2.3311619) q[1];
sx q[1];
rz(2.7644972) q[1];
rz(-pi) q[2];
rz(-1.5756597) q[3];
sx q[3];
rz(-1.8412207) q[3];
sx q[3];
rz(-1.566554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2202806) q[2];
sx q[2];
rz(-1.7582515) q[2];
sx q[2];
rz(-0.55345654) q[2];
rz(-2.4297595) q[3];
sx q[3];
rz(-2.7332941) q[3];
sx q[3];
rz(-2.0994983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9217459) q[0];
sx q[0];
rz(-0.60315673) q[0];
sx q[0];
rz(-3.0043816) q[0];
rz(-0.28221054) q[1];
sx q[1];
rz(-2.3574775) q[1];
sx q[1];
rz(-2.2025253) q[1];
rz(2.0201335) q[2];
sx q[2];
rz(-1.1171787) q[2];
sx q[2];
rz(1.6906307) q[2];
rz(2.739493) q[3];
sx q[3];
rz(-1.3105018) q[3];
sx q[3];
rz(-0.49285938) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
