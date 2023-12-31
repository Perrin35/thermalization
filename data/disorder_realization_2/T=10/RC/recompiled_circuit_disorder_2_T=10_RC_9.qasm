OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6931273) q[0];
sx q[0];
rz(-0.52283302) q[0];
sx q[0];
rz(-0.62358207) q[0];
rz(-2.8514255) q[1];
sx q[1];
rz(-0.71915141) q[1];
sx q[1];
rz(-2.6410988) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5252285) q[0];
sx q[0];
rz(-0.59983569) q[0];
sx q[0];
rz(-3.1347549) q[0];
x q[1];
rz(-1.0153158) q[2];
sx q[2];
rz(-2.9657288) q[2];
sx q[2];
rz(-1.6989087) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9762293) q[1];
sx q[1];
rz(-2.2803218) q[1];
sx q[1];
rz(2.707259) q[1];
x q[2];
rz(1.2455363) q[3];
sx q[3];
rz(-0.45913011) q[3];
sx q[3];
rz(1.9838651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3502675) q[2];
sx q[2];
rz(-1.2056377) q[2];
sx q[2];
rz(-1.2228489) q[2];
rz(-1.6932999) q[3];
sx q[3];
rz(-2.1494614) q[3];
sx q[3];
rz(0.97035113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37110776) q[0];
sx q[0];
rz(-1.6587695) q[0];
sx q[0];
rz(-2.0626542) q[0];
rz(-1.7547912) q[1];
sx q[1];
rz(-2.3290122) q[1];
sx q[1];
rz(-0.66545495) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7121885) q[0];
sx q[0];
rz(-2.5948338) q[0];
sx q[0];
rz(2.5736546) q[0];
x q[1];
rz(2.8198492) q[2];
sx q[2];
rz(-1.4085359) q[2];
sx q[2];
rz(-2.4010047) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7484819) q[1];
sx q[1];
rz(-2.0512274) q[1];
sx q[1];
rz(-1.837681) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1075675) q[3];
sx q[3];
rz(-2.2799387) q[3];
sx q[3];
rz(-0.89108407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.60454303) q[2];
sx q[2];
rz(-1.4031354) q[2];
sx q[2];
rz(0.075142168) q[2];
rz(1.4705307) q[3];
sx q[3];
rz(-2.0139147) q[3];
sx q[3];
rz(2.4501734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5866518) q[0];
sx q[0];
rz(-1.2723158) q[0];
sx q[0];
rz(-2.3828322) q[0];
rz(-1.2930019) q[1];
sx q[1];
rz(-1.2932152) q[1];
sx q[1];
rz(1.4000777) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5132719) q[0];
sx q[0];
rz(-1.8372046) q[0];
sx q[0];
rz(2.7404286) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6252458) q[2];
sx q[2];
rz(-0.22908224) q[2];
sx q[2];
rz(0.89876995) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8850419) q[1];
sx q[1];
rz(-0.98837822) q[1];
sx q[1];
rz(-0.38862733) q[1];
rz(-pi) q[2];
rz(1.8758043) q[3];
sx q[3];
rz(-1.4121571) q[3];
sx q[3];
rz(-2.1624485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.65163461) q[2];
sx q[2];
rz(-2.659446) q[2];
sx q[2];
rz(-0.65762323) q[2];
rz(1.1714606) q[3];
sx q[3];
rz(-1.6572584) q[3];
sx q[3];
rz(-1.7224147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
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
rz(-0.79384971) q[0];
sx q[0];
rz(-2.1934953) q[0];
sx q[0];
rz(-1.6963652) q[0];
rz(1.6943278) q[1];
sx q[1];
rz(-1.6479965) q[1];
sx q[1];
rz(2.7935374) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6570243) q[0];
sx q[0];
rz(-2.1977402) q[0];
sx q[0];
rz(-1.2059962) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1966755) q[2];
sx q[2];
rz(-1.4305563) q[2];
sx q[2];
rz(-1.8082878) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0192249) q[1];
sx q[1];
rz(-0.61673635) q[1];
sx q[1];
rz(1.7255746) q[1];
rz(-pi) q[2];
rz(0.69339852) q[3];
sx q[3];
rz(-2.9326673) q[3];
sx q[3];
rz(2.8914176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.41670123) q[2];
sx q[2];
rz(-1.8323703) q[2];
sx q[2];
rz(-2.7187738) q[2];
rz(-0.73741284) q[3];
sx q[3];
rz(-2.335572) q[3];
sx q[3];
rz(0.084658682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2622862) q[0];
sx q[0];
rz(-1.8286185) q[0];
sx q[0];
rz(2.6111531) q[0];
rz(2.2166705) q[1];
sx q[1];
rz(-1.5735807) q[1];
sx q[1];
rz(1.2984498) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12621524) q[0];
sx q[0];
rz(-1.4091638) q[0];
sx q[0];
rz(2.9956908) q[0];
rz(-pi) q[1];
rz(-0.41579397) q[2];
sx q[2];
rz(-2.4928164) q[2];
sx q[2];
rz(-0.4817889) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7057719) q[1];
sx q[1];
rz(-2.768369) q[1];
sx q[1];
rz(2.865764) q[1];
rz(-2.2523746) q[3];
sx q[3];
rz(-0.71435706) q[3];
sx q[3];
rz(1.0970955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2003145) q[2];
sx q[2];
rz(-1.1555187) q[2];
sx q[2];
rz(0.47362622) q[2];
rz(-0.099362699) q[3];
sx q[3];
rz(-1.2800346) q[3];
sx q[3];
rz(0.84053269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50399238) q[0];
sx q[0];
rz(-1.7853328) q[0];
sx q[0];
rz(2.1437058) q[0];
rz(2.2672794) q[1];
sx q[1];
rz(-1.0214146) q[1];
sx q[1];
rz(0.46674892) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022910718) q[0];
sx q[0];
rz(-2.4372299) q[0];
sx q[0];
rz(-2.8055311) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7919962) q[2];
sx q[2];
rz(-1.9890519) q[2];
sx q[2];
rz(0.43005558) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7044428) q[1];
sx q[1];
rz(-1.721742) q[1];
sx q[1];
rz(1.6187861) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0420226) q[3];
sx q[3];
rz(-0.18897945) q[3];
sx q[3];
rz(2.6721862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.85990396) q[2];
sx q[2];
rz(-0.47912654) q[2];
sx q[2];
rz(-1.5647282) q[2];
rz(-0.62670296) q[3];
sx q[3];
rz(-1.3650711) q[3];
sx q[3];
rz(-0.68157649) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3307813) q[0];
sx q[0];
rz(-2.2450876) q[0];
sx q[0];
rz(2.7217641) q[0];
rz(0.22142521) q[1];
sx q[1];
rz(-0.47859335) q[1];
sx q[1];
rz(-1.5931169) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3395183) q[0];
sx q[0];
rz(-0.087856494) q[0];
sx q[0];
rz(3.0457892) q[0];
x q[1];
rz(0.19303796) q[2];
sx q[2];
rz(-1.8503354) q[2];
sx q[2];
rz(-0.9370196) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6429813) q[1];
sx q[1];
rz(-1.6763408) q[1];
sx q[1];
rz(-1.5555744) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1920131) q[3];
sx q[3];
rz(-0.55286828) q[3];
sx q[3];
rz(-2.051193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.55591136) q[2];
sx q[2];
rz(-2.6101117) q[2];
sx q[2];
rz(-1.4038203) q[2];
rz(-0.79706556) q[3];
sx q[3];
rz(-0.46204391) q[3];
sx q[3];
rz(-1.2169303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7109011) q[0];
sx q[0];
rz(-2.4276908) q[0];
sx q[0];
rz(-0.28924334) q[0];
rz(-0.62492433) q[1];
sx q[1];
rz(-1.1661252) q[1];
sx q[1];
rz(1.3141059) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1159191) q[0];
sx q[0];
rz(-2.325255) q[0];
sx q[0];
rz(1.7220108) q[0];
rz(0.35347519) q[2];
sx q[2];
rz(-2.3620053) q[2];
sx q[2];
rz(-1.0970864) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9574979) q[1];
sx q[1];
rz(-3.0093319) q[1];
sx q[1];
rz(1.1485419) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2460327) q[3];
sx q[3];
rz(-2.357558) q[3];
sx q[3];
rz(1.8079545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4079995) q[2];
sx q[2];
rz(-1.5635798) q[2];
sx q[2];
rz(1.4286263) q[2];
rz(0.96380487) q[3];
sx q[3];
rz(-1.0771841) q[3];
sx q[3];
rz(1.0296286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52255094) q[0];
sx q[0];
rz(-1.6864809) q[0];
sx q[0];
rz(1.2458941) q[0];
rz(3.030581) q[1];
sx q[1];
rz(-1.1975892) q[1];
sx q[1];
rz(-0.54661173) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3205991) q[0];
sx q[0];
rz(-0.59251596) q[0];
sx q[0];
rz(0.86235637) q[0];
rz(-pi) q[1];
x q[1];
rz(0.75066363) q[2];
sx q[2];
rz(-0.51807907) q[2];
sx q[2];
rz(2.2028365) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0134125) q[1];
sx q[1];
rz(-1.9617404) q[1];
sx q[1];
rz(-0.18545111) q[1];
x q[2];
rz(0.48874493) q[3];
sx q[3];
rz(-2.1663323) q[3];
sx q[3];
rz(-2.3084156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0299915) q[2];
sx q[2];
rz(-2.2932055) q[2];
sx q[2];
rz(1.6938422) q[2];
rz(1.1374121) q[3];
sx q[3];
rz(-1.8177989) q[3];
sx q[3];
rz(-2.494273) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85957134) q[0];
sx q[0];
rz(-2.8347926) q[0];
sx q[0];
rz(0.71722537) q[0];
rz(-1.2099129) q[1];
sx q[1];
rz(-2.8094493) q[1];
sx q[1];
rz(-2.4338914) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.087698547) q[0];
sx q[0];
rz(-2.2367034) q[0];
sx q[0];
rz(0.44184394) q[0];
x q[1];
rz(0.14214469) q[2];
sx q[2];
rz(-1.3488349) q[2];
sx q[2];
rz(-0.97505002) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3860491) q[1];
sx q[1];
rz(-1.2554597) q[1];
sx q[1];
rz(-0.83868933) q[1];
rz(-0.95146146) q[3];
sx q[3];
rz(-1.2776432) q[3];
sx q[3];
rz(0.64630634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5132961) q[2];
sx q[2];
rz(-0.6568903) q[2];
sx q[2];
rz(-0.90325242) q[2];
rz(-1.603027) q[3];
sx q[3];
rz(-0.86849803) q[3];
sx q[3];
rz(2.2911151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83508867) q[0];
sx q[0];
rz(-0.36515129) q[0];
sx q[0];
rz(-0.93602244) q[0];
rz(0.8159591) q[1];
sx q[1];
rz(-2.7201256) q[1];
sx q[1];
rz(1.0526007) q[1];
rz(0.46335285) q[2];
sx q[2];
rz(-1.9909161) q[2];
sx q[2];
rz(-1.9009895) q[2];
rz(-1.5296616) q[3];
sx q[3];
rz(-1.516468) q[3];
sx q[3];
rz(-0.80249912) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
