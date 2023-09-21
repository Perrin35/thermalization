OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.49178034) q[0];
sx q[0];
rz(-2.8556813) q[0];
sx q[0];
rz(-0.51529348) q[0];
rz(4.4858785) q[1];
sx q[1];
rz(2.9872515) q[1];
sx q[1];
rz(6.8607688) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89864697) q[0];
sx q[0];
rz(-2.4624914) q[0];
sx q[0];
rz(0.26429096) q[0];
rz(0.53519997) q[2];
sx q[2];
rz(-2.0173965) q[2];
sx q[2];
rz(-1.610178) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.56596245) q[1];
sx q[1];
rz(-1.9735676) q[1];
sx q[1];
rz(-0.29213841) q[1];
rz(-2.3732784) q[3];
sx q[3];
rz(-2.3294805) q[3];
sx q[3];
rz(2.1384359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.704533) q[2];
sx q[2];
rz(-1.5960863) q[2];
sx q[2];
rz(-2.4543767) q[2];
rz(-2.1263188) q[3];
sx q[3];
rz(-1.7679368) q[3];
sx q[3];
rz(-0.12250531) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17094831) q[0];
sx q[0];
rz(-2.0630554) q[0];
sx q[0];
rz(1.8815536) q[0];
rz(-1.0062224) q[1];
sx q[1];
rz(-2.1496014) q[1];
sx q[1];
rz(0.84567436) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.082367912) q[0];
sx q[0];
rz(-2.0031643) q[0];
sx q[0];
rz(2.5655377) q[0];
rz(-pi) q[1];
rz(1.4405865) q[2];
sx q[2];
rz(-1.7034334) q[2];
sx q[2];
rz(2.8108033) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.37296346) q[1];
sx q[1];
rz(-1.0781204) q[1];
sx q[1];
rz(1.3586587) q[1];
x q[2];
rz(2.2250697) q[3];
sx q[3];
rz(-2.407981) q[3];
sx q[3];
rz(-0.77378002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3530897) q[2];
sx q[2];
rz(-0.22558364) q[2];
sx q[2];
rz(2.6611924) q[2];
rz(-1.7885615) q[3];
sx q[3];
rz(-2.0856817) q[3];
sx q[3];
rz(-1.1876748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1903494) q[0];
sx q[0];
rz(-2.9028063) q[0];
sx q[0];
rz(-2.3685266) q[0];
rz(3.0103325) q[1];
sx q[1];
rz(-1.2845598) q[1];
sx q[1];
rz(-2.0551596) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4341136) q[0];
sx q[0];
rz(-1.1093603) q[0];
sx q[0];
rz(0.001860851) q[0];
x q[1];
rz(-1.7227206) q[2];
sx q[2];
rz(-0.65290367) q[2];
sx q[2];
rz(0.34130794) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.8311027) q[1];
sx q[1];
rz(-1.576014) q[1];
sx q[1];
rz(2.8527841) q[1];
x q[2];
rz(-2.9680786) q[3];
sx q[3];
rz(-1.1713235) q[3];
sx q[3];
rz(0.80843335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1290258) q[2];
sx q[2];
rz(-1.4386703) q[2];
sx q[2];
rz(1.3712937) q[2];
rz(2.7584372) q[3];
sx q[3];
rz(-1.2569191) q[3];
sx q[3];
rz(-0.80254054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6894158) q[0];
sx q[0];
rz(-1.8912264) q[0];
sx q[0];
rz(-3.0932328) q[0];
rz(-0.16391779) q[1];
sx q[1];
rz(-2.7719031) q[1];
sx q[1];
rz(-1.6960467) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.018054124) q[0];
sx q[0];
rz(-2.2211383) q[0];
sx q[0];
rz(1.7975848) q[0];
rz(-pi) q[1];
rz(-1.4807329) q[2];
sx q[2];
rz(-1.3552595) q[2];
sx q[2];
rz(2.4243674) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.85019894) q[1];
sx q[1];
rz(-1.7899917) q[1];
sx q[1];
rz(-1.6367957) q[1];
rz(-pi) q[2];
rz(2.0914145) q[3];
sx q[3];
rz(-1.8225267) q[3];
sx q[3];
rz(-2.7057735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.066102862) q[2];
sx q[2];
rz(-1.3572071) q[2];
sx q[2];
rz(-1.0243105) q[2];
rz(-1.5284437) q[3];
sx q[3];
rz(-1.5214835) q[3];
sx q[3];
rz(0.23322341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7016474) q[0];
sx q[0];
rz(-0.82413903) q[0];
sx q[0];
rz(-1.8540927) q[0];
rz(-2.8225186) q[1];
sx q[1];
rz(-1.5417475) q[1];
sx q[1];
rz(-2.2873926) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.094164205) q[0];
sx q[0];
rz(-1.578997) q[0];
sx q[0];
rz(-0.8568944) q[0];
x q[1];
rz(0.87128432) q[2];
sx q[2];
rz(-1.185002) q[2];
sx q[2];
rz(2.3988349) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2367868) q[1];
sx q[1];
rz(-0.85610897) q[1];
sx q[1];
rz(2.3805815) q[1];
rz(2.5913127) q[3];
sx q[3];
rz(-2.0645421) q[3];
sx q[3];
rz(-0.20875904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.95191082) q[2];
sx q[2];
rz(-2.5482735) q[2];
sx q[2];
rz(2.5642776) q[2];
rz(0.50950766) q[3];
sx q[3];
rz(-2.7039492) q[3];
sx q[3];
rz(-0.90665162) q[3];
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
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0034870738) q[0];
sx q[0];
rz(-2.0697937) q[0];
sx q[0];
rz(-0.072120897) q[0];
rz(2.0347319) q[1];
sx q[1];
rz(-0.51263428) q[1];
sx q[1];
rz(-3.0153826) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0059347) q[0];
sx q[0];
rz(-1.6023484) q[0];
sx q[0];
rz(0.21261442) q[0];
x q[1];
rz(1.1926786) q[2];
sx q[2];
rz(-2.3092804) q[2];
sx q[2];
rz(-0.2371012) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.63657657) q[1];
sx q[1];
rz(-2.2762183) q[1];
sx q[1];
rz(2.8935863) q[1];
rz(-pi) q[2];
x q[2];
rz(1.24228) q[3];
sx q[3];
rz(-1.3887172) q[3];
sx q[3];
rz(-1.8329221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.90298992) q[2];
sx q[2];
rz(-2.949252) q[2];
sx q[2];
rz(-2.3664756) q[2];
rz(0.827968) q[3];
sx q[3];
rz(-2.8505846) q[3];
sx q[3];
rz(-1.1221788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59259748) q[0];
sx q[0];
rz(-2.7153375) q[0];
sx q[0];
rz(-0.098408498) q[0];
rz(-1.1920284) q[1];
sx q[1];
rz(-1.3339309) q[1];
sx q[1];
rz(-2.5820406) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3768809) q[0];
sx q[0];
rz(-1.6983713) q[0];
sx q[0];
rz(2.5401831) q[0];
rz(-0.26142188) q[2];
sx q[2];
rz(-1.6209507) q[2];
sx q[2];
rz(-2.9068771) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.93058649) q[1];
sx q[1];
rz(-1.483327) q[1];
sx q[1];
rz(1.4409815) q[1];
rz(-pi) q[2];
rz(-0.26569326) q[3];
sx q[3];
rz(-2.4943647) q[3];
sx q[3];
rz(-0.94097394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6825535) q[2];
sx q[2];
rz(-1.8457396) q[2];
sx q[2];
rz(-2.7977978) q[2];
rz(-2.5750459) q[3];
sx q[3];
rz(-0.44851258) q[3];
sx q[3];
rz(2.6678273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.375181) q[0];
sx q[0];
rz(-1.3324998) q[0];
sx q[0];
rz(2.4108316) q[0];
rz(2.9991951) q[1];
sx q[1];
rz(-1.8715033) q[1];
sx q[1];
rz(-0.87160814) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7155834) q[0];
sx q[0];
rz(-1.4953488) q[0];
sx q[0];
rz(-1.5619318) q[0];
rz(-pi) q[1];
rz(1.4191188) q[2];
sx q[2];
rz(-1.9722087) q[2];
sx q[2];
rz(2.7733208) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.1069113) q[1];
sx q[1];
rz(-2.2195663) q[1];
sx q[1];
rz(1.1120863) q[1];
rz(-0.96846795) q[3];
sx q[3];
rz(-1.6864711) q[3];
sx q[3];
rz(-2.065421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4153851) q[2];
sx q[2];
rz(-1.0691079) q[2];
sx q[2];
rz(-1.139572) q[2];
rz(1.6428044) q[3];
sx q[3];
rz(-2.747624) q[3];
sx q[3];
rz(-2.22877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20294872) q[0];
sx q[0];
rz(-1.4511755) q[0];
sx q[0];
rz(-1.9198445) q[0];
rz(-2.9755759) q[1];
sx q[1];
rz(-1.32042) q[1];
sx q[1];
rz(1.5244012) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.411392) q[0];
sx q[0];
rz(-0.087326614) q[0];
sx q[0];
rz(-1.9090396) q[0];
rz(-pi) q[1];
rz(1.7622856) q[2];
sx q[2];
rz(-2.0404437) q[2];
sx q[2];
rz(0.42275235) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7414654) q[1];
sx q[1];
rz(-1.2418081) q[1];
sx q[1];
rz(-2.9743183) q[1];
rz(2.5025326) q[3];
sx q[3];
rz(-1.3165054) q[3];
sx q[3];
rz(-1.9494848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.021585492) q[2];
sx q[2];
rz(-1.6759796) q[2];
sx q[2];
rz(2.7900556) q[2];
rz(-2.0848138) q[3];
sx q[3];
rz(-2.6119699) q[3];
sx q[3];
rz(-2.3969011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64365023) q[0];
sx q[0];
rz(-0.90181667) q[0];
sx q[0];
rz(1.836401) q[0];
rz(-2.7611043) q[1];
sx q[1];
rz(-1.0419798) q[1];
sx q[1];
rz(2.8881853) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2449269) q[0];
sx q[0];
rz(-1.9650808) q[0];
sx q[0];
rz(-0.11163296) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.16742736) q[2];
sx q[2];
rz(-1.2360459) q[2];
sx q[2];
rz(1.6850922) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7518172) q[1];
sx q[1];
rz(-2.7174065) q[1];
sx q[1];
rz(2.5245689) q[1];
rz(-pi) q[2];
rz(1.1425584) q[3];
sx q[3];
rz(-1.0034475) q[3];
sx q[3];
rz(0.67728562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5499251) q[2];
sx q[2];
rz(-2.2558236) q[2];
sx q[2];
rz(0.5029451) q[2];
rz(2.2425966) q[3];
sx q[3];
rz(-1.2939724) q[3];
sx q[3];
rz(-1.1635273) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9713365) q[0];
sx q[0];
rz(-1.5383056) q[0];
sx q[0];
rz(-2.8785895) q[0];
rz(-0.7111711) q[1];
sx q[1];
rz(-1.0881337) q[1];
sx q[1];
rz(1.7137391) q[1];
rz(1.3118369) q[2];
sx q[2];
rz(-1.8428409) q[2];
sx q[2];
rz(0.98696282) q[2];
rz(-0.9234225) q[3];
sx q[3];
rz(-2.2150726) q[3];
sx q[3];
rz(1.1005145) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];