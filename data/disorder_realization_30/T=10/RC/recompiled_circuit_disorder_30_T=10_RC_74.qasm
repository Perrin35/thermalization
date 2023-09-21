OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.49139872) q[0];
sx q[0];
rz(-0.2645275) q[0];
sx q[0];
rz(-0.39443031) q[0];
rz(0.0061622942) q[1];
sx q[1];
rz(-0.34024629) q[1];
sx q[1];
rz(1.9415829) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29205706) q[0];
sx q[0];
rz(-1.6186876) q[0];
sx q[0];
rz(0.0056664771) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1985487) q[2];
sx q[2];
rz(-2.606751) q[2];
sx q[2];
rz(-2.5778071) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1026417) q[1];
sx q[1];
rz(-1.1167553) q[1];
sx q[1];
rz(0.42788831) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6148189) q[3];
sx q[3];
rz(-1.1929973) q[3];
sx q[3];
rz(1.6954741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3068984) q[2];
sx q[2];
rz(-1.4935741) q[2];
sx q[2];
rz(-0.28960323) q[2];
rz(0.87537193) q[3];
sx q[3];
rz(-2.1353728) q[3];
sx q[3];
rz(-0.059710596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5263379) q[0];
sx q[0];
rz(-2.2741788) q[0];
sx q[0];
rz(0.34399024) q[0];
rz(-0.084331766) q[1];
sx q[1];
rz(-0.66939676) q[1];
sx q[1];
rz(-1.3551691) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7536613) q[0];
sx q[0];
rz(-0.61457115) q[0];
sx q[0];
rz(-1.6835) q[0];
rz(-pi) q[1];
rz(1.6927035) q[2];
sx q[2];
rz(-1.4093471) q[2];
sx q[2];
rz(-3.0485857) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1234839) q[1];
sx q[1];
rz(-1.5068441) q[1];
sx q[1];
rz(1.121184) q[1];
rz(-pi) q[2];
rz(1.387272) q[3];
sx q[3];
rz(-1.4944544) q[3];
sx q[3];
rz(-2.9294088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.29558674) q[2];
sx q[2];
rz(-2.3082374) q[2];
sx q[2];
rz(0.53768349) q[2];
rz(-2.6387571) q[3];
sx q[3];
rz(-0.42481315) q[3];
sx q[3];
rz(-1.0104377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4780592) q[0];
sx q[0];
rz(-2.6038267) q[0];
sx q[0];
rz(-0.64087254) q[0];
rz(-0.74869853) q[1];
sx q[1];
rz(-1.0083895) q[1];
sx q[1];
rz(-1.0650939) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.070806064) q[0];
sx q[0];
rz(-2.7822128) q[0];
sx q[0];
rz(-1.5887567) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2421354) q[2];
sx q[2];
rz(-0.81842917) q[2];
sx q[2];
rz(0.98086548) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0712191) q[1];
sx q[1];
rz(-2.1871236) q[1];
sx q[1];
rz(3.0973869) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.84824003) q[3];
sx q[3];
rz(-1.7602929) q[3];
sx q[3];
rz(-3.0998067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.37725267) q[2];
sx q[2];
rz(-1.179402) q[2];
sx q[2];
rz(-2.4242145) q[2];
rz(-2.453089) q[3];
sx q[3];
rz(-0.62763667) q[3];
sx q[3];
rz(2.8440516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3142969) q[0];
sx q[0];
rz(-1.2001487) q[0];
sx q[0];
rz(0.24969077) q[0];
rz(2.1266134) q[1];
sx q[1];
rz(-2.8453638) q[1];
sx q[1];
rz(-0.011118523) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5549094) q[0];
sx q[0];
rz(-1.1197829) q[0];
sx q[0];
rz(-0.45278544) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.87984933) q[2];
sx q[2];
rz(-0.86266154) q[2];
sx q[2];
rz(-0.81530064) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3938155) q[1];
sx q[1];
rz(-1.2443466) q[1];
sx q[1];
rz(1.0541037) q[1];
rz(2.6075881) q[3];
sx q[3];
rz(-0.36557331) q[3];
sx q[3];
rz(-0.064985736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6461688) q[2];
sx q[2];
rz(-1.255722) q[2];
sx q[2];
rz(0.28309506) q[2];
rz(0.66343534) q[3];
sx q[3];
rz(-0.59316558) q[3];
sx q[3];
rz(0.88808131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11113142) q[0];
sx q[0];
rz(-2.9009394) q[0];
sx q[0];
rz(2.8097613) q[0];
rz(-0.49452531) q[1];
sx q[1];
rz(-1.8437513) q[1];
sx q[1];
rz(-1.8146851) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6700232) q[0];
sx q[0];
rz(-1.7243392) q[0];
sx q[0];
rz(-0.75385401) q[0];
rz(-2.5577776) q[2];
sx q[2];
rz(-0.8875672) q[2];
sx q[2];
rz(-2.3655287) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.27969589) q[1];
sx q[1];
rz(-1.3394636) q[1];
sx q[1];
rz(0.33613236) q[1];
x q[2];
rz(-0.29019659) q[3];
sx q[3];
rz(-1.962933) q[3];
sx q[3];
rz(-0.41662595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1422687) q[2];
sx q[2];
rz(-2.6596255) q[2];
sx q[2];
rz(1.2456606) q[2];
rz(1.8866395) q[3];
sx q[3];
rz(-2.3001223) q[3];
sx q[3];
rz(0.83827034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6923043) q[0];
sx q[0];
rz(-3.1384387) q[0];
sx q[0];
rz(-2.4601049) q[0];
rz(2.9340414) q[1];
sx q[1];
rz(-2.6782268) q[1];
sx q[1];
rz(-2.025827) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.241248) q[0];
sx q[0];
rz(-1.9755409) q[0];
sx q[0];
rz(0.83165283) q[0];
rz(-pi) q[1];
rz(2.7115466) q[2];
sx q[2];
rz(-2.1950245) q[2];
sx q[2];
rz(-0.66832322) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.15397729) q[1];
sx q[1];
rz(-0.9824285) q[1];
sx q[1];
rz(-1.5495367) q[1];
x q[2];
rz(-0.22478215) q[3];
sx q[3];
rz(-0.68192476) q[3];
sx q[3];
rz(-0.080760591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9289124) q[2];
sx q[2];
rz(-1.2968411) q[2];
sx q[2];
rz(-0.35432717) q[2];
rz(-0.3195233) q[3];
sx q[3];
rz(-1.9890246) q[3];
sx q[3];
rz(-0.37187809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1324683) q[0];
sx q[0];
rz(-2.9102944) q[0];
sx q[0];
rz(-0.67434597) q[0];
rz(2.0293503) q[1];
sx q[1];
rz(-0.66450417) q[1];
sx q[1];
rz(0.56232125) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4599265) q[0];
sx q[0];
rz(-0.66837464) q[0];
sx q[0];
rz(1.5468803) q[0];
rz(-pi) q[1];
rz(-0.28568761) q[2];
sx q[2];
rz(-0.34491587) q[2];
sx q[2];
rz(-0.069318511) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6866236) q[1];
sx q[1];
rz(-1.3775871) q[1];
sx q[1];
rz(-1.7186233) q[1];
rz(-pi) q[2];
rz(-2.167422) q[3];
sx q[3];
rz(-1.4985634) q[3];
sx q[3];
rz(0.40856397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.101863) q[2];
sx q[2];
rz(-2.1477284) q[2];
sx q[2];
rz(-2.8015461) q[2];
rz(0.21751054) q[3];
sx q[3];
rz(-2.1931931) q[3];
sx q[3];
rz(0.31869179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44889221) q[0];
sx q[0];
rz(-3.0954439) q[0];
sx q[0];
rz(-0.39644077) q[0];
rz(0.13892826) q[1];
sx q[1];
rz(-0.46008343) q[1];
sx q[1];
rz(1.5213535) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.053145807) q[0];
sx q[0];
rz(-1.7674812) q[0];
sx q[0];
rz(-1.501207) q[0];
x q[1];
rz(1.0993768) q[2];
sx q[2];
rz(-1.2087012) q[2];
sx q[2];
rz(-0.21381703) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7594787) q[1];
sx q[1];
rz(-1.5112875) q[1];
sx q[1];
rz(-0.77854034) q[1];
x q[2];
rz(1.7117386) q[3];
sx q[3];
rz(-0.68416506) q[3];
sx q[3];
rz(-1.14389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5788995) q[2];
sx q[2];
rz(-2.0843299) q[2];
sx q[2];
rz(-2.8472624) q[2];
rz(2.0108022) q[3];
sx q[3];
rz(-1.7629938) q[3];
sx q[3];
rz(2.14595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49333736) q[0];
sx q[0];
rz(-2.2264037) q[0];
sx q[0];
rz(2.7822568) q[0];
rz(2.1954779) q[1];
sx q[1];
rz(-2.7455536) q[1];
sx q[1];
rz(2.8709581) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5198869) q[0];
sx q[0];
rz(-1.6019551) q[0];
sx q[0];
rz(3.1057538) q[0];
rz(-2.052202) q[2];
sx q[2];
rz(-2.319616) q[2];
sx q[2];
rz(0.35441986) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2808025) q[1];
sx q[1];
rz(-1.6165015) q[1];
sx q[1];
rz(1.526282) q[1];
x q[2];
rz(-2.4828033) q[3];
sx q[3];
rz(-0.82595347) q[3];
sx q[3];
rz(-0.39289075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3140807) q[2];
sx q[2];
rz(-2.9120047) q[2];
sx q[2];
rz(-2.6861526) q[2];
rz(0.81196249) q[3];
sx q[3];
rz(-2.2596695) q[3];
sx q[3];
rz(0.5493831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6311326) q[0];
sx q[0];
rz(-1.5083418) q[0];
sx q[0];
rz(-2.4023138) q[0];
rz(-0.23070681) q[1];
sx q[1];
rz(-2.6735327) q[1];
sx q[1];
rz(-2.646692) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1586944) q[0];
sx q[0];
rz(-1.532908) q[0];
sx q[0];
rz(2.057103) q[0];
x q[1];
rz(-1.5779183) q[2];
sx q[2];
rz(-1.2429503) q[2];
sx q[2];
rz(-0.42051007) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7086664) q[1];
sx q[1];
rz(-2.4927944) q[1];
sx q[1];
rz(1.385958) q[1];
rz(-2.2675715) q[3];
sx q[3];
rz(-1.6439983) q[3];
sx q[3];
rz(-2.8332568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0080002) q[2];
sx q[2];
rz(-2.0501037) q[2];
sx q[2];
rz(-0.32785329) q[2];
rz(-2.949529) q[3];
sx q[3];
rz(-0.24774194) q[3];
sx q[3];
rz(1.0333992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1223758) q[0];
sx q[0];
rz(-1.6288971) q[0];
sx q[0];
rz(-1.6678641) q[0];
rz(0.83256759) q[1];
sx q[1];
rz(-1.7201798) q[1];
sx q[1];
rz(-1.7725772) q[1];
rz(1.9033296) q[2];
sx q[2];
rz(-2.5730972) q[2];
sx q[2];
rz(2.2039883) q[2];
rz(-0.26631793) q[3];
sx q[3];
rz(-1.3238293) q[3];
sx q[3];
rz(-2.7662591) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
