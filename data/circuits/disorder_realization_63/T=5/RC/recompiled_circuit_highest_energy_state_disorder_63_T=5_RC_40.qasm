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
rz(-1.0492078) q[0];
sx q[0];
rz(7.0424289) q[0];
sx q[0];
rz(9.0483604) q[0];
rz(1.5578101) q[1];
sx q[1];
rz(-1.0695142) q[1];
sx q[1];
rz(0.95626107) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96891785) q[0];
sx q[0];
rz(-2.4233074) q[0];
sx q[0];
rz(-1.0802313) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3772493) q[2];
sx q[2];
rz(-1.1551394) q[2];
sx q[2];
rz(0.035482835) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.2884644) q[1];
sx q[1];
rz(-0.83544105) q[1];
sx q[1];
rz(0.16448994) q[1];
rz(2.4231367) q[3];
sx q[3];
rz(-2.7877533) q[3];
sx q[3];
rz(-1.2350936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3236397) q[2];
sx q[2];
rz(-0.91922593) q[2];
sx q[2];
rz(1.2032571) q[2];
rz(-1.0407) q[3];
sx q[3];
rz(-2.2668656) q[3];
sx q[3];
rz(0.11999764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0284718) q[0];
sx q[0];
rz(-2.5682243) q[0];
sx q[0];
rz(-2.0290802) q[0];
rz(1.6587229) q[1];
sx q[1];
rz(-2.5506546) q[1];
sx q[1];
rz(-0.15731752) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33799115) q[0];
sx q[0];
rz(-2.4287619) q[0];
sx q[0];
rz(2.8719851) q[0];
rz(-1.4837166) q[2];
sx q[2];
rz(-1.198999) q[2];
sx q[2];
rz(-1.4390989) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1575732) q[1];
sx q[1];
rz(-1.3270151) q[1];
sx q[1];
rz(0.88884377) q[1];
rz(-pi) q[2];
rz(1.6460765) q[3];
sx q[3];
rz(-1.7930248) q[3];
sx q[3];
rz(2.8622975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4979672) q[2];
sx q[2];
rz(-2.6622541) q[2];
sx q[2];
rz(0.41110006) q[2];
rz(2.5655365) q[3];
sx q[3];
rz(-2.1274302) q[3];
sx q[3];
rz(0.99803734) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3882947) q[0];
sx q[0];
rz(-1.0030712) q[0];
sx q[0];
rz(-0.74111795) q[0];
rz(-1.4499715) q[1];
sx q[1];
rz(-1.3131817) q[1];
sx q[1];
rz(2.5695739) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2021671) q[0];
sx q[0];
rz(-1.8644445) q[0];
sx q[0];
rz(2.9333326) q[0];
rz(1.3640312) q[2];
sx q[2];
rz(-1.8305147) q[2];
sx q[2];
rz(-2.5842359) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1831844) q[1];
sx q[1];
rz(-2.5377877) q[1];
sx q[1];
rz(0.81012695) q[1];
rz(-pi) q[2];
rz(-2.1628863) q[3];
sx q[3];
rz(-2.7436069) q[3];
sx q[3];
rz(0.0090816895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.584562) q[2];
sx q[2];
rz(-1.1921554) q[2];
sx q[2];
rz(-2.3438047) q[2];
rz(-1.3095464) q[3];
sx q[3];
rz(-1.8833501) q[3];
sx q[3];
rz(1.4057188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13084594) q[0];
sx q[0];
rz(-2.2444785) q[0];
sx q[0];
rz(2.77453) q[0];
rz(-1.2182073) q[1];
sx q[1];
rz(-1.7299078) q[1];
sx q[1];
rz(2.9201115) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1786792) q[0];
sx q[0];
rz(-1.7811462) q[0];
sx q[0];
rz(1.6400997) q[0];
rz(-pi) q[1];
rz(1.70752) q[2];
sx q[2];
rz(-1.2061937) q[2];
sx q[2];
rz(-0.89791384) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.023497009) q[1];
sx q[1];
rz(-2.6129236) q[1];
sx q[1];
rz(-2.4882727) q[1];
rz(2.6148197) q[3];
sx q[3];
rz(-1.793141) q[3];
sx q[3];
rz(2.6759345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9941142) q[2];
sx q[2];
rz(-1.1722112) q[2];
sx q[2];
rz(-1.3383024) q[2];
rz(-2.639468) q[3];
sx q[3];
rz(-0.10052557) q[3];
sx q[3];
rz(-2.8978735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62622708) q[0];
sx q[0];
rz(-0.60972649) q[0];
sx q[0];
rz(1.6060265) q[0];
rz(3.0961127) q[1];
sx q[1];
rz(-1.3810424) q[1];
sx q[1];
rz(2.6754726) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.461994) q[0];
sx q[0];
rz(-1.3211622) q[0];
sx q[0];
rz(-0.14139463) q[0];
rz(-pi) q[1];
x q[1];
rz(1.999249) q[2];
sx q[2];
rz(-2.4523297) q[2];
sx q[2];
rz(1.9826629) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1093394) q[1];
sx q[1];
rz(-1.4795327) q[1];
sx q[1];
rz(1.5534459) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1640571) q[3];
sx q[3];
rz(-1.007742) q[3];
sx q[3];
rz(-2.7138591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8972299) q[2];
sx q[2];
rz(-1.676061) q[2];
sx q[2];
rz(-0.7828632) q[2];
rz(-3.1252981) q[3];
sx q[3];
rz(-1.8330845) q[3];
sx q[3];
rz(2.079594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8472854) q[0];
sx q[0];
rz(-2.6562302) q[0];
sx q[0];
rz(-2.7622727) q[0];
rz(0.30917057) q[1];
sx q[1];
rz(-1.9425756) q[1];
sx q[1];
rz(-0.22383037) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7302542) q[0];
sx q[0];
rz(-2.3650432) q[0];
sx q[0];
rz(1.1194487) q[0];
x q[1];
rz(1.7621222) q[2];
sx q[2];
rz(-1.680964) q[2];
sx q[2];
rz(-2.3913111) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6675615) q[1];
sx q[1];
rz(-1.4646475) q[1];
sx q[1];
rz(-0.1072704) q[1];
rz(-pi) q[2];
rz(-2.4357314) q[3];
sx q[3];
rz(-1.1339203) q[3];
sx q[3];
rz(-2.049751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0763187) q[2];
sx q[2];
rz(-1.9741917) q[2];
sx q[2];
rz(2.8889612) q[2];
rz(2.1624055) q[3];
sx q[3];
rz(-2.2604209) q[3];
sx q[3];
rz(2.7948936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
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
rz(2.613649) q[0];
sx q[0];
rz(-2.376463) q[0];
sx q[0];
rz(-1.8494404) q[0];
rz(-0.53391236) q[1];
sx q[1];
rz(-1.6894059) q[1];
sx q[1];
rz(0.36010489) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0468217) q[0];
sx q[0];
rz(-1.5015389) q[0];
sx q[0];
rz(-1.9317301) q[0];
rz(-pi) q[1];
rz(2.656237) q[2];
sx q[2];
rz(-1.8543285) q[2];
sx q[2];
rz(-2.527371) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.26038995) q[1];
sx q[1];
rz(-2.5298928) q[1];
sx q[1];
rz(2.7124497) q[1];
rz(-1.6735745) q[3];
sx q[3];
rz(-1.6037827) q[3];
sx q[3];
rz(1.7793293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.99122938) q[2];
sx q[2];
rz(-1.6232619) q[2];
sx q[2];
rz(2.4156127) q[2];
rz(-2.7109072) q[3];
sx q[3];
rz(-2.3613598) q[3];
sx q[3];
rz(-2.6881257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3579269) q[0];
sx q[0];
rz(-1.6538606) q[0];
sx q[0];
rz(2.3554392) q[0];
rz(-2.9642726) q[1];
sx q[1];
rz(-1.0847849) q[1];
sx q[1];
rz(-1.0160944) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35619104) q[0];
sx q[0];
rz(-1.2895661) q[0];
sx q[0];
rz(0.10939557) q[0];
rz(-pi) q[1];
rz(1.7012322) q[2];
sx q[2];
rz(-1.36048) q[2];
sx q[2];
rz(0.73784251) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1397649) q[1];
sx q[1];
rz(-1.5878126) q[1];
sx q[1];
rz(-0.66405343) q[1];
x q[2];
rz(2.3887064) q[3];
sx q[3];
rz(-1.2883923) q[3];
sx q[3];
rz(-1.2644067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5440386) q[2];
sx q[2];
rz(-2.1860213) q[2];
sx q[2];
rz(-0.7698783) q[2];
rz(2.9652413) q[3];
sx q[3];
rz(-1.4498815) q[3];
sx q[3];
rz(-2.8281853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65799323) q[0];
sx q[0];
rz(-0.083881065) q[0];
sx q[0];
rz(1.0461079) q[0];
rz(1.5484035) q[1];
sx q[1];
rz(-1.7240588) q[1];
sx q[1];
rz(-1.2215325) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0964805) q[0];
sx q[0];
rz(-1.4661705) q[0];
sx q[0];
rz(1.2136995) q[0];
rz(-pi) q[1];
x q[1];
rz(0.1510001) q[2];
sx q[2];
rz(-1.2556453) q[2];
sx q[2];
rz(-1.1282819) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0960779) q[1];
sx q[1];
rz(-1.7723262) q[1];
sx q[1];
rz(-1.3526826) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.043552355) q[3];
sx q[3];
rz(-1.5756902) q[3];
sx q[3];
rz(-2.2088449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1239502) q[2];
sx q[2];
rz(-1.7240179) q[2];
sx q[2];
rz(1.1851912) q[2];
rz(2.8017398) q[3];
sx q[3];
rz(-2.161882) q[3];
sx q[3];
rz(3.0237696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-1.3621984) q[0];
sx q[0];
rz(-1.8147991) q[0];
sx q[0];
rz(0.69778824) q[0];
rz(0.70480529) q[1];
sx q[1];
rz(-1.1230527) q[1];
sx q[1];
rz(-2.8885081) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6694723) q[0];
sx q[0];
rz(-1.6447495) q[0];
sx q[0];
rz(-1.8042613) q[0];
x q[1];
rz(0.65480755) q[2];
sx q[2];
rz(-1.1455481) q[2];
sx q[2];
rz(3.0523093) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.608111) q[1];
sx q[1];
rz(-2.2072756) q[1];
sx q[1];
rz(0.57884207) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3138103) q[3];
sx q[3];
rz(-0.8567613) q[3];
sx q[3];
rz(1.0257126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.42085984) q[2];
sx q[2];
rz(-1.3497817) q[2];
sx q[2];
rz(-0.68812686) q[2];
rz(-2.1963035) q[3];
sx q[3];
rz(-1.9663845) q[3];
sx q[3];
rz(0.92760408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5505493) q[0];
sx q[0];
rz(-1.782438) q[0];
sx q[0];
rz(1.8989643) q[0];
rz(-0.91724829) q[1];
sx q[1];
rz(-1.6684253) q[1];
sx q[1];
rz(-2.565276) q[1];
rz(2.5653432) q[2];
sx q[2];
rz(-1.2628308) q[2];
sx q[2];
rz(-2.2028883) q[2];
rz(1.3618906) q[3];
sx q[3];
rz(-2.7068797) q[3];
sx q[3];
rz(-2.5754365) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
