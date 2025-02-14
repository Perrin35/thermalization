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
rz(-2.4969555) q[0];
sx q[0];
rz(-0.58965373) q[0];
sx q[0];
rz(2.0539505) q[0];
rz(-0.015406869) q[1];
sx q[1];
rz(3.5313731) q[1];
sx q[1];
rz(11.402147) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8293267) q[0];
sx q[0];
rz(-3.0115033) q[0];
sx q[0];
rz(-1.6072558) q[0];
x q[1];
rz(-2.4504775) q[2];
sx q[2];
rz(-1.8370312) q[2];
sx q[2];
rz(-0.51048265) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3642295) q[1];
sx q[1];
rz(-2.573488) q[1];
sx q[1];
rz(1.9757759) q[1];
x q[2];
rz(0.075957493) q[3];
sx q[3];
rz(-1.4531368) q[3];
sx q[3];
rz(-2.3880588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2414134) q[2];
sx q[2];
rz(-2.5799077) q[2];
sx q[2];
rz(0.64603311) q[2];
rz(0.25520405) q[3];
sx q[3];
rz(-2.6845158) q[3];
sx q[3];
rz(-1.3946165) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0274886) q[0];
sx q[0];
rz(-2.8059967) q[0];
sx q[0];
rz(-0.80701971) q[0];
rz(1.457816) q[1];
sx q[1];
rz(-2.5648263) q[1];
sx q[1];
rz(-1.797765) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3407077) q[0];
sx q[0];
rz(-1.721764) q[0];
sx q[0];
rz(-0.29125352) q[0];
rz(-0.41092964) q[2];
sx q[2];
rz(-2.0996473) q[2];
sx q[2];
rz(1.628132) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2652774) q[1];
sx q[1];
rz(-0.96926276) q[1];
sx q[1];
rz(2.5391891) q[1];
rz(3.0968118) q[3];
sx q[3];
rz(-0.84242994) q[3];
sx q[3];
rz(-1.7573259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1809711) q[2];
sx q[2];
rz(-1.0701067) q[2];
sx q[2];
rz(0.19372678) q[2];
rz(1.5242029) q[3];
sx q[3];
rz(-2.3834855) q[3];
sx q[3];
rz(-2.962842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5209565) q[0];
sx q[0];
rz(-2.9036324) q[0];
sx q[0];
rz(-0.97610193) q[0];
rz(2.2528265) q[1];
sx q[1];
rz(-0.74405324) q[1];
sx q[1];
rz(-2.9685453) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9219396) q[0];
sx q[0];
rz(-1.5679857) q[0];
sx q[0];
rz(-1.5332743) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6209774) q[2];
sx q[2];
rz(-1.4618502) q[2];
sx q[2];
rz(-1.9828988) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.060932) q[1];
sx q[1];
rz(-2.073604) q[1];
sx q[1];
rz(1.7601556) q[1];
rz(1.6794231) q[3];
sx q[3];
rz(-1.2501688) q[3];
sx q[3];
rz(-2.4214793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0024857) q[2];
sx q[2];
rz(-1.3287013) q[2];
sx q[2];
rz(2.4277182) q[2];
rz(0.21126963) q[3];
sx q[3];
rz(-2.1043089) q[3];
sx q[3];
rz(0.56536388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36770377) q[0];
sx q[0];
rz(-1.1834894) q[0];
sx q[0];
rz(2.0088038) q[0];
rz(2.6702113) q[1];
sx q[1];
rz(-1.2939204) q[1];
sx q[1];
rz(0.43101355) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16551183) q[0];
sx q[0];
rz(-1.3595194) q[0];
sx q[0];
rz(-1.6863281) q[0];
x q[1];
rz(1.9507031) q[2];
sx q[2];
rz(-1.2326733) q[2];
sx q[2];
rz(-2.9924711) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9543332) q[1];
sx q[1];
rz(-1.0403087) q[1];
sx q[1];
rz(-1.9217291) q[1];
x q[2];
rz(-2.1924344) q[3];
sx q[3];
rz(-1.7931058) q[3];
sx q[3];
rz(-2.3883826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6212578) q[2];
sx q[2];
rz(-0.73953491) q[2];
sx q[2];
rz(0.4062824) q[2];
rz(1.2064365) q[3];
sx q[3];
rz(-1.0932357) q[3];
sx q[3];
rz(-2.5549197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63221145) q[0];
sx q[0];
rz(-2.2659232) q[0];
sx q[0];
rz(0.32633728) q[0];
rz(2.1874766) q[1];
sx q[1];
rz(-2.4931144) q[1];
sx q[1];
rz(-0.023524806) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92281728) q[0];
sx q[0];
rz(-1.1975702) q[0];
sx q[0];
rz(1.1362227) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.40074375) q[2];
sx q[2];
rz(-2.4432428) q[2];
sx q[2];
rz(-1.9633479) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6885029) q[1];
sx q[1];
rz(-1.9132691) q[1];
sx q[1];
rz(-0.30330412) q[1];
x q[2];
rz(-2.708528) q[3];
sx q[3];
rz(-0.32363656) q[3];
sx q[3];
rz(1.8214846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.894459) q[2];
sx q[2];
rz(-0.45866141) q[2];
sx q[2];
rz(-0.66391724) q[2];
rz(2.2513921) q[3];
sx q[3];
rz(-1.592344) q[3];
sx q[3];
rz(0.012705407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7909872) q[0];
sx q[0];
rz(-0.43142879) q[0];
sx q[0];
rz(0.33082333) q[0];
rz(-0.36258969) q[1];
sx q[1];
rz(-1.9622842) q[1];
sx q[1];
rz(2.6258452) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73677255) q[0];
sx q[0];
rz(-1.7869925) q[0];
sx q[0];
rz(-1.3275073) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.547566) q[2];
sx q[2];
rz(-1.7286679) q[2];
sx q[2];
rz(1.8094339) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.95553229) q[1];
sx q[1];
rz(-1.6481676) q[1];
sx q[1];
rz(-2.2633228) q[1];
rz(-pi) q[2];
rz(-1.4384076) q[3];
sx q[3];
rz(-2.3623132) q[3];
sx q[3];
rz(0.25151032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.35046878) q[2];
sx q[2];
rz(-1.4501269) q[2];
sx q[2];
rz(0.80751944) q[2];
rz(3.0430072) q[3];
sx q[3];
rz(-0.8980631) q[3];
sx q[3];
rz(-1.0928104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43807855) q[0];
sx q[0];
rz(-0.57323891) q[0];
sx q[0];
rz(2.692063) q[0];
rz(-2.4564157) q[1];
sx q[1];
rz(-2.7339869) q[1];
sx q[1];
rz(2.5221241) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6934366) q[0];
sx q[0];
rz(-1.5991028) q[0];
sx q[0];
rz(-2.8319915) q[0];
rz(-pi) q[1];
x q[1];
rz(0.52823587) q[2];
sx q[2];
rz(-1.6124469) q[2];
sx q[2];
rz(-1.8245175) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0692783) q[1];
sx q[1];
rz(-2.2439146) q[1];
sx q[1];
rz(-2.0050383) q[1];
x q[2];
rz(2.9569217) q[3];
sx q[3];
rz(-1.7725367) q[3];
sx q[3];
rz(-1.4278442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5770136) q[2];
sx q[2];
rz(-1.4489633) q[2];
sx q[2];
rz(-2.9435834) q[2];
rz(1.3624582) q[3];
sx q[3];
rz(-2.8415316) q[3];
sx q[3];
rz(2.4411966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.144416) q[0];
sx q[0];
rz(-0.62189019) q[0];
sx q[0];
rz(1.8560334) q[0];
rz(-2.8649435) q[1];
sx q[1];
rz(-1.3686907) q[1];
sx q[1];
rz(3.0329774) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.016182) q[0];
sx q[0];
rz(-1.6526395) q[0];
sx q[0];
rz(-1.0684408) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0622903) q[2];
sx q[2];
rz(-2.1933043) q[2];
sx q[2];
rz(-2.9626737) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7759571) q[1];
sx q[1];
rz(-1.3564907) q[1];
sx q[1];
rz(-2.2909095) q[1];
rz(-1.2694938) q[3];
sx q[3];
rz(-2.7387894) q[3];
sx q[3];
rz(2.3294152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.30839977) q[2];
sx q[2];
rz(-1.1966285) q[2];
sx q[2];
rz(-1.2742554) q[2];
rz(1.3537004) q[3];
sx q[3];
rz(-2.7022868) q[3];
sx q[3];
rz(-0.86658365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9891147) q[0];
sx q[0];
rz(-0.72565961) q[0];
sx q[0];
rz(-0.12055483) q[0];
rz(0.74147916) q[1];
sx q[1];
rz(-1.9375786) q[1];
sx q[1];
rz(3.0659058) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78223373) q[0];
sx q[0];
rz(-1.7348716) q[0];
sx q[0];
rz(0.30996451) q[0];
rz(-pi) q[1];
rz(-1.7255177) q[2];
sx q[2];
rz(-0.46636367) q[2];
sx q[2];
rz(1.6365964) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.32966742) q[1];
sx q[1];
rz(-0.65663785) q[1];
sx q[1];
rz(0.041402264) q[1];
rz(1.754154) q[3];
sx q[3];
rz(-2.1638738) q[3];
sx q[3];
rz(2.0545404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.072448298) q[2];
sx q[2];
rz(-0.88280237) q[2];
sx q[2];
rz(0.34633386) q[2];
rz(2.196178) q[3];
sx q[3];
rz(-1.8190705) q[3];
sx q[3];
rz(2.1944428) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.866975) q[0];
sx q[0];
rz(-1.9434384) q[0];
sx q[0];
rz(-2.6937038) q[0];
rz(0.92944324) q[1];
sx q[1];
rz(-1.7421236) q[1];
sx q[1];
rz(2.4670752) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68088082) q[0];
sx q[0];
rz(-1.370791) q[0];
sx q[0];
rz(0.77447156) q[0];
rz(-pi) q[1];
rz(0.10402502) q[2];
sx q[2];
rz(-1.3859704) q[2];
sx q[2];
rz(1.8341174) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.062194) q[1];
sx q[1];
rz(-0.17339686) q[1];
sx q[1];
rz(-1.5773415) q[1];
rz(-3.0661656) q[3];
sx q[3];
rz(-1.9884681) q[3];
sx q[3];
rz(3.053523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.21005361) q[2];
sx q[2];
rz(-1.9320107) q[2];
sx q[2];
rz(-2.8741969) q[2];
rz(1.1072655) q[3];
sx q[3];
rz(-2.3591154) q[3];
sx q[3];
rz(0.84899181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4103107) q[0];
sx q[0];
rz(-1.5189497) q[0];
sx q[0];
rz(2.0547163) q[0];
rz(-2.3375753) q[1];
sx q[1];
rz(-1.6537279) q[1];
sx q[1];
rz(-2.0860685) q[1];
rz(-3.0699365) q[2];
sx q[2];
rz(-1.717567) q[2];
sx q[2];
rz(-1.4302413) q[2];
rz(-1.4645529) q[3];
sx q[3];
rz(-2.4709656) q[3];
sx q[3];
rz(-2.7460302) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
