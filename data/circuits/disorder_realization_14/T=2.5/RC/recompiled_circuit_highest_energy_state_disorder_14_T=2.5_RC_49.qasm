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
rz(2.4013588) q[0];
sx q[0];
rz(-1.6594247) q[0];
sx q[0];
rz(-2.8066714) q[0];
rz(0.51796335) q[1];
sx q[1];
rz(5.2809102) q[1];
sx q[1];
rz(10.032293) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.018463919) q[0];
sx q[0];
rz(-2.0573318) q[0];
sx q[0];
rz(-2.8790265) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0614228) q[2];
sx q[2];
rz(-1.6730089) q[2];
sx q[2];
rz(-1.4776023) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.059587) q[1];
sx q[1];
rz(-1.730189) q[1];
sx q[1];
rz(0.48377796) q[1];
x q[2];
rz(-2.257454) q[3];
sx q[3];
rz(-0.55301412) q[3];
sx q[3];
rz(0.64698863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6174378) q[2];
sx q[2];
rz(-1.2158771) q[2];
sx q[2];
rz(-2.4784135) q[2];
rz(-3.0607306) q[3];
sx q[3];
rz(-2.9359449) q[3];
sx q[3];
rz(1.1801571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2422159) q[0];
sx q[0];
rz(-1.7689393) q[0];
sx q[0];
rz(0.85897613) q[0];
rz(1.8513177) q[1];
sx q[1];
rz(-1.4651508) q[1];
sx q[1];
rz(-1.7346409) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0075571) q[0];
sx q[0];
rz(-1.3430183) q[0];
sx q[0];
rz(-2.0640949) q[0];
rz(-pi) q[1];
rz(-0.036769899) q[2];
sx q[2];
rz(-1.2388133) q[2];
sx q[2];
rz(0.22184243) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.21995658) q[1];
sx q[1];
rz(-1.0590648) q[1];
sx q[1];
rz(0.10351609) q[1];
x q[2];
rz(2.5464393) q[3];
sx q[3];
rz(-0.66616026) q[3];
sx q[3];
rz(-2.7056138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0824288) q[2];
sx q[2];
rz(-2.6291206) q[2];
sx q[2];
rz(1.814369) q[2];
rz(1.0446769) q[3];
sx q[3];
rz(-2.4278214) q[3];
sx q[3];
rz(2.7248342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5752207) q[0];
sx q[0];
rz(-0.91082585) q[0];
sx q[0];
rz(2.7161993) q[0];
rz(-1.7644024) q[1];
sx q[1];
rz(-1.6433989) q[1];
sx q[1];
rz(-1.4345217) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1578428) q[0];
sx q[0];
rz(-2.2089777) q[0];
sx q[0];
rz(-1.7179836) q[0];
x q[1];
rz(0.70830958) q[2];
sx q[2];
rz(-1.596144) q[2];
sx q[2];
rz(3.1299394) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.19615281) q[1];
sx q[1];
rz(-1.9827794) q[1];
sx q[1];
rz(0.71228551) q[1];
x q[2];
rz(2.1892912) q[3];
sx q[3];
rz(-1.3188014) q[3];
sx q[3];
rz(-2.9661055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7003358) q[2];
sx q[2];
rz(-1.2563027) q[2];
sx q[2];
rz(-1.2909935) q[2];
rz(1.7839606) q[3];
sx q[3];
rz(-1.6358401) q[3];
sx q[3];
rz(1.6443845) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5212379) q[0];
sx q[0];
rz(-2.1339895) q[0];
sx q[0];
rz(-2.3714016) q[0];
rz(-0.99984804) q[1];
sx q[1];
rz(-2.5412173) q[1];
sx q[1];
rz(-1.3410478) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1091546) q[0];
sx q[0];
rz(-1.0312092) q[0];
sx q[0];
rz(-1.3225609) q[0];
rz(-pi) q[1];
rz(-2.6821939) q[2];
sx q[2];
rz(-0.93743491) q[2];
sx q[2];
rz(-0.49672302) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.25136687) q[1];
sx q[1];
rz(-2.736675) q[1];
sx q[1];
rz(0.7271073) q[1];
x q[2];
rz(3.0422736) q[3];
sx q[3];
rz(-1.2867974) q[3];
sx q[3];
rz(-1.8530396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3406713) q[2];
sx q[2];
rz(-2.0719216) q[2];
sx q[2];
rz(0.7684024) q[2];
rz(-0.10041222) q[3];
sx q[3];
rz(-1.6008335) q[3];
sx q[3];
rz(1.2787308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95555821) q[0];
sx q[0];
rz(-2.8362507) q[0];
sx q[0];
rz(0.22698639) q[0];
rz(1.7598033) q[1];
sx q[1];
rz(-0.58265668) q[1];
sx q[1];
rz(1.7452128) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67302746) q[0];
sx q[0];
rz(-1.2224406) q[0];
sx q[0];
rz(0.38469108) q[0];
rz(-pi) q[1];
rz(-0.41047217) q[2];
sx q[2];
rz(-1.7232401) q[2];
sx q[2];
rz(1.4727915) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0048843) q[1];
sx q[1];
rz(-1.6423823) q[1];
sx q[1];
rz(-2.5077255) q[1];
x q[2];
rz(-1.5377858) q[3];
sx q[3];
rz(-2.3950658) q[3];
sx q[3];
rz(-1.3804264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.22623006) q[2];
sx q[2];
rz(-1.9910944) q[2];
sx q[2];
rz(0.076233141) q[2];
rz(-1.3876312) q[3];
sx q[3];
rz(-2.1662655) q[3];
sx q[3];
rz(1.7414198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(0.48080322) q[0];
sx q[0];
rz(-1.5605518) q[0];
sx q[0];
rz(1.3355108) q[0];
rz(-2.238359) q[1];
sx q[1];
rz(-1.3390373) q[1];
sx q[1];
rz(2.9551771) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53509287) q[0];
sx q[0];
rz(-2.4267174) q[0];
sx q[0];
rz(2.2660648) q[0];
rz(-pi) q[1];
rz(-1.769563) q[2];
sx q[2];
rz(-1.0093401) q[2];
sx q[2];
rz(1.1793062) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.32115667) q[1];
sx q[1];
rz(-1.3296491) q[1];
sx q[1];
rz(-2.2187869) q[1];
rz(-pi) q[2];
rz(0.49390467) q[3];
sx q[3];
rz(-2.0155689) q[3];
sx q[3];
rz(-2.1879823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4621801) q[2];
sx q[2];
rz(-1.9483515) q[2];
sx q[2];
rz(-1.818044) q[2];
rz(1.1421674) q[3];
sx q[3];
rz(-0.78502941) q[3];
sx q[3];
rz(-0.89046684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46397504) q[0];
sx q[0];
rz(-2.9587726) q[0];
sx q[0];
rz(-0.63419813) q[0];
rz(-2.0967261) q[1];
sx q[1];
rz(-2.2506782) q[1];
sx q[1];
rz(-0.96010906) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7564108) q[0];
sx q[0];
rz(-0.39968458) q[0];
sx q[0];
rz(1.5002285) q[0];
rz(-pi) q[1];
x q[1];
rz(0.025770806) q[2];
sx q[2];
rz(-0.30116815) q[2];
sx q[2];
rz(1.1045375) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3893435) q[1];
sx q[1];
rz(-1.5592119) q[1];
sx q[1];
rz(2.0962135) q[1];
rz(-2.3784901) q[3];
sx q[3];
rz(-1.0013442) q[3];
sx q[3];
rz(-2.4214937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1977957) q[2];
sx q[2];
rz(-0.65239492) q[2];
sx q[2];
rz(-1.5459527) q[2];
rz(1.4173077) q[3];
sx q[3];
rz(-0.82481074) q[3];
sx q[3];
rz(-1.5203169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5328131) q[0];
sx q[0];
rz(-2.9488035) q[0];
sx q[0];
rz(0.97484318) q[0];
rz(3.0335562) q[1];
sx q[1];
rz(-1.255475) q[1];
sx q[1];
rz(-1.1688165) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10442142) q[0];
sx q[0];
rz(-1.9871431) q[0];
sx q[0];
rz(-2.092931) q[0];
rz(-1.4270876e-05) q[2];
sx q[2];
rz(-0.76185267) q[2];
sx q[2];
rz(0.95214168) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.16202422) q[1];
sx q[1];
rz(-1.8037533) q[1];
sx q[1];
rz(0.090090171) q[1];
x q[2];
rz(1.5511369) q[3];
sx q[3];
rz(-2.0624614) q[3];
sx q[3];
rz(-1.2439072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.10890266) q[2];
sx q[2];
rz(-0.87778512) q[2];
sx q[2];
rz(-2.4857793) q[2];
rz(-3.0374895) q[3];
sx q[3];
rz(-1.3452353) q[3];
sx q[3];
rz(0.18812215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26466894) q[0];
sx q[0];
rz(-1.8380565) q[0];
sx q[0];
rz(1.4917829) q[0];
rz(-1.0393633) q[1];
sx q[1];
rz(-0.84016687) q[1];
sx q[1];
rz(-2.6731491) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2830848) q[0];
sx q[0];
rz(-1.383184) q[0];
sx q[0];
rz(1.5751189) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8554815) q[2];
sx q[2];
rz(-1.1435025) q[2];
sx q[2];
rz(-0.22186771) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5891582) q[1];
sx q[1];
rz(-2.3367408) q[1];
sx q[1];
rz(0.22487747) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2896721) q[3];
sx q[3];
rz(-1.2996965) q[3];
sx q[3];
rz(-2.2106314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.09482065) q[2];
sx q[2];
rz(-1.580661) q[2];
sx q[2];
rz(-0.92791933) q[2];
rz(-1.3377442) q[3];
sx q[3];
rz(-0.51023054) q[3];
sx q[3];
rz(-1.4891589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0874262) q[0];
sx q[0];
rz(-2.1091643) q[0];
sx q[0];
rz(-1.077865) q[0];
rz(-2.7742591) q[1];
sx q[1];
rz(-1.8108188) q[1];
sx q[1];
rz(-1.2841388) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3340942) q[0];
sx q[0];
rz(-1.1002514) q[0];
sx q[0];
rz(-0.28326359) q[0];
rz(-pi) q[1];
x q[1];
rz(0.11576368) q[2];
sx q[2];
rz(-2.2211005) q[2];
sx q[2];
rz(-2.8011326) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.974677) q[1];
sx q[1];
rz(-1.484718) q[1];
sx q[1];
rz(-0.15847107) q[1];
rz(-pi) q[2];
rz(-0.23350164) q[3];
sx q[3];
rz(-2.3931008) q[3];
sx q[3];
rz(-0.12602636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1401356) q[2];
sx q[2];
rz(-2.6900901) q[2];
sx q[2];
rz(-0.4293116) q[2];
rz(-2.9863827) q[3];
sx q[3];
rz(-2.8838172) q[3];
sx q[3];
rz(-1.3252873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24406381) q[0];
sx q[0];
rz(-1.5499935) q[0];
sx q[0];
rz(1.5503379) q[0];
rz(-2.2776729) q[1];
sx q[1];
rz(-2.7680631) q[1];
sx q[1];
rz(1.6815129) q[1];
rz(2.7511394) q[2];
sx q[2];
rz(-0.57885546) q[2];
sx q[2];
rz(-0.59138966) q[2];
rz(-0.33959099) q[3];
sx q[3];
rz(-0.97151269) q[3];
sx q[3];
rz(1.1480939) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
