OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6150317) q[0];
sx q[0];
rz(-0.57305133) q[0];
sx q[0];
rz(0.84258643) q[0];
rz(2.1057582) q[1];
sx q[1];
rz(-1.0993212) q[1];
sx q[1];
rz(-1.6834747) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6095088) q[0];
sx q[0];
rz(-2.3713787) q[0];
sx q[0];
rz(-0.30755933) q[0];
rz(-pi) q[1];
rz(-0.61383944) q[2];
sx q[2];
rz(-1.5868574) q[2];
sx q[2];
rz(-1.2889372) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2394489) q[1];
sx q[1];
rz(-1.4414756) q[1];
sx q[1];
rz(0.37954482) q[1];
x q[2];
rz(-1.4487212) q[3];
sx q[3];
rz(-0.97994084) q[3];
sx q[3];
rz(-1.5199682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.52790102) q[2];
sx q[2];
rz(-1.0062904) q[2];
sx q[2];
rz(-2.9620985) q[2];
rz(-1.2256631) q[3];
sx q[3];
rz(-1.7951199) q[3];
sx q[3];
rz(2.3195482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3935788) q[0];
sx q[0];
rz(-2.2606235) q[0];
sx q[0];
rz(2.8161312) q[0];
rz(-1.7851967) q[1];
sx q[1];
rz(-2.0929095) q[1];
sx q[1];
rz(-1.9869841) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5886473) q[0];
sx q[0];
rz(-1.5536904) q[0];
sx q[0];
rz(-0.018789142) q[0];
rz(2.7484659) q[2];
sx q[2];
rz(-2.1596585) q[2];
sx q[2];
rz(-1.7413505) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.11101152) q[1];
sx q[1];
rz(-1.6554553) q[1];
sx q[1];
rz(-0.76534033) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96315893) q[3];
sx q[3];
rz(-0.31664407) q[3];
sx q[3];
rz(2.7505927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4521728) q[2];
sx q[2];
rz(-1.8916811) q[2];
sx q[2];
rz(0.88341218) q[2];
rz(-2.6702821) q[3];
sx q[3];
rz(-1.703197) q[3];
sx q[3];
rz(0.78770351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.31323355) q[0];
sx q[0];
rz(-1.6468843) q[0];
sx q[0];
rz(1.5154243) q[0];
rz(-2.5405163) q[1];
sx q[1];
rz(-0.54769146) q[1];
sx q[1];
rz(2.0498958) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1900345) q[0];
sx q[0];
rz(-1.0150195) q[0];
sx q[0];
rz(2.4843198) q[0];
x q[1];
rz(-2.7269084) q[2];
sx q[2];
rz(-2.1846002) q[2];
sx q[2];
rz(0.85418044) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6669238) q[1];
sx q[1];
rz(-0.83819929) q[1];
sx q[1];
rz(-1.0679507) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6352429) q[3];
sx q[3];
rz(-1.6985053) q[3];
sx q[3];
rz(0.32999048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.320257) q[2];
sx q[2];
rz(-2.6358423) q[2];
sx q[2];
rz(2.2606405) q[2];
rz(1.3736003) q[3];
sx q[3];
rz(-1.6146086) q[3];
sx q[3];
rz(-1.0176456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83051935) q[0];
sx q[0];
rz(-1.7493462) q[0];
sx q[0];
rz(-2.7048892) q[0];
rz(-2.9084335) q[1];
sx q[1];
rz(-1.2522839) q[1];
sx q[1];
rz(-0.31035796) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.045517) q[0];
sx q[0];
rz(-1.1928416) q[0];
sx q[0];
rz(1.4287352) q[0];
rz(-0.15375806) q[2];
sx q[2];
rz(-2.4506844) q[2];
sx q[2];
rz(1.5916057) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1102317) q[1];
sx q[1];
rz(-0.35846113) q[1];
sx q[1];
rz(1.5178174) q[1];
x q[2];
rz(-0.035590812) q[3];
sx q[3];
rz(-1.6944052) q[3];
sx q[3];
rz(1.3328758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.13005304) q[2];
sx q[2];
rz(-2.4235642) q[2];
sx q[2];
rz(-2.0641573) q[2];
rz(-0.056190101) q[3];
sx q[3];
rz(-2.5037933) q[3];
sx q[3];
rz(1.594054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9524277) q[0];
sx q[0];
rz(-1.0452894) q[0];
sx q[0];
rz(2.8919343) q[0];
rz(1.5769618) q[1];
sx q[1];
rz(-0.77762929) q[1];
sx q[1];
rz(-2.2713984) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91709671) q[0];
sx q[0];
rz(-1.3670237) q[0];
sx q[0];
rz(-2.8208371) q[0];
rz(-pi) q[1];
rz(1.382803) q[2];
sx q[2];
rz(-0.88289875) q[2];
sx q[2];
rz(0.57304136) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.89228499) q[1];
sx q[1];
rz(-1.2586437) q[1];
sx q[1];
rz(2.409163) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9719475) q[3];
sx q[3];
rz(-2.1262453) q[3];
sx q[3];
rz(0.58562216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.27328086) q[2];
sx q[2];
rz(-1.3262649) q[2];
sx q[2];
rz(-2.4678521) q[2];
rz(2.8379748) q[3];
sx q[3];
rz(-1.2250591) q[3];
sx q[3];
rz(-1.3195066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34981397) q[0];
sx q[0];
rz(-0.93739167) q[0];
sx q[0];
rz(-0.2579903) q[0];
rz(-0.42516431) q[1];
sx q[1];
rz(-2.185967) q[1];
sx q[1];
rz(-1.649883) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0817889) q[0];
sx q[0];
rz(-0.44133082) q[0];
sx q[0];
rz(-2.9102737) q[0];
x q[1];
rz(-2.4620373) q[2];
sx q[2];
rz(-1.2365885) q[2];
sx q[2];
rz(0.92781767) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7266453) q[1];
sx q[1];
rz(-0.84268314) q[1];
sx q[1];
rz(0.77264087) q[1];
rz(3.0400279) q[3];
sx q[3];
rz(-1.933681) q[3];
sx q[3];
rz(-0.18728072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1288746) q[2];
sx q[2];
rz(-2.1779163) q[2];
sx q[2];
rz(0.091726124) q[2];
rz(2.2979459) q[3];
sx q[3];
rz(-2.1618312) q[3];
sx q[3];
rz(0.89404026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3180852) q[0];
sx q[0];
rz(-1.2473236) q[0];
sx q[0];
rz(-2.7303625) q[0];
rz(2.2757018) q[1];
sx q[1];
rz(-0.31232467) q[1];
sx q[1];
rz(0.033989865) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2720374) q[0];
sx q[0];
rz(-0.79622686) q[0];
sx q[0];
rz(1.7982593) q[0];
x q[1];
rz(2.2529644) q[2];
sx q[2];
rz(-0.76105984) q[2];
sx q[2];
rz(1.467848) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6778292) q[1];
sx q[1];
rz(-0.67968183) q[1];
sx q[1];
rz(1.6512647) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4230698) q[3];
sx q[3];
rz(-2.1805694) q[3];
sx q[3];
rz(2.8019398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6035446) q[2];
sx q[2];
rz(-0.59843439) q[2];
sx q[2];
rz(-2.2650488) q[2];
rz(0.34902469) q[3];
sx q[3];
rz(-1.2004431) q[3];
sx q[3];
rz(0.14311895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8975163) q[0];
sx q[0];
rz(-1.4325457) q[0];
sx q[0];
rz(-2.7476655) q[0];
rz(0.36755964) q[1];
sx q[1];
rz(-1.7575248) q[1];
sx q[1];
rz(1.6961018) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96955339) q[0];
sx q[0];
rz(-1.5759828) q[0];
sx q[0];
rz(1.1374377) q[0];
rz(-pi) q[1];
x q[1];
rz(0.58480279) q[2];
sx q[2];
rz(-1.0324761) q[2];
sx q[2];
rz(0.98758299) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3280914) q[1];
sx q[1];
rz(-0.5760759) q[1];
sx q[1];
rz(-2.8495795) q[1];
rz(-pi) q[2];
x q[2];
rz(0.94533841) q[3];
sx q[3];
rz(-2.1350386) q[3];
sx q[3];
rz(-0.64627796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8470856) q[2];
sx q[2];
rz(-0.89035788) q[2];
sx q[2];
rz(2.7344446) q[2];
rz(1.6242705) q[3];
sx q[3];
rz(-1.1573236) q[3];
sx q[3];
rz(-0.24967641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80609926) q[0];
sx q[0];
rz(-2.6265916) q[0];
sx q[0];
rz(1.8898213) q[0];
rz(-2.4720526) q[1];
sx q[1];
rz(-1.957683) q[1];
sx q[1];
rz(-0.30977419) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4315902) q[0];
sx q[0];
rz(-0.54034034) q[0];
sx q[0];
rz(-2.8938328) q[0];
x q[1];
rz(1.6236213) q[2];
sx q[2];
rz(-2.9644358) q[2];
sx q[2];
rz(-2.0668541) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8621091) q[1];
sx q[1];
rz(-2.9276507) q[1];
sx q[1];
rz(-1.2941542) q[1];
rz(-1.314332) q[3];
sx q[3];
rz(-2.1520352) q[3];
sx q[3];
rz(-0.47282156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4391675) q[2];
sx q[2];
rz(-0.71321407) q[2];
sx q[2];
rz(1.9343728) q[2];
rz(-2.1045945) q[3];
sx q[3];
rz(-1.2456649) q[3];
sx q[3];
rz(0.65565482) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0913775) q[0];
sx q[0];
rz(-1.3239048) q[0];
sx q[0];
rz(1.2058831) q[0];
rz(-2.5559015) q[1];
sx q[1];
rz(-2.060545) q[1];
sx q[1];
rz(-1.6419798) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3033894) q[0];
sx q[0];
rz(-1.2801542) q[0];
sx q[0];
rz(-0.10586664) q[0];
rz(-pi) q[1];
rz(1.6386547) q[2];
sx q[2];
rz(-0.80791622) q[2];
sx q[2];
rz(2.9615336) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3824532) q[1];
sx q[1];
rz(-2.0331953) q[1];
sx q[1];
rz(2.9049302) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1390926) q[3];
sx q[3];
rz(-1.7435939) q[3];
sx q[3];
rz(0.97027422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2575834) q[2];
sx q[2];
rz(-1.7929701) q[2];
sx q[2];
rz(-2.1949027) q[2];
rz(2.7729014) q[3];
sx q[3];
rz(-1.5654516) q[3];
sx q[3];
rz(-0.45599109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7832227) q[0];
sx q[0];
rz(-1.1483648) q[0];
sx q[0];
rz(-0.4233465) q[0];
rz(-3.070667) q[1];
sx q[1];
rz(-1.6880886) q[1];
sx q[1];
rz(-0.26500519) q[1];
rz(1.0735738) q[2];
sx q[2];
rz(-1.1560658) q[2];
sx q[2];
rz(1.8552468) q[2];
rz(1.0740888) q[3];
sx q[3];
rz(-2.2278193) q[3];
sx q[3];
rz(-0.57886119) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
