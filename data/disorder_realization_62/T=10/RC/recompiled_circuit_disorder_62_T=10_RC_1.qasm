OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.9010889) q[0];
sx q[0];
rz(-0.17740372) q[0];
sx q[0];
rz(-2.0071964) q[0];
rz(1.1881243) q[1];
sx q[1];
rz(-2.1048574) q[1];
sx q[1];
rz(-0.66361767) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7952607) q[0];
sx q[0];
rz(-1.5471317) q[0];
sx q[0];
rz(-0.46121009) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.3590091) q[2];
sx q[2];
rz(-0.81072545) q[2];
sx q[2];
rz(0.63149482) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2156148) q[1];
sx q[1];
rz(-1.3091334) q[1];
sx q[1];
rz(1.91933) q[1];
rz(1.8144242) q[3];
sx q[3];
rz(-1.4654136) q[3];
sx q[3];
rz(-1.6107314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2628281) q[2];
sx q[2];
rz(-0.44439134) q[2];
sx q[2];
rz(3.0905241) q[2];
rz(-0.55705327) q[3];
sx q[3];
rz(-2.3414108) q[3];
sx q[3];
rz(-1.5867656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58650815) q[0];
sx q[0];
rz(-0.83207911) q[0];
sx q[0];
rz(2.5449975) q[0];
rz(-0.82582981) q[1];
sx q[1];
rz(-1.4412216) q[1];
sx q[1];
rz(1.9155496) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2873043) q[0];
sx q[0];
rz(-2.0137631) q[0];
sx q[0];
rz(-1.7914377) q[0];
x q[1];
rz(-2.3184899) q[2];
sx q[2];
rz(-2.1777993) q[2];
sx q[2];
rz(2.5603103) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.1204685) q[1];
sx q[1];
rz(-1.7906584) q[1];
sx q[1];
rz(-3.1411509) q[1];
rz(0.0094718178) q[3];
sx q[3];
rz(-0.99820271) q[3];
sx q[3];
rz(0.026281683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3423959) q[2];
sx q[2];
rz(-1.1725972) q[2];
sx q[2];
rz(-0.33102316) q[2];
rz(2.3349169) q[3];
sx q[3];
rz(-1.9162063) q[3];
sx q[3];
rz(-1.7139009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1598635) q[0];
sx q[0];
rz(-1.9112497) q[0];
sx q[0];
rz(-1.249041) q[0];
rz(-0.088009134) q[1];
sx q[1];
rz(-1.1227337) q[1];
sx q[1];
rz(1.0294611) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3345966) q[0];
sx q[0];
rz(-2.1801729) q[0];
sx q[0];
rz(-2.8732804) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7715893) q[2];
sx q[2];
rz(-0.68379935) q[2];
sx q[2];
rz(1.5497269) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2468977) q[1];
sx q[1];
rz(-2.9322335) q[1];
sx q[1];
rz(-0.73114242) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.401628) q[3];
sx q[3];
rz(-1.660534) q[3];
sx q[3];
rz(-1.9268074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.007894667) q[2];
sx q[2];
rz(-1.7241314) q[2];
sx q[2];
rz(-2.6339445) q[2];
rz(1.3890022) q[3];
sx q[3];
rz(-2.8116083) q[3];
sx q[3];
rz(-1.9784137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4841109) q[0];
sx q[0];
rz(-0.27814516) q[0];
sx q[0];
rz(-1.5959651) q[0];
rz(-2.0987299) q[1];
sx q[1];
rz(-1.1735801) q[1];
sx q[1];
rz(-1.5159336) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8221178) q[0];
sx q[0];
rz(-1.428077) q[0];
sx q[0];
rz(0.6832173) q[0];
x q[1];
rz(-1.0117202) q[2];
sx q[2];
rz(-1.498073) q[2];
sx q[2];
rz(-0.59414547) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3486466) q[1];
sx q[1];
rz(-0.9325087) q[1];
sx q[1];
rz(-1.2299041) q[1];
x q[2];
rz(1.7901778) q[3];
sx q[3];
rz(-2.1483768) q[3];
sx q[3];
rz(0.95336174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.76379124) q[2];
sx q[2];
rz(-2.2968473) q[2];
sx q[2];
rz(1.8466922) q[2];
rz(3.0002248) q[3];
sx q[3];
rz(-2.5989792) q[3];
sx q[3];
rz(2.1550089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
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
rz(-0.35448733) q[0];
sx q[0];
rz(-1.961345) q[0];
sx q[0];
rz(1.5198583) q[0];
rz(-0.63201085) q[1];
sx q[1];
rz(-2.4193587) q[1];
sx q[1];
rz(0.89486665) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1631854) q[0];
sx q[0];
rz(-2.1911231) q[0];
sx q[0];
rz(-2.083367) q[0];
rz(-pi) q[1];
rz(-0.72563719) q[2];
sx q[2];
rz(-2.0113532) q[2];
sx q[2];
rz(0.63647905) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4836854) q[1];
sx q[1];
rz(-2.2135418) q[1];
sx q[1];
rz(0.55773463) q[1];
rz(-pi) q[2];
rz(-0.28646333) q[3];
sx q[3];
rz(-1.5634007) q[3];
sx q[3];
rz(-1.3163819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.616509) q[2];
sx q[2];
rz(-2.775165) q[2];
sx q[2];
rz(1.1425225) q[2];
rz(2.3948495) q[3];
sx q[3];
rz(-1.8871504) q[3];
sx q[3];
rz(1.07553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.557945) q[0];
sx q[0];
rz(-0.32145158) q[0];
sx q[0];
rz(1.7640132) q[0];
rz(0.47239834) q[1];
sx q[1];
rz(-0.51858416) q[1];
sx q[1];
rz(2.6766052) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12596345) q[0];
sx q[0];
rz(-1.1383801) q[0];
sx q[0];
rz(-2.7727491) q[0];
rz(-pi) q[1];
rz(-1.497252) q[2];
sx q[2];
rz(-1.2300756) q[2];
sx q[2];
rz(1.7433012) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.28593674) q[1];
sx q[1];
rz(-0.93614139) q[1];
sx q[1];
rz(-3.082285) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7808077) q[3];
sx q[3];
rz(-0.87267733) q[3];
sx q[3];
rz(2.7426646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.49089367) q[2];
sx q[2];
rz(-1.2630254) q[2];
sx q[2];
rz(0.4450376) q[2];
rz(-2.2079091) q[3];
sx q[3];
rz(-1.4327587) q[3];
sx q[3];
rz(2.8745108) q[3];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0141107) q[0];
sx q[0];
rz(-1.575379) q[0];
sx q[0];
rz(-1.7215464) q[0];
rz(0.02380112) q[1];
sx q[1];
rz(-2.5271466) q[1];
sx q[1];
rz(-2.9856317) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50990803) q[0];
sx q[0];
rz(-1.8242867) q[0];
sx q[0];
rz(-1.7476728) q[0];
rz(-pi) q[1];
rz(-0.67955741) q[2];
sx q[2];
rz(-2.5153) q[2];
sx q[2];
rz(-1.2223513) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1289039) q[1];
sx q[1];
rz(-1.3959937) q[1];
sx q[1];
rz(-0.054794475) q[1];
rz(2.3994) q[3];
sx q[3];
rz(-2.4626197) q[3];
sx q[3];
rz(1.9051444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0722787) q[2];
sx q[2];
rz(-2.5645655) q[2];
sx q[2];
rz(-2.0689266) q[2];
rz(-0.32564751) q[3];
sx q[3];
rz(-1.1762534) q[3];
sx q[3];
rz(-1.5163039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1919365) q[0];
sx q[0];
rz(-0.40238109) q[0];
sx q[0];
rz(2.8038213) q[0];
rz(-2.0514964) q[1];
sx q[1];
rz(-0.97507674) q[1];
sx q[1];
rz(0.24857323) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2598495) q[0];
sx q[0];
rz(-1.8200841) q[0];
sx q[0];
rz(-0.52815204) q[0];
x q[1];
rz(-0.82069223) q[2];
sx q[2];
rz(-2.4889915) q[2];
sx q[2];
rz(1.5479659) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.12411815) q[1];
sx q[1];
rz(-1.0867456) q[1];
sx q[1];
rz(2.2494621) q[1];
x q[2];
rz(-0.076898889) q[3];
sx q[3];
rz(-1.4014763) q[3];
sx q[3];
rz(0.83991915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8481855) q[2];
sx q[2];
rz(-0.53331393) q[2];
sx q[2];
rz(-0.08671134) q[2];
rz(2.6596206) q[3];
sx q[3];
rz(-0.50589365) q[3];
sx q[3];
rz(-1.988525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50865737) q[0];
sx q[0];
rz(-2.1338699) q[0];
sx q[0];
rz(-2.8588262) q[0];
rz(-0.70156082) q[1];
sx q[1];
rz(-0.82156721) q[1];
sx q[1];
rz(-1.823002) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5134207) q[0];
sx q[0];
rz(-1.4903729) q[0];
sx q[0];
rz(1.1286939) q[0];
rz(-pi) q[1];
rz(0.96005) q[2];
sx q[2];
rz(-1.0078444) q[2];
sx q[2];
rz(-2.8658531) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9642155) q[1];
sx q[1];
rz(-0.69987684) q[1];
sx q[1];
rz(1.3436951) q[1];
rz(2.7301844) q[3];
sx q[3];
rz(-2.3477926) q[3];
sx q[3];
rz(-1.0994764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41137722) q[2];
sx q[2];
rz(-2.8432379) q[2];
sx q[2];
rz(-0.39917699) q[2];
rz(0.88360751) q[3];
sx q[3];
rz(-1.6688321) q[3];
sx q[3];
rz(-1.2214899) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76319641) q[0];
sx q[0];
rz(-2.3662687) q[0];
sx q[0];
rz(-1.5989074) q[0];
rz(-1.0653161) q[1];
sx q[1];
rz(-0.97384802) q[1];
sx q[1];
rz(1.8803966) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48860088) q[0];
sx q[0];
rz(-1.4654667) q[0];
sx q[0];
rz(-2.7306042) q[0];
rz(-pi) q[1];
x q[1];
rz(3.10896) q[2];
sx q[2];
rz(-2.543078) q[2];
sx q[2];
rz(-2.0805217) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.65834261) q[1];
sx q[1];
rz(-1.634942) q[1];
sx q[1];
rz(2.642753) q[1];
rz(-pi) q[2];
x q[2];
rz(2.946978) q[3];
sx q[3];
rz(-2.448304) q[3];
sx q[3];
rz(-0.64893901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.5579055) q[2];
sx q[2];
rz(-0.62290278) q[2];
sx q[2];
rz(-2.5718001) q[2];
rz(1.2184881) q[3];
sx q[3];
rz(-2.4945538) q[3];
sx q[3];
rz(2.5789554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8284843) q[0];
sx q[0];
rz(-2.2611571) q[0];
sx q[0];
rz(-1.8631998) q[0];
rz(2.5333511) q[1];
sx q[1];
rz(-0.47641644) q[1];
sx q[1];
rz(0.48412916) q[1];
rz(2.1424978) q[2];
sx q[2];
rz(-1.5487164) q[2];
sx q[2];
rz(-1.6185417) q[2];
rz(2.5504997) q[3];
sx q[3];
rz(-1.4553087) q[3];
sx q[3];
rz(1.6607264) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
