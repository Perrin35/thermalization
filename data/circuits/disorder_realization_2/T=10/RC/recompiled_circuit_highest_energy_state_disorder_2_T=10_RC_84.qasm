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
rz(-0.18345565) q[0];
sx q[0];
rz(-0.7440716) q[0];
sx q[0];
rz(2.0896572) q[0];
rz(-2.9479041) q[1];
sx q[1];
rz(-2.5084578) q[1];
sx q[1];
rz(2.5685891) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6010701) q[0];
sx q[0];
rz(-2.1946476) q[0];
sx q[0];
rz(0.076657587) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.51197211) q[2];
sx q[2];
rz(-0.7535156) q[2];
sx q[2];
rz(-1.3775502) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6261648) q[1];
sx q[1];
rz(-2.0671131) q[1];
sx q[1];
rz(1.4908916) q[1];
rz(-2.4450847) q[3];
sx q[3];
rz(-1.4918431) q[3];
sx q[3];
rz(-2.1895308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8594592) q[2];
sx q[2];
rz(-2.8318475) q[2];
sx q[2];
rz(2.0852883) q[2];
rz(0.022196444) q[3];
sx q[3];
rz(-2.3789417) q[3];
sx q[3];
rz(1.7215884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9098772) q[0];
sx q[0];
rz(-2.660399) q[0];
sx q[0];
rz(0.43014446) q[0];
rz(-0.12769708) q[1];
sx q[1];
rz(-1.1558497) q[1];
sx q[1];
rz(1.7040303) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24457045) q[0];
sx q[0];
rz(-1.8720227) q[0];
sx q[0];
rz(1.8409077) q[0];
rz(0.050927863) q[2];
sx q[2];
rz(-1.5136592) q[2];
sx q[2];
rz(-2.6661851) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.592652) q[1];
sx q[1];
rz(-2.2363064) q[1];
sx q[1];
rz(-1.9926461) q[1];
rz(-2.7892465) q[3];
sx q[3];
rz(-0.56803136) q[3];
sx q[3];
rz(1.2888792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0251856) q[2];
sx q[2];
rz(-1.4613287) q[2];
sx q[2];
rz(0.44542584) q[2];
rz(0.40924254) q[3];
sx q[3];
rz(-2.0425115) q[3];
sx q[3];
rz(0.87944952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55260783) q[0];
sx q[0];
rz(-3.0486076) q[0];
sx q[0];
rz(2.4267922) q[0];
rz(-1.0572761) q[1];
sx q[1];
rz(-2.7209268) q[1];
sx q[1];
rz(-0.32726273) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3540368) q[0];
sx q[0];
rz(-1.4041413) q[0];
sx q[0];
rz(2.119407) q[0];
rz(1.0017603) q[2];
sx q[2];
rz(-1.941701) q[2];
sx q[2];
rz(0.47874622) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.63673692) q[1];
sx q[1];
rz(-1.550192) q[1];
sx q[1];
rz(1.4192753) q[1];
rz(-pi) q[2];
rz(-0.9483665) q[3];
sx q[3];
rz(-0.41928681) q[3];
sx q[3];
rz(-0.47800999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1383692) q[2];
sx q[2];
rz(-0.97347632) q[2];
sx q[2];
rz(-1.5039911) q[2];
rz(1.2184294) q[3];
sx q[3];
rz(-2.7259493) q[3];
sx q[3];
rz(0.98201069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0686491) q[0];
sx q[0];
rz(-1.2462085) q[0];
sx q[0];
rz(0.15596998) q[0];
rz(0.34863696) q[1];
sx q[1];
rz(-2.5376384) q[1];
sx q[1];
rz(1.9085931) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1240163) q[0];
sx q[0];
rz(-1.6917545) q[0];
sx q[0];
rz(1.8328387) q[0];
rz(0.92278752) q[2];
sx q[2];
rz(-1.7402667) q[2];
sx q[2];
rz(-3.0344935) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.33267) q[1];
sx q[1];
rz(-2.4415209) q[1];
sx q[1];
rz(-2.0168522) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25147049) q[3];
sx q[3];
rz(-1.5850388) q[3];
sx q[3];
rz(-0.79232681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4431346) q[2];
sx q[2];
rz(-0.036018697) q[2];
sx q[2];
rz(-1.6301463) q[2];
rz(-2.002142) q[3];
sx q[3];
rz(-1.3316493) q[3];
sx q[3];
rz(0.99036923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-1.3646506) q[0];
sx q[0];
rz(-2.7237837) q[0];
sx q[0];
rz(1.1530217) q[0];
rz(0.54368377) q[1];
sx q[1];
rz(-1.8749571) q[1];
sx q[1];
rz(-2.8797454) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3735298) q[0];
sx q[0];
rz(-0.090001194) q[0];
sx q[0];
rz(-1.9530208) q[0];
rz(-pi) q[1];
rz(-0.98723094) q[2];
sx q[2];
rz(-1.2886184) q[2];
sx q[2];
rz(-2.5099177) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0386069) q[1];
sx q[1];
rz(-0.9894045) q[1];
sx q[1];
rz(-1.8379148) q[1];
rz(-0.7821726) q[3];
sx q[3];
rz(-1.425079) q[3];
sx q[3];
rz(-2.4268933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5220945) q[2];
sx q[2];
rz(-0.60569373) q[2];
sx q[2];
rz(1.4052793) q[2];
rz(-2.3960579) q[3];
sx q[3];
rz(-1.2657974) q[3];
sx q[3];
rz(2.1982927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.140542) q[0];
sx q[0];
rz(-3.0241835) q[0];
sx q[0];
rz(-2.7897799) q[0];
rz(2.70128) q[1];
sx q[1];
rz(-1.7761209) q[1];
sx q[1];
rz(0.17131677) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7650314) q[0];
sx q[0];
rz(-2.2269571) q[0];
sx q[0];
rz(1.879346) q[0];
x q[1];
rz(-1.8753042) q[2];
sx q[2];
rz(-0.84104462) q[2];
sx q[2];
rz(-0.27308057) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.80444634) q[1];
sx q[1];
rz(-0.69883543) q[1];
sx q[1];
rz(3.1091467) q[1];
rz(-pi) q[2];
rz(-2.4279057) q[3];
sx q[3];
rz(-2.1296576) q[3];
sx q[3];
rz(1.1446468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.064284023) q[2];
sx q[2];
rz(-1.372154) q[2];
sx q[2];
rz(2.1997931) q[2];
rz(-1.9866379) q[3];
sx q[3];
rz(-2.933511) q[3];
sx q[3];
rz(-1.8010767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6963541) q[0];
sx q[0];
rz(-1.1055163) q[0];
sx q[0];
rz(0.0048986991) q[0];
rz(-0.44662961) q[1];
sx q[1];
rz(-0.59372562) q[1];
sx q[1];
rz(-0.68797025) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60058355) q[0];
sx q[0];
rz(-2.2760411) q[0];
sx q[0];
rz(1.7974822) q[0];
rz(-pi) q[1];
rz(-1.2254459) q[2];
sx q[2];
rz(-0.98042578) q[2];
sx q[2];
rz(-1.2638059) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1157284) q[1];
sx q[1];
rz(-1.1394355) q[1];
sx q[1];
rz(-0.23484767) q[1];
rz(-pi) q[2];
rz(0.79598456) q[3];
sx q[3];
rz(-0.86382691) q[3];
sx q[3];
rz(-0.10315264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.61116162) q[2];
sx q[2];
rz(-2.7335584) q[2];
sx q[2];
rz(-0.92998695) q[2];
rz(-1.9200578) q[3];
sx q[3];
rz(-2.0285716) q[3];
sx q[3];
rz(0.59337029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4596443) q[0];
sx q[0];
rz(-1.821803) q[0];
sx q[0];
rz(0.090106877) q[0];
rz(1.7168761) q[1];
sx q[1];
rz(-0.63980353) q[1];
sx q[1];
rz(2.7065014) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6259436) q[0];
sx q[0];
rz(-2.0688147) q[0];
sx q[0];
rz(-1.7367212) q[0];
rz(-2.4736011) q[2];
sx q[2];
rz(-1.7133623) q[2];
sx q[2];
rz(-3.0270456) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7104946) q[1];
sx q[1];
rz(-2.5113547) q[1];
sx q[1];
rz(-1.3517387) q[1];
rz(-pi) q[2];
rz(1.5699205) q[3];
sx q[3];
rz(-1.9951576) q[3];
sx q[3];
rz(-2.3800338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6626176) q[2];
sx q[2];
rz(-1.0650029) q[2];
sx q[2];
rz(-2.5153861) q[2];
rz(-0.19691697) q[3];
sx q[3];
rz(-2.1508689) q[3];
sx q[3];
rz(1.3885434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58186746) q[0];
sx q[0];
rz(-1.4511061) q[0];
sx q[0];
rz(0.054280601) q[0];
rz(0.42298969) q[1];
sx q[1];
rz(-1.3754247) q[1];
sx q[1];
rz(1.8064226) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30935449) q[0];
sx q[0];
rz(-1.0296427) q[0];
sx q[0];
rz(-2.4914129) q[0];
rz(-pi) q[1];
x q[1];
rz(0.35690709) q[2];
sx q[2];
rz(-2.1015328) q[2];
sx q[2];
rz(-0.70009795) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4563417) q[1];
sx q[1];
rz(-2.0944746) q[1];
sx q[1];
rz(-2.1728553) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0306251) q[3];
sx q[3];
rz(-2.9879192) q[3];
sx q[3];
rz(2.3914162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.52428952) q[2];
sx q[2];
rz(-1.6733988) q[2];
sx q[2];
rz(-0.35150251) q[2];
rz(1.7678123) q[3];
sx q[3];
rz(-2.6280554) q[3];
sx q[3];
rz(-2.8721749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3897301) q[0];
sx q[0];
rz(-2.8808012) q[0];
sx q[0];
rz(0.75916284) q[0];
rz(2.0579386) q[1];
sx q[1];
rz(-1.6108797) q[1];
sx q[1];
rz(-1.0677451) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6747492) q[0];
sx q[0];
rz(-2.3567794) q[0];
sx q[0];
rz(-0.87133566) q[0];
x q[1];
rz(1.2933838) q[2];
sx q[2];
rz(-1.8545863) q[2];
sx q[2];
rz(1.0671771) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8651004) q[1];
sx q[1];
rz(-2.4363764) q[1];
sx q[1];
rz(0.92030763) q[1];
x q[2];
rz(-1.178252) q[3];
sx q[3];
rz(-0.7362698) q[3];
sx q[3];
rz(0.53680778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7119673) q[2];
sx q[2];
rz(-1.9316614) q[2];
sx q[2];
rz(2.9313226) q[2];
rz(0.78768864) q[3];
sx q[3];
rz(-1.7588047) q[3];
sx q[3];
rz(-2.6113966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6982211) q[0];
sx q[0];
rz(-2.5826695) q[0];
sx q[0];
rz(-1.2179751) q[0];
rz(-1.5086077) q[1];
sx q[1];
rz(-1.8764381) q[1];
sx q[1];
rz(1.1988342) q[1];
rz(1.0883349) q[2];
sx q[2];
rz(-0.74813048) q[2];
sx q[2];
rz(2.8449423) q[2];
rz(1.0630332) q[3];
sx q[3];
rz(-0.31217839) q[3];
sx q[3];
rz(0.97806539) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
