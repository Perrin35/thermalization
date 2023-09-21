OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.66184008) q[0];
sx q[0];
rz(2.2979484) q[0];
sx q[0];
rz(9.2568682) q[0];
rz(-1.9703938) q[1];
sx q[1];
rz(-0.29532239) q[1];
sx q[1];
rz(-0.056161031) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0060318) q[0];
sx q[0];
rz(-2.6132085) q[0];
sx q[0];
rz(1.1392659) q[0];
rz(-pi) q[1];
rz(-2.9287455) q[2];
sx q[2];
rz(-2.2058862) q[2];
sx q[2];
rz(-2.0095306) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.72370428) q[1];
sx q[1];
rz(-2.1065518) q[1];
sx q[1];
rz(-0.31236155) q[1];
x q[2];
rz(-0.25986259) q[3];
sx q[3];
rz(-1.6136323) q[3];
sx q[3];
rz(-2.6916137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.37796676) q[2];
sx q[2];
rz(-0.28181919) q[2];
sx q[2];
rz(0.4326694) q[2];
rz(-1.1928605) q[3];
sx q[3];
rz(-1.9038707) q[3];
sx q[3];
rz(0.38309923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6137961) q[0];
sx q[0];
rz(-0.48848099) q[0];
sx q[0];
rz(1.312785) q[0];
rz(0.20547543) q[1];
sx q[1];
rz(-0.97646362) q[1];
sx q[1];
rz(1.9899433) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8173556) q[0];
sx q[0];
rz(-1.2756057) q[0];
sx q[0];
rz(-0.8582219) q[0];
rz(-pi) q[1];
rz(-2.5580514) q[2];
sx q[2];
rz(-1.1569287) q[2];
sx q[2];
rz(-1.5838503) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0484867) q[1];
sx q[1];
rz(-1.152532) q[1];
sx q[1];
rz(-2.6622245) q[1];
rz(-1.9563975) q[3];
sx q[3];
rz(-2.0413627) q[3];
sx q[3];
rz(3.1009931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.0097222086) q[2];
sx q[2];
rz(-1.6513609) q[2];
sx q[2];
rz(0.22182626) q[2];
rz(0.37718537) q[3];
sx q[3];
rz(-0.42755869) q[3];
sx q[3];
rz(-2.3691573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8310228) q[0];
sx q[0];
rz(-3.0492058) q[0];
sx q[0];
rz(-3.1047399) q[0];
rz(-0.82551461) q[1];
sx q[1];
rz(-1.3157536) q[1];
sx q[1];
rz(-3.085014) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73233561) q[0];
sx q[0];
rz(-1.1261254) q[0];
sx q[0];
rz(2.1952941) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4372196) q[2];
sx q[2];
rz(-1.8250416) q[2];
sx q[2];
rz(2.996252) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0128855) q[1];
sx q[1];
rz(-1.371908) q[1];
sx q[1];
rz(2.835564) q[1];
rz(-2.4967381) q[3];
sx q[3];
rz(-1.2554902) q[3];
sx q[3];
rz(2.9063318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8686707) q[2];
sx q[2];
rz(-1.6438831) q[2];
sx q[2];
rz(-0.92612129) q[2];
rz(2.5849294) q[3];
sx q[3];
rz(-2.848048) q[3];
sx q[3];
rz(2.0986957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2264003) q[0];
sx q[0];
rz(-0.72015786) q[0];
sx q[0];
rz(0.74209374) q[0];
rz(2.0023951) q[1];
sx q[1];
rz(-0.4793872) q[1];
sx q[1];
rz(-0.46359584) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47953654) q[0];
sx q[0];
rz(-2.2502796) q[0];
sx q[0];
rz(-2.7583073) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6323339) q[2];
sx q[2];
rz(-1.9029641) q[2];
sx q[2];
rz(2.5663944) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1222893) q[1];
sx q[1];
rz(-0.79320723) q[1];
sx q[1];
rz(-1.3293468) q[1];
rz(-pi) q[2];
x q[2];
rz(0.87426825) q[3];
sx q[3];
rz(-0.6797176) q[3];
sx q[3];
rz(0.10248871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7745557) q[2];
sx q[2];
rz(-2.98428) q[2];
sx q[2];
rz(3.0920933) q[2];
rz(0.1285304) q[3];
sx q[3];
rz(-1.5934207) q[3];
sx q[3];
rz(0.11894225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13609919) q[0];
sx q[0];
rz(-0.43731421) q[0];
sx q[0];
rz(-2.8438925) q[0];
rz(0.4822576) q[1];
sx q[1];
rz(-0.75459701) q[1];
sx q[1];
rz(-2.1972426) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1221065) q[0];
sx q[0];
rz(-1.5257611) q[0];
sx q[0];
rz(-1.3615863) q[0];
x q[1];
rz(2.7371251) q[2];
sx q[2];
rz(-2.3239115) q[2];
sx q[2];
rz(-2.5596465) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.22630616) q[1];
sx q[1];
rz(-1.5686791) q[1];
sx q[1];
rz(-1.6318984) q[1];
rz(-pi) q[2];
x q[2];
rz(0.43290187) q[3];
sx q[3];
rz(-1.8042943) q[3];
sx q[3];
rz(-1.3732861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8828316) q[2];
sx q[2];
rz(-2.0488887) q[2];
sx q[2];
rz(3.0333701) q[2];
rz(-3.1392858) q[3];
sx q[3];
rz(-1.6121515) q[3];
sx q[3];
rz(0.32430696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(-2.7047983) q[0];
sx q[0];
rz(-2.695485) q[0];
sx q[0];
rz(0.58445245) q[0];
rz(2.2553518) q[1];
sx q[1];
rz(-0.61683547) q[1];
sx q[1];
rz(-3.086673) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7391101) q[0];
sx q[0];
rz(-1.4882898) q[0];
sx q[0];
rz(-1.301618) q[0];
rz(-pi) q[1];
x q[1];
rz(0.018718406) q[2];
sx q[2];
rz(-1.2476377) q[2];
sx q[2];
rz(-1.2013555) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7683148) q[1];
sx q[1];
rz(-1.0265961) q[1];
sx q[1];
rz(-2.0139704) q[1];
x q[2];
rz(-0.92298569) q[3];
sx q[3];
rz(-0.69159782) q[3];
sx q[3];
rz(1.637961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.75446689) q[2];
sx q[2];
rz(-3.0209164) q[2];
sx q[2];
rz(1.0167271) q[2];
rz(2.5975442) q[3];
sx q[3];
rz(-0.35990158) q[3];
sx q[3];
rz(-1.8090766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.465437) q[0];
sx q[0];
rz(-0.98452079) q[0];
sx q[0];
rz(-0.28453919) q[0];
rz(-2.1971205) q[1];
sx q[1];
rz(-1.1962793) q[1];
sx q[1];
rz(-0.91032666) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7401687) q[0];
sx q[0];
rz(-1.6356042) q[0];
sx q[0];
rz(-3.0868953) q[0];
rz(-pi) q[1];
rz(1.1548642) q[2];
sx q[2];
rz(-1.8850733) q[2];
sx q[2];
rz(-1.2736125) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7200304) q[1];
sx q[1];
rz(-1.3898464) q[1];
sx q[1];
rz(-0.49444316) q[1];
rz(-pi) q[2];
rz(-0.63776871) q[3];
sx q[3];
rz(-0.83003269) q[3];
sx q[3];
rz(0.72731599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3900782) q[2];
sx q[2];
rz(-3.0780767) q[2];
sx q[2];
rz(-2.2195623) q[2];
rz(-0.56728029) q[3];
sx q[3];
rz(-1.7276948) q[3];
sx q[3];
rz(2.1218307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2475125) q[0];
sx q[0];
rz(-0.67665726) q[0];
sx q[0];
rz(-3.0122053) q[0];
rz(2.5091876) q[1];
sx q[1];
rz(-1.0267195) q[1];
sx q[1];
rz(0.30050373) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43113118) q[0];
sx q[0];
rz(-0.80695242) q[0];
sx q[0];
rz(-1.0459082) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1883165) q[2];
sx q[2];
rz(-2.2308908) q[2];
sx q[2];
rz(-0.88027871) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.064425163) q[1];
sx q[1];
rz(-2.7556813) q[1];
sx q[1];
rz(0.47972958) q[1];
rz(-pi) q[2];
x q[2];
rz(0.31605966) q[3];
sx q[3];
rz(-2.3028767) q[3];
sx q[3];
rz(-2.4162606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.58632103) q[2];
sx q[2];
rz(-2.1466612) q[2];
sx q[2];
rz(-2.3596181) q[2];
rz(-2.590495) q[3];
sx q[3];
rz(-1.3827773) q[3];
sx q[3];
rz(-2.9836695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5683811) q[0];
sx q[0];
rz(-1.1567572) q[0];
sx q[0];
rz(0.12776275) q[0];
rz(-0.54221517) q[1];
sx q[1];
rz(-2.1844889) q[1];
sx q[1];
rz(2.382747) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7153873) q[0];
sx q[0];
rz(-1.9949159) q[0];
sx q[0];
rz(-3.12294) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7636289) q[2];
sx q[2];
rz(-0.97059965) q[2];
sx q[2];
rz(2.6892975) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5938877) q[1];
sx q[1];
rz(-0.87915671) q[1];
sx q[1];
rz(1.7734852) q[1];
x q[2];
rz(-2.2364053) q[3];
sx q[3];
rz(-1.5593411) q[3];
sx q[3];
rz(-0.36597914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1252497) q[2];
sx q[2];
rz(-1.3602076) q[2];
sx q[2];
rz(2.6515567) q[2];
rz(1.4222493) q[3];
sx q[3];
rz(-1.9336721) q[3];
sx q[3];
rz(-2.0786044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7816417) q[0];
sx q[0];
rz(-2.521823) q[0];
sx q[0];
rz(0.075335659) q[0];
rz(-0.8967337) q[1];
sx q[1];
rz(-1.1743841) q[1];
sx q[1];
rz(-2.5316701) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2595554) q[0];
sx q[0];
rz(-1.8543188) q[0];
sx q[0];
rz(2.3638704) q[0];
x q[1];
rz(1.2953193) q[2];
sx q[2];
rz(-1.2071949) q[2];
sx q[2];
rz(1.3678577) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.48253179) q[1];
sx q[1];
rz(-0.7809124) q[1];
sx q[1];
rz(-0.34300967) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6033953) q[3];
sx q[3];
rz(-2.3231069) q[3];
sx q[3];
rz(1.324211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.909409) q[2];
sx q[2];
rz(-2.3164618) q[2];
sx q[2];
rz(2.4278736) q[2];
rz(-2.7632726) q[3];
sx q[3];
rz(-0.49351966) q[3];
sx q[3];
rz(2.2617214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.338035) q[0];
sx q[0];
rz(-1.150158) q[0];
sx q[0];
rz(-1.5858571) q[0];
rz(-0.65080416) q[1];
sx q[1];
rz(-1.4918068) q[1];
sx q[1];
rz(3.0202958) q[1];
rz(-1.4177633) q[2];
sx q[2];
rz(-0.19822181) q[2];
sx q[2];
rz(-0.76186686) q[2];
rz(-1.1013423) q[3];
sx q[3];
rz(-2.6228842) q[3];
sx q[3];
rz(1.7026671) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];