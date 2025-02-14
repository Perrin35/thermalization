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
rz(1.3833157) q[0];
sx q[0];
rz(-1.436469) q[0];
sx q[0];
rz(-2.1767148) q[0];
rz(2.3896253) q[1];
sx q[1];
rz(-2.7120092) q[1];
sx q[1];
rz(1.0493976) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2222264) q[0];
sx q[0];
rz(-2.1458186) q[0];
sx q[0];
rz(-2.7336043) q[0];
rz(-pi) q[1];
rz(-1.1683091) q[2];
sx q[2];
rz(-2.7912729) q[2];
sx q[2];
rz(-2.868319) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.2856372) q[1];
sx q[1];
rz(-2.1826943) q[1];
sx q[1];
rz(-2.6735071) q[1];
rz(-pi) q[2];
rz(0.5792867) q[3];
sx q[3];
rz(-1.4930973) q[3];
sx q[3];
rz(-2.169211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.45399484) q[2];
sx q[2];
rz(-1.5291841) q[2];
sx q[2];
rz(0.018608658) q[2];
rz(-2.5971557) q[3];
sx q[3];
rz(-2.8114522) q[3];
sx q[3];
rz(-0.64793599) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3675156) q[0];
sx q[0];
rz(-1.6870966) q[0];
sx q[0];
rz(2.6065705) q[0];
rz(-2.6016443) q[1];
sx q[1];
rz(-2.5060563) q[1];
sx q[1];
rz(-1.7417057) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0468804) q[0];
sx q[0];
rz(-2.0167354) q[0];
sx q[0];
rz(0.42973862) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1715047) q[2];
sx q[2];
rz(-2.0510489) q[2];
sx q[2];
rz(2.9532972) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.78838879) q[1];
sx q[1];
rz(-1.8192768) q[1];
sx q[1];
rz(-0.39239359) q[1];
rz(2.1407645) q[3];
sx q[3];
rz(-1.0335575) q[3];
sx q[3];
rz(1.1972103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0672368) q[2];
sx q[2];
rz(-1.8017733) q[2];
sx q[2];
rz(-0.92607099) q[2];
rz(2.1615084) q[3];
sx q[3];
rz(-2.1772549) q[3];
sx q[3];
rz(-2.4721691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3493018) q[0];
sx q[0];
rz(-1.2783569) q[0];
sx q[0];
rz(-1.5128304) q[0];
rz(2.8577562) q[1];
sx q[1];
rz(-0.92548871) q[1];
sx q[1];
rz(1.1121174) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84217269) q[0];
sx q[0];
rz(-1.5779691) q[0];
sx q[0];
rz(3.1345063) q[0];
rz(-pi) q[1];
rz(1.4623423) q[2];
sx q[2];
rz(-1.6466093) q[2];
sx q[2];
rz(0.52559847) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2123472) q[1];
sx q[1];
rz(-2.622859) q[1];
sx q[1];
rz(-0.38350819) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1792382) q[3];
sx q[3];
rz(-2.7076004) q[3];
sx q[3];
rz(2.2631462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5230368) q[2];
sx q[2];
rz(-1.7967537) q[2];
sx q[2];
rz(1.1019361) q[2];
rz(-2.2526422) q[3];
sx q[3];
rz(-2.6933653) q[3];
sx q[3];
rz(-2.2811208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76982826) q[0];
sx q[0];
rz(-2.9673321) q[0];
sx q[0];
rz(-2.4523822) q[0];
rz(-0.31461942) q[1];
sx q[1];
rz(-0.92051053) q[1];
sx q[1];
rz(-1.4124195) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61050451) q[0];
sx q[0];
rz(-2.304739) q[0];
sx q[0];
rz(2.713301) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4749522) q[2];
sx q[2];
rz(-1.7828701) q[2];
sx q[2];
rz(2.1532358) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8103579) q[1];
sx q[1];
rz(-2.7572933) q[1];
sx q[1];
rz(-0.15366252) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2127498) q[3];
sx q[3];
rz(-2.6302377) q[3];
sx q[3];
rz(2.8062888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.41669258) q[2];
sx q[2];
rz(-2.576454) q[2];
sx q[2];
rz(-0.55595428) q[2];
rz(1.8703095) q[3];
sx q[3];
rz(-1.8794182) q[3];
sx q[3];
rz(0.2909734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81412643) q[0];
sx q[0];
rz(-3.0004582) q[0];
sx q[0];
rz(-0.59447527) q[0];
rz(-2.5796083) q[1];
sx q[1];
rz(-0.85834208) q[1];
sx q[1];
rz(-0.97602239) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41244477) q[0];
sx q[0];
rz(-1.9322104) q[0];
sx q[0];
rz(1.1114208) q[0];
rz(-1.541829) q[2];
sx q[2];
rz(-2.3246362) q[2];
sx q[2];
rz(0.08836937) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4415293) q[1];
sx q[1];
rz(-0.91845998) q[1];
sx q[1];
rz(-1.9487914) q[1];
rz(-1.3537977) q[3];
sx q[3];
rz(-1.2760538) q[3];
sx q[3];
rz(-0.47302305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.48656616) q[2];
sx q[2];
rz(-1.2532633) q[2];
sx q[2];
rz(-0.61069926) q[2];
rz(0.30580172) q[3];
sx q[3];
rz(-0.96547258) q[3];
sx q[3];
rz(1.3293728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82764757) q[0];
sx q[0];
rz(-3.1231472) q[0];
sx q[0];
rz(1.6754643) q[0];
rz(2.3620391) q[1];
sx q[1];
rz(-1.9773217) q[1];
sx q[1];
rz(-2.4028042) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2125428) q[0];
sx q[0];
rz(-1.3669786) q[0];
sx q[0];
rz(-1.2840367) q[0];
rz(-pi) q[1];
rz(0.30679484) q[2];
sx q[2];
rz(-1.5606797) q[2];
sx q[2];
rz(-2.1515357) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.14911095) q[1];
sx q[1];
rz(-0.68435366) q[1];
sx q[1];
rz(-2.1696287) q[1];
x q[2];
rz(0.36130623) q[3];
sx q[3];
rz(-1.2194467) q[3];
sx q[3];
rz(-1.6811937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.61591992) q[2];
sx q[2];
rz(-1.2391261) q[2];
sx q[2];
rz(-0.8626779) q[2];
rz(0.45977965) q[3];
sx q[3];
rz(-1.516187) q[3];
sx q[3];
rz(-1.9988683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.90726844) q[0];
sx q[0];
rz(-0.37721226) q[0];
sx q[0];
rz(1.020485) q[0];
rz(-0.057295784) q[1];
sx q[1];
rz(-1.4833996) q[1];
sx q[1];
rz(-0.78757706) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60815281) q[0];
sx q[0];
rz(-2.449547) q[0];
sx q[0];
rz(0.37779053) q[0];
rz(-1.4992981) q[2];
sx q[2];
rz(-1.8282969) q[2];
sx q[2];
rz(-1.471164) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5912093) q[1];
sx q[1];
rz(-1.0799284) q[1];
sx q[1];
rz(0.20259133) q[1];
rz(1.5191139) q[3];
sx q[3];
rz(-0.58850901) q[3];
sx q[3];
rz(-1.9543598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3202177) q[2];
sx q[2];
rz(-1.8193974) q[2];
sx q[2];
rz(0.58464948) q[2];
rz(0.65230495) q[3];
sx q[3];
rz(-1.2290596) q[3];
sx q[3];
rz(1.8023796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.898191) q[0];
sx q[0];
rz(-1.8011872) q[0];
sx q[0];
rz(2.8133494) q[0];
rz(-0.39168656) q[1];
sx q[1];
rz(-1.990254) q[1];
sx q[1];
rz(1.4788871) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23173702) q[0];
sx q[0];
rz(-1.8340936) q[0];
sx q[0];
rz(1.80491) q[0];
x q[1];
rz(-2.6172892) q[2];
sx q[2];
rz(-1.1534165) q[2];
sx q[2];
rz(-2.5427172) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7110243) q[1];
sx q[1];
rz(-1.4757753) q[1];
sx q[1];
rz(1.8501758) q[1];
x q[2];
rz(1.9401266) q[3];
sx q[3];
rz(-1.7672667) q[3];
sx q[3];
rz(-2.7854491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.11999232) q[2];
sx q[2];
rz(-1.3730717) q[2];
sx q[2];
rz(0.78869406) q[2];
rz(-2.9005519) q[3];
sx q[3];
rz(-0.92625109) q[3];
sx q[3];
rz(-2.327976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64860827) q[0];
sx q[0];
rz(-2.6032175) q[0];
sx q[0];
rz(0.96187821) q[0];
rz(2.9073763) q[1];
sx q[1];
rz(-2.0030256) q[1];
sx q[1];
rz(-0.73807565) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2797425) q[0];
sx q[0];
rz(-1.8075917) q[0];
sx q[0];
rz(-0.40646942) q[0];
x q[1];
rz(-0.63491027) q[2];
sx q[2];
rz(-1.2037841) q[2];
sx q[2];
rz(0.5762595) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9752561) q[1];
sx q[1];
rz(-1.8615906) q[1];
sx q[1];
rz(2.0086121) q[1];
x q[2];
rz(-1.8515737) q[3];
sx q[3];
rz(-2.2178136) q[3];
sx q[3];
rz(3.0386277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.48478475) q[2];
sx q[2];
rz(-1.022555) q[2];
sx q[2];
rz(-1.2602932) q[2];
rz(1.3336746) q[3];
sx q[3];
rz(-1.0013872) q[3];
sx q[3];
rz(0.26237747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8986847) q[0];
sx q[0];
rz(-2.3299291) q[0];
sx q[0];
rz(-2.2743478) q[0];
rz(-1.0137089) q[1];
sx q[1];
rz(-1.3288682) q[1];
sx q[1];
rz(-1.0702466) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7263171) q[0];
sx q[0];
rz(-1.4196102) q[0];
sx q[0];
rz(-0.12355208) q[0];
x q[1];
rz(0.81999166) q[2];
sx q[2];
rz(-0.89897663) q[2];
sx q[2];
rz(1.9209678) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6079191) q[1];
sx q[1];
rz(-1.1864788) q[1];
sx q[1];
rz(0.55748765) q[1];
rz(-pi) q[2];
rz(0.66859122) q[3];
sx q[3];
rz(-0.1771268) q[3];
sx q[3];
rz(1.9884381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9762207) q[2];
sx q[2];
rz(-0.82087159) q[2];
sx q[2];
rz(1.2316068) q[2];
rz(1.5008789) q[3];
sx q[3];
rz(-0.21017635) q[3];
sx q[3];
rz(-2.4098082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83723849) q[0];
sx q[0];
rz(-1.1501034) q[0];
sx q[0];
rz(1.3788086) q[0];
rz(0.41481836) q[1];
sx q[1];
rz(-1.7638313) q[1];
sx q[1];
rz(1.5324963) q[1];
rz(-1.6258705) q[2];
sx q[2];
rz(-2.016042) q[2];
sx q[2];
rz(2.945937) q[2];
rz(2.1555156) q[3];
sx q[3];
rz(-2.232983) q[3];
sx q[3];
rz(1.9017526) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
