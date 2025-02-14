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
rz(1.4662161) q[0];
sx q[0];
rz(5.3386547) q[0];
sx q[0];
rz(9.6261779) q[0];
rz(1.072285) q[1];
sx q[1];
rz(-0.76289248) q[1];
sx q[1];
rz(-1.720517) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.759255) q[0];
sx q[0];
rz(-1.6256285) q[0];
sx q[0];
rz(2.1817529) q[0];
rz(-pi) q[1];
rz(-2.1915053) q[2];
sx q[2];
rz(-0.7751152) q[2];
sx q[2];
rz(-0.24573869) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.78064976) q[1];
sx q[1];
rz(-0.93027481) q[1];
sx q[1];
rz(-1.0761797) q[1];
rz(1.9626612) q[3];
sx q[3];
rz(-1.2614935) q[3];
sx q[3];
rz(0.40180692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1066771) q[2];
sx q[2];
rz(-2.3372529) q[2];
sx q[2];
rz(0.2790645) q[2];
rz(1.2532824) q[3];
sx q[3];
rz(-2.8918355) q[3];
sx q[3];
rz(2.9899924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-1.5308373) q[0];
sx q[0];
rz(-0.84722561) q[0];
sx q[0];
rz(-2.8952059) q[0];
rz(-1.1950182) q[1];
sx q[1];
rz(-2.2476826) q[1];
sx q[1];
rz(1.6349207) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5824597) q[0];
sx q[0];
rz(-1.8163067) q[0];
sx q[0];
rz(-1.1051635) q[0];
rz(-1.4343403) q[2];
sx q[2];
rz(-0.30559691) q[2];
sx q[2];
rz(-1.157925) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3132224) q[1];
sx q[1];
rz(-1.5130318) q[1];
sx q[1];
rz(1.934113) q[1];
x q[2];
rz(2.706932) q[3];
sx q[3];
rz(-2.0512329) q[3];
sx q[3];
rz(-0.91033376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.40111497) q[2];
sx q[2];
rz(-0.97492188) q[2];
sx q[2];
rz(-1.231989) q[2];
rz(1.1022107) q[3];
sx q[3];
rz(-3.1372742) q[3];
sx q[3];
rz(-1.8050885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44322893) q[0];
sx q[0];
rz(-2.4446428) q[0];
sx q[0];
rz(-0.90782905) q[0];
rz(-0.56697956) q[1];
sx q[1];
rz(-0.98273977) q[1];
sx q[1];
rz(3.0388015) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11353569) q[0];
sx q[0];
rz(-2.5350219) q[0];
sx q[0];
rz(-2.7675178) q[0];
rz(-2.2305713) q[2];
sx q[2];
rz(-3.0098923) q[2];
sx q[2];
rz(-2.1610799) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6425848) q[1];
sx q[1];
rz(-0.73621589) q[1];
sx q[1];
rz(-0.6249439) q[1];
x q[2];
rz(-1.0877043) q[3];
sx q[3];
rz(-2.5317041) q[3];
sx q[3];
rz(2.1263378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8583782) q[2];
sx q[2];
rz(-2.3329222) q[2];
sx q[2];
rz(1.4374479) q[2];
rz(-0.39250675) q[3];
sx q[3];
rz(-1.1350574) q[3];
sx q[3];
rz(-2.8431622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(2.0482386) q[0];
sx q[0];
rz(-1.3627351) q[0];
sx q[0];
rz(-1.1887953) q[0];
rz(-2.8554754) q[1];
sx q[1];
rz(-0.61758271) q[1];
sx q[1];
rz(1.6292705) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28933576) q[0];
sx q[0];
rz(-1.8484383) q[0];
sx q[0];
rz(-2.1136978) q[0];
rz(-pi) q[1];
rz(1.9061162) q[2];
sx q[2];
rz(-2.309707) q[2];
sx q[2];
rz(0.12140935) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0551853) q[1];
sx q[1];
rz(-2.1595104) q[1];
sx q[1];
rz(0.22365315) q[1];
x q[2];
rz(-1.2028473) q[3];
sx q[3];
rz(-0.5371437) q[3];
sx q[3];
rz(2.9972671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.35147038) q[2];
sx q[2];
rz(-2.0581547) q[2];
sx q[2];
rz(0.67698395) q[2];
rz(-1.2466768) q[3];
sx q[3];
rz(-1.2723943) q[3];
sx q[3];
rz(1.5164794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-2.2077654) q[0];
sx q[0];
rz(-2.8849869) q[0];
sx q[0];
rz(3.1166792) q[0];
rz(2.086153) q[1];
sx q[1];
rz(-2.1565304) q[1];
sx q[1];
rz(-0.28344646) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.442207) q[0];
sx q[0];
rz(-1.3106579) q[0];
sx q[0];
rz(0.26630638) q[0];
x q[1];
rz(1.3749429) q[2];
sx q[2];
rz(-1.4417267) q[2];
sx q[2];
rz(-1.376525) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.746096) q[1];
sx q[1];
rz(-2.0714134) q[1];
sx q[1];
rz(3.1067645) q[1];
x q[2];
rz(1.9881387) q[3];
sx q[3];
rz(-0.71708369) q[3];
sx q[3];
rz(-2.5337608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8289566) q[2];
sx q[2];
rz(-2.7241311) q[2];
sx q[2];
rz(0.32583315) q[2];
rz(1.2827986) q[3];
sx q[3];
rz(-1.9841586) q[3];
sx q[3];
rz(-1.5475984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40389898) q[0];
sx q[0];
rz(-1.6588545) q[0];
sx q[0];
rz(2.7666336) q[0];
rz(2.6629579) q[1];
sx q[1];
rz(-0.67239434) q[1];
sx q[1];
rz(0.37014827) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2787196) q[0];
sx q[0];
rz(-1.5476883) q[0];
sx q[0];
rz(1.5824759) q[0];
x q[1];
rz(0.92801969) q[2];
sx q[2];
rz(-0.22980873) q[2];
sx q[2];
rz(1.9191051) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3016175) q[1];
sx q[1];
rz(-2.2711922) q[1];
sx q[1];
rz(-0.61995929) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1305976) q[3];
sx q[3];
rz(-0.42783005) q[3];
sx q[3];
rz(1.4233936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.39773539) q[2];
sx q[2];
rz(-2.9501259) q[2];
sx q[2];
rz(-1.3810623) q[2];
rz(1.0487652) q[3];
sx q[3];
rz(-1.7966813) q[3];
sx q[3];
rz(0.63961187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25310707) q[0];
sx q[0];
rz(-2.3024004) q[0];
sx q[0];
rz(2.2174368) q[0];
rz(-0.91668516) q[1];
sx q[1];
rz(-1.0662096) q[1];
sx q[1];
rz(-1.7402657) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7766905) q[0];
sx q[0];
rz(-0.49743891) q[0];
sx q[0];
rz(-0.12125347) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9193814) q[2];
sx q[2];
rz(-0.66721081) q[2];
sx q[2];
rz(-1.8052342) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7622213) q[1];
sx q[1];
rz(-1.7306384) q[1];
sx q[1];
rz(-2.400378) q[1];
rz(-pi) q[2];
rz(-1.0823233) q[3];
sx q[3];
rz(-1.3783749) q[3];
sx q[3];
rz(0.50979641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9132729) q[2];
sx q[2];
rz(-1.6797804) q[2];
sx q[2];
rz(-1.224996) q[2];
rz(-0.0082958881) q[3];
sx q[3];
rz(-3.1305997) q[3];
sx q[3];
rz(-2.2744961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.55815721) q[0];
sx q[0];
rz(-0.83034101) q[0];
sx q[0];
rz(0.48026568) q[0];
rz(2.9971314) q[1];
sx q[1];
rz(-0.90765777) q[1];
sx q[1];
rz(-0.86437782) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9821765) q[0];
sx q[0];
rz(-1.5242531) q[0];
sx q[0];
rz(-1.4070562) q[0];
rz(-2.7988465) q[2];
sx q[2];
rz(-0.10048332) q[2];
sx q[2];
rz(-0.1048987) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.581199) q[1];
sx q[1];
rz(-0.59948363) q[1];
sx q[1];
rz(1.1349212) q[1];
rz(-pi) q[2];
rz(1.8501353) q[3];
sx q[3];
rz(-1.8350826) q[3];
sx q[3];
rz(-0.13388982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9353443) q[2];
sx q[2];
rz(-1.0103005) q[2];
sx q[2];
rz(-2.5065191) q[2];
rz(-0.57279974) q[3];
sx q[3];
rz(-1.3558931) q[3];
sx q[3];
rz(-0.87535453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0315345) q[0];
sx q[0];
rz(-0.77339554) q[0];
sx q[0];
rz(-3.0173259) q[0];
rz(0.36525137) q[1];
sx q[1];
rz(-1.7144014) q[1];
sx q[1];
rz(-1.7100547) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54473684) q[0];
sx q[0];
rz(-1.3787621) q[0];
sx q[0];
rz(1.7512683) q[0];
x q[1];
rz(-0.65406873) q[2];
sx q[2];
rz(-1.1032915) q[2];
sx q[2];
rz(1.2142177) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6811583) q[1];
sx q[1];
rz(-1.2660813) q[1];
sx q[1];
rz(2.4439993) q[1];
rz(2.0084736) q[3];
sx q[3];
rz(-2.056155) q[3];
sx q[3];
rz(-0.35972586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.26620904) q[2];
sx q[2];
rz(-0.92140809) q[2];
sx q[2];
rz(1.9688152) q[2];
rz(1.4069936) q[3];
sx q[3];
rz(-1.809779) q[3];
sx q[3];
rz(1.9269358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59239546) q[0];
sx q[0];
rz(-1.1347436) q[0];
sx q[0];
rz(-0.59463516) q[0];
rz(2.1356964) q[1];
sx q[1];
rz(-1.2024095) q[1];
sx q[1];
rz(-0.88919052) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2834463) q[0];
sx q[0];
rz(-0.22617561) q[0];
sx q[0];
rz(-1.7986138) q[0];
rz(-pi) q[1];
rz(3.0709549) q[2];
sx q[2];
rz(-2.0301691) q[2];
sx q[2];
rz(-0.86462155) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8264849) q[1];
sx q[1];
rz(-1.6651622) q[1];
sx q[1];
rz(3.0815275) q[1];
rz(-pi) q[2];
rz(-3.0746769) q[3];
sx q[3];
rz(-0.51735462) q[3];
sx q[3];
rz(0.87393119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.95120007) q[2];
sx q[2];
rz(-1.1756281) q[2];
sx q[2];
rz(2.5962043) q[2];
rz(-0.96327463) q[3];
sx q[3];
rz(-1.9656209) q[3];
sx q[3];
rz(1.6776599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49190285) q[0];
sx q[0];
rz(-2.2408673) q[0];
sx q[0];
rz(1.3229205) q[0];
rz(-0.84359618) q[1];
sx q[1];
rz(-1.0782764) q[1];
sx q[1];
rz(-1.2596399) q[1];
rz(2.9424473) q[2];
sx q[2];
rz(-2.0024588) q[2];
sx q[2];
rz(-1.3695516) q[2];
rz(1.0641742) q[3];
sx q[3];
rz(-2.0941877) q[3];
sx q[3];
rz(0.82529395) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
