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
rz(-1.6753766) q[0];
sx q[0];
rz(-2.197062) q[0];
sx q[0];
rz(-0.20139995) q[0];
rz(1.072285) q[1];
sx q[1];
rz(-0.76289248) q[1];
sx q[1];
rz(1.4210757) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3823377) q[0];
sx q[0];
rz(-1.5159642) q[0];
sx q[0];
rz(2.1817529) q[0];
x q[1];
rz(-0.95008738) q[2];
sx q[2];
rz(-2.3664775) q[2];
sx q[2];
rz(-0.24573869) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0951251) q[1];
sx q[1];
rz(-0.78739843) q[1];
sx q[1];
rz(-2.5745029) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9626612) q[3];
sx q[3];
rz(-1.8800991) q[3];
sx q[3];
rz(-0.40180692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1066771) q[2];
sx q[2];
rz(-0.80433977) q[2];
sx q[2];
rz(2.8625281) q[2];
rz(1.8883102) q[3];
sx q[3];
rz(-0.24975714) q[3];
sx q[3];
rz(2.9899924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6107553) q[0];
sx q[0];
rz(-0.84722561) q[0];
sx q[0];
rz(2.8952059) q[0];
rz(-1.1950182) q[1];
sx q[1];
rz(-0.89391005) q[1];
sx q[1];
rz(1.5066719) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89013571) q[0];
sx q[0];
rz(-1.1201753) q[0];
sx q[0];
rz(-2.8681953) q[0];
rz(-1.7072524) q[2];
sx q[2];
rz(-0.30559691) q[2];
sx q[2];
rz(1.157925) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8283702) q[1];
sx q[1];
rz(-1.6285609) q[1];
sx q[1];
rz(-1.934113) q[1];
x q[2];
rz(2.706932) q[3];
sx q[3];
rz(-2.0512329) q[3];
sx q[3];
rz(2.2312589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.40111497) q[2];
sx q[2];
rz(-0.97492188) q[2];
sx q[2];
rz(-1.9096036) q[2];
rz(1.1022107) q[3];
sx q[3];
rz(-3.1372742) q[3];
sx q[3];
rz(-1.8050885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44322893) q[0];
sx q[0];
rz(-2.4446428) q[0];
sx q[0];
rz(-2.2337636) q[0];
rz(0.56697956) q[1];
sx q[1];
rz(-0.98273977) q[1];
sx q[1];
rz(-3.0388015) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11353569) q[0];
sx q[0];
rz(-2.5350219) q[0];
sx q[0];
rz(0.3740749) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91102131) q[2];
sx q[2];
rz(-3.0098923) q[2];
sx q[2];
rz(0.98051276) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4145706) q[1];
sx q[1];
rz(-2.1466781) q[1];
sx q[1];
rz(1.0833137) q[1];
rz(-pi) q[2];
rz(1.0166753) q[3];
sx q[3];
rz(-1.8401056) q[3];
sx q[3];
rz(2.179972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8583782) q[2];
sx q[2];
rz(-2.3329222) q[2];
sx q[2];
rz(-1.7041448) q[2];
rz(0.39250675) q[3];
sx q[3];
rz(-1.1350574) q[3];
sx q[3];
rz(2.8431622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0933541) q[0];
sx q[0];
rz(-1.7788576) q[0];
sx q[0];
rz(1.9527973) q[0];
rz(-2.8554754) q[1];
sx q[1];
rz(-0.61758271) q[1];
sx q[1];
rz(1.6292705) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.286519) q[0];
sx q[0];
rz(-2.5382156) q[0];
sx q[0];
rz(1.0666749) q[0];
x q[1];
rz(2.3740923) q[2];
sx q[2];
rz(-1.3250809) q[2];
sx q[2];
rz(-1.9227288) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.69763598) q[1];
sx q[1];
rz(-0.62503615) q[1];
sx q[1];
rz(1.8915063) q[1];
rz(-pi) q[2];
rz(2.9305601) q[3];
sx q[3];
rz(-1.0730626) q[3];
sx q[3];
rz(-0.56609234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7901223) q[2];
sx q[2];
rz(-1.0834379) q[2];
sx q[2];
rz(2.4646087) q[2];
rz(-1.2466768) q[3];
sx q[3];
rz(-1.2723943) q[3];
sx q[3];
rz(1.5164794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2077654) q[0];
sx q[0];
rz(-0.25660577) q[0];
sx q[0];
rz(-0.024913464) q[0];
rz(1.0554396) q[1];
sx q[1];
rz(-0.98506227) q[1];
sx q[1];
rz(2.8581462) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0830529) q[0];
sx q[0];
rz(-1.3136567) q[0];
sx q[0];
rz(-1.8399946) q[0];
rz(-pi) q[1];
rz(-0.9825969) q[2];
sx q[2];
rz(-0.23410205) q[2];
sx q[2];
rz(-2.7603619) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.746096) q[1];
sx q[1];
rz(-1.0701792) q[1];
sx q[1];
rz(0.034828111) q[1];
x q[2];
rz(-1.9881387) q[3];
sx q[3];
rz(-2.424509) q[3];
sx q[3];
rz(-2.5337608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.31263605) q[2];
sx q[2];
rz(-2.7241311) q[2];
sx q[2];
rz(-0.32583315) q[2];
rz(1.2827986) q[3];
sx q[3];
rz(-1.157434) q[3];
sx q[3];
rz(-1.5939943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
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
rz(0.40389898) q[0];
sx q[0];
rz(-1.6588545) q[0];
sx q[0];
rz(-2.7666336) q[0];
rz(-0.47863475) q[1];
sx q[1];
rz(-2.4691983) q[1];
sx q[1];
rz(2.7714444) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2787196) q[0];
sx q[0];
rz(-1.5939043) q[0];
sx q[0];
rz(1.5591168) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7559075) q[2];
sx q[2];
rz(-1.7077669) q[2];
sx q[2];
rz(0.97835195) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3016175) q[1];
sx q[1];
rz(-2.2711922) q[1];
sx q[1];
rz(0.61995929) q[1];
rz(1.5657827) q[3];
sx q[3];
rz(-1.1429938) q[3];
sx q[3];
rz(-1.4354777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7438573) q[2];
sx q[2];
rz(-0.19146679) q[2];
sx q[2];
rz(-1.7605304) q[2];
rz(-1.0487652) q[3];
sx q[3];
rz(-1.3449113) q[3];
sx q[3];
rz(-2.5019808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.8884856) q[0];
sx q[0];
rz(-0.83919224) q[0];
sx q[0];
rz(2.2174368) q[0];
rz(0.91668516) q[1];
sx q[1];
rz(-2.075383) q[1];
sx q[1];
rz(1.4013269) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2271254) q[0];
sx q[0];
rz(-2.0642529) q[0];
sx q[0];
rz(1.5052133) q[0];
rz(-pi) q[1];
rz(1.7426874) q[2];
sx q[2];
rz(-0.92280932) q[2];
sx q[2];
rz(-1.5251708) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7622213) q[1];
sx q[1];
rz(-1.7306384) q[1];
sx q[1];
rz(-0.74121468) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1772778) q[3];
sx q[3];
rz(-2.6194508) q[3];
sx q[3];
rz(-2.4260184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2283198) q[2];
sx q[2];
rz(-1.4618123) q[2];
sx q[2];
rz(-1.9165967) q[2];
rz(0.0082958881) q[3];
sx q[3];
rz(-0.010992916) q[3];
sx q[3];
rz(0.86709658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5834354) q[0];
sx q[0];
rz(-2.3112516) q[0];
sx q[0];
rz(-0.48026568) q[0];
rz(-0.14446124) q[1];
sx q[1];
rz(-0.90765777) q[1];
sx q[1];
rz(2.2772148) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9821765) q[0];
sx q[0];
rz(-1.6173395) q[0];
sx q[0];
rz(-1.4070562) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7988465) q[2];
sx q[2];
rz(-3.0411093) q[2];
sx q[2];
rz(-3.036694) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0675812) q[1];
sx q[1];
rz(-2.1076822) q[1];
sx q[1];
rz(-0.28089653) q[1];
x q[2];
rz(-1.2914574) q[3];
sx q[3];
rz(-1.8350826) q[3];
sx q[3];
rz(3.0077028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.20624837) q[2];
sx q[2];
rz(-2.1312921) q[2];
sx q[2];
rz(-0.63507357) q[2];
rz(-0.57279974) q[3];
sx q[3];
rz(-1.7856995) q[3];
sx q[3];
rz(-2.2662381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0315345) q[0];
sx q[0];
rz(-0.77339554) q[0];
sx q[0];
rz(3.0173259) q[0];
rz(-0.36525137) q[1];
sx q[1];
rz(-1.4271913) q[1];
sx q[1];
rz(-1.7100547) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54473684) q[0];
sx q[0];
rz(-1.3787621) q[0];
sx q[0];
rz(-1.7512683) q[0];
x q[1];
rz(1.0042436) q[2];
sx q[2];
rz(-2.1449617) q[2];
sx q[2];
rz(3.1176709) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4540951) q[1];
sx q[1];
rz(-2.3907067) q[1];
sx q[1];
rz(0.45529699) q[1];
rz(-pi) q[2];
rz(2.6142653) q[3];
sx q[3];
rz(-1.1865215) q[3];
sx q[3];
rz(-1.4260071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8753836) q[2];
sx q[2];
rz(-0.92140809) q[2];
sx q[2];
rz(1.1727775) q[2];
rz(1.734599) q[3];
sx q[3];
rz(-1.809779) q[3];
sx q[3];
rz(-1.9269358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5491972) q[0];
sx q[0];
rz(-1.1347436) q[0];
sx q[0];
rz(2.5469575) q[0];
rz(-1.0058962) q[1];
sx q[1];
rz(-1.2024095) q[1];
sx q[1];
rz(-0.88919052) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62458663) q[0];
sx q[0];
rz(-1.7910302) q[0];
sx q[0];
rz(3.0896679) q[0];
x q[1];
rz(-3.0709549) q[2];
sx q[2];
rz(-1.1114235) q[2];
sx q[2];
rz(2.2769711) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2500221) q[1];
sx q[1];
rz(-1.5109987) q[1];
sx q[1];
rz(1.6653316) q[1];
rz(-3.0746769) q[3];
sx q[3];
rz(-0.51735462) q[3];
sx q[3];
rz(0.87393119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.95120007) q[2];
sx q[2];
rz(-1.9659646) q[2];
sx q[2];
rz(-2.5962043) q[2];
rz(2.178318) q[3];
sx q[3];
rz(-1.9656209) q[3];
sx q[3];
rz(1.6776599) q[3];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6496898) q[0];
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
rz(-0.58335494) q[3];
sx q[3];
rz(-1.1370549) q[3];
sx q[3];
rz(2.1255253) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
