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
rz(-1.7614814) q[0];
sx q[0];
rz(-1.5505646) q[0];
sx q[0];
rz(0.38183364) q[0];
rz(2.7497357) q[1];
sx q[1];
rz(-1.5500103) q[1];
sx q[1];
rz(-0.92441192) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3782717) q[0];
sx q[0];
rz(-2.1594783) q[0];
sx q[0];
rz(0.6529385) q[0];
x q[1];
rz(-1.5509694) q[2];
sx q[2];
rz(-1.2006491) q[2];
sx q[2];
rz(-1.4408979) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0142483) q[1];
sx q[1];
rz(-2.6963628) q[1];
sx q[1];
rz(-0.30502747) q[1];
rz(2.1006868) q[3];
sx q[3];
rz(-1.232551) q[3];
sx q[3];
rz(1.4640946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9367289) q[2];
sx q[2];
rz(-1.8483346) q[2];
sx q[2];
rz(1.5748242) q[2];
rz(-1.9751679) q[3];
sx q[3];
rz(-2.0765442) q[3];
sx q[3];
rz(-0.69275698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4676056) q[0];
sx q[0];
rz(-0.80802149) q[0];
sx q[0];
rz(1.9655193) q[0];
rz(0.067497079) q[1];
sx q[1];
rz(-2.564023) q[1];
sx q[1];
rz(-0.31879058) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97385308) q[0];
sx q[0];
rz(-1.8813846) q[0];
sx q[0];
rz(-2.3936901) q[0];
rz(-pi) q[1];
rz(1.8068878) q[2];
sx q[2];
rz(-1.9892879) q[2];
sx q[2];
rz(-0.22007569) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8242256) q[1];
sx q[1];
rz(-1.1928802) q[1];
sx q[1];
rz(2.1721858) q[1];
rz(-1.2845006) q[3];
sx q[3];
rz(-1.4501374) q[3];
sx q[3];
rz(-0.21721043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0789644) q[2];
sx q[2];
rz(-2.2603409) q[2];
sx q[2];
rz(0.47743615) q[2];
rz(-1.7834974) q[3];
sx q[3];
rz(-2.4886459) q[3];
sx q[3];
rz(3.131598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39182144) q[0];
sx q[0];
rz(-1.3359767) q[0];
sx q[0];
rz(2.0023477) q[0];
rz(-0.40621743) q[1];
sx q[1];
rz(-1.258305) q[1];
sx q[1];
rz(-0.67799813) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0292873) q[0];
sx q[0];
rz(-2.5480518) q[0];
sx q[0];
rz(1.4112006) q[0];
rz(-pi) q[1];
rz(2.1966372) q[2];
sx q[2];
rz(-1.2765108) q[2];
sx q[2];
rz(-1.6950032) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.78732508) q[1];
sx q[1];
rz(-2.4128822) q[1];
sx q[1];
rz(-2.4142152) q[1];
rz(-0.33379995) q[3];
sx q[3];
rz(-1.085328) q[3];
sx q[3];
rz(-0.36241058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4262126) q[2];
sx q[2];
rz(-1.6343626) q[2];
sx q[2];
rz(-2.0646084) q[2];
rz(1.2601323) q[3];
sx q[3];
rz(-1.0271122) q[3];
sx q[3];
rz(-0.25585678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28427163) q[0];
sx q[0];
rz(-1.1742598) q[0];
sx q[0];
rz(-0.76048365) q[0];
rz(-0.35176945) q[1];
sx q[1];
rz(-2.4302509) q[1];
sx q[1];
rz(2.9913091) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1635123) q[0];
sx q[0];
rz(-0.0038359782) q[0];
sx q[0];
rz(-1.0413961) q[0];
rz(0.17268659) q[2];
sx q[2];
rz(-1.193668) q[2];
sx q[2];
rz(-1.0972301) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1898415) q[1];
sx q[1];
rz(-0.69165666) q[1];
sx q[1];
rz(-1.0361978) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.45771367) q[3];
sx q[3];
rz(-1.8031075) q[3];
sx q[3];
rz(-1.4087698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1111697) q[2];
sx q[2];
rz(-0.22746484) q[2];
sx q[2];
rz(-0.58491582) q[2];
rz(3.0758744) q[3];
sx q[3];
rz(-1.1394371) q[3];
sx q[3];
rz(2.0871302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2718662) q[0];
sx q[0];
rz(-1.2249185) q[0];
sx q[0];
rz(-2.4140893) q[0];
rz(-1.2959405) q[1];
sx q[1];
rz(-1.0546874) q[1];
sx q[1];
rz(-0.71472275) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22435221) q[0];
sx q[0];
rz(-0.2324129) q[0];
sx q[0];
rz(0.31347855) q[0];
rz(-pi) q[1];
x q[1];
rz(0.3008414) q[2];
sx q[2];
rz(-1.6812426) q[2];
sx q[2];
rz(-3.0120603) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2421869) q[1];
sx q[1];
rz(-1.9493628) q[1];
sx q[1];
rz(2.942571) q[1];
rz(-pi) q[2];
rz(-0.994579) q[3];
sx q[3];
rz(-1.2354697) q[3];
sx q[3];
rz(-1.9063584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1456445) q[2];
sx q[2];
rz(-1.6137292) q[2];
sx q[2];
rz(-1.145251) q[2];
rz(-0.20259419) q[3];
sx q[3];
rz(-3.0718206) q[3];
sx q[3];
rz(-0.32874671) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2957434) q[0];
sx q[0];
rz(-0.29243094) q[0];
sx q[0];
rz(0.16786815) q[0];
rz(0.56495086) q[1];
sx q[1];
rz(-2.0940557) q[1];
sx q[1];
rz(-1.8126743) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5707626) q[0];
sx q[0];
rz(-2.0432297) q[0];
sx q[0];
rz(-0.60586141) q[0];
rz(-pi) q[1];
rz(2.8064427) q[2];
sx q[2];
rz(-2.9932635) q[2];
sx q[2];
rz(0.89491612) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4327287) q[1];
sx q[1];
rz(-1.8058746) q[1];
sx q[1];
rz(2.8793596) q[1];
x q[2];
rz(0.22340138) q[3];
sx q[3];
rz(-1.3376682) q[3];
sx q[3];
rz(0.65321556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7557175) q[2];
sx q[2];
rz(-1.2473325) q[2];
sx q[2];
rz(0.34783777) q[2];
rz(-2.8268585) q[3];
sx q[3];
rz(-2.0458524) q[3];
sx q[3];
rz(1.1381963) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9921853) q[0];
sx q[0];
rz(-1.6353761) q[0];
sx q[0];
rz(-0.95112479) q[0];
rz(2.8490207) q[1];
sx q[1];
rz(-2.2508299) q[1];
sx q[1];
rz(2.1975885) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4571465) q[0];
sx q[0];
rz(-1.5553885) q[0];
sx q[0];
rz(-1.8788565) q[0];
rz(2.3500309) q[2];
sx q[2];
rz(-2.5063609) q[2];
sx q[2];
rz(-0.61295618) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.71783644) q[1];
sx q[1];
rz(-1.5482117) q[1];
sx q[1];
rz(1.7439688) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9721479) q[3];
sx q[3];
rz(-1.7987393) q[3];
sx q[3];
rz(0.55647382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.49147478) q[2];
sx q[2];
rz(-2.398141) q[2];
sx q[2];
rz(1.4415119) q[2];
rz(0.024638351) q[3];
sx q[3];
rz(-1.8162497) q[3];
sx q[3];
rz(-2.1520481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3330419) q[0];
sx q[0];
rz(-1.4104183) q[0];
sx q[0];
rz(0.1380052) q[0];
rz(2.0637312) q[1];
sx q[1];
rz(-1.6777918) q[1];
sx q[1];
rz(-0.93625751) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2168343) q[0];
sx q[0];
rz(-0.48473293) q[0];
sx q[0];
rz(0.68883552) q[0];
x q[1];
rz(0.9844043) q[2];
sx q[2];
rz(-2.3086562) q[2];
sx q[2];
rz(-3.0393103) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1106134) q[1];
sx q[1];
rz(-0.18337164) q[1];
sx q[1];
rz(1.490209) q[1];
rz(2.8209723) q[3];
sx q[3];
rz(-0.71133876) q[3];
sx q[3];
rz(-2.5742755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1428947) q[2];
sx q[2];
rz(-1.8691209) q[2];
sx q[2];
rz(-1.8227089) q[2];
rz(3.0190492) q[3];
sx q[3];
rz(-0.69869852) q[3];
sx q[3];
rz(-0.43012777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1634624) q[0];
sx q[0];
rz(-0.4438816) q[0];
sx q[0];
rz(2.7769856) q[0];
rz(-1.7568024) q[1];
sx q[1];
rz(-2.2099647) q[1];
sx q[1];
rz(-2.2273831) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96959463) q[0];
sx q[0];
rz(-1.9978317) q[0];
sx q[0];
rz(0.14694218) q[0];
rz(-0.73458521) q[2];
sx q[2];
rz(-1.5365639) q[2];
sx q[2];
rz(-1.5058668) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.89489854) q[1];
sx q[1];
rz(-1.7082038) q[1];
sx q[1];
rz(0.64260428) q[1];
rz(1.4214324) q[3];
sx q[3];
rz(-1.0348667) q[3];
sx q[3];
rz(-1.6370629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8387575) q[2];
sx q[2];
rz(-0.5842394) q[2];
sx q[2];
rz(-1.5128822) q[2];
rz(-2.5527939) q[3];
sx q[3];
rz(-1.9353119) q[3];
sx q[3];
rz(-2.485386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6673679) q[0];
sx q[0];
rz(-0.56741699) q[0];
sx q[0];
rz(-2.8571416) q[0];
rz(2.4868763) q[1];
sx q[1];
rz(-1.7660716) q[1];
sx q[1];
rz(2.7213352) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3271844) q[0];
sx q[0];
rz(-1.1564009) q[0];
sx q[0];
rz(2.9073858) q[0];
rz(2.3510397) q[2];
sx q[2];
rz(-1.2324863) q[2];
sx q[2];
rz(-1.4684791) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0033231) q[1];
sx q[1];
rz(-1.6594634) q[1];
sx q[1];
rz(1.4515463) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5758366) q[3];
sx q[3];
rz(-1.8744191) q[3];
sx q[3];
rz(-0.88241258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1493211) q[2];
sx q[2];
rz(-2.6565266) q[2];
sx q[2];
rz(0.9672271) q[2];
rz(-1.017717) q[3];
sx q[3];
rz(-1.8571721) q[3];
sx q[3];
rz(2.4093936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16348542) q[0];
sx q[0];
rz(-1.0657943) q[0];
sx q[0];
rz(-0.97070538) q[0];
rz(-0.12374395) q[1];
sx q[1];
rz(-0.94656222) q[1];
sx q[1];
rz(1.8306517) q[1];
rz(-0.38519771) q[2];
sx q[2];
rz(-2.7976947) q[2];
sx q[2];
rz(-1.7277389) q[2];
rz(-1.6061892) q[3];
sx q[3];
rz(-1.8826859) q[3];
sx q[3];
rz(-1.7425698) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
