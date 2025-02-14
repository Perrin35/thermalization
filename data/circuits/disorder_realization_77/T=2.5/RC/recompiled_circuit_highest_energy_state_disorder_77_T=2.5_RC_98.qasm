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
rz(-0.4975118) q[0];
sx q[0];
rz(4.4805718) q[0];
sx q[0];
rz(6.5988402) q[0];
rz(1.1701801) q[1];
sx q[1];
rz(-0.38771954) q[1];
sx q[1];
rz(-2.5375836) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7496628) q[0];
sx q[0];
rz(-1.6214341) q[0];
sx q[0];
rz(-1.873436) q[0];
rz(-2.6409615) q[2];
sx q[2];
rz(-0.96304578) q[2];
sx q[2];
rz(-1.0333824) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.042272557) q[1];
sx q[1];
rz(-1.0612773) q[1];
sx q[1];
rz(0.38739844) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8995057) q[3];
sx q[3];
rz(-2.8162873) q[3];
sx q[3];
rz(1.5055552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.151256) q[2];
sx q[2];
rz(-1.1948816) q[2];
sx q[2];
rz(-1.3585496) q[2];
rz(1.547706) q[3];
sx q[3];
rz(-2.2113776) q[3];
sx q[3];
rz(-1.5096629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5556521) q[0];
sx q[0];
rz(-2.5037615) q[0];
sx q[0];
rz(-1.1974539) q[0];
rz(-2.6400631) q[1];
sx q[1];
rz(-1.5781559) q[1];
sx q[1];
rz(1.4362358) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2510808) q[0];
sx q[0];
rz(-2.5174826) q[0];
sx q[0];
rz(1.1680383) q[0];
x q[1];
rz(-3.0244963) q[2];
sx q[2];
rz(-0.70636049) q[2];
sx q[2];
rz(1.7603086) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7926678) q[1];
sx q[1];
rz(-1.7670005) q[1];
sx q[1];
rz(-1.6755107) q[1];
x q[2];
rz(0.45577502) q[3];
sx q[3];
rz(-1.7592128) q[3];
sx q[3];
rz(0.51793232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2531835) q[2];
sx q[2];
rz(-1.9056355) q[2];
sx q[2];
rz(3.1186228) q[2];
rz(-2.2394771) q[3];
sx q[3];
rz(-1.918856) q[3];
sx q[3];
rz(2.7149916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70812923) q[0];
sx q[0];
rz(-1.3490278) q[0];
sx q[0];
rz(-0.9915114) q[0];
rz(2.1082711) q[1];
sx q[1];
rz(-0.28305498) q[1];
sx q[1];
rz(1.2352157) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9231222) q[0];
sx q[0];
rz(-1.8184796) q[0];
sx q[0];
rz(-2.5367303) q[0];
rz(-pi) q[1];
rz(-2.5318145) q[2];
sx q[2];
rz(-1.1705361) q[2];
sx q[2];
rz(2.5230809) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.24656235) q[1];
sx q[1];
rz(-0.8324648) q[1];
sx q[1];
rz(-2.8042996) q[1];
x q[2];
rz(0.83241141) q[3];
sx q[3];
rz(-0.73736546) q[3];
sx q[3];
rz(-2.5016145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.11423763) q[2];
sx q[2];
rz(-2.3894252) q[2];
sx q[2];
rz(2.1503964) q[2];
rz(1.8441955) q[3];
sx q[3];
rz(-1.2605366) q[3];
sx q[3];
rz(-1.4170925) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39795136) q[0];
sx q[0];
rz(-0.93248168) q[0];
sx q[0];
rz(-2.3241296) q[0];
rz(0.3262597) q[1];
sx q[1];
rz(-1.1810818) q[1];
sx q[1];
rz(0.78525966) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24910771) q[0];
sx q[0];
rz(-1.8726702) q[0];
sx q[0];
rz(-2.8623146) q[0];
x q[1];
rz(-0.26008028) q[2];
sx q[2];
rz(-2.4818015) q[2];
sx q[2];
rz(-1.5957956) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1568415) q[1];
sx q[1];
rz(-1.1492711) q[1];
sx q[1];
rz(-2.0578029) q[1];
rz(-pi) q[2];
rz(-2.9296285) q[3];
sx q[3];
rz(-0.3623687) q[3];
sx q[3];
rz(-0.22892287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.12104812) q[2];
sx q[2];
rz(-0.64457568) q[2];
sx q[2];
rz(-0.56309593) q[2];
rz(-2.2480887) q[3];
sx q[3];
rz(-1.1734791) q[3];
sx q[3];
rz(1.9068498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23984443) q[0];
sx q[0];
rz(-0.32308602) q[0];
sx q[0];
rz(1.595994) q[0];
rz(-1.0218989) q[1];
sx q[1];
rz(-1.1232168) q[1];
sx q[1];
rz(0.18128577) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90732987) q[0];
sx q[0];
rz(-1.0815623) q[0];
sx q[0];
rz(0.2184043) q[0];
rz(-pi) q[1];
rz(0.092575707) q[2];
sx q[2];
rz(-0.47738722) q[2];
sx q[2];
rz(0.30711781) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.681057) q[1];
sx q[1];
rz(-1.3486282) q[1];
sx q[1];
rz(0.028789) q[1];
rz(1.7277754) q[3];
sx q[3];
rz(-2.008956) q[3];
sx q[3];
rz(-2.639132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1229317) q[2];
sx q[2];
rz(-0.67492008) q[2];
sx q[2];
rz(1.2459416) q[2];
rz(2.5490226) q[3];
sx q[3];
rz(-1.4453245) q[3];
sx q[3];
rz(2.3004801) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5331921) q[0];
sx q[0];
rz(-0.99219457) q[0];
sx q[0];
rz(0.30884185) q[0];
rz(-2.8187075) q[1];
sx q[1];
rz(-0.37728089) q[1];
sx q[1];
rz(3.0016532) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2353781) q[0];
sx q[0];
rz(-1.8188634) q[0];
sx q[0];
rz(-2.1885022) q[0];
rz(-pi) q[1];
rz(-1.3792737) q[2];
sx q[2];
rz(-2.6783248) q[2];
sx q[2];
rz(0.5199711) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6778826) q[1];
sx q[1];
rz(-1.8720798) q[1];
sx q[1];
rz(0.82656411) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30045565) q[3];
sx q[3];
rz(-2.2078035) q[3];
sx q[3];
rz(-0.14618044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.66190019) q[2];
sx q[2];
rz(-0.65325824) q[2];
sx q[2];
rz(2.8996186) q[2];
rz(0.74609977) q[3];
sx q[3];
rz(-1.3662162) q[3];
sx q[3];
rz(1.411875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11874966) q[0];
sx q[0];
rz(-2.0976837) q[0];
sx q[0];
rz(-1.2710849) q[0];
rz(-2.5785043) q[1];
sx q[1];
rz(-2.2124898) q[1];
sx q[1];
rz(1.7005327) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6800425) q[0];
sx q[0];
rz(-0.6667887) q[0];
sx q[0];
rz(2.0186485) q[0];
x q[1];
rz(1.7092429) q[2];
sx q[2];
rz(-1.9488153) q[2];
sx q[2];
rz(2.0496617) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.043277117) q[1];
sx q[1];
rz(-0.47887427) q[1];
sx q[1];
rz(-0.4988341) q[1];
x q[2];
rz(-2.8981179) q[3];
sx q[3];
rz(-2.4015275) q[3];
sx q[3];
rz(-2.8748663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9219804) q[2];
sx q[2];
rz(-1.8461123) q[2];
sx q[2];
rz(1.316635) q[2];
rz(-2.5804139) q[3];
sx q[3];
rz(-2.1771274) q[3];
sx q[3];
rz(1.6130028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0372666) q[0];
sx q[0];
rz(-2.9236561) q[0];
sx q[0];
rz(-2.2547146) q[0];
rz(-0.2001702) q[1];
sx q[1];
rz(-1.6693516) q[1];
sx q[1];
rz(-2.1194469) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0933233) q[0];
sx q[0];
rz(-1.5453858) q[0];
sx q[0];
rz(1.4396458) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.736945) q[2];
sx q[2];
rz(-0.96074694) q[2];
sx q[2];
rz(1.7231154) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2547438) q[1];
sx q[1];
rz(-1.5690156) q[1];
sx q[1];
rz(0.82474553) q[1];
x q[2];
rz(-1.2194013) q[3];
sx q[3];
rz(-2.1179652) q[3];
sx q[3];
rz(1.6148293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6518121) q[2];
sx q[2];
rz(-0.45280364) q[2];
sx q[2];
rz(-0.81857267) q[2];
rz(0.98322785) q[3];
sx q[3];
rz(-2.1849617) q[3];
sx q[3];
rz(-2.6035068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1061358) q[0];
sx q[0];
rz(-1.9945194) q[0];
sx q[0];
rz(1.5552833) q[0];
rz(-2.4447794) q[1];
sx q[1];
rz(-2.9882444) q[1];
sx q[1];
rz(0.78479016) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0138088) q[0];
sx q[0];
rz(-2.0079566) q[0];
sx q[0];
rz(0.60370914) q[0];
rz(-pi) q[1];
rz(2.9892342) q[2];
sx q[2];
rz(-1.3703642) q[2];
sx q[2];
rz(2.7510425) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5604374) q[1];
sx q[1];
rz(-1.135901) q[1];
sx q[1];
rz(1.5337318) q[1];
x q[2];
rz(-2.574563) q[3];
sx q[3];
rz(-2.4483557) q[3];
sx q[3];
rz(0.45613134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4697504) q[2];
sx q[2];
rz(-1.2776813) q[2];
sx q[2];
rz(2.3308241) q[2];
rz(-1.9226711) q[3];
sx q[3];
rz(-0.54088497) q[3];
sx q[3];
rz(-2.0764652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.78275457) q[0];
sx q[0];
rz(-0.29958075) q[0];
sx q[0];
rz(1.7175571) q[0];
rz(0.77761039) q[1];
sx q[1];
rz(-0.97336665) q[1];
sx q[1];
rz(-0.75540677) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7292405) q[0];
sx q[0];
rz(-1.5753813) q[0];
sx q[0];
rz(-1.3608749) q[0];
x q[1];
rz(-1.5883114) q[2];
sx q[2];
rz(-1.9672158) q[2];
sx q[2];
rz(-2.6902386) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0499784) q[1];
sx q[1];
rz(-2.0777933) q[1];
sx q[1];
rz(1.9732287) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.741119) q[3];
sx q[3];
rz(-0.71610928) q[3];
sx q[3];
rz(0.082883714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5648254) q[2];
sx q[2];
rz(-0.29594031) q[2];
sx q[2];
rz(2.9887065) q[2];
rz(-1.8704002) q[3];
sx q[3];
rz(-1.8891687) q[3];
sx q[3];
rz(-0.66711867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85970238) q[0];
sx q[0];
rz(-2.6489881) q[0];
sx q[0];
rz(-0.76982605) q[0];
rz(-1.9181171) q[1];
sx q[1];
rz(-1.9203095) q[1];
sx q[1];
rz(-0.65345678) q[1];
rz(-2.5443947) q[2];
sx q[2];
rz(-1.9382678) q[2];
sx q[2];
rz(1.1282327) q[2];
rz(2.7889403) q[3];
sx q[3];
rz(-1.9518387) q[3];
sx q[3];
rz(0.52458682) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
