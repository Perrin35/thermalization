OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9392202) q[0];
sx q[0];
rz(-0.4063172) q[0];
sx q[0];
rz(0.82011861) q[0];
rz(2.7804873) q[1];
sx q[1];
rz(-0.63280025) q[1];
sx q[1];
rz(-0.83067218) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2031189) q[0];
sx q[0];
rz(-1.9460558) q[0];
sx q[0];
rz(1.24036) q[0];
rz(-pi) q[1];
rz(2.0482424) q[2];
sx q[2];
rz(-0.89077836) q[2];
sx q[2];
rz(-2.8505461) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0349866) q[1];
sx q[1];
rz(-0.88284661) q[1];
sx q[1];
rz(0.41253849) q[1];
rz(-pi) q[2];
x q[2];
rz(0.72676267) q[3];
sx q[3];
rz(-0.57075497) q[3];
sx q[3];
rz(1.6970413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7044907) q[2];
sx q[2];
rz(-0.4318684) q[2];
sx q[2];
rz(0.1201771) q[2];
rz(1.1581356) q[3];
sx q[3];
rz(-1.7430256) q[3];
sx q[3];
rz(0.91896287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0607818) q[0];
sx q[0];
rz(-1.7827543) q[0];
sx q[0];
rz(-0.91180116) q[0];
rz(-0.78951019) q[1];
sx q[1];
rz(-0.98840886) q[1];
sx q[1];
rz(-2.8149014) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2355736) q[0];
sx q[0];
rz(-0.86750194) q[0];
sx q[0];
rz(0.88615449) q[0];
rz(-pi) q[1];
rz(-2.0147444) q[2];
sx q[2];
rz(-1.4269281) q[2];
sx q[2];
rz(2.4614046) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5888728) q[1];
sx q[1];
rz(-0.96541222) q[1];
sx q[1];
rz(2.6436716) q[1];
rz(0.085751199) q[3];
sx q[3];
rz(-0.88480703) q[3];
sx q[3];
rz(-3.0119197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6313173) q[2];
sx q[2];
rz(-0.77284208) q[2];
sx q[2];
rz(1.2871683) q[2];
rz(3.0316947) q[3];
sx q[3];
rz(-1.4108312) q[3];
sx q[3];
rz(-1.3818285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54863769) q[0];
sx q[0];
rz(-2.4023963) q[0];
sx q[0];
rz(0.32989311) q[0];
rz(2.864481) q[1];
sx q[1];
rz(-1.3169293) q[1];
sx q[1];
rz(-1.057391) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43305106) q[0];
sx q[0];
rz(-1.2305224) q[0];
sx q[0];
rz(0.021854594) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3710652) q[2];
sx q[2];
rz(-1.235851) q[2];
sx q[2];
rz(-1.0945601) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1765269) q[1];
sx q[1];
rz(-2.2245193) q[1];
sx q[1];
rz(0.42018004) q[1];
rz(1.3789165) q[3];
sx q[3];
rz(-1.1550316) q[3];
sx q[3];
rz(0.92431812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.7827591) q[2];
sx q[2];
rz(-2.0262572) q[2];
sx q[2];
rz(1.7791629) q[2];
rz(0.6247012) q[3];
sx q[3];
rz(-1.0995068) q[3];
sx q[3];
rz(-2.3220298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3574922) q[0];
sx q[0];
rz(-0.52629137) q[0];
sx q[0];
rz(-2.6065361) q[0];
rz(-2.0013981) q[1];
sx q[1];
rz(-1.3277206) q[1];
sx q[1];
rz(0.16539703) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9388158) q[0];
sx q[0];
rz(-1.6949777) q[0];
sx q[0];
rz(-1.7569655) q[0];
rz(-pi) q[1];
rz(2.303316) q[2];
sx q[2];
rz(-2.8998313) q[2];
sx q[2];
rz(1.4571112) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4277028) q[1];
sx q[1];
rz(-2.1117003) q[1];
sx q[1];
rz(-1.9546024) q[1];
rz(1.1449279) q[3];
sx q[3];
rz(-0.59755675) q[3];
sx q[3];
rz(-1.4262517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.39607221) q[2];
sx q[2];
rz(-1.3362276) q[2];
sx q[2];
rz(2.9361434) q[2];
rz(2.0139587) q[3];
sx q[3];
rz(-1.1536359) q[3];
sx q[3];
rz(-1.8849461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66185343) q[0];
sx q[0];
rz(-2.2275708) q[0];
sx q[0];
rz(2.7868295) q[0];
rz(1.9873437) q[1];
sx q[1];
rz(-2.216979) q[1];
sx q[1];
rz(0.23194557) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7888768) q[0];
sx q[0];
rz(-2.0577601) q[0];
sx q[0];
rz(0.48077521) q[0];
rz(-pi) q[1];
rz(2.9734128) q[2];
sx q[2];
rz(-0.7025223) q[2];
sx q[2];
rz(-2.0481734) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0143118) q[1];
sx q[1];
rz(-0.31186549) q[1];
sx q[1];
rz(-0.77906268) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5902953) q[3];
sx q[3];
rz(-1.9392506) q[3];
sx q[3];
rz(1.3080477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0014687) q[2];
sx q[2];
rz(-1.3342369) q[2];
sx q[2];
rz(-1.9011964) q[2];
rz(-0.59605789) q[3];
sx q[3];
rz(-1.3052992) q[3];
sx q[3];
rz(1.3034472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.0444788) q[0];
sx q[0];
rz(-3.0713186) q[0];
sx q[0];
rz(-0.20275673) q[0];
rz(-0.98908201) q[1];
sx q[1];
rz(-1.6977856) q[1];
sx q[1];
rz(-0.99745497) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.049012262) q[0];
sx q[0];
rz(-0.870734) q[0];
sx q[0];
rz(-0.17605619) q[0];
rz(-pi) q[1];
rz(2.4100254) q[2];
sx q[2];
rz(-1.8533857) q[2];
sx q[2];
rz(-0.76991316) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.62085405) q[1];
sx q[1];
rz(-1.8100909) q[1];
sx q[1];
rz(-1.1362856) q[1];
x q[2];
rz(0.84647471) q[3];
sx q[3];
rz(-1.9123565) q[3];
sx q[3];
rz(-0.70111707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9138907) q[2];
sx q[2];
rz(-1.1662377) q[2];
sx q[2];
rz(2.690199) q[2];
rz(-0.40870062) q[3];
sx q[3];
rz(-1.6025851) q[3];
sx q[3];
rz(1.2020948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.30615) q[0];
sx q[0];
rz(-0.4040443) q[0];
sx q[0];
rz(0.62414449) q[0];
rz(-1.6250601) q[1];
sx q[1];
rz(-2.8846278) q[1];
sx q[1];
rz(0.61378941) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81543024) q[0];
sx q[0];
rz(-0.76061941) q[0];
sx q[0];
rz(1.8421696) q[0];
rz(-0.28366144) q[2];
sx q[2];
rz(-1.6322871) q[2];
sx q[2];
rz(-2.9279857) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2144199) q[1];
sx q[1];
rz(-1.5121636) q[1];
sx q[1];
rz(0.10417948) q[1];
rz(-pi) q[2];
rz(-1.8011439) q[3];
sx q[3];
rz(-2.6655572) q[3];
sx q[3];
rz(1.7649094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.43341407) q[2];
sx q[2];
rz(-2.2257979) q[2];
sx q[2];
rz(-1.1748574) q[2];
rz(-0.60837778) q[3];
sx q[3];
rz(-1.404168) q[3];
sx q[3];
rz(-1.7181989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2414918) q[0];
sx q[0];
rz(-1.2986203) q[0];
sx q[0];
rz(2.6532145) q[0];
rz(-1.5178559) q[1];
sx q[1];
rz(-1.7428215) q[1];
sx q[1];
rz(0.98446313) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0122089) q[0];
sx q[0];
rz(-1.0249656) q[0];
sx q[0];
rz(1.6145541) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1659053) q[2];
sx q[2];
rz(-2.3003909) q[2];
sx q[2];
rz(-1.0798432) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8967646) q[1];
sx q[1];
rz(-2.7762189) q[1];
sx q[1];
rz(0.016896292) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2582614) q[3];
sx q[3];
rz(-2.0274037) q[3];
sx q[3];
rz(1.6707591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4349334) q[2];
sx q[2];
rz(-1.443053) q[2];
sx q[2];
rz(2.6064176) q[2];
rz(-2.0914071) q[3];
sx q[3];
rz(-1.340056) q[3];
sx q[3];
rz(-1.0092658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.500279) q[0];
sx q[0];
rz(-2.2135493) q[0];
sx q[0];
rz(0.6859268) q[0];
rz(-2.7507239) q[1];
sx q[1];
rz(-2.1978244) q[1];
sx q[1];
rz(2.2156782) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6322964) q[0];
sx q[0];
rz(-2.8863781) q[0];
sx q[0];
rz(2.1584828) q[0];
rz(-pi) q[1];
x q[1];
rz(0.7943031) q[2];
sx q[2];
rz(-1.0216733) q[2];
sx q[2];
rz(1.4505475) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8402108) q[1];
sx q[1];
rz(-1.7173319) q[1];
sx q[1];
rz(2.0135897) q[1];
x q[2];
rz(1.6894475) q[3];
sx q[3];
rz(-1.7735529) q[3];
sx q[3];
rz(0.53393902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9459076) q[2];
sx q[2];
rz(-0.47161272) q[2];
sx q[2];
rz(-2.6055028) q[2];
rz(-2.7219971) q[3];
sx q[3];
rz(-2.5511238) q[3];
sx q[3];
rz(1.6528116) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0855899) q[0];
sx q[0];
rz(-0.35878006) q[0];
sx q[0];
rz(-0.39500239) q[0];
rz(1.5123873) q[1];
sx q[1];
rz(-1.2724266) q[1];
sx q[1];
rz(1.013247) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2241867) q[0];
sx q[0];
rz(-1.7099483) q[0];
sx q[0];
rz(1.805472) q[0];
rz(-3.0229438) q[2];
sx q[2];
rz(-0.35602202) q[2];
sx q[2];
rz(0.3686541) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0638748) q[1];
sx q[1];
rz(-0.43242726) q[1];
sx q[1];
rz(-2.9090911) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2447484) q[3];
sx q[3];
rz(-0.92771155) q[3];
sx q[3];
rz(-2.602592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.44832486) q[2];
sx q[2];
rz(-1.5193181) q[2];
sx q[2];
rz(-0.75941336) q[2];
rz(-1.365472) q[3];
sx q[3];
rz(-2.0335734) q[3];
sx q[3];
rz(0.56994462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7286745) q[0];
sx q[0];
rz(-1.8771794) q[0];
sx q[0];
rz(-1.9369453) q[0];
rz(2.5683174) q[1];
sx q[1];
rz(-1.234006) q[1];
sx q[1];
rz(-1.3201859) q[1];
rz(-1.7309932) q[2];
sx q[2];
rz(-2.6916531) q[2];
sx q[2];
rz(2.0315363) q[2];
rz(-1.6462973) q[3];
sx q[3];
rz(-0.69926881) q[3];
sx q[3];
rz(0.26509501) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
