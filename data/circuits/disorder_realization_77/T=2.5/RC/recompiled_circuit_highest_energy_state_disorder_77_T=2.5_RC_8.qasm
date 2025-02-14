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
rz(2.6440808) q[0];
sx q[0];
rz(-1.3389791) q[0];
sx q[0];
rz(2.8259377) q[0];
rz(1.1701801) q[1];
sx q[1];
rz(-0.38771954) q[1];
sx q[1];
rz(-2.5375836) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7496628) q[0];
sx q[0];
rz(-1.6214341) q[0];
sx q[0];
rz(-1.2681566) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6409615) q[2];
sx q[2];
rz(-2.1785469) q[2];
sx q[2];
rz(-2.1082102) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7249704) q[1];
sx q[1];
rz(-1.23471) q[1];
sx q[1];
rz(-2.1137456) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8798001) q[3];
sx q[3];
rz(-1.4674392) q[3];
sx q[3];
rz(2.7637533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9903367) q[2];
sx q[2];
rz(-1.9467111) q[2];
sx q[2];
rz(-1.3585496) q[2];
rz(1.547706) q[3];
sx q[3];
rz(-0.93021506) q[3];
sx q[3];
rz(1.5096629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(-0.58594054) q[0];
sx q[0];
rz(-2.5037615) q[0];
sx q[0];
rz(1.9441388) q[0];
rz(-0.50152957) q[1];
sx q[1];
rz(-1.5634368) q[1];
sx q[1];
rz(1.4362358) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89051188) q[0];
sx q[0];
rz(-2.5174826) q[0];
sx q[0];
rz(1.1680383) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6701489) q[2];
sx q[2];
rz(-0.87026419) q[2];
sx q[2];
rz(-1.6068899) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9402071) q[1];
sx q[1];
rz(-1.4680982) q[1];
sx q[1];
rz(0.19725712) q[1];
x q[2];
rz(-0.45577502) q[3];
sx q[3];
rz(-1.7592128) q[3];
sx q[3];
rz(2.6236603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2531835) q[2];
sx q[2];
rz(-1.2359572) q[2];
sx q[2];
rz(-3.1186228) q[2];
rz(2.2394771) q[3];
sx q[3];
rz(-1.918856) q[3];
sx q[3];
rz(-2.7149916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-2.4334634) q[0];
sx q[0];
rz(-1.3490278) q[0];
sx q[0];
rz(2.1500812) q[0];
rz(-2.1082711) q[1];
sx q[1];
rz(-0.28305498) q[1];
sx q[1];
rz(1.906377) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9231222) q[0];
sx q[0];
rz(-1.8184796) q[0];
sx q[0];
rz(2.5367303) q[0];
rz(-pi) q[1];
rz(2.5318145) q[2];
sx q[2];
rz(-1.9710566) q[2];
sx q[2];
rz(-0.61851172) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.72690645) q[1];
sx q[1];
rz(-2.3432846) q[1];
sx q[1];
rz(1.2220135) q[1];
rz(2.1622873) q[3];
sx q[3];
rz(-1.1011754) q[3];
sx q[3];
rz(0.33794935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.11423763) q[2];
sx q[2];
rz(-0.75216746) q[2];
sx q[2];
rz(0.99119622) q[2];
rz(1.2973971) q[3];
sx q[3];
rz(-1.2605366) q[3];
sx q[3];
rz(1.4170925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7436413) q[0];
sx q[0];
rz(-0.93248168) q[0];
sx q[0];
rz(2.3241296) q[0];
rz(2.815333) q[1];
sx q[1];
rz(-1.9605109) q[1];
sx q[1];
rz(-2.356333) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8924849) q[0];
sx q[0];
rz(-1.2689225) q[0];
sx q[0];
rz(0.27927804) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7677069) q[2];
sx q[2];
rz(-2.2047538) q[2];
sx q[2];
rz(-1.2209354) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2438812) q[1];
sx q[1];
rz(-0.63279018) q[1];
sx q[1];
rz(-0.80676078) q[1];
rz(-pi) q[2];
rz(2.9296285) q[3];
sx q[3];
rz(-2.779224) q[3];
sx q[3];
rz(2.9126698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0205445) q[2];
sx q[2];
rz(-2.497017) q[2];
sx q[2];
rz(-0.56309593) q[2];
rz(-2.2480887) q[3];
sx q[3];
rz(-1.9681135) q[3];
sx q[3];
rz(1.2347429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9017482) q[0];
sx q[0];
rz(-0.32308602) q[0];
sx q[0];
rz(1.5455986) q[0];
rz(2.1196938) q[1];
sx q[1];
rz(-1.1232168) q[1];
sx q[1];
rz(0.18128577) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.793042) q[0];
sx q[0];
rz(-2.6094545) q[0];
sx q[0];
rz(1.1842968) q[0];
rz(-pi) q[1];
rz(-0.092575707) q[2];
sx q[2];
rz(-2.6642054) q[2];
sx q[2];
rz(0.30711781) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5904894) q[1];
sx q[1];
rz(-0.22399513) q[1];
sx q[1];
rz(1.444054) q[1];
x q[2];
rz(-1.4138172) q[3];
sx q[3];
rz(-1.1326367) q[3];
sx q[3];
rz(2.639132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.01866092) q[2];
sx q[2];
rz(-0.67492008) q[2];
sx q[2];
rz(-1.8956511) q[2];
rz(-2.5490226) q[3];
sx q[3];
rz(-1.4453245) q[3];
sx q[3];
rz(-2.3004801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5331921) q[0];
sx q[0];
rz(-2.1493981) q[0];
sx q[0];
rz(2.8327508) q[0];
rz(-0.32288512) q[1];
sx q[1];
rz(-0.37728089) q[1];
sx q[1];
rz(0.13993941) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16271834) q[0];
sx q[0];
rz(-2.1669183) q[0];
sx q[0];
rz(-2.8403502) q[0];
rz(-pi) q[1];
rz(1.1148648) q[2];
sx q[2];
rz(-1.6559634) q[2];
sx q[2];
rz(-2.2625201) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4637101) q[1];
sx q[1];
rz(-1.8720798) q[1];
sx q[1];
rz(0.82656411) q[1];
rz(1.9512964) q[3];
sx q[3];
rz(-2.4463013) q[3];
sx q[3];
rz(-0.33405802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4796925) q[2];
sx q[2];
rz(-2.4883344) q[2];
sx q[2];
rz(2.8996186) q[2];
rz(2.3954929) q[3];
sx q[3];
rz(-1.7753764) q[3];
sx q[3];
rz(1.411875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11874966) q[0];
sx q[0];
rz(-1.043909) q[0];
sx q[0];
rz(-1.2710849) q[0];
rz(0.56308833) q[1];
sx q[1];
rz(-2.2124898) q[1];
sx q[1];
rz(1.7005327) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4615501) q[0];
sx q[0];
rz(-0.6667887) q[0];
sx q[0];
rz(2.0186485) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7602642) q[2];
sx q[2];
rz(-1.4421808) q[2];
sx q[2];
rz(-2.6113457) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.043277117) q[1];
sx q[1];
rz(-0.47887427) q[1];
sx q[1];
rz(-0.4988341) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7874927) q[3];
sx q[3];
rz(-0.85734493) q[3];
sx q[3];
rz(2.5503623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9219804) q[2];
sx q[2];
rz(-1.2954804) q[2];
sx q[2];
rz(1.316635) q[2];
rz(-0.5611788) q[3];
sx q[3];
rz(-0.96446529) q[3];
sx q[3];
rz(-1.5285899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.104326) q[0];
sx q[0];
rz(-2.9236561) q[0];
sx q[0];
rz(-2.2547146) q[0];
rz(-2.9414224) q[1];
sx q[1];
rz(-1.472241) q[1];
sx q[1];
rz(1.0221457) q[1];
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
rz(-2.3273507) q[2];
sx q[2];
rz(-0.98746429) q[2];
sx q[2];
rz(2.5098206) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8868489) q[1];
sx q[1];
rz(-1.5690156) q[1];
sx q[1];
rz(0.82474553) q[1];
rz(-pi) q[2];
rz(0.57557801) q[3];
sx q[3];
rz(-1.8691571) q[3];
sx q[3];
rz(-0.23250599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.48978051) q[2];
sx q[2];
rz(-2.688789) q[2];
sx q[2];
rz(2.32302) q[2];
rz(2.1583648) q[3];
sx q[3];
rz(-0.95663095) q[3];
sx q[3];
rz(-2.6035068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.035456903) q[0];
sx q[0];
rz(-1.1470733) q[0];
sx q[0];
rz(1.5863093) q[0];
rz(-0.69681329) q[1];
sx q[1];
rz(-0.15334829) q[1];
sx q[1];
rz(0.78479016) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12778388) q[0];
sx q[0];
rz(-1.133636) q[0];
sx q[0];
rz(0.60370914) q[0];
rz(0.15235849) q[2];
sx q[2];
rz(-1.7712284) q[2];
sx q[2];
rz(-0.39055017) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6482249) q[1];
sx q[1];
rz(-0.43637143) q[1];
sx q[1];
rz(-3.0619951) q[1];
x q[2];
rz(-2.574563) q[3];
sx q[3];
rz(-2.4483557) q[3];
sx q[3];
rz(-2.6854613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.67184225) q[2];
sx q[2];
rz(-1.2776813) q[2];
sx q[2];
rz(-0.81076852) q[2];
rz(-1.9226711) q[3];
sx q[3];
rz(-0.54088497) q[3];
sx q[3];
rz(1.0651275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78275457) q[0];
sx q[0];
rz(-2.8420119) q[0];
sx q[0];
rz(1.7175571) q[0];
rz(2.3639823) q[1];
sx q[1];
rz(-2.168226) q[1];
sx q[1];
rz(-0.75540677) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1799602) q[0];
sx q[0];
rz(-2.9316219) q[0];
sx q[0];
rz(-1.5487973) q[0];
x q[1];
rz(0.39647409) q[2];
sx q[2];
rz(-1.5546397) q[2];
sx q[2];
rz(-2.0153869) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2753693) q[1];
sx q[1];
rz(-1.2213372) q[1];
sx q[1];
rz(2.5984955) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4660268) q[3];
sx q[3];
rz(-1.3119952) q[3];
sx q[3];
rz(1.7970366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5648254) q[2];
sx q[2];
rz(-0.29594031) q[2];
sx q[2];
rz(2.9887065) q[2];
rz(-1.2711924) q[3];
sx q[3];
rz(-1.252424) q[3];
sx q[3];
rz(-0.66711867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2818903) q[0];
sx q[0];
rz(-0.49260456) q[0];
sx q[0];
rz(2.3717666) q[0];
rz(-1.2234756) q[1];
sx q[1];
rz(-1.2212831) q[1];
sx q[1];
rz(2.4881359) q[1];
rz(-0.59719795) q[2];
sx q[2];
rz(-1.2033249) q[2];
sx q[2];
rz(-2.0133599) q[2];
rz(0.35265233) q[3];
sx q[3];
rz(-1.1897539) q[3];
sx q[3];
rz(-2.6170058) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
