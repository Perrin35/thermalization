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
rz(-0.31565491) q[0];
rz(1.1701801) q[1];
sx q[1];
rz(-0.38771954) q[1];
sx q[1];
rz(0.60400909) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9785288) q[0];
sx q[0];
rz(-1.8730358) q[0];
sx q[0];
rz(3.0885484) q[0];
rz(-2.1748105) q[2];
sx q[2];
rz(-0.76672115) q[2];
sx q[2];
rz(1.3441835) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4166222) q[1];
sx q[1];
rz(-1.23471) q[1];
sx q[1];
rz(1.027847) q[1];
rz(-pi) q[2];
rz(3.0331361) q[3];
sx q[3];
rz(-1.2634957) q[3];
sx q[3];
rz(1.1600329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9903367) q[2];
sx q[2];
rz(-1.9467111) q[2];
sx q[2];
rz(-1.3585496) q[2];
rz(-1.5938866) q[3];
sx q[3];
rz(-2.2113776) q[3];
sx q[3];
rz(-1.5096629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5556521) q[0];
sx q[0];
rz(-0.63783115) q[0];
sx q[0];
rz(1.1974539) q[0];
rz(-0.50152957) q[1];
sx q[1];
rz(-1.5781559) q[1];
sx q[1];
rz(-1.4362358) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.013151) q[0];
sx q[0];
rz(-1.3396946) q[0];
sx q[0];
rz(-0.9856772) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6701489) q[2];
sx q[2];
rz(-2.2713285) q[2];
sx q[2];
rz(-1.5347028) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.20138559) q[1];
sx q[1];
rz(-1.4680982) q[1];
sx q[1];
rz(0.19725712) q[1];
rz(-pi) q[2];
x q[2];
rz(0.40879876) q[3];
sx q[3];
rz(-2.650947) q[3];
sx q[3];
rz(1.7237494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.88840914) q[2];
sx q[2];
rz(-1.2359572) q[2];
sx q[2];
rz(-3.1186228) q[2];
rz(-2.2394771) q[3];
sx q[3];
rz(-1.2227367) q[3];
sx q[3];
rz(-2.7149916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70812923) q[0];
sx q[0];
rz(-1.3490278) q[0];
sx q[0];
rz(0.9915114) q[0];
rz(-2.1082711) q[1];
sx q[1];
rz(-0.28305498) q[1];
sx q[1];
rz(1.906377) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2184705) q[0];
sx q[0];
rz(-1.8184796) q[0];
sx q[0];
rz(0.60486233) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5052983) q[2];
sx q[2];
rz(-2.4264196) q[2];
sx q[2];
rz(1.6805122) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8950303) q[1];
sx q[1];
rz(-2.3091279) q[1];
sx q[1];
rz(-2.8042996) q[1];
rz(0.54872437) q[3];
sx q[3];
rz(-1.0503891) q[3];
sx q[3];
rz(1.5279552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.027355) q[2];
sx q[2];
rz(-2.3894252) q[2];
sx q[2];
rz(2.1503964) q[2];
rz(-1.2973971) q[3];
sx q[3];
rz(-1.2605366) q[3];
sx q[3];
rz(1.7245002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39795136) q[0];
sx q[0];
rz(-0.93248168) q[0];
sx q[0];
rz(0.81746307) q[0];
rz(2.815333) q[1];
sx q[1];
rz(-1.9605109) q[1];
sx q[1];
rz(0.78525966) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2366339) q[0];
sx q[0];
rz(-1.8371305) q[0];
sx q[0];
rz(1.2575218) q[0];
rz(-2.4982959) q[2];
sx q[2];
rz(-1.4125106) q[2];
sx q[2];
rz(2.9093551) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1568415) q[1];
sx q[1];
rz(-1.9923216) q[1];
sx q[1];
rz(-2.0578029) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.35492949) q[3];
sx q[3];
rz(-1.4961492) q[3];
sx q[3];
rz(-1.1433034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0205445) q[2];
sx q[2];
rz(-0.64457568) q[2];
sx q[2];
rz(2.5784967) q[2];
rz(-0.8935039) q[3];
sx q[3];
rz(-1.9681135) q[3];
sx q[3];
rz(-1.2347429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23984443) q[0];
sx q[0];
rz(-2.8185066) q[0];
sx q[0];
rz(-1.595994) q[0];
rz(-1.0218989) q[1];
sx q[1];
rz(-1.1232168) q[1];
sx q[1];
rz(0.18128577) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5820532) q[0];
sx q[0];
rz(-1.7632369) q[0];
sx q[0];
rz(-2.0700702) q[0];
x q[1];
rz(-1.6185804) q[2];
sx q[2];
rz(-1.0956229) q[2];
sx q[2];
rz(-2.938624) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.1039155) q[1];
sx q[1];
rz(-1.5427151) q[1];
sx q[1];
rz(-1.3485391) q[1];
rz(-pi) q[2];
rz(0.44293176) q[3];
sx q[3];
rz(-1.4287523) q[3];
sx q[3];
rz(2.1403109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.01866092) q[2];
sx q[2];
rz(-2.4666726) q[2];
sx q[2];
rz(-1.2459416) q[2];
rz(2.5490226) q[3];
sx q[3];
rz(-1.4453245) q[3];
sx q[3];
rz(-0.84111253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6084006) q[0];
sx q[0];
rz(-2.1493981) q[0];
sx q[0];
rz(-0.30884185) q[0];
rz(0.32288512) q[1];
sx q[1];
rz(-2.7643118) q[1];
sx q[1];
rz(0.13993941) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9788743) q[0];
sx q[0];
rz(-2.1669183) q[0];
sx q[0];
rz(0.30124248) q[0];
rz(-0.094801589) q[2];
sx q[2];
rz(-2.0249512) q[2];
sx q[2];
rz(-2.4081782) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7227712) q[1];
sx q[1];
rz(-0.79192894) q[1];
sx q[1];
rz(-2.0008816) q[1];
rz(-pi) q[2];
x q[2];
rz(2.841137) q[3];
sx q[3];
rz(-0.93378919) q[3];
sx q[3];
rz(-0.14618044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4796925) q[2];
sx q[2];
rz(-2.4883344) q[2];
sx q[2];
rz(-2.8996186) q[2];
rz(0.74609977) q[3];
sx q[3];
rz(-1.7753764) q[3];
sx q[3];
rz(1.7297176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11874966) q[0];
sx q[0];
rz(-2.0976837) q[0];
sx q[0];
rz(-1.2710849) q[0];
rz(2.5785043) q[1];
sx q[1];
rz(-0.92910281) q[1];
sx q[1];
rz(-1.44106) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4615501) q[0];
sx q[0];
rz(-2.474804) q[0];
sx q[0];
rz(-2.0186485) q[0];
x q[1];
rz(-0.38132847) q[2];
sx q[2];
rz(-1.4421808) q[2];
sx q[2];
rz(2.6113457) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0771516) q[1];
sx q[1];
rz(-1.3485326) q[1];
sx q[1];
rz(0.42776107) q[1];
x q[2];
rz(-0.7251803) q[3];
sx q[3];
rz(-1.7340875) q[3];
sx q[3];
rz(1.1226599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9219804) q[2];
sx q[2];
rz(-1.8461123) q[2];
sx q[2];
rz(-1.8249576) q[2];
rz(-0.5611788) q[3];
sx q[3];
rz(-0.96446529) q[3];
sx q[3];
rz(-1.5285899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.104326) q[0];
sx q[0];
rz(-0.21793652) q[0];
sx q[0];
rz(2.2547146) q[0];
rz(-0.2001702) q[1];
sx q[1];
rz(-1.6693516) q[1];
sx q[1];
rz(1.0221457) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5258785) q[0];
sx q[0];
rz(-1.4396884) q[0];
sx q[0];
rz(3.1159621) q[0];
x q[1];
rz(-0.81424197) q[2];
sx q[2];
rz(-2.1541284) q[2];
sx q[2];
rz(-0.63177201) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.31769842) q[1];
sx q[1];
rz(-2.3168457) q[1];
sx q[1];
rz(-3.1391678) q[1];
rz(2.5660146) q[3];
sx q[3];
rz(-1.8691571) q[3];
sx q[3];
rz(-2.9090867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.48978051) q[2];
sx q[2];
rz(-0.45280364) q[2];
sx q[2];
rz(2.32302) q[2];
rz(0.98322785) q[3];
sx q[3];
rz(-0.95663095) q[3];
sx q[3];
rz(2.6035068) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1061358) q[0];
sx q[0];
rz(-1.1470733) q[0];
sx q[0];
rz(-1.5552833) q[0];
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
rz(-1.9826384) q[0];
sx q[0];
rz(-2.1110015) q[0];
sx q[0];
rz(-1.0544974) q[0];
rz(0.92918877) q[2];
sx q[2];
rz(-0.25114533) q[2];
sx q[2];
rz(-0.26640688) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6482249) q[1];
sx q[1];
rz(-0.43637143) q[1];
sx q[1];
rz(-0.079597537) q[1];
rz(2.5303477) q[3];
sx q[3];
rz(-1.2204303) q[3];
sx q[3];
rz(-0.65920748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.67184225) q[2];
sx q[2];
rz(-1.2776813) q[2];
sx q[2];
rz(2.3308241) q[2];
rz(-1.9226711) q[3];
sx q[3];
rz(-2.6007077) q[3];
sx q[3];
rz(-1.0651275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78275457) q[0];
sx q[0];
rz(-2.8420119) q[0];
sx q[0];
rz(1.7175571) q[0];
rz(0.77761039) q[1];
sx q[1];
rz(-2.168226) q[1];
sx q[1];
rz(-2.3861859) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7292405) q[0];
sx q[0];
rz(-1.5753813) q[0];
sx q[0];
rz(1.3608749) q[0];
rz(-pi) q[1];
rz(0.39647409) q[2];
sx q[2];
rz(-1.586953) q[2];
sx q[2];
rz(2.0153869) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3302119) q[1];
sx q[1];
rz(-2.5054058) q[1];
sx q[1];
rz(-0.61417555) q[1];
rz(-pi) q[2];
rz(1.2437263) q[3];
sx q[3];
rz(-2.2199515) q[3];
sx q[3];
rz(-2.7130733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5767673) q[2];
sx q[2];
rz(-0.29594031) q[2];
sx q[2];
rz(-0.15288615) q[2];
rz(1.8704002) q[3];
sx q[3];
rz(-1.8891687) q[3];
sx q[3];
rz(-2.474474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2818903) q[0];
sx q[0];
rz(-0.49260456) q[0];
sx q[0];
rz(2.3717666) q[0];
rz(-1.9181171) q[1];
sx q[1];
rz(-1.9203095) q[1];
sx q[1];
rz(-0.65345678) q[1];
rz(-0.6003004) q[2];
sx q[2];
rz(-0.68927286) q[2];
sx q[2];
rz(2.2130028) q[2];
rz(-2.2822904) q[3];
sx q[3];
rz(-0.51325428) q[3];
sx q[3];
rz(-0.2556066) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
