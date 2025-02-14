OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.34831369) q[0];
sx q[0];
rz(-1.4077185) q[0];
sx q[0];
rz(-0.044535927) q[0];
rz(-1.1107923) q[1];
sx q[1];
rz(5.0587237) q[1];
sx q[1];
rz(8.6375477) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6605741) q[0];
sx q[0];
rz(-0.74587599) q[0];
sx q[0];
rz(-0.77156271) q[0];
x q[1];
rz(-1.6881053) q[2];
sx q[2];
rz(-0.32735887) q[2];
sx q[2];
rz(0.028520564) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.41809717) q[1];
sx q[1];
rz(-0.58291328) q[1];
sx q[1];
rz(-1.6969796) q[1];
rz(-pi) q[2];
x q[2];
rz(0.98712943) q[3];
sx q[3];
rz(-0.41582022) q[3];
sx q[3];
rz(-2.4087931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.74701509) q[2];
sx q[2];
rz(-2.0614823) q[2];
sx q[2];
rz(-0.75418312) q[2];
rz(1.4150103) q[3];
sx q[3];
rz(-2.6455046) q[3];
sx q[3];
rz(-0.32895857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0542145) q[0];
sx q[0];
rz(-0.35254756) q[0];
sx q[0];
rz(1.3341599) q[0];
rz(1.8502024) q[1];
sx q[1];
rz(-2.1033557) q[1];
sx q[1];
rz(0.87563595) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19568811) q[0];
sx q[0];
rz(-0.92600805) q[0];
sx q[0];
rz(-2.2553205) q[0];
x q[1];
rz(1.4562143) q[2];
sx q[2];
rz(-2.6422524) q[2];
sx q[2];
rz(0.90081604) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8671682) q[1];
sx q[1];
rz(-2.6058785) q[1];
sx q[1];
rz(-2.1650326) q[1];
rz(-pi) q[2];
x q[2];
rz(0.29107901) q[3];
sx q[3];
rz(-1.7699695) q[3];
sx q[3];
rz(1.8739623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4448173) q[2];
sx q[2];
rz(-2.4397218) q[2];
sx q[2];
rz(-0.87265054) q[2];
rz(0.78898346) q[3];
sx q[3];
rz(-0.9011457) q[3];
sx q[3];
rz(0.22855973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92662421) q[0];
sx q[0];
rz(-1.7890395) q[0];
sx q[0];
rz(-3.1335926) q[0];
rz(2.3232715) q[1];
sx q[1];
rz(-0.74940959) q[1];
sx q[1];
rz(1.9047033) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43046666) q[0];
sx q[0];
rz(-0.9014117) q[0];
sx q[0];
rz(-0.41394625) q[0];
rz(-0.90571442) q[2];
sx q[2];
rz(-0.53230282) q[2];
sx q[2];
rz(-1.5168845) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.8244918) q[1];
sx q[1];
rz(-0.6380837) q[1];
sx q[1];
rz(-0.53148766) q[1];
x q[2];
rz(2.6872356) q[3];
sx q[3];
rz(-0.26232346) q[3];
sx q[3];
rz(1.232805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7311953) q[2];
sx q[2];
rz(-3.0072913) q[2];
sx q[2];
rz(2.2006939) q[2];
rz(-1.1040374) q[3];
sx q[3];
rz(-1.7887812) q[3];
sx q[3];
rz(1.997939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0333772) q[0];
sx q[0];
rz(-0.5431076) q[0];
sx q[0];
rz(0.23750842) q[0];
rz(-2.4820651) q[1];
sx q[1];
rz(-2.567465) q[1];
sx q[1];
rz(3.1245756) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4537246) q[0];
sx q[0];
rz(-0.17565082) q[0];
sx q[0];
rz(-1.2916628) q[0];
x q[1];
rz(-0.41827664) q[2];
sx q[2];
rz(-2.661663) q[2];
sx q[2];
rz(1.0470225) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2260277) q[1];
sx q[1];
rz(-1.6751958) q[1];
sx q[1];
rz(0.56609367) q[1];
rz(-pi) q[2];
x q[2];
rz(0.7155101) q[3];
sx q[3];
rz(-2.1984716) q[3];
sx q[3];
rz(2.024141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.90618769) q[2];
sx q[2];
rz(-1.2013288) q[2];
sx q[2];
rz(-2.2271633) q[2];
rz(0.36214456) q[3];
sx q[3];
rz(-2.3320964) q[3];
sx q[3];
rz(0.58524281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2528766) q[0];
sx q[0];
rz(-1.7676366) q[0];
sx q[0];
rz(-0.1121029) q[0];
rz(-1.0653227) q[1];
sx q[1];
rz(-2.0707097) q[1];
sx q[1];
rz(1.0692474) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1493268) q[0];
sx q[0];
rz(-2.7130824) q[0];
sx q[0];
rz(0.70121588) q[0];
rz(-pi) q[1];
x q[1];
rz(0.71938558) q[2];
sx q[2];
rz(-0.66443887) q[2];
sx q[2];
rz(-1.4321355) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.18373173) q[1];
sx q[1];
rz(-1.2535447) q[1];
sx q[1];
rz(-0.31362335) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0085722) q[3];
sx q[3];
rz(-2.1658796) q[3];
sx q[3];
rz(-1.8436335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7142882) q[2];
sx q[2];
rz(-2.4236743) q[2];
sx q[2];
rz(2.60738) q[2];
rz(2.838375) q[3];
sx q[3];
rz(-1.5568308) q[3];
sx q[3];
rz(1.6259954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41820207) q[0];
sx q[0];
rz(-1.284282) q[0];
sx q[0];
rz(-0.89163017) q[0];
rz(1.760969) q[1];
sx q[1];
rz(-1.9386407) q[1];
sx q[1];
rz(-2.2183529) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1311125) q[0];
sx q[0];
rz(-1.4725871) q[0];
sx q[0];
rz(-1.6232071) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60439749) q[2];
sx q[2];
rz(-0.45995074) q[2];
sx q[2];
rz(2.9876171) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1404851) q[1];
sx q[1];
rz(-1.342089) q[1];
sx q[1];
rz(-1.2910991) q[1];
rz(2.4999077) q[3];
sx q[3];
rz(-1.2417792) q[3];
sx q[3];
rz(-1.6365479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.096006958) q[2];
sx q[2];
rz(-0.43232375) q[2];
sx q[2];
rz(1.2621137) q[2];
rz(-1.1612085) q[3];
sx q[3];
rz(-1.592417) q[3];
sx q[3];
rz(-1.2459292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68814174) q[0];
sx q[0];
rz(-2.3095486) q[0];
sx q[0];
rz(-1.6081109) q[0];
rz(-0.43117943) q[1];
sx q[1];
rz(-0.9175514) q[1];
sx q[1];
rz(1.4088438) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.726175) q[0];
sx q[0];
rz(-1.8573995) q[0];
sx q[0];
rz(-0.64985259) q[0];
rz(-pi) q[1];
rz(-0.088690745) q[2];
sx q[2];
rz(-1.9541395) q[2];
sx q[2];
rz(0.92244043) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.81804336) q[1];
sx q[1];
rz(-2.9406639) q[1];
sx q[1];
rz(2.4646321) q[1];
x q[2];
rz(1.3471782) q[3];
sx q[3];
rz(-2.6287492) q[3];
sx q[3];
rz(-1.4033069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2048753) q[2];
sx q[2];
rz(-1.1254346) q[2];
sx q[2];
rz(-1.3196866) q[2];
rz(0.034505757) q[3];
sx q[3];
rz(-1.6104108) q[3];
sx q[3];
rz(-1.7721133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40783229) q[0];
sx q[0];
rz(-2.5875081) q[0];
sx q[0];
rz(-1.2169417) q[0];
rz(-0.92562428) q[1];
sx q[1];
rz(-1.3652912) q[1];
sx q[1];
rz(1.9409174) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.901163) q[0];
sx q[0];
rz(-1.0679185) q[0];
sx q[0];
rz(0.85977413) q[0];
x q[1];
rz(-1.9782135) q[2];
sx q[2];
rz(-1.0098741) q[2];
sx q[2];
rz(-1.7307868) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0266704) q[1];
sx q[1];
rz(-1.5296796) q[1];
sx q[1];
rz(2.369472) q[1];
rz(0.36603277) q[3];
sx q[3];
rz(-1.9683553) q[3];
sx q[3];
rz(-1.1432858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.22244464) q[2];
sx q[2];
rz(-1.5347975) q[2];
sx q[2];
rz(2.0929125) q[2];
rz(2.9564296) q[3];
sx q[3];
rz(-1.984963) q[3];
sx q[3];
rz(0.10704253) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1222526) q[0];
sx q[0];
rz(-0.66146079) q[0];
sx q[0];
rz(0.81073236) q[0];
rz(0.73075378) q[1];
sx q[1];
rz(-2.0850875) q[1];
sx q[1];
rz(0.96053851) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7030684) q[0];
sx q[0];
rz(-1.5674066) q[0];
sx q[0];
rz(-1.5698213) q[0];
rz(-pi) q[1];
rz(2.075718) q[2];
sx q[2];
rz(-1.9532579) q[2];
sx q[2];
rz(-1.3825716) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9593138) q[1];
sx q[1];
rz(-2.5393705) q[1];
sx q[1];
rz(-2.0456362) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1436362) q[3];
sx q[3];
rz(-1.33516) q[3];
sx q[3];
rz(2.4164661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1824128) q[2];
sx q[2];
rz(-1.7824087) q[2];
sx q[2];
rz(1.5394999) q[2];
rz(1.0473853) q[3];
sx q[3];
rz(-1.3771219) q[3];
sx q[3];
rz(1.876095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78028107) q[0];
sx q[0];
rz(-0.86532101) q[0];
sx q[0];
rz(-2.5541232) q[0];
rz(2.7777708) q[1];
sx q[1];
rz(-1.1452585) q[1];
sx q[1];
rz(-0.90715539) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.912187) q[0];
sx q[0];
rz(-2.7813781) q[0];
sx q[0];
rz(0.2304669) q[0];
x q[1];
rz(-1.9586708) q[2];
sx q[2];
rz(-1.1287612) q[2];
sx q[2];
rz(1.7427874) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.75395155) q[1];
sx q[1];
rz(-2.6389737) q[1];
sx q[1];
rz(-2.6941264) q[1];
rz(1.9349707) q[3];
sx q[3];
rz(-1.4910021) q[3];
sx q[3];
rz(-1.4353115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5017447) q[2];
sx q[2];
rz(-0.088968337) q[2];
sx q[2];
rz(-0.42195827) q[2];
rz(-0.0095327775) q[3];
sx q[3];
rz(-1.5671174) q[3];
sx q[3];
rz(-2.6154521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52787732) q[0];
sx q[0];
rz(-1.8462702) q[0];
sx q[0];
rz(0.12311002) q[0];
rz(0.70115024) q[1];
sx q[1];
rz(-1.9155365) q[1];
sx q[1];
rz(1.6446) q[1];
rz(0.43177615) q[2];
sx q[2];
rz(-2.3974621) q[2];
sx q[2];
rz(1.5324788) q[2];
rz(2.8641239) q[3];
sx q[3];
rz(-1.1744842) q[3];
sx q[3];
rz(-0.13406772) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
