OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.45286173) q[0];
sx q[0];
rz(2.9714669) q[0];
sx q[0];
rz(10.210769) q[0];
rz(-2.535948) q[1];
sx q[1];
rz(-0.65094596) q[1];
sx q[1];
rz(0.63408607) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85672985) q[0];
sx q[0];
rz(-2.3861109) q[0];
sx q[0];
rz(1.612624) q[0];
x q[1];
rz(-1.6526821) q[2];
sx q[2];
rz(-2.4856644) q[2];
sx q[2];
rz(-1.8413078) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4418728) q[1];
sx q[1];
rz(-0.77077121) q[1];
sx q[1];
rz(2.2295879) q[1];
rz(0.4001873) q[3];
sx q[3];
rz(-1.549198) q[3];
sx q[3];
rz(0.17822972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9156076) q[2];
sx q[2];
rz(-0.51708022) q[2];
sx q[2];
rz(1.263164) q[2];
rz(-1.2849215) q[3];
sx q[3];
rz(-1.6804755) q[3];
sx q[3];
rz(0.29299709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9830575) q[0];
sx q[0];
rz(-1.2936658) q[0];
sx q[0];
rz(-0.43757004) q[0];
rz(0.63105398) q[1];
sx q[1];
rz(-2.8129306) q[1];
sx q[1];
rz(-0.22110573) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3723345) q[0];
sx q[0];
rz(-0.49494574) q[0];
sx q[0];
rz(1.1713722) q[0];
rz(1.919585) q[2];
sx q[2];
rz(-1.4663327) q[2];
sx q[2];
rz(1.2806569) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4202538) q[1];
sx q[1];
rz(-1.6920648) q[1];
sx q[1];
rz(0.46055693) q[1];
x q[2];
rz(-0.50237327) q[3];
sx q[3];
rz(-1.6564192) q[3];
sx q[3];
rz(-1.3115713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.21800403) q[2];
sx q[2];
rz(-1.709334) q[2];
sx q[2];
rz(-0.58829266) q[2];
rz(2.6925987) q[3];
sx q[3];
rz(-2.7136927) q[3];
sx q[3];
rz(-2.0935521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5730729) q[0];
sx q[0];
rz(-0.097671106) q[0];
sx q[0];
rz(1.0004689) q[0];
rz(0.72552848) q[1];
sx q[1];
rz(-2.0670481) q[1];
sx q[1];
rz(-2.3838938) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.476186) q[0];
sx q[0];
rz(-1.2520257) q[0];
sx q[0];
rz(0.55253367) q[0];
rz(-pi) q[1];
rz(2.9341142) q[2];
sx q[2];
rz(-1.7277272) q[2];
sx q[2];
rz(0.66397882) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.11848395) q[1];
sx q[1];
rz(-0.99220905) q[1];
sx q[1];
rz(2.5793377) q[1];
rz(-1.682231) q[3];
sx q[3];
rz(-2.206344) q[3];
sx q[3];
rz(-3.0524658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9329325) q[2];
sx q[2];
rz(-0.22298403) q[2];
sx q[2];
rz(2.3051252) q[2];
rz(-1.4423192) q[3];
sx q[3];
rz(-1.4740372) q[3];
sx q[3];
rz(-2.9377655) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.02012415) q[0];
sx q[0];
rz(-3.0614873) q[0];
sx q[0];
rz(-0.51112038) q[0];
rz(-0.077443667) q[1];
sx q[1];
rz(-2.7192392) q[1];
sx q[1];
rz(-1.8130594) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1140095) q[0];
sx q[0];
rz(-0.264835) q[0];
sx q[0];
rz(0.94075216) q[0];
rz(-pi) q[1];
rz(0.51283299) q[2];
sx q[2];
rz(-2.3094258) q[2];
sx q[2];
rz(1.2733449) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1991785) q[1];
sx q[1];
rz(-1.5713099) q[1];
sx q[1];
rz(-1.3039939) q[1];
rz(-pi) q[2];
rz(-2.02416) q[3];
sx q[3];
rz(-2.8957267) q[3];
sx q[3];
rz(-0.21594957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0478583) q[2];
sx q[2];
rz(-1.9231984) q[2];
sx q[2];
rz(2.034534) q[2];
rz(2.6691061) q[3];
sx q[3];
rz(-1.6146336) q[3];
sx q[3];
rz(-2.2842177) q[3];
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
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040314019) q[0];
sx q[0];
rz(-2.2450228) q[0];
sx q[0];
rz(-0.3381981) q[0];
rz(-1.2942554) q[1];
sx q[1];
rz(-1.5816403) q[1];
sx q[1];
rz(-2.6370874) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.757526) q[0];
sx q[0];
rz(-2.5528918) q[0];
sx q[0];
rz(1.9996044) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1930824) q[2];
sx q[2];
rz(-1.1337122) q[2];
sx q[2];
rz(1.2010241) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5711172) q[1];
sx q[1];
rz(-1.332453) q[1];
sx q[1];
rz(-1.2586602) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3584064) q[3];
sx q[3];
rz(-0.83055701) q[3];
sx q[3];
rz(0.98379788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.111104) q[2];
sx q[2];
rz(-2.3031074) q[2];
sx q[2];
rz(1.6476691) q[2];
rz(1.4533639) q[3];
sx q[3];
rz(-0.93795347) q[3];
sx q[3];
rz(1.3822314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21022739) q[0];
sx q[0];
rz(-0.87451044) q[0];
sx q[0];
rz(-2.0507623) q[0];
rz(0.53972721) q[1];
sx q[1];
rz(-1.9878186) q[1];
sx q[1];
rz(-0.18879034) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15121962) q[0];
sx q[0];
rz(-1.2114721) q[0];
sx q[0];
rz(-1.3034526) q[0];
rz(-1.9300869) q[2];
sx q[2];
rz(-2.852716) q[2];
sx q[2];
rz(-1.2668244) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.530045) q[1];
sx q[1];
rz(-1.5049184) q[1];
sx q[1];
rz(-0.24286119) q[1];
rz(2.4981899) q[3];
sx q[3];
rz(-2.8426369) q[3];
sx q[3];
rz(-1.9770196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.171689) q[2];
sx q[2];
rz(-1.0287372) q[2];
sx q[2];
rz(-1.3151273) q[2];
rz(-1.5054437) q[3];
sx q[3];
rz(-0.35741487) q[3];
sx q[3];
rz(1.858254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0514907) q[0];
sx q[0];
rz(-2.3110456) q[0];
sx q[0];
rz(2.8175957) q[0];
rz(1.8404768) q[1];
sx q[1];
rz(-2.3062861) q[1];
sx q[1];
rz(1.3791929) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24810476) q[0];
sx q[0];
rz(-1.4890492) q[0];
sx q[0];
rz(-1.8155314) q[0];
rz(0.46160134) q[2];
sx q[2];
rz(-1.9399376) q[2];
sx q[2];
rz(0.75418562) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0566237) q[1];
sx q[1];
rz(-1.8693722) q[1];
sx q[1];
rz(-1.741239) q[1];
rz(-pi) q[2];
rz(2.2961388) q[3];
sx q[3];
rz(-0.77973706) q[3];
sx q[3];
rz(-0.23507915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1088915) q[2];
sx q[2];
rz(-1.6872493) q[2];
sx q[2];
rz(-0.94318715) q[2];
rz(-0.33106783) q[3];
sx q[3];
rz(-1.7739242) q[3];
sx q[3];
rz(2.846472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11809764) q[0];
sx q[0];
rz(-2.5546615) q[0];
sx q[0];
rz(2.3773637) q[0];
rz(0.14097342) q[1];
sx q[1];
rz(-2.7455175) q[1];
sx q[1];
rz(1.8483298) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8338776) q[0];
sx q[0];
rz(-0.47979646) q[0];
sx q[0];
rz(1.1307554) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6663315) q[2];
sx q[2];
rz(-0.97180688) q[2];
sx q[2];
rz(0.73087382) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5544719) q[1];
sx q[1];
rz(-1.9738102) q[1];
sx q[1];
rz(2.9906669) q[1];
x q[2];
rz(2.7420298) q[3];
sx q[3];
rz(-1.7559397) q[3];
sx q[3];
rz(-0.012133908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3747037) q[2];
sx q[2];
rz(-1.039144) q[2];
sx q[2];
rz(0.58132201) q[2];
rz(0.86822048) q[3];
sx q[3];
rz(-2.7243082) q[3];
sx q[3];
rz(-1.3114312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.593489) q[0];
sx q[0];
rz(-0.5797525) q[0];
sx q[0];
rz(1.2506437) q[0];
rz(-2.1447694) q[1];
sx q[1];
rz(-0.88368982) q[1];
sx q[1];
rz(-1.8539799) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33289136) q[0];
sx q[0];
rz(-1.5455855) q[0];
sx q[0];
rz(-1.0813792) q[0];
rz(-0.31918819) q[2];
sx q[2];
rz(-2.6337998) q[2];
sx q[2];
rz(-2.9920981) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.3567644) q[1];
sx q[1];
rz(-0.67912662) q[1];
sx q[1];
rz(-1.9041512) q[1];
rz(-pi) q[2];
rz(0.82448126) q[3];
sx q[3];
rz(-1.609625) q[3];
sx q[3];
rz(2.3076434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.320497) q[2];
sx q[2];
rz(-2.7536776) q[2];
sx q[2];
rz(-1.256475) q[2];
rz(0.16658941) q[3];
sx q[3];
rz(-1.5649256) q[3];
sx q[3];
rz(-2.0690209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88084108) q[0];
sx q[0];
rz(-0.80417019) q[0];
sx q[0];
rz(2.8163731) q[0];
rz(1.1351769) q[1];
sx q[1];
rz(-0.67712855) q[1];
sx q[1];
rz(2.7744055) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31736483) q[0];
sx q[0];
rz(-1.6019078) q[0];
sx q[0];
rz(0.048006417) q[0];
x q[1];
rz(-0.37819241) q[2];
sx q[2];
rz(-2.4792255) q[2];
sx q[2];
rz(-2.3364002) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.90875188) q[1];
sx q[1];
rz(-2.4526261) q[1];
sx q[1];
rz(-0.31713756) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1807947) q[3];
sx q[3];
rz(-0.73878091) q[3];
sx q[3];
rz(1.184954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8964768) q[2];
sx q[2];
rz(-0.8090691) q[2];
sx q[2];
rz(-1.0661351) q[2];
rz(-3.0623479) q[3];
sx q[3];
rz(-2.5388122) q[3];
sx q[3];
rz(1.665303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29522482) q[0];
sx q[0];
rz(-1.5867148) q[0];
sx q[0];
rz(-1.5750194) q[0];
rz(2.8425343) q[1];
sx q[1];
rz(-0.60332861) q[1];
sx q[1];
rz(0.52287846) q[1];
rz(1.6499741) q[2];
sx q[2];
rz(-2.3186602) q[2];
sx q[2];
rz(3.0849948) q[2];
rz(-0.84898938) q[3];
sx q[3];
rz(-1.6860387) q[3];
sx q[3];
rz(0.54308346) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
