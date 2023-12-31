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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68361002) q[0];
sx q[0];
rz(-1.542122) q[0];
sx q[0];
rz(0.81575127) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6526821) q[2];
sx q[2];
rz(-2.4856644) q[2];
sx q[2];
rz(-1.3002849) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0178535) q[1];
sx q[1];
rz(-2.1542319) q[1];
sx q[1];
rz(2.6052193) q[1];
rz(-pi) q[2];
rz(0.4001873) q[3];
sx q[3];
rz(-1.549198) q[3];
sx q[3];
rz(-2.9633629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2259851) q[2];
sx q[2];
rz(-0.51708022) q[2];
sx q[2];
rz(-1.263164) q[2];
rz(1.8566711) q[3];
sx q[3];
rz(-1.6804755) q[3];
sx q[3];
rz(-2.8485956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.15853515) q[0];
sx q[0];
rz(-1.2936658) q[0];
sx q[0];
rz(-2.7040226) q[0];
rz(0.63105398) q[1];
sx q[1];
rz(-0.32866207) q[1];
sx q[1];
rz(0.22110573) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1572004) q[0];
sx q[0];
rz(-1.7565787) q[0];
sx q[0];
rz(2.0322582) q[0];
x q[1];
rz(0.11110335) q[2];
sx q[2];
rz(-1.9176033) q[2];
sx q[2];
rz(-0.32804104) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4202538) q[1];
sx q[1];
rz(-1.6920648) q[1];
sx q[1];
rz(0.46055693) q[1];
x q[2];
rz(-2.6392194) q[3];
sx q[3];
rz(-1.6564192) q[3];
sx q[3];
rz(-1.8300213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.21800403) q[2];
sx q[2];
rz(-1.709334) q[2];
sx q[2];
rz(-2.5533) q[2];
rz(0.44899392) q[3];
sx q[3];
rz(-0.42789999) q[3];
sx q[3];
rz(1.0480405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56851971) q[0];
sx q[0];
rz(-3.0439215) q[0];
sx q[0];
rz(2.1411238) q[0];
rz(2.4160642) q[1];
sx q[1];
rz(-2.0670481) q[1];
sx q[1];
rz(-0.75769889) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.476186) q[0];
sx q[0];
rz(-1.889567) q[0];
sx q[0];
rz(2.589059) q[0];
rz(-pi) q[1];
rz(1.7311086) q[2];
sx q[2];
rz(-1.3659039) q[2];
sx q[2];
rz(0.93970539) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.12049) q[1];
sx q[1];
rz(-1.1081401) q[1];
sx q[1];
rz(-0.91336577) q[1];
x q[2];
rz(2.9919639) q[3];
sx q[3];
rz(-2.4976839) q[3];
sx q[3];
rz(0.097188918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2086601) q[2];
sx q[2];
rz(-2.9186086) q[2];
sx q[2];
rz(2.3051252) q[2];
rz(1.4423192) q[3];
sx q[3];
rz(-1.6675555) q[3];
sx q[3];
rz(-2.9377655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1214685) q[0];
sx q[0];
rz(-0.080105372) q[0];
sx q[0];
rz(0.51112038) q[0];
rz(-0.077443667) q[1];
sx q[1];
rz(-0.42235342) q[1];
sx q[1];
rz(1.8130594) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3805566) q[0];
sx q[0];
rz(-1.3576926) q[0];
sx q[0];
rz(0.15844945) q[0];
rz(-1.0765692) q[2];
sx q[2];
rz(-0.87072125) q[2];
sx q[2];
rz(-1.1717403) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7718539) q[1];
sx q[1];
rz(-2.8747897) q[1];
sx q[1];
rz(1.5688483) q[1];
x q[2];
rz(-2.02416) q[3];
sx q[3];
rz(-2.8957267) q[3];
sx q[3];
rz(-0.21594957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.093734309) q[2];
sx q[2];
rz(-1.2183943) q[2];
sx q[2];
rz(2.034534) q[2];
rz(-0.47248653) q[3];
sx q[3];
rz(-1.6146336) q[3];
sx q[3];
rz(-2.2842177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1012786) q[0];
sx q[0];
rz(-2.2450228) q[0];
sx q[0];
rz(-2.8033946) q[0];
rz(-1.2942554) q[1];
sx q[1];
rz(-1.5816403) q[1];
sx q[1];
rz(-2.6370874) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.757526) q[0];
sx q[0];
rz(-0.58870089) q[0];
sx q[0];
rz(1.9996044) q[0];
rz(-1.1930824) q[2];
sx q[2];
rz(-1.1337122) q[2];
sx q[2];
rz(-1.9405685) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0658768) q[1];
sx q[1];
rz(-1.8738235) q[1];
sx q[1];
rz(-0.24995835) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2284331) q[3];
sx q[3];
rz(-1.0201766) q[3];
sx q[3];
rz(-0.0084358128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0304886) q[2];
sx q[2];
rz(-2.3031074) q[2];
sx q[2];
rz(-1.6476691) q[2];
rz(1.6882287) q[3];
sx q[3];
rz(-0.93795347) q[3];
sx q[3];
rz(-1.3822314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21022739) q[0];
sx q[0];
rz(-0.87451044) q[0];
sx q[0];
rz(1.0908303) q[0];
rz(2.6018654) q[1];
sx q[1];
rz(-1.153774) q[1];
sx q[1];
rz(-0.18879034) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.323558) q[0];
sx q[0];
rz(-1.3209045) q[0];
sx q[0];
rz(-0.37139335) q[0];
rz(1.9300869) q[2];
sx q[2];
rz(-0.28887666) q[2];
sx q[2];
rz(1.8747683) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.69953883) q[1];
sx q[1];
rz(-0.25146723) q[1];
sx q[1];
rz(2.8738408) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3879697) q[3];
sx q[3];
rz(-1.8086686) q[3];
sx q[3];
rz(0.49926234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.171689) q[2];
sx q[2];
rz(-1.0287372) q[2];
sx q[2];
rz(-1.8264654) q[2];
rz(1.6361489) q[3];
sx q[3];
rz(-0.35741487) q[3];
sx q[3];
rz(1.858254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
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
rz(-1.0514907) q[0];
sx q[0];
rz(-0.83054709) q[0];
sx q[0];
rz(-2.8175957) q[0];
rz(-1.3011159) q[1];
sx q[1];
rz(-2.3062861) q[1];
sx q[1];
rz(1.3791929) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24810476) q[0];
sx q[0];
rz(-1.6525434) q[0];
sx q[0];
rz(-1.8155314) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1629282) q[2];
sx q[2];
rz(-1.9991572) q[2];
sx q[2];
rz(-0.63901627) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.56475793) q[1];
sx q[1];
rz(-1.733629) q[1];
sx q[1];
rz(2.8388883) q[1];
rz(-pi) q[2];
rz(-0.84545387) q[3];
sx q[3];
rz(-2.3618556) q[3];
sx q[3];
rz(-2.9065135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.032701187) q[2];
sx q[2];
rz(-1.6872493) q[2];
sx q[2];
rz(-2.1984055) q[2];
rz(0.33106783) q[3];
sx q[3];
rz(-1.7739242) q[3];
sx q[3];
rz(-2.846472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.023495) q[0];
sx q[0];
rz(-2.5546615) q[0];
sx q[0];
rz(-2.3773637) q[0];
rz(-3.0006192) q[1];
sx q[1];
rz(-0.39607513) q[1];
sx q[1];
rz(-1.8483298) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4828669) q[0];
sx q[0];
rz(-1.7687161) q[0];
sx q[0];
rz(-1.1307964) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6663315) q[2];
sx q[2];
rz(-2.1697858) q[2];
sx q[2];
rz(0.73087382) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.075899374) q[1];
sx q[1];
rz(-1.7095487) q[1];
sx q[1];
rz(-1.9779512) q[1];
rz(-pi) q[2];
rz(-0.44869081) q[3];
sx q[3];
rz(-0.43827) q[3];
sx q[3];
rz(1.1718307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.76688898) q[2];
sx q[2];
rz(-1.039144) q[2];
sx q[2];
rz(0.58132201) q[2];
rz(-0.86822048) q[3];
sx q[3];
rz(-0.41728443) q[3];
sx q[3];
rz(-1.3114312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54810369) q[0];
sx q[0];
rz(-0.5797525) q[0];
sx q[0];
rz(-1.8909489) q[0];
rz(2.1447694) q[1];
sx q[1];
rz(-2.2579028) q[1];
sx q[1];
rz(1.2876127) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8563961) q[0];
sx q[0];
rz(-0.49001339) q[0];
sx q[0];
rz(1.6243837) q[0];
rz(-2.655517) q[2];
sx q[2];
rz(-1.7239778) q[2];
sx q[2];
rz(1.7024405) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.3567644) q[1];
sx q[1];
rz(-0.67912662) q[1];
sx q[1];
rz(1.9041512) q[1];
rz(-pi) q[2];
rz(0.052863315) q[3];
sx q[3];
rz(-0.82517805) q[3];
sx q[3];
rz(-0.70096522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8210956) q[2];
sx q[2];
rz(-2.7536776) q[2];
sx q[2];
rz(-1.256475) q[2];
rz(-0.16658941) q[3];
sx q[3];
rz(-1.5649256) q[3];
sx q[3];
rz(2.0690209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88084108) q[0];
sx q[0];
rz(-0.80417019) q[0];
sx q[0];
rz(2.8163731) q[0];
rz(2.0064158) q[1];
sx q[1];
rz(-0.67712855) q[1];
sx q[1];
rz(0.36718711) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8242278) q[0];
sx q[0];
rz(-1.6019078) q[0];
sx q[0];
rz(3.0935862) q[0];
x q[1];
rz(-2.514421) q[2];
sx q[2];
rz(-1.3417202) q[2];
sx q[2];
rz(2.0723745) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.91017427) q[1];
sx q[1];
rz(-1.3712198) q[1];
sx q[1];
rz(2.477596) q[1];
x q[2];
rz(-0.48093421) q[3];
sx q[3];
rz(-2.1554865) q[3];
sx q[3];
rz(1.9422873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.24511589) q[2];
sx q[2];
rz(-2.3325236) q[2];
sx q[2];
rz(1.0661351) q[2];
rz(0.079244763) q[3];
sx q[3];
rz(-0.60278046) q[3];
sx q[3];
rz(-1.665303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
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
rz(2.3921641) q[2];
sx q[2];
rz(-1.6288169) q[2];
sx q[2];
rz(-1.5734869) q[2];
rz(-1.3973665) q[3];
sx q[3];
rz(-2.4122824) q[3];
sx q[3];
rz(-1.1576049) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
