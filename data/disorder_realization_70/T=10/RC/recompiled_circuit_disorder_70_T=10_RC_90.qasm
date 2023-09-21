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
rz(-0.17012574) q[0];
sx q[0];
rz(2.3556019) q[0];
rz(0.6056447) q[1];
sx q[1];
rz(3.7925386) q[1];
sx q[1];
rz(8.7906919) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91416042) q[0];
sx q[0];
rz(-2.3254546) q[0];
sx q[0];
rz(-0.039365191) q[0];
rz(-pi) q[1];
rz(-0.91648957) q[2];
sx q[2];
rz(-1.620703) q[2];
sx q[2];
rz(-2.8061342) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7635599) q[1];
sx q[1];
rz(-1.130192) q[1];
sx q[1];
rz(0.9159169) q[1];
rz(1.5942469) q[3];
sx q[3];
rz(-1.970885) q[3];
sx q[3];
rz(-1.401702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2259851) q[2];
sx q[2];
rz(-0.51708022) q[2];
sx q[2];
rz(-1.263164) q[2];
rz(1.8566711) q[3];
sx q[3];
rz(-1.6804755) q[3];
sx q[3];
rz(0.29299709) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9830575) q[0];
sx q[0];
rz(-1.8479269) q[0];
sx q[0];
rz(-0.43757004) q[0];
rz(-0.63105398) q[1];
sx q[1];
rz(-2.8129306) q[1];
sx q[1];
rz(-2.9204869) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76925812) q[0];
sx q[0];
rz(-2.6466469) q[0];
sx q[0];
rz(1.1713722) q[0];
x q[1];
rz(-0.11110335) q[2];
sx q[2];
rz(-1.2239893) q[2];
sx q[2];
rz(2.8135516) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0529777) q[1];
sx q[1];
rz(-0.47514519) q[1];
sx q[1];
rz(2.8739724) q[1];
rz(0.17640555) q[3];
sx q[3];
rz(-0.50900148) q[3];
sx q[3];
rz(-2.727946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9235886) q[2];
sx q[2];
rz(-1.709334) q[2];
sx q[2];
rz(-2.5533) q[2];
rz(-0.44899392) q[3];
sx q[3];
rz(-0.42789999) q[3];
sx q[3];
rz(2.0935521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56851971) q[0];
sx q[0];
rz(-3.0439215) q[0];
sx q[0];
rz(-1.0004689) q[0];
rz(2.4160642) q[1];
sx q[1];
rz(-2.0670481) q[1];
sx q[1];
rz(2.3838938) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.476186) q[0];
sx q[0];
rz(-1.2520257) q[0];
sx q[0];
rz(0.55253367) q[0];
rz(-pi) q[1];
x q[1];
rz(0.65501113) q[2];
sx q[2];
rz(-0.25946028) q[2];
sx q[2];
rz(-0.26817817) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1670907) q[1];
sx q[1];
rz(-2.3579512) q[1];
sx q[1];
rz(-2.2553315) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9919639) q[3];
sx q[3];
rz(-0.64390874) q[3];
sx q[3];
rz(-3.0444037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9329325) q[2];
sx q[2];
rz(-0.22298403) q[2];
sx q[2];
rz(2.3051252) q[2];
rz(-1.6992735) q[3];
sx q[3];
rz(-1.6675555) q[3];
sx q[3];
rz(0.20382717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-3.1214685) q[0];
sx q[0];
rz(-3.0614873) q[0];
sx q[0];
rz(0.51112038) q[0];
rz(3.064149) q[1];
sx q[1];
rz(-2.7192392) q[1];
sx q[1];
rz(-1.8130594) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9851345) q[0];
sx q[0];
rz(-1.4159604) q[0];
sx q[0];
rz(1.7865208) q[0];
x q[1];
rz(-2.0650234) q[2];
sx q[2];
rz(-0.87072125) q[2];
sx q[2];
rz(-1.9698524) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7718539) q[1];
sx q[1];
rz(-2.8747897) q[1];
sx q[1];
rz(-1.5688483) q[1];
x q[2];
rz(1.7926746) q[3];
sx q[3];
rz(-1.6776049) q[3];
sx q[3];
rz(-0.91339236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.093734309) q[2];
sx q[2];
rz(-1.9231984) q[2];
sx q[2];
rz(-2.034534) q[2];
rz(-2.6691061) q[3];
sx q[3];
rz(-1.5269591) q[3];
sx q[3];
rz(-2.2842177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.040314019) q[0];
sx q[0];
rz(-2.2450228) q[0];
sx q[0];
rz(0.3381981) q[0];
rz(-1.2942554) q[1];
sx q[1];
rz(-1.5599524) q[1];
sx q[1];
rz(-0.50450528) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38406661) q[0];
sx q[0];
rz(-2.5528918) q[0];
sx q[0];
rz(-1.1419883) q[0];
rz(-pi) q[1];
x q[1];
rz(0.46576969) q[2];
sx q[2];
rz(-1.911474) q[2];
sx q[2];
rz(-0.53616947) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.63197631) q[1];
sx q[1];
rz(-0.39034931) q[1];
sx q[1];
rz(0.901464) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4818146) q[3];
sx q[3];
rz(-1.0228844) q[3];
sx q[3];
rz(-1.1783311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0304886) q[2];
sx q[2];
rz(-2.3031074) q[2];
sx q[2];
rz(-1.4939235) q[2];
rz(1.4533639) q[3];
sx q[3];
rz(-2.2036392) q[3];
sx q[3];
rz(1.7593613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9313653) q[0];
sx q[0];
rz(-0.87451044) q[0];
sx q[0];
rz(1.0908303) q[0];
rz(2.6018654) q[1];
sx q[1];
rz(-1.9878186) q[1];
sx q[1];
rz(0.18879034) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.328619) q[0];
sx q[0];
rz(-0.44437528) q[0];
sx q[0];
rz(-2.5286753) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9300869) q[2];
sx q[2];
rz(-0.28887666) q[2];
sx q[2];
rz(1.8747683) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.69953883) q[1];
sx q[1];
rz(-2.8901254) q[1];
sx q[1];
rz(-2.8738408) q[1];
x q[2];
rz(2.4981899) q[3];
sx q[3];
rz(-2.8426369) q[3];
sx q[3];
rz(1.1645731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9699036) q[2];
sx q[2];
rz(-2.1128555) q[2];
sx q[2];
rz(1.8264654) q[2];
rz(1.6361489) q[3];
sx q[3];
rz(-2.7841778) q[3];
sx q[3];
rz(-1.858254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0514907) q[0];
sx q[0];
rz(-2.3110456) q[0];
sx q[0];
rz(2.8175957) q[0];
rz(-1.3011159) q[1];
sx q[1];
rz(-2.3062861) q[1];
sx q[1];
rz(-1.7623998) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1349072) q[0];
sx q[0];
rz(-0.25776699) q[0];
sx q[0];
rz(-1.8968614) q[0];
rz(1.1629282) q[2];
sx q[2];
rz(-1.1424354) q[2];
sx q[2];
rz(-2.5025764) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6143601) q[1];
sx q[1];
rz(-2.7990606) q[1];
sx q[1];
rz(2.637898) q[1];
x q[2];
rz(2.561065) q[3];
sx q[3];
rz(-2.1248098) q[3];
sx q[3];
rz(1.1298657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.032701187) q[2];
sx q[2];
rz(-1.6872493) q[2];
sx q[2];
rz(0.94318715) q[2];
rz(0.33106783) q[3];
sx q[3];
rz(-1.7739242) q[3];
sx q[3];
rz(-2.846472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.023495) q[0];
sx q[0];
rz(-2.5546615) q[0];
sx q[0];
rz(2.3773637) q[0];
rz(-0.14097342) q[1];
sx q[1];
rz(-2.7455175) q[1];
sx q[1];
rz(-1.8483298) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4828669) q[0];
sx q[0];
rz(-1.7687161) q[0];
sx q[0];
rz(-2.0107962) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0027577) q[2];
sx q[2];
rz(-2.5359557) q[2];
sx q[2];
rz(-2.5790737) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1845491) q[1];
sx q[1];
rz(-0.42889412) q[1];
sx q[1];
rz(-1.2317608) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.771365) q[3];
sx q[3];
rz(-1.9631533) q[3];
sx q[3];
rz(1.6605103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.76688898) q[2];
sx q[2];
rz(-2.1024487) q[2];
sx q[2];
rz(-2.5602706) q[2];
rz(2.2733722) q[3];
sx q[3];
rz(-2.7243082) q[3];
sx q[3];
rz(-1.8301615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54810369) q[0];
sx q[0];
rz(-2.5618401) q[0];
sx q[0];
rz(1.8909489) q[0];
rz(2.1447694) q[1];
sx q[1];
rz(-0.88368982) q[1];
sx q[1];
rz(1.8539799) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8087013) q[0];
sx q[0];
rz(-1.5455855) q[0];
sx q[0];
rz(-1.0813792) q[0];
rz(-pi) q[1];
x q[1];
rz(2.655517) q[2];
sx q[2];
rz(-1.7239778) q[2];
sx q[2];
rz(-1.7024405) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.3567644) q[1];
sx q[1];
rz(-2.462466) q[1];
sx q[1];
rz(-1.2374415) q[1];
rz(3.0887293) q[3];
sx q[3];
rz(-2.3164146) q[3];
sx q[3];
rz(-0.70096522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8210956) q[2];
sx q[2];
rz(-0.38791502) q[2];
sx q[2];
rz(-1.256475) q[2];
rz(2.9750032) q[3];
sx q[3];
rz(-1.5649256) q[3];
sx q[3];
rz(-1.0725718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88084108) q[0];
sx q[0];
rz(-0.80417019) q[0];
sx q[0];
rz(-2.8163731) q[0];
rz(-1.1351769) q[1];
sx q[1];
rz(-2.4644641) q[1];
sx q[1];
rz(-0.36718711) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8866667) q[0];
sx q[0];
rz(-1.5228132) q[0];
sx q[0];
rz(1.539649) q[0];
rz(-2.7634002) q[2];
sx q[2];
rz(-2.4792255) q[2];
sx q[2];
rz(-0.80519245) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.90875188) q[1];
sx q[1];
rz(-0.68896657) q[1];
sx q[1];
rz(-2.8244551) q[1];
rz(-0.96079798) q[3];
sx q[3];
rz(-0.73878091) q[3];
sx q[3];
rz(-1.184954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.24511589) q[2];
sx q[2];
rz(-0.8090691) q[2];
sx q[2];
rz(-2.0754576) q[2];
rz(3.0623479) q[3];
sx q[3];
rz(-2.5388122) q[3];
sx q[3];
rz(-1.665303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8463678) q[0];
sx q[0];
rz(-1.5867148) q[0];
sx q[0];
rz(-1.5750194) q[0];
rz(-2.8425343) q[1];
sx q[1];
rz(-2.538264) q[1];
sx q[1];
rz(-2.6187142) q[1];
rz(-1.6499741) q[2];
sx q[2];
rz(-0.82293246) q[2];
sx q[2];
rz(-0.056597829) q[2];
rz(-2.2926033) q[3];
sx q[3];
rz(-1.4555539) q[3];
sx q[3];
rz(-2.5985092) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];