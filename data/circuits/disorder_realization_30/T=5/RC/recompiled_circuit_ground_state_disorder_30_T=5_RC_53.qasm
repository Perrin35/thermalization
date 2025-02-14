OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8024017) q[0];
sx q[0];
rz(-1.985476) q[0];
sx q[0];
rz(2.6240786) q[0];
rz(1.327688) q[1];
sx q[1];
rz(-2.2496536) q[1];
sx q[1];
rz(0.54203027) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5919246) q[0];
sx q[0];
rz(-0.79093638) q[0];
sx q[0];
rz(2.8998081) q[0];
x q[1];
rz(-1.8250685) q[2];
sx q[2];
rz(-0.82077998) q[2];
sx q[2];
rz(-2.9645142) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0412647) q[1];
sx q[1];
rz(-2.3243743) q[1];
sx q[1];
rz(0.27268091) q[1];
rz(-pi) q[2];
rz(1.3969421) q[3];
sx q[3];
rz(-1.3411634) q[3];
sx q[3];
rz(1.053912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.0077670495) q[2];
sx q[2];
rz(-1.0453036) q[2];
sx q[2];
rz(2.2110151) q[2];
rz(-3.0104356) q[3];
sx q[3];
rz(-2.7635062) q[3];
sx q[3];
rz(1.4906496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5432878) q[0];
sx q[0];
rz(-2.2853993) q[0];
sx q[0];
rz(1.8970066) q[0];
rz(-0.92567956) q[1];
sx q[1];
rz(-1.8857748) q[1];
sx q[1];
rz(1.6474887) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7667395) q[0];
sx q[0];
rz(-1.8962066) q[0];
sx q[0];
rz(-0.84393566) q[0];
rz(-pi) q[1];
rz(-0.41183128) q[2];
sx q[2];
rz(-0.71626012) q[2];
sx q[2];
rz(2.0873983) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8204262) q[1];
sx q[1];
rz(-2.0211206) q[1];
sx q[1];
rz(1.9558332) q[1];
rz(-3.1140987) q[3];
sx q[3];
rz(-2.8373233) q[3];
sx q[3];
rz(0.19434838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1458448) q[2];
sx q[2];
rz(-0.44439849) q[2];
sx q[2];
rz(2.2399529) q[2];
rz(1.7835279) q[3];
sx q[3];
rz(-1.8243022) q[3];
sx q[3];
rz(0.54734126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3667592) q[0];
sx q[0];
rz(-1.9623373) q[0];
sx q[0];
rz(-1.1460079) q[0];
rz(-0.92578069) q[1];
sx q[1];
rz(-2.130276) q[1];
sx q[1];
rz(2.944223) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28846962) q[0];
sx q[0];
rz(-0.98308167) q[0];
sx q[0];
rz(-0.16492407) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.80855753) q[2];
sx q[2];
rz(-2.7312091) q[2];
sx q[2];
rz(1.1754328) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1038405) q[1];
sx q[1];
rz(-1.7985797) q[1];
sx q[1];
rz(0.040018602) q[1];
rz(-pi) q[2];
rz(0.32549627) q[3];
sx q[3];
rz(-2.5739658) q[3];
sx q[3];
rz(-0.63206712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4027412) q[2];
sx q[2];
rz(-0.98029843) q[2];
sx q[2];
rz(2.951238) q[2];
rz(0.061577408) q[3];
sx q[3];
rz(-2.172251) q[3];
sx q[3];
rz(2.0891345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15678081) q[0];
sx q[0];
rz(-1.7233912) q[0];
sx q[0];
rz(-2.5132827) q[0];
rz(-1.5208987) q[1];
sx q[1];
rz(-1.9338806) q[1];
sx q[1];
rz(0.77888387) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4951866) q[0];
sx q[0];
rz(-0.71346006) q[0];
sx q[0];
rz(-0.97795217) q[0];
rz(2.4905048) q[2];
sx q[2];
rz(-0.44241239) q[2];
sx q[2];
rz(-0.96786066) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4877704) q[1];
sx q[1];
rz(-1.5573335) q[1];
sx q[1];
rz(-3.0938992) q[1];
x q[2];
rz(-0.86872673) q[3];
sx q[3];
rz(-1.0208289) q[3];
sx q[3];
rz(1.7977253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4817619) q[2];
sx q[2];
rz(-1.0783106) q[2];
sx q[2];
rz(4.3241186e-05) q[2];
rz(-1.0846042) q[3];
sx q[3];
rz(-1.1559887) q[3];
sx q[3];
rz(2.675975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(0.67149177) q[0];
sx q[0];
rz(-1.4441613) q[0];
sx q[0];
rz(1.2985562) q[0];
rz(-2.5647054) q[1];
sx q[1];
rz(-1.8678317) q[1];
sx q[1];
rz(0.44688046) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53696886) q[0];
sx q[0];
rz(-2.4659221) q[0];
sx q[0];
rz(2.8506822) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2715029) q[2];
sx q[2];
rz(-1.6422049) q[2];
sx q[2];
rz(-1.7091319) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8557103) q[1];
sx q[1];
rz(-1.1709164) q[1];
sx q[1];
rz(0.22004308) q[1];
rz(-pi) q[2];
rz(-0.3077363) q[3];
sx q[3];
rz(-2.947181) q[3];
sx q[3];
rz(0.82378635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.4345066) q[2];
sx q[2];
rz(-2.862317) q[2];
sx q[2];
rz(2.9217829) q[2];
rz(1.48742) q[3];
sx q[3];
rz(-1.4498962) q[3];
sx q[3];
rz(-2.4421104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(0.18911067) q[0];
sx q[0];
rz(-2.7087152) q[0];
sx q[0];
rz(-0.89749807) q[0];
rz(-1.3709566) q[1];
sx q[1];
rz(-2.1015344) q[1];
sx q[1];
rz(2.2460361) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8093256) q[0];
sx q[0];
rz(-0.73609645) q[0];
sx q[0];
rz(0.21047896) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.95889556) q[2];
sx q[2];
rz(-2.09225) q[2];
sx q[2];
rz(0.34470169) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.31187801) q[1];
sx q[1];
rz(-0.51110454) q[1];
sx q[1];
rz(-2.4097777) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0895133) q[3];
sx q[3];
rz(-1.9949504) q[3];
sx q[3];
rz(-1.7103575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.36876496) q[2];
sx q[2];
rz(-2.0389098) q[2];
sx q[2];
rz(-0.15711288) q[2];
rz(-2.469192) q[3];
sx q[3];
rz(-1.7417358) q[3];
sx q[3];
rz(-3.0290643) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2812578) q[0];
sx q[0];
rz(-0.66547886) q[0];
sx q[0];
rz(1.459664) q[0];
rz(2.7809987) q[1];
sx q[1];
rz(-1.4692042) q[1];
sx q[1];
rz(-1.8416539) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4091638) q[0];
sx q[0];
rz(-2.7061908) q[0];
sx q[0];
rz(2.8485427) q[0];
rz(-pi) q[1];
rz(0.62076135) q[2];
sx q[2];
rz(-1.1859535) q[2];
sx q[2];
rz(0.2523202) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5977138) q[1];
sx q[1];
rz(-1.8671486) q[1];
sx q[1];
rz(1.6719901) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1645369) q[3];
sx q[3];
rz(-0.99525577) q[3];
sx q[3];
rz(1.5062603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2099057) q[2];
sx q[2];
rz(-1.0301215) q[2];
sx q[2];
rz(2.9295861) q[2];
rz(1.0509521) q[3];
sx q[3];
rz(-1.7651599) q[3];
sx q[3];
rz(-0.088137805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5365005) q[0];
sx q[0];
rz(-2.3516646) q[0];
sx q[0];
rz(-2.9686046) q[0];
rz(-1.1146924) q[1];
sx q[1];
rz(-1.0429017) q[1];
sx q[1];
rz(-3.0317422) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2762624) q[0];
sx q[0];
rz(-2.1234351) q[0];
sx q[0];
rz(-1.7981973) q[0];
rz(-pi) q[1];
rz(1.2908247) q[2];
sx q[2];
rz(-2.3738573) q[2];
sx q[2];
rz(2.6565832) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.122095) q[1];
sx q[1];
rz(-1.8260584) q[1];
sx q[1];
rz(2.9649516) q[1];
rz(-pi) q[2];
x q[2];
rz(0.55446378) q[3];
sx q[3];
rz(-1.4687456) q[3];
sx q[3];
rz(-2.3481675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7383808) q[2];
sx q[2];
rz(-1.9387559) q[2];
sx q[2];
rz(0.66486764) q[2];
rz(2.6730149) q[3];
sx q[3];
rz(-1.6178308) q[3];
sx q[3];
rz(2.328228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.453603) q[0];
sx q[0];
rz(-2.5316694) q[0];
sx q[0];
rz(2.4000121) q[0];
rz(-1.954151) q[1];
sx q[1];
rz(-1.5267742) q[1];
sx q[1];
rz(-0.56914079) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5461972) q[0];
sx q[0];
rz(-0.44537195) q[0];
sx q[0];
rz(-2.4085042) q[0];
x q[1];
rz(2.9958565) q[2];
sx q[2];
rz(-0.45782858) q[2];
sx q[2];
rz(2.482058) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5829825) q[1];
sx q[1];
rz(-0.64839586) q[1];
sx q[1];
rz(2.3984539) q[1];
x q[2];
rz(0.83888104) q[3];
sx q[3];
rz(-2.122022) q[3];
sx q[3];
rz(-1.4428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6471214) q[2];
sx q[2];
rz(-0.98860604) q[2];
sx q[2];
rz(-0.76623255) q[2];
rz(-1.1726441) q[3];
sx q[3];
rz(-1.5394883) q[3];
sx q[3];
rz(-2.9197781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93429339) q[0];
sx q[0];
rz(-0.3485637) q[0];
sx q[0];
rz(2.9497414) q[0];
rz(0.81709298) q[1];
sx q[1];
rz(-0.25661883) q[1];
sx q[1];
rz(0.70768913) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95078844) q[0];
sx q[0];
rz(-0.30256264) q[0];
sx q[0];
rz(-2.5214419) q[0];
rz(-pi) q[1];
rz(0.3336256) q[2];
sx q[2];
rz(-1.0051491) q[2];
sx q[2];
rz(-2.8666292) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0470815) q[1];
sx q[1];
rz(-2.1489177) q[1];
sx q[1];
rz(-0.44433997) q[1];
x q[2];
rz(0.97174208) q[3];
sx q[3];
rz(-2.0721367) q[3];
sx q[3];
rz(2.3399692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2231458) q[2];
sx q[2];
rz(-0.94474363) q[2];
sx q[2];
rz(-0.72500149) q[2];
rz(0.21507344) q[3];
sx q[3];
rz(-0.30078617) q[3];
sx q[3];
rz(-0.71896368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9216777) q[0];
sx q[0];
rz(-2.324993) q[0];
sx q[0];
rz(-2.8391229) q[0];
rz(-2.0401781) q[1];
sx q[1];
rz(-1.9093724) q[1];
sx q[1];
rz(-2.2473635) q[1];
rz(-0.11753043) q[2];
sx q[2];
rz(-1.0385658) q[2];
sx q[2];
rz(-0.081296878) q[2];
rz(-1.5433031) q[3];
sx q[3];
rz(-2.0381654) q[3];
sx q[3];
rz(0.5034133) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
