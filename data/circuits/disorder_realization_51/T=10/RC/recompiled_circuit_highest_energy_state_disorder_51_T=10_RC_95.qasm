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
rz(0.32938862) q[0];
sx q[0];
rz(3.7731054) q[0];
sx q[0];
rz(9.4256529) q[0];
rz(-0.64088351) q[1];
sx q[1];
rz(-1.0008608) q[1];
sx q[1];
rz(0.34520087) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5986818) q[0];
sx q[0];
rz(-0.094176725) q[0];
sx q[0];
rz(-1.8369294) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7183185) q[2];
sx q[2];
rz(-1.194724) q[2];
sx q[2];
rz(-3.0126115) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.36269128) q[1];
sx q[1];
rz(-2.3786015) q[1];
sx q[1];
rz(0.94339006) q[1];
rz(3.1253417) q[3];
sx q[3];
rz(-1.4964087) q[3];
sx q[3];
rz(-0.95916884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.22902809) q[2];
sx q[2];
rz(-1.3338858) q[2];
sx q[2];
rz(0.20797569) q[2];
rz(0.29911706) q[3];
sx q[3];
rz(-2.5695473) q[3];
sx q[3];
rz(1.1408898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8107373) q[0];
sx q[0];
rz(-2.9006697) q[0];
sx q[0];
rz(-0.99288565) q[0];
rz(1.8244686) q[1];
sx q[1];
rz(-0.31610745) q[1];
sx q[1];
rz(0.77450007) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74906384) q[0];
sx q[0];
rz(-1.7490024) q[0];
sx q[0];
rz(3.130413) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0335017) q[2];
sx q[2];
rz(-2.3024493) q[2];
sx q[2];
rz(-3.0947859) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6933813) q[1];
sx q[1];
rz(-0.44084545) q[1];
sx q[1];
rz(0.54852672) q[1];
x q[2];
rz(-0.22529545) q[3];
sx q[3];
rz(-1.1114128) q[3];
sx q[3];
rz(-2.6090906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.896686) q[2];
sx q[2];
rz(-1.8550355) q[2];
sx q[2];
rz(-1.0150821) q[2];
rz(-0.24909881) q[3];
sx q[3];
rz(-0.86779147) q[3];
sx q[3];
rz(2.7740313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2772813) q[0];
sx q[0];
rz(-0.69121498) q[0];
sx q[0];
rz(0.15750289) q[0];
rz(1.0109673) q[1];
sx q[1];
rz(-2.8087661) q[1];
sx q[1];
rz(0.13883042) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75051266) q[0];
sx q[0];
rz(-1.2475999) q[0];
sx q[0];
rz(-1.1884407) q[0];
x q[1];
rz(-0.71535297) q[2];
sx q[2];
rz(-1.2497447) q[2];
sx q[2];
rz(0.80594826) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5032387) q[1];
sx q[1];
rz(-1.279083) q[1];
sx q[1];
rz(2.9364763) q[1];
rz(-pi) q[2];
rz(-0.28335684) q[3];
sx q[3];
rz(-1.4025926) q[3];
sx q[3];
rz(3.0137872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.46131721) q[2];
sx q[2];
rz(-0.91247827) q[2];
sx q[2];
rz(0.3581363) q[2];
rz(2.4628468) q[3];
sx q[3];
rz(-2.1923784) q[3];
sx q[3];
rz(-2.1024735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.045227483) q[0];
sx q[0];
rz(-2.1684833) q[0];
sx q[0];
rz(-0.3048234) q[0];
rz(2.7335956) q[1];
sx q[1];
rz(-1.7034737) q[1];
sx q[1];
rz(-2.2136484) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5261032) q[0];
sx q[0];
rz(-2.9080221) q[0];
sx q[0];
rz(2.156267) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5408542) q[2];
sx q[2];
rz(-2.1081807) q[2];
sx q[2];
rz(-1.5058277) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.881262) q[1];
sx q[1];
rz(-1.4196383) q[1];
sx q[1];
rz(0.69044729) q[1];
rz(0.66189142) q[3];
sx q[3];
rz(-2.3958979) q[3];
sx q[3];
rz(0.5058561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1912332) q[2];
sx q[2];
rz(-0.57098907) q[2];
sx q[2];
rz(2.9534269) q[2];
rz(-1.865271) q[3];
sx q[3];
rz(-1.3151582) q[3];
sx q[3];
rz(1.4307384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6998049) q[0];
sx q[0];
rz(-0.34407523) q[0];
sx q[0];
rz(-3.1222043) q[0];
rz(2.0806606) q[1];
sx q[1];
rz(-1.7273936) q[1];
sx q[1];
rz(-1.9020938) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0415619) q[0];
sx q[0];
rz(-1.2912371) q[0];
sx q[0];
rz(0.97645219) q[0];
rz(2.2197228) q[2];
sx q[2];
rz(-2.1989294) q[2];
sx q[2];
rz(-2.3482196) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0026999) q[1];
sx q[1];
rz(-2.5475916) q[1];
sx q[1];
rz(-2.6245489) q[1];
x q[2];
rz(-0.40941671) q[3];
sx q[3];
rz(-2.0338661) q[3];
sx q[3];
rz(-0.45209979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1218607) q[2];
sx q[2];
rz(-1.3132361) q[2];
sx q[2];
rz(1.8149553) q[2];
rz(-2.895368) q[3];
sx q[3];
rz(-1.1028057) q[3];
sx q[3];
rz(0.8518014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5354079) q[0];
sx q[0];
rz(-2.1562205) q[0];
sx q[0];
rz(-0.73915172) q[0];
rz(-0.34573653) q[1];
sx q[1];
rz(-1.3290936) q[1];
sx q[1];
rz(-2.2703222) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3987253) q[0];
sx q[0];
rz(-2.350432) q[0];
sx q[0];
rz(-1.6683116) q[0];
rz(-pi) q[1];
rz(-0.99314697) q[2];
sx q[2];
rz(-2.2340074) q[2];
sx q[2];
rz(-0.011485966) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4193486) q[1];
sx q[1];
rz(-0.39545317) q[1];
sx q[1];
rz(1.761318) q[1];
x q[2];
rz(2.5090747) q[3];
sx q[3];
rz(-0.7154724) q[3];
sx q[3];
rz(-2.4533437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.31200108) q[2];
sx q[2];
rz(-2.5025554) q[2];
sx q[2];
rz(2.269022) q[2];
rz(1.5729337) q[3];
sx q[3];
rz(-2.5376153) q[3];
sx q[3];
rz(2.2085371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0912112) q[0];
sx q[0];
rz(-0.25930431) q[0];
sx q[0];
rz(0.76164371) q[0];
rz(-2.7161982) q[1];
sx q[1];
rz(-2.0014706) q[1];
sx q[1];
rz(2.2147307) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2970718) q[0];
sx q[0];
rz(-1.335134) q[0];
sx q[0];
rz(-1.8051487) q[0];
rz(-pi) q[1];
rz(2.4270646) q[2];
sx q[2];
rz(-2.0863669) q[2];
sx q[2];
rz(0.58802468) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.35173479) q[1];
sx q[1];
rz(-2.7119293) q[1];
sx q[1];
rz(-2.1273145) q[1];
rz(-0.048153444) q[3];
sx q[3];
rz(-2.1934436) q[3];
sx q[3];
rz(2.4751055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0285792) q[2];
sx q[2];
rz(-2.4455652) q[2];
sx q[2];
rz(-0.87116233) q[2];
rz(0.29843676) q[3];
sx q[3];
rz(-1.8961743) q[3];
sx q[3];
rz(-1.3309853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7262064) q[0];
sx q[0];
rz(-2.6415249) q[0];
sx q[0];
rz(2.723208) q[0];
rz(0.2977953) q[1];
sx q[1];
rz(-1.9458385) q[1];
sx q[1];
rz(3.0072838) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90420049) q[0];
sx q[0];
rz(-0.89590329) q[0];
sx q[0];
rz(2.2740721) q[0];
rz(1.1654794) q[2];
sx q[2];
rz(-1.3994872) q[2];
sx q[2];
rz(1.5660945) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2063576) q[1];
sx q[1];
rz(-0.22769732) q[1];
sx q[1];
rz(0.63657324) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0928602) q[3];
sx q[3];
rz(-0.52937859) q[3];
sx q[3];
rz(-3.0335308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7316651) q[2];
sx q[2];
rz(-1.2890559) q[2];
sx q[2];
rz(0.64201391) q[2];
rz(1.0632473) q[3];
sx q[3];
rz(-0.65170538) q[3];
sx q[3];
rz(-2.1003907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6006271) q[0];
sx q[0];
rz(-0.0890812) q[0];
sx q[0];
rz(-0.11216057) q[0];
rz(1.8070096) q[1];
sx q[1];
rz(-0.59252512) q[1];
sx q[1];
rz(1.0293915) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5350069) q[0];
sx q[0];
rz(-2.2748053) q[0];
sx q[0];
rz(-3.0993942) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5310982) q[2];
sx q[2];
rz(-1.8744812) q[2];
sx q[2];
rz(3.1162709) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2160714) q[1];
sx q[1];
rz(-0.23434429) q[1];
sx q[1];
rz(-1.3804803) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6736027) q[3];
sx q[3];
rz(-1.6525998) q[3];
sx q[3];
rz(1.8379267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7182497) q[2];
sx q[2];
rz(-0.074904718) q[2];
sx q[2];
rz(-0.038012803) q[2];
rz(-2.5293317) q[3];
sx q[3];
rz(-2.2050048) q[3];
sx q[3];
rz(-0.14122252) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24899471) q[0];
sx q[0];
rz(-1.6706415) q[0];
sx q[0];
rz(-0.41879642) q[0];
rz(1.2576125) q[1];
sx q[1];
rz(-1.5092756) q[1];
sx q[1];
rz(0.14258252) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6194942) q[0];
sx q[0];
rz(-1.5987248) q[0];
sx q[0];
rz(-0.15837196) q[0];
rz(-pi) q[1];
rz(-1.293574) q[2];
sx q[2];
rz(-0.96053329) q[2];
sx q[2];
rz(1.1278314) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.31227641) q[1];
sx q[1];
rz(-1.3061151) q[1];
sx q[1];
rz(-2.1618202) q[1];
rz(-pi) q[2];
rz(0.3440109) q[3];
sx q[3];
rz(-1.1340525) q[3];
sx q[3];
rz(-0.94319447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3451781) q[2];
sx q[2];
rz(-0.08064457) q[2];
sx q[2];
rz(-2.8734015) q[2];
rz(1.7372519) q[3];
sx q[3];
rz(-1.0920478) q[3];
sx q[3];
rz(2.0973189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0161229) q[0];
sx q[0];
rz(-0.37624993) q[0];
sx q[0];
rz(-2.7952623) q[0];
rz(0.12915962) q[1];
sx q[1];
rz(-1.2551413) q[1];
sx q[1];
rz(-1.7358949) q[1];
rz(1.9228946) q[2];
sx q[2];
rz(-2.0205971) q[2];
sx q[2];
rz(-1.481075) q[2];
rz(1.7862475) q[3];
sx q[3];
rz(-2.8402495) q[3];
sx q[3];
rz(0.80339669) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
