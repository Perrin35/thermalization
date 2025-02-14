OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.793279) q[0];
sx q[0];
rz(-1.7338742) q[0];
sx q[0];
rz(-3.0970567) q[0];
rz(-1.1107923) q[1];
sx q[1];
rz(5.0587237) q[1];
sx q[1];
rz(8.6375477) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4314595) q[0];
sx q[0];
rz(-1.077912) q[0];
sx q[0];
rz(-2.5566275) q[0];
rz(1.2455315) q[2];
sx q[2];
rz(-1.6084387) q[2];
sx q[2];
rz(1.4311439) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2582142) q[1];
sx q[1];
rz(-1.6401263) q[1];
sx q[1];
rz(-0.99154559) q[1];
x q[2];
rz(-0.23870339) q[3];
sx q[3];
rz(-1.9145962) q[3];
sx q[3];
rz(1.3581004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.74701509) q[2];
sx q[2];
rz(-2.0614823) q[2];
sx q[2];
rz(-2.3874095) q[2];
rz(-1.7265823) q[3];
sx q[3];
rz(-0.49608803) q[3];
sx q[3];
rz(-2.8126341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0542145) q[0];
sx q[0];
rz(-2.7890451) q[0];
sx q[0];
rz(-1.8074328) q[0];
rz(-1.2913903) q[1];
sx q[1];
rz(-2.1033557) q[1];
sx q[1];
rz(-2.2659567) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9459045) q[0];
sx q[0];
rz(-2.2155846) q[0];
sx q[0];
rz(-0.88627215) q[0];
rz(-pi) q[1];
rz(1.6853784) q[2];
sx q[2];
rz(-2.6422524) q[2];
sx q[2];
rz(-0.90081604) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7499979) q[1];
sx q[1];
rz(-1.1340967) q[1];
sx q[1];
rz(2.8207247) q[1];
rz(2.528627) q[3];
sx q[3];
rz(-0.35109441) q[3];
sx q[3];
rz(2.8610817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.69677532) q[2];
sx q[2];
rz(-0.70187086) q[2];
sx q[2];
rz(-2.2689421) q[2];
rz(-2.3526092) q[3];
sx q[3];
rz(-2.240447) q[3];
sx q[3];
rz(-0.22855973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2149684) q[0];
sx q[0];
rz(-1.3525532) q[0];
sx q[0];
rz(0.0080000814) q[0];
rz(-2.3232715) q[1];
sx q[1];
rz(-2.3921831) q[1];
sx q[1];
rz(-1.2368894) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.711126) q[0];
sx q[0];
rz(-2.240181) q[0];
sx q[0];
rz(0.41394625) q[0];
rz(-pi) q[1];
rz(1.1367927) q[2];
sx q[2];
rz(-1.2522337) q[2];
sx q[2];
rz(0.54036507) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9489386) q[1];
sx q[1];
rz(-1.0315507) q[1];
sx q[1];
rz(1.2112898) q[1];
rz(1.6881006) q[3];
sx q[3];
rz(-1.8059732) q[3];
sx q[3];
rz(-1.7010353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7311953) q[2];
sx q[2];
rz(-3.0072913) q[2];
sx q[2];
rz(-2.2006939) q[2];
rz(1.1040374) q[3];
sx q[3];
rz(-1.3528115) q[3];
sx q[3];
rz(1.997939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10821548) q[0];
sx q[0];
rz(-2.5984851) q[0];
sx q[0];
rz(-2.9040842) q[0];
rz(-2.4820651) q[1];
sx q[1];
rz(-0.57412761) q[1];
sx q[1];
rz(-3.1245756) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68786808) q[0];
sx q[0];
rz(-0.17565082) q[0];
sx q[0];
rz(1.2916628) q[0];
rz(2.723316) q[2];
sx q[2];
rz(-2.661663) q[2];
sx q[2];
rz(-2.0945702) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5524872) q[1];
sx q[1];
rz(-1.0081588) q[1];
sx q[1];
rz(1.6943114) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4260826) q[3];
sx q[3];
rz(-0.9431211) q[3];
sx q[3];
rz(2.024141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.235405) q[2];
sx q[2];
rz(-1.2013288) q[2];
sx q[2];
rz(0.91442937) q[2];
rz(0.36214456) q[3];
sx q[3];
rz(-2.3320964) q[3];
sx q[3];
rz(0.58524281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.888716) q[0];
sx q[0];
rz(-1.7676366) q[0];
sx q[0];
rz(0.1121029) q[0];
rz(2.0762699) q[1];
sx q[1];
rz(-2.0707097) q[1];
sx q[1];
rz(-2.0723453) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1493268) q[0];
sx q[0];
rz(-2.7130824) q[0];
sx q[0];
rz(2.4403768) q[0];
rz(-pi) q[1];
rz(-0.5324131) q[2];
sx q[2];
rz(-1.1523917) q[2];
sx q[2];
rz(-2.6765228) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.62110626) q[1];
sx q[1];
rz(-2.6992528) q[1];
sx q[1];
rz(0.81653313) q[1];
rz(-pi) q[2];
rz(0.97158708) q[3];
sx q[3];
rz(-1.6808482) q[3];
sx q[3];
rz(-0.19796619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4273044) q[2];
sx q[2];
rz(-2.4236743) q[2];
sx q[2];
rz(-0.53421268) q[2];
rz(-0.30321768) q[3];
sx q[3];
rz(-1.5568308) q[3];
sx q[3];
rz(1.6259954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7233906) q[0];
sx q[0];
rz(-1.8573107) q[0];
sx q[0];
rz(0.89163017) q[0];
rz(1.760969) q[1];
sx q[1];
rz(-1.9386407) q[1];
sx q[1];
rz(-2.2183529) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5017446) q[0];
sx q[0];
rz(-3.0303133) q[0];
sx q[0];
rz(-2.6529045) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7545287) q[2];
sx q[2];
rz(-1.825807) q[2];
sx q[2];
rz(-0.86266359) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.504661) q[1];
sx q[1];
rz(-1.8430222) q[1];
sx q[1];
rz(-0.23761655) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4999077) q[3];
sx q[3];
rz(-1.2417792) q[3];
sx q[3];
rz(-1.6365479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0455857) q[2];
sx q[2];
rz(-0.43232375) q[2];
sx q[2];
rz(-1.879479) q[2];
rz(1.9803842) q[3];
sx q[3];
rz(-1.592417) q[3];
sx q[3];
rz(-1.2459292) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4534509) q[0];
sx q[0];
rz(-2.3095486) q[0];
sx q[0];
rz(-1.6081109) q[0];
rz(0.43117943) q[1];
sx q[1];
rz(-2.2240413) q[1];
sx q[1];
rz(-1.7327488) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20074318) q[0];
sx q[0];
rz(-2.4398167) q[0];
sx q[0];
rz(-0.45325847) q[0];
x q[1];
rz(1.3546014) q[2];
sx q[2];
rz(-0.39297418) q[2];
sx q[2];
rz(-2.4525688) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5049487) q[1];
sx q[1];
rz(-1.4145936) q[1];
sx q[1];
rz(-1.4438932) q[1];
x q[2];
rz(0.124229) q[3];
sx q[3];
rz(-1.0719187) q[3];
sx q[3];
rz(1.6586097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.93671736) q[2];
sx q[2];
rz(-1.1254346) q[2];
sx q[2];
rz(-1.8219061) q[2];
rz(-0.034505757) q[3];
sx q[3];
rz(-1.6104108) q[3];
sx q[3];
rz(1.7721133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40783229) q[0];
sx q[0];
rz(-0.5540846) q[0];
sx q[0];
rz(1.2169417) q[0];
rz(-2.2159684) q[1];
sx q[1];
rz(-1.3652912) q[1];
sx q[1];
rz(1.2006753) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72382765) q[0];
sx q[0];
rz(-2.1795131) q[0];
sx q[0];
rz(-0.62792741) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.56270069) q[2];
sx q[2];
rz(-0.68018736) q[2];
sx q[2];
rz(-2.0923751) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6457121) q[1];
sx q[1];
rz(-2.3420942) q[1];
sx q[1];
rz(-1.6281566) q[1];
x q[2];
rz(-1.993409) q[3];
sx q[3];
rz(-1.9070996) q[3];
sx q[3];
rz(-2.5667532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.22244464) q[2];
sx q[2];
rz(-1.5347975) q[2];
sx q[2];
rz(-1.0486802) q[2];
rz(2.9564296) q[3];
sx q[3];
rz(-1.1566297) q[3];
sx q[3];
rz(3.0345501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1222526) q[0];
sx q[0];
rz(-2.4801319) q[0];
sx q[0];
rz(2.3308603) q[0];
rz(0.73075378) q[1];
sx q[1];
rz(-1.0565051) q[1];
sx q[1];
rz(-0.96053851) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86772462) q[0];
sx q[0];
rz(-1.5698213) q[0];
sx q[0];
rz(-0.0033897059) q[0];
x q[1];
rz(2.075718) q[2];
sx q[2];
rz(-1.9532579) q[2];
sx q[2];
rz(1.7590211) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5170578) q[1];
sx q[1];
rz(-1.0427999) q[1];
sx q[1];
rz(2.8370884) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0455564) q[3];
sx q[3];
rz(-0.48431319) q[3];
sx q[3];
rz(-2.7698295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1824128) q[2];
sx q[2];
rz(-1.3591839) q[2];
sx q[2];
rz(1.6020927) q[2];
rz(-2.0942073) q[3];
sx q[3];
rz(-1.7644707) q[3];
sx q[3];
rz(-1.876095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3613116) q[0];
sx q[0];
rz(-0.86532101) q[0];
sx q[0];
rz(-0.58746946) q[0];
rz(-2.7777708) q[1];
sx q[1];
rz(-1.1452585) q[1];
sx q[1];
rz(-2.2344373) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.016257441) q[0];
sx q[0];
rz(-1.9210707) q[0];
sx q[0];
rz(-1.6566234) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47260471) q[2];
sx q[2];
rz(-1.9197004) q[2];
sx q[2];
rz(0.34502236) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.75395155) q[1];
sx q[1];
rz(-2.6389737) q[1];
sx q[1];
rz(2.6941264) q[1];
rz(-pi) q[2];
x q[2];
rz(1.206622) q[3];
sx q[3];
rz(-1.4910021) q[3];
sx q[3];
rz(1.4353115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5017447) q[2];
sx q[2];
rz(-0.088968337) q[2];
sx q[2];
rz(2.7196344) q[2];
rz(-3.1320599) q[3];
sx q[3];
rz(-1.5744753) q[3];
sx q[3];
rz(0.52614051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.52787732) q[0];
sx q[0];
rz(-1.8462702) q[0];
sx q[0];
rz(0.12311002) q[0];
rz(2.4404424) q[1];
sx q[1];
rz(-1.2260561) q[1];
sx q[1];
rz(-1.4969926) q[1];
rz(0.43177615) q[2];
sx q[2];
rz(-2.3974621) q[2];
sx q[2];
rz(1.5324788) q[2];
rz(1.9811859) q[3];
sx q[3];
rz(-1.3153362) q[3];
sx q[3];
rz(1.5462331) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
