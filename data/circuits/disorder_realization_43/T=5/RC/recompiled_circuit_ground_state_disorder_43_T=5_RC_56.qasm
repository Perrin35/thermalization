OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2744098) q[0];
sx q[0];
rz(2.6743968) q[0];
sx q[0];
rz(13.462486) q[0];
rz(2.3236302) q[1];
sx q[1];
rz(-1.5939413) q[1];
sx q[1];
rz(0.98734468) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5817669) q[0];
sx q[0];
rz(-0.94115138) q[0];
sx q[0];
rz(-2.6658076) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1842264) q[2];
sx q[2];
rz(-2.2920458) q[2];
sx q[2];
rz(2.3867875) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7973229) q[1];
sx q[1];
rz(-1.2943648) q[1];
sx q[1];
rz(-1.396342) q[1];
x q[2];
rz(1.4127068) q[3];
sx q[3];
rz(-1.2088547) q[3];
sx q[3];
rz(0.87495041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2200615) q[2];
sx q[2];
rz(-2.107373) q[2];
sx q[2];
rz(-0.49911487) q[2];
rz(0.81579298) q[3];
sx q[3];
rz(-2.54839) q[3];
sx q[3];
rz(-2.7844586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65207425) q[0];
sx q[0];
rz(-0.3638142) q[0];
sx q[0];
rz(0.2336842) q[0];
rz(3.0417327) q[1];
sx q[1];
rz(-1.654511) q[1];
sx q[1];
rz(-1.608009) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31029345) q[0];
sx q[0];
rz(-1.4742032) q[0];
sx q[0];
rz(-1.6929168) q[0];
rz(-pi) q[1];
x q[1];
rz(0.38012114) q[2];
sx q[2];
rz(-0.46762662) q[2];
sx q[2];
rz(-1.9531701) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.69791693) q[1];
sx q[1];
rz(-1.6622826) q[1];
sx q[1];
rz(2.1771624) q[1];
rz(2.1183999) q[3];
sx q[3];
rz(-0.55063081) q[3];
sx q[3];
rz(0.44704244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0443772) q[2];
sx q[2];
rz(-2.2558687) q[2];
sx q[2];
rz(-0.054917939) q[2];
rz(2.4954259) q[3];
sx q[3];
rz(-1.61444) q[3];
sx q[3];
rz(0.072877876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9570479) q[0];
sx q[0];
rz(-2.8962729) q[0];
sx q[0];
rz(2.3967337) q[0];
rz(-1.8648196) q[1];
sx q[1];
rz(-1.0294015) q[1];
sx q[1];
rz(-0.863711) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.058216393) q[0];
sx q[0];
rz(-1.5124707) q[0];
sx q[0];
rz(-2.7013679) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0625819) q[2];
sx q[2];
rz(-0.84948601) q[2];
sx q[2];
rz(-0.53141293) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3491213) q[1];
sx q[1];
rz(-2.0466261) q[1];
sx q[1];
rz(-0.70648593) q[1];
x q[2];
rz(-2.5287348) q[3];
sx q[3];
rz(-2.3817059) q[3];
sx q[3];
rz(2.5155544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.471571) q[2];
sx q[2];
rz(-1.5315285) q[2];
sx q[2];
rz(2.7460597) q[2];
rz(2.1192571) q[3];
sx q[3];
rz(-0.46415713) q[3];
sx q[3];
rz(2.6497604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1763497) q[0];
sx q[0];
rz(-1.1818161) q[0];
sx q[0];
rz(3.0856207) q[0];
rz(2.5296027) q[1];
sx q[1];
rz(-1.5925708) q[1];
sx q[1];
rz(2.2956119) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6406785) q[0];
sx q[0];
rz(-2.0862824) q[0];
sx q[0];
rz(0.2978931) q[0];
rz(-pi) q[1];
rz(0.34807713) q[2];
sx q[2];
rz(-0.98053369) q[2];
sx q[2];
rz(1.8463617) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.31797925) q[1];
sx q[1];
rz(-0.3682963) q[1];
sx q[1];
rz(0.55693926) q[1];
x q[2];
rz(-2.0457129) q[3];
sx q[3];
rz(-1.9088512) q[3];
sx q[3];
rz(-0.89382225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1299639) q[2];
sx q[2];
rz(-1.446898) q[2];
sx q[2];
rz(2.4212627) q[2];
rz(-2.7885041) q[3];
sx q[3];
rz(-1.947764) q[3];
sx q[3];
rz(0.24875719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7406152) q[0];
sx q[0];
rz(-1.4530797) q[0];
sx q[0];
rz(0.98689669) q[0];
rz(1.997021) q[1];
sx q[1];
rz(-0.83895504) q[1];
sx q[1];
rz(-0.95485895) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.389176) q[0];
sx q[0];
rz(-1.5446394) q[0];
sx q[0];
rz(0.90361066) q[0];
rz(-pi) q[1];
rz(2.0667162) q[2];
sx q[2];
rz(-0.62037797) q[2];
sx q[2];
rz(-1.7091144) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.92063078) q[1];
sx q[1];
rz(-2.536533) q[1];
sx q[1];
rz(-0.69843881) q[1];
x q[2];
rz(1.8243039) q[3];
sx q[3];
rz(-1.714332) q[3];
sx q[3];
rz(1.1689127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5130875) q[2];
sx q[2];
rz(-2.2272765) q[2];
sx q[2];
rz(-1.0465735) q[2];
rz(-0.70932499) q[3];
sx q[3];
rz(-1.9966634) q[3];
sx q[3];
rz(-2.5078702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3777622) q[0];
sx q[0];
rz(-1.7338294) q[0];
sx q[0];
rz(2.8705226) q[0];
rz(-0.079004869) q[1];
sx q[1];
rz(-0.43402356) q[1];
sx q[1];
rz(1.431042) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97320172) q[0];
sx q[0];
rz(-1.7152837) q[0];
sx q[0];
rz(2.9299749) q[0];
x q[1];
rz(0.56366326) q[2];
sx q[2];
rz(-1.0978292) q[2];
sx q[2];
rz(-1.2207292) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8003098) q[1];
sx q[1];
rz(-0.73004111) q[1];
sx q[1];
rz(-2.7679539) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9892788) q[3];
sx q[3];
rz(-1.4238723) q[3];
sx q[3];
rz(0.72673405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.048910353) q[2];
sx q[2];
rz(-1.9581257) q[2];
sx q[2];
rz(-0.023822039) q[2];
rz(1.9134936) q[3];
sx q[3];
rz(-0.3796328) q[3];
sx q[3];
rz(-0.14321271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.989885) q[0];
sx q[0];
rz(-0.30549529) q[0];
sx q[0];
rz(0.09224961) q[0];
rz(2.8406738) q[1];
sx q[1];
rz(-1.9124799) q[1];
sx q[1];
rz(-1.0801962) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5841519) q[0];
sx q[0];
rz(-1.6655465) q[0];
sx q[0];
rz(1.2242203) q[0];
rz(-pi) q[1];
rz(-1.9930196) q[2];
sx q[2];
rz(-1.9034837) q[2];
sx q[2];
rz(-2.0300558) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2094592) q[1];
sx q[1];
rz(-1.5725296) q[1];
sx q[1];
rz(-0.14793795) q[1];
rz(0.37410847) q[3];
sx q[3];
rz(-0.7302993) q[3];
sx q[3];
rz(0.41159831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.13176189) q[2];
sx q[2];
rz(-2.7907382) q[2];
sx q[2];
rz(-0.62823137) q[2];
rz(2.175711) q[3];
sx q[3];
rz(-1.5647669) q[3];
sx q[3];
rz(-2.8204744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1107165) q[0];
sx q[0];
rz(-0.70699152) q[0];
sx q[0];
rz(1.734717) q[0];
rz(-0.12786099) q[1];
sx q[1];
rz(-1.4297994) q[1];
sx q[1];
rz(1.1258639) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0681852) q[0];
sx q[0];
rz(-0.74755423) q[0];
sx q[0];
rz(-2.8697467) q[0];
x q[1];
rz(1.5324513) q[2];
sx q[2];
rz(-0.45203096) q[2];
sx q[2];
rz(-2.0386809) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8496823) q[1];
sx q[1];
rz(-2.9314802) q[1];
sx q[1];
rz(-0.30847524) q[1];
rz(-0.58784318) q[3];
sx q[3];
rz(-0.70826642) q[3];
sx q[3];
rz(-1.721134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.13059482) q[2];
sx q[2];
rz(-1.1329634) q[2];
sx q[2];
rz(3.0478743) q[2];
rz(-1.4334009) q[3];
sx q[3];
rz(-2.8580557) q[3];
sx q[3];
rz(1.3933498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4505287) q[0];
sx q[0];
rz(-2.0557623) q[0];
sx q[0];
rz(-1.0040671) q[0];
rz(-1.1035236) q[1];
sx q[1];
rz(-0.48201489) q[1];
sx q[1];
rz(1.4080661) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.133014) q[0];
sx q[0];
rz(-2.6264418) q[0];
sx q[0];
rz(-3.1293082) q[0];
rz(-1.682071) q[2];
sx q[2];
rz(-1.2603052) q[2];
sx q[2];
rz(1.3729915) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1939331) q[1];
sx q[1];
rz(-1.1378189) q[1];
sx q[1];
rz(2.9440332) q[1];
x q[2];
rz(-2.6226165) q[3];
sx q[3];
rz(-1.5277315) q[3];
sx q[3];
rz(-1.2534634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7496926) q[2];
sx q[2];
rz(-1.5211279) q[2];
sx q[2];
rz(-2.5938972) q[2];
rz(-0.31217602) q[3];
sx q[3];
rz(-1.7269937) q[3];
sx q[3];
rz(1.7267905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3373229) q[0];
sx q[0];
rz(-1.6707358) q[0];
sx q[0];
rz(-2.3626589) q[0];
rz(0.3745105) q[1];
sx q[1];
rz(-1.6284527) q[1];
sx q[1];
rz(2.3096854) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4171281) q[0];
sx q[0];
rz(-0.75675557) q[0];
sx q[0];
rz(-2.7555076) q[0];
x q[1];
rz(-0.7124165) q[2];
sx q[2];
rz(-2.2210026) q[2];
sx q[2];
rz(-2.4622816) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.9495957) q[1];
sx q[1];
rz(-1.0284541) q[1];
sx q[1];
rz(2.3344912) q[1];
rz(-1.5396773) q[3];
sx q[3];
rz(-0.77431576) q[3];
sx q[3];
rz(-1.4207135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1067918) q[2];
sx q[2];
rz(-2.782244) q[2];
sx q[2];
rz(0.0084776004) q[2];
rz(-0.40555412) q[3];
sx q[3];
rz(-1.4529934) q[3];
sx q[3];
rz(-1.4439553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1494898) q[0];
sx q[0];
rz(-1.533951) q[0];
sx q[0];
rz(0.78716192) q[0];
rz(-1.0943195) q[1];
sx q[1];
rz(-0.37847395) q[1];
sx q[1];
rz(2.7056221) q[1];
rz(1.0395861) q[2];
sx q[2];
rz(-0.24460228) q[2];
sx q[2];
rz(-1.2096418) q[2];
rz(1.575749) q[3];
sx q[3];
rz(-1.6650143) q[3];
sx q[3];
rz(-3.0190013) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
