OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.9353256) q[0];
sx q[0];
rz(-1.5201818) q[0];
sx q[0];
rz(-1.4186463) q[0];
rz(-3.0980134) q[1];
sx q[1];
rz(-1.8585304) q[1];
sx q[1];
rz(-0.15951523) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43082008) q[0];
sx q[0];
rz(-2.9732467) q[0];
sx q[0];
rz(-1.7414067) q[0];
rz(2.5297727) q[2];
sx q[2];
rz(-1.3587591) q[2];
sx q[2];
rz(-3.0481047) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.72580569) q[1];
sx q[1];
rz(-1.09607) q[1];
sx q[1];
rz(2.9863547) q[1];
x q[2];
rz(0.84955022) q[3];
sx q[3];
rz(-1.8462984) q[3];
sx q[3];
rz(1.7003789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.688545) q[2];
sx q[2];
rz(-1.9828372) q[2];
sx q[2];
rz(0.090252074) q[2];
rz(0.074617535) q[3];
sx q[3];
rz(-1.6736232) q[3];
sx q[3];
rz(-2.727865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2523044) q[0];
sx q[0];
rz(-1.9828718) q[0];
sx q[0];
rz(1.4537551) q[0];
rz(2.3989035) q[1];
sx q[1];
rz(-1.366726) q[1];
sx q[1];
rz(0.76505605) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9738282) q[0];
sx q[0];
rz(-0.59685464) q[0];
sx q[0];
rz(-3.0227227) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3397331) q[2];
sx q[2];
rz(-0.86881402) q[2];
sx q[2];
rz(0.02846708) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.504952) q[1];
sx q[1];
rz(-0.8209629) q[1];
sx q[1];
rz(2.0269462) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5935947) q[3];
sx q[3];
rz(-1.573085) q[3];
sx q[3];
rz(-1.8016004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.226977) q[2];
sx q[2];
rz(-1.0116297) q[2];
sx q[2];
rz(-1.6870618) q[2];
rz(2.4948273) q[3];
sx q[3];
rz(-2.5819467) q[3];
sx q[3];
rz(1.8618934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0525368) q[0];
sx q[0];
rz(-1.4707969) q[0];
sx q[0];
rz(-3.0964417) q[0];
rz(-1.9106983) q[1];
sx q[1];
rz(-1.3094333) q[1];
sx q[1];
rz(-2.0848354) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8674744) q[0];
sx q[0];
rz(-1.8961268) q[0];
sx q[0];
rz(-2.1528457) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.64521328) q[2];
sx q[2];
rz(-1.1629379) q[2];
sx q[2];
rz(-2.4347664) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5910198) q[1];
sx q[1];
rz(-1.824784) q[1];
sx q[1];
rz(-2.2588782) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8523434) q[3];
sx q[3];
rz(-0.70812884) q[3];
sx q[3];
rz(0.81048465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7354273) q[2];
sx q[2];
rz(-2.1801517) q[2];
sx q[2];
rz(2.4889634) q[2];
rz(-1.8485707) q[3];
sx q[3];
rz(-2.3725464) q[3];
sx q[3];
rz(-0.098171083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80766135) q[0];
sx q[0];
rz(-0.46285358) q[0];
sx q[0];
rz(-1.7474784) q[0];
rz(1.8380503) q[1];
sx q[1];
rz(-1.9775672) q[1];
sx q[1];
rz(1.0882783) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1980285) q[0];
sx q[0];
rz(-1.8450415) q[0];
sx q[0];
rz(-2.345577) q[0];
rz(-0.91371961) q[2];
sx q[2];
rz(-0.83763257) q[2];
sx q[2];
rz(-2.4108771) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7409121) q[1];
sx q[1];
rz(-2.0517382) q[1];
sx q[1];
rz(-2.4467642) q[1];
x q[2];
rz(-0.49934629) q[3];
sx q[3];
rz(-1.5245228) q[3];
sx q[3];
rz(2.3114227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3950222) q[2];
sx q[2];
rz(-1.7638548) q[2];
sx q[2];
rz(2.6611888) q[2];
rz(0.26614842) q[3];
sx q[3];
rz(-1.1689309) q[3];
sx q[3];
rz(-0.84421414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10923037) q[0];
sx q[0];
rz(-0.6210331) q[0];
sx q[0];
rz(1.2048298) q[0];
rz(-0.042639848) q[1];
sx q[1];
rz(-1.7667021) q[1];
sx q[1];
rz(-0.94243324) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5980127) q[0];
sx q[0];
rz(-0.44853285) q[0];
sx q[0];
rz(-1.2537812) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3187014) q[2];
sx q[2];
rz(-2.3510755) q[2];
sx q[2];
rz(-2.2728777) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9246848) q[1];
sx q[1];
rz(-2.496301) q[1];
sx q[1];
rz(2.2297736) q[1];
x q[2];
rz(0.63905255) q[3];
sx q[3];
rz(-1.4830695) q[3];
sx q[3];
rz(1.7114342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0254536) q[2];
sx q[2];
rz(-2.4567273) q[2];
sx q[2];
rz(-1.2241414) q[2];
rz(-0.72743607) q[3];
sx q[3];
rz(-1.4801315) q[3];
sx q[3];
rz(-3.091403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8023476) q[0];
sx q[0];
rz(-2.6363063) q[0];
sx q[0];
rz(-2.8832866) q[0];
rz(-0.91336617) q[1];
sx q[1];
rz(-2.2759627) q[1];
sx q[1];
rz(0.9679274) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3334078) q[0];
sx q[0];
rz(-0.70751941) q[0];
sx q[0];
rz(2.524441) q[0];
rz(-1.9759693) q[2];
sx q[2];
rz(-2.3771667) q[2];
sx q[2];
rz(-2.7583099) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.55019864) q[1];
sx q[1];
rz(-2.0481143) q[1];
sx q[1];
rz(1.5309019) q[1];
rz(-2.4658396) q[3];
sx q[3];
rz(-0.76774137) q[3];
sx q[3];
rz(-2.2217285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4919058) q[2];
sx q[2];
rz(-2.457149) q[2];
sx q[2];
rz(0.41927949) q[2];
rz(-2.3573917) q[3];
sx q[3];
rz(-2.1214285) q[3];
sx q[3];
rz(-1.7694582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(3.075031) q[0];
sx q[0];
rz(-3.0176268) q[0];
sx q[0];
rz(-0.01874622) q[0];
rz(3.0491414) q[1];
sx q[1];
rz(-1.3321313) q[1];
sx q[1];
rz(-3.0466373) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0844011) q[0];
sx q[0];
rz(-2.4782628) q[0];
sx q[0];
rz(-0.54121547) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5592741) q[2];
sx q[2];
rz(-2.4597801) q[2];
sx q[2];
rz(2.5043629) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.590789) q[1];
sx q[1];
rz(-1.1491286) q[1];
sx q[1];
rz(1.292839) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9704814) q[3];
sx q[3];
rz(-1.7198472) q[3];
sx q[3];
rz(0.51250729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2350601) q[2];
sx q[2];
rz(-1.9573213) q[2];
sx q[2];
rz(2.2825784) q[2];
rz(0.61339316) q[3];
sx q[3];
rz(-2.9412013) q[3];
sx q[3];
rz(-1.332351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5253946) q[0];
sx q[0];
rz(-1.2735676) q[0];
sx q[0];
rz(-1.6599576) q[0];
rz(2.0741277) q[1];
sx q[1];
rz(-2.4165202) q[1];
sx q[1];
rz(-1.3190837) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7049371) q[0];
sx q[0];
rz(-2.4709513) q[0];
sx q[0];
rz(1.4192307) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.51609938) q[2];
sx q[2];
rz(-1.2428428) q[2];
sx q[2];
rz(-0.60521561) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6665277) q[1];
sx q[1];
rz(-1.3366146) q[1];
sx q[1];
rz(-0.98824309) q[1];
rz(-pi) q[2];
x q[2];
rz(2.784801) q[3];
sx q[3];
rz(-1.6902552) q[3];
sx q[3];
rz(0.91537634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8760425) q[2];
sx q[2];
rz(-0.37028131) q[2];
sx q[2];
rz(0.041157095) q[2];
rz(-2.0187078) q[3];
sx q[3];
rz(-1.4833781) q[3];
sx q[3];
rz(-0.52367228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6599643) q[0];
sx q[0];
rz(-2.1631503) q[0];
sx q[0];
rz(1.6690669) q[0];
rz(-0.72498471) q[1];
sx q[1];
rz(-1.9703194) q[1];
sx q[1];
rz(-2.5670126) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.816088) q[0];
sx q[0];
rz(-2.8003152) q[0];
sx q[0];
rz(-2.9983872) q[0];
x q[1];
rz(-2.8594744) q[2];
sx q[2];
rz(-1.1638255) q[2];
sx q[2];
rz(-2.5438683) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2255937) q[1];
sx q[1];
rz(-2.288842) q[1];
sx q[1];
rz(-1.6101933) q[1];
rz(2.9933917) q[3];
sx q[3];
rz(-2.253014) q[3];
sx q[3];
rz(-2.2489376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.64728111) q[2];
sx q[2];
rz(-0.24792555) q[2];
sx q[2];
rz(-2.1809123) q[2];
rz(-0.48323768) q[3];
sx q[3];
rz(-1.5471231) q[3];
sx q[3];
rz(-1.5620935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1305337) q[0];
sx q[0];
rz(-0.49376765) q[0];
sx q[0];
rz(-0.44179398) q[0];
rz(-0.68253016) q[1];
sx q[1];
rz(-1.6923169) q[1];
sx q[1];
rz(-2.0880879) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5614612) q[0];
sx q[0];
rz(-2.2765358) q[0];
sx q[0];
rz(0.34039008) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3467265) q[2];
sx q[2];
rz(-2.7599979) q[2];
sx q[2];
rz(-0.14118186) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7451585) q[1];
sx q[1];
rz(-0.72275439) q[1];
sx q[1];
rz(2.0682343) q[1];
x q[2];
rz(0.39393306) q[3];
sx q[3];
rz(-1.5026758) q[3];
sx q[3];
rz(1.7914655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0423476) q[2];
sx q[2];
rz(-1.2688364) q[2];
sx q[2];
rz(-2.919) q[2];
rz(-1.3709204) q[3];
sx q[3];
rz(-1.7226847) q[3];
sx q[3];
rz(0.62209779) q[3];
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
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2400146) q[0];
sx q[0];
rz(-1.0116901) q[0];
sx q[0];
rz(-2.1160175) q[0];
rz(-1.1657794) q[1];
sx q[1];
rz(-2.5318601) q[1];
sx q[1];
rz(-0.21808521) q[1];
rz(2.4161584) q[2];
sx q[2];
rz(-2.465783) q[2];
sx q[2];
rz(-3.0259544) q[2];
rz(1.5722269) q[3];
sx q[3];
rz(-0.32162255) q[3];
sx q[3];
rz(-2.376775) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
