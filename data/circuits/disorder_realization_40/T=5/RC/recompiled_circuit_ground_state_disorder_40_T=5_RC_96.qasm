OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.52387828) q[0];
sx q[0];
rz(3.7359306) q[0];
sx q[0];
rz(9.115968) q[0];
rz(-2.8438957) q[1];
sx q[1];
rz(-1.2574137) q[1];
sx q[1];
rz(2.5296192) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64498211) q[0];
sx q[0];
rz(-2.5685446) q[0];
sx q[0];
rz(1.7368421) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0711028) q[2];
sx q[2];
rz(-0.32499748) q[2];
sx q[2];
rz(0.18735841) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6943629) q[1];
sx q[1];
rz(-1.9026533) q[1];
sx q[1];
rz(2.993078) q[1];
rz(-pi) q[2];
rz(0.6925005) q[3];
sx q[3];
rz(-1.6095265) q[3];
sx q[3];
rz(0.2693143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1674126) q[2];
sx q[2];
rz(-1.9638502) q[2];
sx q[2];
rz(-0.56498945) q[2];
rz(2.6307093) q[3];
sx q[3];
rz(-2.9027945) q[3];
sx q[3];
rz(1.6142982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6086455) q[0];
sx q[0];
rz(-0.2810418) q[0];
sx q[0];
rz(0.99579048) q[0];
rz(-0.58798724) q[1];
sx q[1];
rz(-0.43991393) q[1];
sx q[1];
rz(-1.8221375) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3036766) q[0];
sx q[0];
rz(-0.48978031) q[0];
sx q[0];
rz(1.0465001) q[0];
rz(0.075602268) q[2];
sx q[2];
rz(-1.2443064) q[2];
sx q[2];
rz(-1.16815) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7403789) q[1];
sx q[1];
rz(-1.7771562) q[1];
sx q[1];
rz(-2.3310082) q[1];
rz(-2.6590517) q[3];
sx q[3];
rz(-1.6785033) q[3];
sx q[3];
rz(1.5272702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1002645) q[2];
sx q[2];
rz(-2.1375956) q[2];
sx q[2];
rz(-0.022493258) q[2];
rz(3.0371173) q[3];
sx q[3];
rz(-1.6040809) q[3];
sx q[3];
rz(0.87048602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.81258881) q[0];
sx q[0];
rz(-2.1109695) q[0];
sx q[0];
rz(-1.6382244) q[0];
rz(-0.080987856) q[1];
sx q[1];
rz(-2.4826725) q[1];
sx q[1];
rz(1.9714877) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92540298) q[0];
sx q[0];
rz(-1.7316375) q[0];
sx q[0];
rz(0.58165929) q[0];
x q[1];
rz(-0.64441893) q[2];
sx q[2];
rz(-1.7303409) q[2];
sx q[2];
rz(1.5089515) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2418211) q[1];
sx q[1];
rz(-1.542143) q[1];
sx q[1];
rz(-1.0706399) q[1];
x q[2];
rz(-2.0017573) q[3];
sx q[3];
rz(-0.54910481) q[3];
sx q[3];
rz(0.74326754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.66204232) q[2];
sx q[2];
rz(-0.68632555) q[2];
sx q[2];
rz(0.043206841) q[2];
rz(-2.9442545) q[3];
sx q[3];
rz(-1.03136) q[3];
sx q[3];
rz(-2.7092547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2546805) q[0];
sx q[0];
rz(-1.0858902) q[0];
sx q[0];
rz(0.32549724) q[0];
rz(0.85572851) q[1];
sx q[1];
rz(-0.41494644) q[1];
sx q[1];
rz(-1.0861446) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7522858) q[0];
sx q[0];
rz(-2.6323937) q[0];
sx q[0];
rz(0.26637116) q[0];
rz(-pi) q[1];
rz(2.3620783) q[2];
sx q[2];
rz(-0.91060591) q[2];
sx q[2];
rz(-2.9110031) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.74198308) q[1];
sx q[1];
rz(-2.0860414) q[1];
sx q[1];
rz(-2.9523729) q[1];
rz(-1.7138163) q[3];
sx q[3];
rz(-1.656683) q[3];
sx q[3];
rz(-0.68616223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.93422741) q[2];
sx q[2];
rz(-0.97565979) q[2];
sx q[2];
rz(0.18386851) q[2];
rz(-1.3217226) q[3];
sx q[3];
rz(-0.70122856) q[3];
sx q[3];
rz(-0.13300657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70165271) q[0];
sx q[0];
rz(-2.8337605) q[0];
sx q[0];
rz(-0.93691784) q[0];
rz(-2.2266455) q[1];
sx q[1];
rz(-2.2888384) q[1];
sx q[1];
rz(1.459704) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5640083) q[0];
sx q[0];
rz(-1.47808) q[0];
sx q[0];
rz(-0.87323879) q[0];
rz(0.95090805) q[2];
sx q[2];
rz(-2.0712598) q[2];
sx q[2];
rz(-1.9096979) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8847757) q[1];
sx q[1];
rz(-2.3814572) q[1];
sx q[1];
rz(1.4291309) q[1];
x q[2];
rz(-2.567148) q[3];
sx q[3];
rz(-0.46009053) q[3];
sx q[3];
rz(-2.121832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0978284) q[2];
sx q[2];
rz(-1.3358668) q[2];
sx q[2];
rz(2.9975927) q[2];
rz(-2.0740267) q[3];
sx q[3];
rz(-0.34393603) q[3];
sx q[3];
rz(2.4501154) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3095793) q[0];
sx q[0];
rz(-0.80736512) q[0];
sx q[0];
rz(2.373234) q[0];
rz(-2.2840624) q[1];
sx q[1];
rz(-1.0464959) q[1];
sx q[1];
rz(2.5879477) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2087649) q[0];
sx q[0];
rz(-2.0220482) q[0];
sx q[0];
rz(2.8594467) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3403999) q[2];
sx q[2];
rz(-1.4095708) q[2];
sx q[2];
rz(0.025242199) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.91398134) q[1];
sx q[1];
rz(-2.1030594) q[1];
sx q[1];
rz(2.6460939) q[1];
x q[2];
rz(-0.97691128) q[3];
sx q[3];
rz(-2.2008228) q[3];
sx q[3];
rz(-1.6225369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7553317) q[2];
sx q[2];
rz(-0.43060455) q[2];
sx q[2];
rz(0.44580305) q[2];
rz(1.2465994) q[3];
sx q[3];
rz(-1.3695025) q[3];
sx q[3];
rz(2.2915452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.069139473) q[0];
sx q[0];
rz(-2.1816165) q[0];
sx q[0];
rz(3.0015216) q[0];
rz(2.2135997) q[1];
sx q[1];
rz(-1.3368139) q[1];
sx q[1];
rz(1.6915406) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1221913) q[0];
sx q[0];
rz(-0.1941351) q[0];
sx q[0];
rz(-0.96786626) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23132236) q[2];
sx q[2];
rz(-1.9440821) q[2];
sx q[2];
rz(1.1333677) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.90277744) q[1];
sx q[1];
rz(-2.0014781) q[1];
sx q[1];
rz(2.4611453) q[1];
rz(-2.9599056) q[3];
sx q[3];
rz(-1.7227001) q[3];
sx q[3];
rz(-3.0345598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.037584) q[2];
sx q[2];
rz(-0.17825492) q[2];
sx q[2];
rz(-0.083871052) q[2];
rz(0.9907848) q[3];
sx q[3];
rz(-1.1889941) q[3];
sx q[3];
rz(-1.8535463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.8952119) q[0];
sx q[0];
rz(-1.8354494) q[0];
sx q[0];
rz(2.7847248) q[0];
rz(-1.6995947) q[1];
sx q[1];
rz(-2.5394963) q[1];
sx q[1];
rz(3.1325565) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5166695) q[0];
sx q[0];
rz(-2.8213266) q[0];
sx q[0];
rz(2.3413168) q[0];
rz(1.1844877) q[2];
sx q[2];
rz(-2.2335459) q[2];
sx q[2];
rz(0.3370924) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9408843) q[1];
sx q[1];
rz(-1.4794596) q[1];
sx q[1];
rz(-1.2979085) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5311386) q[3];
sx q[3];
rz(-2.0137827) q[3];
sx q[3];
rz(-2.7970527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.25702) q[2];
sx q[2];
rz(-1.7947861) q[2];
sx q[2];
rz(-0.70551562) q[2];
rz(0.26933119) q[3];
sx q[3];
rz(-2.0651385) q[3];
sx q[3];
rz(-0.96946019) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89328289) q[0];
sx q[0];
rz(-0.80535424) q[0];
sx q[0];
rz(2.4321108) q[0];
rz(2.4995038) q[1];
sx q[1];
rz(-2.700192) q[1];
sx q[1];
rz(0.4253687) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.294153) q[0];
sx q[0];
rz(-0.80753875) q[0];
sx q[0];
rz(-0.60955131) q[0];
x q[1];
rz(-0.9832731) q[2];
sx q[2];
rz(-2.4725178) q[2];
sx q[2];
rz(-2.4226505) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.44008419) q[1];
sx q[1];
rz(-0.45508859) q[1];
sx q[1];
rz(0.20642682) q[1];
x q[2];
rz(2.5675699) q[3];
sx q[3];
rz(-1.8282601) q[3];
sx q[3];
rz(-2.4935745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8713259) q[2];
sx q[2];
rz(-0.17921236) q[2];
sx q[2];
rz(-0.47879177) q[2];
rz(-2.791413) q[3];
sx q[3];
rz(-1.9637354) q[3];
sx q[3];
rz(-2.71463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66788524) q[0];
sx q[0];
rz(-2.844664) q[0];
sx q[0];
rz(-0.83734751) q[0];
rz(0.26652023) q[1];
sx q[1];
rz(-1.3198677) q[1];
sx q[1];
rz(-1.3574319) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0650554) q[0];
sx q[0];
rz(-1.6572857) q[0];
sx q[0];
rz(3.1381395) q[0];
x q[1];
rz(-0.063284782) q[2];
sx q[2];
rz(-1.4522168) q[2];
sx q[2];
rz(2.0381387) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3597922) q[1];
sx q[1];
rz(-2.4891653) q[1];
sx q[1];
rz(3.1210207) q[1];
rz(-pi) q[2];
x q[2];
rz(0.97748791) q[3];
sx q[3];
rz(-2.5989669) q[3];
sx q[3];
rz(0.092433905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.16167656) q[2];
sx q[2];
rz(-0.53407532) q[2];
sx q[2];
rz(-2.7034289) q[2];
rz(-0.91580373) q[3];
sx q[3];
rz(-1.6180429) q[3];
sx q[3];
rz(0.80317909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.5220779) q[0];
sx q[0];
rz(-1.8997471) q[0];
sx q[0];
rz(2.1258623) q[0];
rz(0.10454128) q[1];
sx q[1];
rz(-1.8395945) q[1];
sx q[1];
rz(1.3855388) q[1];
rz(0.06212297) q[2];
sx q[2];
rz(-1.9005513) q[2];
sx q[2];
rz(-1.8951436) q[2];
rz(-0.87540353) q[3];
sx q[3];
rz(-0.84423748) q[3];
sx q[3];
rz(0.68035833) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
