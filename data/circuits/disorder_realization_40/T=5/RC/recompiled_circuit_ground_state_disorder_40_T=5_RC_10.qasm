OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6177144) q[0];
sx q[0];
rz(-0.59433794) q[0];
sx q[0];
rz(-2.8327827) q[0];
rz(-2.8438957) q[1];
sx q[1];
rz(-1.2574137) q[1];
sx q[1];
rz(2.5296192) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4966105) q[0];
sx q[0];
rz(-2.5685446) q[0];
sx q[0];
rz(-1.7368421) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9813479) q[2];
sx q[2];
rz(-1.2868243) q[2];
sx q[2];
rz(-0.71061963) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2638611) q[1];
sx q[1];
rz(-2.7791443) q[1];
sx q[1];
rz(1.1652105) q[1];
rz(0.060617491) q[3];
sx q[3];
rz(-2.4481886) q[3];
sx q[3];
rz(-1.3481247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1674126) q[2];
sx q[2];
rz(-1.9638502) q[2];
sx q[2];
rz(-2.5766032) q[2];
rz(2.6307093) q[3];
sx q[3];
rz(-0.23879819) q[3];
sx q[3];
rz(1.5272944) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6086455) q[0];
sx q[0];
rz(-2.8605509) q[0];
sx q[0];
rz(-0.99579048) q[0];
rz(0.58798724) q[1];
sx q[1];
rz(-0.43991393) q[1];
sx q[1];
rz(1.8221375) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72356862) q[0];
sx q[0];
rz(-1.9902285) q[0];
sx q[0];
rz(2.8807958) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3513279) q[2];
sx q[2];
rz(-0.33482345) q[2];
sx q[2];
rz(1.741516) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.118757) q[1];
sx q[1];
rz(-0.83056994) q[1];
sx q[1];
rz(2.8603795) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.22892816) q[3];
sx q[3];
rz(-0.49348885) q[3];
sx q[3];
rz(2.8956653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1002645) q[2];
sx q[2];
rz(-2.1375956) q[2];
sx q[2];
rz(-3.1190994) q[2];
rz(-0.10447539) q[3];
sx q[3];
rz(-1.6040809) q[3];
sx q[3];
rz(-2.2711066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81258881) q[0];
sx q[0];
rz(-1.0306232) q[0];
sx q[0];
rz(-1.5033683) q[0];
rz(-0.080987856) q[1];
sx q[1];
rz(-2.4826725) q[1];
sx q[1];
rz(1.9714877) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54047548) q[0];
sx q[0];
rz(-0.99760054) q[0];
sx q[0];
rz(1.3790087) q[0];
rz(-2.8798772) q[2];
sx q[2];
rz(-0.66113512) q[2];
sx q[2];
rz(-2.9950855) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8997716) q[1];
sx q[1];
rz(-1.5994497) q[1];
sx q[1];
rz(2.0709527) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0634138) q[3];
sx q[3];
rz(-1.3510002) q[3];
sx q[3];
rz(-2.6878217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4795503) q[2];
sx q[2];
rz(-2.4552671) q[2];
sx q[2];
rz(0.043206841) q[2];
rz(-0.19733812) q[3];
sx q[3];
rz(-1.03136) q[3];
sx q[3];
rz(2.7092547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8869121) q[0];
sx q[0];
rz(-2.0557025) q[0];
sx q[0];
rz(2.8160954) q[0];
rz(2.2858641) q[1];
sx q[1];
rz(-0.41494644) q[1];
sx q[1];
rz(1.0861446) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.08642149) q[0];
sx q[0];
rz(-1.0811816) q[0];
sx q[0];
rz(1.7167164) q[0];
x q[1];
rz(-0.77951435) q[2];
sx q[2];
rz(-0.91060591) q[2];
sx q[2];
rz(0.23058952) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3996096) q[1];
sx q[1];
rz(-2.0860414) q[1];
sx q[1];
rz(-2.9523729) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4277763) q[3];
sx q[3];
rz(-1.4849097) q[3];
sx q[3];
rz(2.4554304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2073652) q[2];
sx q[2];
rz(-2.1659329) q[2];
sx q[2];
rz(0.18386851) q[2];
rz(-1.3217226) q[3];
sx q[3];
rz(-2.4403641) q[3];
sx q[3];
rz(-3.0085861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70165271) q[0];
sx q[0];
rz(-0.30783215) q[0];
sx q[0];
rz(-2.2046748) q[0];
rz(-2.2266455) q[1];
sx q[1];
rz(-2.2888384) q[1];
sx q[1];
rz(-1.6818887) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5640083) q[0];
sx q[0];
rz(-1.47808) q[0];
sx q[0];
rz(0.87323879) q[0];
rz(-pi) q[1];
rz(2.5499623) q[2];
sx q[2];
rz(-2.1056294) q[2];
sx q[2];
rz(-2.4727351) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8847757) q[1];
sx q[1];
rz(-2.3814572) q[1];
sx q[1];
rz(-1.4291309) q[1];
rz(1.3077626) q[3];
sx q[3];
rz(-1.9527779) q[3];
sx q[3];
rz(-0.39417496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0437643) q[2];
sx q[2];
rz(-1.3358668) q[2];
sx q[2];
rz(-2.9975927) q[2];
rz(-1.0675659) q[3];
sx q[3];
rz(-0.34393603) q[3];
sx q[3];
rz(-2.4501154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-1.8320134) q[0];
sx q[0];
rz(-2.3342275) q[0];
sx q[0];
rz(0.76835865) q[0];
rz(0.8575303) q[1];
sx q[1];
rz(-1.0464959) q[1];
sx q[1];
rz(-0.55364496) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93282774) q[0];
sx q[0];
rz(-1.1195445) q[0];
sx q[0];
rz(-0.28214595) q[0];
rz(-pi) q[1];
rz(-2.3403999) q[2];
sx q[2];
rz(-1.7320219) q[2];
sx q[2];
rz(3.1163505) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2276113) q[1];
sx q[1];
rz(-1.0385333) q[1];
sx q[1];
rz(0.49549876) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.72154273) q[3];
sx q[3];
rz(-1.1016204) q[3];
sx q[3];
rz(2.814722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.386261) q[2];
sx q[2];
rz(-2.7109881) q[2];
sx q[2];
rz(-0.44580305) q[2];
rz(1.8949932) q[3];
sx q[3];
rz(-1.7720902) q[3];
sx q[3];
rz(-0.8500475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.069139473) q[0];
sx q[0];
rz(-2.1816165) q[0];
sx q[0];
rz(3.0015216) q[0];
rz(-0.92799294) q[1];
sx q[1];
rz(-1.8047787) q[1];
sx q[1];
rz(-1.6915406) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0194013) q[0];
sx q[0];
rz(-0.1941351) q[0];
sx q[0];
rz(0.96786626) q[0];
x q[1];
rz(1.0411925) q[2];
sx q[2];
rz(-0.43627377) q[2];
sx q[2];
rz(-2.5817007) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2388152) q[1];
sx q[1];
rz(-1.1401145) q[1];
sx q[1];
rz(2.4611453) q[1];
x q[2];
rz(-2.9599056) q[3];
sx q[3];
rz(-1.7227001) q[3];
sx q[3];
rz(0.10703281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.037584) q[2];
sx q[2];
rz(-2.9633377) q[2];
sx q[2];
rz(-0.083871052) q[2];
rz(-2.1508079) q[3];
sx q[3];
rz(-1.1889941) q[3];
sx q[3];
rz(1.2880464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24638076) q[0];
sx q[0];
rz(-1.3061433) q[0];
sx q[0];
rz(-2.7847248) q[0];
rz(1.4419979) q[1];
sx q[1];
rz(-2.5394963) q[1];
sx q[1];
rz(3.1325565) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9401902) q[0];
sx q[0];
rz(-1.7918669) q[0];
sx q[0];
rz(1.3371435) q[0];
rz(-0.44972723) q[2];
sx q[2];
rz(-0.75216659) q[2];
sx q[2];
rz(-0.2470242) q[2];
sx q[3];
rz(-pi/2) q[3];
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
x q[2];
rz(1.5311386) q[3];
sx q[3];
rz(-1.1278099) q[3];
sx q[3];
rz(-2.7970527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8845727) q[2];
sx q[2];
rz(-1.7947861) q[2];
sx q[2];
rz(-0.70551562) q[2];
rz(-2.8722615) q[3];
sx q[3];
rz(-1.0764542) q[3];
sx q[3];
rz(0.96946019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89328289) q[0];
sx q[0];
rz(-2.3362384) q[0];
sx q[0];
rz(-2.4321108) q[0];
rz(-2.4995038) q[1];
sx q[1];
rz(-2.700192) q[1];
sx q[1];
rz(-0.4253687) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50373721) q[0];
sx q[0];
rz(-2.2048973) q[0];
sx q[0];
rz(1.0315328) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.413093) q[2];
sx q[2];
rz(-1.0283264) q[2];
sx q[2];
rz(1.7187207) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.21103745) q[1];
sx q[1];
rz(-2.0155219) q[1];
sx q[1];
rz(-1.4708323) q[1];
x q[2];
rz(-2.5675699) q[3];
sx q[3];
rz(-1.3133326) q[3];
sx q[3];
rz(-2.4935745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2702668) q[2];
sx q[2];
rz(-2.9623803) q[2];
sx q[2];
rz(-0.47879177) q[2];
rz(-2.791413) q[3];
sx q[3];
rz(-1.9637354) q[3];
sx q[3];
rz(0.42696264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4737074) q[0];
sx q[0];
rz(-2.844664) q[0];
sx q[0];
rz(2.3042451) q[0];
rz(-2.8750724) q[1];
sx q[1];
rz(-1.8217249) q[1];
sx q[1];
rz(-1.7841608) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11649179) q[0];
sx q[0];
rz(-3.0550346) q[0];
sx q[0];
rz(1.5309912) q[0];
x q[1];
rz(-3.0783079) q[2];
sx q[2];
rz(-1.6893759) q[2];
sx q[2];
rz(-1.103454) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7818004) q[1];
sx q[1];
rz(-2.4891653) q[1];
sx q[1];
rz(-3.1210207) q[1];
rz(-pi) q[2];
rz(0.97748791) q[3];
sx q[3];
rz(-2.5989669) q[3];
sx q[3];
rz(0.092433905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9799161) q[2];
sx q[2];
rz(-0.53407532) q[2];
sx q[2];
rz(-0.43816379) q[2];
rz(2.2257889) q[3];
sx q[3];
rz(-1.5235498) q[3];
sx q[3];
rz(2.3384136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61951471) q[0];
sx q[0];
rz(-1.8997471) q[0];
sx q[0];
rz(2.1258623) q[0];
rz(-3.0370514) q[1];
sx q[1];
rz(-1.8395945) q[1];
sx q[1];
rz(1.3855388) q[1];
rz(1.3913515) q[2];
sx q[2];
rz(-2.806247) q[2];
sx q[2];
rz(1.4362337) q[2];
rz(-0.62461169) q[3];
sx q[3];
rz(-0.95952989) q[3];
sx q[3];
rz(2.923514) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
