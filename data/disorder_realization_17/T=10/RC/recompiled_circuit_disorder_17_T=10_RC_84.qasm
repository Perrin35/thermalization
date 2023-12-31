OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2330115) q[0];
sx q[0];
rz(-1.1865948) q[0];
sx q[0];
rz(-1.9957805) q[0];
rz(0.97283483) q[1];
sx q[1];
rz(4.8117074) q[1];
sx q[1];
rz(9.7172578) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8669517) q[0];
sx q[0];
rz(-2.7434218) q[0];
sx q[0];
rz(1.7142332) q[0];
x q[1];
rz(-1.5195261) q[2];
sx q[2];
rz(-2.3495557) q[2];
sx q[2];
rz(-0.37026065) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8279743) q[1];
sx q[1];
rz(-0.59578124) q[1];
sx q[1];
rz(1.6525364) q[1];
rz(-pi) q[2];
rz(-1.5670033) q[3];
sx q[3];
rz(-0.47382254) q[3];
sx q[3];
rz(0.68236085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.73074377) q[2];
sx q[2];
rz(-2.0480672) q[2];
sx q[2];
rz(-0.4494108) q[2];
rz(-2.4959026) q[3];
sx q[3];
rz(-0.68104762) q[3];
sx q[3];
rz(1.8507563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21355024) q[0];
sx q[0];
rz(-2.2556861) q[0];
sx q[0];
rz(0.74203062) q[0];
rz(1.6702601) q[1];
sx q[1];
rz(-2.4447618) q[1];
sx q[1];
rz(-1.3630294) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78995401) q[0];
sx q[0];
rz(-0.75766701) q[0];
sx q[0];
rz(-0.58654465) q[0];
rz(-pi) q[1];
x q[1];
rz(0.46911932) q[2];
sx q[2];
rz(-1.3781147) q[2];
sx q[2];
rz(2.0872781) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.084761707) q[1];
sx q[1];
rz(-2.8014604) q[1];
sx q[1];
rz(3.0840193) q[1];
rz(2.9491049) q[3];
sx q[3];
rz(-1.8322332) q[3];
sx q[3];
rz(0.9769494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.30119511) q[2];
sx q[2];
rz(-2.6448554) q[2];
sx q[2];
rz(-1.2724686) q[2];
rz(-0.33723351) q[3];
sx q[3];
rz(-1.971258) q[3];
sx q[3];
rz(0.01005323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.067588016) q[0];
sx q[0];
rz(-0.33960605) q[0];
sx q[0];
rz(0.2628251) q[0];
rz(-2.7858531) q[1];
sx q[1];
rz(-1.6825312) q[1];
sx q[1];
rz(1.9062818) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10997406) q[0];
sx q[0];
rz(-1.5556766) q[0];
sx q[0];
rz(-1.581024) q[0];
rz(-pi) q[1];
rz(-3.1149408) q[2];
sx q[2];
rz(-1.8111374) q[2];
sx q[2];
rz(-1.196256) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3706627) q[1];
sx q[1];
rz(-1.2491033) q[1];
sx q[1];
rz(1.4598203) q[1];
rz(-pi) q[2];
x q[2];
rz(0.87576207) q[3];
sx q[3];
rz(-1.2998298) q[3];
sx q[3];
rz(0.7489487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1245023) q[2];
sx q[2];
rz(-1.3829145) q[2];
sx q[2];
rz(-2.9612605) q[2];
rz(-0.63594937) q[3];
sx q[3];
rz(-1.162581) q[3];
sx q[3];
rz(-0.37160555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76698774) q[0];
sx q[0];
rz(-0.42414442) q[0];
sx q[0];
rz(0.71907991) q[0];
rz(-2.6285697) q[1];
sx q[1];
rz(-1.9388371) q[1];
sx q[1];
rz(2.1069353) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35095222) q[0];
sx q[0];
rz(-1.5119023) q[0];
sx q[0];
rz(1.5509997) q[0];
rz(-0.047495202) q[2];
sx q[2];
rz(-1.4357899) q[2];
sx q[2];
rz(-2.6035655) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3116341) q[1];
sx q[1];
rz(-1.2882075) q[1];
sx q[1];
rz(-0.71210536) q[1];
x q[2];
rz(1.3961117) q[3];
sx q[3];
rz(-1.1038728) q[3];
sx q[3];
rz(-0.13388453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.93306142) q[2];
sx q[2];
rz(-0.26991093) q[2];
sx q[2];
rz(2.9053524) q[2];
rz(-1.158372) q[3];
sx q[3];
rz(-1.002243) q[3];
sx q[3];
rz(2.2896144) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.029595705) q[0];
sx q[0];
rz(-2.1313666) q[0];
sx q[0];
rz(-0.90233666) q[0];
rz(-2.2205655) q[1];
sx q[1];
rz(-0.61704707) q[1];
sx q[1];
rz(-0.62228084) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.92684) q[0];
sx q[0];
rz(-1.2878969) q[0];
sx q[0];
rz(-2.9515285) q[0];
rz(-pi) q[1];
rz(-0.12139608) q[2];
sx q[2];
rz(-2.4897794) q[2];
sx q[2];
rz(0.35169841) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8497148) q[1];
sx q[1];
rz(-1.7397225) q[1];
sx q[1];
rz(-0.090976322) q[1];
x q[2];
rz(-0.673224) q[3];
sx q[3];
rz(-2.554318) q[3];
sx q[3];
rz(-2.4622963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2174012) q[2];
sx q[2];
rz(-1.3240341) q[2];
sx q[2];
rz(0.021281555) q[2];
rz(1.4422669) q[3];
sx q[3];
rz(-2.4271624) q[3];
sx q[3];
rz(-2.2309979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
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
rz(-1.4845881) q[0];
sx q[0];
rz(-2.5452884) q[0];
sx q[0];
rz(-3.0969627) q[0];
rz(2.1394829) q[1];
sx q[1];
rz(-2.1410746) q[1];
sx q[1];
rz(3.1071641) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8196626) q[0];
sx q[0];
rz(-0.18580431) q[0];
sx q[0];
rz(1.4968605) q[0];
rz(2.0738515) q[2];
sx q[2];
rz(-1.6776553) q[2];
sx q[2];
rz(-1.6018794) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5163102) q[1];
sx q[1];
rz(-2.1045661) q[1];
sx q[1];
rz(-2.7620402) q[1];
rz(-1.7354292) q[3];
sx q[3];
rz(-2.2693686) q[3];
sx q[3];
rz(2.5005831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.78952152) q[2];
sx q[2];
rz(-1.2723203) q[2];
sx q[2];
rz(0.96674031) q[2];
rz(3.1294075) q[3];
sx q[3];
rz(-0.92958486) q[3];
sx q[3];
rz(-3.0325082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1352585) q[0];
sx q[0];
rz(-2.874458) q[0];
sx q[0];
rz(2.1475041) q[0];
rz(3.0864691) q[1];
sx q[1];
rz(-1.8693285) q[1];
sx q[1];
rz(-0.73928839) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2544884) q[0];
sx q[0];
rz(-1.6657889) q[0];
sx q[0];
rz(-1.7083005) q[0];
rz(-pi) q[1];
rz(0.67989345) q[2];
sx q[2];
rz(-0.69603622) q[2];
sx q[2];
rz(1.807715) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.97923219) q[1];
sx q[1];
rz(-1.3812961) q[1];
sx q[1];
rz(-0.53342553) q[1];
rz(1.8592632) q[3];
sx q[3];
rz(-0.97086421) q[3];
sx q[3];
rz(-2.695042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.58549515) q[2];
sx q[2];
rz(-1.408351) q[2];
sx q[2];
rz(-0.56376702) q[2];
rz(0.24027763) q[3];
sx q[3];
rz(-2.6085745) q[3];
sx q[3];
rz(-1.3747922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69650841) q[0];
sx q[0];
rz(-0.17906469) q[0];
sx q[0];
rz(0.58404303) q[0];
rz(-0.42075992) q[1];
sx q[1];
rz(-2.1003484) q[1];
sx q[1];
rz(-0.27841321) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4415092) q[0];
sx q[0];
rz(-1.004305) q[0];
sx q[0];
rz(2.9928815) q[0];
rz(2.7483447) q[2];
sx q[2];
rz(-1.3967782) q[2];
sx q[2];
rz(2.4207123) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.091374) q[1];
sx q[1];
rz(-2.723657) q[1];
sx q[1];
rz(1.0068847) q[1];
x q[2];
rz(-1.1823468) q[3];
sx q[3];
rz(-2.5591345) q[3];
sx q[3];
rz(2.9782481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2723096) q[2];
sx q[2];
rz(-0.81988207) q[2];
sx q[2];
rz(-2.753624) q[2];
rz(1.7913473) q[3];
sx q[3];
rz(-0.46348059) q[3];
sx q[3];
rz(0.2750245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85912722) q[0];
sx q[0];
rz(-0.21331856) q[0];
sx q[0];
rz(0.7318837) q[0];
rz(0.16648079) q[1];
sx q[1];
rz(-0.44547588) q[1];
sx q[1];
rz(-3.0678715) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1151704) q[0];
sx q[0];
rz(-2.0048884) q[0];
sx q[0];
rz(0.11154453) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6905626) q[2];
sx q[2];
rz(-1.4904067) q[2];
sx q[2];
rz(2.5188841) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2697061) q[1];
sx q[1];
rz(-1.5333813) q[1];
sx q[1];
rz(0.58362959) q[1];
rz(1.4599667) q[3];
sx q[3];
rz(-1.0004527) q[3];
sx q[3];
rz(-3.0487206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5471197) q[2];
sx q[2];
rz(-2.9056845) q[2];
sx q[2];
rz(-2.7129042) q[2];
rz(1.8244913) q[3];
sx q[3];
rz(-1.7842112) q[3];
sx q[3];
rz(-1.930442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73228943) q[0];
sx q[0];
rz(-1.4855054) q[0];
sx q[0];
rz(1.0428585) q[0];
rz(1.6304784) q[1];
sx q[1];
rz(-1.7621367) q[1];
sx q[1];
rz(1.3483378) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82768517) q[0];
sx q[0];
rz(-1.5965009) q[0];
sx q[0];
rz(-1.4575507) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6595608) q[2];
sx q[2];
rz(-0.73854337) q[2];
sx q[2];
rz(0.86631394) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7235782) q[1];
sx q[1];
rz(-1.2472767) q[1];
sx q[1];
rz(0.07467204) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3799558) q[3];
sx q[3];
rz(-1.8075426) q[3];
sx q[3];
rz(-2.5773347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.80031359) q[2];
sx q[2];
rz(-0.51395243) q[2];
sx q[2];
rz(-3.0467765) q[2];
rz(1.9684277) q[3];
sx q[3];
rz(-1.7404107) q[3];
sx q[3];
rz(-0.15869424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6288347) q[0];
sx q[0];
rz(-2.6678968) q[0];
sx q[0];
rz(-2.0976023) q[0];
rz(-1.5402773) q[1];
sx q[1];
rz(-1.5834783) q[1];
sx q[1];
rz(-1.6316354) q[1];
rz(-3.0586254) q[2];
sx q[2];
rz(-1.7832179) q[2];
sx q[2];
rz(-0.28945343) q[2];
rz(-2.3306866) q[3];
sx q[3];
rz(-1.1699642) q[3];
sx q[3];
rz(-1.5311833) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
