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
rz(-0.6119734) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64498211) q[0];
sx q[0];
rz(-2.5685446) q[0];
sx q[0];
rz(1.7368421) q[0];
rz(1.0704899) q[2];
sx q[2];
rz(-2.8165952) q[2];
sx q[2];
rz(-0.18735841) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8777316) q[1];
sx q[1];
rz(-2.7791443) q[1];
sx q[1];
rz(-1.1652105) q[1];
rz(-pi) q[2];
x q[2];
rz(0.060617491) q[3];
sx q[3];
rz(-2.4481886) q[3];
sx q[3];
rz(1.7934679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.97418) q[2];
sx q[2];
rz(-1.1777425) q[2];
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
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6086455) q[0];
sx q[0];
rz(-2.8605509) q[0];
sx q[0];
rz(2.1458022) q[0];
rz(-2.5536054) q[1];
sx q[1];
rz(-0.43991393) q[1];
sx q[1];
rz(-1.3194552) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73896995) q[0];
sx q[0];
rz(-1.8085294) q[0];
sx q[0];
rz(2.0031702) q[0];
rz(-pi) q[1];
x q[1];
rz(0.075602268) q[2];
sx q[2];
rz(-1.8972862) q[2];
sx q[2];
rz(1.16815) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4012137) q[1];
sx q[1];
rz(-1.7771562) q[1];
sx q[1];
rz(2.3310082) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6922582) q[3];
sx q[3];
rz(-1.0912885) q[3];
sx q[3];
rz(3.1288655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1002645) q[2];
sx q[2];
rz(-2.1375956) q[2];
sx q[2];
rz(-3.1190994) q[2];
rz(0.10447539) q[3];
sx q[3];
rz(-1.5375117) q[3];
sx q[3];
rz(-2.2711066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3290038) q[0];
sx q[0];
rz(-2.1109695) q[0];
sx q[0];
rz(-1.5033683) q[0];
rz(0.080987856) q[1];
sx q[1];
rz(-0.65892017) q[1];
sx q[1];
rz(1.9714877) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2161897) q[0];
sx q[0];
rz(-1.7316375) q[0];
sx q[0];
rz(-0.58165929) q[0];
rz(-pi) q[1];
rz(1.3721714) q[2];
sx q[2];
rz(-2.2057057) q[2];
sx q[2];
rz(-2.9609307) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7602424) q[1];
sx q[1];
rz(-2.6406857) q[1];
sx q[1];
rz(1.5111021) q[1];
x q[2];
rz(1.0634138) q[3];
sx q[3];
rz(-1.7905924) q[3];
sx q[3];
rz(0.45377094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.66204232) q[2];
sx q[2];
rz(-0.68632555) q[2];
sx q[2];
rz(0.043206841) q[2];
rz(-0.19733812) q[3];
sx q[3];
rz(-2.1102326) q[3];
sx q[3];
rz(0.432338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8869121) q[0];
sx q[0];
rz(-1.0858902) q[0];
sx q[0];
rz(-2.8160954) q[0];
rz(-0.85572851) q[1];
sx q[1];
rz(-0.41494644) q[1];
sx q[1];
rz(1.0861446) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.08642149) q[0];
sx q[0];
rz(-2.060411) q[0];
sx q[0];
rz(1.7167164) q[0];
rz(2.3065673) q[2];
sx q[2];
rz(-0.97423282) q[2];
sx q[2];
rz(1.2460097) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1126744) q[1];
sx q[1];
rz(-0.54593819) q[1];
sx q[1];
rz(1.250099) q[1];
rz(-pi) q[2];
rz(1.4277763) q[3];
sx q[3];
rz(-1.656683) q[3];
sx q[3];
rz(-0.68616223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2073652) q[2];
sx q[2];
rz(-2.1659329) q[2];
sx q[2];
rz(-0.18386851) q[2];
rz(1.81987) q[3];
sx q[3];
rz(-2.4403641) q[3];
sx q[3];
rz(0.13300657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2258218) q[0];
sx q[0];
rz(-2.2647595) q[0];
sx q[0];
rz(-3.0208605) q[0];
rz(-pi) q[1];
rz(0.59163036) q[2];
sx q[2];
rz(-1.0359633) q[2];
sx q[2];
rz(-2.4727351) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7246157) q[1];
sx q[1];
rz(-1.6682345) q[1];
sx q[1];
rz(2.3259054) q[1];
rz(-pi) q[2];
rz(1.3077626) q[3];
sx q[3];
rz(-1.9527779) q[3];
sx q[3];
rz(-0.39417496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0978284) q[2];
sx q[2];
rz(-1.8057258) q[2];
sx q[2];
rz(2.9975927) q[2];
rz(2.0740267) q[3];
sx q[3];
rz(-2.7976566) q[3];
sx q[3];
rz(2.4501154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3095793) q[0];
sx q[0];
rz(-2.3342275) q[0];
sx q[0];
rz(2.373234) q[0];
rz(0.8575303) q[1];
sx q[1];
rz(-2.0950967) q[1];
sx q[1];
rz(-2.5879477) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2087649) q[0];
sx q[0];
rz(-2.0220482) q[0];
sx q[0];
rz(0.28214595) q[0];
rz(-1.800399) q[2];
sx q[2];
rz(-0.78289778) q[2];
sx q[2];
rz(-1.7097434) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2276113) q[1];
sx q[1];
rz(-1.0385333) q[1];
sx q[1];
rz(-2.6460939) q[1];
rz(-pi) q[2];
x q[2];
rz(0.72154273) q[3];
sx q[3];
rz(-2.0399722) q[3];
sx q[3];
rz(-0.32687068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7553317) q[2];
sx q[2];
rz(-0.43060455) q[2];
sx q[2];
rz(2.6957896) q[2];
rz(-1.2465994) q[3];
sx q[3];
rz(-1.7720902) q[3];
sx q[3];
rz(2.2915452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
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
rz(0.069139473) q[0];
sx q[0];
rz(-2.1816165) q[0];
sx q[0];
rz(0.14007105) q[0];
rz(-2.2135997) q[1];
sx q[1];
rz(-1.8047787) q[1];
sx q[1];
rz(1.6915406) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1221913) q[0];
sx q[0];
rz(-2.9474576) q[0];
sx q[0];
rz(2.1737264) q[0];
x q[1];
rz(2.9102703) q[2];
sx q[2];
rz(-1.1975106) q[2];
sx q[2];
rz(-1.1333677) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2388152) q[1];
sx q[1];
rz(-2.0014781) q[1];
sx q[1];
rz(0.68044739) q[1];
rz(-pi) q[2];
x q[2];
rz(0.70287649) q[3];
sx q[3];
rz(-2.9053087) q[3];
sx q[3];
rz(0.98881665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1040087) q[2];
sx q[2];
rz(-0.17825492) q[2];
sx q[2];
rz(3.0577216) q[2];
rz(0.9907848) q[3];
sx q[3];
rz(-1.9525986) q[3];
sx q[3];
rz(1.8535463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
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
rz(1.4419979) q[1];
sx q[1];
rz(-0.60209638) q[1];
sx q[1];
rz(-3.1325565) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4215338) q[0];
sx q[0];
rz(-1.3429317) q[0];
sx q[0];
rz(2.9145545) q[0];
rz(-pi) q[1];
rz(0.70019146) q[2];
sx q[2];
rz(-1.2692361) q[2];
sx q[2];
rz(1.4789326) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9408843) q[1];
sx q[1];
rz(-1.6621331) q[1];
sx q[1];
rz(-1.8436842) q[1];
rz(-pi) q[2];
rz(0.44329109) q[3];
sx q[3];
rz(-1.6066243) q[3];
sx q[3];
rz(1.2092502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8845727) q[2];
sx q[2];
rz(-1.3468065) q[2];
sx q[2];
rz(2.436077) q[2];
rz(2.8722615) q[3];
sx q[3];
rz(-2.0651385) q[3];
sx q[3];
rz(-2.1721325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(2.2483098) q[0];
sx q[0];
rz(-2.3362384) q[0];
sx q[0];
rz(2.4321108) q[0];
rz(0.64208883) q[1];
sx q[1];
rz(-2.700192) q[1];
sx q[1];
rz(-0.4253687) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50373721) q[0];
sx q[0];
rz(-2.2048973) q[0];
sx q[0];
rz(-1.0315328) q[0];
rz(2.7284996) q[2];
sx q[2];
rz(-2.1132662) q[2];
sx q[2];
rz(-1.7187207) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9305552) q[1];
sx q[1];
rz(-1.1260707) q[1];
sx q[1];
rz(-1.6707604) q[1];
rz(1.2669417) q[3];
sx q[3];
rz(-1.0179449) q[3];
sx q[3];
rz(-1.0859717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8713259) q[2];
sx q[2];
rz(-0.17921236) q[2];
sx q[2];
rz(-2.6628009) q[2];
rz(0.35017961) q[3];
sx q[3];
rz(-1.9637354) q[3];
sx q[3];
rz(-2.71463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(2.4737074) q[0];
sx q[0];
rz(-0.29692867) q[0];
sx q[0];
rz(0.83734751) q[0];
rz(-2.8750724) q[1];
sx q[1];
rz(-1.8217249) q[1];
sx q[1];
rz(-1.7841608) q[1];
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
rz(-1.6893759) q[2];
sx q[2];
rz(-2.0381387) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7818004) q[1];
sx q[1];
rz(-2.4891653) q[1];
sx q[1];
rz(3.1210207) q[1];
x q[2];
rz(2.0344026) q[3];
sx q[3];
rz(-1.2779116) q[3];
sx q[3];
rz(1.1395265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.16167656) q[2];
sx q[2];
rz(-2.6075173) q[2];
sx q[2];
rz(-2.7034289) q[2];
rz(2.2257889) q[3];
sx q[3];
rz(-1.5235498) q[3];
sx q[3];
rz(2.3384136) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5220779) q[0];
sx q[0];
rz(-1.8997471) q[0];
sx q[0];
rz(2.1258623) q[0];
rz(3.0370514) q[1];
sx q[1];
rz(-1.3019982) q[1];
sx q[1];
rz(-1.7560538) q[1];
rz(-1.3913515) q[2];
sx q[2];
rz(-0.33534563) q[2];
sx q[2];
rz(-1.705359) q[2];
rz(0.87540353) q[3];
sx q[3];
rz(-2.2973552) q[3];
sx q[3];
rz(-2.4612343) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
