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
rz(-0.62701464) q[0];
sx q[0];
rz(3.7778683) q[0];
sx q[0];
rz(10.502622) q[0];
rz(1.0765422) q[1];
sx q[1];
rz(-0.62907469) q[1];
sx q[1];
rz(2.10973) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92562308) q[0];
sx q[0];
rz(-1.5703989) q[0];
sx q[0];
rz(-0.0041603869) q[0];
x q[1];
rz(0.59649845) q[2];
sx q[2];
rz(-2.7257127) q[2];
sx q[2];
rz(0.88230995) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2638322) q[1];
sx q[1];
rz(-1.9139557) q[1];
sx q[1];
rz(-2.7105217) q[1];
x q[2];
rz(2.8058654) q[3];
sx q[3];
rz(-1.6435267) q[3];
sx q[3];
rz(2.5617395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2304113) q[2];
sx q[2];
rz(-1.1831256) q[2];
sx q[2];
rz(2.6739056) q[2];
rz(-0.45105252) q[3];
sx q[3];
rz(-2.7428198) q[3];
sx q[3];
rz(2.8635645) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69538799) q[0];
sx q[0];
rz(-2.9187293) q[0];
sx q[0];
rz(-0.15431246) q[0];
rz(-0.34126869) q[1];
sx q[1];
rz(-2.7879265) q[1];
sx q[1];
rz(2.7214859) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.882928) q[0];
sx q[0];
rz(-2.1999584) q[0];
sx q[0];
rz(1.8101519) q[0];
x q[1];
rz(-3.0947761) q[2];
sx q[2];
rz(-2.280378) q[2];
sx q[2];
rz(0.094120838) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1359033) q[1];
sx q[1];
rz(-1.5694251) q[1];
sx q[1];
rz(0.47079177) q[1];
x q[2];
rz(0.16815925) q[3];
sx q[3];
rz(-2.3871867) q[3];
sx q[3];
rz(1.2446294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.61098617) q[2];
sx q[2];
rz(-2.2425118) q[2];
sx q[2];
rz(-0.77303028) q[2];
rz(0.84345877) q[3];
sx q[3];
rz(-1.5777028) q[3];
sx q[3];
rz(2.5696866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17324363) q[0];
sx q[0];
rz(-1.0272212) q[0];
sx q[0];
rz(-2.9395043) q[0];
rz(2.659722) q[1];
sx q[1];
rz(-0.20129573) q[1];
sx q[1];
rz(0.906382) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2807686) q[0];
sx q[0];
rz(-1.1204136) q[0];
sx q[0];
rz(-2.3748939) q[0];
rz(-pi) q[1];
rz(-0.34445539) q[2];
sx q[2];
rz(-2.0018501) q[2];
sx q[2];
rz(1.5004917) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3292092) q[1];
sx q[1];
rz(-2.1076729) q[1];
sx q[1];
rz(-1.3841793) q[1];
rz(2.0790948) q[3];
sx q[3];
rz(-0.33514272) q[3];
sx q[3];
rz(-0.80860521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0704982) q[2];
sx q[2];
rz(-2.0281894) q[2];
sx q[2];
rz(2.5615198) q[2];
rz(-1.5602559) q[3];
sx q[3];
rz(-1.7507078) q[3];
sx q[3];
rz(1.7177379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.4105014) q[0];
sx q[0];
rz(-3.0803362) q[0];
sx q[0];
rz(-3.0938003) q[0];
rz(-1.8675249) q[1];
sx q[1];
rz(-1.2893226) q[1];
sx q[1];
rz(-0.045225708) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2936989) q[0];
sx q[0];
rz(-1.4220823) q[0];
sx q[0];
rz(-0.92647628) q[0];
x q[1];
rz(-2.8130262) q[2];
sx q[2];
rz(-2.1801278) q[2];
sx q[2];
rz(-0.48664618) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3056335) q[1];
sx q[1];
rz(-0.0063414185) q[1];
sx q[1];
rz(2.0672227) q[1];
x q[2];
rz(2.1158695) q[3];
sx q[3];
rz(-1.1909606) q[3];
sx q[3];
rz(-1.7120509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.73146003) q[2];
sx q[2];
rz(-1.2612017) q[2];
sx q[2];
rz(0.88978466) q[2];
rz(-0.54640031) q[3];
sx q[3];
rz(-2.140464) q[3];
sx q[3];
rz(-0.31118292) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50362098) q[0];
sx q[0];
rz(-1.6872971) q[0];
sx q[0];
rz(1.7244435) q[0];
rz(1.1196989) q[1];
sx q[1];
rz(-1.3816625) q[1];
sx q[1];
rz(-1.1823357) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6851824) q[0];
sx q[0];
rz(-1.6401924) q[0];
sx q[0];
rz(-0.38763898) q[0];
rz(-pi) q[1];
rz(2.3619217) q[2];
sx q[2];
rz(-2.5354249) q[2];
sx q[2];
rz(-3.0749644) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.97302283) q[1];
sx q[1];
rz(-0.72271361) q[1];
sx q[1];
rz(0.074345868) q[1];
rz(-pi) q[2];
rz(2.3919258) q[3];
sx q[3];
rz(-1.2359859) q[3];
sx q[3];
rz(2.1129169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.917439) q[2];
sx q[2];
rz(-0.76442337) q[2];
sx q[2];
rz(2.7177641) q[2];
rz(-1.1118927) q[3];
sx q[3];
rz(-2.6424776) q[3];
sx q[3];
rz(0.65139884) q[3];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3424585) q[0];
sx q[0];
rz(-0.0722216) q[0];
sx q[0];
rz(-2.1891731) q[0];
rz(0.77596387) q[1];
sx q[1];
rz(-0.9351848) q[1];
sx q[1];
rz(1.5663358) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.090388894) q[0];
sx q[0];
rz(-3.0984833) q[0];
sx q[0];
rz(0.48844047) q[0];
rz(0.03869073) q[2];
sx q[2];
rz(-0.96444791) q[2];
sx q[2];
rz(2.4232466) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1426865) q[1];
sx q[1];
rz(-0.64780462) q[1];
sx q[1];
rz(0.36869375) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9569472) q[3];
sx q[3];
rz(-1.0870516) q[3];
sx q[3];
rz(0.32102206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7427407) q[2];
sx q[2];
rz(-0.80790085) q[2];
sx q[2];
rz(-2.3952132) q[2];
rz(1.0910723) q[3];
sx q[3];
rz(-1.7079791) q[3];
sx q[3];
rz(2.5892042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2090476) q[0];
sx q[0];
rz(-0.046367558) q[0];
sx q[0];
rz(-2.5982502) q[0];
rz(-0.53687334) q[1];
sx q[1];
rz(-0.32506341) q[1];
sx q[1];
rz(-0.12319014) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9730102) q[0];
sx q[0];
rz(-1.6014528) q[0];
sx q[0];
rz(0.93498303) q[0];
rz(1.1511717) q[2];
sx q[2];
rz(-1.0977572) q[2];
sx q[2];
rz(-1.1993739) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1363594) q[1];
sx q[1];
rz(-2.2670304) q[1];
sx q[1];
rz(-0.74058786) q[1];
x q[2];
rz(-2.3402392) q[3];
sx q[3];
rz(-2.1983302) q[3];
sx q[3];
rz(0.31454861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8366375) q[2];
sx q[2];
rz(-1.7540437) q[2];
sx q[2];
rz(2.6564964) q[2];
rz(-1.1525611) q[3];
sx q[3];
rz(-1.3644812) q[3];
sx q[3];
rz(-3.0936354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7570067) q[0];
sx q[0];
rz(-2.204019) q[0];
sx q[0];
rz(-0.49569976) q[0];
rz(-0.063830201) q[1];
sx q[1];
rz(-0.7494691) q[1];
sx q[1];
rz(-3.0705423) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9793689) q[0];
sx q[0];
rz(-2.6155229) q[0];
sx q[0];
rz(2.6499676) q[0];
rz(0.60511968) q[2];
sx q[2];
rz(-1.7150884) q[2];
sx q[2];
rz(-0.52478204) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7531705) q[1];
sx q[1];
rz(-0.28825295) q[1];
sx q[1];
rz(0.99797319) q[1];
x q[2];
rz(-1.4859753) q[3];
sx q[3];
rz(-2.0165947) q[3];
sx q[3];
rz(-2.5100207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.12585982) q[2];
sx q[2];
rz(-1.4733529) q[2];
sx q[2];
rz(0.86137548) q[2];
rz(-1.3043978) q[3];
sx q[3];
rz(-0.574489) q[3];
sx q[3];
rz(1.5931574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9651589) q[0];
sx q[0];
rz(-2.6506944) q[0];
sx q[0];
rz(1.7023671) q[0];
rz(3.0610415) q[1];
sx q[1];
rz(-1.6702024) q[1];
sx q[1];
rz(1.390994) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4347067) q[0];
sx q[0];
rz(-1.1483743) q[0];
sx q[0];
rz(2.1377853) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8755781) q[2];
sx q[2];
rz(-1.551318) q[2];
sx q[2];
rz(-1.8754885) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.94165588) q[1];
sx q[1];
rz(-0.98683954) q[1];
sx q[1];
rz(-2.5879509) q[1];
rz(-pi) q[2];
rz(0.79896121) q[3];
sx q[3];
rz(-2.6332601) q[3];
sx q[3];
rz(1.5647183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.78802687) q[2];
sx q[2];
rz(-1.8466419) q[2];
sx q[2];
rz(-3.0143152) q[2];
rz(2.2150529) q[3];
sx q[3];
rz(-2.8997731) q[3];
sx q[3];
rz(1.4827137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38630286) q[0];
sx q[0];
rz(-2.4918064) q[0];
sx q[0];
rz(1.6408828) q[0];
rz(0.81037784) q[1];
sx q[1];
rz(-2.2893298) q[1];
sx q[1];
rz(0.77218974) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0579073) q[0];
sx q[0];
rz(-0.9389239) q[0];
sx q[0];
rz(-2.4993012) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8462662) q[2];
sx q[2];
rz(-1.3662896) q[2];
sx q[2];
rz(2.6596065) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4160847) q[1];
sx q[1];
rz(-0.5495175) q[1];
sx q[1];
rz(1.5931507) q[1];
x q[2];
rz(0.051372929) q[3];
sx q[3];
rz(-0.89150611) q[3];
sx q[3];
rz(1.5254586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7605674) q[2];
sx q[2];
rz(-0.80485359) q[2];
sx q[2];
rz(-2.5123361) q[2];
rz(-2.4106846) q[3];
sx q[3];
rz(-0.84354246) q[3];
sx q[3];
rz(1.4430911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-2.6998941) q[0];
sx q[0];
rz(-0.48342539) q[0];
sx q[0];
rz(-0.40406686) q[0];
rz(2.7370257) q[1];
sx q[1];
rz(-1.7351983) q[1];
sx q[1];
rz(2.4529967) q[1];
rz(-1.5028904) q[2];
sx q[2];
rz(-1.2896027) q[2];
sx q[2];
rz(2.2647535) q[2];
rz(-0.47687198) q[3];
sx q[3];
rz(-0.58860368) q[3];
sx q[3];
rz(0.5640201) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
