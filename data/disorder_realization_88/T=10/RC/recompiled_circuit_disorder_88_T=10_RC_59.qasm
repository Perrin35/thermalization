OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8945583) q[0];
sx q[0];
rz(-1.2556827) q[0];
sx q[0];
rz(-0.32796252) q[0];
rz(-0.2344996) q[1];
sx q[1];
rz(3.3426715) q[1];
sx q[1];
rz(9.3333416) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8592005) q[0];
sx q[0];
rz(-2.1452367) q[0];
sx q[0];
rz(-2.0496856) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7056866) q[2];
sx q[2];
rz(-1.8276668) q[2];
sx q[2];
rz(2.169662) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.92802231) q[1];
sx q[1];
rz(-1.3334647) q[1];
sx q[1];
rz(-0.084753239) q[1];
rz(-pi) q[2];
x q[2];
rz(2.923008) q[3];
sx q[3];
rz(-0.9871452) q[3];
sx q[3];
rz(-2.1368795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.66449195) q[2];
sx q[2];
rz(-0.97350073) q[2];
sx q[2];
rz(2.0155902) q[2];
rz(0.27515718) q[3];
sx q[3];
rz(-2.5313009) q[3];
sx q[3];
rz(2.2385105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1317516) q[0];
sx q[0];
rz(-0.69350243) q[0];
sx q[0];
rz(2.4480208) q[0];
rz(1.0961078) q[1];
sx q[1];
rz(-2.1577436) q[1];
sx q[1];
rz(2.9512761) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41102558) q[0];
sx q[0];
rz(-1.1682296) q[0];
sx q[0];
rz(2.9257724) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61043592) q[2];
sx q[2];
rz(-2.5154841) q[2];
sx q[2];
rz(-0.49952835) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7479334) q[1];
sx q[1];
rz(-0.84888443) q[1];
sx q[1];
rz(0.02971239) q[1];
rz(-pi) q[2];
rz(-0.5204366) q[3];
sx q[3];
rz(-2.3428829) q[3];
sx q[3];
rz(-2.6526116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.50743121) q[2];
sx q[2];
rz(-2.4298411) q[2];
sx q[2];
rz(-1.2197536) q[2];
rz(2.9988585) q[3];
sx q[3];
rz(-2.1295363) q[3];
sx q[3];
rz(0.90863168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.398657) q[0];
sx q[0];
rz(-1.0899028) q[0];
sx q[0];
rz(0.8202585) q[0];
rz(-2.8495158) q[1];
sx q[1];
rz(-1.074011) q[1];
sx q[1];
rz(-1.8935727) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6662212) q[0];
sx q[0];
rz(-2.1682122) q[0];
sx q[0];
rz(-1.4587547) q[0];
x q[1];
rz(1.6349995) q[2];
sx q[2];
rz(-1.3639796) q[2];
sx q[2];
rz(-1.1484255) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.51860147) q[1];
sx q[1];
rz(-0.90497436) q[1];
sx q[1];
rz(-1.4856505) q[1];
x q[2];
rz(-1.9565342) q[3];
sx q[3];
rz(-0.68813656) q[3];
sx q[3];
rz(0.08337534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.32039207) q[2];
sx q[2];
rz(-2.0911066) q[2];
sx q[2];
rz(-1.2134264) q[2];
rz(0.16472566) q[3];
sx q[3];
rz(-2.2475593) q[3];
sx q[3];
rz(-3.059982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40959013) q[0];
sx q[0];
rz(-1.859917) q[0];
sx q[0];
rz(0.18606342) q[0];
rz(-2.9371254) q[1];
sx q[1];
rz(-0.47195131) q[1];
sx q[1];
rz(-1.2971372) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0766749) q[0];
sx q[0];
rz(-3.0041822) q[0];
sx q[0];
rz(1.2491559) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9739431) q[2];
sx q[2];
rz(-1.8561346) q[2];
sx q[2];
rz(-2.167785) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.101673) q[1];
sx q[1];
rz(-0.97517255) q[1];
sx q[1];
rz(0.69570978) q[1];
x q[2];
rz(-1.3912195) q[3];
sx q[3];
rz(-1.5571801) q[3];
sx q[3];
rz(0.40819685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2356448) q[2];
sx q[2];
rz(-2.3489958) q[2];
sx q[2];
rz(-2.2223991) q[2];
rz(-0.32133189) q[3];
sx q[3];
rz(-2.0733757) q[3];
sx q[3];
rz(1.8937768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17386757) q[0];
sx q[0];
rz(-0.50650948) q[0];
sx q[0];
rz(0.87819535) q[0];
rz(1.325266) q[1];
sx q[1];
rz(-1.7293431) q[1];
sx q[1];
rz(1.7153046) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8910687) q[0];
sx q[0];
rz(-0.91618012) q[0];
sx q[0];
rz(3.0380681) q[0];
rz(-0.23767383) q[2];
sx q[2];
rz(-1.2665247) q[2];
sx q[2];
rz(-1.301847) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.95477415) q[1];
sx q[1];
rz(-0.84328077) q[1];
sx q[1];
rz(0.14528841) q[1];
x q[2];
rz(0.32894965) q[3];
sx q[3];
rz(-2.3429686) q[3];
sx q[3];
rz(0.22508276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2720126) q[2];
sx q[2];
rz(-2*pi/13) q[2];
sx q[2];
rz(-0.041794725) q[2];
rz(-0.061491866) q[3];
sx q[3];
rz(-1.0890591) q[3];
sx q[3];
rz(-2.7048236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7866661) q[0];
sx q[0];
rz(-3.0556398) q[0];
sx q[0];
rz(1.0304931) q[0];
rz(0.73973918) q[1];
sx q[1];
rz(-1.5286427) q[1];
sx q[1];
rz(-0.57156634) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3540928) q[0];
sx q[0];
rz(-1.0618853) q[0];
sx q[0];
rz(-0.66977588) q[0];
rz(1.4095441) q[2];
sx q[2];
rz(-2.469392) q[2];
sx q[2];
rz(1.6598998) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7476269) q[1];
sx q[1];
rz(-0.52226258) q[1];
sx q[1];
rz(2.0678492) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3011175) q[3];
sx q[3];
rz(-2.2395036) q[3];
sx q[3];
rz(0.83276444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.968154) q[2];
sx q[2];
rz(-2.3996694) q[2];
sx q[2];
rz(-2.5773933) q[2];
rz(-3.0155904) q[3];
sx q[3];
rz(-1.6832451) q[3];
sx q[3];
rz(2.8469767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0903704) q[0];
sx q[0];
rz(-0.19989656) q[0];
sx q[0];
rz(-0.71227658) q[0];
rz(-2.6157216) q[1];
sx q[1];
rz(-2.7253175) q[1];
sx q[1];
rz(2.4760822) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3367046) q[0];
sx q[0];
rz(-0.24222736) q[0];
sx q[0];
rz(-1.5664943) q[0];
rz(-1.7475142) q[2];
sx q[2];
rz(-0.70509796) q[2];
sx q[2];
rz(-1.9213898) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8477518) q[1];
sx q[1];
rz(-1.7898702) q[1];
sx q[1];
rz(2.0391383) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1036759) q[3];
sx q[3];
rz(-0.27927342) q[3];
sx q[3];
rz(2.6168407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9843288) q[2];
sx q[2];
rz(-0.66005808) q[2];
sx q[2];
rz(-1.4228014) q[2];
rz(3.026399) q[3];
sx q[3];
rz(-0.51518232) q[3];
sx q[3];
rz(-0.19259024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.049906235) q[0];
sx q[0];
rz(-2.005907) q[0];
sx q[0];
rz(-2.9507622) q[0];
rz(2.514839) q[1];
sx q[1];
rz(-2.125506) q[1];
sx q[1];
rz(-0.33871067) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3074293) q[0];
sx q[0];
rz(-1.2273664) q[0];
sx q[0];
rz(-0.15983454) q[0];
x q[1];
rz(0.45192265) q[2];
sx q[2];
rz(-1.3616614) q[2];
sx q[2];
rz(-2.506633) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.72570669) q[1];
sx q[1];
rz(-0.48735122) q[1];
sx q[1];
rz(-1.332167) q[1];
rz(0.2770822) q[3];
sx q[3];
rz(-1.2556561) q[3];
sx q[3];
rz(1.806123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.64615858) q[2];
sx q[2];
rz(-0.63445264) q[2];
sx q[2];
rz(-0.70043606) q[2];
rz(-2.2436079) q[3];
sx q[3];
rz(-1.8549517) q[3];
sx q[3];
rz(-2.9680796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53196466) q[0];
sx q[0];
rz(-2.6599929) q[0];
sx q[0];
rz(-2.4832446) q[0];
rz(2.530653) q[1];
sx q[1];
rz(-1.8860203) q[1];
sx q[1];
rz(-0.13959612) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7394373) q[0];
sx q[0];
rz(-1.6318775) q[0];
sx q[0];
rz(0.020045965) q[0];
rz(-pi) q[1];
rz(1.5084969) q[2];
sx q[2];
rz(-0.75714105) q[2];
sx q[2];
rz(1.8119259) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.67006754) q[1];
sx q[1];
rz(-1.0499665) q[1];
sx q[1];
rz(-1.4937558) q[1];
x q[2];
rz(-0.052501909) q[3];
sx q[3];
rz(-2.2146261) q[3];
sx q[3];
rz(-0.98464636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6167986) q[2];
sx q[2];
rz(-1.0107661) q[2];
sx q[2];
rz(-2.810478) q[2];
rz(-0.75774276) q[3];
sx q[3];
rz(-2.7527633) q[3];
sx q[3];
rz(3.0537135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8452633) q[0];
sx q[0];
rz(-2.1305278) q[0];
sx q[0];
rz(2.9560126) q[0];
rz(-1.0962076) q[1];
sx q[1];
rz(-0.21462333) q[1];
sx q[1];
rz(1.6569482) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97666868) q[0];
sx q[0];
rz(-0.52545588) q[0];
sx q[0];
rz(-2.1303961) q[0];
x q[1];
rz(1.2946285) q[2];
sx q[2];
rz(-1.7880926) q[2];
sx q[2];
rz(1.5751788) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9975002) q[1];
sx q[1];
rz(-1.1966685) q[1];
sx q[1];
rz(3.0348026) q[1];
rz(-pi) q[2];
x q[2];
rz(0.86587571) q[3];
sx q[3];
rz(-1.5670334) q[3];
sx q[3];
rz(-0.83434425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.93402702) q[2];
sx q[2];
rz(-2.2521314) q[2];
sx q[2];
rz(2.5893842) q[2];
rz(-2.3637555) q[3];
sx q[3];
rz(-0.85770291) q[3];
sx q[3];
rz(2.623693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8778397) q[0];
sx q[0];
rz(-1.5120266) q[0];
sx q[0];
rz(1.7396447) q[0];
rz(0.18763018) q[1];
sx q[1];
rz(-1.7788806) q[1];
sx q[1];
rz(2.3685041) q[1];
rz(-2.8616703) q[2];
sx q[2];
rz(-1.4301849) q[2];
sx q[2];
rz(2.4591597) q[2];
rz(-2.3579303) q[3];
sx q[3];
rz(-2.2829934) q[3];
sx q[3];
rz(-0.34096277) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
