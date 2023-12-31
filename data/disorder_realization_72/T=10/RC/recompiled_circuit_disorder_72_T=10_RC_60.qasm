OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.62087286) q[0];
sx q[0];
rz(1.7680661) q[0];
sx q[0];
rz(11.058523) q[0];
rz(0.047343407) q[1];
sx q[1];
rz(-2.3634057) q[1];
sx q[1];
rz(0.49931061) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7802785) q[0];
sx q[0];
rz(-1.8457992) q[0];
sx q[0];
rz(-1.6367903) q[0];
x q[1];
rz(-1.7668031) q[2];
sx q[2];
rz(-0.73050806) q[2];
sx q[2];
rz(1.0175878) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3589188) q[1];
sx q[1];
rz(-2.3896857) q[1];
sx q[1];
rz(2.9208675) q[1];
x q[2];
rz(-1.3189089) q[3];
sx q[3];
rz(-0.10676521) q[3];
sx q[3];
rz(1.3085384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.22380655) q[2];
sx q[2];
rz(-2.1710158) q[2];
sx q[2];
rz(-1.0144368) q[2];
rz(-0.23400083) q[3];
sx q[3];
rz(-2.6205385) q[3];
sx q[3];
rz(-2.8570989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4734128) q[0];
sx q[0];
rz(-1.6750591) q[0];
sx q[0];
rz(2.1372674) q[0];
rz(-1.6218119) q[1];
sx q[1];
rz(-0.92679778) q[1];
sx q[1];
rz(1.0027592) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6111345) q[0];
sx q[0];
rz(-2.3728752) q[0];
sx q[0];
rz(2.4918633) q[0];
rz(-2.3147896) q[2];
sx q[2];
rz(-2.8892592) q[2];
sx q[2];
rz(1.2834872) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4304632) q[1];
sx q[1];
rz(-1.7257479) q[1];
sx q[1];
rz(0.98053812) q[1];
rz(0.49591222) q[3];
sx q[3];
rz(-0.65973385) q[3];
sx q[3];
rz(-2.0663313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.26943794) q[2];
sx q[2];
rz(-0.99439159) q[2];
sx q[2];
rz(-0.15094748) q[2];
rz(0.41444591) q[3];
sx q[3];
rz(-0.60025418) q[3];
sx q[3];
rz(-3.0533561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8931483) q[0];
sx q[0];
rz(-1.3131161) q[0];
sx q[0];
rz(0.77899581) q[0];
rz(-2.7422854) q[1];
sx q[1];
rz(-1.2483968) q[1];
sx q[1];
rz(2.2580106) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8944091) q[0];
sx q[0];
rz(-2.7054225) q[0];
sx q[0];
rz(0.39788525) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9687037) q[2];
sx q[2];
rz(-1.298561) q[2];
sx q[2];
rz(0.18323252) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5501432) q[1];
sx q[1];
rz(-1.5482229) q[1];
sx q[1];
rz(-1.2858461) q[1];
rz(-pi) q[2];
x q[2];
rz(0.15150841) q[3];
sx q[3];
rz(-2.476859) q[3];
sx q[3];
rz(-2.1745149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.3016004) q[2];
sx q[2];
rz(-1.2568544) q[2];
sx q[2];
rz(0.42923129) q[2];
rz(2.1515576) q[3];
sx q[3];
rz(-2.0638549) q[3];
sx q[3];
rz(-0.86301962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(0.4317959) q[0];
sx q[0];
rz(-1.3732095) q[0];
sx q[0];
rz(0.61169949) q[0];
rz(-2.0344095) q[1];
sx q[1];
rz(-2.2829843) q[1];
sx q[1];
rz(-2.591419) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.600425) q[0];
sx q[0];
rz(-2.4006002) q[0];
sx q[0];
rz(-3.0279972) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8407211) q[2];
sx q[2];
rz(-0.30756018) q[2];
sx q[2];
rz(0.077066271) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.95851129) q[1];
sx q[1];
rz(-1.4554879) q[1];
sx q[1];
rz(0.99460852) q[1];
rz(1.6226193) q[3];
sx q[3];
rz(-1.341815) q[3];
sx q[3];
rz(0.94829544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.52175534) q[2];
sx q[2];
rz(-2.6553314) q[2];
sx q[2];
rz(-2.8660529) q[2];
rz(0.11166212) q[3];
sx q[3];
rz(-1.9409981) q[3];
sx q[3];
rz(-0.47079852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4438641) q[0];
sx q[0];
rz(-1.0794909) q[0];
sx q[0];
rz(0.87669796) q[0];
rz(0.69119167) q[1];
sx q[1];
rz(-2.2677939) q[1];
sx q[1];
rz(0.91526389) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4782151) q[0];
sx q[0];
rz(-1.6179248) q[0];
sx q[0];
rz(2.1992654) q[0];
rz(-pi) q[1];
rz(3.0643164) q[2];
sx q[2];
rz(-2.6151267) q[2];
sx q[2];
rz(0.72409814) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.535714) q[1];
sx q[1];
rz(-1.5389581) q[1];
sx q[1];
rz(-1.72623) q[1];
rz(2.3239273) q[3];
sx q[3];
rz(-2.2125803) q[3];
sx q[3];
rz(-1.6587917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9203732) q[2];
sx q[2];
rz(-1.0125786) q[2];
sx q[2];
rz(2.6780224) q[2];
rz(0.56435895) q[3];
sx q[3];
rz(-2.1488991) q[3];
sx q[3];
rz(-0.92818964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1148465) q[0];
sx q[0];
rz(-2.5214054) q[0];
sx q[0];
rz(-2.0625431) q[0];
rz(-2.5462529) q[1];
sx q[1];
rz(-2.3880312) q[1];
sx q[1];
rz(0.39658305) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.840344) q[0];
sx q[0];
rz(-0.52919555) q[0];
sx q[0];
rz(-1.0521786) q[0];
x q[1];
rz(2.6655212) q[2];
sx q[2];
rz(-1.8487612) q[2];
sx q[2];
rz(-1.4897886) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.10627667) q[1];
sx q[1];
rz(-0.66674399) q[1];
sx q[1];
rz(1.7229401) q[1];
x q[2];
rz(-1.3271689) q[3];
sx q[3];
rz(-1.0001839) q[3];
sx q[3];
rz(2.7057196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.98809272) q[2];
sx q[2];
rz(-2.2160539) q[2];
sx q[2];
rz(1.5552103) q[2];
rz(-1.68613) q[3];
sx q[3];
rz(-0.60417914) q[3];
sx q[3];
rz(0.82715183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39847386) q[0];
sx q[0];
rz(-1.159659) q[0];
sx q[0];
rz(-2.356785) q[0];
rz(1.8709042) q[1];
sx q[1];
rz(-1.7653468) q[1];
sx q[1];
rz(-1.126359) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2369909) q[0];
sx q[0];
rz(-2.8303574) q[0];
sx q[0];
rz(2.8471332) q[0];
rz(-pi) q[1];
rz(-2.0790714) q[2];
sx q[2];
rz(-1.695343) q[2];
sx q[2];
rz(1.6342271) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2396444) q[1];
sx q[1];
rz(-0.48929729) q[1];
sx q[1];
rz(-0.40892618) q[1];
rz(-pi) q[2];
rz(-1.8282312) q[3];
sx q[3];
rz(-0.35850393) q[3];
sx q[3];
rz(3.0646216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.209098) q[2];
sx q[2];
rz(-1.1065437) q[2];
sx q[2];
rz(-0.15360019) q[2];
rz(-0.30512729) q[3];
sx q[3];
rz(-2.025445) q[3];
sx q[3];
rz(-1.37384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7711733) q[0];
sx q[0];
rz(-2.6035247) q[0];
sx q[0];
rz(0.11288189) q[0];
rz(-1.0007292) q[1];
sx q[1];
rz(-0.74644867) q[1];
sx q[1];
rz(-3.1088366) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2643124) q[0];
sx q[0];
rz(-2.7240629) q[0];
sx q[0];
rz(2.2420922) q[0];
rz(2.2135127) q[2];
sx q[2];
rz(-0.66814458) q[2];
sx q[2];
rz(-1.1754787) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7777562) q[1];
sx q[1];
rz(-1.5118074) q[1];
sx q[1];
rz(1.9883518) q[1];
x q[2];
rz(1.1912187) q[3];
sx q[3];
rz(-1.4879585) q[3];
sx q[3];
rz(2.6759202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.96413606) q[2];
sx q[2];
rz(-0.63825858) q[2];
sx q[2];
rz(-0.62210554) q[2];
rz(-1.9744251) q[3];
sx q[3];
rz(-1.111235) q[3];
sx q[3];
rz(0.39045236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.099982925) q[0];
sx q[0];
rz(-2.6384625) q[0];
sx q[0];
rz(-1.5266248) q[0];
rz(-0.73293066) q[1];
sx q[1];
rz(-2.4885978) q[1];
sx q[1];
rz(-0.48535767) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63472414) q[0];
sx q[0];
rz(-1.6244495) q[0];
sx q[0];
rz(3.0070478) q[0];
rz(-pi) q[1];
rz(0.59818563) q[2];
sx q[2];
rz(-1.4387812) q[2];
sx q[2];
rz(2.8796632) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.24385142) q[1];
sx q[1];
rz(-0.85695367) q[1];
sx q[1];
rz(2.7276911) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.63391179) q[3];
sx q[3];
rz(-0.68448193) q[3];
sx q[3];
rz(-0.67542911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0316524) q[2];
sx q[2];
rz(-1.8376708) q[2];
sx q[2];
rz(-0.1594485) q[2];
rz(-1.4032646) q[3];
sx q[3];
rz(-2.7519029) q[3];
sx q[3];
rz(-3.1260417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6195246) q[0];
sx q[0];
rz(-0.63729006) q[0];
sx q[0];
rz(-1.7074701) q[0];
rz(-1.8822949) q[1];
sx q[1];
rz(-2.1914296) q[1];
sx q[1];
rz(2.5433345) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78865096) q[0];
sx q[0];
rz(-1.3539679) q[0];
sx q[0];
rz(-1.7822595) q[0];
x q[1];
rz(-1.422545) q[2];
sx q[2];
rz(-1.3607927) q[2];
sx q[2];
rz(-1.6809747) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.88565247) q[1];
sx q[1];
rz(-1.9347895) q[1];
sx q[1];
rz(0.10140681) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.80501276) q[3];
sx q[3];
rz(-2.9152438) q[3];
sx q[3];
rz(-1.5387907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0593807) q[2];
sx q[2];
rz(-1.2450612) q[2];
sx q[2];
rz(0.65199488) q[2];
rz(-0.59371289) q[3];
sx q[3];
rz(-1.9605325) q[3];
sx q[3];
rz(-2.0098856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-1.839529) q[0];
sx q[0];
rz(-0.65757127) q[0];
sx q[0];
rz(0.99304637) q[0];
rz(1.2364173) q[1];
sx q[1];
rz(-2.1724783) q[1];
sx q[1];
rz(1.9289) q[1];
rz(-0.17157208) q[2];
sx q[2];
rz(-0.62660672) q[2];
sx q[2];
rz(-1.3852711) q[2];
rz(-2.7491309) q[3];
sx q[3];
rz(-2.2069208) q[3];
sx q[3];
rz(1.8451286) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
