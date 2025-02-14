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
rz(-2.4515406) q[0];
sx q[0];
rz(-2.2909988) q[0];
sx q[0];
rz(2.9414862) q[0];
rz(2.9432358) q[1];
sx q[1];
rz(-2.618572) q[1];
sx q[1];
rz(-1.3318292) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49591767) q[0];
sx q[0];
rz(-1.1326651) q[0];
sx q[0];
rz(-2.8854831) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4780689) q[2];
sx q[2];
rz(-2.4276456) q[2];
sx q[2];
rz(0.51807846) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6864904) q[1];
sx q[1];
rz(-2.2481866) q[1];
sx q[1];
rz(1.5090921) q[1];
rz(-pi) q[2];
rz(-2.0468007) q[3];
sx q[3];
rz(-1.4972357) q[3];
sx q[3];
rz(-1.5293763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0693822) q[2];
sx q[2];
rz(-0.36082265) q[2];
sx q[2];
rz(-1.8270095) q[2];
rz(2.2383111) q[3];
sx q[3];
rz(-1.3638146) q[3];
sx q[3];
rz(2.7956853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8697934) q[0];
sx q[0];
rz(-1.485774) q[0];
sx q[0];
rz(-0.85579175) q[0];
rz(-0.77404147) q[1];
sx q[1];
rz(-1.1639405) q[1];
sx q[1];
rz(0.78786293) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1163569) q[0];
sx q[0];
rz(-1.5486002) q[0];
sx q[0];
rz(1.8216013) q[0];
rz(-pi) q[1];
rz(1.5283952) q[2];
sx q[2];
rz(-1.9147434) q[2];
sx q[2];
rz(-2.6478772) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7447506) q[1];
sx q[1];
rz(-1.1424658) q[1];
sx q[1];
rz(2.7672) q[1];
x q[2];
rz(1.3608906) q[3];
sx q[3];
rz(-0.91448254) q[3];
sx q[3];
rz(-2.8950952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0756125) q[2];
sx q[2];
rz(-1.9115261) q[2];
sx q[2];
rz(0.39609972) q[2];
rz(-1.5710477) q[3];
sx q[3];
rz(-1.5038265) q[3];
sx q[3];
rz(-1.3554696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1876672) q[0];
sx q[0];
rz(-1.4952156) q[0];
sx q[0];
rz(3.0809825) q[0];
rz(-2.4568779) q[1];
sx q[1];
rz(-1.4002607) q[1];
sx q[1];
rz(2.8025467) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1364839) q[0];
sx q[0];
rz(-1.5915131) q[0];
sx q[0];
rz(-2.5584975) q[0];
x q[1];
rz(-1.6940368) q[2];
sx q[2];
rz(-1.7584821) q[2];
sx q[2];
rz(-1.0836243) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3587957) q[1];
sx q[1];
rz(-2.7227245) q[1];
sx q[1];
rz(2.2936506) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5825082) q[3];
sx q[3];
rz(-2.3475966) q[3];
sx q[3];
rz(-2.7744966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6591349) q[2];
sx q[2];
rz(-1.277801) q[2];
sx q[2];
rz(0.47046146) q[2];
rz(2.2792234) q[3];
sx q[3];
rz(-0.44895288) q[3];
sx q[3];
rz(1.5615777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3622793) q[0];
sx q[0];
rz(-2.0231495) q[0];
sx q[0];
rz(-2.733574) q[0];
rz(-1.7904003) q[1];
sx q[1];
rz(-1.2974757) q[1];
sx q[1];
rz(1.321235) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7119638) q[0];
sx q[0];
rz(-2.5416059) q[0];
sx q[0];
rz(-2.5866051) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2259813) q[2];
sx q[2];
rz(-1.639099) q[2];
sx q[2];
rz(-3.0223522) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3258354) q[1];
sx q[1];
rz(-1.6920631) q[1];
sx q[1];
rz(-1.6440637) q[1];
rz(-pi) q[2];
rz(0.77797555) q[3];
sx q[3];
rz(-0.56853349) q[3];
sx q[3];
rz(-1.83873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2230175) q[2];
sx q[2];
rz(-1.7097946) q[2];
sx q[2];
rz(-2.4187386) q[2];
rz(-0.84960788) q[3];
sx q[3];
rz(-1.1725715) q[3];
sx q[3];
rz(1.0378708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9549114) q[0];
sx q[0];
rz(-2.6246922) q[0];
sx q[0];
rz(-2.5816259) q[0];
rz(0.2050744) q[1];
sx q[1];
rz(-1.8849025) q[1];
sx q[1];
rz(0.29108873) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81924651) q[0];
sx q[0];
rz(-2.8221606) q[0];
sx q[0];
rz(0.56630212) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.79925691) q[2];
sx q[2];
rz(-1.5267258) q[2];
sx q[2];
rz(1.26233) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.12778681) q[1];
sx q[1];
rz(-1.3723137) q[1];
sx q[1];
rz(0.86680331) q[1];
rz(2.970333) q[3];
sx q[3];
rz(-0.75920877) q[3];
sx q[3];
rz(1.1472792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9186972) q[2];
sx q[2];
rz(-1.6632068) q[2];
sx q[2];
rz(2.843294) q[2];
rz(1.9722624) q[3];
sx q[3];
rz(-0.9669286) q[3];
sx q[3];
rz(1.5112618) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30763141) q[0];
sx q[0];
rz(-1.4606322) q[0];
sx q[0];
rz(2.5391915) q[0];
rz(-2.7700453) q[1];
sx q[1];
rz(-1.5448152) q[1];
sx q[1];
rz(2.1098302) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96249798) q[0];
sx q[0];
rz(-2.966605) q[0];
sx q[0];
rz(1.4797158) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7784836) q[2];
sx q[2];
rz(-1.9676932) q[2];
sx q[2];
rz(-2.7092421) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6556485) q[1];
sx q[1];
rz(-1.3820501) q[1];
sx q[1];
rz(-2.6552378) q[1];
rz(-0.57547456) q[3];
sx q[3];
rz(-2.1361975) q[3];
sx q[3];
rz(-2.7421212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.80456698) q[2];
sx q[2];
rz(-2.0402543) q[2];
sx q[2];
rz(-1.6229013) q[2];
rz(0.8484146) q[3];
sx q[3];
rz(-2.266326) q[3];
sx q[3];
rz(0.039610473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0254211) q[0];
sx q[0];
rz(-3.0905368) q[0];
sx q[0];
rz(0.99789944) q[0];
rz(2.8669224) q[1];
sx q[1];
rz(-0.88013595) q[1];
sx q[1];
rz(1.3708699) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28579636) q[0];
sx q[0];
rz(-1.5573475) q[0];
sx q[0];
rz(-0.6310668) q[0];
x q[1];
rz(-2.9878658) q[2];
sx q[2];
rz(-1.9228553) q[2];
sx q[2];
rz(-1.597126) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.22860195) q[1];
sx q[1];
rz(-2.1739042) q[1];
sx q[1];
rz(1.8643537) q[1];
rz(-pi) q[2];
x q[2];
rz(0.73565817) q[3];
sx q[3];
rz(-2.1950985) q[3];
sx q[3];
rz(-1.5746547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2844598) q[2];
sx q[2];
rz(-2.1993115) q[2];
sx q[2];
rz(3.1119463) q[2];
rz(2.5263785) q[3];
sx q[3];
rz(-1.7447724) q[3];
sx q[3];
rz(-0.342338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-3.0520332) q[0];
sx q[0];
rz(-2.679306) q[0];
sx q[0];
rz(-2.3585368) q[0];
rz(2.0217333) q[1];
sx q[1];
rz(-2.4744787) q[1];
sx q[1];
rz(1.3979744) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.072294839) q[0];
sx q[0];
rz(-2.0136478) q[0];
sx q[0];
rz(-2.686727) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2863761) q[2];
sx q[2];
rz(-1.3462974) q[2];
sx q[2];
rz(-0.38849354) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.82247558) q[1];
sx q[1];
rz(-0.90293316) q[1];
sx q[1];
rz(-3.130359) q[1];
rz(-pi) q[2];
rz(-0.78432958) q[3];
sx q[3];
rz(-1.722986) q[3];
sx q[3];
rz(-1.0206211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.29584259) q[2];
sx q[2];
rz(-0.58640277) q[2];
sx q[2];
rz(1.3207377) q[2];
rz(-2.1218421) q[3];
sx q[3];
rz(-1.7176065) q[3];
sx q[3];
rz(1.7216871) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19972292) q[0];
sx q[0];
rz(-2.0482735) q[0];
sx q[0];
rz(0.95440188) q[0];
rz(1.154254) q[1];
sx q[1];
rz(-2.4955165) q[1];
sx q[1];
rz(2.0571713) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0823776) q[0];
sx q[0];
rz(-2.431916) q[0];
sx q[0];
rz(1.6355455) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8516638) q[2];
sx q[2];
rz(-1.8164338) q[2];
sx q[2];
rz(2.7410067) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.22218765) q[1];
sx q[1];
rz(-2.0424941) q[1];
sx q[1];
rz(-3.0604048) q[1];
x q[2];
rz(0.74862759) q[3];
sx q[3];
rz(-1.4350865) q[3];
sx q[3];
rz(-0.77881294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.083100975) q[2];
sx q[2];
rz(-2.9497171) q[2];
sx q[2];
rz(-2.7044738) q[2];
rz(1.7982177) q[3];
sx q[3];
rz(-1.3513887) q[3];
sx q[3];
rz(-0.04537151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.042353543) q[0];
sx q[0];
rz(-2.1411242) q[0];
sx q[0];
rz(-1.4105256) q[0];
rz(-0.22077416) q[1];
sx q[1];
rz(-0.53123728) q[1];
sx q[1];
rz(-0.9332307) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2483053) q[0];
sx q[0];
rz(-2.3961005) q[0];
sx q[0];
rz(-3.1214691) q[0];
rz(-pi) q[1];
x q[1];
rz(0.48671203) q[2];
sx q[2];
rz(-0.67811869) q[2];
sx q[2];
rz(-2.202452) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.5367914) q[1];
sx q[1];
rz(-1.4385421) q[1];
sx q[1];
rz(2.3459646) q[1];
x q[2];
rz(1.7228863) q[3];
sx q[3];
rz(-2.503431) q[3];
sx q[3];
rz(1.9877951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.431388) q[2];
sx q[2];
rz(-0.34222558) q[2];
sx q[2];
rz(2.5101275) q[2];
rz(2.4225875) q[3];
sx q[3];
rz(-1.3471666) q[3];
sx q[3];
rz(0.24148153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63856335) q[0];
sx q[0];
rz(-1.6521492) q[0];
sx q[0];
rz(-1.8304992) q[0];
rz(0.44166625) q[1];
sx q[1];
rz(-1.8001945) q[1];
sx q[1];
rz(1.1460907) q[1];
rz(1.1900525) q[2];
sx q[2];
rz(-1.6002009) q[2];
sx q[2];
rz(2.1717351) q[2];
rz(1.6735703) q[3];
sx q[3];
rz(-2.1379875) q[3];
sx q[3];
rz(0.47982346) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
