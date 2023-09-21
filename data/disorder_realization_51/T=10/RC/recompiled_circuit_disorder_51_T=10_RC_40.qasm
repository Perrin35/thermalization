OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.59453073) q[0];
sx q[0];
rz(-1.1214331) q[0];
sx q[0];
rz(0.18145951) q[0];
rz(2.060086) q[1];
sx q[1];
rz(5.6097538) q[1];
sx q[1];
rz(14.654832) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5232789) q[0];
sx q[0];
rz(-1.7821454) q[0];
sx q[0];
rz(2.0413415) q[0];
x q[1];
rz(-0.62280957) q[2];
sx q[2];
rz(-1.8138759) q[2];
sx q[2];
rz(0.18261766) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.44838201) q[1];
sx q[1];
rz(-0.795185) q[1];
sx q[1];
rz(-0.27547142) q[1];
rz(-pi) q[2];
rz(-2.9795322) q[3];
sx q[3];
rz(-1.2352304) q[3];
sx q[3];
rz(3.0367362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9782605) q[2];
sx q[2];
rz(-1.1085199) q[2];
sx q[2];
rz(-1.367761) q[2];
rz(-2.1286428) q[3];
sx q[3];
rz(-2.2949341) q[3];
sx q[3];
rz(-3.0701385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(2.0435836) q[0];
sx q[0];
rz(-1.4494891) q[0];
sx q[0];
rz(0.59511551) q[0];
rz(-2.0960506) q[1];
sx q[1];
rz(-1.4141934) q[1];
sx q[1];
rz(-1.6275303) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94905108) q[0];
sx q[0];
rz(-2.4727614) q[0];
sx q[0];
rz(2.2632954) q[0];
x q[1];
rz(1.6018277) q[2];
sx q[2];
rz(-2.4040262) q[2];
sx q[2];
rz(-0.34027983) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.82799998) q[1];
sx q[1];
rz(-0.26121059) q[1];
sx q[1];
rz(-2.2952495) q[1];
x q[2];
rz(3.0916445) q[3];
sx q[3];
rz(-1.6729828) q[3];
sx q[3];
rz(3.0199043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4864768) q[2];
sx q[2];
rz(-2.1003508) q[2];
sx q[2];
rz(-0.79616037) q[2];
rz(-2.1697309) q[3];
sx q[3];
rz(-0.704851) q[3];
sx q[3];
rz(-2.9698353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7650836) q[0];
sx q[0];
rz(-2.3628545) q[0];
sx q[0];
rz(0.064095108) q[0];
rz(-0.31072101) q[1];
sx q[1];
rz(-1.6711845) q[1];
sx q[1];
rz(-1.6832738) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4401535) q[0];
sx q[0];
rz(-2.9868786) q[0];
sx q[0];
rz(-1.2604777) q[0];
x q[1];
rz(-0.65285039) q[2];
sx q[2];
rz(-1.0703501) q[2];
sx q[2];
rz(2.1476538) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8628143) q[1];
sx q[1];
rz(-1.9796951) q[1];
sx q[1];
rz(1.5120974) q[1];
x q[2];
rz(-1.5536669) q[3];
sx q[3];
rz(-1.8855842) q[3];
sx q[3];
rz(2.0103612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1594499) q[2];
sx q[2];
rz(-0.84656707) q[2];
sx q[2];
rz(2.4625051) q[2];
rz(-1.3726161) q[3];
sx q[3];
rz(-1.2981828) q[3];
sx q[3];
rz(2.1931271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8452334) q[0];
sx q[0];
rz(-0.56931749) q[0];
sx q[0];
rz(1.1244208) q[0];
rz(-1.2202948) q[1];
sx q[1];
rz(-1.1923469) q[1];
sx q[1];
rz(-1.4845928) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4666045) q[0];
sx q[0];
rz(-3.0825442) q[0];
sx q[0];
rz(2.1112291) q[0];
rz(-pi) q[1];
rz(-2.0882656) q[2];
sx q[2];
rz(-1.4033068) q[2];
sx q[2];
rz(2.6094764) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2034519) q[1];
sx q[1];
rz(-1.005299) q[1];
sx q[1];
rz(0.52312619) q[1];
rz(0.72317601) q[3];
sx q[3];
rz(-1.3190862) q[3];
sx q[3];
rz(2.6808006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0174039) q[2];
sx q[2];
rz(-1.2839395) q[2];
sx q[2];
rz(1.0162639) q[2];
rz(-1.4034363) q[3];
sx q[3];
rz(-1.6442464) q[3];
sx q[3];
rz(-2.0531634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-0.50773412) q[0];
sx q[0];
rz(-1.8554747) q[0];
sx q[0];
rz(-2.741709) q[0];
rz(-1.9790861) q[1];
sx q[1];
rz(-1.8116654) q[1];
sx q[1];
rz(-2.9679325) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0608873) q[0];
sx q[0];
rz(-1.4982491) q[0];
sx q[0];
rz(1.5425496) q[0];
rz(-pi) q[1];
rz(0.47541754) q[2];
sx q[2];
rz(-2.7042275) q[2];
sx q[2];
rz(-1.4471444) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.793321) q[1];
sx q[1];
rz(-1.4655359) q[1];
sx q[1];
rz(0.028491032) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7791622) q[3];
sx q[3];
rz(-1.2456129) q[3];
sx q[3];
rz(1.1680101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.48012039) q[2];
sx q[2];
rz(-1.1602594) q[2];
sx q[2];
rz(-0.76812569) q[2];
rz(2.2875732) q[3];
sx q[3];
rz(-1.4203527) q[3];
sx q[3];
rz(0.97222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0356692) q[0];
sx q[0];
rz(-0.24374715) q[0];
sx q[0];
rz(-1.6712028) q[0];
rz(-2.629783) q[1];
sx q[1];
rz(-0.51135951) q[1];
sx q[1];
rz(-1.1434198) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91771942) q[0];
sx q[0];
rz(-1.6544764) q[0];
sx q[0];
rz(2.9907945) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.60322275) q[2];
sx q[2];
rz(-1.2957186) q[2];
sx q[2];
rz(-0.075866931) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1089576) q[1];
sx q[1];
rz(-2.7320478) q[1];
sx q[1];
rz(2.0202548) q[1];
x q[2];
rz(2.1120758) q[3];
sx q[3];
rz(-1.3068849) q[3];
sx q[3];
rz(1.2332066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.64289552) q[2];
sx q[2];
rz(-0.84474793) q[2];
sx q[2];
rz(1.5303622) q[2];
rz(-1.4536084) q[3];
sx q[3];
rz(-1.0691103) q[3];
sx q[3];
rz(-0.24916515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11944184) q[0];
sx q[0];
rz(-2.3972153) q[0];
sx q[0];
rz(1.7013593) q[0];
rz(2.4123736) q[1];
sx q[1];
rz(-1.9897285) q[1];
sx q[1];
rz(1.1332606) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0798577) q[0];
sx q[0];
rz(-1.2485866) q[0];
sx q[0];
rz(1.6155924) q[0];
rz(-0.71224833) q[2];
sx q[2];
rz(-1.5812261) q[2];
sx q[2];
rz(-2.1005837) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.82517351) q[1];
sx q[1];
rz(-2.2146041) q[1];
sx q[1];
rz(2.7979147) q[1];
rz(-pi) q[2];
rz(0.53448581) q[3];
sx q[3];
rz(-2.839698) q[3];
sx q[3];
rz(1.2045977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7363654) q[2];
sx q[2];
rz(-1.1232802) q[2];
sx q[2];
rz(-2.6531632) q[2];
rz(-1.3119665) q[3];
sx q[3];
rz(-0.043881504) q[3];
sx q[3];
rz(-2.5002938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064780386) q[0];
sx q[0];
rz(-1.8490054) q[0];
sx q[0];
rz(-0.07117614) q[0];
rz(-3.1094303) q[1];
sx q[1];
rz(-1.3379438) q[1];
sx q[1];
rz(1.2088998) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6541518) q[0];
sx q[0];
rz(-1.7983266) q[0];
sx q[0];
rz(0.80765101) q[0];
rz(2.9572763) q[2];
sx q[2];
rz(-1.3375998) q[2];
sx q[2];
rz(1.2554982) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3370918) q[1];
sx q[1];
rz(-0.13131222) q[1];
sx q[1];
rz(1.4485703) q[1];
rz(-0.62187059) q[3];
sx q[3];
rz(-1.2688046) q[3];
sx q[3];
rz(0.57884502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9341087) q[2];
sx q[2];
rz(-2.948163) q[2];
sx q[2];
rz(-1.570638) q[2];
rz(2.2682244) q[3];
sx q[3];
rz(-1.4040754) q[3];
sx q[3];
rz(-0.35203716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97380012) q[0];
sx q[0];
rz(-1.5252824) q[0];
sx q[0];
rz(-2.8299676) q[0];
rz(-2.3198126) q[1];
sx q[1];
rz(-2.5506134) q[1];
sx q[1];
rz(1.5100381) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25699297) q[0];
sx q[0];
rz(-1.25367) q[0];
sx q[0];
rz(-1.3604128) q[0];
x q[1];
rz(1.2276149) q[2];
sx q[2];
rz(-1.7556264) q[2];
sx q[2];
rz(-0.96175805) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9540625) q[1];
sx q[1];
rz(-0.64779753) q[1];
sx q[1];
rz(2.8026583) q[1];
rz(-pi) q[2];
rz(0.23267965) q[3];
sx q[3];
rz(-1.7494697) q[3];
sx q[3];
rz(1.0089547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8682378) q[2];
sx q[2];
rz(-0.53871012) q[2];
sx q[2];
rz(2.3256425) q[2];
rz(-0.50968918) q[3];
sx q[3];
rz(-0.78521252) q[3];
sx q[3];
rz(1.2379439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2508535) q[0];
sx q[0];
rz(-1.7128523) q[0];
sx q[0];
rz(-2.9113286) q[0];
rz(2.5157805) q[1];
sx q[1];
rz(-0.9451378) q[1];
sx q[1];
rz(2.4831916) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2021159) q[0];
sx q[0];
rz(-0.80047031) q[0];
sx q[0];
rz(0.70864422) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4670829) q[2];
sx q[2];
rz(-1.837338) q[2];
sx q[2];
rz(-0.36905497) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8734332) q[1];
sx q[1];
rz(-0.50926103) q[1];
sx q[1];
rz(-1.2251653) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6055466) q[3];
sx q[3];
rz(-1.8744161) q[3];
sx q[3];
rz(2.924078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.6752424) q[2];
sx q[2];
rz(-2.3190976) q[2];
sx q[2];
rz(-0.33774439) q[2];
rz(1.0234458) q[3];
sx q[3];
rz(-1.3600391) q[3];
sx q[3];
rz(-1.2158998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.52453775) q[0];
sx q[0];
rz(-2.239997) q[0];
sx q[0];
rz(-0.064185113) q[0];
rz(1.9032003) q[1];
sx q[1];
rz(-1.4307784) q[1];
sx q[1];
rz(1.4684114) q[1];
rz(-1.3255618) q[2];
sx q[2];
rz(-2.500848) q[2];
sx q[2];
rz(1.4304888) q[2];
rz(-0.099048793) q[3];
sx q[3];
rz(-2.5669813) q[3];
sx q[3];
rz(3.0909227) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];