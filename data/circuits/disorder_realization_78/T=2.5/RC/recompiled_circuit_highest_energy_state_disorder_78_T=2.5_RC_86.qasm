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
rz(0.084963381) q[0];
sx q[0];
rz(3.4439937) q[0];
sx q[0];
rz(6.2511282) q[0];
rz(3.1337466) q[1];
sx q[1];
rz(3.6454522) q[1];
sx q[1];
rz(7.03581) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0467211) q[0];
sx q[0];
rz(-1.7716452) q[0];
sx q[0];
rz(0.3263536) q[0];
x q[1];
rz(1.9448124) q[2];
sx q[2];
rz(-1.2311282) q[2];
sx q[2];
rz(2.938478) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6258273) q[1];
sx q[1];
rz(-1.2316362) q[1];
sx q[1];
rz(-1.304342) q[1];
x q[2];
rz(-2.0977375) q[3];
sx q[3];
rz(-2.733272) q[3];
sx q[3];
rz(0.645831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.96693119) q[2];
sx q[2];
rz(-1.9661247) q[2];
sx q[2];
rz(-0.82007972) q[2];
rz(-0.44102937) q[3];
sx q[3];
rz(-1.1038154) q[3];
sx q[3];
rz(-0.57344121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.012861982) q[0];
sx q[0];
rz(-0.91330376) q[0];
sx q[0];
rz(-2.7313857) q[0];
rz(2.7606616) q[1];
sx q[1];
rz(-1.610264) q[1];
sx q[1];
rz(-1.358323) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8976645) q[0];
sx q[0];
rz(-1.9635734) q[0];
sx q[0];
rz(-0.94148366) q[0];
rz(-0.83186457) q[2];
sx q[2];
rz(-0.68507776) q[2];
sx q[2];
rz(-2.5545504) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7100923) q[1];
sx q[1];
rz(-2.9572801) q[1];
sx q[1];
rz(-1.9955562) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93102867) q[3];
sx q[3];
rz(-1.4461541) q[3];
sx q[3];
rz(2.3384936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7831948) q[2];
sx q[2];
rz(-2.7326549) q[2];
sx q[2];
rz(0.50296339) q[2];
rz(-1.8424235) q[3];
sx q[3];
rz(-0.75853577) q[3];
sx q[3];
rz(-2.9097596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4101039) q[0];
sx q[0];
rz(-2.6298611) q[0];
sx q[0];
rz(-0.9758392) q[0];
rz(-2.2052235) q[1];
sx q[1];
rz(-1.5543289) q[1];
sx q[1];
rz(0.13793129) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16743037) q[0];
sx q[0];
rz(-2.5241521) q[0];
sx q[0];
rz(-0.90888826) q[0];
rz(-2.7439026) q[2];
sx q[2];
rz(-1.6776591) q[2];
sx q[2];
rz(1.7025089) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8616069) q[1];
sx q[1];
rz(-0.87738887) q[1];
sx q[1];
rz(-0.1711636) q[1];
rz(-pi) q[2];
x q[2];
rz(0.46584399) q[3];
sx q[3];
rz(-1.3261822) q[3];
sx q[3];
rz(-0.96058339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.58933538) q[2];
sx q[2];
rz(-1.6103123) q[2];
sx q[2];
rz(1.4790347) q[2];
rz(2.3432483) q[3];
sx q[3];
rz(-2.2975497) q[3];
sx q[3];
rz(-0.020309694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8524356) q[0];
sx q[0];
rz(-3.030179) q[0];
sx q[0];
rz(1.0668466) q[0];
rz(-1.5486708) q[1];
sx q[1];
rz(-1.8916062) q[1];
sx q[1];
rz(0.41935316) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3187554) q[0];
sx q[0];
rz(-1.7298123) q[0];
sx q[0];
rz(1.8029638) q[0];
rz(1.661395) q[2];
sx q[2];
rz(-1.0444469) q[2];
sx q[2];
rz(-1.3535997) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.6904628) q[1];
sx q[1];
rz(-1.2138543) q[1];
sx q[1];
rz(-3.073077) q[1];
rz(-pi) q[2];
rz(2.7953496) q[3];
sx q[3];
rz(-1.1107363) q[3];
sx q[3];
rz(-1.8025043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.35820094) q[2];
sx q[2];
rz(-2.6322067) q[2];
sx q[2];
rz(-1.5979213) q[2];
rz(1.915043) q[3];
sx q[3];
rz(-1.2164601) q[3];
sx q[3];
rz(-1.8515057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.876038) q[0];
sx q[0];
rz(-2.4198678) q[0];
sx q[0];
rz(0.75620404) q[0];
rz(-1.2528231) q[1];
sx q[1];
rz(-0.93106657) q[1];
sx q[1];
rz(-1.4160215) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.904793) q[0];
sx q[0];
rz(-2.5052245) q[0];
sx q[0];
rz(2.7704253) q[0];
rz(-1.6966693) q[2];
sx q[2];
rz(-1.2294266) q[2];
sx q[2];
rz(1.204513) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.49633138) q[1];
sx q[1];
rz(-1.0320391) q[1];
sx q[1];
rz(-1.5625728) q[1];
x q[2];
rz(2.6719499) q[3];
sx q[3];
rz(-2.7386463) q[3];
sx q[3];
rz(0.5136036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.34117928) q[2];
sx q[2];
rz(-0.74113733) q[2];
sx q[2];
rz(-2.9608012) q[2];
rz(-1.6527269) q[3];
sx q[3];
rz(-1.3988262) q[3];
sx q[3];
rz(-1.6256049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4329231) q[0];
sx q[0];
rz(-1.0754508) q[0];
sx q[0];
rz(-0.80023009) q[0];
rz(-0.284614) q[1];
sx q[1];
rz(-2.0814643) q[1];
sx q[1];
rz(3.0175993) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.894688) q[0];
sx q[0];
rz(-0.13291153) q[0];
sx q[0];
rz(-1.7869496) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0818865) q[2];
sx q[2];
rz(-1.8054031) q[2];
sx q[2];
rz(1.4690746) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.89470321) q[1];
sx q[1];
rz(-1.8754301) q[1];
sx q[1];
rz(1.4851634) q[1];
x q[2];
rz(-2.0923397) q[3];
sx q[3];
rz(-0.66415706) q[3];
sx q[3];
rz(0.60598999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2064712) q[2];
sx q[2];
rz(-1.7266885) q[2];
sx q[2];
rz(-3.0653811) q[2];
rz(-1.8501836) q[3];
sx q[3];
rz(-2.0468057) q[3];
sx q[3];
rz(-2.5230303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9718219) q[0];
sx q[0];
rz(-1.5795647) q[0];
sx q[0];
rz(2.3845657) q[0];
rz(2.3454759) q[1];
sx q[1];
rz(-1.7346104) q[1];
sx q[1];
rz(0.98181358) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9011089) q[0];
sx q[0];
rz(-2.3496858) q[0];
sx q[0];
rz(-1.1674985) q[0];
rz(0.99346353) q[2];
sx q[2];
rz(-0.12778035) q[2];
sx q[2];
rz(-1.8005231) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.031739) q[1];
sx q[1];
rz(-0.64424911) q[1];
sx q[1];
rz(2.2883313) q[1];
rz(-pi) q[2];
x q[2];
rz(0.14431868) q[3];
sx q[3];
rz(-1.2929521) q[3];
sx q[3];
rz(1.8622189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7299399) q[2];
sx q[2];
rz(-0.81092683) q[2];
sx q[2];
rz(-0.5438424) q[2];
rz(-0.020921556) q[3];
sx q[3];
rz(-1.7048416) q[3];
sx q[3];
rz(-0.81834832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92886096) q[0];
sx q[0];
rz(-0.91101557) q[0];
sx q[0];
rz(0.62335706) q[0];
rz(0.32132545) q[1];
sx q[1];
rz(-2.5265381) q[1];
sx q[1];
rz(-2.8772433) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0118222) q[0];
sx q[0];
rz(-1.7826102) q[0];
sx q[0];
rz(2.2925633) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2351722) q[2];
sx q[2];
rz(-1.610029) q[2];
sx q[2];
rz(0.36612636) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.26786727) q[1];
sx q[1];
rz(-1.8494938) q[1];
sx q[1];
rz(-2.3257117) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8885625) q[3];
sx q[3];
rz(-1.4444286) q[3];
sx q[3];
rz(-1.7333836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.136772) q[2];
sx q[2];
rz(-2.4592168) q[2];
sx q[2];
rz(2.6668059) q[2];
rz(-2.7759806) q[3];
sx q[3];
rz(-0.48706278) q[3];
sx q[3];
rz(0.18317187) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8356165) q[0];
sx q[0];
rz(-0.37537471) q[0];
sx q[0];
rz(2.6557652) q[0];
rz(-1.0659418) q[1];
sx q[1];
rz(-1.424788) q[1];
sx q[1];
rz(-0.33445439) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9850824) q[0];
sx q[0];
rz(-1.2691536) q[0];
sx q[0];
rz(1.7443875) q[0];
rz(-pi) q[1];
x q[1];
rz(0.43844079) q[2];
sx q[2];
rz(-1.4435569) q[2];
sx q[2];
rz(-1.3194989) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.124467) q[1];
sx q[1];
rz(-1.7199868) q[1];
sx q[1];
rz(-2.6890432) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.79631898) q[3];
sx q[3];
rz(-1.9069172) q[3];
sx q[3];
rz(-1.4848061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7740384) q[2];
sx q[2];
rz(-2.2944821) q[2];
sx q[2];
rz(-0.96662194) q[2];
rz(-1.4878081) q[3];
sx q[3];
rz(-2.8130468) q[3];
sx q[3];
rz(0.070913471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-1.3286572) q[0];
sx q[0];
rz(-2.2870977) q[0];
sx q[0];
rz(0.27035126) q[0];
rz(-1.5125754) q[1];
sx q[1];
rz(-0.83643475) q[1];
sx q[1];
rz(2.5152452) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26890818) q[0];
sx q[0];
rz(-1.88694) q[0];
sx q[0];
rz(2.288398) q[0];
rz(-pi) q[1];
x q[1];
rz(0.43105189) q[2];
sx q[2];
rz(-0.97857514) q[2];
sx q[2];
rz(2.2811802) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8423564) q[1];
sx q[1];
rz(-1.6116214) q[1];
sx q[1];
rz(-3.051332) q[1];
rz(-pi) q[2];
rz(2.7553431) q[3];
sx q[3];
rz(-1.874141) q[3];
sx q[3];
rz(-1.7084427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.35499972) q[2];
sx q[2];
rz(-2.2248416) q[2];
sx q[2];
rz(1.5691441) q[2];
rz(-0.7507503) q[3];
sx q[3];
rz(-0.34199491) q[3];
sx q[3];
rz(-2.2641838) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5759721) q[0];
sx q[0];
rz(-1.2979869) q[0];
sx q[0];
rz(-0.57814231) q[0];
rz(1.3427973) q[1];
sx q[1];
rz(-2.8248351) q[1];
sx q[1];
rz(0.14541365) q[1];
rz(-0.00040309488) q[2];
sx q[2];
rz(-3.0192791) q[2];
sx q[2];
rz(3.064439) q[2];
rz(-0.25788767) q[3];
sx q[3];
rz(-1.5189239) q[3];
sx q[3];
rz(2.2699395) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
