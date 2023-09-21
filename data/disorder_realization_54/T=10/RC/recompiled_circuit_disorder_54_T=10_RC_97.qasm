OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.78046507) q[0];
sx q[0];
rz(3.830885) q[0];
sx q[0];
rz(9.7552714) q[0];
rz(-2.8117872) q[1];
sx q[1];
rz(-2.2916315) q[1];
sx q[1];
rz(-0.70911521) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7666727) q[0];
sx q[0];
rz(-1.2776889) q[0];
sx q[0];
rz(-0.93107443) q[0];
x q[1];
rz(1.4719226) q[2];
sx q[2];
rz(-2.8266202) q[2];
sx q[2];
rz(-1.0613943) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6735437) q[1];
sx q[1];
rz(-1.1357726) q[1];
sx q[1];
rz(2.166964) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1756644) q[3];
sx q[3];
rz(-0.60797193) q[3];
sx q[3];
rz(2.0917497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7314529) q[2];
sx q[2];
rz(-2.7316766) q[2];
sx q[2];
rz(-1.5343792) q[2];
rz(0.93506995) q[3];
sx q[3];
rz(-1.8362703) q[3];
sx q[3];
rz(2.3527761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-1.4193029) q[0];
sx q[0];
rz(-0.062148217) q[0];
sx q[0];
rz(0.62227917) q[0];
rz(2.9653446) q[1];
sx q[1];
rz(-1.9269678) q[1];
sx q[1];
rz(2.2252749) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57740649) q[0];
sx q[0];
rz(-0.86125492) q[0];
sx q[0];
rz(-0.90155154) q[0];
x q[1];
rz(-1.4257405) q[2];
sx q[2];
rz(-2.3999891) q[2];
sx q[2];
rz(1.2765826) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2434477) q[1];
sx q[1];
rz(-1.2444278) q[1];
sx q[1];
rz(1.2136202) q[1];
rz(-pi) q[2];
rz(1.968859) q[3];
sx q[3];
rz(-1.5916087) q[3];
sx q[3];
rz(-1.2956937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9006485) q[2];
sx q[2];
rz(-0.99664656) q[2];
sx q[2];
rz(2.3201578) q[2];
rz(-0.017283043) q[3];
sx q[3];
rz(-2.343943) q[3];
sx q[3];
rz(2.3582874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18773742) q[0];
sx q[0];
rz(-1.563235) q[0];
sx q[0];
rz(-2.1333372) q[0];
rz(-3.1058274) q[1];
sx q[1];
rz(-1.6556031) q[1];
sx q[1];
rz(0.52454138) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5265822) q[0];
sx q[0];
rz(-2.9752762) q[0];
sx q[0];
rz(1.7517356) q[0];
x q[1];
rz(-2.4630765) q[2];
sx q[2];
rz(-1.3635474) q[2];
sx q[2];
rz(3.0845272) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4265392) q[1];
sx q[1];
rz(-1.5476777) q[1];
sx q[1];
rz(-0.061736488) q[1];
rz(1.2479765) q[3];
sx q[3];
rz(-2.7016692) q[3];
sx q[3];
rz(-0.80252121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.77164578) q[2];
sx q[2];
rz(-1.5180072) q[2];
sx q[2];
rz(1.4952205) q[2];
rz(1.3211936) q[3];
sx q[3];
rz(-1.4097872) q[3];
sx q[3];
rz(2.7041919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8111073) q[0];
sx q[0];
rz(-2.2225668) q[0];
sx q[0];
rz(3.0526429) q[0];
rz(2.6308909) q[1];
sx q[1];
rz(-2.316244) q[1];
sx q[1];
rz(-0.68960062) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21259637) q[0];
sx q[0];
rz(-1.4082452) q[0];
sx q[0];
rz(-2.8430804) q[0];
rz(-pi) q[1];
x q[1];
rz(0.56324048) q[2];
sx q[2];
rz(-13/(3*pi)) q[2];
sx q[2];
rz(-2.9134977) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.888962) q[1];
sx q[1];
rz(-1.2446212) q[1];
sx q[1];
rz(2.3701282) q[1];
rz(-pi) q[2];
rz(-2.6570286) q[3];
sx q[3];
rz(-1.6996517) q[3];
sx q[3];
rz(1.7948922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4758063) q[2];
sx q[2];
rz(-2.4643354) q[2];
sx q[2];
rz(-0.62292567) q[2];
rz(1.1359435) q[3];
sx q[3];
rz(-0.1853075) q[3];
sx q[3];
rz(-0.46666551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85161197) q[0];
sx q[0];
rz(-1.8988134) q[0];
sx q[0];
rz(2.5581397) q[0];
rz(1.1460229) q[1];
sx q[1];
rz(-1.6505516) q[1];
sx q[1];
rz(1.4978283) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72901112) q[0];
sx q[0];
rz(-0.81775613) q[0];
sx q[0];
rz(1.9304995) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6696817) q[2];
sx q[2];
rz(-0.77557287) q[2];
sx q[2];
rz(2.0688187) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.46206611) q[1];
sx q[1];
rz(-2.4375009) q[1];
sx q[1];
rz(0.042298869) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4086401) q[3];
sx q[3];
rz(-1.0609396) q[3];
sx q[3];
rz(1.7061403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4862165) q[2];
sx q[2];
rz(-1.5950173) q[2];
sx q[2];
rz(1.5779457) q[2];
rz(-2.2359713) q[3];
sx q[3];
rz(-0.27035299) q[3];
sx q[3];
rz(0.63846987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6435796) q[0];
sx q[0];
rz(-1.6828515) q[0];
sx q[0];
rz(0.046982732) q[0];
rz(-0.14818305) q[1];
sx q[1];
rz(-2.2427185) q[1];
sx q[1];
rz(-1.7061589) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5934138) q[0];
sx q[0];
rz(-0.81028623) q[0];
sx q[0];
rz(-0.20472783) q[0];
x q[1];
rz(3.0827423) q[2];
sx q[2];
rz(-2.119679) q[2];
sx q[2];
rz(-0.6097874) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.50373494) q[1];
sx q[1];
rz(-1.4459472) q[1];
sx q[1];
rz(1.9810852) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.027408882) q[3];
sx q[3];
rz(-1.3331183) q[3];
sx q[3];
rz(-2.4159367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0662213) q[2];
sx q[2];
rz(-1.0152638) q[2];
sx q[2];
rz(-3.0701239) q[2];
rz(-1.6890769) q[3];
sx q[3];
rz(-1.444933) q[3];
sx q[3];
rz(2.215109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8554095) q[0];
sx q[0];
rz(-1.3278642) q[0];
sx q[0];
rz(-0.57762161) q[0];
rz(-1.8619934) q[1];
sx q[1];
rz(-1.7956693) q[1];
sx q[1];
rz(-2.1320027) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5492064) q[0];
sx q[0];
rz(-1.2388595) q[0];
sx q[0];
rz(-0.50360002) q[0];
rz(-2.6206714) q[2];
sx q[2];
rz(-2.705057) q[2];
sx q[2];
rz(1.5182564) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2254667) q[1];
sx q[1];
rz(-1.2327317) q[1];
sx q[1];
rz(2.613693) q[1];
rz(-pi) q[2];
rz(-2.3791802) q[3];
sx q[3];
rz(-1.1165285) q[3];
sx q[3];
rz(-1.4095969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0916831) q[2];
sx q[2];
rz(-0.35353264) q[2];
sx q[2];
rz(-1.9419149) q[2];
rz(2.6692634) q[3];
sx q[3];
rz(-1.6475369) q[3];
sx q[3];
rz(-2.1388617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49884477) q[0];
sx q[0];
rz(-1.5002748) q[0];
sx q[0];
rz(-2.902466) q[0];
rz(-0.7827951) q[1];
sx q[1];
rz(-2.6233964) q[1];
sx q[1];
rz(2.696864) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.037576588) q[0];
sx q[0];
rz(-0.74066478) q[0];
sx q[0];
rz(0.91233493) q[0];
rz(-pi) q[1];
rz(2.4773981) q[2];
sx q[2];
rz(-1.0401298) q[2];
sx q[2];
rz(1.9265837) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4152894) q[1];
sx q[1];
rz(-2.2518034) q[1];
sx q[1];
rz(-0.16420941) q[1];
x q[2];
rz(-1.4374251) q[3];
sx q[3];
rz(-1.1629472) q[3];
sx q[3];
rz(-2.9123902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7198221) q[2];
sx q[2];
rz(-2.3193216) q[2];
sx q[2];
rz(-0.81531173) q[2];
rz(-2.1067965) q[3];
sx q[3];
rz(-1.5765604) q[3];
sx q[3];
rz(-1.8243272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8274882) q[0];
sx q[0];
rz(-1.0359456) q[0];
sx q[0];
rz(-0.28717336) q[0];
rz(0.18889591) q[1];
sx q[1];
rz(-0.72383988) q[1];
sx q[1];
rz(-2.8093991) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0567386) q[0];
sx q[0];
rz(-2.2915654) q[0];
sx q[0];
rz(0.42215729) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.8719445) q[2];
sx q[2];
rz(-1.5351864) q[2];
sx q[2];
rz(1.2138838) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0041973) q[1];
sx q[1];
rz(-1.9837556) q[1];
sx q[1];
rz(1.8505627) q[1];
rz(-pi) q[2];
rz(-2.513047) q[3];
sx q[3];
rz(-2.5124031) q[3];
sx q[3];
rz(2.7043846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.93280783) q[2];
sx q[2];
rz(-2.1366426) q[2];
sx q[2];
rz(-0.83958158) q[2];
rz(1.8509289) q[3];
sx q[3];
rz(-1.3495812) q[3];
sx q[3];
rz(0.65657842) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0861417) q[0];
sx q[0];
rz(-2.3738326) q[0];
sx q[0];
rz(1.0593876) q[0];
rz(2.1620031) q[1];
sx q[1];
rz(-1.9444119) q[1];
sx q[1];
rz(-0.25451452) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1232676) q[0];
sx q[0];
rz(-1.3402481) q[0];
sx q[0];
rz(0.075320764) q[0];
x q[1];
rz(2.7024686) q[2];
sx q[2];
rz(-0.52999485) q[2];
sx q[2];
rz(-2.8119171) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3855615) q[1];
sx q[1];
rz(-0.33547151) q[1];
sx q[1];
rz(2.8940593) q[1];
rz(-pi) q[2];
rz(-0.20930807) q[3];
sx q[3];
rz(-2.6946687) q[3];
sx q[3];
rz(-0.79093864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0467726) q[2];
sx q[2];
rz(-2.5929055) q[2];
sx q[2];
rz(-2.0937031) q[2];
rz(-1.7808328) q[3];
sx q[3];
rz(-2.4436617) q[3];
sx q[3];
rz(-0.60539436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1142674) q[0];
sx q[0];
rz(-2.0226759) q[0];
sx q[0];
rz(-0.080060536) q[0];
rz(0.36021532) q[1];
sx q[1];
rz(-1.6811014) q[1];
sx q[1];
rz(-0.95492687) q[1];
rz(2.0965626) q[2];
sx q[2];
rz(-1.3617931) q[2];
sx q[2];
rz(1.1522273) q[2];
rz(1.8923106) q[3];
sx q[3];
rz(-2.5026863) q[3];
sx q[3];
rz(1.7659059) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
