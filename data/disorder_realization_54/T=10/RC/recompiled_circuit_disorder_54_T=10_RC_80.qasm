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
rz(-2.4523003) q[0];
sx q[0];
rz(0.33049345) q[0];
rz(0.32980546) q[1];
sx q[1];
rz(5.4332241) q[1];
sx q[1];
rz(10.133893) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7666727) q[0];
sx q[0];
rz(-1.2776889) q[0];
sx q[0];
rz(0.93107443) q[0];
rz(-pi) q[1];
rz(1.8843295) q[2];
sx q[2];
rz(-1.6013813) q[2];
sx q[2];
rz(-0.41536301) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3173649) q[1];
sx q[1];
rz(-2.1050276) q[1];
sx q[1];
rz(-0.51170106) q[1];
rz(-pi) q[2];
rz(2.0907008) q[3];
sx q[3];
rz(-1.2399779) q[3];
sx q[3];
rz(2.1045121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7314529) q[2];
sx q[2];
rz(-2.7316766) q[2];
sx q[2];
rz(-1.6072134) q[2];
rz(-0.93506995) q[3];
sx q[3];
rz(-1.3053223) q[3];
sx q[3];
rz(-0.7888166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7222897) q[0];
sx q[0];
rz(-0.062148217) q[0];
sx q[0];
rz(2.5193135) q[0];
rz(-0.17624804) q[1];
sx q[1];
rz(-1.2146249) q[1];
sx q[1];
rz(0.91631779) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8371671) q[0];
sx q[0];
rz(-0.93351782) q[0];
sx q[0];
rz(2.5159555) q[0];
rz(2.3071438) q[2];
sx q[2];
rz(-1.6685899) q[2];
sx q[2];
rz(-0.18690878) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1044836) q[1];
sx q[1];
rz(-0.47905211) q[1];
sx q[1];
rz(-2.3399809) q[1];
rz(3.1190156) q[3];
sx q[3];
rz(-1.1728247) q[3];
sx q[3];
rz(2.857739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9006485) q[2];
sx q[2];
rz(-0.99664656) q[2];
sx q[2];
rz(0.82143482) q[2];
rz(-3.1243096) q[3];
sx q[3];
rz(-2.343943) q[3];
sx q[3];
rz(-2.3582874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9538552) q[0];
sx q[0];
rz(-1.5783577) q[0];
sx q[0];
rz(-2.1333372) q[0];
rz(3.1058274) q[1];
sx q[1];
rz(-1.6556031) q[1];
sx q[1];
rz(-0.52454138) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4315955) q[0];
sx q[0];
rz(-1.734373) q[0];
sx q[0];
rz(0.030199108) q[0];
x q[1];
rz(2.8183297) q[2];
sx q[2];
rz(-0.70463902) q[2];
sx q[2];
rz(1.3779674) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7150535) q[1];
sx q[1];
rz(-1.5939149) q[1];
sx q[1];
rz(0.061736488) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9906304) q[3];
sx q[3];
rz(-1.7063147) q[3];
sx q[3];
rz(0.47437048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.77164578) q[2];
sx q[2];
rz(-1.5180072) q[2];
sx q[2];
rz(1.4952205) q[2];
rz(-1.8203991) q[3];
sx q[3];
rz(-1.4097872) q[3];
sx q[3];
rz(-0.43740073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8111073) q[0];
sx q[0];
rz(-0.91902584) q[0];
sx q[0];
rz(3.0526429) q[0];
rz(-2.6308909) q[1];
sx q[1];
rz(-0.82534868) q[1];
sx q[1];
rz(-0.68960062) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8331497) q[0];
sx q[0];
rz(-1.8652548) q[0];
sx q[0];
rz(-1.4008646) q[0];
rz(-1.7961411) q[2];
sx q[2];
rz(-1.0190522) q[2];
sx q[2];
rz(-1.4622886) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1414813) q[1];
sx q[1];
rz(-2.3173216) q[1];
sx q[1];
rz(-0.4517171) q[1];
rz(-pi) q[2];
rz(0.48456405) q[3];
sx q[3];
rz(-1.441941) q[3];
sx q[3];
rz(-1.7948922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.66578635) q[2];
sx q[2];
rz(-2.4643354) q[2];
sx q[2];
rz(-2.518667) q[2];
rz(-1.1359435) q[3];
sx q[3];
rz(-2.9562852) q[3];
sx q[3];
rz(2.6749271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.2899807) q[0];
sx q[0];
rz(-1.8988134) q[0];
sx q[0];
rz(2.5581397) q[0];
rz(-1.9955697) q[1];
sx q[1];
rz(-1.4910411) q[1];
sx q[1];
rz(1.6437644) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72901112) q[0];
sx q[0];
rz(-2.3238365) q[0];
sx q[0];
rz(1.9304995) q[0];
rz(-pi) q[1];
rz(-0.096502467) q[2];
sx q[2];
rz(-2.3415903) q[2];
sx q[2];
rz(-0.9347136) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0651107) q[1];
sx q[1];
rz(-1.5434192) q[1];
sx q[1];
rz(-0.70365023) q[1];
x q[2];
rz(2.6260914) q[3];
sx q[3];
rz(-1.4294129) q[3];
sx q[3];
rz(-3.0859214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4862165) q[2];
sx q[2];
rz(-1.5465753) q[2];
sx q[2];
rz(1.5636469) q[2];
rz(0.90562138) q[3];
sx q[3];
rz(-2.8712397) q[3];
sx q[3];
rz(-0.63846987) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49801302) q[0];
sx q[0];
rz(-1.4587412) q[0];
sx q[0];
rz(3.0946099) q[0];
rz(-2.9934096) q[1];
sx q[1];
rz(-0.89887416) q[1];
sx q[1];
rz(1.4354338) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5934138) q[0];
sx q[0];
rz(-0.81028623) q[0];
sx q[0];
rz(0.20472783) q[0];
x q[1];
rz(-1.6666744) q[2];
sx q[2];
rz(-2.5898858) q[2];
sx q[2];
rz(-2.4193537) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.50373494) q[1];
sx q[1];
rz(-1.4459472) q[1];
sx q[1];
rz(1.1605074) q[1];
rz(-0.027408882) q[3];
sx q[3];
rz(-1.8084744) q[3];
sx q[3];
rz(2.4159367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0753714) q[2];
sx q[2];
rz(-2.1263289) q[2];
sx q[2];
rz(0.071468778) q[2];
rz(1.6890769) q[3];
sx q[3];
rz(-1.444933) q[3];
sx q[3];
rz(0.92648363) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28618318) q[0];
sx q[0];
rz(-1.8137285) q[0];
sx q[0];
rz(-0.57762161) q[0];
rz(1.8619934) q[1];
sx q[1];
rz(-1.3459233) q[1];
sx q[1];
rz(-2.1320027) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9855232) q[0];
sx q[0];
rz(-2.0445604) q[0];
sx q[0];
rz(1.1958634) q[0];
rz(2.7570653) q[2];
sx q[2];
rz(-1.782801) q[2];
sx q[2];
rz(-0.42696135) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.91612591) q[1];
sx q[1];
rz(-1.2327317) q[1];
sx q[1];
rz(-2.613693) q[1];
x q[2];
rz(2.3791802) q[3];
sx q[3];
rz(-1.1165285) q[3];
sx q[3];
rz(-1.7319958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0499095) q[2];
sx q[2];
rz(-2.78806) q[2];
sx q[2];
rz(1.1996777) q[2];
rz(-2.6692634) q[3];
sx q[3];
rz(-1.6475369) q[3];
sx q[3];
rz(2.1388617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6427479) q[0];
sx q[0];
rz(-1.5002748) q[0];
sx q[0];
rz(-2.902466) q[0];
rz(-0.7827951) q[1];
sx q[1];
rz(-0.51819623) q[1];
sx q[1];
rz(-2.696864) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3702104) q[0];
sx q[0];
rz(-2.133773) q[0];
sx q[0];
rz(0.51008205) q[0];
x q[1];
rz(-0.66419454) q[2];
sx q[2];
rz(-1.0401298) q[2];
sx q[2];
rz(1.9265837) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4689444) q[1];
sx q[1];
rz(-0.69744195) q[1];
sx q[1];
rz(1.3717321) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.41110699) q[3];
sx q[3];
rz(-1.4484222) q[3];
sx q[3];
rz(-1.3947595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7198221) q[2];
sx q[2];
rz(-2.3193216) q[2];
sx q[2];
rz(-0.81531173) q[2];
rz(-2.1067965) q[3];
sx q[3];
rz(-1.5650322) q[3];
sx q[3];
rz(-1.3172654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8274882) q[0];
sx q[0];
rz(-2.1056471) q[0];
sx q[0];
rz(-0.28717336) q[0];
rz(-2.9526967) q[1];
sx q[1];
rz(-2.4177528) q[1];
sx q[1];
rz(-0.33219355) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9438292) q[0];
sx q[0];
rz(-1.257886) q[0];
sx q[0];
rz(-2.3373332) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.095093) q[2];
sx q[2];
rz(-2.2691155) q[2];
sx q[2];
rz(2.7547714) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5933983) q[1];
sx q[1];
rz(-1.3151004) q[1];
sx q[1];
rz(-0.42773186) q[1];
rz(-pi) q[2];
x q[2];
rz(1.166415) q[3];
sx q[3];
rz(-1.0746733) q[3];
sx q[3];
rz(-0.29508428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2087848) q[2];
sx q[2];
rz(-1.00495) q[2];
sx q[2];
rz(0.83958158) q[2];
rz(-1.8509289) q[3];
sx q[3];
rz(-1.7920114) q[3];
sx q[3];
rz(-2.4850142) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.055450913) q[0];
sx q[0];
rz(-2.3738326) q[0];
sx q[0];
rz(-1.0593876) q[0];
rz(2.1620031) q[1];
sx q[1];
rz(-1.9444119) q[1];
sx q[1];
rz(-0.25451452) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6063639) q[0];
sx q[0];
rz(-1.6441206) q[0];
sx q[0];
rz(1.3396157) q[0];
rz(-1.326667) q[2];
sx q[2];
rz(-1.0955398) q[2];
sx q[2];
rz(0.82820669) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6470681) q[1];
sx q[1];
rz(-1.8956603) q[1];
sx q[1];
rz(1.6560133) q[1];
x q[2];
rz(-2.9322846) q[3];
sx q[3];
rz(-0.44692398) q[3];
sx q[3];
rz(-0.79093864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0467726) q[2];
sx q[2];
rz(-0.54868713) q[2];
sx q[2];
rz(-2.0937031) q[2];
rz(1.3607599) q[3];
sx q[3];
rz(-2.4436617) q[3];
sx q[3];
rz(-0.60539436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.027325252) q[0];
sx q[0];
rz(-1.1189168) q[0];
sx q[0];
rz(3.0615321) q[0];
rz(-0.36021532) q[1];
sx q[1];
rz(-1.4604912) q[1];
sx q[1];
rz(2.1866658) q[1];
rz(0.24047492) q[2];
sx q[2];
rz(-2.0839811) q[2];
sx q[2];
rz(2.6032084) q[2];
rz(-0.23056728) q[3];
sx q[3];
rz(-2.1721526) q[3];
sx q[3];
rz(-0.9823907) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
