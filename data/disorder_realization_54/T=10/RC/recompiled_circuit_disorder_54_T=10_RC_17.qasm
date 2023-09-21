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
rz(2.4324774) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1256975) q[0];
sx q[0];
rz(-2.1791434) q[0];
sx q[0];
rz(-2.7817821) q[0];
x q[1];
rz(-1.8843295) q[2];
sx q[2];
rz(-1.6013813) q[2];
sx q[2];
rz(-2.7262296) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6735437) q[1];
sx q[1];
rz(-2.00582) q[1];
sx q[1];
rz(-2.166964) q[1];
rz(1.0508918) q[3];
sx q[3];
rz(-1.9016148) q[3];
sx q[3];
rz(2.1045121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7314529) q[2];
sx q[2];
rz(-2.7316766) q[2];
sx q[2];
rz(1.6072134) q[2];
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
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4193029) q[0];
sx q[0];
rz(-3.0794444) q[0];
sx q[0];
rz(2.5193135) q[0];
rz(-2.9653446) q[1];
sx q[1];
rz(-1.9269678) q[1];
sx q[1];
rz(0.91631779) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5641862) q[0];
sx q[0];
rz(-0.86125492) q[0];
sx q[0];
rz(-0.90155154) q[0];
rz(-1.7158521) q[2];
sx q[2];
rz(-0.74160355) q[2];
sx q[2];
rz(-1.8650101) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.037109) q[1];
sx q[1];
rz(-2.6625405) q[1];
sx q[1];
rz(2.3399809) q[1];
rz(-pi) q[2];
rz(-1.968859) q[3];
sx q[3];
rz(-1.549984) q[3];
sx q[3];
rz(-1.2956937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9006485) q[2];
sx q[2];
rz(-2.1449461) q[2];
sx q[2];
rz(-2.3201578) q[2];
rz(-0.017283043) q[3];
sx q[3];
rz(-2.343943) q[3];
sx q[3];
rz(-0.78330529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9538552) q[0];
sx q[0];
rz(-1.5783577) q[0];
sx q[0];
rz(2.1333372) q[0];
rz(-3.1058274) q[1];
sx q[1];
rz(-1.4859896) q[1];
sx q[1];
rz(2.6170513) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5265822) q[0];
sx q[0];
rz(-2.9752762) q[0];
sx q[0];
rz(1.7517356) q[0];
rz(-pi) q[1];
rz(1.8345941) q[2];
sx q[2];
rz(-2.2321777) q[2];
sx q[2];
rz(-1.7922572) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4265392) q[1];
sx q[1];
rz(-1.5476777) q[1];
sx q[1];
rz(-0.061736488) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
rz(-2.3699469) q[2];
sx q[2];
rz(-1.5180072) q[2];
sx q[2];
rz(1.4952205) q[2];
rz(-1.8203991) q[3];
sx q[3];
rz(-1.7318055) q[3];
sx q[3];
rz(-2.7041919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8111073) q[0];
sx q[0];
rz(-2.2225668) q[0];
sx q[0];
rz(-0.088949732) q[0];
rz(0.51070172) q[1];
sx q[1];
rz(-2.316244) q[1];
sx q[1];
rz(-2.451992) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8423858) q[0];
sx q[0];
rz(-0.33873522) q[0];
sx q[0];
rz(2.6329106) q[0];
x q[1];
rz(-2.7933649) q[2];
sx q[2];
rz(-2.5500482) q[2];
sx q[2];
rz(-2.0915742) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1414813) q[1];
sx q[1];
rz(-2.3173216) q[1];
sx q[1];
rz(-2.6898756) q[1];
rz(2.6570286) q[3];
sx q[3];
rz(-1.441941) q[3];
sx q[3];
rz(1.7948922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4758063) q[2];
sx q[2];
rz(-0.67725724) q[2];
sx q[2];
rz(-0.62292567) q[2];
rz(-2.0056491) q[3];
sx q[3];
rz(-0.1853075) q[3];
sx q[3];
rz(-0.46666551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85161197) q[0];
sx q[0];
rz(-1.2427793) q[0];
sx q[0];
rz(-0.583453) q[0];
rz(-1.9955697) q[1];
sx q[1];
rz(-1.4910411) q[1];
sx q[1];
rz(1.6437644) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0480872) q[0];
sx q[0];
rz(-1.3110647) q[0];
sx q[0];
rz(0.78608677) q[0];
rz(-pi) q[1];
x q[1];
rz(0.79767144) q[2];
sx q[2];
rz(-1.6399709) q[2];
sx q[2];
rz(-0.56874146) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.46206611) q[1];
sx q[1];
rz(-0.70409173) q[1];
sx q[1];
rz(3.0992938) q[1];
rz(2.8605117) q[3];
sx q[3];
rz(-2.6087458) q[3];
sx q[3];
rz(-1.3827152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4862165) q[2];
sx q[2];
rz(-1.5950173) q[2];
sx q[2];
rz(-1.5779457) q[2];
rz(-0.90562138) q[3];
sx q[3];
rz(-0.27035299) q[3];
sx q[3];
rz(-0.63846987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49801302) q[0];
sx q[0];
rz(-1.4587412) q[0];
sx q[0];
rz(-3.0946099) q[0];
rz(-0.14818305) q[1];
sx q[1];
rz(-2.2427185) q[1];
sx q[1];
rz(1.4354338) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2611321) q[0];
sx q[0];
rz(-1.7186223) q[0];
sx q[0];
rz(0.79974215) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0211421) q[2];
sx q[2];
rz(-1.5205985) q[2];
sx q[2];
rz(-0.93027885) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.50373494) q[1];
sx q[1];
rz(-1.6956455) q[1];
sx q[1];
rz(1.9810852) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8085603) q[3];
sx q[3];
rz(-1.5974345) q[3];
sx q[3];
rz(2.3029072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0662213) q[2];
sx q[2];
rz(-1.0152638) q[2];
sx q[2];
rz(3.0701239) q[2];
rz(1.4525157) q[3];
sx q[3];
rz(-1.444933) q[3];
sx q[3];
rz(2.215109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8554095) q[0];
sx q[0];
rz(-1.8137285) q[0];
sx q[0];
rz(-0.57762161) q[0];
rz(-1.2795992) q[1];
sx q[1];
rz(-1.7956693) q[1];
sx q[1];
rz(2.1320027) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5492064) q[0];
sx q[0];
rz(-1.2388595) q[0];
sx q[0];
rz(0.50360002) q[0];
x q[1];
rz(1.7989484) q[2];
sx q[2];
rz(-1.9462799) q[2];
sx q[2];
rz(1.0588888) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.84570388) q[1];
sx q[1];
rz(-1.0755952) q[1];
sx q[1];
rz(1.1843029) q[1];
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
rz(-2.0916831) q[2];
sx q[2];
rz(-0.35353264) q[2];
sx q[2];
rz(1.9419149) q[2];
rz(2.6692634) q[3];
sx q[3];
rz(-1.4940558) q[3];
sx q[3];
rz(-1.002731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6427479) q[0];
sx q[0];
rz(-1.5002748) q[0];
sx q[0];
rz(0.2391267) q[0];
rz(-2.3587976) q[1];
sx q[1];
rz(-0.51819623) q[1];
sx q[1];
rz(2.696864) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0520154) q[0];
sx q[0];
rz(-1.9964295) q[0];
sx q[0];
rz(2.1968967) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66419454) q[2];
sx q[2];
rz(-1.0401298) q[2];
sx q[2];
rz(1.2150089) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.051561292) q[1];
sx q[1];
rz(-1.4434442) q[1];
sx q[1];
rz(-2.2584372) q[1];
rz(2.8430311) q[3];
sx q[3];
rz(-0.42793722) q[3];
sx q[3];
rz(3.0446133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.42177054) q[2];
sx q[2];
rz(-0.82227102) q[2];
sx q[2];
rz(-2.3262809) q[2];
rz(-2.1067965) q[3];
sx q[3];
rz(-1.5650322) q[3];
sx q[3];
rz(1.8243272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3141044) q[0];
sx q[0];
rz(-1.0359456) q[0];
sx q[0];
rz(-0.28717336) q[0];
rz(-0.18889591) q[1];
sx q[1];
rz(-2.4177528) q[1];
sx q[1];
rz(-2.8093991) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0848541) q[0];
sx q[0];
rz(-0.85002725) q[0];
sx q[0];
rz(-0.42215729) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5154778) q[2];
sx q[2];
rz(-0.69960591) q[2];
sx q[2];
rz(0.31457065) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5481944) q[1];
sx q[1];
rz(-1.3151004) q[1];
sx q[1];
rz(2.7138608) q[1];
rz(-pi) q[2];
rz(-2.6094749) q[3];
sx q[3];
rz(-1.9241153) q[3];
sx q[3];
rz(1.6649099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2087848) q[2];
sx q[2];
rz(-1.00495) q[2];
sx q[2];
rz(2.3020111) q[2];
rz(1.8509289) q[3];
sx q[3];
rz(-1.3495812) q[3];
sx q[3];
rz(0.65657842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(3.0861417) q[0];
sx q[0];
rz(-0.76776004) q[0];
sx q[0];
rz(1.0593876) q[0];
rz(0.97958952) q[1];
sx q[1];
rz(-1.9444119) q[1];
sx q[1];
rz(-2.8870781) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6063639) q[0];
sx q[0];
rz(-1.6441206) q[0];
sx q[0];
rz(-1.801977) q[0];
rz(-pi) q[1];
rz(2.6539831) q[2];
sx q[2];
rz(-1.3541823) q[2];
sx q[2];
rz(-0.85607869) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0490129) q[1];
sx q[1];
rz(-1.6515459) q[1];
sx q[1];
rz(0.32596522) q[1];
rz(1.6700527) q[3];
sx q[3];
rz(-2.0072862) q[3];
sx q[3];
rz(-1.0222767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.09482) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1142674) q[0];
sx q[0];
rz(-1.1189168) q[0];
sx q[0];
rz(3.0615321) q[0];
rz(-0.36021532) q[1];
sx q[1];
rz(-1.4604912) q[1];
sx q[1];
rz(2.1866658) q[1];
rz(2.0965626) q[2];
sx q[2];
rz(-1.3617931) q[2];
sx q[2];
rz(1.1522273) q[2];
rz(0.23056728) q[3];
sx q[3];
rz(-0.96944001) q[3];
sx q[3];
rz(2.159202) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];