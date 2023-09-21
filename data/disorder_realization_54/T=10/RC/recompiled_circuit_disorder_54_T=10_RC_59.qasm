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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.015895122) q[0];
sx q[0];
rz(-0.96244922) q[0];
sx q[0];
rz(-2.7817821) q[0];
x q[1];
rz(1.66967) q[2];
sx q[2];
rz(-2.8266202) q[2];
sx q[2];
rz(1.0613943) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6735437) q[1];
sx q[1];
rz(-1.1357726) q[1];
sx q[1];
rz(-0.97462868) q[1];
rz(-1.0508918) q[3];
sx q[3];
rz(-1.9016148) q[3];
sx q[3];
rz(-2.1045121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7314529) q[2];
sx q[2];
rz(-2.7316766) q[2];
sx q[2];
rz(1.6072134) q[2];
rz(0.93506995) q[3];
sx q[3];
rz(-1.3053223) q[3];
sx q[3];
rz(-2.3527761) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4193029) q[0];
sx q[0];
rz(-0.062148217) q[0];
sx q[0];
rz(-2.5193135) q[0];
rz(-0.17624804) q[1];
sx q[1];
rz(-1.2146249) q[1];
sx q[1];
rz(-2.2252749) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57740649) q[0];
sx q[0];
rz(-0.86125492) q[0];
sx q[0];
rz(-2.2400411) q[0];
x q[1];
rz(-0.83444886) q[2];
sx q[2];
rz(-1.6685899) q[2];
sx q[2];
rz(-0.18690878) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9333222) q[1];
sx q[1];
rz(-1.2332488) q[1];
sx q[1];
rz(2.7949105) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1727337) q[3];
sx q[3];
rz(-1.5916087) q[3];
sx q[3];
rz(1.845899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.24094412) q[2];
sx q[2];
rz(-0.99664656) q[2];
sx q[2];
rz(0.82143482) q[2];
rz(-0.017283043) q[3];
sx q[3];
rz(-0.79764962) q[3];
sx q[3];
rz(-2.3582874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18773742) q[0];
sx q[0];
rz(-1.5783577) q[0];
sx q[0];
rz(-2.1333372) q[0];
rz(0.035765212) q[1];
sx q[1];
rz(-1.4859896) q[1];
sx q[1];
rz(-0.52454138) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0073111) q[0];
sx q[0];
rz(-1.6005922) q[0];
sx q[0];
rz(1.4071464) q[0];
x q[1];
rz(0.67851615) q[2];
sx q[2];
rz(-1.3635474) q[2];
sx q[2];
rz(3.0845272) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.14568612) q[1];
sx q[1];
rz(-1.5090764) q[1];
sx q[1];
rz(-1.593959) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9933661) q[3];
sx q[3];
rz(-1.1550511) q[3];
sx q[3];
rz(-1.1566597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.77164578) q[2];
sx q[2];
rz(-1.6235855) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
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
rz(-0.82534868) q[1];
sx q[1];
rz(0.68960062) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21259637) q[0];
sx q[0];
rz(-1.7333475) q[0];
sx q[0];
rz(0.29851229) q[0];
x q[1];
rz(2.7933649) q[2];
sx q[2];
rz(-2.5500482) q[2];
sx q[2];
rz(2.0915742) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.888962) q[1];
sx q[1];
rz(-1.2446212) q[1];
sx q[1];
rz(0.77146448) q[1];
x q[2];
rz(2.6570286) q[3];
sx q[3];
rz(-1.441941) q[3];
sx q[3];
rz(1.7948922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.66578635) q[2];
sx q[2];
rz(-2.4643354) q[2];
sx q[2];
rz(0.62292567) q[2];
rz(-2.0056491) q[3];
sx q[3];
rz(-0.1853075) q[3];
sx q[3];
rz(-0.46666551) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2899807) q[0];
sx q[0];
rz(-1.2427793) q[0];
sx q[0];
rz(-2.5581397) q[0];
rz(-1.1460229) q[1];
sx q[1];
rz(-1.4910411) q[1];
sx q[1];
rz(1.4978283) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0935055) q[0];
sx q[0];
rz(-1.830528) q[0];
sx q[0];
rz(-0.78608677) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3439212) q[2];
sx q[2];
rz(-1.5016218) q[2];
sx q[2];
rz(0.56874146) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.51753804) q[1];
sx q[1];
rz(-2.2741286) q[1];
sx q[1];
rz(-1.6066949) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7329526) q[3];
sx q[3];
rz(-2.080653) q[3];
sx q[3];
rz(1.7061403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4862165) q[2];
sx q[2];
rz(-1.5950173) q[2];
sx q[2];
rz(1.5779457) q[2];
rz(0.90562138) q[3];
sx q[3];
rz(-0.27035299) q[3];
sx q[3];
rz(-2.5031228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49801302) q[0];
sx q[0];
rz(-1.6828515) q[0];
sx q[0];
rz(-0.046982732) q[0];
rz(2.9934096) q[1];
sx q[1];
rz(-0.89887416) q[1];
sx q[1];
rz(-1.4354338) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5481789) q[0];
sx q[0];
rz(-0.81028623) q[0];
sx q[0];
rz(-0.20472783) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6666744) q[2];
sx q[2];
rz(-0.55170689) q[2];
sx q[2];
rz(0.72223896) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1286436) q[1];
sx q[1];
rz(-1.9777021) q[1];
sx q[1];
rz(-3.0055771) q[1];
rz(-3.1141838) q[3];
sx q[3];
rz(-1.8084744) q[3];
sx q[3];
rz(0.72565597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0753714) q[2];
sx q[2];
rz(-1.0152638) q[2];
sx q[2];
rz(0.071468778) q[2];
rz(-1.6890769) q[3];
sx q[3];
rz(-1.6966597) q[3];
sx q[3];
rz(0.92648363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8554095) q[0];
sx q[0];
rz(-1.3278642) q[0];
sx q[0];
rz(-2.563971) q[0];
rz(-1.8619934) q[1];
sx q[1];
rz(-1.7956693) q[1];
sx q[1];
rz(-2.1320027) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1560695) q[0];
sx q[0];
rz(-1.0970322) q[0];
sx q[0];
rz(-1.9457293) q[0];
x q[1];
rz(-0.52092123) q[2];
sx q[2];
rz(-2.705057) q[2];
sx q[2];
rz(1.6233363) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2254667) q[1];
sx q[1];
rz(-1.2327317) q[1];
sx q[1];
rz(2.613693) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5261577) q[3];
sx q[3];
rz(-0.86343599) q[3];
sx q[3];
rz(-2.5497041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0499095) q[2];
sx q[2];
rz(-2.78806) q[2];
sx q[2];
rz(1.9419149) q[2];
rz(-0.47232929) q[3];
sx q[3];
rz(-1.6475369) q[3];
sx q[3];
rz(1.002731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49884477) q[0];
sx q[0];
rz(-1.5002748) q[0];
sx q[0];
rz(2.902466) q[0];
rz(0.7827951) q[1];
sx q[1];
rz(-2.6233964) q[1];
sx q[1];
rz(-2.696864) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77138222) q[0];
sx q[0];
rz(-1.0078197) q[0];
sx q[0];
rz(0.51008205) q[0];
rz(-2.3808001) q[2];
sx q[2];
rz(-0.82423254) q[2];
sx q[2];
rz(0.21812083) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.051561292) q[1];
sx q[1];
rz(-1.6981484) q[1];
sx q[1];
rz(2.2584372) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8430311) q[3];
sx q[3];
rz(-0.42793722) q[3];
sx q[3];
rz(-3.0446133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7198221) q[2];
sx q[2];
rz(-0.82227102) q[2];
sx q[2];
rz(-2.3262809) q[2];
rz(2.1067965) q[3];
sx q[3];
rz(-1.5765604) q[3];
sx q[3];
rz(1.8243272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3141044) q[0];
sx q[0];
rz(-2.1056471) q[0];
sx q[0];
rz(0.28717336) q[0];
rz(-0.18889591) q[1];
sx q[1];
rz(-0.72383988) q[1];
sx q[1];
rz(-0.33219355) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0848541) q[0];
sx q[0];
rz(-0.85002725) q[0];
sx q[0];
rz(-0.42215729) q[0];
x q[1];
rz(-3.095093) q[2];
sx q[2];
rz(-2.2691155) q[2];
sx q[2];
rz(2.7547714) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5933983) q[1];
sx q[1];
rz(-1.8264923) q[1];
sx q[1];
rz(0.42773186) q[1];
rz(-pi) q[2];
rz(-0.53211777) q[3];
sx q[3];
rz(-1.9241153) q[3];
sx q[3];
rz(-1.6649099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.93280783) q[2];
sx q[2];
rz(-2.1366426) q[2];
sx q[2];
rz(0.83958158) q[2];
rz(1.2906637) q[3];
sx q[3];
rz(-1.7920114) q[3];
sx q[3];
rz(-2.4850142) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0861417) q[0];
sx q[0];
rz(-0.76776004) q[0];
sx q[0];
rz(1.0593876) q[0];
rz(-0.97958952) q[1];
sx q[1];
rz(-1.9444119) q[1];
sx q[1];
rz(2.8870781) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5352288) q[0];
sx q[0];
rz(-1.6441206) q[0];
sx q[0];
rz(-1.3396157) q[0];
rz(0.48760957) q[2];
sx q[2];
rz(-1.7874103) q[2];
sx q[2];
rz(-0.85607869) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0925797) q[1];
sx q[1];
rz(-1.6515459) q[1];
sx q[1];
rz(2.8156274) q[1];
rz(-2.9322846) q[3];
sx q[3];
rz(-0.44692398) q[3];
sx q[3];
rz(2.350654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.09482) q[2];
sx q[2];
rz(-0.54868713) q[2];
sx q[2];
rz(2.0937031) q[2];
rz(1.3607599) q[3];
sx q[3];
rz(-0.69793099) q[3];
sx q[3];
rz(-2.5361983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.027325252) q[0];
sx q[0];
rz(-1.1189168) q[0];
sx q[0];
rz(3.0615321) q[0];
rz(0.36021532) q[1];
sx q[1];
rz(-1.6811014) q[1];
sx q[1];
rz(-0.95492687) q[1];
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