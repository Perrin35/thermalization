OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3611276) q[0];
sx q[0];
rz(-0.68929231) q[0];
sx q[0];
rz(-0.33049345) q[0];
rz(0.32980546) q[1];
sx q[1];
rz(-0.84996119) q[1];
sx q[1];
rz(0.70911521) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5753484) q[0];
sx q[0];
rz(-2.4465804) q[0];
sx q[0];
rz(-2.038875) q[0];
rz(-pi) q[1];
rz(-1.2572631) q[2];
sx q[2];
rz(-1.6013813) q[2];
sx q[2];
rz(-0.41536301) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.658537) q[1];
sx q[1];
rz(-0.72209789) q[1];
sx q[1];
rz(-2.2621821) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0508918) q[3];
sx q[3];
rz(-1.2399779) q[3];
sx q[3];
rz(-2.1045121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7314529) q[2];
sx q[2];
rz(-0.4099161) q[2];
sx q[2];
rz(1.6072134) q[2];
rz(-2.2065227) q[3];
sx q[3];
rz(-1.3053223) q[3];
sx q[3];
rz(-2.3527761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7222897) q[0];
sx q[0];
rz(-0.062148217) q[0];
sx q[0];
rz(0.62227917) q[0];
rz(-2.9653446) q[1];
sx q[1];
rz(-1.2146249) q[1];
sx q[1];
rz(2.2252749) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5641862) q[0];
sx q[0];
rz(-2.2803377) q[0];
sx q[0];
rz(0.90155154) q[0];
rz(-pi) q[1];
rz(3.0099478) q[2];
sx q[2];
rz(-2.3028214) q[2];
sx q[2];
rz(-1.6694348) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.898145) q[1];
sx q[1];
rz(-1.2444278) q[1];
sx q[1];
rz(-1.2136202) q[1];
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
rz(-0.24094412) q[2];
sx q[2];
rz(-2.1449461) q[2];
sx q[2];
rz(2.3201578) q[2];
rz(0.017283043) q[3];
sx q[3];
rz(-0.79764962) q[3];
sx q[3];
rz(-0.78330529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18773742) q[0];
sx q[0];
rz(-1.563235) q[0];
sx q[0];
rz(-1.0082555) q[0];
rz(-0.035765212) q[1];
sx q[1];
rz(-1.4859896) q[1];
sx q[1];
rz(0.52454138) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6150104) q[0];
sx q[0];
rz(-2.9752762) q[0];
sx q[0];
rz(1.389857) q[0];
rz(-1.8345941) q[2];
sx q[2];
rz(-2.2321777) q[2];
sx q[2];
rz(-1.3493354) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.14568612) q[1];
sx q[1];
rz(-1.6325163) q[1];
sx q[1];
rz(1.5476336) q[1];
rz(-1.9906304) q[3];
sx q[3];
rz(-1.435278) q[3];
sx q[3];
rz(-0.47437048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3699469) q[2];
sx q[2];
rz(-1.5180072) q[2];
sx q[2];
rz(-1.6463722) q[2];
rz(1.8203991) q[3];
sx q[3];
rz(-1.4097872) q[3];
sx q[3];
rz(0.43740073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8111073) q[0];
sx q[0];
rz(-2.2225668) q[0];
sx q[0];
rz(0.088949732) q[0];
rz(-0.51070172) q[1];
sx q[1];
rz(-2.316244) q[1];
sx q[1];
rz(-0.68960062) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2992069) q[0];
sx q[0];
rz(-0.33873522) q[0];
sx q[0];
rz(2.6329106) q[0];
x q[1];
rz(2.7933649) q[2];
sx q[2];
rz(-0.59154445) q[2];
sx q[2];
rz(1.0500184) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5213485) q[1];
sx q[1];
rz(-0.84940956) q[1];
sx q[1];
rz(1.1299302) q[1];
rz(0.27130178) q[3];
sx q[3];
rz(-0.50008431) q[3];
sx q[3];
rz(0.015319583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4758063) q[2];
sx q[2];
rz(-0.67725724) q[2];
sx q[2];
rz(-2.518667) q[2];
rz(2.0056491) q[3];
sx q[3];
rz(-2.9562852) q[3];
sx q[3];
rz(2.6749271) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85161197) q[0];
sx q[0];
rz(-1.8988134) q[0];
sx q[0];
rz(-0.583453) q[0];
rz(1.1460229) q[1];
sx q[1];
rz(-1.6505516) q[1];
sx q[1];
rz(-1.6437644) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0935055) q[0];
sx q[0];
rz(-1.830528) q[0];
sx q[0];
rz(-2.3555059) q[0];
rz(3.0450902) q[2];
sx q[2];
rz(-2.3415903) q[2];
sx q[2];
rz(-0.9347136) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0651107) q[1];
sx q[1];
rz(-1.5981734) q[1];
sx q[1];
rz(-0.70365023) q[1];
x q[2];
rz(-0.28108092) q[3];
sx q[3];
rz(-2.6087458) q[3];
sx q[3];
rz(-1.3827152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4862165) q[2];
sx q[2];
rz(-1.5950173) q[2];
sx q[2];
rz(1.5636469) q[2];
rz(-0.90562138) q[3];
sx q[3];
rz(-0.27035299) q[3];
sx q[3];
rz(-0.63846987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49801302) q[0];
sx q[0];
rz(-1.4587412) q[0];
sx q[0];
rz(-0.046982732) q[0];
rz(-2.9934096) q[1];
sx q[1];
rz(-2.2427185) q[1];
sx q[1];
rz(-1.4354338) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3008266) q[0];
sx q[0];
rz(-2.3593785) q[0];
sx q[0];
rz(-1.7813111) q[0];
rz(-pi) q[1];
rz(0.058850364) q[2];
sx q[2];
rz(-1.0219136) q[2];
sx q[2];
rz(-0.6097874) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.012949) q[1];
sx q[1];
rz(-1.1638906) q[1];
sx q[1];
rz(-0.13601555) q[1];
rz(-pi) q[2];
rz(-1.3330323) q[3];
sx q[3];
rz(-1.5974345) q[3];
sx q[3];
rz(-2.3029072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0662213) q[2];
sx q[2];
rz(-1.0152638) q[2];
sx q[2];
rz(0.071468778) q[2];
rz(1.6890769) q[3];
sx q[3];
rz(-1.444933) q[3];
sx q[3];
rz(-2.215109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8554095) q[0];
sx q[0];
rz(-1.8137285) q[0];
sx q[0];
rz(2.563971) q[0];
rz(-1.8619934) q[1];
sx q[1];
rz(-1.7956693) q[1];
sx q[1];
rz(1.0095899) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9855232) q[0];
sx q[0];
rz(-1.0970322) q[0];
sx q[0];
rz(-1.1958634) q[0];
x q[1];
rz(-0.38452734) q[2];
sx q[2];
rz(-1.782801) q[2];
sx q[2];
rz(2.7146313) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.84570388) q[1];
sx q[1];
rz(-1.0755952) q[1];
sx q[1];
rz(1.9572898) q[1];
x q[2];
rz(0.76241242) q[3];
sx q[3];
rz(-2.0250642) q[3];
sx q[3];
rz(1.4095969) q[3];
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
rz(0.47232929) q[3];
sx q[3];
rz(-1.4940558) q[3];
sx q[3];
rz(-2.1388617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6427479) q[0];
sx q[0];
rz(-1.6413178) q[0];
sx q[0];
rz(0.2391267) q[0];
rz(0.7827951) q[1];
sx q[1];
rz(-2.6233964) q[1];
sx q[1];
rz(0.4447287) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1040161) q[0];
sx q[0];
rz(-0.74066478) q[0];
sx q[0];
rz(-0.91233493) q[0];
x q[1];
rz(-0.66419454) q[2];
sx q[2];
rz(-1.0401298) q[2];
sx q[2];
rz(-1.2150089) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7263033) q[1];
sx q[1];
rz(-2.2518034) q[1];
sx q[1];
rz(-0.16420941) q[1];
rz(-0.41110699) q[3];
sx q[3];
rz(-1.6931705) q[3];
sx q[3];
rz(-1.7468332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7198221) q[2];
sx q[2];
rz(-2.3193216) q[2];
sx q[2];
rz(2.3262809) q[2];
rz(-1.0347962) q[3];
sx q[3];
rz(-1.5765604) q[3];
sx q[3];
rz(-1.3172654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3141044) q[0];
sx q[0];
rz(-1.0359456) q[0];
sx q[0];
rz(-0.28717336) q[0];
rz(2.9526967) q[1];
sx q[1];
rz(-2.4177528) q[1];
sx q[1];
rz(-2.8093991) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0848541) q[0];
sx q[0];
rz(-0.85002725) q[0];
sx q[0];
rz(2.7194354) q[0];
x q[1];
rz(-1.6261149) q[2];
sx q[2];
rz(-2.4419867) q[2];
sx q[2];
rz(2.827022) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5481944) q[1];
sx q[1];
rz(-1.3151004) q[1];
sx q[1];
rz(0.42773186) q[1];
x q[2];
rz(0.53211777) q[3];
sx q[3];
rz(-1.2174774) q[3];
sx q[3];
rz(1.4766828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2087848) q[2];
sx q[2];
rz(-1.00495) q[2];
sx q[2];
rz(-2.3020111) q[2];
rz(1.2906637) q[3];
sx q[3];
rz(-1.3495812) q[3];
sx q[3];
rz(2.4850142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
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
rz(-1.1971808) q[1];
sx q[1];
rz(0.25451452) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8043038) q[0];
sx q[0];
rz(-0.24233195) q[0];
sx q[0];
rz(1.2605577) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.43912402) q[2];
sx q[2];
rz(-2.6115978) q[2];
sx q[2];
rz(2.8119171) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6470681) q[1];
sx q[1];
rz(-1.8956603) q[1];
sx q[1];
rz(1.4855794) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.43838318) q[3];
sx q[3];
rz(-1.66072) q[3];
sx q[3];
rz(-2.5509978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0467726) q[2];
sx q[2];
rz(-2.5929055) q[2];
sx q[2];
rz(1.0478896) q[2];
rz(-1.7808328) q[3];
sx q[3];
rz(-2.4436617) q[3];
sx q[3];
rz(-0.60539436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1142674) q[0];
sx q[0];
rz(-1.1189168) q[0];
sx q[0];
rz(3.0615321) q[0];
rz(0.36021532) q[1];
sx q[1];
rz(-1.6811014) q[1];
sx q[1];
rz(-0.95492687) q[1];
rz(-1.1709471) q[2];
sx q[2];
rz(-0.56213899) q[2];
sx q[2];
rz(3.0664372) q[2];
rz(-0.9568692) q[3];
sx q[3];
rz(-1.7603684) q[3];
sx q[3];
rz(-2.6852222) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];