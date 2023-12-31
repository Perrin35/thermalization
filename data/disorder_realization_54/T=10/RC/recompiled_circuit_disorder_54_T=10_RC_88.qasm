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
rz(-0.84996119) q[1];
sx q[1];
rz(0.70911521) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1256975) q[0];
sx q[0];
rz(-2.1791434) q[0];
sx q[0];
rz(2.7817821) q[0];
rz(-1.4719226) q[2];
sx q[2];
rz(-2.8266202) q[2];
sx q[2];
rz(1.0613943) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.658537) q[1];
sx q[1];
rz(-2.4194948) q[1];
sx q[1];
rz(-0.87941054) q[1];
x q[2];
rz(-2.1756644) q[3];
sx q[3];
rz(-0.60797193) q[3];
sx q[3];
rz(-1.0498429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4101397) q[2];
sx q[2];
rz(-0.4099161) q[2];
sx q[2];
rz(-1.5343792) q[2];
rz(0.93506995) q[3];
sx q[3];
rz(-1.8362703) q[3];
sx q[3];
rz(-0.7888166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4193029) q[0];
sx q[0];
rz(-3.0794444) q[0];
sx q[0];
rz(2.5193135) q[0];
rz(2.9653446) q[1];
sx q[1];
rz(-1.9269678) q[1];
sx q[1];
rz(-0.91631779) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8371671) q[0];
sx q[0];
rz(-0.93351782) q[0];
sx q[0];
rz(0.62563719) q[0];
rz(-pi) q[1];
rz(2.3071438) q[2];
sx q[2];
rz(-1.6685899) q[2];
sx q[2];
rz(-0.18690878) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.037109) q[1];
sx q[1];
rz(-0.47905211) q[1];
sx q[1];
rz(2.3399809) q[1];
rz(-pi) q[2];
rz(-0.022577062) q[3];
sx q[3];
rz(-1.1728247) q[3];
sx q[3];
rz(-0.28385362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.24094412) q[2];
sx q[2];
rz(-0.99664656) q[2];
sx q[2];
rz(-2.3201578) q[2];
rz(-3.1243096) q[3];
sx q[3];
rz(-2.343943) q[3];
sx q[3];
rz(0.78330529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9538552) q[0];
sx q[0];
rz(-1.5783577) q[0];
sx q[0];
rz(-2.1333372) q[0];
rz(-0.035765212) q[1];
sx q[1];
rz(-1.6556031) q[1];
sx q[1];
rz(2.6170513) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7099972) q[0];
sx q[0];
rz(-1.4072197) q[0];
sx q[0];
rz(3.1113935) q[0];
rz(-pi) q[1];
rz(2.4630765) q[2];
sx q[2];
rz(-1.3635474) q[2];
sx q[2];
rz(-3.0845272) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9959065) q[1];
sx q[1];
rz(-1.5090764) q[1];
sx q[1];
rz(1.5476336) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1509622) q[3];
sx q[3];
rz(-1.7063147) q[3];
sx q[3];
rz(-2.6672222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3699469) q[2];
sx q[2];
rz(-1.6235855) q[2];
sx q[2];
rz(1.6463722) q[2];
rz(1.3211936) q[3];
sx q[3];
rz(-1.4097872) q[3];
sx q[3];
rz(2.7041919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.8111073) q[0];
sx q[0];
rz(-2.2225668) q[0];
sx q[0];
rz(3.0526429) q[0];
rz(2.6308909) q[1];
sx q[1];
rz(-2.316244) q[1];
sx q[1];
rz(2.451992) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9289963) q[0];
sx q[0];
rz(-1.7333475) q[0];
sx q[0];
rz(2.8430804) q[0];
rz(-0.56324048) q[2];
sx q[2];
rz(-13/(3*pi)) q[2];
sx q[2];
rz(2.9134977) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2526306) q[1];
sx q[1];
rz(-1.8969715) q[1];
sx q[1];
rz(0.77146448) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8702909) q[3];
sx q[3];
rz(-0.50008431) q[3];
sx q[3];
rz(3.1262731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.66578635) q[2];
sx q[2];
rz(-2.4643354) q[2];
sx q[2];
rz(-0.62292567) q[2];
rz(1.1359435) q[3];
sx q[3];
rz(-2.9562852) q[3];
sx q[3];
rz(0.46666551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85161197) q[0];
sx q[0];
rz(-1.2427793) q[0];
sx q[0];
rz(0.583453) q[0];
rz(1.9955697) q[1];
sx q[1];
rz(-1.6505516) q[1];
sx q[1];
rz(-1.4978283) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0935055) q[0];
sx q[0];
rz(-1.830528) q[0];
sx q[0];
rz(0.78608677) q[0];
rz(-pi) q[1];
rz(-1.4719109) q[2];
sx q[2];
rz(-2.3660198) q[2];
sx q[2];
rz(2.0688187) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6240546) q[1];
sx q[1];
rz(-2.2741286) q[1];
sx q[1];
rz(-1.6066949) q[1];
rz(1.4086401) q[3];
sx q[3];
rz(-2.080653) q[3];
sx q[3];
rz(-1.4354524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.65537611) q[2];
sx q[2];
rz(-1.5950173) q[2];
sx q[2];
rz(-1.5779457) q[2];
rz(-2.2359713) q[3];
sx q[3];
rz(-0.27035299) q[3];
sx q[3];
rz(-2.5031228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49801302) q[0];
sx q[0];
rz(-1.4587412) q[0];
sx q[0];
rz(0.046982732) q[0];
rz(-0.14818305) q[1];
sx q[1];
rz(-2.2427185) q[1];
sx q[1];
rz(-1.7061589) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2611321) q[0];
sx q[0];
rz(-1.4229703) q[0];
sx q[0];
rz(-0.79974215) q[0];
x q[1];
rz(-0.058850364) q[2];
sx q[2];
rz(-1.0219136) q[2];
sx q[2];
rz(0.6097874) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.012949) q[1];
sx q[1];
rz(-1.9777021) q[1];
sx q[1];
rz(-0.13601555) q[1];
x q[2];
rz(-0.027408882) q[3];
sx q[3];
rz(-1.3331183) q[3];
sx q[3];
rz(0.72565597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0753714) q[2];
sx q[2];
rz(-1.0152638) q[2];
sx q[2];
rz(-0.071468778) q[2];
rz(-1.6890769) q[3];
sx q[3];
rz(-1.6966597) q[3];
sx q[3];
rz(0.92648363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8554095) q[0];
sx q[0];
rz(-1.3278642) q[0];
sx q[0];
rz(0.57762161) q[0];
rz(-1.8619934) q[1];
sx q[1];
rz(-1.7956693) q[1];
sx q[1];
rz(-2.1320027) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59238626) q[0];
sx q[0];
rz(-1.9027332) q[0];
sx q[0];
rz(-0.50360002) q[0];
x q[1];
rz(-1.3426443) q[2];
sx q[2];
rz(-1.1953127) q[2];
sx q[2];
rz(2.0827039) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2254667) q[1];
sx q[1];
rz(-1.2327317) q[1];
sx q[1];
rz(0.52789968) q[1];
x q[2];
rz(0.97686751) q[3];
sx q[3];
rz(-0.90126029) q[3];
sx q[3];
rz(-0.23564786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0499095) q[2];
sx q[2];
rz(-0.35353264) q[2];
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
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6427479) q[0];
sx q[0];
rz(-1.6413178) q[0];
sx q[0];
rz(0.2391267) q[0];
rz(2.3587976) q[1];
sx q[1];
rz(-2.6233964) q[1];
sx q[1];
rz(2.696864) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1040161) q[0];
sx q[0];
rz(-2.4009279) q[0];
sx q[0];
rz(-0.91233493) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66419454) q[2];
sx q[2];
rz(-1.0401298) q[2];
sx q[2];
rz(-1.9265837) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6726482) q[1];
sx q[1];
rz(-0.69744195) q[1];
sx q[1];
rz(1.3717321) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4374251) q[3];
sx q[3];
rz(-1.9786454) q[3];
sx q[3];
rz(-2.9123902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7198221) q[2];
sx q[2];
rz(-0.82227102) q[2];
sx q[2];
rz(-0.81531173) q[2];
rz(-2.1067965) q[3];
sx q[3];
rz(-1.5650322) q[3];
sx q[3];
rz(1.8243272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3141044) q[0];
sx q[0];
rz(-1.0359456) q[0];
sx q[0];
rz(0.28717336) q[0];
rz(-0.18889591) q[1];
sx q[1];
rz(-0.72383988) q[1];
sx q[1];
rz(2.8093991) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19776343) q[0];
sx q[0];
rz(-1.257886) q[0];
sx q[0];
rz(2.3373332) q[0];
rz(-pi) q[1];
rz(-1.6261149) q[2];
sx q[2];
rz(-2.4419867) q[2];
sx q[2];
rz(-0.31457065) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5481944) q[1];
sx q[1];
rz(-1.3151004) q[1];
sx q[1];
rz(-2.7138608) q[1];
rz(-pi) q[2];
rz(0.62854564) q[3];
sx q[3];
rz(-0.62918951) q[3];
sx q[3];
rz(0.43720804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2087848) q[2];
sx q[2];
rz(-2.1366426) q[2];
sx q[2];
rz(-2.3020111) q[2];
rz(-1.8509289) q[3];
sx q[3];
rz(-1.3495812) q[3];
sx q[3];
rz(-0.65657842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-3.0861417) q[0];
sx q[0];
rz(-2.3738326) q[0];
sx q[0];
rz(1.0593876) q[0];
rz(2.1620031) q[1];
sx q[1];
rz(-1.1971808) q[1];
sx q[1];
rz(-2.8870781) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33728889) q[0];
sx q[0];
rz(-2.8992607) q[0];
sx q[0];
rz(-1.8810349) q[0];
x q[1];
rz(-0.43912402) q[2];
sx q[2];
rz(-0.52999485) q[2];
sx q[2];
rz(-2.8119171) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6470681) q[1];
sx q[1];
rz(-1.8956603) q[1];
sx q[1];
rz(-1.4855794) q[1];
rz(-pi) q[2];
rz(-1.47154) q[3];
sx q[3];
rz(-2.0072862) q[3];
sx q[3];
rz(2.1193159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0467726) q[2];
sx q[2];
rz(-0.54868713) q[2];
sx q[2];
rz(2.0937031) q[2];
rz(-1.7808328) q[3];
sx q[3];
rz(-0.69793099) q[3];
sx q[3];
rz(0.60539436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-3.1142674) q[0];
sx q[0];
rz(-2.0226759) q[0];
sx q[0];
rz(-0.080060536) q[0];
rz(-2.7813773) q[1];
sx q[1];
rz(-1.6811014) q[1];
sx q[1];
rz(-0.95492687) q[1];
rz(0.24047492) q[2];
sx q[2];
rz(-2.0839811) q[2];
sx q[2];
rz(2.6032084) q[2];
rz(-1.249282) q[3];
sx q[3];
rz(-2.5026863) q[3];
sx q[3];
rz(1.7659059) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
