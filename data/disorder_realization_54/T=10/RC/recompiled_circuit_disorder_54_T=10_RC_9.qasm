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
rz(0.32980546) q[1];
sx q[1];
rz(5.4332241) q[1];
sx q[1];
rz(10.133893) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.015895122) q[0];
sx q[0];
rz(-2.1791434) q[0];
sx q[0];
rz(2.7817821) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.032151392) q[2];
sx q[2];
rz(-1.8841779) q[2];
sx q[2];
rz(-1.9762447) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6735437) q[1];
sx q[1];
rz(-1.1357726) q[1];
sx q[1];
rz(2.166964) q[1];
x q[2];
rz(-2.764774) q[3];
sx q[3];
rz(-2.0599277) q[3];
sx q[3];
rz(-2.7917142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4101397) q[2];
sx q[2];
rz(-0.4099161) q[2];
sx q[2];
rz(-1.6072134) q[2];
rz(2.2065227) q[3];
sx q[3];
rz(-1.3053223) q[3];
sx q[3];
rz(2.3527761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.4193029) q[0];
sx q[0];
rz(-3.0794444) q[0];
sx q[0];
rz(-2.5193135) q[0];
rz(0.17624804) q[1];
sx q[1];
rz(-1.9269678) q[1];
sx q[1];
rz(0.91631779) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.672357) q[0];
sx q[0];
rz(-1.0807481) q[0];
sx q[0];
rz(0.83067466) q[0];
rz(1.7158521) q[2];
sx q[2];
rz(-0.74160355) q[2];
sx q[2];
rz(-1.2765826) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1044836) q[1];
sx q[1];
rz(-0.47905211) q[1];
sx q[1];
rz(2.3399809) q[1];
rz(1.1727337) q[3];
sx q[3];
rz(-1.549984) q[3];
sx q[3];
rz(1.845899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.24094412) q[2];
sx q[2];
rz(-2.1449461) q[2];
sx q[2];
rz(2.3201578) q[2];
rz(-0.017283043) q[3];
sx q[3];
rz(-0.79764962) q[3];
sx q[3];
rz(-2.3582874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9538552) q[0];
sx q[0];
rz(-1.563235) q[0];
sx q[0];
rz(1.0082555) q[0];
rz(-0.035765212) q[1];
sx q[1];
rz(-1.6556031) q[1];
sx q[1];
rz(2.6170513) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4315955) q[0];
sx q[0];
rz(-1.734373) q[0];
sx q[0];
rz(0.030199108) q[0];
x q[1];
rz(-0.67851615) q[2];
sx q[2];
rz(-1.3635474) q[2];
sx q[2];
rz(0.057065406) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7150535) q[1];
sx q[1];
rz(-1.5476777) q[1];
sx q[1];
rz(-0.061736488) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1509622) q[3];
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
rz(0.77164578) q[2];
sx q[2];
rz(-1.5180072) q[2];
sx q[2];
rz(-1.6463722) q[2];
rz(1.3211936) q[3];
sx q[3];
rz(-1.7318055) q[3];
sx q[3];
rz(-2.7041919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-0.33048531) q[0];
sx q[0];
rz(-0.91902584) q[0];
sx q[0];
rz(-3.0526429) q[0];
rz(0.51070172) q[1];
sx q[1];
rz(-0.82534868) q[1];
sx q[1];
rz(2.451992) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21259637) q[0];
sx q[0];
rz(-1.7333475) q[0];
sx q[0];
rz(0.29851229) q[0];
rz(-pi) q[1];
rz(-0.34822779) q[2];
sx q[2];
rz(-2.5500482) q[2];
sx q[2];
rz(-1.0500184) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1414813) q[1];
sx q[1];
rz(-0.82427102) q[1];
sx q[1];
rz(0.4517171) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7161937) q[3];
sx q[3];
rz(-1.0905915) q[3];
sx q[3];
rz(2.8499545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.66578635) q[2];
sx q[2];
rz(-2.4643354) q[2];
sx q[2];
rz(-0.62292567) q[2];
rz(2.0056491) q[3];
sx q[3];
rz(-2.9562852) q[3];
sx q[3];
rz(2.6749271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2899807) q[0];
sx q[0];
rz(-1.2427793) q[0];
sx q[0];
rz(0.583453) q[0];
rz(-1.9955697) q[1];
sx q[1];
rz(-1.6505516) q[1];
sx q[1];
rz(-1.6437644) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4125815) q[0];
sx q[0];
rz(-2.3238365) q[0];
sx q[0];
rz(-1.2110932) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0450902) q[2];
sx q[2];
rz(-2.3415903) q[2];
sx q[2];
rz(2.2068791) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0651107) q[1];
sx q[1];
rz(-1.5981734) q[1];
sx q[1];
rz(-2.4379424) q[1];
x q[2];
rz(-1.7329526) q[3];
sx q[3];
rz(-2.080653) q[3];
sx q[3];
rz(-1.4354524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4862165) q[2];
sx q[2];
rz(-1.5465753) q[2];
sx q[2];
rz(1.5636469) q[2];
rz(0.90562138) q[3];
sx q[3];
rz(-0.27035299) q[3];
sx q[3];
rz(0.63846987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6435796) q[0];
sx q[0];
rz(-1.4587412) q[0];
sx q[0];
rz(-3.0946099) q[0];
rz(0.14818305) q[1];
sx q[1];
rz(-2.2427185) q[1];
sx q[1];
rz(1.7061589) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88046056) q[0];
sx q[0];
rz(-1.7186223) q[0];
sx q[0];
rz(-0.79974215) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4749182) q[2];
sx q[2];
rz(-0.55170689) q[2];
sx q[2];
rz(-2.4193537) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1286436) q[1];
sx q[1];
rz(-1.1638906) q[1];
sx q[1];
rz(3.0055771) q[1];
x q[2];
rz(1.8085603) q[3];
sx q[3];
rz(-1.5974345) q[3];
sx q[3];
rz(-2.3029072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0753714) q[2];
sx q[2];
rz(-1.0152638) q[2];
sx q[2];
rz(-3.0701239) q[2];
rz(-1.4525157) q[3];
sx q[3];
rz(-1.444933) q[3];
sx q[3];
rz(0.92648363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(1.0095899) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1560695) q[0];
sx q[0];
rz(-1.0970322) q[0];
sx q[0];
rz(-1.9457293) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6206714) q[2];
sx q[2];
rz(-0.43653566) q[2];
sx q[2];
rz(1.6233363) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.84570388) q[1];
sx q[1];
rz(-2.0659975) q[1];
sx q[1];
rz(-1.9572898) q[1];
rz(2.3791802) q[3];
sx q[3];
rz(-1.1165285) q[3];
sx q[3];
rz(-1.7319958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0499095) q[2];
sx q[2];
rz(-2.78806) q[2];
sx q[2];
rz(-1.1996777) q[2];
rz(-2.6692634) q[3];
sx q[3];
rz(-1.4940558) q[3];
sx q[3];
rz(-2.1388617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6427479) q[0];
sx q[0];
rz(-1.6413178) q[0];
sx q[0];
rz(-2.902466) q[0];
rz(-0.7827951) q[1];
sx q[1];
rz(-0.51819623) q[1];
sx q[1];
rz(0.4447287) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1040161) q[0];
sx q[0];
rz(-2.4009279) q[0];
sx q[0];
rz(-2.2292577) q[0];
x q[1];
rz(-2.2112446) q[2];
sx q[2];
rz(-1.010251) q[2];
sx q[2];
rz(2.408574) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7263033) q[1];
sx q[1];
rz(-0.88978926) q[1];
sx q[1];
rz(-0.16420941) q[1];
rz(-2.8430311) q[3];
sx q[3];
rz(-2.7136554) q[3];
sx q[3];
rz(-0.096979389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7198221) q[2];
sx q[2];
rz(-2.3193216) q[2];
sx q[2];
rz(2.3262809) q[2];
rz(-2.1067965) q[3];
sx q[3];
rz(-1.5650322) q[3];
sx q[3];
rz(1.8243272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3141044) q[0];
sx q[0];
rz(-1.0359456) q[0];
sx q[0];
rz(0.28717336) q[0];
rz(0.18889591) q[1];
sx q[1];
rz(-2.4177528) q[1];
sx q[1];
rz(-0.33219355) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6824326) q[0];
sx q[0];
rz(-0.81572616) q[0];
sx q[0];
rz(2.0072323) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2696482) q[2];
sx q[2];
rz(-1.6064062) q[2];
sx q[2];
rz(-1.2138838) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5481944) q[1];
sx q[1];
rz(-1.3151004) q[1];
sx q[1];
rz(0.42773186) q[1];
rz(-pi) q[2];
rz(-0.53211777) q[3];
sx q[3];
rz(-1.2174774) q[3];
sx q[3];
rz(1.6649099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.93280783) q[2];
sx q[2];
rz(-1.00495) q[2];
sx q[2];
rz(-0.83958158) q[2];
rz(-1.2906637) q[3];
sx q[3];
rz(-1.7920114) q[3];
sx q[3];
rz(-0.65657842) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.055450913) q[0];
sx q[0];
rz(-0.76776004) q[0];
sx q[0];
rz(1.0593876) q[0];
rz(0.97958952) q[1];
sx q[1];
rz(-1.9444119) q[1];
sx q[1];
rz(0.25451452) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1232676) q[0];
sx q[0];
rz(-1.3402481) q[0];
sx q[0];
rz(-0.075320764) q[0];
x q[1];
rz(-1.326667) q[2];
sx q[2];
rz(-1.0955398) q[2];
sx q[2];
rz(-2.313386) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6470681) q[1];
sx q[1];
rz(-1.2459323) q[1];
sx q[1];
rz(-1.6560133) q[1];
rz(2.7032095) q[3];
sx q[3];
rz(-1.66072) q[3];
sx q[3];
rz(-2.5509978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.09482) q[2];
sx q[2];
rz(-0.54868713) q[2];
sx q[2];
rz(-1.0478896) q[2];
rz(1.3607599) q[3];
sx q[3];
rz(-0.69793099) q[3];
sx q[3];
rz(-2.5361983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.027325252) q[0];
sx q[0];
rz(-2.0226759) q[0];
sx q[0];
rz(-0.080060536) q[0];
rz(2.7813773) q[1];
sx q[1];
rz(-1.4604912) q[1];
sx q[1];
rz(2.1866658) q[1];
rz(0.24047492) q[2];
sx q[2];
rz(-2.0839811) q[2];
sx q[2];
rz(2.6032084) q[2];
rz(2.1847235) q[3];
sx q[3];
rz(-1.7603684) q[3];
sx q[3];
rz(-2.6852222) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
