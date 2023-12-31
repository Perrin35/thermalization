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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7666727) q[0];
sx q[0];
rz(-1.2776889) q[0];
sx q[0];
rz(-2.2105182) q[0];
rz(1.4719226) q[2];
sx q[2];
rz(-2.8266202) q[2];
sx q[2];
rz(2.0801983) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3173649) q[1];
sx q[1];
rz(-1.0365651) q[1];
sx q[1];
rz(2.6298916) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1756644) q[3];
sx q[3];
rz(-2.5336207) q[3];
sx q[3];
rz(-1.0498429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4101397) q[2];
sx q[2];
rz(-0.4099161) q[2];
sx q[2];
rz(-1.6072134) q[2];
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
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7222897) q[0];
sx q[0];
rz(-0.062148217) q[0];
sx q[0];
rz(-2.5193135) q[0];
rz(-2.9653446) q[1];
sx q[1];
rz(-1.9269678) q[1];
sx q[1];
rz(0.91631779) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30442552) q[0];
sx q[0];
rz(-2.2080748) q[0];
sx q[0];
rz(2.5159555) q[0];
rz(-0.13164481) q[2];
sx q[2];
rz(-2.3028214) q[2];
sx q[2];
rz(-1.6694348) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1044836) q[1];
sx q[1];
rz(-0.47905211) q[1];
sx q[1];
rz(0.80161174) q[1];
x q[2];
rz(1.5171492) q[3];
sx q[3];
rz(-0.39857736) q[3];
sx q[3];
rz(0.22565354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9006485) q[2];
sx q[2];
rz(-2.1449461) q[2];
sx q[2];
rz(2.3201578) q[2];
rz(3.1243096) q[3];
sx q[3];
rz(-0.79764962) q[3];
sx q[3];
rz(0.78330529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9538552) q[0];
sx q[0];
rz(-1.563235) q[0];
sx q[0];
rz(2.1333372) q[0];
rz(3.1058274) q[1];
sx q[1];
rz(-1.6556031) q[1];
sx q[1];
rz(-0.52454138) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6150104) q[0];
sx q[0];
rz(-2.9752762) q[0];
sx q[0];
rz(1.7517356) q[0];
rz(-pi) q[1];
rz(1.3069986) q[2];
sx q[2];
rz(-0.909415) q[2];
sx q[2];
rz(1.3493354) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4265392) q[1];
sx q[1];
rz(-1.5939149) q[1];
sx q[1];
rz(0.061736488) q[1];
rz(-pi) q[2];
rz(-1.9906304) q[3];
sx q[3];
rz(-1.7063147) q[3];
sx q[3];
rz(-2.6672222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3699469) q[2];
sx q[2];
rz(-1.6235855) q[2];
sx q[2];
rz(-1.4952205) q[2];
rz(-1.8203991) q[3];
sx q[3];
rz(-1.7318055) q[3];
sx q[3];
rz(0.43740073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2992069) q[0];
sx q[0];
rz(-2.8028574) q[0];
sx q[0];
rz(-2.6329106) q[0];
x q[1];
rz(1.3454516) q[2];
sx q[2];
rz(-1.0190522) q[2];
sx q[2];
rz(-1.4622886) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.888962) q[1];
sx q[1];
rz(-1.8969715) q[1];
sx q[1];
rz(-0.77146448) q[1];
rz(0.27130178) q[3];
sx q[3];
rz(-0.50008431) q[3];
sx q[3];
rz(0.015319583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2899807) q[0];
sx q[0];
rz(-1.8988134) q[0];
sx q[0];
rz(0.583453) q[0];
rz(-1.1460229) q[1];
sx q[1];
rz(-1.6505516) q[1];
sx q[1];
rz(1.6437644) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22623423) q[0];
sx q[0];
rz(-2.3225473) q[0];
sx q[0];
rz(0.35924964) q[0];
rz(0.79767144) q[2];
sx q[2];
rz(-1.5016218) q[2];
sx q[2];
rz(-2.5728512) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6240546) q[1];
sx q[1];
rz(-0.86746403) q[1];
sx q[1];
rz(1.6066949) q[1];
rz(0.28108092) q[3];
sx q[3];
rz(-0.53284684) q[3];
sx q[3];
rz(-1.3827152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.65537611) q[2];
sx q[2];
rz(-1.5465753) q[2];
sx q[2];
rz(-1.5779457) q[2];
rz(2.2359713) q[3];
sx q[3];
rz(-2.8712397) q[3];
sx q[3];
rz(0.63846987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6435796) q[0];
sx q[0];
rz(-1.6828515) q[0];
sx q[0];
rz(-3.0946099) q[0];
rz(2.9934096) q[1];
sx q[1];
rz(-0.89887416) q[1];
sx q[1];
rz(-1.4354338) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2611321) q[0];
sx q[0];
rz(-1.7186223) q[0];
sx q[0];
rz(2.3418505) q[0];
rz(-pi) q[1];
rz(2.1204505) q[2];
sx q[2];
rz(-1.6209941) q[2];
sx q[2];
rz(0.93027885) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3458851) q[1];
sx q[1];
rz(-2.7137623) q[1];
sx q[1];
rz(1.2659628) q[1];
rz(-pi) q[2];
rz(-1.8085603) q[3];
sx q[3];
rz(-1.5441582) q[3];
sx q[3];
rz(0.8386855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0753714) q[2];
sx q[2];
rz(-2.1263289) q[2];
sx q[2];
rz(3.0701239) q[2];
rz(1.6890769) q[3];
sx q[3];
rz(-1.444933) q[3];
sx q[3];
rz(0.92648363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28618318) q[0];
sx q[0];
rz(-1.8137285) q[0];
sx q[0];
rz(-2.563971) q[0];
rz(-1.8619934) q[1];
sx q[1];
rz(-1.7956693) q[1];
sx q[1];
rz(-2.1320027) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5492064) q[0];
sx q[0];
rz(-1.2388595) q[0];
sx q[0];
rz(-2.6379926) q[0];
rz(-pi) q[1];
rz(0.38452734) q[2];
sx q[2];
rz(-1.3587917) q[2];
sx q[2];
rz(2.7146313) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0040782) q[1];
sx q[1];
rz(-0.61811781) q[1];
sx q[1];
rz(0.60933463) q[1];
rz(-pi) q[2];
x q[2];
rz(0.97686751) q[3];
sx q[3];
rz(-0.90126029) q[3];
sx q[3];
rz(-0.23564786) q[3];
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
rz(1.9419149) q[2];
rz(2.6692634) q[3];
sx q[3];
rz(-1.4940558) q[3];
sx q[3];
rz(2.1388617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-0.49884477) q[0];
sx q[0];
rz(-1.5002748) q[0];
sx q[0];
rz(-2.902466) q[0];
rz(2.3587976) q[1];
sx q[1];
rz(-0.51819623) q[1];
sx q[1];
rz(-2.696864) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3702104) q[0];
sx q[0];
rz(-2.133773) q[0];
sx q[0];
rz(-2.6315106) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2112446) q[2];
sx q[2];
rz(-1.010251) q[2];
sx q[2];
rz(0.7330187) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0900314) q[1];
sx q[1];
rz(-1.6981484) q[1];
sx q[1];
rz(-0.8831555) q[1];
x q[2];
rz(-1.7041676) q[3];
sx q[3];
rz(-1.9786454) q[3];
sx q[3];
rz(-2.9123902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7198221) q[2];
sx q[2];
rz(-0.82227102) q[2];
sx q[2];
rz(-2.3262809) q[2];
rz(-2.1067965) q[3];
sx q[3];
rz(-1.5765604) q[3];
sx q[3];
rz(-1.8243272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3141044) q[0];
sx q[0];
rz(-1.0359456) q[0];
sx q[0];
rz(-0.28717336) q[0];
rz(0.18889591) q[1];
sx q[1];
rz(-0.72383988) q[1];
sx q[1];
rz(-2.8093991) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4591601) q[0];
sx q[0];
rz(-0.81572616) q[0];
sx q[0];
rz(2.0072323) q[0];
rz(1.6261149) q[2];
sx q[2];
rz(-0.69960591) q[2];
sx q[2];
rz(2.827022) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.13739535) q[1];
sx q[1];
rz(-1.157837) q[1];
sx q[1];
rz(1.8505627) q[1];
rz(-pi) q[2];
rz(-0.62854564) q[3];
sx q[3];
rz(-2.5124031) q[3];
sx q[3];
rz(-2.7043846) q[3];
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
rz(-1.8509289) q[3];
sx q[3];
rz(-1.7920114) q[3];
sx q[3];
rz(-2.4850142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.055450913) q[0];
sx q[0];
rz(-2.3738326) q[0];
sx q[0];
rz(2.0822051) q[0];
rz(0.97958952) q[1];
sx q[1];
rz(-1.9444119) q[1];
sx q[1];
rz(0.25451452) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.018325018) q[0];
sx q[0];
rz(-1.8013445) q[0];
sx q[0];
rz(-0.075320764) q[0];
x q[1];
rz(-1.8149257) q[2];
sx q[2];
rz(-2.0460528) q[2];
sx q[2];
rz(-2.313386) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
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
x q[2];
rz(1.47154) q[3];
sx q[3];
rz(-1.1343065) q[3];
sx q[3];
rz(2.1193159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0467726) q[2];
sx q[2];
rz(-2.5929055) q[2];
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
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.027325252) q[0];
sx q[0];
rz(-1.1189168) q[0];
sx q[0];
rz(3.0615321) q[0];
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
