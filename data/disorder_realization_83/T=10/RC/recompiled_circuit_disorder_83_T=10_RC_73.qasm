OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5542334) q[0];
sx q[0];
rz(-2.1589307) q[0];
sx q[0];
rz(0.76173705) q[0];
rz(0.76454437) q[1];
sx q[1];
rz(-2.0643056) q[1];
sx q[1];
rz(0.74365562) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7997768) q[0];
sx q[0];
rz(-1.3306949) q[0];
sx q[0];
rz(0.067702985) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15723575) q[2];
sx q[2];
rz(-0.6586282) q[2];
sx q[2];
rz(3.0008891) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.88120645) q[1];
sx q[1];
rz(-1.2938061) q[1];
sx q[1];
rz(-1.9828412) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3841964) q[3];
sx q[3];
rz(-0.36986923) q[3];
sx q[3];
rz(-2.6922525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7068229) q[2];
sx q[2];
rz(-2.7396024) q[2];
sx q[2];
rz(-3.0337231) q[2];
rz(-2.9989631) q[3];
sx q[3];
rz(-1.7248036) q[3];
sx q[3];
rz(0.67255783) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07664872) q[0];
sx q[0];
rz(-0.80524421) q[0];
sx q[0];
rz(-0.23072492) q[0];
rz(-1.327286) q[1];
sx q[1];
rz(-0.67115152) q[1];
sx q[1];
rz(0.040963106) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36386585) q[0];
sx q[0];
rz(-1.4085318) q[0];
sx q[0];
rz(-1.8198387) q[0];
rz(-pi) q[1];
x q[1];
rz(0.83611739) q[2];
sx q[2];
rz(-0.40514075) q[2];
sx q[2];
rz(2.1737614) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7296655) q[1];
sx q[1];
rz(-1.1150868) q[1];
sx q[1];
rz(0.63974849) q[1];
rz(-pi) q[2];
rz(-1.9713692) q[3];
sx q[3];
rz(-1.922732) q[3];
sx q[3];
rz(-2.7924201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.162398) q[2];
sx q[2];
rz(-1.4731864) q[2];
sx q[2];
rz(-2.5734148) q[2];
rz(2.6290821) q[3];
sx q[3];
rz(-0.5957225) q[3];
sx q[3];
rz(0.56186831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7115962) q[0];
sx q[0];
rz(-1.5189518) q[0];
sx q[0];
rz(-0.45021737) q[0];
rz(-1.846116) q[1];
sx q[1];
rz(-1.9772915) q[1];
sx q[1];
rz(0.67726642) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60904658) q[0];
sx q[0];
rz(-0.17621528) q[0];
sx q[0];
rz(2.8410068) q[0];
rz(-pi) q[1];
rz(1.5718939) q[2];
sx q[2];
rz(-1.5771616) q[2];
sx q[2];
rz(2.9874143) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9159569) q[1];
sx q[1];
rz(-1.7491873) q[1];
sx q[1];
rz(-0.13326463) q[1];
rz(-pi) q[2];
rz(-1.9618109) q[3];
sx q[3];
rz(-2.0174332) q[3];
sx q[3];
rz(-3.1004268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.21851097) q[2];
sx q[2];
rz(-1.9873025) q[2];
sx q[2];
rz(-0.046860524) q[2];
rz(-0.81165195) q[3];
sx q[3];
rz(-0.67957687) q[3];
sx q[3];
rz(-0.17015447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99455225) q[0];
sx q[0];
rz(-0.11163286) q[0];
sx q[0];
rz(0.051483367) q[0];
rz(-0.4908081) q[1];
sx q[1];
rz(-0.97508109) q[1];
sx q[1];
rz(-1.9690008) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8797982) q[0];
sx q[0];
rz(-2.2248631) q[0];
sx q[0];
rz(1.7394702) q[0];
rz(-pi) q[1];
rz(1.555221) q[2];
sx q[2];
rz(-0.96484631) q[2];
sx q[2];
rz(1.3941927) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.40201515) q[1];
sx q[1];
rz(-0.47164729) q[1];
sx q[1];
rz(-2.2775843) q[1];
rz(-pi) q[2];
rz(-1.4106393) q[3];
sx q[3];
rz(-0.43678108) q[3];
sx q[3];
rz(0.76524759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2525758) q[2];
sx q[2];
rz(-2.2106407) q[2];
sx q[2];
rz(-0.48689294) q[2];
rz(-1.1934818) q[3];
sx q[3];
rz(-0.32961696) q[3];
sx q[3];
rz(3.1161599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.837773) q[0];
sx q[0];
rz(-1.0205512) q[0];
sx q[0];
rz(1.3254962) q[0];
rz(-0.59108132) q[1];
sx q[1];
rz(-0.68060827) q[1];
sx q[1];
rz(0.65471929) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7703169) q[0];
sx q[0];
rz(-1.9195942) q[0];
sx q[0];
rz(2.6692997) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7227313) q[2];
sx q[2];
rz(-1.3992116) q[2];
sx q[2];
rz(-0.78531839) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3760738) q[1];
sx q[1];
rz(-2.3507833) q[1];
sx q[1];
rz(-2.5942624) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4773024) q[3];
sx q[3];
rz(-1.9545385) q[3];
sx q[3];
rz(1.547471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6902265) q[2];
sx q[2];
rz(-1.2008685) q[2];
sx q[2];
rz(2.5704685) q[2];
rz(-0.58978224) q[3];
sx q[3];
rz(-0.46001205) q[3];
sx q[3];
rz(-0.067972876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(0.94474435) q[0];
sx q[0];
rz(-1.1905043) q[0];
sx q[0];
rz(-0.29904547) q[0];
rz(-1.8213182) q[1];
sx q[1];
rz(-0.25779217) q[1];
sx q[1];
rz(-1.6437795) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4371571) q[0];
sx q[0];
rz(-0.91059443) q[0];
sx q[0];
rz(3.0755088) q[0];
rz(-pi) q[1];
x q[1];
rz(0.57967474) q[2];
sx q[2];
rz(-1.056676) q[2];
sx q[2];
rz(-2.3953953) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7945929) q[1];
sx q[1];
rz(-1.1113249) q[1];
sx q[1];
rz(-0.32230349) q[1];
rz(-2.8619814) q[3];
sx q[3];
rz(-1.8108484) q[3];
sx q[3];
rz(-0.9595426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.51612878) q[2];
sx q[2];
rz(-1.4279782) q[2];
sx q[2];
rz(-0.027475474) q[2];
rz(-2.6190858) q[3];
sx q[3];
rz(-2.3404739) q[3];
sx q[3];
rz(0.81645042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41247535) q[0];
sx q[0];
rz(-1.1183879) q[0];
sx q[0];
rz(-3.0274042) q[0];
rz(2.1633637) q[1];
sx q[1];
rz(-2.5949635) q[1];
sx q[1];
rz(-0.79089975) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7366911) q[0];
sx q[0];
rz(-1.4762523) q[0];
sx q[0];
rz(0.94876429) q[0];
rz(-1.2355455) q[2];
sx q[2];
rz(-2.2860048) q[2];
sx q[2];
rz(0.4244948) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0485718) q[1];
sx q[1];
rz(-1.2067817) q[1];
sx q[1];
rz(-0.084333468) q[1];
rz(-pi) q[2];
rz(-2.2579455) q[3];
sx q[3];
rz(-1.2847319) q[3];
sx q[3];
rz(1.4253695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2699282) q[2];
sx q[2];
rz(-1.1130788) q[2];
sx q[2];
rz(-0.86501914) q[2];
rz(0.67251742) q[3];
sx q[3];
rz(-1.1285684) q[3];
sx q[3];
rz(3.045936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4928116) q[0];
sx q[0];
rz(-0.5195986) q[0];
sx q[0];
rz(-0.46736026) q[0];
rz(2.6043747) q[1];
sx q[1];
rz(-2.1610114) q[1];
sx q[1];
rz(0.25407243) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92914903) q[0];
sx q[0];
rz(-1.375276) q[0];
sx q[0];
rz(-1.6385727) q[0];
x q[1];
rz(1.5723096) q[2];
sx q[2];
rz(-1.5760033) q[2];
sx q[2];
rz(1.7092012) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.2368187) q[1];
sx q[1];
rz(-0.33729759) q[1];
sx q[1];
rz(-0.39223139) q[1];
rz(-pi) q[2];
rz(-2.5423126) q[3];
sx q[3];
rz(-1.5044971) q[3];
sx q[3];
rz(1.9950206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.49999985) q[2];
sx q[2];
rz(-0.77074146) q[2];
sx q[2];
rz(0.33561486) q[2];
rz(-2.8149758) q[3];
sx q[3];
rz(-0.89544046) q[3];
sx q[3];
rz(-1.4183104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0582054) q[0];
sx q[0];
rz(-2.931262) q[0];
sx q[0];
rz(-1.0539508) q[0];
rz(-2.9846233) q[1];
sx q[1];
rz(-1.7228246) q[1];
sx q[1];
rz(-2.1597247) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0197524) q[0];
sx q[0];
rz(-1.4593908) q[0];
sx q[0];
rz(0.41282546) q[0];
rz(-pi) q[1];
rz(-1.4023196) q[2];
sx q[2];
rz(-0.57541621) q[2];
sx q[2];
rz(-1.9926496) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.30584221) q[1];
sx q[1];
rz(-1.3853449) q[1];
sx q[1];
rz(-0.058538392) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8544159) q[3];
sx q[3];
rz(-2.4099726) q[3];
sx q[3];
rz(-1.2958131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.92423576) q[2];
sx q[2];
rz(-2.0841667) q[2];
sx q[2];
rz(-0.99739933) q[2];
rz(-2.635397) q[3];
sx q[3];
rz(-2.180407) q[3];
sx q[3];
rz(-3.1072646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2845594) q[0];
sx q[0];
rz(-0.71145809) q[0];
sx q[0];
rz(-2.6249264) q[0];
rz(0.38756469) q[1];
sx q[1];
rz(-1.060408) q[1];
sx q[1];
rz(-2.8732252) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4386913) q[0];
sx q[0];
rz(-0.99150204) q[0];
sx q[0];
rz(-2.2772574) q[0];
x q[1];
rz(0.95558138) q[2];
sx q[2];
rz(-0.68253839) q[2];
sx q[2];
rz(0.62894422) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8908773) q[1];
sx q[1];
rz(-2.4836575) q[1];
sx q[1];
rz(-0.4472181) q[1];
x q[2];
rz(0.99276944) q[3];
sx q[3];
rz(-1.3475932) q[3];
sx q[3];
rz(2.2900667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9364075) q[2];
sx q[2];
rz(-2.3582017) q[2];
sx q[2];
rz(-0.56383413) q[2];
rz(-1.9514203) q[3];
sx q[3];
rz(-0.95919132) q[3];
sx q[3];
rz(0.64175516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.664809) q[0];
sx q[0];
rz(-1.2510779) q[0];
sx q[0];
rz(-1.0673987) q[0];
rz(-1.8019567) q[1];
sx q[1];
rz(-1.4383153) q[1];
sx q[1];
rz(-1.7972606) q[1];
rz(1.8742758) q[2];
sx q[2];
rz(-1.0721285) q[2];
sx q[2];
rz(-0.68821651) q[2];
rz(-0.31233882) q[3];
sx q[3];
rz(-1.3795508) q[3];
sx q[3];
rz(3.0637904) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
