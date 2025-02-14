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
rz(1.4242564) q[0];
sx q[0];
rz(1.9219226) q[0];
sx q[0];
rz(9.9998247) q[0];
rz(0.14248928) q[1];
sx q[1];
rz(-1.9228851) q[1];
sx q[1];
rz(-0.86427468) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2144379) q[0];
sx q[0];
rz(-1.7933284) q[0];
sx q[0];
rz(0.19630918) q[0];
rz(-pi) q[1];
rz(-1.7792542) q[2];
sx q[2];
rz(-0.64178665) q[2];
sx q[2];
rz(1.889495) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.083347224) q[1];
sx q[1];
rz(-1.5687404) q[1];
sx q[1];
rz(0.53002341) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5745886) q[3];
sx q[3];
rz(-2.0006913) q[3];
sx q[3];
rz(-2.5741732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.064726) q[2];
sx q[2];
rz(-1.4996424) q[2];
sx q[2];
rz(-1.0092674) q[2];
rz(0.81389728) q[3];
sx q[3];
rz(-0.32865694) q[3];
sx q[3];
rz(0.3391268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-0.086804248) q[0];
sx q[0];
rz(-1.4084933) q[0];
sx q[0];
rz(-2.8993697) q[0];
rz(1.9107266) q[1];
sx q[1];
rz(-0.63175285) q[1];
sx q[1];
rz(2.3435074) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.674192) q[0];
sx q[0];
rz(-0.62819203) q[0];
sx q[0];
rz(2.9487361) q[0];
rz(-pi) q[1];
rz(-1.961444) q[2];
sx q[2];
rz(-0.36875781) q[2];
sx q[2];
rz(-2.7386896) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.68951571) q[1];
sx q[1];
rz(-0.59714666) q[1];
sx q[1];
rz(-2.4680016) q[1];
rz(2.734058) q[3];
sx q[3];
rz(-0.38713249) q[3];
sx q[3];
rz(2.3705496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.65608394) q[2];
sx q[2];
rz(-1.069243) q[2];
sx q[2];
rz(2.3089224) q[2];
rz(0.63140702) q[3];
sx q[3];
rz(-2.4000945) q[3];
sx q[3];
rz(-0.00028636534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0176508) q[0];
sx q[0];
rz(-0.90921062) q[0];
sx q[0];
rz(-0.18774524) q[0];
rz(-0.84367696) q[1];
sx q[1];
rz(-1.2742821) q[1];
sx q[1];
rz(-0.63293308) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5143941) q[0];
sx q[0];
rz(-0.36982515) q[0];
sx q[0];
rz(-1.373795) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.723188) q[2];
sx q[2];
rz(-2.4103571) q[2];
sx q[2];
rz(-2.3894073) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.46309457) q[1];
sx q[1];
rz(-0.60623525) q[1];
sx q[1];
rz(0.54529584) q[1];
rz(-pi) q[2];
rz(0.12944451) q[3];
sx q[3];
rz(-1.2264381) q[3];
sx q[3];
rz(-0.5705006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8741499) q[2];
sx q[2];
rz(-0.94620693) q[2];
sx q[2];
rz(-1.2713185) q[2];
rz(1.6455796) q[3];
sx q[3];
rz(-1.6451903) q[3];
sx q[3];
rz(-2.1700844) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2520168) q[0];
sx q[0];
rz(-0.75279623) q[0];
sx q[0];
rz(0.65943199) q[0];
rz(-1.8239498) q[1];
sx q[1];
rz(-1.8023856) q[1];
sx q[1];
rz(-0.79016322) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0280608) q[0];
sx q[0];
rz(-1.4469997) q[0];
sx q[0];
rz(-2.8040299) q[0];
rz(-1.6035346) q[2];
sx q[2];
rz(-1.590609) q[2];
sx q[2];
rz(-1.1184426) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9875264) q[1];
sx q[1];
rz(-0.48168698) q[1];
sx q[1];
rz(-0.14194004) q[1];
rz(-pi) q[2];
rz(-2.5874917) q[3];
sx q[3];
rz(-1.3375207) q[3];
sx q[3];
rz(1.3194039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.059375199) q[2];
sx q[2];
rz(-1.6148115) q[2];
sx q[2];
rz(2.8455632) q[2];
rz(-0.56240231) q[3];
sx q[3];
rz(-2.2735169) q[3];
sx q[3];
rz(-0.11399046) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.247308) q[0];
sx q[0];
rz(-1.1881275) q[0];
sx q[0];
rz(0.2071912) q[0];
rz(-2.1098792) q[1];
sx q[1];
rz(-1.9887911) q[1];
sx q[1];
rz(-1.656104) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9821859) q[0];
sx q[0];
rz(-0.80538087) q[0];
sx q[0];
rz(-0.55434395) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1274952) q[2];
sx q[2];
rz(-2.5336899) q[2];
sx q[2];
rz(0.098092484) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.19073701) q[1];
sx q[1];
rz(-1.2225635) q[1];
sx q[1];
rz(-1.7067043) q[1];
rz(0.84019227) q[3];
sx q[3];
rz(-1.6512135) q[3];
sx q[3];
rz(1.059746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1032054) q[2];
sx q[2];
rz(-2.3296671) q[2];
sx q[2];
rz(-1.1078328) q[2];
rz(-2.2191018) q[3];
sx q[3];
rz(-1.4100217) q[3];
sx q[3];
rz(0.24423519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6237727) q[0];
sx q[0];
rz(-0.51384846) q[0];
sx q[0];
rz(0.22860953) q[0];
rz(2.8672583) q[1];
sx q[1];
rz(-0.81293303) q[1];
sx q[1];
rz(0.45023227) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0666445) q[0];
sx q[0];
rz(-2.3111176) q[0];
sx q[0];
rz(1.226783) q[0];
rz(-pi) q[1];
rz(-1.6041555) q[2];
sx q[2];
rz(-0.22936996) q[2];
sx q[2];
rz(-2.0663886) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6630312) q[1];
sx q[1];
rz(-1.3197462) q[1];
sx q[1];
rz(-0.51086225) q[1];
rz(-pi) q[2];
rz(0.26900503) q[3];
sx q[3];
rz(-1.9719567) q[3];
sx q[3];
rz(-1.3022193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.36394128) q[2];
sx q[2];
rz(-2.5886017) q[2];
sx q[2];
rz(-0.027776329) q[2];
rz(-0.76644301) q[3];
sx q[3];
rz(-1.1602297) q[3];
sx q[3];
rz(-2.4465223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22208333) q[0];
sx q[0];
rz(-1.6351901) q[0];
sx q[0];
rz(2.5191504) q[0];
rz(-2.4355603) q[1];
sx q[1];
rz(-2.249554) q[1];
sx q[1];
rz(2.5546254) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3496124) q[0];
sx q[0];
rz(-2.0087578) q[0];
sx q[0];
rz(-0.95923202) q[0];
x q[1];
rz(0.47759159) q[2];
sx q[2];
rz(-0.53672635) q[2];
sx q[2];
rz(2.0411185) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1026969) q[1];
sx q[1];
rz(-1.4998815) q[1];
sx q[1];
rz(0.50693477) q[1];
x q[2];
rz(3.1312084) q[3];
sx q[3];
rz(-2.2828013) q[3];
sx q[3];
rz(-0.66085789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.075835) q[2];
sx q[2];
rz(-1.4328052) q[2];
sx q[2];
rz(-0.59448057) q[2];
rz(-0.25932702) q[3];
sx q[3];
rz(-1.2649053) q[3];
sx q[3];
rz(1.9243141) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1031621) q[0];
sx q[0];
rz(-2.2362464) q[0];
sx q[0];
rz(2.5627947) q[0];
rz(-1.3102866) q[1];
sx q[1];
rz(-1.2625887) q[1];
sx q[1];
rz(2.328918) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90234251) q[0];
sx q[0];
rz(-0.50241155) q[0];
sx q[0];
rz(0.30186304) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9695639) q[2];
sx q[2];
rz(-0.59894717) q[2];
sx q[2];
rz(-2.2676165) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.26749565) q[1];
sx q[1];
rz(-2.2363642) q[1];
sx q[1];
rz(-1.1049306) q[1];
x q[2];
rz(-0.50847362) q[3];
sx q[3];
rz(-2.0805854) q[3];
sx q[3];
rz(-2.2126861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7909214) q[2];
sx q[2];
rz(-1.8427883) q[2];
sx q[2];
rz(0.92362967) q[2];
rz(-1.0203429) q[3];
sx q[3];
rz(-0.89710051) q[3];
sx q[3];
rz(-0.11788192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0807121) q[0];
sx q[0];
rz(-1.7055644) q[0];
sx q[0];
rz(-1.4870148) q[0];
rz(-0.05038536) q[1];
sx q[1];
rz(-0.81704187) q[1];
sx q[1];
rz(2.494716) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2912783) q[0];
sx q[0];
rz(-2.3731557) q[0];
sx q[0];
rz(-0.9356229) q[0];
rz(-2.0503309) q[2];
sx q[2];
rz(-2.1796436) q[2];
sx q[2];
rz(3.0225168) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8343492) q[1];
sx q[1];
rz(-1.3615047) q[1];
sx q[1];
rz(2.7268975) q[1];
rz(-pi) q[2];
rz(-0.96109747) q[3];
sx q[3];
rz(-2.8831579) q[3];
sx q[3];
rz(-1.0040545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.50463027) q[2];
sx q[2];
rz(-1.5550104) q[2];
sx q[2];
rz(0.75616765) q[2];
rz(1.6715096) q[3];
sx q[3];
rz(-0.78592891) q[3];
sx q[3];
rz(2.3362931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.037755448) q[0];
sx q[0];
rz(-1.2498195) q[0];
sx q[0];
rz(2.0334429) q[0];
rz(0.26421079) q[1];
sx q[1];
rz(-1.6725531) q[1];
sx q[1];
rz(2.5392551) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7450934) q[0];
sx q[0];
rz(-1.4317997) q[0];
sx q[0];
rz(-1.8462015) q[0];
x q[1];
rz(0.71508821) q[2];
sx q[2];
rz(-2.9093766) q[2];
sx q[2];
rz(-2.2975936) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2341799) q[1];
sx q[1];
rz(-0.8846332) q[1];
sx q[1];
rz(-3.0790672) q[1];
x q[2];
rz(2.4906795) q[3];
sx q[3];
rz(-1.7554566) q[3];
sx q[3];
rz(2.1391275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5084874) q[2];
sx q[2];
rz(-2.5999531) q[2];
sx q[2];
rz(-0.87149054) q[2];
rz(-1.2095215) q[3];
sx q[3];
rz(-0.28035527) q[3];
sx q[3];
rz(2.5476294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7021983) q[0];
sx q[0];
rz(-1.7560503) q[0];
sx q[0];
rz(-0.51881292) q[0];
rz(-2.5595472) q[1];
sx q[1];
rz(-2.0379635) q[1];
sx q[1];
rz(1.4720974) q[1];
rz(1.9459361) q[2];
sx q[2];
rz(-2.2715501) q[2];
sx q[2];
rz(3.0987433) q[2];
rz(-2.4529759) q[3];
sx q[3];
rz(-1.4388765) q[3];
sx q[3];
rz(2.8394113) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
