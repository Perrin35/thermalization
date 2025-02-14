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
rz(-1.9788111) q[0];
sx q[0];
rz(-2.2012308) q[0];
sx q[0];
rz(2.9475687) q[0];
rz(-3.6858978) q[1];
sx q[1];
rz(1.5679918) q[1];
sx q[1];
rz(12.369649) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17795632) q[0];
sx q[0];
rz(-2.5546409) q[0];
sx q[0];
rz(1.3474083) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9241184) q[2];
sx q[2];
rz(-1.5494692) q[2];
sx q[2];
rz(2.3623717) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2620723) q[1];
sx q[1];
rz(-1.889887) q[1];
sx q[1];
rz(1.2882922) q[1];
rz(-pi) q[2];
rz(1.3338577) q[3];
sx q[3];
rz(-0.58962017) q[3];
sx q[3];
rz(-0.20649621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0804312) q[2];
sx q[2];
rz(-1.4208527) q[2];
sx q[2];
rz(-0.013193456) q[2];
rz(2.0773928) q[3];
sx q[3];
rz(-1.0366169) q[3];
sx q[3];
rz(-2.9707151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5203633) q[0];
sx q[0];
rz(-2.5907488) q[0];
sx q[0];
rz(-0.4826104) q[0];
rz(0.47166011) q[1];
sx q[1];
rz(-2.1444247) q[1];
sx q[1];
rz(-0.68798033) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4783966) q[0];
sx q[0];
rz(-1.8747878) q[0];
sx q[0];
rz(-1.3711956) q[0];
rz(0.28765388) q[2];
sx q[2];
rz(-2.322933) q[2];
sx q[2];
rz(2.5729736) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5291366) q[1];
sx q[1];
rz(-1.0998526) q[1];
sx q[1];
rz(-1.6199084) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.08556871) q[3];
sx q[3];
rz(-1.7734756) q[3];
sx q[3];
rz(0.50093716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1241577) q[2];
sx q[2];
rz(-2.0686801) q[2];
sx q[2];
rz(-2.9166481) q[2];
rz(-1.0791091) q[3];
sx q[3];
rz(-2.2033117) q[3];
sx q[3];
rz(0.45043954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9709622) q[0];
sx q[0];
rz(-2.5522975) q[0];
sx q[0];
rz(1.2901837) q[0];
rz(1.2782485) q[1];
sx q[1];
rz(-1.1561013) q[1];
sx q[1];
rz(1.8589171) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3564094) q[0];
sx q[0];
rz(-0.98726596) q[0];
sx q[0];
rz(-0.43278354) q[0];
x q[1];
rz(3.01802) q[2];
sx q[2];
rz(-1.4915492) q[2];
sx q[2];
rz(1.0956956) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8508529) q[1];
sx q[1];
rz(-3.1147482) q[1];
sx q[1];
rz(2.31183) q[1];
rz(2.3219789) q[3];
sx q[3];
rz(-0.96246877) q[3];
sx q[3];
rz(1.5386594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4625385) q[2];
sx q[2];
rz(-1.5306229) q[2];
sx q[2];
rz(-2.6666717) q[2];
rz(1.9654407) q[3];
sx q[3];
rz(-0.85634309) q[3];
sx q[3];
rz(2.8600051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3731641) q[0];
sx q[0];
rz(-1.7730862) q[0];
sx q[0];
rz(2.9737293) q[0];
rz(-2.8236875) q[1];
sx q[1];
rz(-1.2891506) q[1];
sx q[1];
rz(1.0344523) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8694591) q[0];
sx q[0];
rz(-1.3087166) q[0];
sx q[0];
rz(0.036152187) q[0];
rz(-pi) q[1];
x q[1];
rz(2.690763) q[2];
sx q[2];
rz(-0.73720142) q[2];
sx q[2];
rz(-2.1429755) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2047836) q[1];
sx q[1];
rz(-2.2310871) q[1];
sx q[1];
rz(-2.5461063) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.529783) q[3];
sx q[3];
rz(-1.6535078) q[3];
sx q[3];
rz(1.3349319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1075403) q[2];
sx q[2];
rz(-0.62959051) q[2];
sx q[2];
rz(0.26629392) q[2];
rz(2.9722424) q[3];
sx q[3];
rz(-2.0784056) q[3];
sx q[3];
rz(0.17899409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.87194815) q[0];
sx q[0];
rz(-0.18505159) q[0];
sx q[0];
rz(-1.8484775) q[0];
rz(0.070501892) q[1];
sx q[1];
rz(-0.87424707) q[1];
sx q[1];
rz(-0.806113) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9490813) q[0];
sx q[0];
rz(-1.7777788) q[0];
sx q[0];
rz(3.0144318) q[0];
rz(1.9367123) q[2];
sx q[2];
rz(-0.15934773) q[2];
sx q[2];
rz(1.510335) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.55743417) q[1];
sx q[1];
rz(-2.5100252) q[1];
sx q[1];
rz(-0.33051349) q[1];
rz(-pi) q[2];
rz(0.3573017) q[3];
sx q[3];
rz(-0.56824486) q[3];
sx q[3];
rz(-0.33333353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.21141323) q[2];
sx q[2];
rz(-1.6665919) q[2];
sx q[2];
rz(0.30964568) q[2];
rz(-0.51112255) q[3];
sx q[3];
rz(-0.57356727) q[3];
sx q[3];
rz(0.82093325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9822134) q[0];
sx q[0];
rz(-0.62726227) q[0];
sx q[0];
rz(-2.7045767) q[0];
rz(2.6868195) q[1];
sx q[1];
rz(-0.79650703) q[1];
sx q[1];
rz(0.88266596) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2734786) q[0];
sx q[0];
rz(-1.4820288) q[0];
sx q[0];
rz(-1.8411536) q[0];
rz(-pi) q[1];
rz(3.109819) q[2];
sx q[2];
rz(-1.1924044) q[2];
sx q[2];
rz(0.057003958) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5875579) q[1];
sx q[1];
rz(-1.2393426) q[1];
sx q[1];
rz(1.0629088) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.168247) q[3];
sx q[3];
rz(-0.69434598) q[3];
sx q[3];
rz(-2.2996817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.63558811) q[2];
sx q[2];
rz(-1.4281861) q[2];
sx q[2];
rz(-0.44090718) q[2];
rz(2.0697557) q[3];
sx q[3];
rz(-0.25522885) q[3];
sx q[3];
rz(0.60060874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1275682) q[0];
sx q[0];
rz(-1.2034282) q[0];
sx q[0];
rz(-1.2822275) q[0];
rz(1.0857238) q[1];
sx q[1];
rz(-1.5241357) q[1];
sx q[1];
rz(-1.6283183) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0659637) q[0];
sx q[0];
rz(-2.3018078) q[0];
sx q[0];
rz(-0.15878598) q[0];
x q[1];
rz(0.70019987) q[2];
sx q[2];
rz(-1.321903) q[2];
sx q[2];
rz(-0.85485103) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3320184) q[1];
sx q[1];
rz(-0.93832131) q[1];
sx q[1];
rz(-2.6676086) q[1];
rz(-0.75700883) q[3];
sx q[3];
rz(-1.2881713) q[3];
sx q[3];
rz(0.69863897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6699803) q[2];
sx q[2];
rz(-2.4358304) q[2];
sx q[2];
rz(-0.82028779) q[2];
rz(2.7465076) q[3];
sx q[3];
rz(-0.79212752) q[3];
sx q[3];
rz(-0.61409605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.26246) q[0];
sx q[0];
rz(-1.4985871) q[0];
sx q[0];
rz(-2.5347128) q[0];
rz(-2.795769) q[1];
sx q[1];
rz(-1.1864097) q[1];
sx q[1];
rz(2.3343559) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31921747) q[0];
sx q[0];
rz(-1.9142879) q[0];
sx q[0];
rz(0.10256711) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.26793311) q[2];
sx q[2];
rz(-1.9606646) q[2];
sx q[2];
rz(-0.34696503) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9620558) q[1];
sx q[1];
rz(-2.7946094) q[1];
sx q[1];
rz(-0.68106236) q[1];
x q[2];
rz(-0.59698481) q[3];
sx q[3];
rz(-1.0506949) q[3];
sx q[3];
rz(-0.4515243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.81749934) q[2];
sx q[2];
rz(-1.9976511) q[2];
sx q[2];
rz(2.6452046) q[2];
rz(-0.084970623) q[3];
sx q[3];
rz(-0.43046633) q[3];
sx q[3];
rz(0.5808723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9896511) q[0];
sx q[0];
rz(-1.4284416) q[0];
sx q[0];
rz(-0.38452837) q[0];
rz(2.8893068) q[1];
sx q[1];
rz(-1.3832046) q[1];
sx q[1];
rz(-1.5581473) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91004086) q[0];
sx q[0];
rz(-1.1052026) q[0];
sx q[0];
rz(1.5814073) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6884638) q[2];
sx q[2];
rz(-1.8431547) q[2];
sx q[2];
rz(2.5982369) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6391771) q[1];
sx q[1];
rz(-2.1599837) q[1];
sx q[1];
rz(-3.1261954) q[1];
rz(-pi) q[2];
rz(-2.2966398) q[3];
sx q[3];
rz(-0.87741565) q[3];
sx q[3];
rz(-1.4385731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.79534507) q[2];
sx q[2];
rz(-1.2582422) q[2];
sx q[2];
rz(-2.0015049) q[2];
rz(1.7831066) q[3];
sx q[3];
rz(-1.8890817) q[3];
sx q[3];
rz(-0.31806773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7110905) q[0];
sx q[0];
rz(-1.4416114) q[0];
sx q[0];
rz(0.39518133) q[0];
rz(-1.9197865) q[1];
sx q[1];
rz(-2.2876574) q[1];
sx q[1];
rz(-0.88776678) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7956314) q[0];
sx q[0];
rz(-2.2610616) q[0];
sx q[0];
rz(1.2119341) q[0];
x q[1];
rz(-3.1380464) q[2];
sx q[2];
rz(-1.0958091) q[2];
sx q[2];
rz(-1.3925683) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.50248442) q[1];
sx q[1];
rz(-0.35668761) q[1];
sx q[1];
rz(-1.4977895) q[1];
rz(-pi) q[2];
rz(1.2597144) q[3];
sx q[3];
rz(-1.4691989) q[3];
sx q[3];
rz(-1.7465308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.44783529) q[2];
sx q[2];
rz(-1.3112023) q[2];
sx q[2];
rz(2.5698404) q[2];
rz(-1.6311496) q[3];
sx q[3];
rz(-1.4331199) q[3];
sx q[3];
rz(1.800764) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0157264) q[0];
sx q[0];
rz(-1.8342352) q[0];
sx q[0];
rz(-2.9974708) q[0];
rz(1.9563328) q[1];
sx q[1];
rz(-0.85694255) q[1];
sx q[1];
rz(-1.8640929) q[1];
rz(0.38707564) q[2];
sx q[2];
rz(-1.7995926) q[2];
sx q[2];
rz(2.8298701) q[2];
rz(-1.2492467) q[3];
sx q[3];
rz(-0.95148181) q[3];
sx q[3];
rz(-0.8368809) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
