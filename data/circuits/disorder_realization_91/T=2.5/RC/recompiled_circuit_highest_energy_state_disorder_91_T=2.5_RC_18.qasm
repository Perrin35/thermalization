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
rz(0.36747992) q[0];
sx q[0];
rz(-1.3290661) q[0];
sx q[0];
rz(1.4867866) q[0];
rz(1.5005255) q[1];
sx q[1];
rz(-2.5379116) q[1];
sx q[1];
rz(-0.85890213) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6971815) q[0];
sx q[0];
rz(-1.0134122) q[0];
sx q[0];
rz(-1.0943221) q[0];
x q[1];
rz(2.3043465) q[2];
sx q[2];
rz(-2.4105802) q[2];
sx q[2];
rz(-1.8164509) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2341237) q[1];
sx q[1];
rz(-2.4501743) q[1];
sx q[1];
rz(1.3983634) q[1];
x q[2];
rz(-1.2626507) q[3];
sx q[3];
rz(-1.7761782) q[3];
sx q[3];
rz(2.8440203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9438802) q[2];
sx q[2];
rz(-0.79215017) q[2];
sx q[2];
rz(1.0666749) q[2];
rz(-0.094430447) q[3];
sx q[3];
rz(-0.61757278) q[3];
sx q[3];
rz(-2.7134231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18148947) q[0];
sx q[0];
rz(-2.8747989) q[0];
sx q[0];
rz(-1.2130523) q[0];
rz(-2.6890697) q[1];
sx q[1];
rz(-2.8892543) q[1];
sx q[1];
rz(2.4620893) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0328355) q[0];
sx q[0];
rz(-1.8201314) q[0];
sx q[0];
rz(-1.6052206) q[0];
rz(1.7361197) q[2];
sx q[2];
rz(-1.2836873) q[2];
sx q[2];
rz(1.5218671) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.22917381) q[1];
sx q[1];
rz(-2.3245735) q[1];
sx q[1];
rz(0.55025834) q[1];
x q[2];
rz(0.68778681) q[3];
sx q[3];
rz(-1.7264328) q[3];
sx q[3];
rz(2.9693672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.54000336) q[2];
sx q[2];
rz(-2.3175779) q[2];
sx q[2];
rz(-0.69671112) q[2];
rz(-1.2434897) q[3];
sx q[3];
rz(-1.1949298) q[3];
sx q[3];
rz(-2.5534326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8884647) q[0];
sx q[0];
rz(-1.0364113) q[0];
sx q[0];
rz(2.3179407) q[0];
rz(0.15692391) q[1];
sx q[1];
rz(-1.7053968) q[1];
sx q[1];
rz(2.7395111) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7347468) q[0];
sx q[0];
rz(-1.9663133) q[0];
sx q[0];
rz(2.9989373) q[0];
rz(-1.7355315) q[2];
sx q[2];
rz(-0.28898063) q[2];
sx q[2];
rz(-1.8584205) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7618708) q[1];
sx q[1];
rz(-2.5609697) q[1];
sx q[1];
rz(-2.1942744) q[1];
x q[2];
rz(-1.1810494) q[3];
sx q[3];
rz(-1.7754835) q[3];
sx q[3];
rz(-1.2677416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0154401) q[2];
sx q[2];
rz(-0.93706477) q[2];
sx q[2];
rz(0.13119571) q[2];
rz(-2.0845856) q[3];
sx q[3];
rz(-1.9481235) q[3];
sx q[3];
rz(1.2098275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7159202) q[0];
sx q[0];
rz(-1.965006) q[0];
sx q[0];
rz(-1.7536989) q[0];
rz(1.362494) q[1];
sx q[1];
rz(-1.8905996) q[1];
sx q[1];
rz(1.9349792) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8598452) q[0];
sx q[0];
rz(-2.3424405) q[0];
sx q[0];
rz(2.4837982) q[0];
rz(2.0974656) q[2];
sx q[2];
rz(-0.9894254) q[2];
sx q[2];
rz(-2.2033221) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8959103) q[1];
sx q[1];
rz(-0.31202641) q[1];
sx q[1];
rz(-2.8883977) q[1];
x q[2];
rz(2.016417) q[3];
sx q[3];
rz(-0.88332159) q[3];
sx q[3];
rz(-0.690122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6917307) q[2];
sx q[2];
rz(-1.3351771) q[2];
sx q[2];
rz(-0.94592363) q[2];
rz(-2.6876884) q[3];
sx q[3];
rz(-2.4920521) q[3];
sx q[3];
rz(-0.18979931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7396963) q[0];
sx q[0];
rz(-1.3119768) q[0];
sx q[0];
rz(-2.358118) q[0];
rz(-2.2354194) q[1];
sx q[1];
rz(-2.2515191) q[1];
sx q[1];
rz(-2.5772212) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63097341) q[0];
sx q[0];
rz(-1.4518018) q[0];
sx q[0];
rz(-1.4257159) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4323813) q[2];
sx q[2];
rz(-0.40988806) q[2];
sx q[2];
rz(0.60333911) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.61144343) q[1];
sx q[1];
rz(-0.93031973) q[1];
sx q[1];
rz(2.3534858) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7310745) q[3];
sx q[3];
rz(-2.2791282) q[3];
sx q[3];
rz(1.1089408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.924661) q[2];
sx q[2];
rz(-1.7380119) q[2];
sx q[2];
rz(-3.0972287) q[2];
rz(1.4102604) q[3];
sx q[3];
rz(-0.20874617) q[3];
sx q[3];
rz(-1.1948208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5656972) q[0];
sx q[0];
rz(-0.79224753) q[0];
sx q[0];
rz(-0.76882291) q[0];
rz(-2.5104751) q[1];
sx q[1];
rz(-1.2080071) q[1];
sx q[1];
rz(0.15359503) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054718356) q[0];
sx q[0];
rz(-1.2445893) q[0];
sx q[0];
rz(-0.98814641) q[0];
rz(1.8472438) q[2];
sx q[2];
rz(-0.9381367) q[2];
sx q[2];
rz(1.1142434) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.61443084) q[1];
sx q[1];
rz(-1.8519743) q[1];
sx q[1];
rz(1.6826732) q[1];
x q[2];
rz(1.6669964) q[3];
sx q[3];
rz(-2.3398388) q[3];
sx q[3];
rz(2.8303162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0742053) q[2];
sx q[2];
rz(-1.5437443) q[2];
sx q[2];
rz(0.81573168) q[2];
rz(-3.0873114) q[3];
sx q[3];
rz(-1.7547601) q[3];
sx q[3];
rz(1.774971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.5197007) q[0];
sx q[0];
rz(-2.0760355) q[0];
sx q[0];
rz(-0.52072293) q[0];
rz(-2.0479274) q[1];
sx q[1];
rz(-0.92898527) q[1];
sx q[1];
rz(1.7589794) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4224515) q[0];
sx q[0];
rz(-2.1364742) q[0];
sx q[0];
rz(-2.9702787) q[0];
x q[1];
rz(-0.6925631) q[2];
sx q[2];
rz(-0.64856883) q[2];
sx q[2];
rz(2.8798333) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3471038) q[1];
sx q[1];
rz(-2.1784049) q[1];
sx q[1];
rz(0.5083063) q[1];
x q[2];
rz(-0.46926559) q[3];
sx q[3];
rz(-1.7466144) q[3];
sx q[3];
rz(-1.3202536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9426721) q[2];
sx q[2];
rz(-0.3519381) q[2];
sx q[2];
rz(1.675763) q[2];
rz(0.28907019) q[3];
sx q[3];
rz(-1.7710268) q[3];
sx q[3];
rz(1.8606961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11956231) q[0];
sx q[0];
rz(-2.1896095) q[0];
sx q[0];
rz(-2.6493678) q[0];
rz(2.4960663) q[1];
sx q[1];
rz(-1.5954834) q[1];
sx q[1];
rz(-0.92799497) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6715901) q[0];
sx q[0];
rz(-2.8264464) q[0];
sx q[0];
rz(-2.0285839) q[0];
x q[1];
rz(2.3113584) q[2];
sx q[2];
rz(-2.6343759) q[2];
sx q[2];
rz(0.00032256034) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7240844) q[1];
sx q[1];
rz(-0.88005304) q[1];
sx q[1];
rz(-2.7907985) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2008529) q[3];
sx q[3];
rz(-1.2725012) q[3];
sx q[3];
rz(1.6443133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.32138985) q[2];
sx q[2];
rz(-1.5265744) q[2];
sx q[2];
rz(-0.75902933) q[2];
rz(1.8697033) q[3];
sx q[3];
rz(-1.3262409) q[3];
sx q[3];
rz(2.429764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(1.7976545) q[0];
sx q[0];
rz(-2.5929218) q[0];
sx q[0];
rz(0.77242533) q[0];
rz(1.3683569) q[1];
sx q[1];
rz(-1.1027579) q[1];
sx q[1];
rz(-0.30430421) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40201449) q[0];
sx q[0];
rz(-2.9106044) q[0];
sx q[0];
rz(-0.84257479) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7965806) q[2];
sx q[2];
rz(-1.5743557) q[2];
sx q[2];
rz(2.4709216) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.22812477) q[1];
sx q[1];
rz(-2.2992837) q[1];
sx q[1];
rz(-1.3291275) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7469962) q[3];
sx q[3];
rz(-1.6871042) q[3];
sx q[3];
rz(-1.1726086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7336537) q[2];
sx q[2];
rz(-2.759178) q[2];
sx q[2];
rz(1.1747053) q[2];
rz(-0.63742739) q[3];
sx q[3];
rz(-1.8791684) q[3];
sx q[3];
rz(-1.2229961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8680854) q[0];
sx q[0];
rz(-0.33184505) q[0];
sx q[0];
rz(2.3688431) q[0];
rz(2.6840774) q[1];
sx q[1];
rz(-1.3950149) q[1];
sx q[1];
rz(1.7083907) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49559418) q[0];
sx q[0];
rz(-2.3353205) q[0];
sx q[0];
rz(0.92732556) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9338701) q[2];
sx q[2];
rz(-1.6350897) q[2];
sx q[2];
rz(0.63426547) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6683176) q[1];
sx q[1];
rz(-2.6087481) q[1];
sx q[1];
rz(-2.1168296) q[1];
rz(0.97416157) q[3];
sx q[3];
rz(-1.0785558) q[3];
sx q[3];
rz(-2.8288768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.92051238) q[2];
sx q[2];
rz(-0.94197333) q[2];
sx q[2];
rz(2.4648049) q[2];
rz(-1.6109198) q[3];
sx q[3];
rz(-2.8407606) q[3];
sx q[3];
rz(-1.0845186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21657011) q[0];
sx q[0];
rz(-1.3022447) q[0];
sx q[0];
rz(-1.362823) q[0];
rz(-2.2577747) q[1];
sx q[1];
rz(-2.4088036) q[1];
sx q[1];
rz(2.1688681) q[1];
rz(-1.0013318) q[2];
sx q[2];
rz(-1.9591667) q[2];
sx q[2];
rz(-2.9137076) q[2];
rz(-1.8098214) q[3];
sx q[3];
rz(-2.0229133) q[3];
sx q[3];
rz(-2.3537707) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
