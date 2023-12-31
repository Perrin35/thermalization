OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7005641) q[0];
sx q[0];
rz(-1.9987885) q[0];
sx q[0];
rz(-1.9300652) q[0];
rz(6.05655) q[1];
sx q[1];
rz(1.5645138) q[1];
sx q[1];
rz(5.984879) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.053781833) q[0];
sx q[0];
rz(-2.2349173) q[0];
sx q[0];
rz(-1.640663) q[0];
x q[1];
rz(-1.9036129) q[2];
sx q[2];
rz(-1.4822072) q[2];
sx q[2];
rz(0.36117902) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8611697) q[1];
sx q[1];
rz(-2.0457595) q[1];
sx q[1];
rz(0.95649398) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.88793036) q[3];
sx q[3];
rz(-0.13266064) q[3];
sx q[3];
rz(-0.56548972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6478708) q[2];
sx q[2];
rz(-1.2080668) q[2];
sx q[2];
rz(-0.99386627) q[2];
rz(-2.1422051) q[3];
sx q[3];
rz(-1.9013654) q[3];
sx q[3];
rz(2.5527111) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4988929) q[0];
sx q[0];
rz(-1.2058586) q[0];
sx q[0];
rz(3.1233741) q[0];
rz(-0.81623626) q[1];
sx q[1];
rz(-1.0304334) q[1];
sx q[1];
rz(0.47168628) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1175849) q[0];
sx q[0];
rz(-1.4837259) q[0];
sx q[0];
rz(2.9136806) q[0];
rz(2.6366028) q[2];
sx q[2];
rz(-1.8171176) q[2];
sx q[2];
rz(1.1438952) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.7447409) q[1];
sx q[1];
rz(-0.62780118) q[1];
sx q[1];
rz(2.6746034) q[1];
rz(-pi) q[2];
rz(-1.3470115) q[3];
sx q[3];
rz(-2.5219678) q[3];
sx q[3];
rz(2.1686045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.54962426) q[2];
sx q[2];
rz(-1.9886118) q[2];
sx q[2];
rz(-2.1726051) q[2];
rz(2.5668868) q[3];
sx q[3];
rz(-0.55137268) q[3];
sx q[3];
rz(2.1000752) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4330924) q[0];
sx q[0];
rz(-2.0860465) q[0];
sx q[0];
rz(2.95978) q[0];
rz(1.1026985) q[1];
sx q[1];
rz(-1.6405374) q[1];
sx q[1];
rz(1.4556494) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91711125) q[0];
sx q[0];
rz(-0.79229504) q[0];
sx q[0];
rz(2.2729421) q[0];
rz(-pi) q[1];
rz(1.1481029) q[2];
sx q[2];
rz(-2.2659677) q[2];
sx q[2];
rz(-2.9425651) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1959343) q[1];
sx q[1];
rz(-0.12609005) q[1];
sx q[1];
rz(-0.083421589) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6772179) q[3];
sx q[3];
rz(-2.1933746) q[3];
sx q[3];
rz(-2.4604083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.87749798) q[2];
sx q[2];
rz(-1.654518) q[2];
sx q[2];
rz(-0.96763119) q[2];
rz(-2.4140221) q[3];
sx q[3];
rz(-1.8811767) q[3];
sx q[3];
rz(2.9038866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2847292) q[0];
sx q[0];
rz(-0.52607042) q[0];
sx q[0];
rz(2.5033584) q[0];
rz(1.1278641) q[1];
sx q[1];
rz(-0.82740873) q[1];
sx q[1];
rz(-1.2329873) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1283135) q[0];
sx q[0];
rz(-2.0648801) q[0];
sx q[0];
rz(-1.4923151) q[0];
rz(-pi) q[1];
rz(2.4013176) q[2];
sx q[2];
rz(-0.7358272) q[2];
sx q[2];
rz(2.1528113) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4302664) q[1];
sx q[1];
rz(-1.3923936) q[1];
sx q[1];
rz(-3.1049411) q[1];
rz(-pi) q[2];
rz(0.22709417) q[3];
sx q[3];
rz(-0.89082754) q[3];
sx q[3];
rz(2.1122776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2234852) q[2];
sx q[2];
rz(-1.8635668) q[2];
sx q[2];
rz(2.3275862) q[2];
rz(1.043184) q[3];
sx q[3];
rz(-2.510575) q[3];
sx q[3];
rz(1.1842747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2994613) q[0];
sx q[0];
rz(-1.786754) q[0];
sx q[0];
rz(-2.2498851) q[0];
rz(1.8978329) q[1];
sx q[1];
rz(-1.7638821) q[1];
sx q[1];
rz(2.9290501) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.606013) q[0];
sx q[0];
rz(-1.5545462) q[0];
sx q[0];
rz(-1.5910801) q[0];
rz(-pi) q[1];
rz(-1.5151305) q[2];
sx q[2];
rz(-0.60534436) q[2];
sx q[2];
rz(2.2948613) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.37413874) q[1];
sx q[1];
rz(-2.1530188) q[1];
sx q[1];
rz(-2.2351082) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0605293) q[3];
sx q[3];
rz(-0.79499309) q[3];
sx q[3];
rz(2.0612962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1084958) q[2];
sx q[2];
rz(-2.1704845) q[2];
sx q[2];
rz(-2.6203716) q[2];
rz(-1.7565578) q[3];
sx q[3];
rz(-1.9135467) q[3];
sx q[3];
rz(-1.2683755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8591156) q[0];
sx q[0];
rz(-0.92882597) q[0];
sx q[0];
rz(-0.22512063) q[0];
rz(1.7865932) q[1];
sx q[1];
rz(-1.0083818) q[1];
sx q[1];
rz(-0.37757847) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7901944) q[0];
sx q[0];
rz(-1.5488008) q[0];
sx q[0];
rz(1.3047332) q[0];
rz(-pi) q[1];
rz(0.2874561) q[2];
sx q[2];
rz(-1.6171347) q[2];
sx q[2];
rz(-0.031678274) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7403455) q[1];
sx q[1];
rz(-1.1402854) q[1];
sx q[1];
rz(3.0262448) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6013006) q[3];
sx q[3];
rz(-0.77611938) q[3];
sx q[3];
rz(-0.12150773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1569415) q[2];
sx q[2];
rz(-2.3539383) q[2];
sx q[2];
rz(-0.48103508) q[2];
rz(2.7379819) q[3];
sx q[3];
rz(-2.1026881) q[3];
sx q[3];
rz(0.31479442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0890546) q[0];
sx q[0];
rz(-0.60482329) q[0];
sx q[0];
rz(-2.9470434) q[0];
rz(0.21952195) q[1];
sx q[1];
rz(-1.4621282) q[1];
sx q[1];
rz(0.25442466) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5146778) q[0];
sx q[0];
rz(-1.8906381) q[0];
sx q[0];
rz(-3.0266648) q[0];
rz(-pi) q[1];
rz(1.3527855) q[2];
sx q[2];
rz(-2.5180452) q[2];
sx q[2];
rz(-0.64507285) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9421778) q[1];
sx q[1];
rz(-0.53655469) q[1];
sx q[1];
rz(0.55120991) q[1];
rz(0.80612225) q[3];
sx q[3];
rz(-0.61180173) q[3];
sx q[3];
rz(-0.34477371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.97757942) q[2];
sx q[2];
rz(-1.7501202) q[2];
sx q[2];
rz(-1.8257726) q[2];
rz(0.87604648) q[3];
sx q[3];
rz(-0.13893572) q[3];
sx q[3];
rz(2.1379437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9309689) q[0];
sx q[0];
rz(-2.7656778) q[0];
sx q[0];
rz(1.6865431) q[0];
rz(-0.82398206) q[1];
sx q[1];
rz(-2.1897557) q[1];
sx q[1];
rz(-1.5751858) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4644509) q[0];
sx q[0];
rz(-1.8832708) q[0];
sx q[0];
rz(0.72014767) q[0];
rz(1.0303866) q[2];
sx q[2];
rz(-1.4166797) q[2];
sx q[2];
rz(-0.87481462) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.52289256) q[1];
sx q[1];
rz(-2.61781) q[1];
sx q[1];
rz(-0.79843847) q[1];
x q[2];
rz(2.6353587) q[3];
sx q[3];
rz(-1.1971548) q[3];
sx q[3];
rz(-0.95844275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2660797) q[2];
sx q[2];
rz(-1.2337039) q[2];
sx q[2];
rz(-2.8640462) q[2];
rz(-1.6905486) q[3];
sx q[3];
rz(-0.45193672) q[3];
sx q[3];
rz(-2.1267166) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9797416) q[0];
sx q[0];
rz(-0.45270544) q[0];
sx q[0];
rz(1.45654) q[0];
rz(0.62943554) q[1];
sx q[1];
rz(-1.1673735) q[1];
sx q[1];
rz(1.1368407) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0757383) q[0];
sx q[0];
rz(-2.542001) q[0];
sx q[0];
rz(-1.0435186) q[0];
rz(0.85176977) q[2];
sx q[2];
rz(-0.53817828) q[2];
sx q[2];
rz(0.9966419) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6362308) q[1];
sx q[1];
rz(-0.33126918) q[1];
sx q[1];
rz(-1.865571) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4056899) q[3];
sx q[3];
rz(-1.2888442) q[3];
sx q[3];
rz(-1.499093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4920766) q[2];
sx q[2];
rz(-1.7121544) q[2];
sx q[2];
rz(-2.1006404) q[2];
rz(3.1395636) q[3];
sx q[3];
rz(-0.40015951) q[3];
sx q[3];
rz(0.31203976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3025538) q[0];
sx q[0];
rz(-0.22452393) q[0];
sx q[0];
rz(-0.94605207) q[0];
rz(-2.229915) q[1];
sx q[1];
rz(-1.2152351) q[1];
sx q[1];
rz(2.5295703) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48869041) q[0];
sx q[0];
rz(-0.86581794) q[0];
sx q[0];
rz(2.5253354) q[0];
rz(-pi) q[1];
rz(1.3304747) q[2];
sx q[2];
rz(-1.2346621) q[2];
sx q[2];
rz(-2.0517595) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.40572383) q[1];
sx q[1];
rz(-1.4475665) q[1];
sx q[1];
rz(0.65995364) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9480115) q[3];
sx q[3];
rz(-1.0672788) q[3];
sx q[3];
rz(0.52770381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0845906) q[2];
sx q[2];
rz(-2.4995063) q[2];
sx q[2];
rz(-2.4882312) q[2];
rz(0.35081321) q[3];
sx q[3];
rz(-1.5272798) q[3];
sx q[3];
rz(-2.4408834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54031298) q[0];
sx q[0];
rz(-1.5355587) q[0];
sx q[0];
rz(-3.0328947) q[0];
rz(-0.75469771) q[1];
sx q[1];
rz(-1.3194059) q[1];
sx q[1];
rz(-1.5059765) q[1];
rz(-1.0036219) q[2];
sx q[2];
rz(-2.1688609) q[2];
sx q[2];
rz(-1.4458956) q[2];
rz(-1.8176953) q[3];
sx q[3];
rz(-0.3331475) q[3];
sx q[3];
rz(-3.0039136) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
