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
rz(-2.7741127) q[0];
sx q[0];
rz(-1.8125266) q[0];
sx q[0];
rz(-1.4867866) q[0];
rz(-1.6410671) q[1];
sx q[1];
rz(2.5379116) q[1];
sx q[1];
rz(11.707468) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2175332) q[0];
sx q[0];
rz(-2.4251063) q[0];
sx q[0];
rz(-0.63436809) q[0];
x q[1];
rz(-2.1584081) q[2];
sx q[2];
rz(-1.1074083) q[2];
sx q[2];
rz(-2.7962229) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6851226) q[1];
sx q[1];
rz(-0.89159144) q[1];
sx q[1];
rz(-0.14107708) q[1];
rz(-pi) q[2];
rz(-1.8789419) q[3];
sx q[3];
rz(-1.7761782) q[3];
sx q[3];
rz(0.29757231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9438802) q[2];
sx q[2];
rz(-0.79215017) q[2];
sx q[2];
rz(1.0666749) q[2];
rz(0.094430447) q[3];
sx q[3];
rz(-2.5240199) q[3];
sx q[3];
rz(-2.7134231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9601032) q[0];
sx q[0];
rz(-2.8747989) q[0];
sx q[0];
rz(1.9285404) q[0];
rz(-2.6890697) q[1];
sx q[1];
rz(-2.8892543) q[1];
sx q[1];
rz(2.4620893) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54645854) q[0];
sx q[0];
rz(-1.537437) q[0];
sx q[0];
rz(-0.24947687) q[0];
x q[1];
rz(-1.7361197) q[2];
sx q[2];
rz(-1.8579053) q[2];
sx q[2];
rz(-1.6197255) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.22917381) q[1];
sx q[1];
rz(-2.3245735) q[1];
sx q[1];
rz(2.5913343) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3704471) q[3];
sx q[3];
rz(-0.89290038) q[3];
sx q[3];
rz(1.6163449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.54000336) q[2];
sx q[2];
rz(-0.82401472) q[2];
sx q[2];
rz(-0.69671112) q[2];
rz(-1.8981029) q[3];
sx q[3];
rz(-1.1949298) q[3];
sx q[3];
rz(2.5534326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
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
rz(2.9846687) q[1];
sx q[1];
rz(-1.7053968) q[1];
sx q[1];
rz(-2.7395111) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4068459) q[0];
sx q[0];
rz(-1.9663133) q[0];
sx q[0];
rz(0.14265539) q[0];
rz(1.8560748) q[2];
sx q[2];
rz(-1.6175468) q[2];
sx q[2];
rz(-0.12960427) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3797219) q[1];
sx q[1];
rz(-0.58062299) q[1];
sx q[1];
rz(0.94731829) q[1];
x q[2];
rz(-1.1810494) q[3];
sx q[3];
rz(-1.7754835) q[3];
sx q[3];
rz(1.8738511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0154401) q[2];
sx q[2];
rz(-0.93706477) q[2];
sx q[2];
rz(-3.0103969) q[2];
rz(1.0570071) q[3];
sx q[3];
rz(-1.9481235) q[3];
sx q[3];
rz(-1.9317651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4256725) q[0];
sx q[0];
rz(-1.965006) q[0];
sx q[0];
rz(-1.7536989) q[0];
rz(1.7790986) q[1];
sx q[1];
rz(-1.250993) q[1];
sx q[1];
rz(1.9349792) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8598452) q[0];
sx q[0];
rz(-2.3424405) q[0];
sx q[0];
rz(2.4837982) q[0];
rz(0.64995857) q[2];
sx q[2];
rz(-1.1372677) q[2];
sx q[2];
rz(-2.1999846) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5111419) q[1];
sx q[1];
rz(-1.2690374) q[1];
sx q[1];
rz(-1.4901699) q[1];
rz(-pi) q[2];
rz(-2.6582022) q[3];
sx q[3];
rz(-2.3425274) q[3];
sx q[3];
rz(1.8062325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6917307) q[2];
sx q[2];
rz(-1.8064156) q[2];
sx q[2];
rz(0.94592363) q[2];
rz(-2.6876884) q[3];
sx q[3];
rz(-2.4920521) q[3];
sx q[3];
rz(-0.18979931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7396963) q[0];
sx q[0];
rz(-1.8296158) q[0];
sx q[0];
rz(0.78347462) q[0];
rz(2.2354194) q[1];
sx q[1];
rz(-0.89007354) q[1];
sx q[1];
rz(0.56437147) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6221478) q[0];
sx q[0];
rz(-0.18737245) q[0];
sx q[0];
rz(-0.87981422) q[0];
rz(-pi) q[1];
rz(-0.059877574) q[2];
sx q[2];
rz(-1.9765325) q[2];
sx q[2];
rz(-0.45260383) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7184938) q[1];
sx q[1];
rz(-0.96984399) q[1];
sx q[1];
rz(2.3312631) q[1];
rz(-pi) q[2];
rz(2.9574017) q[3];
sx q[3];
rz(-2.4184265) q[3];
sx q[3];
rz(1.3525055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2169317) q[2];
sx q[2];
rz(-1.7380119) q[2];
sx q[2];
rz(0.04436392) q[2];
rz(-1.7313322) q[3];
sx q[3];
rz(-2.9328465) q[3];
sx q[3];
rz(1.1948208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5656972) q[0];
sx q[0];
rz(-0.79224753) q[0];
sx q[0];
rz(0.76882291) q[0];
rz(-0.63111758) q[1];
sx q[1];
rz(-1.2080071) q[1];
sx q[1];
rz(-0.15359503) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7241868) q[0];
sx q[0];
rz(-2.119078) q[0];
sx q[0];
rz(0.38492695) q[0];
rz(-1.8472438) q[2];
sx q[2];
rz(-2.203456) q[2];
sx q[2];
rz(1.1142434) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.61443084) q[1];
sx q[1];
rz(-1.2896184) q[1];
sx q[1];
rz(1.6826732) q[1];
rz(-pi) q[2];
x q[2];
rz(0.098922313) q[3];
sx q[3];
rz(-2.3677845) q[3];
sx q[3];
rz(0.1734002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.067387335) q[2];
sx q[2];
rz(-1.5437443) q[2];
sx q[2];
rz(-0.81573168) q[2];
rz(3.0873114) q[3];
sx q[3];
rz(-1.7547601) q[3];
sx q[3];
rz(1.3666216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.62189198) q[0];
sx q[0];
rz(-1.0655572) q[0];
sx q[0];
rz(0.52072293) q[0];
rz(1.0936652) q[1];
sx q[1];
rz(-2.2126074) q[1];
sx q[1];
rz(-1.7589794) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1102229) q[0];
sx q[0];
rz(-0.58833226) q[0];
sx q[0];
rz(-1.3084685) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.6925631) q[2];
sx q[2];
rz(-0.64856883) q[2];
sx q[2];
rz(2.8798333) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7944889) q[1];
sx q[1];
rz(-0.96318775) q[1];
sx q[1];
rz(0.5083063) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3741854) q[3];
sx q[3];
rz(-2.0322606) q[3];
sx q[3];
rz(-2.9795071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9426721) q[2];
sx q[2];
rz(-0.3519381) q[2];
sx q[2];
rz(1.675763) q[2];
rz(-2.8525225) q[3];
sx q[3];
rz(-1.3705658) q[3];
sx q[3];
rz(-1.8606961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11956231) q[0];
sx q[0];
rz(-2.1896095) q[0];
sx q[0];
rz(0.49222487) q[0];
rz(-0.64552632) q[1];
sx q[1];
rz(-1.5954834) q[1];
sx q[1];
rz(2.2135977) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6715901) q[0];
sx q[0];
rz(-0.3151463) q[0];
sx q[0];
rz(-1.1130087) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3113584) q[2];
sx q[2];
rz(-2.6343759) q[2];
sx q[2];
rz(-0.00032256034) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.89615738) q[1];
sx q[1];
rz(-0.76145524) q[1];
sx q[1];
rz(1.9647962) q[1];
x q[2];
rz(-0.31854872) q[3];
sx q[3];
rz(-1.9236698) q[3];
sx q[3];
rz(0.039963756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.32138985) q[2];
sx q[2];
rz(-1.6150183) q[2];
sx q[2];
rz(0.75902933) q[2];
rz(-1.8697033) q[3];
sx q[3];
rz(-1.3262409) q[3];
sx q[3];
rz(-2.429764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7976545) q[0];
sx q[0];
rz(-0.54867083) q[0];
sx q[0];
rz(-2.3691673) q[0];
rz(-1.7732357) q[1];
sx q[1];
rz(-1.1027579) q[1];
sx q[1];
rz(-0.30430421) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9979663) q[0];
sx q[0];
rz(-1.7425101) q[0];
sx q[0];
rz(2.986326) q[0];
x q[1];
rz(0.0036520521) q[2];
sx q[2];
rz(-1.7965791) q[2];
sx q[2];
rz(0.89930764) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.22812477) q[1];
sx q[1];
rz(-0.84230891) q[1];
sx q[1];
rz(-1.8124652) q[1];
rz(-0.39459644) q[3];
sx q[3];
rz(-1.4544885) q[3];
sx q[3];
rz(1.1726086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.407939) q[2];
sx q[2];
rz(-0.38241461) q[2];
sx q[2];
rz(1.9668874) q[2];
rz(0.63742739) q[3];
sx q[3];
rz(-1.2624242) q[3];
sx q[3];
rz(-1.2229961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-1.8680854) q[0];
sx q[0];
rz(-0.33184505) q[0];
sx q[0];
rz(2.3688431) q[0];
rz(-2.6840774) q[1];
sx q[1];
rz(-1.7465778) q[1];
sx q[1];
rz(1.7083907) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59638176) q[0];
sx q[0];
rz(-1.1229674) q[0];
sx q[0];
rz(-0.87558545) q[0];
rz(-pi) q[1];
rz(-1.2077226) q[2];
sx q[2];
rz(-1.506503) q[2];
sx q[2];
rz(0.63426547) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0876574) q[1];
sx q[1];
rz(-2.0198576) q[1];
sx q[1];
rz(2.8444071) q[1];
x q[2];
rz(0.57515842) q[3];
sx q[3];
rz(-2.0888101) q[3];
sx q[3];
rz(0.94747551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2210803) q[2];
sx q[2];
rz(-2.1996193) q[2];
sx q[2];
rz(-2.4648049) q[2];
rz(-1.5306728) q[3];
sx q[3];
rz(-2.8407606) q[3];
sx q[3];
rz(1.0845186) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21657011) q[0];
sx q[0];
rz(-1.8393479) q[0];
sx q[0];
rz(1.7787697) q[0];
rz(2.2577747) q[1];
sx q[1];
rz(-0.73278905) q[1];
sx q[1];
rz(-0.97272452) q[1];
rz(0.45223805) q[2];
sx q[2];
rz(-1.0483208) q[2];
sx q[2];
rz(2.0365289) q[2];
rz(-0.46355214) q[3];
sx q[3];
rz(-1.3561854) q[3];
sx q[3];
rz(-0.88903313) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
