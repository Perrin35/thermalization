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
rz(1.3951294) q[0];
sx q[0];
rz(3.8962235) q[0];
sx q[0];
rz(8.3409283) q[0];
rz(-0.87767449) q[1];
sx q[1];
rz(-1.2134774) q[1];
sx q[1];
rz(-0.47668639) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9636757) q[0];
sx q[0];
rz(-2.6612894) q[0];
sx q[0];
rz(-2.1558574) q[0];
rz(-pi) q[1];
rz(-1.1083015) q[2];
sx q[2];
rz(-1.1574189) q[2];
sx q[2];
rz(0.55962901) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.649316) q[1];
sx q[1];
rz(-1.5536123) q[1];
sx q[1];
rz(-1.3130929) q[1];
rz(2.6891119) q[3];
sx q[3];
rz(-0.47305952) q[3];
sx q[3];
rz(-0.61247952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.013926355) q[2];
sx q[2];
rz(-1.650859) q[2];
sx q[2];
rz(2.9909383) q[2];
rz(-1.6023747) q[3];
sx q[3];
rz(-2.7645002) q[3];
sx q[3];
rz(2.8902875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.056331228) q[0];
sx q[0];
rz(-1.3717317) q[0];
sx q[0];
rz(2.7714609) q[0];
rz(-2.6689957) q[1];
sx q[1];
rz(-0.31817803) q[1];
sx q[1];
rz(-0.95742375) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2070302) q[0];
sx q[0];
rz(-2.0683204) q[0];
sx q[0];
rz(-3.050404) q[0];
rz(-1.9198869) q[2];
sx q[2];
rz(-0.47609509) q[2];
sx q[2];
rz(-0.27175909) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.53471662) q[1];
sx q[1];
rz(-0.83321111) q[1];
sx q[1];
rz(-1.4556769) q[1];
x q[2];
rz(1.9022835) q[3];
sx q[3];
rz(-0.9210862) q[3];
sx q[3];
rz(3.0420764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.0051673278) q[2];
sx q[2];
rz(-2.397126) q[2];
sx q[2];
rz(2.1180507) q[2];
rz(-1.7510022) q[3];
sx q[3];
rz(-1.7460881) q[3];
sx q[3];
rz(1.8243779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.7904974) q[0];
sx q[0];
rz(-2.1320765) q[0];
sx q[0];
rz(0.56280953) q[0];
rz(-1.6273181) q[1];
sx q[1];
rz(-1.4311675) q[1];
sx q[1];
rz(3.0624342) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3252727) q[0];
sx q[0];
rz(-2.5421341) q[0];
sx q[0];
rz(0.58167235) q[0];
rz(-pi) q[1];
rz(1.6116033) q[2];
sx q[2];
rz(-2.0611603) q[2];
sx q[2];
rz(-0.23863579) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.2217949) q[1];
sx q[1];
rz(-2.256003) q[1];
sx q[1];
rz(2.5849839) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8757365) q[3];
sx q[3];
rz(-0.99627021) q[3];
sx q[3];
rz(2.4158626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2341653) q[2];
sx q[2];
rz(-1.5831466) q[2];
sx q[2];
rz(-1.3321715) q[2];
rz(-2.7857156) q[3];
sx q[3];
rz(-1.9477113) q[3];
sx q[3];
rz(2.1607384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6818162) q[0];
sx q[0];
rz(-1.1695319) q[0];
sx q[0];
rz(-2.0303149) q[0];
rz(-0.92012826) q[1];
sx q[1];
rz(-1.5524813) q[1];
sx q[1];
rz(-0.2535325) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4944015) q[0];
sx q[0];
rz(-1.8888942) q[0];
sx q[0];
rz(0.53089106) q[0];
x q[1];
rz(-2.5229613) q[2];
sx q[2];
rz(-0.19605532) q[2];
sx q[2];
rz(-0.679099) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.82099709) q[1];
sx q[1];
rz(-1.7554469) q[1];
sx q[1];
rz(0.76652938) q[1];
rz(1.0456123) q[3];
sx q[3];
rz(-2.5730657) q[3];
sx q[3];
rz(0.15318682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.8186875) q[2];
sx q[2];
rz(-2.5762168) q[2];
sx q[2];
rz(-1.6820071) q[2];
rz(-0.5736351) q[3];
sx q[3];
rz(-1.2021659) q[3];
sx q[3];
rz(-0.044376317) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6673073) q[0];
sx q[0];
rz(-2.0720338) q[0];
sx q[0];
rz(1.3473508) q[0];
rz(0.59066331) q[1];
sx q[1];
rz(-1.9213516) q[1];
sx q[1];
rz(-1.4264872) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2329946) q[0];
sx q[0];
rz(-0.54366606) q[0];
sx q[0];
rz(-2.3999016) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9587742) q[2];
sx q[2];
rz(-0.96599865) q[2];
sx q[2];
rz(-1.6590349) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.062078309) q[1];
sx q[1];
rz(-0.2582363) q[1];
sx q[1];
rz(0.18194992) q[1];
x q[2];
rz(2.9699202) q[3];
sx q[3];
rz(-0.35997691) q[3];
sx q[3];
rz(1.6115007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9023989) q[2];
sx q[2];
rz(-2.0812483) q[2];
sx q[2];
rz(-1.1725461) q[2];
rz(3.0555365) q[3];
sx q[3];
rz(-0.9525758) q[3];
sx q[3];
rz(-0.15611592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0758783) q[0];
sx q[0];
rz(-1.5944163) q[0];
sx q[0];
rz(-2.2789047) q[0];
rz(2.8942096) q[1];
sx q[1];
rz(-1.8377973) q[1];
sx q[1];
rz(2.6695796) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87403546) q[0];
sx q[0];
rz(-0.78473259) q[0];
sx q[0];
rz(2.6536056) q[0];
rz(-pi) q[1];
rz(-1.0649933) q[2];
sx q[2];
rz(-1.7031295) q[2];
sx q[2];
rz(1.9605317) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0473049) q[1];
sx q[1];
rz(-1.3940812) q[1];
sx q[1];
rz(-0.30728886) q[1];
rz(2.9293044) q[3];
sx q[3];
rz(-0.95930305) q[3];
sx q[3];
rz(1.6595728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.90999675) q[2];
sx q[2];
rz(-0.96698499) q[2];
sx q[2];
rz(1.3989353) q[2];
rz(-1.4553962) q[3];
sx q[3];
rz(-2.3536436) q[3];
sx q[3];
rz(0.56070352) q[3];
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
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48651925) q[0];
sx q[0];
rz(-1.1460679) q[0];
sx q[0];
rz(-0.51862496) q[0];
rz(-0.71867603) q[1];
sx q[1];
rz(-0.51281723) q[1];
sx q[1];
rz(1.3291043) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3603921) q[0];
sx q[0];
rz(-2.5787163) q[0];
sx q[0];
rz(-0.84873523) q[0];
x q[1];
rz(-3.091264) q[2];
sx q[2];
rz(-1.8476881) q[2];
sx q[2];
rz(-0.30093873) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.74187255) q[1];
sx q[1];
rz(-0.73556566) q[1];
sx q[1];
rz(-1.5217785) q[1];
rz(1.6006599) q[3];
sx q[3];
rz(-0.58664188) q[3];
sx q[3];
rz(-1.7778974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.79080498) q[2];
sx q[2];
rz(-0.69956508) q[2];
sx q[2];
rz(3.1034071) q[2];
rz(-0.83032483) q[3];
sx q[3];
rz(-1.7540951) q[3];
sx q[3];
rz(-0.047688095) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6787978) q[0];
sx q[0];
rz(-2.9252453) q[0];
sx q[0];
rz(-1.7713254) q[0];
rz(2.7554152) q[1];
sx q[1];
rz(-1.4048856) q[1];
sx q[1];
rz(-0.6699627) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44035916) q[0];
sx q[0];
rz(-2.0645077) q[0];
sx q[0];
rz(-0.53892737) q[0];
rz(-pi) q[1];
rz(0.092548056) q[2];
sx q[2];
rz(-2.4227648) q[2];
sx q[2];
rz(0.75564849) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.48310623) q[1];
sx q[1];
rz(-1.6838594) q[1];
sx q[1];
rz(-2.1674564) q[1];
rz(-1.7427722) q[3];
sx q[3];
rz(-2.9265518) q[3];
sx q[3];
rz(0.57143927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2204444) q[2];
sx q[2];
rz(-2.8577652) q[2];
sx q[2];
rz(2.0849126) q[2];
rz(1.7250666) q[3];
sx q[3];
rz(-1.2456015) q[3];
sx q[3];
rz(-1.2904803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.073762745) q[0];
sx q[0];
rz(-2.1011598) q[0];
sx q[0];
rz(2.1347866) q[0];
rz(0.49268588) q[1];
sx q[1];
rz(-1.410781) q[1];
sx q[1];
rz(-0.83008343) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7573038) q[0];
sx q[0];
rz(-0.28097935) q[0];
sx q[0];
rz(0.12019867) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7858429) q[2];
sx q[2];
rz(-2.1587662) q[2];
sx q[2];
rz(2.7933592) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0574903) q[1];
sx q[1];
rz(-1.3907281) q[1];
sx q[1];
rz(1.7258641) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1287635) q[3];
sx q[3];
rz(-2.0137824) q[3];
sx q[3];
rz(-1.6010546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.81862187) q[2];
sx q[2];
rz(-1.3016204) q[2];
sx q[2];
rz(1.0969561) q[2];
rz(1.8638301) q[3];
sx q[3];
rz(-0.68456972) q[3];
sx q[3];
rz(2.4975615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85035664) q[0];
sx q[0];
rz(-1.1024029) q[0];
sx q[0];
rz(-0.37503234) q[0];
rz(0.11861435) q[1];
sx q[1];
rz(-1.7053441) q[1];
sx q[1];
rz(-0.92432252) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18436954) q[0];
sx q[0];
rz(-0.65378377) q[0];
sx q[0];
rz(-0.48982805) q[0];
x q[1];
rz(2.5293328) q[2];
sx q[2];
rz(-2.3697402) q[2];
sx q[2];
rz(-1.5651694) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.19434793) q[1];
sx q[1];
rz(-1.2100265) q[1];
sx q[1];
rz(1.6014577) q[1];
rz(-pi) q[2];
rz(-2.6931346) q[3];
sx q[3];
rz(-1.8844386) q[3];
sx q[3];
rz(2.3109316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2529926) q[2];
sx q[2];
rz(-1.000095) q[2];
sx q[2];
rz(2.2484153) q[2];
rz(-1.1183974) q[3];
sx q[3];
rz(-1.018254) q[3];
sx q[3];
rz(0.65677381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7669582) q[0];
sx q[0];
rz(-2.0603016) q[0];
sx q[0];
rz(1.7745071) q[0];
rz(-2.959666) q[1];
sx q[1];
rz(-0.55768273) q[1];
sx q[1];
rz(2.2348977) q[1];
rz(-1.2557658) q[2];
sx q[2];
rz(-1.4834822) q[2];
sx q[2];
rz(-0.33374141) q[2];
rz(2.6245194) q[3];
sx q[3];
rz(-1.2096249) q[3];
sx q[3];
rz(0.67740868) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
