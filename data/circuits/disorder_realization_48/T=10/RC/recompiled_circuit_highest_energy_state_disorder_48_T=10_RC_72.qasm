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
rz(-1.7464632) q[0];
sx q[0];
rz(-0.75463086) q[0];
sx q[0];
rz(1.0838497) q[0];
rz(-0.87767449) q[1];
sx q[1];
rz(-1.2134774) q[1];
sx q[1];
rz(2.6649063) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81945588) q[0];
sx q[0];
rz(-1.9662204) q[0];
sx q[0];
rz(-0.28015341) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0332912) q[2];
sx q[2];
rz(-1.1574189) q[2];
sx q[2];
rz(-0.55962901) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4922766) q[1];
sx q[1];
rz(-1.5879803) q[1];
sx q[1];
rz(1.3130929) q[1];
rz(-pi) q[2];
rz(-0.4313978) q[3];
sx q[3];
rz(-1.7713303) q[3];
sx q[3];
rz(1.7748743) q[3];
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
rz(1.6023747) q[3];
sx q[3];
rz(-2.7645002) q[3];
sx q[3];
rz(0.25130513) q[3];
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
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0852614) q[0];
sx q[0];
rz(-1.3717317) q[0];
sx q[0];
rz(-2.7714609) q[0];
rz(-2.6689957) q[1];
sx q[1];
rz(-2.8234146) q[1];
sx q[1];
rz(0.95742375) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32015285) q[0];
sx q[0];
rz(-1.6509045) q[0];
sx q[0];
rz(1.0715241) q[0];
rz(1.1195807) q[2];
sx q[2];
rz(-1.4133845) q[2];
sx q[2];
rz(-1.5296641) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1831187) q[1];
sx q[1];
rz(-1.4856824) q[1];
sx q[1];
rz(0.74089153) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7368746) q[3];
sx q[3];
rz(-2.4232695) q[3];
sx q[3];
rz(-0.61678929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.0051673278) q[2];
sx q[2];
rz(-0.74446669) q[2];
sx q[2];
rz(1.0235419) q[2];
rz(1.3905904) q[3];
sx q[3];
rz(-1.3955045) q[3];
sx q[3];
rz(1.3172147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35109529) q[0];
sx q[0];
rz(-1.0095162) q[0];
sx q[0];
rz(0.56280953) q[0];
rz(1.6273181) q[1];
sx q[1];
rz(-1.7104251) q[1];
sx q[1];
rz(-0.079158457) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.889733) q[0];
sx q[0];
rz(-1.2556228) q[0];
sx q[0];
rz(-0.51879518) q[0];
x q[1];
rz(2.6508826) q[2];
sx q[2];
rz(-1.5348002) q[2];
sx q[2];
rz(-1.3513868) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9197977) q[1];
sx q[1];
rz(-0.88558965) q[1];
sx q[1];
rz(-2.5849839) q[1];
rz(0.78000809) q[3];
sx q[3];
rz(-2.2714408) q[3];
sx q[3];
rz(2.8740945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9074273) q[2];
sx q[2];
rz(-1.558446) q[2];
sx q[2];
rz(1.8094212) q[2];
rz(2.7857156) q[3];
sx q[3];
rz(-1.1938813) q[3];
sx q[3];
rz(2.1607384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.45977649) q[0];
sx q[0];
rz(-1.9720607) q[0];
sx q[0];
rz(-2.0303149) q[0];
rz(-2.2214644) q[1];
sx q[1];
rz(-1.5524813) q[1];
sx q[1];
rz(-2.8880602) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1052083) q[0];
sx q[0];
rz(-1.0691133) q[0];
sx q[0];
rz(1.2060449) q[0];
rz(-pi) q[1];
rz(0.61863135) q[2];
sx q[2];
rz(-0.19605532) q[2];
sx q[2];
rz(2.4624937) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.56139466) q[1];
sx q[1];
rz(-0.78402482) q[1];
sx q[1];
rz(2.8785588) q[1];
rz(2.0959804) q[3];
sx q[3];
rz(-2.5730657) q[3];
sx q[3];
rz(2.9884058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.8186875) q[2];
sx q[2];
rz(-2.5762168) q[2];
sx q[2];
rz(1.6820071) q[2];
rz(0.5736351) q[3];
sx q[3];
rz(-1.2021659) q[3];
sx q[3];
rz(-3.0972163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(0.4742853) q[0];
sx q[0];
rz(-2.0720338) q[0];
sx q[0];
rz(1.3473508) q[0];
rz(2.5509293) q[1];
sx q[1];
rz(-1.2202411) q[1];
sx q[1];
rz(-1.4264872) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0027699) q[0];
sx q[0];
rz(-1.213824) q[0];
sx q[0];
rz(-2.7223552) q[0];
rz(-1.9587742) q[2];
sx q[2];
rz(-2.175594) q[2];
sx q[2];
rz(-1.4825578) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0795143) q[1];
sx q[1];
rz(-2.8833564) q[1];
sx q[1];
rz(-0.18194992) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.17167244) q[3];
sx q[3];
rz(-0.35997691) q[3];
sx q[3];
rz(1.6115007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.23919375) q[2];
sx q[2];
rz(-2.0812483) q[2];
sx q[2];
rz(1.1725461) q[2];
rz(-0.086056195) q[3];
sx q[3];
rz(-2.1890169) q[3];
sx q[3];
rz(0.15611592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(1.0657144) q[0];
sx q[0];
rz(-1.5944163) q[0];
sx q[0];
rz(-2.2789047) q[0];
rz(2.8942096) q[1];
sx q[1];
rz(-1.8377973) q[1];
sx q[1];
rz(-0.47201306) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0855477) q[0];
sx q[0];
rz(-1.9084832) q[0];
sx q[0];
rz(-2.4187536) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8389147) q[2];
sx q[2];
rz(-2.6202218) q[2];
sx q[2];
rz(-2.9857295) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6093543) q[1];
sx q[1];
rz(-1.2684457) q[1];
sx q[1];
rz(-1.7560033) q[1];
rz(-pi) q[2];
x q[2];
rz(0.94859) q[3];
sx q[3];
rz(-1.3974117) q[3];
sx q[3];
rz(-0.21188785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2315959) q[2];
sx q[2];
rz(-0.96698499) q[2];
sx q[2];
rz(1.3989353) q[2];
rz(1.4553962) q[3];
sx q[3];
rz(-2.3536436) q[3];
sx q[3];
rz(2.5808891) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48651925) q[0];
sx q[0];
rz(-1.1460679) q[0];
sx q[0];
rz(0.51862496) q[0];
rz(-2.4229166) q[1];
sx q[1];
rz(-2.6287754) q[1];
sx q[1];
rz(-1.8124883) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78120056) q[0];
sx q[0];
rz(-0.5628764) q[0];
sx q[0];
rz(-0.84873523) q[0];
rz(-pi) q[1];
rz(-1.2935712) q[2];
sx q[2];
rz(-1.5223862) q[2];
sx q[2];
rz(1.2836266) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3997201) q[1];
sx q[1];
rz(-0.73556566) q[1];
sx q[1];
rz(-1.6198141) q[1];
rz(-pi) q[2];
rz(-1.6006599) q[3];
sx q[3];
rz(-0.58664188) q[3];
sx q[3];
rz(-1.3636953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.79080498) q[2];
sx q[2];
rz(-2.4420276) q[2];
sx q[2];
rz(-0.038185509) q[2];
rz(0.83032483) q[3];
sx q[3];
rz(-1.3874976) q[3];
sx q[3];
rz(3.0939046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4627948) q[0];
sx q[0];
rz(-2.9252453) q[0];
sx q[0];
rz(-1.3702673) q[0];
rz(-0.38617745) q[1];
sx q[1];
rz(-1.4048856) q[1];
sx q[1];
rz(2.4716299) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46028194) q[0];
sx q[0];
rz(-2.42762) q[0];
sx q[0];
rz(0.80912368) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.71670312) q[2];
sx q[2];
rz(-1.6316902) q[2];
sx q[2];
rz(-2.3961807) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.130387) q[1];
sx q[1];
rz(-0.97846675) q[1];
sx q[1];
rz(3.0051809) q[1];
rz(-pi) q[2];
rz(-1.7427722) q[3];
sx q[3];
rz(-0.21504083) q[3];
sx q[3];
rz(-0.57143927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2204444) q[2];
sx q[2];
rz(-0.28382742) q[2];
sx q[2];
rz(1.05668) q[2];
rz(-1.7250666) q[3];
sx q[3];
rz(-1.2456015) q[3];
sx q[3];
rz(1.2904803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0678299) q[0];
sx q[0];
rz(-1.0404328) q[0];
sx q[0];
rz(1.006806) q[0];
rz(-2.6489068) q[1];
sx q[1];
rz(-1.7308116) q[1];
sx q[1];
rz(0.83008343) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7573038) q[0];
sx q[0];
rz(-2.8606133) q[0];
sx q[0];
rz(3.021394) q[0];
rz(-0.30980457) q[2];
sx q[2];
rz(-2.5199157) q[2];
sx q[2];
rz(0.72337389) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.600297) q[1];
sx q[1];
rz(-1.4182555) q[1];
sx q[1];
rz(-2.9593854) q[1];
x q[2];
rz(-2.0138143) q[3];
sx q[3];
rz(-1.5823871) q[3];
sx q[3];
rz(-0.024758967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3229708) q[2];
sx q[2];
rz(-1.3016204) q[2];
sx q[2];
rz(-1.0969561) q[2];
rz(-1.2777626) q[3];
sx q[3];
rz(-2.4570229) q[3];
sx q[3];
rz(0.6440312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
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
rz(-0.11861435) q[1];
sx q[1];
rz(-1.4362486) q[1];
sx q[1];
rz(2.2172701) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7345006) q[0];
sx q[0];
rz(-1.0043) q[0];
sx q[0];
rz(1.9167711) q[0];
rz(-pi) q[1];
rz(2.4690041) q[2];
sx q[2];
rz(-1.158357) q[2];
sx q[2];
rz(-2.6809566) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9472447) q[1];
sx q[1];
rz(-1.2100265) q[1];
sx q[1];
rz(-1.5401349) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.64226867) q[3];
sx q[3];
rz(-0.54107252) q[3];
sx q[3];
rz(0.16995811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2529926) q[2];
sx q[2];
rz(-2.1414976) q[2];
sx q[2];
rz(-2.2484153) q[2];
rz(-1.1183974) q[3];
sx q[3];
rz(-1.018254) q[3];
sx q[3];
rz(-2.4848188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7669582) q[0];
sx q[0];
rz(-2.0603016) q[0];
sx q[0];
rz(1.7745071) q[0];
rz(-0.18192667) q[1];
sx q[1];
rz(-2.5839099) q[1];
sx q[1];
rz(-0.90669496) q[1];
rz(3.0497839) q[2];
sx q[2];
rz(-1.257007) q[2];
sx q[2];
rz(-1.9329482) q[2];
rz(-0.51707324) q[3];
sx q[3];
rz(-1.2096249) q[3];
sx q[3];
rz(0.67740868) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
