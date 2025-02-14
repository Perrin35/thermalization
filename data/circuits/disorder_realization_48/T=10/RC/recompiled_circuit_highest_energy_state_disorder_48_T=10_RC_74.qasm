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
rz(2.2639182) q[1];
sx q[1];
rz(-1.9281153) q[1];
sx q[1];
rz(0.47668639) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3221368) q[0];
sx q[0];
rz(-1.9662204) q[0];
sx q[0];
rz(2.8614392) q[0];
x q[1];
rz(-2.0332912) q[2];
sx q[2];
rz(-1.1574189) q[2];
sx q[2];
rz(-0.55962901) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4922766) q[1];
sx q[1];
rz(-1.5879803) q[1];
sx q[1];
rz(-1.3130929) q[1];
rz(-pi) q[2];
x q[2];
rz(0.4313978) q[3];
sx q[3];
rz(-1.3702624) q[3];
sx q[3];
rz(1.7748743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1276663) q[2];
sx q[2];
rz(-1.4907336) q[2];
sx q[2];
rz(-0.15065436) q[2];
rz(1.6023747) q[3];
sx q[3];
rz(-2.7645002) q[3];
sx q[3];
rz(-2.8902875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0852614) q[0];
sx q[0];
rz(-1.3717317) q[0];
sx q[0];
rz(-0.37013176) q[0];
rz(-0.47259694) q[1];
sx q[1];
rz(-2.8234146) q[1];
sx q[1];
rz(-0.95742375) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2070302) q[0];
sx q[0];
rz(-1.0732722) q[0];
sx q[0];
rz(3.050404) q[0];
rz(-1.2217058) q[2];
sx q[2];
rz(-0.47609509) q[2];
sx q[2];
rz(0.27175909) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.606876) q[1];
sx q[1];
rz(-0.83321111) q[1];
sx q[1];
rz(1.4556769) q[1];
rz(0.40471802) q[3];
sx q[3];
rz(-2.4232695) q[3];
sx q[3];
rz(-0.61678929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.0051673278) q[2];
sx q[2];
rz(-0.74446669) q[2];
sx q[2];
rz(-2.1180507) q[2];
rz(1.3905904) q[3];
sx q[3];
rz(-1.3955045) q[3];
sx q[3];
rz(1.3172147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7904974) q[0];
sx q[0];
rz(-1.0095162) q[0];
sx q[0];
rz(-2.5787831) q[0];
rz(-1.5142745) q[1];
sx q[1];
rz(-1.4311675) q[1];
sx q[1];
rz(-3.0624342) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25185967) q[0];
sx q[0];
rz(-1.8859698) q[0];
sx q[0];
rz(-0.51879518) q[0];
rz(-pi) q[1];
rz(-1.5299893) q[2];
sx q[2];
rz(-1.0804324) q[2];
sx q[2];
rz(-2.9029569) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9197977) q[1];
sx q[1];
rz(-2.256003) q[1];
sx q[1];
rz(2.5849839) q[1];
x q[2];
rz(-2.4412254) q[3];
sx q[3];
rz(-1.0031962) q[3];
sx q[3];
rz(1.2705453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9074273) q[2];
sx q[2];
rz(-1.558446) q[2];
sx q[2];
rz(1.8094212) q[2];
rz(-2.7857156) q[3];
sx q[3];
rz(-1.9477113) q[3];
sx q[3];
rz(-0.98085421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(-2.6818162) q[0];
sx q[0];
rz(-1.9720607) q[0];
sx q[0];
rz(2.0303149) q[0];
rz(0.92012826) q[1];
sx q[1];
rz(-1.5524813) q[1];
sx q[1];
rz(0.2535325) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4944015) q[0];
sx q[0];
rz(-1.8888942) q[0];
sx q[0];
rz(-2.6107016) q[0];
x q[1];
rz(0.61863135) q[2];
sx q[2];
rz(-2.9455373) q[2];
sx q[2];
rz(0.679099) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.580198) q[1];
sx q[1];
rz(-0.78402482) q[1];
sx q[1];
rz(2.8785588) q[1];
rz(2.8315968) q[3];
sx q[3];
rz(-1.0862203) q[3];
sx q[3];
rz(2.692401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.8186875) q[2];
sx q[2];
rz(-2.5762168) q[2];
sx q[2];
rz(1.6820071) q[2];
rz(-0.5736351) q[3];
sx q[3];
rz(-1.2021659) q[3];
sx q[3];
rz(-0.044376317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6673073) q[0];
sx q[0];
rz(-1.0695589) q[0];
sx q[0];
rz(-1.7942418) q[0];
rz(-2.5509293) q[1];
sx q[1];
rz(-1.9213516) q[1];
sx q[1];
rz(1.7151054) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7280557) q[0];
sx q[0];
rz(-1.1794834) q[0];
sx q[0];
rz(1.1831468) q[0];
x q[1];
rz(2.6408004) q[2];
sx q[2];
rz(-2.4363849) q[2];
sx q[2];
rz(2.1057745) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.062078309) q[1];
sx q[1];
rz(-0.2582363) q[1];
sx q[1];
rz(-2.9596427) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9699202) q[3];
sx q[3];
rz(-0.35997691) q[3];
sx q[3];
rz(-1.5300919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.23919375) q[2];
sx q[2];
rz(-1.0603443) q[2];
sx q[2];
rz(-1.9690465) q[2];
rz(-0.086056195) q[3];
sx q[3];
rz(-0.9525758) q[3];
sx q[3];
rz(-0.15611592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0657144) q[0];
sx q[0];
rz(-1.5944163) q[0];
sx q[0];
rz(-0.86268798) q[0];
rz(-0.24738303) q[1];
sx q[1];
rz(-1.3037953) q[1];
sx q[1];
rz(0.47201306) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2675572) q[0];
sx q[0];
rz(-2.3568601) q[0];
sx q[0];
rz(0.48798706) q[0];
rz(-1.3026779) q[2];
sx q[2];
rz(-2.6202218) q[2];
sx q[2];
rz(0.15586317) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.029376349) q[1];
sx q[1];
rz(-2.7885155) q[1];
sx q[1];
rz(2.608271) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8627164) q[3];
sx q[3];
rz(-2.4987767) q[3];
sx q[3];
rz(2.0187279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.90999675) q[2];
sx q[2];
rz(-0.96698499) q[2];
sx q[2];
rz(1.7426573) q[2];
rz(-1.4553962) q[3];
sx q[3];
rz(-0.78794909) q[3];
sx q[3];
rz(-0.56070352) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6550734) q[0];
sx q[0];
rz(-1.1460679) q[0];
sx q[0];
rz(-2.6229677) q[0];
rz(-2.4229166) q[1];
sx q[1];
rz(-0.51281723) q[1];
sx q[1];
rz(-1.3291043) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3603921) q[0];
sx q[0];
rz(-2.5787163) q[0];
sx q[0];
rz(-0.84873523) q[0];
rz(-pi) q[1];
rz(3.091264) q[2];
sx q[2];
rz(-1.8476881) q[2];
sx q[2];
rz(0.30093873) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.74187255) q[1];
sx q[1];
rz(-2.406027) q[1];
sx q[1];
rz(-1.5217785) q[1];
x q[2];
rz(1.6006599) q[3];
sx q[3];
rz(-2.5549508) q[3];
sx q[3];
rz(1.7778974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.79080498) q[2];
sx q[2];
rz(-0.69956508) q[2];
sx q[2];
rz(-3.1034071) q[2];
rz(-2.3112678) q[3];
sx q[3];
rz(-1.3874976) q[3];
sx q[3];
rz(3.0939046) q[3];
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
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6787978) q[0];
sx q[0];
rz(-0.21634732) q[0];
sx q[0];
rz(-1.7713254) q[0];
rz(2.7554152) q[1];
sx q[1];
rz(-1.4048856) q[1];
sx q[1];
rz(2.4716299) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46028194) q[0];
sx q[0];
rz(-2.42762) q[0];
sx q[0];
rz(-0.80912368) q[0];
rz(2.4248895) q[2];
sx q[2];
rz(-1.6316902) q[2];
sx q[2];
rz(0.74541192) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0112056) q[1];
sx q[1];
rz(-0.97846675) q[1];
sx q[1];
rz(0.13641178) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3588328) q[3];
sx q[3];
rz(-1.6073213) q[3];
sx q[3];
rz(1.9741457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.92114821) q[2];
sx q[2];
rz(-0.28382742) q[2];
sx q[2];
rz(1.05668) q[2];
rz(-1.4165261) q[3];
sx q[3];
rz(-1.2456015) q[3];
sx q[3];
rz(-1.2904803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
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
rz(-2.6489068) q[1];
sx q[1];
rz(-1.7308116) q[1];
sx q[1];
rz(-2.3115092) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5093436) q[0];
sx q[0];
rz(-1.8496939) q[0];
sx q[0];
rz(-1.6053902) q[0];
rz(-pi) q[1];
rz(-0.5987723) q[2];
sx q[2];
rz(-1.7493003) q[2];
sx q[2];
rz(2.0395961) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.3398847) q[1];
sx q[1];
rz(-2.9045068) q[1];
sx q[1];
rz(2.4380142) q[1];
x q[2];
rz(-3.1287635) q[3];
sx q[3];
rz(-1.1278102) q[3];
sx q[3];
rz(1.6010546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.81862187) q[2];
sx q[2];
rz(-1.3016204) q[2];
sx q[2];
rz(-2.0446365) q[2];
rz(1.2777626) q[3];
sx q[3];
rz(-0.68456972) q[3];
sx q[3];
rz(-2.4975615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85035664) q[0];
sx q[0];
rz(-2.0391897) q[0];
sx q[0];
rz(-0.37503234) q[0];
rz(-0.11861435) q[1];
sx q[1];
rz(-1.4362486) q[1];
sx q[1];
rz(2.2172701) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3548022) q[0];
sx q[0];
rz(-1.2806007) q[0];
sx q[0];
rz(2.5470887) q[0];
rz(1.0607988) q[2];
sx q[2];
rz(-0.96335232) q[2];
sx q[2];
rz(-2.3404672) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0339151) q[1];
sx q[1];
rz(-0.36201358) q[1];
sx q[1];
rz(0.081078366) q[1];
rz(-pi) q[2];
rz(-0.64226867) q[3];
sx q[3];
rz(-2.6005201) q[3];
sx q[3];
rz(2.9716345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.88860005) q[2];
sx q[2];
rz(-2.1414976) q[2];
sx q[2];
rz(-0.89317733) q[2];
rz(-2.0231953) q[3];
sx q[3];
rz(-2.1233386) q[3];
sx q[3];
rz(-2.4848188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7669582) q[0];
sx q[0];
rz(-2.0603016) q[0];
sx q[0];
rz(1.7745071) q[0];
rz(0.18192667) q[1];
sx q[1];
rz(-0.55768273) q[1];
sx q[1];
rz(2.2348977) q[1];
rz(0.091808783) q[2];
sx q[2];
rz(-1.8845857) q[2];
sx q[2];
rz(1.2086445) q[2];
rz(-2.4891067) q[3];
sx q[3];
rz(-2.5204044) q[3];
sx q[3];
rz(1.6922097) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
