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
rz(-2.1674018) q[0];
sx q[0];
rz(-0.84889698) q[0];
sx q[0];
rz(1.4155686) q[0];
rz(2.2805136) q[1];
sx q[1];
rz(-1.5134892) q[1];
sx q[1];
rz(2.3754062) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9250037) q[0];
sx q[0];
rz(-1.6356409) q[0];
sx q[0];
rz(1.5349887) q[0];
x q[1];
rz(0.033896565) q[2];
sx q[2];
rz(-1.6914512) q[2];
sx q[2];
rz(1.1519037) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.79957818) q[1];
sx q[1];
rz(-1.7740846) q[1];
sx q[1];
rz(-1.5631097) q[1];
rz(-pi) q[2];
rz(-1.1580519) q[3];
sx q[3];
rz(-1.6696602) q[3];
sx q[3];
rz(-3.1237941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8908995) q[2];
sx q[2];
rz(-2.3919899) q[2];
sx q[2];
rz(3.0787943) q[2];
rz(1.4248258) q[3];
sx q[3];
rz(-1.3751605) q[3];
sx q[3];
rz(-2.4019901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0733136) q[0];
sx q[0];
rz(-2.8977019) q[0];
sx q[0];
rz(1.1163611) q[0];
rz(1.5721389) q[1];
sx q[1];
rz(-2.8469323) q[1];
sx q[1];
rz(-2.0223298) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37518445) q[0];
sx q[0];
rz(-1.5699727) q[0];
sx q[0];
rz(-2.4341325) q[0];
x q[1];
rz(2.767973) q[2];
sx q[2];
rz(-1.116215) q[2];
sx q[2];
rz(1.0475243) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5183603) q[1];
sx q[1];
rz(-2.6122852) q[1];
sx q[1];
rz(0.12959403) q[1];
x q[2];
rz(-2.4517566) q[3];
sx q[3];
rz(-2.3202826) q[3];
sx q[3];
rz(2.6765055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0548627) q[2];
sx q[2];
rz(-0.77069288) q[2];
sx q[2];
rz(-1.9697624) q[2];
rz(2.150676) q[3];
sx q[3];
rz(-1.0057697) q[3];
sx q[3];
rz(-1.2539554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1427926) q[0];
sx q[0];
rz(-0.87738377) q[0];
sx q[0];
rz(-0.073609322) q[0];
rz(-0.40120801) q[1];
sx q[1];
rz(-0.3708655) q[1];
sx q[1];
rz(2.4295095) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82862604) q[0];
sx q[0];
rz(-3.0168128) q[0];
sx q[0];
rz(1.9931884) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4855494) q[2];
sx q[2];
rz(-1.5696161) q[2];
sx q[2];
rz(-0.16108433) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0910511) q[1];
sx q[1];
rz(-1.8767722) q[1];
sx q[1];
rz(-1.6832965) q[1];
rz(-pi) q[2];
rz(-2.1041342) q[3];
sx q[3];
rz(-1.6210307) q[3];
sx q[3];
rz(-2.2118357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3965343) q[2];
sx q[2];
rz(-1.859954) q[2];
sx q[2];
rz(0.38953951) q[2];
rz(-0.39858308) q[3];
sx q[3];
rz(-1.4094718) q[3];
sx q[3];
rz(1.7194974) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013537708) q[0];
sx q[0];
rz(-0.062352926) q[0];
sx q[0];
rz(2.7259977) q[0];
rz(-1.2526814) q[1];
sx q[1];
rz(-2.0038192) q[1];
sx q[1];
rz(2.2223053) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.141826) q[0];
sx q[0];
rz(-1.6117038) q[0];
sx q[0];
rz(1.7342939) q[0];
rz(-1.5482727) q[2];
sx q[2];
rz(-1.833758) q[2];
sx q[2];
rz(-2.7561273) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.63863149) q[1];
sx q[1];
rz(-0.98010862) q[1];
sx q[1];
rz(-0.69815127) q[1];
x q[2];
rz(2.0218798) q[3];
sx q[3];
rz(-0.85441426) q[3];
sx q[3];
rz(0.35600397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2933942) q[2];
sx q[2];
rz(-1.5092809) q[2];
sx q[2];
rz(-0.97818127) q[2];
rz(0.50655043) q[3];
sx q[3];
rz(-1.4975486) q[3];
sx q[3];
rz(-2.1250561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.070234805) q[0];
sx q[0];
rz(-2.4429584) q[0];
sx q[0];
rz(-0.51682669) q[0];
rz(1.0454073) q[1];
sx q[1];
rz(-1.005859) q[1];
sx q[1];
rz(-2.7396835) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4025426) q[0];
sx q[0];
rz(-2.2687904) q[0];
sx q[0];
rz(-0.9902467) q[0];
rz(2.972225) q[2];
sx q[2];
rz(-1.8529839) q[2];
sx q[2];
rz(2.3991632) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4116844) q[1];
sx q[1];
rz(-2.9422816) q[1];
sx q[1];
rz(1.5491376) q[1];
rz(-pi) q[2];
rz(-3.0203781) q[3];
sx q[3];
rz(-2.182721) q[3];
sx q[3];
rz(-0.90288375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7261293) q[2];
sx q[2];
rz(-1.1795421) q[2];
sx q[2];
rz(2.9187875) q[2];
rz(2.2807617) q[3];
sx q[3];
rz(-1.6313044) q[3];
sx q[3];
rz(-1.8676602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083542682) q[0];
sx q[0];
rz(-1.8996176) q[0];
sx q[0];
rz(-2.2831931) q[0];
rz(-2.3421613) q[1];
sx q[1];
rz(-1.0170499) q[1];
sx q[1];
rz(-1.7564836) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1027062) q[0];
sx q[0];
rz(-1.8441851) q[0];
sx q[0];
rz(0.1929184) q[0];
rz(-pi) q[1];
rz(0.054692312) q[2];
sx q[2];
rz(-1.8963976) q[2];
sx q[2];
rz(-0.97132746) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.68011682) q[1];
sx q[1];
rz(-1.450074) q[1];
sx q[1];
rz(-2.5402447) q[1];
rz(1.0812182) q[3];
sx q[3];
rz(-1.7822232) q[3];
sx q[3];
rz(-1.2871773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4110306) q[2];
sx q[2];
rz(-0.73341122) q[2];
sx q[2];
rz(-0.15156558) q[2];
rz(-1.1847121) q[3];
sx q[3];
rz(-1.1368753) q[3];
sx q[3];
rz(-0.90004164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2043692) q[0];
sx q[0];
rz(-1.3441514) q[0];
sx q[0];
rz(2.3580661) q[0];
rz(-0.20768684) q[1];
sx q[1];
rz(-2.2238104) q[1];
sx q[1];
rz(1.9327365) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0328114) q[0];
sx q[0];
rz(-2.0447123) q[0];
sx q[0];
rz(1.117871) q[0];
x q[1];
rz(-2.1561106) q[2];
sx q[2];
rz(-0.83067465) q[2];
sx q[2];
rz(-2.5565846) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.085371278) q[1];
sx q[1];
rz(-1.610145) q[1];
sx q[1];
rz(0.84169547) q[1];
rz(-2.9840488) q[3];
sx q[3];
rz(-1.9374401) q[3];
sx q[3];
rz(3.1160924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8565159) q[2];
sx q[2];
rz(-0.74399844) q[2];
sx q[2];
rz(-3.0169955) q[2];
rz(1.3397269) q[3];
sx q[3];
rz(-2.7281269) q[3];
sx q[3];
rz(1.9477897) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1272005) q[0];
sx q[0];
rz(-1.7008282) q[0];
sx q[0];
rz(0.097323962) q[0];
rz(2.4864181) q[1];
sx q[1];
rz(-0.56071463) q[1];
sx q[1];
rz(-2.4526144) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1611947) q[0];
sx q[0];
rz(-1.4965048) q[0];
sx q[0];
rz(1.0832537) q[0];
x q[1];
rz(1.5083049) q[2];
sx q[2];
rz(-0.35258503) q[2];
sx q[2];
rz(3.0931713) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2804451) q[1];
sx q[1];
rz(-0.89830924) q[1];
sx q[1];
rz(1.4073744) q[1];
rz(-pi) q[2];
rz(-0.68920691) q[3];
sx q[3];
rz(-1.6968932) q[3];
sx q[3];
rz(1.4339303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9002306) q[2];
sx q[2];
rz(-0.9598383) q[2];
sx q[2];
rz(-2.0755365) q[2];
rz(-2.7212932) q[3];
sx q[3];
rz(-1.2371141) q[3];
sx q[3];
rz(-1.6828116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083753839) q[0];
sx q[0];
rz(-0.76863113) q[0];
sx q[0];
rz(2.2572857) q[0];
rz(0.090944313) q[1];
sx q[1];
rz(-0.93946409) q[1];
sx q[1];
rz(-1.1788064) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5835253) q[0];
sx q[0];
rz(-1.5060165) q[0];
sx q[0];
rz(0.16048783) q[0];
x q[1];
rz(-0.0021807533) q[2];
sx q[2];
rz(-0.90490393) q[2];
sx q[2];
rz(-1.3058174) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7024421) q[1];
sx q[1];
rz(-1.0156173) q[1];
sx q[1];
rz(1.7722284) q[1];
x q[2];
rz(0.20060819) q[3];
sx q[3];
rz(-1.6880013) q[3];
sx q[3];
rz(2.8497118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.74943298) q[2];
sx q[2];
rz(-1.9955111) q[2];
sx q[2];
rz(2.667099) q[2];
rz(-0.042304603) q[3];
sx q[3];
rz(-1.5174805) q[3];
sx q[3];
rz(-1.8208985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6831191) q[0];
sx q[0];
rz(-1.1011769) q[0];
sx q[0];
rz(1.839757) q[0];
rz(-0.63109541) q[1];
sx q[1];
rz(-1.0429691) q[1];
sx q[1];
rz(-2.1249117) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23970579) q[0];
sx q[0];
rz(-1.282987) q[0];
sx q[0];
rz(0.91602703) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.053875845) q[2];
sx q[2];
rz(-1.4720671) q[2];
sx q[2];
rz(-1.3786664) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.8606372) q[1];
sx q[1];
rz(-1.5796575) q[1];
sx q[1];
rz(-2.289516) q[1];
rz(-pi) q[2];
x q[2];
rz(2.437463) q[3];
sx q[3];
rz(-1.0350772) q[3];
sx q[3];
rz(-2.2975722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3116698) q[2];
sx q[2];
rz(-2.416478) q[2];
sx q[2];
rz(2.1666918) q[2];
rz(-1.4026862) q[3];
sx q[3];
rz(-0.97167531) q[3];
sx q[3];
rz(1.5909125) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1350031) q[0];
sx q[0];
rz(-2.4000744) q[0];
sx q[0];
rz(0.97493521) q[0];
rz(-1.7368175) q[1];
sx q[1];
rz(-1.4092696) q[1];
sx q[1];
rz(-0.41913941) q[1];
rz(2.1640833) q[2];
sx q[2];
rz(-1.1326264) q[2];
sx q[2];
rz(-2.626194) q[2];
rz(1.4584803) q[3];
sx q[3];
rz(-1.9989804) q[3];
sx q[3];
rz(1.2811666) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
