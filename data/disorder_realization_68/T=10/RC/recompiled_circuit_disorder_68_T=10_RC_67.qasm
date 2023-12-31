OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.39419898) q[0];
sx q[0];
rz(2.6497901) q[0];
sx q[0];
rz(9.2368035) q[0];
rz(-1.1176874) q[1];
sx q[1];
rz(-1.517065) q[1];
sx q[1];
rz(2.7741073) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7029019) q[0];
sx q[0];
rz(-2.5950948) q[0];
sx q[0];
rz(-1.7706857) q[0];
rz(-pi) q[1];
rz(2.0619225) q[2];
sx q[2];
rz(-1.4516146) q[2];
sx q[2];
rz(-1.0246547) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.130598) q[1];
sx q[1];
rz(-1.8124609) q[1];
sx q[1];
rz(0.17250891) q[1];
rz(-0.34100451) q[3];
sx q[3];
rz(-1.4031646) q[3];
sx q[3];
rz(-2.117702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.964103) q[2];
sx q[2];
rz(-2.6314645) q[2];
sx q[2];
rz(-0.5509848) q[2];
rz(1.3059113) q[3];
sx q[3];
rz(-1.6492313) q[3];
sx q[3];
rz(1.8252385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47857639) q[0];
sx q[0];
rz(-2.1415648) q[0];
sx q[0];
rz(-2.6696894) q[0];
rz(2.7117803) q[1];
sx q[1];
rz(-1.2496354) q[1];
sx q[1];
rz(0.93634161) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9443003) q[0];
sx q[0];
rz(-0.24646491) q[0];
sx q[0];
rz(-2.775807) q[0];
x q[1];
rz(0.41994862) q[2];
sx q[2];
rz(-0.83050767) q[2];
sx q[2];
rz(-0.74479693) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.71948235) q[1];
sx q[1];
rz(-1.6750095) q[1];
sx q[1];
rz(2.4476493) q[1];
rz(-pi) q[2];
rz(-3.1182321) q[3];
sx q[3];
rz(-1.0682032) q[3];
sx q[3];
rz(2.392829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.77461809) q[2];
sx q[2];
rz(-2.8145511) q[2];
sx q[2];
rz(-2.7152087) q[2];
rz(1.2373699) q[3];
sx q[3];
rz(-2.5137413) q[3];
sx q[3];
rz(0.0330851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8957829) q[0];
sx q[0];
rz(-1.824546) q[0];
sx q[0];
rz(2.202503) q[0];
rz(2.242873) q[1];
sx q[1];
rz(-0.4788613) q[1];
sx q[1];
rz(-0.59392196) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9868601) q[0];
sx q[0];
rz(-1.8911456) q[0];
sx q[0];
rz(2.8264168) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5005641) q[2];
sx q[2];
rz(-0.68968455) q[2];
sx q[2];
rz(0.94674142) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.18332874) q[1];
sx q[1];
rz(-1.7839583) q[1];
sx q[1];
rz(-2.0002685) q[1];
rz(0.61641683) q[3];
sx q[3];
rz(-0.55509242) q[3];
sx q[3];
rz(-0.68716955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.64017355) q[2];
sx q[2];
rz(-0.72945014) q[2];
sx q[2];
rz(1.7017986) q[2];
rz(-0.38763186) q[3];
sx q[3];
rz(-1.6250316) q[3];
sx q[3];
rz(-2.7534527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7664292) q[0];
sx q[0];
rz(-1.5901934) q[0];
sx q[0];
rz(0.50278062) q[0];
rz(0.76820961) q[1];
sx q[1];
rz(-0.50351024) q[1];
sx q[1];
rz(0.75685135) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4810281) q[0];
sx q[0];
rz(-1.2111944) q[0];
sx q[0];
rz(-1.9268376) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1050622) q[2];
sx q[2];
rz(-1.5823936) q[2];
sx q[2];
rz(-2.2955745) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.13285747) q[1];
sx q[1];
rz(-2.1454304) q[1];
sx q[1];
rz(2.7112714) q[1];
rz(-pi) q[2];
rz(-2.5960856) q[3];
sx q[3];
rz(-2.3295998) q[3];
sx q[3];
rz(3.0363887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7148774) q[2];
sx q[2];
rz(-1.2299512) q[2];
sx q[2];
rz(1.654401) q[2];
rz(0.58250827) q[3];
sx q[3];
rz(-2.0472066) q[3];
sx q[3];
rz(2.5845161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8476167) q[0];
sx q[0];
rz(-1.0852381) q[0];
sx q[0];
rz(0.75772444) q[0];
rz(-1.853653) q[1];
sx q[1];
rz(-2.2133591) q[1];
sx q[1];
rz(1.0505189) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30848635) q[0];
sx q[0];
rz(-1.8252488) q[0];
sx q[0];
rz(1.893505) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6987223) q[2];
sx q[2];
rz(-1.0204698) q[2];
sx q[2];
rz(1.0843104) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1861021) q[1];
sx q[1];
rz(-2.0932066) q[1];
sx q[1];
rz(2.2239457) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6090622) q[3];
sx q[3];
rz(-0.88056394) q[3];
sx q[3];
rz(-2.0257476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.918255) q[2];
sx q[2];
rz(-0.35327062) q[2];
sx q[2];
rz(2.5081432) q[2];
rz(-1.9472306) q[3];
sx q[3];
rz(-1.4712237) q[3];
sx q[3];
rz(0.71715322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.441992) q[0];
sx q[0];
rz(-0.49015912) q[0];
sx q[0];
rz(-2.8884086) q[0];
rz(1.6075915) q[1];
sx q[1];
rz(-1.7065159) q[1];
sx q[1];
rz(1.6794499) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77308649) q[0];
sx q[0];
rz(-0.22531548) q[0];
sx q[0];
rz(2.7789475) q[0];
rz(2.2199549) q[2];
sx q[2];
rz(-1.7643133) q[2];
sx q[2];
rz(-0.18093872) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6219382) q[1];
sx q[1];
rz(-0.42236537) q[1];
sx q[1];
rz(2.7154891) q[1];
x q[2];
rz(0.31520321) q[3];
sx q[3];
rz(-1.8990371) q[3];
sx q[3];
rz(2.4344276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9399461) q[2];
sx q[2];
rz(-2.3942409) q[2];
sx q[2];
rz(0.80491006) q[2];
rz(-1.1770052) q[3];
sx q[3];
rz(-2.1046808) q[3];
sx q[3];
rz(-1.0837519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8686304) q[0];
sx q[0];
rz(-1.0694163) q[0];
sx q[0];
rz(2.2139363) q[0];
rz(1.0246798) q[1];
sx q[1];
rz(-1.506348) q[1];
sx q[1];
rz(-2.129508) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0976809) q[0];
sx q[0];
rz(-2.0872396) q[0];
sx q[0];
rz(-2.5994356) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7094678) q[2];
sx q[2];
rz(-1.9562634) q[2];
sx q[2];
rz(3.1388381) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.183179) q[1];
sx q[1];
rz(-1.5378386) q[1];
sx q[1];
rz(0.54380137) q[1];
x q[2];
rz(-1.4885694) q[3];
sx q[3];
rz(-1.0136908) q[3];
sx q[3];
rz(2.0092056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.53081375) q[2];
sx q[2];
rz(-1.6616219) q[2];
sx q[2];
rz(2.7116595) q[2];
rz(2.1271465) q[3];
sx q[3];
rz(-2.7323664) q[3];
sx q[3];
rz(-0.53340069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.6124509) q[0];
sx q[0];
rz(-2.2021459) q[0];
sx q[0];
rz(0.21417831) q[0];
rz(-1.051349) q[1];
sx q[1];
rz(-2.9290757) q[1];
sx q[1];
rz(-2.8578551) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8026233) q[0];
sx q[0];
rz(-1.6049275) q[0];
sx q[0];
rz(0.3011093) q[0];
rz(-pi) q[1];
rz(2.0869414) q[2];
sx q[2];
rz(-0.37545855) q[2];
sx q[2];
rz(1.6106538) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9323862) q[1];
sx q[1];
rz(-1.3750409) q[1];
sx q[1];
rz(2.6372361) q[1];
rz(-pi) q[2];
x q[2];
rz(0.65225668) q[3];
sx q[3];
rz(-1.0271003) q[3];
sx q[3];
rz(0.039747681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6909137) q[2];
sx q[2];
rz(-0.51270715) q[2];
sx q[2];
rz(-1.4452176) q[2];
rz(1.5971659) q[3];
sx q[3];
rz(-1.4404567) q[3];
sx q[3];
rz(-0.33932313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.504869) q[0];
sx q[0];
rz(-3.0915785) q[0];
sx q[0];
rz(-0.069256393) q[0];
rz(-1.6537369) q[1];
sx q[1];
rz(-1.8585049) q[1];
sx q[1];
rz(-1.5725296) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78281392) q[0];
sx q[0];
rz(-0.94952119) q[0];
sx q[0];
rz(2.5675979) q[0];
rz(2.9022129) q[2];
sx q[2];
rz(-2.8068672) q[2];
sx q[2];
rz(-0.3604381) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3410586) q[1];
sx q[1];
rz(-1.8039928) q[1];
sx q[1];
rz(1.484616) q[1];
rz(1.2980372) q[3];
sx q[3];
rz(-2.4747304) q[3];
sx q[3];
rz(2.42815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9562324) q[2];
sx q[2];
rz(-2.9113443) q[2];
sx q[2];
rz(-3.0017079) q[2];
rz(2.774003) q[3];
sx q[3];
rz(-1.9544173) q[3];
sx q[3];
rz(-2.1504413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1763828) q[0];
sx q[0];
rz(-0.39127025) q[0];
sx q[0];
rz(-0.64176732) q[0];
rz(1.2311252) q[1];
sx q[1];
rz(-1.9893913) q[1];
sx q[1];
rz(0.26783255) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7627015) q[0];
sx q[0];
rz(-2.2305616) q[0];
sx q[0];
rz(2.1428109) q[0];
rz(-1.7540625) q[2];
sx q[2];
rz(-1.6943309) q[2];
sx q[2];
rz(-1.087041) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.99991998) q[1];
sx q[1];
rz(-1.8598286) q[1];
sx q[1];
rz(2.8029867) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0249428) q[3];
sx q[3];
rz(-0.59026679) q[3];
sx q[3];
rz(-2.2772307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3315167) q[2];
sx q[2];
rz(-1.6167567) q[2];
sx q[2];
rz(-1.4769185) q[2];
rz(-0.26633513) q[3];
sx q[3];
rz(-0.24644066) q[3];
sx q[3];
rz(2.5951071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1289566) q[0];
sx q[0];
rz(-0.92100443) q[0];
sx q[0];
rz(-2.2367649) q[0];
rz(-0.77990445) q[1];
sx q[1];
rz(-0.48702469) q[1];
sx q[1];
rz(-1.3866966) q[1];
rz(-0.98942479) q[2];
sx q[2];
rz(-2.1944254) q[2];
sx q[2];
rz(2.2133322) q[2];
rz(-2.4640502) q[3];
sx q[3];
rz(-1.4031706) q[3];
sx q[3];
rz(-1.8085898) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
