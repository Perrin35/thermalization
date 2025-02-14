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
rz(1.9396012) q[0];
sx q[0];
rz(-0.48299462) q[0];
sx q[0];
rz(-1.510409) q[0];
rz(-0.13934879) q[1];
sx q[1];
rz(-0.55958334) q[1];
sx q[1];
rz(-0.72996563) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8529634) q[0];
sx q[0];
rz(-1.4242111) q[0];
sx q[0];
rz(2.1781871) q[0];
rz(-pi) q[1];
rz(-0.59487409) q[2];
sx q[2];
rz(-1.383179) q[2];
sx q[2];
rz(-1.7930195) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0180359) q[1];
sx q[1];
rz(-0.58096051) q[1];
sx q[1];
rz(-0.058079795) q[1];
rz(-0.22878756) q[3];
sx q[3];
rz(-1.8513894) q[3];
sx q[3];
rz(2.5737263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0240747) q[2];
sx q[2];
rz(-2.8705609) q[2];
sx q[2];
rz(0.81895858) q[2];
rz(-0.0062746127) q[3];
sx q[3];
rz(-1.2333906) q[3];
sx q[3];
rz(0.98244572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5821424) q[0];
sx q[0];
rz(-1.5903951) q[0];
sx q[0];
rz(-2.5421802) q[0];
rz(-0.86743152) q[1];
sx q[1];
rz(-1.0299094) q[1];
sx q[1];
rz(2.3960466) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9052542) q[0];
sx q[0];
rz(-2.8239692) q[0];
sx q[0];
rz(2.7228217) q[0];
x q[1];
rz(-1.8944055) q[2];
sx q[2];
rz(-0.83728803) q[2];
sx q[2];
rz(0.74786438) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.47059575) q[1];
sx q[1];
rz(-1.7299011) q[1];
sx q[1];
rz(2.5039423) q[1];
x q[2];
rz(-0.87941951) q[3];
sx q[3];
rz(-1.9718277) q[3];
sx q[3];
rz(1.629231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6856689) q[2];
sx q[2];
rz(-0.29752877) q[2];
sx q[2];
rz(1.4073184) q[2];
rz(0.84434858) q[3];
sx q[3];
rz(-2.307939) q[3];
sx q[3];
rz(-1.0367397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5582964) q[0];
sx q[0];
rz(-1.946227) q[0];
sx q[0];
rz(-1.012828) q[0];
rz(0.60802513) q[1];
sx q[1];
rz(-1.5826179) q[1];
sx q[1];
rz(1.8720522) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0657261) q[0];
sx q[0];
rz(-2.009543) q[0];
sx q[0];
rz(0.42515305) q[0];
rz(-0.082398675) q[2];
sx q[2];
rz(-2.2401056) q[2];
sx q[2];
rz(-1.5562039) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3976058) q[1];
sx q[1];
rz(-2.3922046) q[1];
sx q[1];
rz(2.8515062) q[1];
rz(-0.49130398) q[3];
sx q[3];
rz(-1.2761371) q[3];
sx q[3];
rz(3.0719985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6262007) q[2];
sx q[2];
rz(-1.0582558) q[2];
sx q[2];
rz(-1.8499648) q[2];
rz(-0.0083222566) q[3];
sx q[3];
rz(-0.95299995) q[3];
sx q[3];
rz(2.9279809) q[3];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13852791) q[0];
sx q[0];
rz(-0.55532885) q[0];
sx q[0];
rz(1.7648765) q[0];
rz(1.5273013) q[1];
sx q[1];
rz(-1.9381899) q[1];
sx q[1];
rz(2.6208904) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4462858) q[0];
sx q[0];
rz(-0.55342544) q[0];
sx q[0];
rz(2.2880716) q[0];
rz(-2.7815656) q[2];
sx q[2];
rz(-2.0930778) q[2];
sx q[2];
rz(-1.479666) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.59078539) q[1];
sx q[1];
rz(-2.0001162) q[1];
sx q[1];
rz(-0.6145668) q[1];
rz(0.52526955) q[3];
sx q[3];
rz(-2.161918) q[3];
sx q[3];
rz(-2.1256465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.5248096) q[2];
sx q[2];
rz(-0.7889792) q[2];
sx q[2];
rz(1.3516124) q[2];
rz(1.0194408) q[3];
sx q[3];
rz(-2.5717058) q[3];
sx q[3];
rz(1.1714237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5919507) q[0];
sx q[0];
rz(-1.6787981) q[0];
sx q[0];
rz(-0.098966448) q[0];
rz(-2.4348266) q[1];
sx q[1];
rz(-0.87044972) q[1];
sx q[1];
rz(0.4695355) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5851368) q[0];
sx q[0];
rz(-0.13988189) q[0];
sx q[0];
rz(1.3993457) q[0];
rz(-2.7623873) q[2];
sx q[2];
rz(-1.0673041) q[2];
sx q[2];
rz(2.6852754) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8626735) q[1];
sx q[1];
rz(-1.9544744) q[1];
sx q[1];
rz(0.63273276) q[1];
rz(1.7170226) q[3];
sx q[3];
rz(-1.0195882) q[3];
sx q[3];
rz(2.220038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2609451) q[2];
sx q[2];
rz(-1.8069043) q[2];
sx q[2];
rz(1.8355231) q[2];
rz(0.4246873) q[3];
sx q[3];
rz(-1.3771907) q[3];
sx q[3];
rz(2.2556321) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0223618) q[0];
sx q[0];
rz(-0.62218085) q[0];
sx q[0];
rz(1.8192044) q[0];
rz(-1.127683) q[1];
sx q[1];
rz(-1.5195945) q[1];
sx q[1];
rz(1.3816396) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9499672) q[0];
sx q[0];
rz(-2.3268496) q[0];
sx q[0];
rz(-3.0227468) q[0];
rz(-2.1476949) q[2];
sx q[2];
rz(-1.111278) q[2];
sx q[2];
rz(1.3764868) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0387011) q[1];
sx q[1];
rz(-0.82566092) q[1];
sx q[1];
rz(1.8243415) q[1];
x q[2];
rz(2.840191) q[3];
sx q[3];
rz(-1.6414789) q[3];
sx q[3];
rz(1.7168728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7606925) q[2];
sx q[2];
rz(-1.5194632) q[2];
sx q[2];
rz(-2.6603928) q[2];
rz(2.6324658) q[3];
sx q[3];
rz(-2.9042518) q[3];
sx q[3];
rz(0.38207644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5612438) q[0];
sx q[0];
rz(-0.48155293) q[0];
sx q[0];
rz(3.0522108) q[0];
rz(1.0786169) q[1];
sx q[1];
rz(-1.6856472) q[1];
sx q[1];
rz(2.1481029) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.908113) q[0];
sx q[0];
rz(-0.78266875) q[0];
sx q[0];
rz(3.0872869) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3553576) q[2];
sx q[2];
rz(-1.0362384) q[2];
sx q[2];
rz(-1.3649808) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1878933) q[1];
sx q[1];
rz(-2.7479798) q[1];
sx q[1];
rz(-1.1792437) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5105623) q[3];
sx q[3];
rz(-1.1388121) q[3];
sx q[3];
rz(-0.61842266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3420458) q[2];
sx q[2];
rz(-2.5265103) q[2];
sx q[2];
rz(-0.40204027) q[2];
rz(-1.0953995) q[3];
sx q[3];
rz(-1.9270555) q[3];
sx q[3];
rz(-2.5636165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47356975) q[0];
sx q[0];
rz(-0.80900017) q[0];
sx q[0];
rz(-1.8079669) q[0];
rz(-2.2857621) q[1];
sx q[1];
rz(-1.4060833) q[1];
sx q[1];
rz(2.2241101) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79865341) q[0];
sx q[0];
rz(-0.22428939) q[0];
sx q[0];
rz(-2.0464315) q[0];
rz(0.11882527) q[2];
sx q[2];
rz(-0.95044327) q[2];
sx q[2];
rz(-2.9147749) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.92261295) q[1];
sx q[1];
rz(-0.41480468) q[1];
sx q[1];
rz(2.1475683) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0621655) q[3];
sx q[3];
rz(-0.98782238) q[3];
sx q[3];
rz(0.82750083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.63460073) q[2];
sx q[2];
rz(-1.4054106) q[2];
sx q[2];
rz(-1.290192) q[2];
rz(-0.90977943) q[3];
sx q[3];
rz(-1.6981643) q[3];
sx q[3];
rz(3.1006052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23583394) q[0];
sx q[0];
rz(-1.1474778) q[0];
sx q[0];
rz(2.8644417) q[0];
rz(1.9174891) q[1];
sx q[1];
rz(-1.6033019) q[1];
sx q[1];
rz(-1.7701497) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4215721) q[0];
sx q[0];
rz(-0.69584457) q[0];
sx q[0];
rz(2.8369342) q[0];
rz(-pi) q[1];
rz(-0.54746898) q[2];
sx q[2];
rz(-1.9488504) q[2];
sx q[2];
rz(0.34502703) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.692457) q[1];
sx q[1];
rz(-2.4583011) q[1];
sx q[1];
rz(1.9016674) q[1];
rz(-0.17454608) q[3];
sx q[3];
rz(-0.8564328) q[3];
sx q[3];
rz(2.3256231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.93854967) q[2];
sx q[2];
rz(-2.2800192) q[2];
sx q[2];
rz(-2.4533563) q[2];
rz(2.8271683) q[3];
sx q[3];
rz(-0.36948547) q[3];
sx q[3];
rz(-1.2615874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37453434) q[0];
sx q[0];
rz(-2.2735167) q[0];
sx q[0];
rz(2.5700997) q[0];
rz(2.4608965) q[1];
sx q[1];
rz(-1.9711767) q[1];
sx q[1];
rz(2.4933955) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2637973) q[0];
sx q[0];
rz(-1.5476942) q[0];
sx q[0];
rz(-0.013851555) q[0];
x q[1];
rz(2.2806281) q[2];
sx q[2];
rz(-2.2662244) q[2];
sx q[2];
rz(-0.26317138) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1614589) q[1];
sx q[1];
rz(-2.4422283) q[1];
sx q[1];
rz(3.1070263) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.33721029) q[3];
sx q[3];
rz(-0.54223947) q[3];
sx q[3];
rz(-0.48619871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.28006831) q[2];
sx q[2];
rz(-1.0736977) q[2];
sx q[2];
rz(-1.6746707) q[2];
rz(-2.9649949) q[3];
sx q[3];
rz(-2.4956775) q[3];
sx q[3];
rz(-1.8471898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27547729) q[0];
sx q[0];
rz(-1.3963516) q[0];
sx q[0];
rz(1.8657952) q[0];
rz(-0.38446174) q[1];
sx q[1];
rz(-1.6356331) q[1];
sx q[1];
rz(2.5148139) q[1];
rz(-2.7914417) q[2];
sx q[2];
rz(-1.2731009) q[2];
sx q[2];
rz(0.44558744) q[2];
rz(-2.6077059) q[3];
sx q[3];
rz(-2.0697099) q[3];
sx q[3];
rz(2.8768215) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
