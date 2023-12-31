OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.49139872) q[0];
sx q[0];
rz(2.8770652) q[0];
sx q[0];
rz(9.8192083) q[0];
rz(0.0061622942) q[1];
sx q[1];
rz(-0.34024629) q[1];
sx q[1];
rz(1.9415829) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8631247) q[0];
sx q[0];
rz(-1.5651363) q[0];
sx q[0];
rz(1.5229043) q[0];
rz(-pi) q[1];
rz(-2.1985487) q[2];
sx q[2];
rz(-0.53484166) q[2];
sx q[2];
rz(0.56378555) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.29772273) q[1];
sx q[1];
rz(-0.6134609) q[1];
sx q[1];
rz(-0.86617275) q[1];
rz(-pi) q[2];
rz(-0.37813152) q[3];
sx q[3];
rz(-1.6117125) q[3];
sx q[3];
rz(3.0331628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8346943) q[2];
sx q[2];
rz(-1.6480185) q[2];
sx q[2];
rz(-2.8519894) q[2];
rz(-0.87537193) q[3];
sx q[3];
rz(-2.1353728) q[3];
sx q[3];
rz(-3.0818821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5263379) q[0];
sx q[0];
rz(-2.2741788) q[0];
sx q[0];
rz(0.34399024) q[0];
rz(-0.084331766) q[1];
sx q[1];
rz(-0.66939676) q[1];
sx q[1];
rz(1.7864236) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2750759) q[0];
sx q[0];
rz(-1.6356902) q[0];
sx q[0];
rz(2.1823723) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.64135735) q[2];
sx q[2];
rz(-0.20198447) q[2];
sx q[2];
rz(-0.55822492) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5835186) q[1];
sx q[1];
rz(-2.0194224) q[1];
sx q[1];
rz(-3.0706057) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7543206) q[3];
sx q[3];
rz(-1.4944544) q[3];
sx q[3];
rz(-0.21218382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8460059) q[2];
sx q[2];
rz(-2.3082374) q[2];
sx q[2];
rz(-2.6039092) q[2];
rz(2.6387571) q[3];
sx q[3];
rz(-0.42481315) q[3];
sx q[3];
rz(-2.1311549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66353345) q[0];
sx q[0];
rz(-0.53776598) q[0];
sx q[0];
rz(2.5007201) q[0];
rz(-2.3928941) q[1];
sx q[1];
rz(-1.0083895) q[1];
sx q[1];
rz(-2.0764988) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6247892) q[0];
sx q[0];
rz(-1.5771126) q[0];
sx q[0];
rz(-1.2114695) q[0];
rz(-pi) q[1];
rz(2.2674019) q[2];
sx q[2];
rz(-2.0421931) q[2];
sx q[2];
rz(-0.092560571) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.47485891) q[1];
sx q[1];
rz(-1.6068646) q[1];
sx q[1];
rz(-0.95400793) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2933526) q[3];
sx q[3];
rz(-1.3812997) q[3];
sx q[3];
rz(3.0998067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.37725267) q[2];
sx q[2];
rz(-1.9621907) q[2];
sx q[2];
rz(-2.4242145) q[2];
rz(2.453089) q[3];
sx q[3];
rz(-2.513956) q[3];
sx q[3];
rz(-0.29754105) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82729572) q[0];
sx q[0];
rz(-1.941444) q[0];
sx q[0];
rz(-0.24969077) q[0];
rz(-1.0149792) q[1];
sx q[1];
rz(-0.29622886) q[1];
sx q[1];
rz(0.011118523) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.714689) q[0];
sx q[0];
rz(-2.5137797) q[0];
sx q[0];
rz(0.83616242) q[0];
rz(-pi) q[1];
rz(-0.83799329) q[2];
sx q[2];
rz(-1.0655155) q[2];
sx q[2];
rz(-2.879564) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0032469) q[1];
sx q[1];
rz(-2.0577288) q[1];
sx q[1];
rz(-0.37133118) q[1];
x q[2];
rz(1.7632145) q[3];
sx q[3];
rz(-1.8835861) q[3];
sx q[3];
rz(2.6421412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.49542385) q[2];
sx q[2];
rz(-1.8858706) q[2];
sx q[2];
rz(-0.28309506) q[2];
rz(-0.66343534) q[3];
sx q[3];
rz(-2.5484271) q[3];
sx q[3];
rz(-2.2535113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11113142) q[0];
sx q[0];
rz(-0.24065329) q[0];
sx q[0];
rz(0.33183137) q[0];
rz(2.6470673) q[1];
sx q[1];
rz(-1.8437513) q[1];
sx q[1];
rz(1.3269075) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3834284) q[0];
sx q[0];
rz(-2.313662) q[0];
sx q[0];
rz(1.3616256) q[0];
rz(-pi) q[1];
rz(-0.97557108) q[2];
sx q[2];
rz(-2.274548) q[2];
sx q[2];
rz(3.1095568) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2697619) q[1];
sx q[1];
rz(-0.40553906) q[1];
sx q[1];
rz(0.62015066) q[1];
rz(0.29019659) q[3];
sx q[3];
rz(-1.1786596) q[3];
sx q[3];
rz(-0.41662595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.999324) q[2];
sx q[2];
rz(-0.48196718) q[2];
sx q[2];
rz(-1.8959321) q[2];
rz(-1.2549531) q[3];
sx q[3];
rz(-0.84147036) q[3];
sx q[3];
rz(-0.83827034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6923043) q[0];
sx q[0];
rz(-3.1384387) q[0];
sx q[0];
rz(-2.4601049) q[0];
rz(-0.20755126) q[1];
sx q[1];
rz(-0.46336585) q[1];
sx q[1];
rz(2.025827) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0633495) q[0];
sx q[0];
rz(-0.82394281) q[0];
sx q[0];
rz(-2.1372165) q[0];
x q[1];
rz(0.43004604) q[2];
sx q[2];
rz(-2.1950245) q[2];
sx q[2];
rz(-2.4732694) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7129732) q[1];
sx q[1];
rz(-1.5884807) q[1];
sx q[1];
rz(-0.58847217) q[1];
rz(0.66949087) q[3];
sx q[3];
rz(-1.4298425) q[3];
sx q[3];
rz(1.475875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9289124) q[2];
sx q[2];
rz(-1.2968411) q[2];
sx q[2];
rz(-2.7872655) q[2];
rz(-2.8220693) q[3];
sx q[3];
rz(-1.1525681) q[3];
sx q[3];
rz(-0.37187809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1324683) q[0];
sx q[0];
rz(-2.9102944) q[0];
sx q[0];
rz(-0.67434597) q[0];
rz(-1.1122423) q[1];
sx q[1];
rz(-0.66450417) q[1];
sx q[1];
rz(-2.5792714) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4599265) q[0];
sx q[0];
rz(-0.66837464) q[0];
sx q[0];
rz(-1.5947123) q[0];
x q[1];
rz(1.4698896) q[2];
sx q[2];
rz(-1.9011874) q[2];
sx q[2];
rz(0.23320564) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6866236) q[1];
sx q[1];
rz(-1.7640055) q[1];
sx q[1];
rz(-1.7186233) q[1];
x q[2];
rz(-0.9741707) q[3];
sx q[3];
rz(-1.4985634) q[3];
sx q[3];
rz(2.7330287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.101863) q[2];
sx q[2];
rz(-0.9938643) q[2];
sx q[2];
rz(2.8015461) q[2];
rz(2.9240821) q[3];
sx q[3];
rz(-0.9483996) q[3];
sx q[3];
rz(0.31869179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6927004) q[0];
sx q[0];
rz(-0.046148766) q[0];
sx q[0];
rz(-0.39644077) q[0];
rz(-0.13892826) q[1];
sx q[1];
rz(-0.46008343) q[1];
sx q[1];
rz(1.6202392) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0884468) q[0];
sx q[0];
rz(-1.3741115) q[0];
sx q[0];
rz(-1.501207) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2659726) q[2];
sx q[2];
rz(-2.5555829) q[2];
sx q[2];
rz(-0.74967521) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0115067) q[1];
sx q[1];
rz(-2.3475921) q[1];
sx q[1];
rz(-1.6542875) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7117386) q[3];
sx q[3];
rz(-0.68416506) q[3];
sx q[3];
rz(-1.14389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5788995) q[2];
sx q[2];
rz(-1.0572628) q[2];
sx q[2];
rz(0.29433027) q[2];
rz(2.0108022) q[3];
sx q[3];
rz(-1.7629938) q[3];
sx q[3];
rz(2.14595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49333736) q[0];
sx q[0];
rz(-0.91518891) q[0];
sx q[0];
rz(-2.7822568) q[0];
rz(0.94611478) q[1];
sx q[1];
rz(-0.39603907) q[1];
sx q[1];
rz(2.8709581) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5198869) q[0];
sx q[0];
rz(-1.5396376) q[0];
sx q[0];
rz(3.1057538) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6793849) q[2];
sx q[2];
rz(-0.86420977) q[2];
sx q[2];
rz(-0.30009899) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.7120413) q[1];
sx q[1];
rz(-1.6152641) q[1];
sx q[1];
rz(-0.045750381) q[1];
rz(-pi) q[2];
rz(-2.1569096) q[3];
sx q[3];
rz(-0.95041785) q[3];
sx q[3];
rz(-1.2445205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3140807) q[2];
sx q[2];
rz(-2.9120047) q[2];
sx q[2];
rz(2.6861526) q[2];
rz(-2.3296302) q[3];
sx q[3];
rz(-0.88192314) q[3];
sx q[3];
rz(2.5922095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6311326) q[0];
sx q[0];
rz(-1.5083418) q[0];
sx q[0];
rz(0.73927885) q[0];
rz(2.9108858) q[1];
sx q[1];
rz(-2.6735327) q[1];
sx q[1];
rz(-2.646692) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56787581) q[0];
sx q[0];
rz(-2.0567237) q[0];
sx q[0];
rz(-3.0987415) q[0];
x q[1];
rz(-2.8137389) q[2];
sx q[2];
rz(-1.577539) q[2];
sx q[2];
rz(-1.1525796) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4329263) q[1];
sx q[1];
rz(-0.64879829) q[1];
sx q[1];
rz(1.385958) q[1];
rz(-pi) q[2];
rz(1.4570191) q[3];
sx q[3];
rz(-2.4416231) q[3];
sx q[3];
rz(1.9663119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0080002) q[2];
sx q[2];
rz(-1.091489) q[2];
sx q[2];
rz(2.8137394) q[2];
rz(-2.949529) q[3];
sx q[3];
rz(-2.8938507) q[3];
sx q[3];
rz(2.1081934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0192169) q[0];
sx q[0];
rz(-1.6288971) q[0];
sx q[0];
rz(-1.6678641) q[0];
rz(-2.3090251) q[1];
sx q[1];
rz(-1.7201798) q[1];
sx q[1];
rz(-1.7725772) q[1];
rz(1.9033296) q[2];
sx q[2];
rz(-2.5730972) q[2];
sx q[2];
rz(2.2039883) q[2];
rz(0.26631793) q[3];
sx q[3];
rz(-1.8177633) q[3];
sx q[3];
rz(0.37533356) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
