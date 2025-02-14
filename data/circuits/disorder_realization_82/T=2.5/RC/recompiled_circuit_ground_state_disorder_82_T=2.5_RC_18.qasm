OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.3759484) q[0];
sx q[0];
rz(-2.7274237) q[0];
sx q[0];
rz(0.8015269) q[0];
rz(-0.0030567788) q[1];
sx q[1];
rz(-0.86441511) q[1];
sx q[1];
rz(-0.094223082) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4163001) q[0];
sx q[0];
rz(-2.2803118) q[0];
sx q[0];
rz(2.6363591) q[0];
x q[1];
rz(1.8434486) q[2];
sx q[2];
rz(-1.5412207) q[2];
sx q[2];
rz(-2.3675413) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7677418) q[1];
sx q[1];
rz(-0.5597502) q[1];
sx q[1];
rz(-2.1327726) q[1];
x q[2];
rz(2.9889936) q[3];
sx q[3];
rz(-1.5729701) q[3];
sx q[3];
rz(-1.847037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6039383) q[2];
sx q[2];
rz(-2.5014169) q[2];
sx q[2];
rz(1.5113277) q[2];
rz(2.9553735) q[3];
sx q[3];
rz(-0.43034601) q[3];
sx q[3];
rz(-2.490624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6179825) q[0];
sx q[0];
rz(-1.9600493) q[0];
sx q[0];
rz(-0.13667662) q[0];
rz(-0.094820529) q[1];
sx q[1];
rz(-0.46351981) q[1];
sx q[1];
rz(-0.76900855) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5154607) q[0];
sx q[0];
rz(-0.40053408) q[0];
sx q[0];
rz(-0.35450165) q[0];
x q[1];
rz(-0.42119512) q[2];
sx q[2];
rz(-2.4026818) q[2];
sx q[2];
rz(-0.74904672) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8254346) q[1];
sx q[1];
rz(-0.76091754) q[1];
sx q[1];
rz(-2.9050211) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7058825) q[3];
sx q[3];
rz(-1.6052433) q[3];
sx q[3];
rz(-2.5273539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.39917055) q[2];
sx q[2];
rz(-2.6593282) q[2];
sx q[2];
rz(-1.7102309) q[2];
rz(-2.6509905) q[3];
sx q[3];
rz(-0.49121818) q[3];
sx q[3];
rz(-0.77559364) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2630149) q[0];
sx q[0];
rz(-0.37913015) q[0];
sx q[0];
rz(-1.0947134) q[0];
rz(1.3021944) q[1];
sx q[1];
rz(-2.1523988) q[1];
sx q[1];
rz(-0.0040231752) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1510295) q[0];
sx q[0];
rz(-0.73712611) q[0];
sx q[0];
rz(-1.3746136) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.268774) q[2];
sx q[2];
rz(-2.1790163) q[2];
sx q[2];
rz(-1.5814511) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0931172) q[1];
sx q[1];
rz(-2.8826536) q[1];
sx q[1];
rz(-3.0228651) q[1];
rz(-pi) q[2];
rz(1.8173006) q[3];
sx q[3];
rz(-0.87216264) q[3];
sx q[3];
rz(2.7072631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5912938) q[2];
sx q[2];
rz(-2.55547) q[2];
sx q[2];
rz(0.44814056) q[2];
rz(2.6552933) q[3];
sx q[3];
rz(-0.86217642) q[3];
sx q[3];
rz(-2.015131) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9124741) q[0];
sx q[0];
rz(-1.4294701) q[0];
sx q[0];
rz(0.66977704) q[0];
rz(-1.9122596) q[1];
sx q[1];
rz(-2.6826617) q[1];
sx q[1];
rz(-2.5965447) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16825039) q[0];
sx q[0];
rz(-1.1118044) q[0];
sx q[0];
rz(0.58990546) q[0];
rz(2.3178905) q[2];
sx q[2];
rz(-2.3266226) q[2];
sx q[2];
rz(-0.65836775) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.74074827) q[1];
sx q[1];
rz(-1.5103589) q[1];
sx q[1];
rz(-1.7857741) q[1];
rz(-pi) q[2];
rz(2.6757984) q[3];
sx q[3];
rz(-0.85448962) q[3];
sx q[3];
rz(-0.6443302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.58686078) q[2];
sx q[2];
rz(-1.6394337) q[2];
sx q[2];
rz(0.16793212) q[2];
rz(-1.2083017) q[3];
sx q[3];
rz(-0.30253634) q[3];
sx q[3];
rz(-1.6944073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.55030441) q[0];
sx q[0];
rz(-1.8809603) q[0];
sx q[0];
rz(0.0038797832) q[0];
rz(-0.30715352) q[1];
sx q[1];
rz(-0.75633621) q[1];
sx q[1];
rz(-2.602813) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8910485) q[0];
sx q[0];
rz(-1.4293423) q[0];
sx q[0];
rz(-0.23141872) q[0];
rz(-1.7096552) q[2];
sx q[2];
rz(-1.1449688) q[2];
sx q[2];
rz(-2.6151997) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5153577) q[1];
sx q[1];
rz(-1.9222676) q[1];
sx q[1];
rz(2.5280747) q[1];
rz(-pi) q[2];
rz(-1.1191551) q[3];
sx q[3];
rz(-0.8913826) q[3];
sx q[3];
rz(1.7070626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6946081) q[2];
sx q[2];
rz(-1.0993404) q[2];
sx q[2];
rz(-1.2896607) q[2];
rz(-1.5223711) q[3];
sx q[3];
rz(-0.74519849) q[3];
sx q[3];
rz(1.4815909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
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
rz(2.8412142) q[0];
sx q[0];
rz(-0.9391681) q[0];
sx q[0];
rz(-3.0389431) q[0];
rz(2.4094021) q[1];
sx q[1];
rz(-2.3468572) q[1];
sx q[1];
rz(0.57714677) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10023663) q[0];
sx q[0];
rz(-1.5994342) q[0];
sx q[0];
rz(-3.1379382) q[0];
rz(-pi) q[1];
x q[1];
rz(0.58697613) q[2];
sx q[2];
rz(-1.7086626) q[2];
sx q[2];
rz(-1.1816813) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0967193) q[1];
sx q[1];
rz(-1.4920939) q[1];
sx q[1];
rz(-3.031879) q[1];
x q[2];
rz(-1.243337) q[3];
sx q[3];
rz(-2.2994882) q[3];
sx q[3];
rz(3.1407243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.010043667) q[2];
sx q[2];
rz(-2.4250344) q[2];
sx q[2];
rz(-1.5232167) q[2];
rz(0.067954436) q[3];
sx q[3];
rz(-2.8165635) q[3];
sx q[3];
rz(0.84097356) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1286569) q[0];
sx q[0];
rz(-1.034863) q[0];
sx q[0];
rz(-2.3701684) q[0];
rz(2.6029288) q[1];
sx q[1];
rz(-1.3464876) q[1];
sx q[1];
rz(2.0753986) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5911884) q[0];
sx q[0];
rz(-0.26445358) q[0];
sx q[0];
rz(1.868684) q[0];
rz(-pi) q[1];
rz(0.23287878) q[2];
sx q[2];
rz(-2.072273) q[2];
sx q[2];
rz(1.6362783) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.33246189) q[1];
sx q[1];
rz(-1.5839424) q[1];
sx q[1];
rz(-0.21717182) q[1];
x q[2];
rz(-2.8071219) q[3];
sx q[3];
rz(-1.3190509) q[3];
sx q[3];
rz(-1.2125963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6923339) q[2];
sx q[2];
rz(-2.3350495) q[2];
sx q[2];
rz(2.6991357) q[2];
rz(0.19065204) q[3];
sx q[3];
rz(-0.72398829) q[3];
sx q[3];
rz(-0.75907069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2451179) q[0];
sx q[0];
rz(-2.9560095) q[0];
sx q[0];
rz(1.3099439) q[0];
rz(-0.39778057) q[1];
sx q[1];
rz(-2.3186627) q[1];
sx q[1];
rz(1.3936874) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8373972) q[0];
sx q[0];
rz(-2.1050354) q[0];
sx q[0];
rz(-0.5128575) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9056121) q[2];
sx q[2];
rz(-1.8644635) q[2];
sx q[2];
rz(-2.8666509) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9164394) q[1];
sx q[1];
rz(-1.272456) q[1];
sx q[1];
rz(-0.48709309) q[1];
x q[2];
rz(1.3721136) q[3];
sx q[3];
rz(-1.4666838) q[3];
sx q[3];
rz(-2.7840691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1132249) q[2];
sx q[2];
rz(-0.71030474) q[2];
sx q[2];
rz(-0.073471546) q[2];
rz(-2.4469589) q[3];
sx q[3];
rz(-1.2725384) q[3];
sx q[3];
rz(-1.0276444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6756814) q[0];
sx q[0];
rz(-0.2008734) q[0];
sx q[0];
rz(-0.95941108) q[0];
rz(2.361946) q[1];
sx q[1];
rz(-2.1387073) q[1];
sx q[1];
rz(-0.27228212) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.773943) q[0];
sx q[0];
rz(-3.127357) q[0];
sx q[0];
rz(1.9444501) q[0];
rz(-0.30961664) q[2];
sx q[2];
rz(-2.2520492) q[2];
sx q[2];
rz(0.75971425) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.12676621) q[1];
sx q[1];
rz(-2.0527168) q[1];
sx q[1];
rz(2.6792106) q[1];
rz(0.7471146) q[3];
sx q[3];
rz(-2.0339587) q[3];
sx q[3];
rz(-0.85827352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4043364) q[2];
sx q[2];
rz(-1.5055089) q[2];
sx q[2];
rz(-2.9340202) q[2];
rz(-0.10457822) q[3];
sx q[3];
rz(-0.50555491) q[3];
sx q[3];
rz(-1.526621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.680147) q[0];
sx q[0];
rz(-1.9141645) q[0];
sx q[0];
rz(-2.0848059) q[0];
rz(-0.54569221) q[1];
sx q[1];
rz(-0.97815198) q[1];
sx q[1];
rz(-0.32870764) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.064081505) q[0];
sx q[0];
rz(-1.2405378) q[0];
sx q[0];
rz(-2.0547377) q[0];
x q[1];
rz(1.516009) q[2];
sx q[2];
rz(-1.4984926) q[2];
sx q[2];
rz(-2.876578) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.88966072) q[1];
sx q[1];
rz(-1.228782) q[1];
sx q[1];
rz(0.017707822) q[1];
x q[2];
rz(-2.8601951) q[3];
sx q[3];
rz(-2.3499422) q[3];
sx q[3];
rz(0.36330128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1591961) q[2];
sx q[2];
rz(-0.13267645) q[2];
sx q[2];
rz(1.1245493) q[2];
rz(-3.0647965) q[3];
sx q[3];
rz(-0.43282893) q[3];
sx q[3];
rz(-2.2176149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9818253) q[0];
sx q[0];
rz(-1.8851017) q[0];
sx q[0];
rz(-1.5633352) q[0];
rz(2.3253597) q[1];
sx q[1];
rz(-1.7385794) q[1];
sx q[1];
rz(2.0508918) q[1];
rz(2.0682206) q[2];
sx q[2];
rz(-0.88928992) q[2];
sx q[2];
rz(2.5273565) q[2];
rz(-3.1327457) q[3];
sx q[3];
rz(-0.4053517) q[3];
sx q[3];
rz(2.3978618) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
