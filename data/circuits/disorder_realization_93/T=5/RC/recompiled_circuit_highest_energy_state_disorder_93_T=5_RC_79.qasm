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
rz(-0.45577249) q[0];
sx q[0];
rz(-1.1205641) q[0];
sx q[0];
rz(1.4135345) q[0];
rz(-2.7780374) q[1];
sx q[1];
rz(-1.8283365) q[1];
sx q[1];
rz(0.29120905) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4158299) q[0];
sx q[0];
rz(-0.32074577) q[0];
sx q[0];
rz(1.5070391) q[0];
x q[1];
rz(1.9102664) q[2];
sx q[2];
rz(-1.6874773) q[2];
sx q[2];
rz(1.9872023) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0800228) q[1];
sx q[1];
rz(-2.1407619) q[1];
sx q[1];
rz(1.7504196) q[1];
rz(-pi) q[2];
rz(-2.324027) q[3];
sx q[3];
rz(-2.052784) q[3];
sx q[3];
rz(0.40460247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.0061965813) q[2];
sx q[2];
rz(-1.5634894) q[2];
sx q[2];
rz(0.023999365) q[2];
rz(-0.35605797) q[3];
sx q[3];
rz(-2.4672716) q[3];
sx q[3];
rz(-3.1168028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.085670797) q[0];
sx q[0];
rz(-0.69127685) q[0];
sx q[0];
rz(-0.63496494) q[0];
rz(2.3786646) q[1];
sx q[1];
rz(-1.678391) q[1];
sx q[1];
rz(-0.91711226) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2208495) q[0];
sx q[0];
rz(-1.5524327) q[0];
sx q[0];
rz(1.5361274) q[0];
rz(-pi) q[1];
rz(2.2516052) q[2];
sx q[2];
rz(-1.1900848) q[2];
sx q[2];
rz(-0.057304545) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.32815427) q[1];
sx q[1];
rz(-2.0106689) q[1];
sx q[1];
rz(-0.085170345) q[1];
rz(1.8648326) q[3];
sx q[3];
rz(-2.5007479) q[3];
sx q[3];
rz(-1.8185735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.94747535) q[2];
sx q[2];
rz(-2.6185991) q[2];
sx q[2];
rz(-0.73327649) q[2];
rz(2.0745847) q[3];
sx q[3];
rz(-1.4209483) q[3];
sx q[3];
rz(2.9616621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5697524) q[0];
sx q[0];
rz(-1.952992) q[0];
sx q[0];
rz(-0.52823129) q[0];
rz(-0.86604467) q[1];
sx q[1];
rz(-1.3287611) q[1];
sx q[1];
rz(-2.6285062) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.00947) q[0];
sx q[0];
rz(-2.7010601) q[0];
sx q[0];
rz(1.4709573) q[0];
rz(-pi) q[1];
rz(0.31103525) q[2];
sx q[2];
rz(-2.1967109) q[2];
sx q[2];
rz(-3.1300822) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5410903) q[1];
sx q[1];
rz(-1.9185431) q[1];
sx q[1];
rz(3.0060124) q[1];
x q[2];
rz(-1.5422899) q[3];
sx q[3];
rz(-1.0712726) q[3];
sx q[3];
rz(1.0699501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.036006007) q[2];
sx q[2];
rz(-1.0893818) q[2];
sx q[2];
rz(-0.038724381) q[2];
rz(0.45768467) q[3];
sx q[3];
rz(-0.80515146) q[3];
sx q[3];
rz(1.3409486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7078581) q[0];
sx q[0];
rz(-0.41446328) q[0];
sx q[0];
rz(-1.1210972) q[0];
rz(-0.7577678) q[1];
sx q[1];
rz(-1.6242177) q[1];
sx q[1];
rz(0.68414348) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0476889) q[0];
sx q[0];
rz(-0.8465811) q[0];
sx q[0];
rz(-2.8709477) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5779547) q[2];
sx q[2];
rz(-1.5809007) q[2];
sx q[2];
rz(-1.9310538) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.18312081) q[1];
sx q[1];
rz(-2.712912) q[1];
sx q[1];
rz(2.3752992) q[1];
x q[2];
rz(0.97400333) q[3];
sx q[3];
rz(-2.0066378) q[3];
sx q[3];
rz(1.7929994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.4844369) q[2];
sx q[2];
rz(-1.5269205) q[2];
sx q[2];
rz(0.26051513) q[2];
rz(-0.83812964) q[3];
sx q[3];
rz(-1.9939491) q[3];
sx q[3];
rz(-3.1269791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.218006) q[0];
sx q[0];
rz(-1.3777233) q[0];
sx q[0];
rz(0.51061428) q[0];
rz(-1.249373) q[1];
sx q[1];
rz(-1.9793972) q[1];
sx q[1];
rz(-1.4543264) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.175486) q[0];
sx q[0];
rz(-1.5990055) q[0];
sx q[0];
rz(2.9278446) q[0];
rz(-1.6188443) q[2];
sx q[2];
rz(-1.4312051) q[2];
sx q[2];
rz(-2.0968116) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0328503) q[1];
sx q[1];
rz(-1.3128442) q[1];
sx q[1];
rz(-2.5775419) q[1];
rz(-pi) q[2];
x q[2];
rz(0.84564836) q[3];
sx q[3];
rz(-2.925784) q[3];
sx q[3];
rz(-0.25419054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5877567) q[2];
sx q[2];
rz(-2.7601056) q[2];
sx q[2];
rz(-0.39056632) q[2];
rz(0.25911123) q[3];
sx q[3];
rz(-1.6389537) q[3];
sx q[3];
rz(0.34671569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7873586) q[0];
sx q[0];
rz(-2.7095095) q[0];
sx q[0];
rz(-2.4969192) q[0];
rz(1.3020172) q[1];
sx q[1];
rz(-1.6141067) q[1];
sx q[1];
rz(2.7929746) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88953274) q[0];
sx q[0];
rz(-2.4675998) q[0];
sx q[0];
rz(-0.78788449) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.63747823) q[2];
sx q[2];
rz(-2.5036252) q[2];
sx q[2];
rz(1.5359985) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.4069566) q[1];
sx q[1];
rz(-2.6183081) q[1];
sx q[1];
rz(-0.95013817) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.11696924) q[3];
sx q[3];
rz(-1.6050964) q[3];
sx q[3];
rz(3.0142865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0384486) q[2];
sx q[2];
rz(-0.95547533) q[2];
sx q[2];
rz(0.038334282) q[2];
rz(0.96758715) q[3];
sx q[3];
rz(-1.4097694) q[3];
sx q[3];
rz(0.35965317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86647812) q[0];
sx q[0];
rz(-2.622128) q[0];
sx q[0];
rz(-3.0679829) q[0];
rz(1.4069125) q[1];
sx q[1];
rz(-0.55714566) q[1];
sx q[1];
rz(-0.73307347) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.063696472) q[0];
sx q[0];
rz(-2.1025582) q[0];
sx q[0];
rz(-2.1666906) q[0];
rz(-pi) q[1];
rz(-2.2233811) q[2];
sx q[2];
rz(-2.2951153) q[2];
sx q[2];
rz(-1.0670964) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9736874) q[1];
sx q[1];
rz(-1.2406772) q[1];
sx q[1];
rz(-0.37574212) q[1];
x q[2];
rz(-0.81871512) q[3];
sx q[3];
rz(-0.33740852) q[3];
sx q[3];
rz(-0.81160566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1072032) q[2];
sx q[2];
rz(-1.4185536) q[2];
sx q[2];
rz(3.1305195) q[2];
rz(0.30558807) q[3];
sx q[3];
rz(-2.0162851) q[3];
sx q[3];
rz(-1.2529469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8082064) q[0];
sx q[0];
rz(-2.5337063) q[0];
sx q[0];
rz(-2.407684) q[0];
rz(0.58087307) q[1];
sx q[1];
rz(-1.0092694) q[1];
sx q[1];
rz(0.098800585) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0721314) q[0];
sx q[0];
rz(-2.0939576) q[0];
sx q[0];
rz(-2.4510372) q[0];
x q[1];
rz(0.71008981) q[2];
sx q[2];
rz(-0.81560336) q[2];
sx q[2];
rz(0.85271533) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.72291127) q[1];
sx q[1];
rz(-2.0127814) q[1];
sx q[1];
rz(1.8479455) q[1];
rz(-pi) q[2];
rz(-3.0116073) q[3];
sx q[3];
rz(-1.49125) q[3];
sx q[3];
rz(1.6787488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8431479) q[2];
sx q[2];
rz(-1.8567905) q[2];
sx q[2];
rz(-2.0419545) q[2];
rz(0.034865033) q[3];
sx q[3];
rz(-0.74368447) q[3];
sx q[3];
rz(2.5449424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.027503969) q[0];
sx q[0];
rz(-1.1570258) q[0];
sx q[0];
rz(1.1771359) q[0];
rz(-1.6674532) q[1];
sx q[1];
rz(-0.93282229) q[1];
sx q[1];
rz(1.4449545) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5367573) q[0];
sx q[0];
rz(-0.49237456) q[0];
sx q[0];
rz(2.5144666) q[0];
rz(-pi) q[1];
rz(-2.3722052) q[2];
sx q[2];
rz(-1.1379998) q[2];
sx q[2];
rz(2.3355049) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.18981537) q[1];
sx q[1];
rz(-0.51945247) q[1];
sx q[1];
rz(1.7096667) q[1];
rz(-0.27650303) q[3];
sx q[3];
rz(-2.2275337) q[3];
sx q[3];
rz(-1.2402676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.16704796) q[2];
sx q[2];
rz(-1.5844774) q[2];
sx q[2];
rz(-2.9696999) q[2];
rz(0.78957549) q[3];
sx q[3];
rz(-2.004345) q[3];
sx q[3];
rz(1.3353698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0563141) q[0];
sx q[0];
rz(-0.21509898) q[0];
sx q[0];
rz(-0.54157448) q[0];
rz(-1.1594634) q[1];
sx q[1];
rz(-1.3075202) q[1];
sx q[1];
rz(2.3084739) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7689) q[0];
sx q[0];
rz(-1.295255) q[0];
sx q[0];
rz(-2.5252456) q[0];
rz(-pi) q[1];
x q[1];
rz(0.53785777) q[2];
sx q[2];
rz(-0.61213697) q[2];
sx q[2];
rz(-2.9878841) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.45634746) q[1];
sx q[1];
rz(-0.73328555) q[1];
sx q[1];
rz(1.9128146) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5686137) q[3];
sx q[3];
rz(-0.72764054) q[3];
sx q[3];
rz(1.997987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.021576015) q[2];
sx q[2];
rz(-0.77777255) q[2];
sx q[2];
rz(2.1982101) q[2];
rz(-2.0639482) q[3];
sx q[3];
rz(-1.5543944) q[3];
sx q[3];
rz(-0.32850346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42644603) q[0];
sx q[0];
rz(-2.1052512) q[0];
sx q[0];
rz(-0.25682009) q[0];
rz(0.23964755) q[1];
sx q[1];
rz(-2.2365204) q[1];
sx q[1];
rz(2.0047275) q[1];
rz(-0.93076651) q[2];
sx q[2];
rz(-0.98251828) q[2];
sx q[2];
rz(0.3457014) q[2];
rz(-1.8244459) q[3];
sx q[3];
rz(-1.8719049) q[3];
sx q[3];
rz(-0.20897839) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
