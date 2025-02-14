OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.2849046) q[0];
sx q[0];
rz(-1.8104799) q[0];
sx q[0];
rz(-0.80944219) q[0];
rz(-2.5419905) q[1];
sx q[1];
rz(-1.124958) q[1];
sx q[1];
rz(0.52176276) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7441352) q[0];
sx q[0];
rz(-1.7108727) q[0];
sx q[0];
rz(3.1404353) q[0];
rz(-pi) q[1];
x q[1];
rz(0.56057616) q[2];
sx q[2];
rz(-1.3180705) q[2];
sx q[2];
rz(-1.8552903) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.062816) q[1];
sx q[1];
rz(-2.605781) q[1];
sx q[1];
rz(1.0173372) q[1];
rz(0.14484804) q[3];
sx q[3];
rz(-2.600935) q[3];
sx q[3];
rz(0.40504211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1521505) q[2];
sx q[2];
rz(-2.4336954) q[2];
sx q[2];
rz(-1.36261) q[2];
rz(-1.4710434) q[3];
sx q[3];
rz(-1.5444642) q[3];
sx q[3];
rz(2.7267) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7489557) q[0];
sx q[0];
rz(-1.7935268) q[0];
sx q[0];
rz(-2.0718527) q[0];
rz(2.3906129) q[1];
sx q[1];
rz(-2.2396478) q[1];
sx q[1];
rz(2.353207) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25229821) q[0];
sx q[0];
rz(-1.7579334) q[0];
sx q[0];
rz(2.5146913) q[0];
rz(-pi) q[1];
rz(2.4761202) q[2];
sx q[2];
rz(-1.62684) q[2];
sx q[2];
rz(-1.9201345) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.8395811) q[1];
sx q[1];
rz(-2.5835648) q[1];
sx q[1];
rz(1.1990947) q[1];
rz(-pi) q[2];
rz(-0.69425868) q[3];
sx q[3];
rz(-2.5115847) q[3];
sx q[3];
rz(-0.97351551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.02701935) q[2];
sx q[2];
rz(-2.7832649) q[2];
sx q[2];
rz(-0.9066073) q[2];
rz(-1.00057) q[3];
sx q[3];
rz(-1.118719) q[3];
sx q[3];
rz(-2.6226543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24666102) q[0];
sx q[0];
rz(-2.7592359) q[0];
sx q[0];
rz(1.3038127) q[0];
rz(-1.1600102) q[1];
sx q[1];
rz(-1.8142533) q[1];
sx q[1];
rz(1.7283641) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6400063) q[0];
sx q[0];
rz(-1.7118071) q[0];
sx q[0];
rz(1.4835351) q[0];
rz(-pi) q[1];
rz(-2.5557466) q[2];
sx q[2];
rz(-2.4163462) q[2];
sx q[2];
rz(-2.9278529) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1294413) q[1];
sx q[1];
rz(-1.1124764) q[1];
sx q[1];
rz(-1.767551) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.71536) q[3];
sx q[3];
rz(-0.72025296) q[3];
sx q[3];
rz(1.6584058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.383519) q[2];
sx q[2];
rz(-1.588409) q[2];
sx q[2];
rz(-0.12913945) q[2];
rz(-0.50390759) q[3];
sx q[3];
rz(-1.7241071) q[3];
sx q[3];
rz(2.5655139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1103766) q[0];
sx q[0];
rz(-1.2126558) q[0];
sx q[0];
rz(0.84621286) q[0];
rz(2.5695678) q[1];
sx q[1];
rz(-1.6366448) q[1];
sx q[1];
rz(1.2342854) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9548428) q[0];
sx q[0];
rz(-0.34663793) q[0];
sx q[0];
rz(-0.28916547) q[0];
rz(-1.1115121) q[2];
sx q[2];
rz(-1.827335) q[2];
sx q[2];
rz(0.81748325) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7031786) q[1];
sx q[1];
rz(-1.4902189) q[1];
sx q[1];
rz(0.2097278) q[1];
rz(-pi) q[2];
rz(-0.79674308) q[3];
sx q[3];
rz(-1.9588184) q[3];
sx q[3];
rz(1.4772367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7531551) q[2];
sx q[2];
rz(-2.2074102) q[2];
sx q[2];
rz(1.4446806) q[2];
rz(-1.2706903) q[3];
sx q[3];
rz(-1.5710187) q[3];
sx q[3];
rz(1.9522033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5020849) q[0];
sx q[0];
rz(-1.4690228) q[0];
sx q[0];
rz(2.0462659) q[0];
rz(-0.0630088) q[1];
sx q[1];
rz(-1.4728225) q[1];
sx q[1];
rz(0.94620401) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5216833) q[0];
sx q[0];
rz(-2.5961726) q[0];
sx q[0];
rz(1.5167461) q[0];
x q[1];
rz(-2.6465552) q[2];
sx q[2];
rz(-2.4433141) q[2];
sx q[2];
rz(-2.3586522) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7530147) q[1];
sx q[1];
rz(-1.565897) q[1];
sx q[1];
rz(-0.51112643) q[1];
rz(1.9619476) q[3];
sx q[3];
rz(-1.7877276) q[3];
sx q[3];
rz(0.85238487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5064064) q[2];
sx q[2];
rz(-2.0621641) q[2];
sx q[2];
rz(-1.5080473) q[2];
rz(-2.9888198) q[3];
sx q[3];
rz(-1.907405) q[3];
sx q[3];
rz(-2.2085786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2372811) q[0];
sx q[0];
rz(-2.5560684) q[0];
sx q[0];
rz(3.0329419) q[0];
rz(-0.12734224) q[1];
sx q[1];
rz(-1.2510977) q[1];
sx q[1];
rz(2.2551575) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21289794) q[0];
sx q[0];
rz(-1.587364) q[0];
sx q[0];
rz(2.6221656) q[0];
x q[1];
rz(1.3662874) q[2];
sx q[2];
rz(-1.5123774) q[2];
sx q[2];
rz(1.2893334) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.163609) q[1];
sx q[1];
rz(-0.22630461) q[1];
sx q[1];
rz(-1.8288906) q[1];
rz(-pi) q[2];
rz(-1.4599019) q[3];
sx q[3];
rz(-1.4023047) q[3];
sx q[3];
rz(0.061010188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2966557) q[2];
sx q[2];
rz(-2.3462494) q[2];
sx q[2];
rz(-2.8080158) q[2];
rz(-0.9345471) q[3];
sx q[3];
rz(-2.3881674) q[3];
sx q[3];
rz(-0.88753382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4535256) q[0];
sx q[0];
rz(-1.4997361) q[0];
sx q[0];
rz(2.8330579) q[0];
rz(0.60641369) q[1];
sx q[1];
rz(-2.2589222) q[1];
sx q[1];
rz(-0.44953129) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.318293) q[0];
sx q[0];
rz(-1.6345981) q[0];
sx q[0];
rz(-2.6266111) q[0];
rz(2.6753898) q[2];
sx q[2];
rz(-2.297819) q[2];
sx q[2];
rz(-2.1920993) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3492185) q[1];
sx q[1];
rz(-1.6728396) q[1];
sx q[1];
rz(-2.1898063) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5964162) q[3];
sx q[3];
rz(-2.1948493) q[3];
sx q[3];
rz(1.6069309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.061607925) q[2];
sx q[2];
rz(-1.4915165) q[2];
sx q[2];
rz(1.9645346) q[2];
rz(2.4401149) q[3];
sx q[3];
rz(-0.25291118) q[3];
sx q[3];
rz(-2.285932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45706448) q[0];
sx q[0];
rz(-0.51790154) q[0];
sx q[0];
rz(1.49217) q[0];
rz(1.0555142) q[1];
sx q[1];
rz(-1.5558473) q[1];
sx q[1];
rz(-2.7141056) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1221177) q[0];
sx q[0];
rz(-2.1699598) q[0];
sx q[0];
rz(-2.504306) q[0];
x q[1];
rz(-2.7107377) q[2];
sx q[2];
rz(-1.4240032) q[2];
sx q[2];
rz(2.1374747) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.52054469) q[1];
sx q[1];
rz(-1.2579903) q[1];
sx q[1];
rz(-1.5942082) q[1];
rz(-2.8264846) q[3];
sx q[3];
rz(-1.3118532) q[3];
sx q[3];
rz(2.0729947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.18275729) q[2];
sx q[2];
rz(-1.3311102) q[2];
sx q[2];
rz(-1.2048362) q[2];
rz(3.0086573) q[3];
sx q[3];
rz(-0.098630579) q[3];
sx q[3];
rz(-1.0814063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8590915) q[0];
sx q[0];
rz(-1.4097255) q[0];
sx q[0];
rz(-0.45243636) q[0];
rz(0.85894194) q[1];
sx q[1];
rz(-0.57933885) q[1];
sx q[1];
rz(1.6890242) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6159264) q[0];
sx q[0];
rz(-0.99352628) q[0];
sx q[0];
rz(2.1494715) q[0];
rz(2.0797009) q[2];
sx q[2];
rz(-2.5246722) q[2];
sx q[2];
rz(-1.2764608) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7048129) q[1];
sx q[1];
rz(-1.5734993) q[1];
sx q[1];
rz(2.0767434) q[1];
x q[2];
rz(-2.3708569) q[3];
sx q[3];
rz(-0.49388921) q[3];
sx q[3];
rz(2.1300742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9026044) q[2];
sx q[2];
rz(-1.8391823) q[2];
sx q[2];
rz(0.39546173) q[2];
rz(-1.0836733) q[3];
sx q[3];
rz(-0.62267059) q[3];
sx q[3];
rz(-1.4130939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85635066) q[0];
sx q[0];
rz(-0.1150035) q[0];
sx q[0];
rz(1.2913936) q[0];
rz(-0.77766386) q[1];
sx q[1];
rz(-1.5592557) q[1];
sx q[1];
rz(1.8596328) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1384821) q[0];
sx q[0];
rz(-1.3457451) q[0];
sx q[0];
rz(-2.9182994) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6474495) q[2];
sx q[2];
rz(-1.3927336) q[2];
sx q[2];
rz(-0.21221976) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.79445542) q[1];
sx q[1];
rz(-1.0442808) q[1];
sx q[1];
rz(1.8361676) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.090023) q[3];
sx q[3];
rz(-1.6656205) q[3];
sx q[3];
rz(-2.7213396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.164244) q[2];
sx q[2];
rz(-2.7898596) q[2];
sx q[2];
rz(0.39883167) q[2];
rz(-0.93419689) q[3];
sx q[3];
rz(-0.7312921) q[3];
sx q[3];
rz(-3.0493693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6964523) q[0];
sx q[0];
rz(-2.2284989) q[0];
sx q[0];
rz(-3.1365119) q[0];
rz(-2.483881) q[1];
sx q[1];
rz(-1.4124159) q[1];
sx q[1];
rz(1.1884069) q[1];
rz(-2.3376113) q[2];
sx q[2];
rz(-1.8009427) q[2];
sx q[2];
rz(1.367205) q[2];
rz(-3.0791238) q[3];
sx q[3];
rz(-2.1556839) q[3];
sx q[3];
rz(1.5869303) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
