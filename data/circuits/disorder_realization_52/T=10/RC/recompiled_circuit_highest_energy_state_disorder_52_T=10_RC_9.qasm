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
rz(0.723847) q[0];
sx q[0];
rz(-1.2011733) q[0];
sx q[0];
rz(0.49402753) q[0];
rz(1.5552893) q[1];
sx q[1];
rz(2.2145693) q[1];
sx q[1];
rz(8.5742843) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14334003) q[0];
sx q[0];
rz(-1.7368642) q[0];
sx q[0];
rz(-2.2448426) q[0];
rz(-pi) q[1];
rz(0.65324983) q[2];
sx q[2];
rz(-1.946665) q[2];
sx q[2];
rz(-2.1819161) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.96877484) q[1];
sx q[1];
rz(-2.7451395) q[1];
sx q[1];
rz(0.60703599) q[1];
rz(-0.1876202) q[3];
sx q[3];
rz(-2.7696262) q[3];
sx q[3];
rz(0.99727977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4832619) q[2];
sx q[2];
rz(-2.5265145) q[2];
sx q[2];
rz(-2.1966546) q[2];
rz(3.1350709) q[3];
sx q[3];
rz(-0.76144731) q[3];
sx q[3];
rz(0.25750461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1061123) q[0];
sx q[0];
rz(-2.1529614) q[0];
sx q[0];
rz(2.535787) q[0];
rz(0.93217355) q[1];
sx q[1];
rz(-1.4235539) q[1];
sx q[1];
rz(-0.1056284) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86565844) q[0];
sx q[0];
rz(-1.4167042) q[0];
sx q[0];
rz(0.6886512) q[0];
x q[1];
rz(-1.0829686) q[2];
sx q[2];
rz(-1.8886856) q[2];
sx q[2];
rz(-2.6018104) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9739237) q[1];
sx q[1];
rz(-2.5906256) q[1];
sx q[1];
rz(2.5552804) q[1];
rz(-2.3096931) q[3];
sx q[3];
rz(-0.53153342) q[3];
sx q[3];
rz(-0.040415045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2848844) q[2];
sx q[2];
rz(-1.363089) q[2];
sx q[2];
rz(-1.523783) q[2];
rz(0.81426042) q[3];
sx q[3];
rz(-1.7819449) q[3];
sx q[3];
rz(2.2566569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.043269) q[0];
sx q[0];
rz(-2.3461778) q[0];
sx q[0];
rz(-1.5765618) q[0];
rz(-2.1463429) q[1];
sx q[1];
rz(-0.97882706) q[1];
sx q[1];
rz(-0.78688041) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9691447) q[0];
sx q[0];
rz(-2.8993164) q[0];
sx q[0];
rz(-2.1680225) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6549395) q[2];
sx q[2];
rz(-1.9168789) q[2];
sx q[2];
rz(3.0402407) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0643493) q[1];
sx q[1];
rz(-1.2573544) q[1];
sx q[1];
rz(-2.7948492) q[1];
rz(-0.059642369) q[3];
sx q[3];
rz(-1.630405) q[3];
sx q[3];
rz(0.27935057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3028822) q[2];
sx q[2];
rz(-1.6532712) q[2];
sx q[2];
rz(-2.5035456) q[2];
rz(0.4979411) q[3];
sx q[3];
rz(-2.0697856) q[3];
sx q[3];
rz(-2.0518484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-1.2447253) q[0];
sx q[0];
rz(-1.7843972) q[0];
sx q[0];
rz(-3.0778399) q[0];
rz(1.6925192) q[1];
sx q[1];
rz(-1.3056825) q[1];
sx q[1];
rz(1.8316899) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5777912) q[0];
sx q[0];
rz(-3.0534857) q[0];
sx q[0];
rz(0.43803517) q[0];
rz(-pi) q[1];
rz(-0.74964995) q[2];
sx q[2];
rz(-0.26826619) q[2];
sx q[2];
rz(-1.1531354) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.64778544) q[1];
sx q[1];
rz(-1.6946548) q[1];
sx q[1];
rz(0.74614831) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2579262) q[3];
sx q[3];
rz(-1.851463) q[3];
sx q[3];
rz(1.3786045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0706851) q[2];
sx q[2];
rz(-2.0347774) q[2];
sx q[2];
rz(0.10565383) q[2];
rz(1.1834772) q[3];
sx q[3];
rz(-2.2453997) q[3];
sx q[3];
rz(0.67160523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.79528177) q[0];
sx q[0];
rz(-2.2310937) q[0];
sx q[0];
rz(-1.7972535) q[0];
rz(0.67289871) q[1];
sx q[1];
rz(-1.3105323) q[1];
sx q[1];
rz(3.1281298) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5199842) q[0];
sx q[0];
rz(-0.88112393) q[0];
sx q[0];
rz(-3.1075928) q[0];
rz(-pi) q[1];
rz(0.014217397) q[2];
sx q[2];
rz(-1.9849615) q[2];
sx q[2];
rz(2.3384475) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2695737) q[1];
sx q[1];
rz(-2.7260352) q[1];
sx q[1];
rz(3.1088105) q[1];
rz(-pi) q[2];
rz(1.1290523) q[3];
sx q[3];
rz(-0.75296445) q[3];
sx q[3];
rz(0.66679614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.54231918) q[2];
sx q[2];
rz(-2.4794674) q[2];
sx q[2];
rz(1.2467965) q[2];
rz(-3.0689734) q[3];
sx q[3];
rz(-0.19425546) q[3];
sx q[3];
rz(0.3092002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4319864) q[0];
sx q[0];
rz(-2.015634) q[0];
sx q[0];
rz(-0.11269888) q[0];
rz(-2.9365149) q[1];
sx q[1];
rz(-1.3321184) q[1];
sx q[1];
rz(0.90788666) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4604707) q[0];
sx q[0];
rz(-2.2325071) q[0];
sx q[0];
rz(-0.76028334) q[0];
x q[1];
rz(2.5285883) q[2];
sx q[2];
rz(-1.2132436) q[2];
sx q[2];
rz(1.845971) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.0081387719) q[1];
sx q[1];
rz(-2.5027103) q[1];
sx q[1];
rz(-0.15365803) q[1];
rz(-1.3566293) q[3];
sx q[3];
rz(-0.91717824) q[3];
sx q[3];
rz(-0.71398338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9997361) q[2];
sx q[2];
rz(-0.33679589) q[2];
sx q[2];
rz(0.067961819) q[2];
rz(1.147602) q[3];
sx q[3];
rz(-1.2679029) q[3];
sx q[3];
rz(-1.8806774) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1934018) q[0];
sx q[0];
rz(-1.8302487) q[0];
sx q[0];
rz(-0.46352682) q[0];
rz(-0.14687982) q[1];
sx q[1];
rz(-1.905922) q[1];
sx q[1];
rz(2.1489977) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8604763) q[0];
sx q[0];
rz(-2.2669811) q[0];
sx q[0];
rz(1.244522) q[0];
rz(1.8568749) q[2];
sx q[2];
rz(-1.6353288) q[2];
sx q[2];
rz(-0.71071029) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.9001683) q[1];
sx q[1];
rz(-1.1196152) q[1];
sx q[1];
rz(-2.3996949) q[1];
rz(-pi) q[2];
x q[2];
rz(0.49673523) q[3];
sx q[3];
rz(-2.3443446) q[3];
sx q[3];
rz(2.0172142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2831882) q[2];
sx q[2];
rz(-1.7121199) q[2];
sx q[2];
rz(2.6962213) q[2];
rz(-1.8177659) q[3];
sx q[3];
rz(-1.0704853) q[3];
sx q[3];
rz(-2.0557192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0367947) q[0];
sx q[0];
rz(-3.1322271) q[0];
sx q[0];
rz(-2.985756) q[0];
rz(1.8644631) q[1];
sx q[1];
rz(-1.3469478) q[1];
sx q[1];
rz(3.1256622) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6496507) q[0];
sx q[0];
rz(-1.3306405) q[0];
sx q[0];
rz(-2.35046) q[0];
x q[1];
rz(0.35959776) q[2];
sx q[2];
rz(-0.82530515) q[2];
sx q[2];
rz(-1.2957089) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5585352) q[1];
sx q[1];
rz(-2.0777521) q[1];
sx q[1];
rz(2.3304549) q[1];
rz(-pi) q[2];
rz(-2.5032477) q[3];
sx q[3];
rz(-1.2223635) q[3];
sx q[3];
rz(1.750647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5126123) q[2];
sx q[2];
rz(-0.11233687) q[2];
sx q[2];
rz(-3.0158499) q[2];
rz(0.9564774) q[3];
sx q[3];
rz(-1.4230909) q[3];
sx q[3];
rz(-0.3392578) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9823031) q[0];
sx q[0];
rz(-0.24022261) q[0];
sx q[0];
rz(-2.6851728) q[0];
rz(1.5688815) q[1];
sx q[1];
rz(-2.1379037) q[1];
sx q[1];
rz(2.454954) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4149041) q[0];
sx q[0];
rz(-0.81302887) q[0];
sx q[0];
rz(-2.4343879) q[0];
rz(-pi) q[1];
rz(-1.3977658) q[2];
sx q[2];
rz(-1.4464207) q[2];
sx q[2];
rz(1.591452) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7788199) q[1];
sx q[1];
rz(-1.6906926) q[1];
sx q[1];
rz(2.9935548) q[1];
rz(-pi) q[2];
x q[2];
rz(0.17466361) q[3];
sx q[3];
rz(-1.5015937) q[3];
sx q[3];
rz(-0.098002794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3829284) q[2];
sx q[2];
rz(-1.1229346) q[2];
sx q[2];
rz(2.0056966) q[2];
rz(2.1618333) q[3];
sx q[3];
rz(-1.3330678) q[3];
sx q[3];
rz(-2.4433344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6887688) q[0];
sx q[0];
rz(-1.3249506) q[0];
sx q[0];
rz(-2.685637) q[0];
rz(-3.0988354) q[1];
sx q[1];
rz(-1.9330934) q[1];
sx q[1];
rz(-0.6932238) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42825481) q[0];
sx q[0];
rz(-1.9549668) q[0];
sx q[0];
rz(-0.1066271) q[0];
rz(-pi) q[1];
rz(-1.9749538) q[2];
sx q[2];
rz(-1.9459138) q[2];
sx q[2];
rz(2.1585495) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0574368) q[1];
sx q[1];
rz(-0.098271253) q[1];
sx q[1];
rz(-1.4986503) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59488036) q[3];
sx q[3];
rz(-0.52973807) q[3];
sx q[3];
rz(1.6486419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2263055) q[2];
sx q[2];
rz(-1.8733571) q[2];
sx q[2];
rz(0.66429663) q[2];
rz(1.9281467) q[3];
sx q[3];
rz(-0.72187859) q[3];
sx q[3];
rz(-2.2693995) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.450347) q[0];
sx q[0];
rz(-0.84083122) q[0];
sx q[0];
rz(-1.3837411) q[0];
rz(-1.4585523) q[1];
sx q[1];
rz(-2.8589307) q[1];
sx q[1];
rz(0.9160441) q[1];
rz(3.0677879) q[2];
sx q[2];
rz(-2.7701785) q[2];
sx q[2];
rz(-2.2899173) q[2];
rz(-2.5057143) q[3];
sx q[3];
rz(-2.3273031) q[3];
sx q[3];
rz(0.81958156) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
