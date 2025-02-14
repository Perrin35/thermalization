OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.779939) q[0];
sx q[0];
rz(-0.10804636) q[0];
sx q[0];
rz(-1.2471696) q[0];
rz(-0.78020686) q[1];
sx q[1];
rz(-1.127004) q[1];
sx q[1];
rz(2.3957774) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5597668) q[0];
sx q[0];
rz(-0.39224658) q[0];
sx q[0];
rz(2.0779209) q[0];
x q[1];
rz(1.5175061) q[2];
sx q[2];
rz(-0.16096965) q[2];
sx q[2];
rz(-0.56426891) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8565897) q[1];
sx q[1];
rz(-1.3599281) q[1];
sx q[1];
rz(-0.1182571) q[1];
rz(-2.8549544) q[3];
sx q[3];
rz(-0.53813808) q[3];
sx q[3];
rz(1.9236444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9790393) q[2];
sx q[2];
rz(-0.59430846) q[2];
sx q[2];
rz(1.7621367) q[2];
rz(0.72921324) q[3];
sx q[3];
rz(-0.76494923) q[3];
sx q[3];
rz(-0.929207) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58296975) q[0];
sx q[0];
rz(-2.0780777) q[0];
sx q[0];
rz(-2.9003918) q[0];
rz(-2.8855715) q[1];
sx q[1];
rz(-2.119901) q[1];
sx q[1];
rz(-2.3604438) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.985551) q[0];
sx q[0];
rz(-0.80901481) q[0];
sx q[0];
rz(2.3963905) q[0];
rz(-pi) q[1];
rz(0.56480572) q[2];
sx q[2];
rz(-1.5283268) q[2];
sx q[2];
rz(-0.04893411) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.95805146) q[1];
sx q[1];
rz(-2.1528917) q[1];
sx q[1];
rz(2.2954825) q[1];
rz(-1.7082907) q[3];
sx q[3];
rz(-1.526727) q[3];
sx q[3];
rz(-1.6018683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8476734) q[2];
sx q[2];
rz(-1.5896229) q[2];
sx q[2];
rz(1.5839362) q[2];
rz(3.1368351) q[3];
sx q[3];
rz(-0.59499756) q[3];
sx q[3];
rz(-2.7958272) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7541589) q[0];
sx q[0];
rz(-1.6922981) q[0];
sx q[0];
rz(-0.19101983) q[0];
rz(2.7905131) q[1];
sx q[1];
rz(-2.7103238) q[1];
sx q[1];
rz(1.692159) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3905778) q[0];
sx q[0];
rz(-1.6216941) q[0];
sx q[0];
rz(2.8021332) q[0];
rz(-1.1507785) q[2];
sx q[2];
rz(-1.1598829) q[2];
sx q[2];
rz(1.0622298) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7731157) q[1];
sx q[1];
rz(-2.0252899) q[1];
sx q[1];
rz(1.2746432) q[1];
rz(0.48630096) q[3];
sx q[3];
rz(-2.4913906) q[3];
sx q[3];
rz(-0.40806684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8423975) q[2];
sx q[2];
rz(-1.7123875) q[2];
sx q[2];
rz(-0.12796417) q[2];
rz(1.1658824) q[3];
sx q[3];
rz(-0.88763014) q[3];
sx q[3];
rz(0.33795801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8562451) q[0];
sx q[0];
rz(-2.0576394) q[0];
sx q[0];
rz(-1.4803084) q[0];
rz(3.0604494) q[1];
sx q[1];
rz(-2.0442043) q[1];
sx q[1];
rz(-2.2611484) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1937516) q[0];
sx q[0];
rz(-1.546085) q[0];
sx q[0];
rz(1.5642883) q[0];
x q[1];
rz(3.092406) q[2];
sx q[2];
rz(-2.1890759) q[2];
sx q[2];
rz(-0.82727369) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2374946) q[1];
sx q[1];
rz(-1.6818871) q[1];
sx q[1];
rz(2.6059125) q[1];
rz(-3.1034558) q[3];
sx q[3];
rz(-0.55805579) q[3];
sx q[3];
rz(1.3636774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5389898) q[2];
sx q[2];
rz(-0.67402855) q[2];
sx q[2];
rz(1.7146141) q[2];
rz(-1.9849298) q[3];
sx q[3];
rz(-1.4256698) q[3];
sx q[3];
rz(0.15779933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2222774) q[0];
sx q[0];
rz(-2.3265525) q[0];
sx q[0];
rz(2.2015233) q[0];
rz(2.6433511) q[1];
sx q[1];
rz(-2.2254641) q[1];
sx q[1];
rz(-2.0169651) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8064708) q[0];
sx q[0];
rz(-2.2735519) q[0];
sx q[0];
rz(0.051856144) q[0];
x q[1];
rz(-1.1612879) q[2];
sx q[2];
rz(-1.3826269) q[2];
sx q[2];
rz(1.7766603) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.64943759) q[1];
sx q[1];
rz(-0.78178015) q[1];
sx q[1];
rz(0.63654727) q[1];
rz(3.1203007) q[3];
sx q[3];
rz(-0.7856889) q[3];
sx q[3];
rz(1.1710492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.70790946) q[2];
sx q[2];
rz(-0.83432546) q[2];
sx q[2];
rz(1.4617807) q[2];
rz(2.8953569) q[3];
sx q[3];
rz(-0.88053954) q[3];
sx q[3];
rz(2.06854) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4435302) q[0];
sx q[0];
rz(-1.227523) q[0];
sx q[0];
rz(-2.6865633) q[0];
rz(0.21806923) q[1];
sx q[1];
rz(-2.2360305) q[1];
sx q[1];
rz(-2.6630482) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25898257) q[0];
sx q[0];
rz(-0.6713258) q[0];
sx q[0];
rz(2.9837926) q[0];
rz(-pi) q[1];
rz(1.215082) q[2];
sx q[2];
rz(-2.1797769) q[2];
sx q[2];
rz(-0.019542309) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.45286059) q[1];
sx q[1];
rz(-1.11519) q[1];
sx q[1];
rz(-0.36700008) q[1];
rz(-2.4535937) q[3];
sx q[3];
rz(-2.0562226) q[3];
sx q[3];
rz(-1.7958876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9559418) q[2];
sx q[2];
rz(-1.4943244) q[2];
sx q[2];
rz(2.4556124) q[2];
rz(1.8020804) q[3];
sx q[3];
rz(-1.1687665) q[3];
sx q[3];
rz(1.3919938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5550391) q[0];
sx q[0];
rz(-1.3864484) q[0];
sx q[0];
rz(-2.6626124) q[0];
rz(-0.85887495) q[1];
sx q[1];
rz(-2.5996467) q[1];
sx q[1];
rz(3.0368793) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0901439) q[0];
sx q[0];
rz(-1.3810754) q[0];
sx q[0];
rz(2.7499697) q[0];
rz(-pi) q[1];
rz(-0.48583651) q[2];
sx q[2];
rz(-1.2434894) q[2];
sx q[2];
rz(-0.32020928) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5471902) q[1];
sx q[1];
rz(-2.5195751) q[1];
sx q[1];
rz(-1.1941431) q[1];
x q[2];
rz(2.0488304) q[3];
sx q[3];
rz(-1.5198738) q[3];
sx q[3];
rz(-1.4560521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0844448) q[2];
sx q[2];
rz(-1.8859325) q[2];
sx q[2];
rz(0.88195938) q[2];
rz(-3.0706578) q[3];
sx q[3];
rz(-1.2627914) q[3];
sx q[3];
rz(0.96496636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98081812) q[0];
sx q[0];
rz(-0.47176281) q[0];
sx q[0];
rz(1.2534575) q[0];
rz(2.4313633) q[1];
sx q[1];
rz(-1.1980779) q[1];
sx q[1];
rz(2.0517147) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7701845) q[0];
sx q[0];
rz(-2.29378) q[0];
sx q[0];
rz(0.56130479) q[0];
x q[1];
rz(0.63927357) q[2];
sx q[2];
rz(-0.54080039) q[2];
sx q[2];
rz(1.5098454) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0347669) q[1];
sx q[1];
rz(-2.2112101) q[1];
sx q[1];
rz(0.57042112) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9097611) q[3];
sx q[3];
rz(-0.58887312) q[3];
sx q[3];
rz(-2.8804121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.39500427) q[2];
sx q[2];
rz(-0.85902625) q[2];
sx q[2];
rz(1.0224226) q[2];
rz(-0.59631452) q[3];
sx q[3];
rz(-0.78323451) q[3];
sx q[3];
rz(2.7885126) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31350598) q[0];
sx q[0];
rz(-2.6441898) q[0];
sx q[0];
rz(1.5561546) q[0];
rz(1.6515139) q[1];
sx q[1];
rz(-0.45526344) q[1];
sx q[1];
rz(2.1447287) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5307025) q[0];
sx q[0];
rz(-0.35323745) q[0];
sx q[0];
rz(-0.20521407) q[0];
x q[1];
rz(-1.562122) q[2];
sx q[2];
rz(-2.9563835) q[2];
sx q[2];
rz(2.8767074) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2311923) q[1];
sx q[1];
rz(-2.2072189) q[1];
sx q[1];
rz(-2.5144122) q[1];
rz(-2.6410651) q[3];
sx q[3];
rz(-2.1033035) q[3];
sx q[3];
rz(-0.73126572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.31522125) q[2];
sx q[2];
rz(-2.4085277) q[2];
sx q[2];
rz(0.23615393) q[2];
rz(-2.835623) q[3];
sx q[3];
rz(-1.9491842) q[3];
sx q[3];
rz(1.4845622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21208256) q[0];
sx q[0];
rz(-2.4903553) q[0];
sx q[0];
rz(2.4834852) q[0];
rz(-1.8923538) q[1];
sx q[1];
rz(-0.58964261) q[1];
sx q[1];
rz(-0.49159893) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99557331) q[0];
sx q[0];
rz(-2.4065912) q[0];
sx q[0];
rz(-2.4400178) q[0];
rz(-1.3479606) q[2];
sx q[2];
rz(-1.0165983) q[2];
sx q[2];
rz(-1.8775307) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.39855865) q[1];
sx q[1];
rz(-1.5614206) q[1];
sx q[1];
rz(0.27573632) q[1];
rz(0.082178311) q[3];
sx q[3];
rz(-2.058147) q[3];
sx q[3];
rz(3.0370219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7813985) q[2];
sx q[2];
rz(-2.7935544) q[2];
sx q[2];
rz(-1.4813102) q[2];
rz(0.71634746) q[3];
sx q[3];
rz(-0.64645386) q[3];
sx q[3];
rz(2.9334478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.79138712) q[0];
sx q[0];
rz(-1.5443784) q[0];
sx q[0];
rz(-2.7271893) q[0];
rz(-2.1898337) q[1];
sx q[1];
rz(-1.8487683) q[1];
sx q[1];
rz(1.9602736) q[1];
rz(2.4368481) q[2];
sx q[2];
rz(-0.8316883) q[2];
sx q[2];
rz(-1.4744454) q[2];
rz(2.2496102) q[3];
sx q[3];
rz(-1.4374229) q[3];
sx q[3];
rz(-2.0363804) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
