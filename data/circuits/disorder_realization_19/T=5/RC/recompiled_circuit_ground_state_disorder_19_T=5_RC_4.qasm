OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8813397) q[0];
sx q[0];
rz(-0.94085675) q[0];
sx q[0];
rz(-0.22766222) q[0];
rz(-2.8582299) q[1];
sx q[1];
rz(-0.41937399) q[1];
sx q[1];
rz(1.286932) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44449466) q[0];
sx q[0];
rz(-2.6711406) q[0];
sx q[0];
rz(1.4509499) q[0];
rz(-pi) q[1];
rz(-2.3782733) q[2];
sx q[2];
rz(-0.60930353) q[2];
sx q[2];
rz(1.0805902) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.70482774) q[1];
sx q[1];
rz(-2.4840762) q[1];
sx q[1];
rz(2.2566811) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8836796) q[3];
sx q[3];
rz(-2.1863424) q[3];
sx q[3];
rz(-2.2039977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8925573) q[2];
sx q[2];
rz(-2.1302569) q[2];
sx q[2];
rz(2.2300143) q[2];
rz(0.75561953) q[3];
sx q[3];
rz(-2.8204462) q[3];
sx q[3];
rz(0.49629456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3590473) q[0];
sx q[0];
rz(-1.3107212) q[0];
sx q[0];
rz(-2.0027335) q[0];
rz(-2.5149939) q[1];
sx q[1];
rz(-0.36830026) q[1];
sx q[1];
rz(1.0777333) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2127578) q[0];
sx q[0];
rz(-0.99883119) q[0];
sx q[0];
rz(-1.1642745) q[0];
rz(2.0135278) q[2];
sx q[2];
rz(-1.2473543) q[2];
sx q[2];
rz(1.5087104) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.30451286) q[1];
sx q[1];
rz(-2.8690352) q[1];
sx q[1];
rz(-2.8902092) q[1];
rz(-pi) q[2];
rz(0.84799453) q[3];
sx q[3];
rz(-1.2395596) q[3];
sx q[3];
rz(0.02324638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6935912) q[2];
sx q[2];
rz(-2.3841264) q[2];
sx q[2];
rz(2.9288911) q[2];
rz(1.5851783) q[3];
sx q[3];
rz(-2.3978265) q[3];
sx q[3];
rz(-2.0785544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14146516) q[0];
sx q[0];
rz(-0.26697049) q[0];
sx q[0];
rz(2.2947327) q[0];
rz(2.8894539) q[1];
sx q[1];
rz(-0.82770258) q[1];
sx q[1];
rz(-2.491378) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2144521) q[0];
sx q[0];
rz(-0.34159476) q[0];
sx q[0];
rz(0.13169698) q[0];
rz(-0.79651259) q[2];
sx q[2];
rz(-1.0625417) q[2];
sx q[2];
rz(0.69895335) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6628428) q[1];
sx q[1];
rz(-1.2744941) q[1];
sx q[1];
rz(-3.0257312) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.43127755) q[3];
sx q[3];
rz(-0.85429885) q[3];
sx q[3];
rz(0.34654472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.83027679) q[2];
sx q[2];
rz(-2.4429784) q[2];
sx q[2];
rz(2.1777731) q[2];
rz(-1.3612932) q[3];
sx q[3];
rz(-2.6908974) q[3];
sx q[3];
rz(-0.042958766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2394102) q[0];
sx q[0];
rz(-0.22265156) q[0];
sx q[0];
rz(-1.0182925) q[0];
rz(2.2564383) q[1];
sx q[1];
rz(-2.8926909) q[1];
sx q[1];
rz(2.8525066) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24482306) q[0];
sx q[0];
rz(-0.45216783) q[0];
sx q[0];
rz(2.3954579) q[0];
rz(0.54479213) q[2];
sx q[2];
rz(-0.71342403) q[2];
sx q[2];
rz(-2.1868844) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4319181) q[1];
sx q[1];
rz(-2.5920715) q[1];
sx q[1];
rz(1.2195682) q[1];
rz(0.21896514) q[3];
sx q[3];
rz(-1.4462399) q[3];
sx q[3];
rz(1.7809465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.79124147) q[2];
sx q[2];
rz(-1.6356607) q[2];
sx q[2];
rz(1.2887456) q[2];
rz(-0.57981235) q[3];
sx q[3];
rz(-2.6668187) q[3];
sx q[3];
rz(0.60540664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8346005) q[0];
sx q[0];
rz(-2.6267316) q[0];
sx q[0];
rz(-0.32421625) q[0];
rz(3.1027555) q[1];
sx q[1];
rz(-2.4034998) q[1];
sx q[1];
rz(1.0968346) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7997579) q[0];
sx q[0];
rz(-1.2909162) q[0];
sx q[0];
rz(1.3169692) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7261502) q[2];
sx q[2];
rz(-2.5553779) q[2];
sx q[2];
rz(-1.4582576) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.13584863) q[1];
sx q[1];
rz(-1.6127819) q[1];
sx q[1];
rz(-1.2688925) q[1];
rz(1.9911258) q[3];
sx q[3];
rz(-1.4864576) q[3];
sx q[3];
rz(-1.4140417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.066112) q[2];
sx q[2];
rz(-2.8754063) q[2];
sx q[2];
rz(-3.10293) q[2];
rz(-1.4085116) q[3];
sx q[3];
rz(-1.446529) q[3];
sx q[3];
rz(-0.2814289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-2.6283145) q[0];
sx q[0];
rz(-2.3230041) q[0];
sx q[0];
rz(1.0191089) q[0];
rz(0.45711532) q[1];
sx q[1];
rz(-0.26908427) q[1];
sx q[1];
rz(1.4310744) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5810878) q[0];
sx q[0];
rz(-1.771534) q[0];
sx q[0];
rz(1.5439347) q[0];
rz(-0.98177399) q[2];
sx q[2];
rz(-1.7830689) q[2];
sx q[2];
rz(-2.6350217) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7841815) q[1];
sx q[1];
rz(-2.7704382) q[1];
sx q[1];
rz(1.3264912) q[1];
rz(-pi) q[2];
rz(0.39737153) q[3];
sx q[3];
rz(-2.2512967) q[3];
sx q[3];
rz(1.3272663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9223601) q[2];
sx q[2];
rz(-1.7722426) q[2];
sx q[2];
rz(2.5617981) q[2];
rz(-0.29948768) q[3];
sx q[3];
rz(-0.53286415) q[3];
sx q[3];
rz(1.2286435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2806468) q[0];
sx q[0];
rz(-1.6609284) q[0];
sx q[0];
rz(-3.0378367) q[0];
rz(-2.5170028) q[1];
sx q[1];
rz(-2.2207405) q[1];
sx q[1];
rz(1.7594899) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1555903) q[0];
sx q[0];
rz(-2.3891923) q[0];
sx q[0];
rz(1.8852021) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29989179) q[2];
sx q[2];
rz(-1.6768528) q[2];
sx q[2];
rz(2.7404355) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.15006062) q[1];
sx q[1];
rz(-1.0673017) q[1];
sx q[1];
rz(0.30725422) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4806466) q[3];
sx q[3];
rz(-1.7978284) q[3];
sx q[3];
rz(1.1477136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3524807) q[2];
sx q[2];
rz(-1.666297) q[2];
sx q[2];
rz(0.86461198) q[2];
rz(-0.23884808) q[3];
sx q[3];
rz(-0.75967234) q[3];
sx q[3];
rz(0.71340942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9637941) q[0];
sx q[0];
rz(-2.5803784) q[0];
sx q[0];
rz(-0.23502769) q[0];
rz(-2.4759953) q[1];
sx q[1];
rz(-2.0033629) q[1];
sx q[1];
rz(-1.6682909) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.145133) q[0];
sx q[0];
rz(-2.1204824) q[0];
sx q[0];
rz(0.24497801) q[0];
rz(-1.4334045) q[2];
sx q[2];
rz(-1.1072888) q[2];
sx q[2];
rz(-2.3799294) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1531495) q[1];
sx q[1];
rz(-1.4979939) q[1];
sx q[1];
rz(-0.18393391) q[1];
x q[2];
rz(-1.6157916) q[3];
sx q[3];
rz(-1.1200841) q[3];
sx q[3];
rz(1.352528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.16234806) q[2];
sx q[2];
rz(-2.9408216) q[2];
sx q[2];
rz(0.46094224) q[2];
rz(2.0567242) q[3];
sx q[3];
rz(-0.74551398) q[3];
sx q[3];
rz(-2.357024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.064086176) q[0];
sx q[0];
rz(-0.53090799) q[0];
sx q[0];
rz(2.532646) q[0];
rz(-0.56600904) q[1];
sx q[1];
rz(-1.9809664) q[1];
sx q[1];
rz(1.9401248) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4964551) q[0];
sx q[0];
rz(-1.7063396) q[0];
sx q[0];
rz(0.0082502967) q[0];
rz(-1.9502152) q[2];
sx q[2];
rz(-2.622329) q[2];
sx q[2];
rz(1.0644703) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2347327) q[1];
sx q[1];
rz(-0.77488778) q[1];
sx q[1];
rz(-3.0084684) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5646832) q[3];
sx q[3];
rz(-2.6196777) q[3];
sx q[3];
rz(-2.1186704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.10298771) q[2];
sx q[2];
rz(-1.0706341) q[2];
sx q[2];
rz(-1.9496244) q[2];
rz(-2.1206756) q[3];
sx q[3];
rz(-2.9396785) q[3];
sx q[3];
rz(1.6010223) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9182619) q[0];
sx q[0];
rz(-0.20063618) q[0];
sx q[0];
rz(0.9675135) q[0];
rz(-1.4406904) q[1];
sx q[1];
rz(-0.43165019) q[1];
sx q[1];
rz(-0.38223019) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8865693) q[0];
sx q[0];
rz(-2.5394727) q[0];
sx q[0];
rz(1.894879) q[0];
x q[1];
rz(-2.380314) q[2];
sx q[2];
rz(-2.571311) q[2];
sx q[2];
rz(2.7340661) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6401419) q[1];
sx q[1];
rz(-2.3489408) q[1];
sx q[1];
rz(-0.16736302) q[1];
x q[2];
rz(-1.554364) q[3];
sx q[3];
rz(-2.1167618) q[3];
sx q[3];
rz(0.16758238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.050015673) q[2];
sx q[2];
rz(-2.2735368) q[2];
sx q[2];
rz(-2.3559605) q[2];
rz(3.1146289) q[3];
sx q[3];
rz(-1.6227159) q[3];
sx q[3];
rz(-0.54076076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80326573) q[0];
sx q[0];
rz(-1.6852408) q[0];
sx q[0];
rz(2.2483873) q[0];
rz(2.4772353) q[1];
sx q[1];
rz(-1.2394445) q[1];
sx q[1];
rz(2.209421) q[1];
rz(1.8705838) q[2];
sx q[2];
rz(-2.2769417) q[2];
sx q[2];
rz(-1.4852038) q[2];
rz(-1.5803278) q[3];
sx q[3];
rz(-0.96612488) q[3];
sx q[3];
rz(0.92074367) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
