OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.93513918) q[0];
sx q[0];
rz(-2.3606665) q[0];
sx q[0];
rz(0.20679064) q[0];
rz(0.38987723) q[1];
sx q[1];
rz(-1.0607399) q[1];
sx q[1];
rz(-0.18016711) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59158303) q[0];
sx q[0];
rz(-1.9073434) q[0];
sx q[0];
rz(-2.6495289) q[0];
rz(-1.4397058) q[2];
sx q[2];
rz(-1.0616454) q[2];
sx q[2];
rz(-2.5703562) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.7000007) q[1];
sx q[1];
rz(-2.4543426) q[1];
sx q[1];
rz(1.6874466) q[1];
x q[2];
rz(1.3869242) q[3];
sx q[3];
rz(-0.89727083) q[3];
sx q[3];
rz(2.9415188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7303598) q[2];
sx q[2];
rz(-2.2683472) q[2];
sx q[2];
rz(-1.8475378) q[2];
rz(-2.7358352) q[3];
sx q[3];
rz(-1.5016705) q[3];
sx q[3];
rz(-2.7348203) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2186573) q[0];
sx q[0];
rz(-1.0359534) q[0];
sx q[0];
rz(0.12582114) q[0];
rz(-0.80548349) q[1];
sx q[1];
rz(-0.80635726) q[1];
sx q[1];
rz(1.3719826) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67324154) q[0];
sx q[0];
rz(-1.5090319) q[0];
sx q[0];
rz(-2.4657235) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4470909) q[2];
sx q[2];
rz(-2.1134085) q[2];
sx q[2];
rz(-1.6306842) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6276715) q[1];
sx q[1];
rz(-1.5253592) q[1];
sx q[1];
rz(0.25000484) q[1];
x q[2];
rz(-2.6729229) q[3];
sx q[3];
rz(-1.3523003) q[3];
sx q[3];
rz(-3.1359429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8877318) q[2];
sx q[2];
rz(-2.4527206) q[2];
sx q[2];
rz(0.99622336) q[2];
rz(2.0455202) q[3];
sx q[3];
rz(-1.6888065) q[3];
sx q[3];
rz(0.85038275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4662194) q[0];
sx q[0];
rz(-0.96005625) q[0];
sx q[0];
rz(2.0825785) q[0];
rz(-1.1478708) q[1];
sx q[1];
rz(-2.4025326) q[1];
sx q[1];
rz(-1.0645197) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6305884) q[0];
sx q[0];
rz(-0.9398886) q[0];
sx q[0];
rz(2.0677807) q[0];
rz(-1.1946482) q[2];
sx q[2];
rz(-2.7333439) q[2];
sx q[2];
rz(-0.86366913) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3036084) q[1];
sx q[1];
rz(-1.260699) q[1];
sx q[1];
rz(3.1293948) q[1];
rz(-pi) q[2];
rz(-0.84633175) q[3];
sx q[3];
rz(-1.6958106) q[3];
sx q[3];
rz(2.4785329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0718096) q[2];
sx q[2];
rz(-1.1867384) q[2];
sx q[2];
rz(0.26322571) q[2];
rz(1.1188544) q[3];
sx q[3];
rz(-1.8656732) q[3];
sx q[3];
rz(0.81937218) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5082821) q[0];
sx q[0];
rz(-0.044697035) q[0];
sx q[0];
rz(1.7472965) q[0];
rz(2.1030203) q[1];
sx q[1];
rz(-1.4515406) q[1];
sx q[1];
rz(-2.0746453) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2571714) q[0];
sx q[0];
rz(-2.1002249) q[0];
sx q[0];
rz(3.0966395) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6520086) q[2];
sx q[2];
rz(-1.8394543) q[2];
sx q[2];
rz(0.072222885) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3921515) q[1];
sx q[1];
rz(-1.1679808) q[1];
sx q[1];
rz(-0.72026003) q[1];
rz(-pi) q[2];
rz(-1.0889441) q[3];
sx q[3];
rz(-1.370016) q[3];
sx q[3];
rz(2.5364385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2944494) q[2];
sx q[2];
rz(-1.7913982) q[2];
sx q[2];
rz(2.8520544) q[2];
rz(-2.6211522) q[3];
sx q[3];
rz(-1.3402904) q[3];
sx q[3];
rz(2.382544) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6510058) q[0];
sx q[0];
rz(-3.1327972) q[0];
sx q[0];
rz(-2.3068413) q[0];
rz(-1.897215) q[1];
sx q[1];
rz(-1.8908187) q[1];
sx q[1];
rz(1.7117737) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1236876) q[0];
sx q[0];
rz(-3.078853) q[0];
sx q[0];
rz(-0.040266589) q[0];
x q[1];
rz(-1.3313815) q[2];
sx q[2];
rz(-1.8481701) q[2];
sx q[2];
rz(-1.1053567) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5606219) q[1];
sx q[1];
rz(-2.3498658) q[1];
sx q[1];
rz(-0.03246275) q[1];
x q[2];
rz(-0.17794869) q[3];
sx q[3];
rz(-2.2214409) q[3];
sx q[3];
rz(-1.4072756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.66118583) q[2];
sx q[2];
rz(-1.0057665) q[2];
sx q[2];
rz(-1.8583813) q[2];
rz(3.0094106) q[3];
sx q[3];
rz(-0.32871267) q[3];
sx q[3];
rz(-2.8018518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4955687) q[0];
sx q[0];
rz(-3.0806354) q[0];
sx q[0];
rz(-2.6640889) q[0];
rz(1.5006789) q[1];
sx q[1];
rz(-1.532282) q[1];
sx q[1];
rz(2.5710411) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7084224) q[0];
sx q[0];
rz(-2.6714984) q[0];
sx q[0];
rz(0.8870468) q[0];
x q[1];
rz(0.678755) q[2];
sx q[2];
rz(-0.55870134) q[2];
sx q[2];
rz(0.15685454) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4482566) q[1];
sx q[1];
rz(-0.92474557) q[1];
sx q[1];
rz(2.4962884) q[1];
rz(-1.5930575) q[3];
sx q[3];
rz(-0.92261693) q[3];
sx q[3];
rz(-2.9472586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.70696124) q[2];
sx q[2];
rz(-1.0501477) q[2];
sx q[2];
rz(2.3449507) q[2];
rz(-2.8213275) q[3];
sx q[3];
rz(-2.0551149) q[3];
sx q[3];
rz(-1.2789352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0657848) q[0];
sx q[0];
rz(-1.0362352) q[0];
sx q[0];
rz(2.7835223) q[0];
rz(0.2886731) q[1];
sx q[1];
rz(-2.6117548) q[1];
sx q[1];
rz(-1.8310865) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3863556) q[0];
sx q[0];
rz(-2.1023395) q[0];
sx q[0];
rz(1.3528354) q[0];
rz(0.99408044) q[2];
sx q[2];
rz(-2.1834282) q[2];
sx q[2];
rz(2.6660369) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9999034) q[1];
sx q[1];
rz(-1.1232166) q[1];
sx q[1];
rz(-2.5738641) q[1];
rz(-pi) q[2];
rz(0.86705039) q[3];
sx q[3];
rz(-2.3208445) q[3];
sx q[3];
rz(-1.5914608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.33401176) q[2];
sx q[2];
rz(-0.58471218) q[2];
sx q[2];
rz(-2.1441148) q[2];
rz(0.57972646) q[3];
sx q[3];
rz(-2.4119191) q[3];
sx q[3];
rz(1.4355481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8758133) q[0];
sx q[0];
rz(-1.6515323) q[0];
sx q[0];
rz(0.34061256) q[0];
rz(1.2365201) q[1];
sx q[1];
rz(-0.75729901) q[1];
sx q[1];
rz(1.1901201) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0252359) q[0];
sx q[0];
rz(-1.3394757) q[0];
sx q[0];
rz(3.0781834) q[0];
rz(-pi) q[1];
rz(-0.068433381) q[2];
sx q[2];
rz(-1.6514196) q[2];
sx q[2];
rz(1.2863408) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3548673) q[1];
sx q[1];
rz(-0.39199542) q[1];
sx q[1];
rz(-0.74704945) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.45300071) q[3];
sx q[3];
rz(-2.0815597) q[3];
sx q[3];
rz(2.9651027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3999346) q[2];
sx q[2];
rz(-0.88230336) q[2];
sx q[2];
rz(2.7271872) q[2];
rz(1.7587781) q[3];
sx q[3];
rz(-1.4004935) q[3];
sx q[3];
rz(-2.8167021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25093108) q[0];
sx q[0];
rz(-1.0216167) q[0];
sx q[0];
rz(-0.1517621) q[0];
rz(-1.7157308) q[1];
sx q[1];
rz(-1.3769923) q[1];
sx q[1];
rz(2.192416) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.855809) q[0];
sx q[0];
rz(-1.4306418) q[0];
sx q[0];
rz(-2.6967718) q[0];
rz(-pi) q[1];
rz(-2.8434535) q[2];
sx q[2];
rz(-2.4900644) q[2];
sx q[2];
rz(2.8704314) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6712499) q[1];
sx q[1];
rz(-1.3312695) q[1];
sx q[1];
rz(2.4688979) q[1];
rz(-pi) q[2];
rz(0.98032326) q[3];
sx q[3];
rz(-1.983641) q[3];
sx q[3];
rz(-1.0124029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40733797) q[2];
sx q[2];
rz(-2.7463425) q[2];
sx q[2];
rz(-0.43241832) q[2];
rz(-1.5405103) q[3];
sx q[3];
rz(-1.5099022) q[3];
sx q[3];
rz(0.66175118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0703053) q[0];
sx q[0];
rz(-2.8685331) q[0];
sx q[0];
rz(2.7684257) q[0];
rz(2.7846653) q[1];
sx q[1];
rz(-2.280805) q[1];
sx q[1];
rz(2.4180791) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0821738) q[0];
sx q[0];
rz(-2.2954303) q[0];
sx q[0];
rz(1.7436149) q[0];
x q[1];
rz(1.6610442) q[2];
sx q[2];
rz(-0.89333488) q[2];
sx q[2];
rz(-1.8842763) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3620421) q[1];
sx q[1];
rz(-2.3099646) q[1];
sx q[1];
rz(-1.9401624) q[1];
rz(0.017541842) q[3];
sx q[3];
rz(-2.8711257) q[3];
sx q[3];
rz(1.5568352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2202806) q[2];
sx q[2];
rz(-1.7582515) q[2];
sx q[2];
rz(2.5881361) q[2];
rz(-2.4297595) q[3];
sx q[3];
rz(-2.7332941) q[3];
sx q[3];
rz(-2.0994983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.2198467) q[0];
sx q[0];
rz(-2.5384359) q[0];
sx q[0];
rz(0.13721101) q[0];
rz(2.8593821) q[1];
sx q[1];
rz(-2.3574775) q[1];
sx q[1];
rz(-2.2025253) q[1];
rz(-2.6454906) q[2];
sx q[2];
rz(-1.9719057) q[2];
sx q[2];
rz(3.0531648) q[2];
rz(-0.59755748) q[3];
sx q[3];
rz(-0.47511027) q[3];
sx q[3];
rz(-2.6078754) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
