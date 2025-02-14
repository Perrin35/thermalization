OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.423288) q[0];
sx q[0];
rz(-2.365132) q[0];
sx q[0];
rz(-2.7756696) q[0];
rz(-2.794682) q[1];
sx q[1];
rz(-2.8568755) q[1];
sx q[1];
rz(1.3887583) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.552504) q[0];
sx q[0];
rz(-2.9285746) q[0];
sx q[0];
rz(-2.8211843) q[0];
rz(0.69268815) q[2];
sx q[2];
rz(-1.8028304) q[2];
sx q[2];
rz(-1.9197861) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5947111) q[1];
sx q[1];
rz(-1.1818019) q[1];
sx q[1];
rz(2.0531897) q[1];
x q[2];
rz(-0.59361122) q[3];
sx q[3];
rz(-2.0214012) q[3];
sx q[3];
rz(-3.1385147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.95691386) q[2];
sx q[2];
rz(-0.92014402) q[2];
sx q[2];
rz(-2.206395) q[2];
rz(2.8397371) q[3];
sx q[3];
rz(-1.5993092) q[3];
sx q[3];
rz(-2.7479318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
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
rz(-0.97717706) q[0];
sx q[0];
rz(-1.2657607) q[0];
sx q[0];
rz(2.6432977) q[0];
rz(-1.7138819) q[1];
sx q[1];
rz(-1.2063113) q[1];
sx q[1];
rz(2.0548342) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0181429) q[0];
sx q[0];
rz(-0.88844127) q[0];
sx q[0];
rz(1.2641973) q[0];
rz(-pi) q[1];
rz(2.9478023) q[2];
sx q[2];
rz(-1.8326129) q[2];
sx q[2];
rz(1.7270391) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.77176911) q[1];
sx q[1];
rz(-0.46716094) q[1];
sx q[1];
rz(-2.7792756) q[1];
rz(-1.5828176) q[3];
sx q[3];
rz(-1.88899) q[3];
sx q[3];
rz(0.98230045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.19865092) q[2];
sx q[2];
rz(-2.8309839) q[2];
sx q[2];
rz(-1.3322213) q[2];
rz(1.9940935) q[3];
sx q[3];
rz(-1.56286) q[3];
sx q[3];
rz(1.8089627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-0.20951095) q[0];
sx q[0];
rz(-1.4028343) q[0];
sx q[0];
rz(0.15723666) q[0];
rz(1.9137742) q[1];
sx q[1];
rz(-0.20918748) q[1];
sx q[1];
rz(-0.42207178) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7868766) q[0];
sx q[0];
rz(-1.7484957) q[0];
sx q[0];
rz(3.0502158) q[0];
x q[1];
rz(0.55251212) q[2];
sx q[2];
rz(-0.32304155) q[2];
sx q[2];
rz(-0.99977899) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6032519) q[1];
sx q[1];
rz(-1.1900239) q[1];
sx q[1];
rz(-0.095007665) q[1];
rz(-1.8881599) q[3];
sx q[3];
rz(-0.79447132) q[3];
sx q[3];
rz(-1.0004071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.813039) q[2];
sx q[2];
rz(-2.1812794) q[2];
sx q[2];
rz(-0.73406827) q[2];
rz(2.0640533) q[3];
sx q[3];
rz(-2.4534093) q[3];
sx q[3];
rz(-3.0663826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5828736) q[0];
sx q[0];
rz(-2.131077) q[0];
sx q[0];
rz(-3.1251113) q[0];
rz(-0.074542848) q[1];
sx q[1];
rz(-0.45769474) q[1];
sx q[1];
rz(-0.25513908) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6612084) q[0];
sx q[0];
rz(-1.0912278) q[0];
sx q[0];
rz(-1.40856) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8066977) q[2];
sx q[2];
rz(-2.1112295) q[2];
sx q[2];
rz(1.188736) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.29683477) q[1];
sx q[1];
rz(-0.95176178) q[1];
sx q[1];
rz(-2.1059787) q[1];
rz(-1.9202407) q[3];
sx q[3];
rz(-1.6779052) q[3];
sx q[3];
rz(-0.28863321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5269346) q[2];
sx q[2];
rz(-0.69710985) q[2];
sx q[2];
rz(-2.441791) q[2];
rz(-3.0497293) q[3];
sx q[3];
rz(-1.9242761) q[3];
sx q[3];
rz(-0.45472586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5192473) q[0];
sx q[0];
rz(-0.82038251) q[0];
sx q[0];
rz(2.0297594) q[0];
rz(-2.9513997) q[1];
sx q[1];
rz(-0.88514248) q[1];
sx q[1];
rz(-1.9817188) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1721508) q[0];
sx q[0];
rz(-1.6052725) q[0];
sx q[0];
rz(1.9146754) q[0];
rz(0.93244035) q[2];
sx q[2];
rz(-1.4529422) q[2];
sx q[2];
rz(-1.7913851) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.24007356) q[1];
sx q[1];
rz(-1.0283955) q[1];
sx q[1];
rz(1.8611004) q[1];
rz(2.0926863) q[3];
sx q[3];
rz(-1.6097704) q[3];
sx q[3];
rz(3.1119973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5708892) q[2];
sx q[2];
rz(-0.8231701) q[2];
sx q[2];
rz(-2.9111351) q[2];
rz(2.0795836) q[3];
sx q[3];
rz(-1.2758723) q[3];
sx q[3];
rz(1.6331858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(2.688711) q[0];
sx q[0];
rz(-2.0827561) q[0];
sx q[0];
rz(-1.3558615) q[0];
rz(0.52976766) q[1];
sx q[1];
rz(-2.3593088) q[1];
sx q[1];
rz(-1.9897602) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7811711) q[0];
sx q[0];
rz(-2.7716405) q[0];
sx q[0];
rz(-2.8656703) q[0];
rz(-pi) q[1];
rz(0.095345796) q[2];
sx q[2];
rz(-2.110425) q[2];
sx q[2];
rz(-2.6971872) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5024912) q[1];
sx q[1];
rz(-1.7646953) q[1];
sx q[1];
rz(0.37096937) q[1];
rz(-0.21422503) q[3];
sx q[3];
rz(-1.4559295) q[3];
sx q[3];
rz(-1.3319931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4749703) q[2];
sx q[2];
rz(-0.64133659) q[2];
sx q[2];
rz(2.6744911) q[2];
rz(-2.1304255) q[3];
sx q[3];
rz(-2.7979388) q[3];
sx q[3];
rz(-1.3694192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2924627) q[0];
sx q[0];
rz(-2.9865773) q[0];
sx q[0];
rz(-2.9169061) q[0];
rz(0.16381964) q[1];
sx q[1];
rz(-2.8024709) q[1];
sx q[1];
rz(-2.4952369) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3098329) q[0];
sx q[0];
rz(-1.0467967) q[0];
sx q[0];
rz(2.8519467) q[0];
rz(-pi) q[1];
rz(-0.21804131) q[2];
sx q[2];
rz(-1.8865117) q[2];
sx q[2];
rz(1.9670847) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.028370628) q[1];
sx q[1];
rz(-2.592855) q[1];
sx q[1];
rz(-2.996925) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.32352738) q[3];
sx q[3];
rz(-2.7351742) q[3];
sx q[3];
rz(-0.42662963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9859163) q[2];
sx q[2];
rz(-1.4480269) q[2];
sx q[2];
rz(-2.7316366) q[2];
rz(2.3112467) q[3];
sx q[3];
rz(-0.94873077) q[3];
sx q[3];
rz(1.4478987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53772563) q[0];
sx q[0];
rz(-0.69320885) q[0];
sx q[0];
rz(-2.895288) q[0];
rz(-2.314997) q[1];
sx q[1];
rz(-1.7431424) q[1];
sx q[1];
rz(-0.26396096) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0279585) q[0];
sx q[0];
rz(-1.6662681) q[0];
sx q[0];
rz(-0.0072633538) q[0];
rz(-2.6648681) q[2];
sx q[2];
rz(-2.2988604) q[2];
sx q[2];
rz(-1.6260894) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5689478) q[1];
sx q[1];
rz(-1.8518475) q[1];
sx q[1];
rz(0.87139178) q[1];
rz(-pi) q[2];
rz(-0.58707763) q[3];
sx q[3];
rz(-1.2012606) q[3];
sx q[3];
rz(1.0472421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0366514) q[2];
sx q[2];
rz(-1.2678601) q[2];
sx q[2];
rz(-0.71511739) q[2];
rz(-3.0702843) q[3];
sx q[3];
rz(-1.5846059) q[3];
sx q[3];
rz(-0.74688545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0095373) q[0];
sx q[0];
rz(-1.6522464) q[0];
sx q[0];
rz(-0.43553964) q[0];
rz(0.14777331) q[1];
sx q[1];
rz(-1.7099893) q[1];
sx q[1];
rz(1.4685644) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4184121) q[0];
sx q[0];
rz(-1.7505129) q[0];
sx q[0];
rz(-2.7030858) q[0];
rz(1.4917489) q[2];
sx q[2];
rz(-2.2666396) q[2];
sx q[2];
rz(-0.33338132) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5295386) q[1];
sx q[1];
rz(-1.4040605) q[1];
sx q[1];
rz(-1.177854) q[1];
x q[2];
rz(2.6665523) q[3];
sx q[3];
rz(-2.6938525) q[3];
sx q[3];
rz(-2.077075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.396951) q[2];
sx q[2];
rz(-1.3775185) q[2];
sx q[2];
rz(-0.89293876) q[2];
rz(-1.2785771) q[3];
sx q[3];
rz(-1.1568926) q[3];
sx q[3];
rz(1.4634092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8466012) q[0];
sx q[0];
rz(-2.8104267) q[0];
sx q[0];
rz(2.5355329) q[0];
rz(0.010295708) q[1];
sx q[1];
rz(-0.98774397) q[1];
sx q[1];
rz(-0.079924718) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.851136) q[0];
sx q[0];
rz(-2.5154167) q[0];
sx q[0];
rz(2.0100223) q[0];
x q[1];
rz(-2.299231) q[2];
sx q[2];
rz(-1.3536705) q[2];
sx q[2];
rz(-0.44936839) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8070302) q[1];
sx q[1];
rz(-1.1638146) q[1];
sx q[1];
rz(-1.5466067) q[1];
x q[2];
rz(-2.2977423) q[3];
sx q[3];
rz(-2.8371713) q[3];
sx q[3];
rz(-2.3850887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.14110485) q[2];
sx q[2];
rz(-0.92685574) q[2];
sx q[2];
rz(1.0151939) q[2];
rz(0.60123932) q[3];
sx q[3];
rz(-0.70178086) q[3];
sx q[3];
rz(-2.0950441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64423185) q[0];
sx q[0];
rz(-1.624122) q[0];
sx q[0];
rz(-2.4540785) q[0];
rz(2.5524706) q[1];
sx q[1];
rz(-2.4644869) q[1];
sx q[1];
rz(-0.45225515) q[1];
rz(1.6670139) q[2];
sx q[2];
rz(-2.0141891) q[2];
sx q[2];
rz(-0.63465848) q[2];
rz(1.6986871) q[3];
sx q[3];
rz(-2.3913132) q[3];
sx q[3];
rz(-2.0146418) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
