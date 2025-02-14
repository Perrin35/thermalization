OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.91650668) q[0];
sx q[0];
rz(-2.9986311) q[0];
sx q[0];
rz(1.4885055) q[0];
rz(-2.0090964) q[1];
sx q[1];
rz(-0.43192616) q[1];
sx q[1];
rz(-2.5052524) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7677732) q[0];
sx q[0];
rz(-2.1507906) q[0];
sx q[0];
rz(-2.679002) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1875528) q[2];
sx q[2];
rz(-1.9343209) q[2];
sx q[2];
rz(2.5853047) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.26633599) q[1];
sx q[1];
rz(-0.41558973) q[1];
sx q[1];
rz(-1.8824091) q[1];
x q[2];
rz(-2.2872055) q[3];
sx q[3];
rz(-2.1463089) q[3];
sx q[3];
rz(-2.7137386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1119614) q[2];
sx q[2];
rz(-1.1250857) q[2];
sx q[2];
rz(-0.65043989) q[2];
rz(-1.316831) q[3];
sx q[3];
rz(-1.7951169) q[3];
sx q[3];
rz(-0.10183798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0541075) q[0];
sx q[0];
rz(-2.5547145) q[0];
sx q[0];
rz(0.45463872) q[0];
rz(-3.1184323) q[1];
sx q[1];
rz(-1.4440447) q[1];
sx q[1];
rz(2.815411) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9232193) q[0];
sx q[0];
rz(-2.120864) q[0];
sx q[0];
rz(-1.108842) q[0];
rz(-1.9942547) q[2];
sx q[2];
rz(-1.8427765) q[2];
sx q[2];
rz(2.6620797) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4224976) q[1];
sx q[1];
rz(-1.0651673) q[1];
sx q[1];
rz(1.0984332) q[1];
rz(2.2315027) q[3];
sx q[3];
rz(-1.6636563) q[3];
sx q[3];
rz(1.6805436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0252016) q[2];
sx q[2];
rz(-1.2025183) q[2];
sx q[2];
rz(-0.95428673) q[2];
rz(-1.8918234) q[3];
sx q[3];
rz(-0.3229177) q[3];
sx q[3];
rz(3.1402816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2713476) q[0];
sx q[0];
rz(-1.988669) q[0];
sx q[0];
rz(-1.9648319) q[0];
rz(2.7765043) q[1];
sx q[1];
rz(-1.3026214) q[1];
sx q[1];
rz(1.6129859) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1592401) q[0];
sx q[0];
rz(-1.591658) q[0];
sx q[0];
rz(-3.1283035) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.079054376) q[2];
sx q[2];
rz(-1.1284661) q[2];
sx q[2];
rz(2.4238264) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.075767081) q[1];
sx q[1];
rz(-1.5403374) q[1];
sx q[1];
rz(1.7225313) q[1];
rz(-0.068644329) q[3];
sx q[3];
rz(-1.6283096) q[3];
sx q[3];
rz(-0.36348331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.24885808) q[2];
sx q[2];
rz(-2.5881519) q[2];
sx q[2];
rz(-1.152732) q[2];
rz(-1.24498) q[3];
sx q[3];
rz(-0.98074073) q[3];
sx q[3];
rz(0.28321701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1588441) q[0];
sx q[0];
rz(-11*pi/12) q[0];
sx q[0];
rz(0.68956476) q[0];
rz(0.025029643) q[1];
sx q[1];
rz(-1.6233147) q[1];
sx q[1];
rz(1.4208581) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.384149) q[0];
sx q[0];
rz(-1.4673646) q[0];
sx q[0];
rz(-0.62532702) q[0];
x q[1];
rz(-1.9474073) q[2];
sx q[2];
rz(-0.45832299) q[2];
sx q[2];
rz(-1.6018181) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5031227) q[1];
sx q[1];
rz(-1.5431738) q[1];
sx q[1];
rz(0.50325985) q[1];
rz(-0.33735621) q[3];
sx q[3];
rz(-1.3135664) q[3];
sx q[3];
rz(-1.6599865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.20874061) q[2];
sx q[2];
rz(-0.17540652) q[2];
sx q[2];
rz(-1.1669195) q[2];
rz(2.6973727) q[3];
sx q[3];
rz(-1.2362213) q[3];
sx q[3];
rz(0.3096295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.2883478) q[0];
sx q[0];
rz(-1.7720368) q[0];
sx q[0];
rz(-0.41611588) q[0];
rz(2.6314645) q[1];
sx q[1];
rz(-0.77113873) q[1];
sx q[1];
rz(2.3663734) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32019553) q[0];
sx q[0];
rz(-0.96548432) q[0];
sx q[0];
rz(1.3277675) q[0];
rz(-pi) q[1];
x q[1];
rz(0.40480308) q[2];
sx q[2];
rz(-0.43663132) q[2];
sx q[2];
rz(-0.2704119) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3415592) q[1];
sx q[1];
rz(-1.6960891) q[1];
sx q[1];
rz(3.028246) q[1];
x q[2];
rz(-2.6121077) q[3];
sx q[3];
rz(-0.66505243) q[3];
sx q[3];
rz(0.82529008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.2148718) q[2];
sx q[2];
rz(-0.99433172) q[2];
sx q[2];
rz(1.035824) q[2];
rz(-2.9521247) q[3];
sx q[3];
rz(-2.6697956) q[3];
sx q[3];
rz(2.6389879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.254461) q[0];
sx q[0];
rz(-0.04763617) q[0];
sx q[0];
rz(-0.3558085) q[0];
rz(1.6461146) q[1];
sx q[1];
rz(-2.1864086) q[1];
sx q[1];
rz(-1.1411508) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0063067181) q[0];
sx q[0];
rz(-1.7809488) q[0];
sx q[0];
rz(0.89583136) q[0];
rz(-pi) q[1];
rz(3.1325794) q[2];
sx q[2];
rz(-0.95925943) q[2];
sx q[2];
rz(-2.1076941) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.43737632) q[1];
sx q[1];
rz(-2.9659038) q[1];
sx q[1];
rz(-0.64196569) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0047046) q[3];
sx q[3];
rz(-1.5597789) q[3];
sx q[3];
rz(2.9316255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.63153875) q[2];
sx q[2];
rz(-0.71439356) q[2];
sx q[2];
rz(1.4368524) q[2];
rz(-1.0531744) q[3];
sx q[3];
rz(-1.5318233) q[3];
sx q[3];
rz(-0.50301445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.961504) q[0];
sx q[0];
rz(-2.8453974) q[0];
sx q[0];
rz(-0.12251138) q[0];
rz(0.65613121) q[1];
sx q[1];
rz(-1.2686814) q[1];
sx q[1];
rz(1.8001385) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65402664) q[0];
sx q[0];
rz(-1.6161421) q[0];
sx q[0];
rz(1.5349755) q[0];
rz(-pi) q[1];
rz(-1.6788922) q[2];
sx q[2];
rz(-1.8524395) q[2];
sx q[2];
rz(-0.070527129) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4402561) q[1];
sx q[1];
rz(-2.1783268) q[1];
sx q[1];
rz(-2.8154147) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4571601) q[3];
sx q[3];
rz(-1.1294147) q[3];
sx q[3];
rz(-1.4809004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.573367) q[2];
sx q[2];
rz(-1.919701) q[2];
sx q[2];
rz(1.4637671) q[2];
rz(-0.73417869) q[3];
sx q[3];
rz(-0.28407431) q[3];
sx q[3];
rz(1.3045788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-1.7952591) q[0];
sx q[0];
rz(-1.7743552) q[0];
sx q[0];
rz(-2.6829868) q[0];
rz(2.1268225) q[1];
sx q[1];
rz(-1.9472803) q[1];
sx q[1];
rz(1.0626622) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6342446) q[0];
sx q[0];
rz(-1.5131823) q[0];
sx q[0];
rz(-2.9751124) q[0];
rz(-pi) q[1];
rz(0.80842473) q[2];
sx q[2];
rz(-1.9701517) q[2];
sx q[2];
rz(-2.6724919) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.52673756) q[1];
sx q[1];
rz(-1.4397845) q[1];
sx q[1];
rz(-3.0796771) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8187232) q[3];
sx q[3];
rz(-1.0738651) q[3];
sx q[3];
rz(-0.54847713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.25166544) q[2];
sx q[2];
rz(-1.7532316) q[2];
sx q[2];
rz(-1.7903222) q[2];
rz(-1.9225559) q[3];
sx q[3];
rz(-2.7128897) q[3];
sx q[3];
rz(2.3082699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0625921) q[0];
sx q[0];
rz(-1.3840249) q[0];
sx q[0];
rz(0.90743995) q[0];
rz(0.22706789) q[1];
sx q[1];
rz(-1.0457958) q[1];
sx q[1];
rz(-0.89920941) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.718841) q[0];
sx q[0];
rz(-2.1678574) q[0];
sx q[0];
rz(0.62946749) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4883409) q[2];
sx q[2];
rz(-1.6489779) q[2];
sx q[2];
rz(-0.19985914) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.21630741) q[1];
sx q[1];
rz(-1.7561963) q[1];
sx q[1];
rz(-1.0550631) q[1];
x q[2];
rz(-0.82269086) q[3];
sx q[3];
rz(-1.5334276) q[3];
sx q[3];
rz(2.5189375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9565309) q[2];
sx q[2];
rz(-1.1659634) q[2];
sx q[2];
rz(-2.5409307) q[2];
rz(2.0447842) q[3];
sx q[3];
rz(-0.59022248) q[3];
sx q[3];
rz(-2.2685331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3908865) q[0];
sx q[0];
rz(-2.0426671) q[0];
sx q[0];
rz(-0.98439687) q[0];
rz(-2.6495433) q[1];
sx q[1];
rz(-1.7285408) q[1];
sx q[1];
rz(-1.9459928) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9162684) q[0];
sx q[0];
rz(-1.9231503) q[0];
sx q[0];
rz(-2.5114162) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.49585988) q[2];
sx q[2];
rz(-0.95520077) q[2];
sx q[2];
rz(2.1834508) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3477563) q[1];
sx q[1];
rz(-1.4897921) q[1];
sx q[1];
rz(1.1872613) q[1];
rz(-pi) q[2];
rz(-0.46903761) q[3];
sx q[3];
rz(-0.82995633) q[3];
sx q[3];
rz(1.1560464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.055858) q[2];
sx q[2];
rz(-1.5208289) q[2];
sx q[2];
rz(0.98304191) q[2];
rz(-0.25162697) q[3];
sx q[3];
rz(-2.7139137) q[3];
sx q[3];
rz(-1.0035286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6428103) q[0];
sx q[0];
rz(-0.97220535) q[0];
sx q[0];
rz(-2.4844949) q[0];
rz(1.8819173) q[1];
sx q[1];
rz(-1.7301662) q[1];
sx q[1];
rz(0.87283254) q[1];
rz(0.21435195) q[2];
sx q[2];
rz(-0.8690693) q[2];
sx q[2];
rz(-2.2590841) q[2];
rz(0.76446492) q[3];
sx q[3];
rz(-0.44776147) q[3];
sx q[3];
rz(0.54816435) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
