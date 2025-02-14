OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.79787624) q[0];
sx q[0];
rz(2.1581082) q[0];
sx q[0];
rz(9.6165514) q[0];
rz(-2.6311488) q[1];
sx q[1];
rz(-1.3671083) q[1];
sx q[1];
rz(-1.7459315) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5001427) q[0];
sx q[0];
rz(-3.1063188) q[0];
sx q[0];
rz(-0.39536898) q[0];
rz(-0.22966603) q[2];
sx q[2];
rz(-1.3896835) q[2];
sx q[2];
rz(2.119198) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6916445) q[1];
sx q[1];
rz(-2.3385923) q[1];
sx q[1];
rz(1.2951485) q[1];
rz(-1.1825425) q[3];
sx q[3];
rz(-0.40273977) q[3];
sx q[3];
rz(2.5975063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6338966) q[2];
sx q[2];
rz(-0.32749978) q[2];
sx q[2];
rz(-0.92602777) q[2];
rz(1.6294468) q[3];
sx q[3];
rz(-0.67073268) q[3];
sx q[3];
rz(1.1431471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40874261) q[0];
sx q[0];
rz(-0.73246211) q[0];
sx q[0];
rz(0.23250411) q[0];
rz(-1.1393503) q[1];
sx q[1];
rz(-1.573223) q[1];
sx q[1];
rz(-2.2154636) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8208198) q[0];
sx q[0];
rz(-2.8348587) q[0];
sx q[0];
rz(-0.83267468) q[0];
rz(-pi) q[1];
rz(-1.4732292) q[2];
sx q[2];
rz(-1.2446523) q[2];
sx q[2];
rz(-0.17352428) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.88484832) q[1];
sx q[1];
rz(-1.3646646) q[1];
sx q[1];
rz(0.54289674) q[1];
rz(1.1180463) q[3];
sx q[3];
rz(-1.0501852) q[3];
sx q[3];
rz(-1.3829766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6050379) q[2];
sx q[2];
rz(-1.6190517) q[2];
sx q[2];
rz(-0.79263318) q[2];
rz(-2.7301181) q[3];
sx q[3];
rz(-0.94235197) q[3];
sx q[3];
rz(-0.80504942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3181535) q[0];
sx q[0];
rz(-1.0018438) q[0];
sx q[0];
rz(0.26082984) q[0];
rz(-2.8584495) q[1];
sx q[1];
rz(-1.6915551) q[1];
sx q[1];
rz(-2.0859437) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84427858) q[0];
sx q[0];
rz(-1.7417272) q[0];
sx q[0];
rz(2.2173082) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.63313578) q[2];
sx q[2];
rz(-1.5738259) q[2];
sx q[2];
rz(2.5315447) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4312517) q[1];
sx q[1];
rz(-1.6762937) q[1];
sx q[1];
rz(0.72781095) q[1];
rz(1.0751455) q[3];
sx q[3];
rz(-1.5167674) q[3];
sx q[3];
rz(0.13755218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5259214) q[2];
sx q[2];
rz(-2.0522223) q[2];
sx q[2];
rz(-1.2582568) q[2];
rz(1.0387756) q[3];
sx q[3];
rz(-0.98826161) q[3];
sx q[3];
rz(-1.4163777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0106169) q[0];
sx q[0];
rz(-1.6058141) q[0];
sx q[0];
rz(0.11949874) q[0];
rz(2.9833228) q[1];
sx q[1];
rz(-0.4522849) q[1];
sx q[1];
rz(-1.7012885) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8058469) q[0];
sx q[0];
rz(-2.3290538) q[0];
sx q[0];
rz(-2.0772019) q[0];
x q[1];
rz(-1.2358627) q[2];
sx q[2];
rz(-1.4960297) q[2];
sx q[2];
rz(-0.78943816) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5938737) q[1];
sx q[1];
rz(-0.55127326) q[1];
sx q[1];
rz(-0.51993315) q[1];
rz(-0.97109183) q[3];
sx q[3];
rz(-2.5673742) q[3];
sx q[3];
rz(0.80050877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4629472) q[2];
sx q[2];
rz(-3.0061649) q[2];
sx q[2];
rz(-0.43369183) q[2];
rz(0.67484754) q[3];
sx q[3];
rz(-2.0516472) q[3];
sx q[3];
rz(-1.1564144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7334412) q[0];
sx q[0];
rz(-1.0557405) q[0];
sx q[0];
rz(0.59930402) q[0];
rz(-0.76404461) q[1];
sx q[1];
rz(-1.0300449) q[1];
sx q[1];
rz(-2.5291671) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0611506) q[0];
sx q[0];
rz(-2.168405) q[0];
sx q[0];
rz(-0.1784647) q[0];
rz(-pi) q[1];
rz(3.1407498) q[2];
sx q[2];
rz(-0.71686059) q[2];
sx q[2];
rz(1.5886024) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1445771) q[1];
sx q[1];
rz(-1.738228) q[1];
sx q[1];
rz(2.5232878) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3633943) q[3];
sx q[3];
rz(-0.68420568) q[3];
sx q[3];
rz(-2.9915031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4446438) q[2];
sx q[2];
rz(-0.79978839) q[2];
sx q[2];
rz(-0.4701699) q[2];
rz(-2.7759975) q[3];
sx q[3];
rz(-1.9496893) q[3];
sx q[3];
rz(-2.5308334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.325901) q[0];
sx q[0];
rz(-2.3453562) q[0];
sx q[0];
rz(2.4796487) q[0];
rz(3.0603307) q[1];
sx q[1];
rz(-2.6632402) q[1];
sx q[1];
rz(-0.14019664) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88410035) q[0];
sx q[0];
rz(-2.0292768) q[0];
sx q[0];
rz(-1.8152587) q[0];
x q[1];
rz(0.95023167) q[2];
sx q[2];
rz(-1.3060313) q[2];
sx q[2];
rz(1.0222514) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.473147) q[1];
sx q[1];
rz(-0.47159099) q[1];
sx q[1];
rz(-1.5525329) q[1];
rz(1.7532888) q[3];
sx q[3];
rz(-1.2267707) q[3];
sx q[3];
rz(-1.9545912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.36416546) q[2];
sx q[2];
rz(-1.3157996) q[2];
sx q[2];
rz(2.8167456) q[2];
rz(2.9835564) q[3];
sx q[3];
rz(-2.7285748) q[3];
sx q[3];
rz(1.0154999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3196816) q[0];
sx q[0];
rz(-0.076229036) q[0];
sx q[0];
rz(2.2538189) q[0];
rz(-0.38331389) q[1];
sx q[1];
rz(-2.2593468) q[1];
sx q[1];
rz(-0.15957889) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0903069) q[0];
sx q[0];
rz(-1.2538099) q[0];
sx q[0];
rz(-0.27329926) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.331686) q[2];
sx q[2];
rz(-2.4121373) q[2];
sx q[2];
rz(0.73939656) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.23160203) q[1];
sx q[1];
rz(-1.9456777) q[1];
sx q[1];
rz(-1.7085307) q[1];
x q[2];
rz(0.51377138) q[3];
sx q[3];
rz(-0.32960609) q[3];
sx q[3];
rz(0.20121516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5923656) q[2];
sx q[2];
rz(-1.1575907) q[2];
sx q[2];
rz(1.4166191) q[2];
rz(1.8727411) q[3];
sx q[3];
rz(-2.3609991) q[3];
sx q[3];
rz(-2.1835073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45605993) q[0];
sx q[0];
rz(-1.1153509) q[0];
sx q[0];
rz(-2.2400895) q[0];
rz(2.447336) q[1];
sx q[1];
rz(-1.6777104) q[1];
sx q[1];
rz(1.136397) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6601453) q[0];
sx q[0];
rz(-1.5390696) q[0];
sx q[0];
rz(1.7313787) q[0];
rz(-0.55113422) q[2];
sx q[2];
rz(-1.2493842) q[2];
sx q[2];
rz(-1.2642191) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9688837) q[1];
sx q[1];
rz(-2.4447021) q[1];
sx q[1];
rz(1.4994166) q[1];
x q[2];
rz(-2.4122167) q[3];
sx q[3];
rz(-1.9751901) q[3];
sx q[3];
rz(2.803363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0006492) q[2];
sx q[2];
rz(-1.8724172) q[2];
sx q[2];
rz(2.3285274) q[2];
rz(0.55417577) q[3];
sx q[3];
rz(-1.412609) q[3];
sx q[3];
rz(-0.36177844) q[3];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56162214) q[0];
sx q[0];
rz(-2.0273835) q[0];
sx q[0];
rz(-0.17350523) q[0];
rz(-2.7165727) q[1];
sx q[1];
rz(-1.5360906) q[1];
sx q[1];
rz(-0.82957155) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7873142) q[0];
sx q[0];
rz(-1.4256546) q[0];
sx q[0];
rz(2.4046242) q[0];
x q[1];
rz(2.4194952) q[2];
sx q[2];
rz(-0.061826454) q[2];
sx q[2];
rz(-2.823148) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0922488) q[1];
sx q[1];
rz(-2.4929765) q[1];
sx q[1];
rz(-0.34347024) q[1];
rz(-pi) q[2];
rz(1.5061257) q[3];
sx q[3];
rz(-1.2351278) q[3];
sx q[3];
rz(0.27529374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7912264) q[2];
sx q[2];
rz(-1.841505) q[2];
sx q[2];
rz(1.680797) q[2];
rz(-2.1469877) q[3];
sx q[3];
rz(-0.9959144) q[3];
sx q[3];
rz(-0.88609707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9643672) q[0];
sx q[0];
rz(-2.565964) q[0];
sx q[0];
rz(3.0395569) q[0];
rz(2.521934) q[1];
sx q[1];
rz(-1.4980059) q[1];
sx q[1];
rz(-2.5443351) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12349081) q[0];
sx q[0];
rz(-2.1649556) q[0];
sx q[0];
rz(1.4090562) q[0];
x q[1];
rz(-1.4240392) q[2];
sx q[2];
rz(-1.9224615) q[2];
sx q[2];
rz(0.41837012) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6957132) q[1];
sx q[1];
rz(-1.5339601) q[1];
sx q[1];
rz(-1.0050773) q[1];
rz(-pi) q[2];
x q[2];
rz(0.064205211) q[3];
sx q[3];
rz(-2.6693025) q[3];
sx q[3];
rz(2.9407168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.20798802) q[2];
sx q[2];
rz(-0.36078578) q[2];
sx q[2];
rz(2.2130845) q[2];
rz(-1.9814631) q[3];
sx q[3];
rz(-1.4312294) q[3];
sx q[3];
rz(0.76735705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3930897) q[0];
sx q[0];
rz(-0.7001644) q[0];
sx q[0];
rz(-2.9099332) q[0];
rz(-2.0139991) q[1];
sx q[1];
rz(-2.6078106) q[1];
sx q[1];
rz(-0.003905205) q[1];
rz(-0.071795287) q[2];
sx q[2];
rz(-2.2459002) q[2];
sx q[2];
rz(-0.33159524) q[2];
rz(1.9740022) q[3];
sx q[3];
rz(-0.73763631) q[3];
sx q[3];
rz(0.49745001) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
