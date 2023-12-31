OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.54685932) q[0];
sx q[0];
rz(4.6580553) q[0];
sx q[0];
rz(9.1604995) q[0];
rz(-0.9737941) q[1];
sx q[1];
rz(5.073054) q[1];
sx q[1];
rz(10.160025) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7967448) q[0];
sx q[0];
rz(-1.42294) q[0];
sx q[0];
rz(2.4792838) q[0];
rz(-pi) q[1];
rz(-2.0779607) q[2];
sx q[2];
rz(-0.9557561) q[2];
sx q[2];
rz(-2.087649) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.211261) q[1];
sx q[1];
rz(-2.8144565) q[1];
sx q[1];
rz(-1.0717908) q[1];
x q[2];
rz(3.0987708) q[3];
sx q[3];
rz(-0.58086568) q[3];
sx q[3];
rz(-2.7659741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4720817) q[2];
sx q[2];
rz(-1.8410204) q[2];
sx q[2];
rz(-2.0377339) q[2];
rz(-1.2708698) q[3];
sx q[3];
rz(-1.2277675) q[3];
sx q[3];
rz(0.27403533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0274149) q[0];
sx q[0];
rz(-0.52863055) q[0];
sx q[0];
rz(0.43637481) q[0];
rz(2.6787058) q[1];
sx q[1];
rz(-1.0375689) q[1];
sx q[1];
rz(-2.8754821) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3811831) q[0];
sx q[0];
rz(-1.1927483) q[0];
sx q[0];
rz(-0.80947431) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.16337784) q[2];
sx q[2];
rz(-1.0014357) q[2];
sx q[2];
rz(0.53158224) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4437342) q[1];
sx q[1];
rz(-2.2599972) q[1];
sx q[1];
rz(-0.84390784) q[1];
rz(-pi) q[2];
rz(2.9588685) q[3];
sx q[3];
rz(-2.1696739) q[3];
sx q[3];
rz(-0.11174186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6767072) q[2];
sx q[2];
rz(-1.8647944) q[2];
sx q[2];
rz(-0.51149386) q[2];
rz(-2.3320847) q[3];
sx q[3];
rz(-1.5313238) q[3];
sx q[3];
rz(2.8539343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4354316) q[0];
sx q[0];
rz(-1.5478739) q[0];
sx q[0];
rz(-2.2128552) q[0];
rz(1.7354895) q[1];
sx q[1];
rz(-0.7000674) q[1];
sx q[1];
rz(-1.7944638) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12362326) q[0];
sx q[0];
rz(-0.53640134) q[0];
sx q[0];
rz(2.6111952) q[0];
x q[1];
rz(1.2855661) q[2];
sx q[2];
rz(-1.4189548) q[2];
sx q[2];
rz(-2.1053932) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.82145065) q[1];
sx q[1];
rz(-0.67678932) q[1];
sx q[1];
rz(-2.0435145) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8637572) q[3];
sx q[3];
rz(-1.45544) q[3];
sx q[3];
rz(-2.399721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9946263) q[2];
sx q[2];
rz(-1.7549843) q[2];
sx q[2];
rz(1.6195126) q[2];
rz(-0.26432031) q[3];
sx q[3];
rz(-2.1051354) q[3];
sx q[3];
rz(0.49595293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3494444) q[0];
sx q[0];
rz(-1.148372) q[0];
sx q[0];
rz(0.96570063) q[0];
rz(-2.4194338) q[1];
sx q[1];
rz(-1.5042217) q[1];
sx q[1];
rz(-0.55975634) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4013195) q[0];
sx q[0];
rz(-1.6519321) q[0];
sx q[0];
rz(-2.7149537) q[0];
rz(2.3279834) q[2];
sx q[2];
rz(-1.0126197) q[2];
sx q[2];
rz(-1.8166325) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8430082) q[1];
sx q[1];
rz(-2.1857939) q[1];
sx q[1];
rz(3.0857012) q[1];
x q[2];
rz(0.23493725) q[3];
sx q[3];
rz(-2.3384691) q[3];
sx q[3];
rz(2.0149751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7136148) q[2];
sx q[2];
rz(-1.6327991) q[2];
sx q[2];
rz(-1.4245865) q[2];
rz(2.8811841) q[3];
sx q[3];
rz(-1.3924761) q[3];
sx q[3];
rz(-0.4310472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99073064) q[0];
sx q[0];
rz(-1.4980415) q[0];
sx q[0];
rz(0.36636233) q[0];
rz(-1.5953966) q[1];
sx q[1];
rz(-2.5876744) q[1];
sx q[1];
rz(0.34367925) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77728358) q[0];
sx q[0];
rz(-1.0316348) q[0];
sx q[0];
rz(-0.76215141) q[0];
rz(-pi) q[1];
rz(2.2720488) q[2];
sx q[2];
rz(-1.1079259) q[2];
sx q[2];
rz(-2.9733544) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.74663631) q[1];
sx q[1];
rz(-0.50368217) q[1];
sx q[1];
rz(1.4175182) q[1];
rz(-1.5424535) q[3];
sx q[3];
rz(-2.4975371) q[3];
sx q[3];
rz(-1.5202265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7498103) q[2];
sx q[2];
rz(-2.2822773) q[2];
sx q[2];
rz(3.0878477) q[2];
rz(1.404445) q[3];
sx q[3];
rz(-2.6440547) q[3];
sx q[3];
rz(2.8500309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6546201) q[0];
sx q[0];
rz(-0.64990652) q[0];
sx q[0];
rz(-2.3858331) q[0];
rz(3.1164363) q[1];
sx q[1];
rz(-0.92725602) q[1];
sx q[1];
rz(-0.25973928) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4070968) q[0];
sx q[0];
rz(-1.8844814) q[0];
sx q[0];
rz(1.19019) q[0];
rz(0.01979205) q[2];
sx q[2];
rz(-1.9070101) q[2];
sx q[2];
rz(2.3078231) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.427634) q[1];
sx q[1];
rz(-2.7379588) q[1];
sx q[1];
rz(0.9637109) q[1];
rz(-pi) q[2];
rz(-3.1236157) q[3];
sx q[3];
rz(-2.020105) q[3];
sx q[3];
rz(0.22389212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.48866895) q[2];
sx q[2];
rz(-2.3355464) q[2];
sx q[2];
rz(3.0409813) q[2];
rz(0.18209022) q[3];
sx q[3];
rz(-0.87825769) q[3];
sx q[3];
rz(1.7939059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9942193) q[0];
sx q[0];
rz(-1.4202776) q[0];
sx q[0];
rz(0.31016645) q[0];
rz(-0.50225964) q[1];
sx q[1];
rz(-2.6175833) q[1];
sx q[1];
rz(2.535634) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80752559) q[0];
sx q[0];
rz(-1.345732) q[0];
sx q[0];
rz(-2.5146757) q[0];
rz(-2.8499243) q[2];
sx q[2];
rz(-1.3783611) q[2];
sx q[2];
rz(0.33484909) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.76304945) q[1];
sx q[1];
rz(-2.0705283) q[1];
sx q[1];
rz(-0.39079697) q[1];
x q[2];
rz(-2.1171655) q[3];
sx q[3];
rz(-1.1064648) q[3];
sx q[3];
rz(-1.6950316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.24017748) q[2];
sx q[2];
rz(-1.1837974) q[2];
sx q[2];
rz(0.85285464) q[2];
rz(-1.7715706) q[3];
sx q[3];
rz(-1.4586689) q[3];
sx q[3];
rz(-0.20496932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2296427) q[0];
sx q[0];
rz(-0.59597534) q[0];
sx q[0];
rz(1.6802616) q[0];
rz(-1.7386859) q[1];
sx q[1];
rz(-2.1673514) q[1];
sx q[1];
rz(0.064037474) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2647576) q[0];
sx q[0];
rz(-1.3834073) q[0];
sx q[0];
rz(0.95215709) q[0];
x q[1];
rz(-2.089558) q[2];
sx q[2];
rz(-2.5737737) q[2];
sx q[2];
rz(-1.0713112) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.774051) q[1];
sx q[1];
rz(-1.6072175) q[1];
sx q[1];
rz(-0.35518412) q[1];
rz(-1.2588345) q[3];
sx q[3];
rz(-1.6083816) q[3];
sx q[3];
rz(-0.037548693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0503851) q[2];
sx q[2];
rz(-2.51077) q[2];
sx q[2];
rz(-1.5861661) q[2];
rz(0.88820109) q[3];
sx q[3];
rz(-1.9017838) q[3];
sx q[3];
rz(-0.92938882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14426194) q[0];
sx q[0];
rz(-2.0781131) q[0];
sx q[0];
rz(1.7096747) q[0];
rz(2.572708) q[1];
sx q[1];
rz(-0.535393) q[1];
sx q[1];
rz(1.127839) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.226798) q[0];
sx q[0];
rz(-2.9920122) q[0];
sx q[0];
rz(-3.0339255) q[0];
x q[1];
rz(-2.608689) q[2];
sx q[2];
rz(-1.0120631) q[2];
sx q[2];
rz(-1.2677873) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4081501) q[1];
sx q[1];
rz(-0.24194716) q[1];
sx q[1];
rz(-1.5223632) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1209675) q[3];
sx q[3];
rz(-0.63311011) q[3];
sx q[3];
rz(1.8819295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.52788064) q[2];
sx q[2];
rz(-2.2884559) q[2];
sx q[2];
rz(0.17364994) q[2];
rz(-2.8052143) q[3];
sx q[3];
rz(-1.9189546) q[3];
sx q[3];
rz(-2.7500847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7062475) q[0];
sx q[0];
rz(-0.58723891) q[0];
sx q[0];
rz(-1.6760814) q[0];
rz(0.82410518) q[1];
sx q[1];
rz(-1.5287639) q[1];
sx q[1];
rz(-0.5724268) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49003285) q[0];
sx q[0];
rz(-0.76793725) q[0];
sx q[0];
rz(-2.4134273) q[0];
x q[1];
rz(-2.0881537) q[2];
sx q[2];
rz(-0.51689076) q[2];
sx q[2];
rz(-3.0992532) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4932369) q[1];
sx q[1];
rz(-1.3322468) q[1];
sx q[1];
rz(1.7931213) q[1];
rz(-2.8045373) q[3];
sx q[3];
rz(-0.73738499) q[3];
sx q[3];
rz(-2.9205703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.77999014) q[2];
sx q[2];
rz(-2.6670691) q[2];
sx q[2];
rz(-0.53722107) q[2];
rz(-1.0572664) q[3];
sx q[3];
rz(-0.89151645) q[3];
sx q[3];
rz(-2.4479772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2873516) q[0];
sx q[0];
rz(-1.1653405) q[0];
sx q[0];
rz(-1.5821138) q[0];
rz(-0.46335012) q[1];
sx q[1];
rz(-0.87711038) q[1];
sx q[1];
rz(-1.6323485) q[1];
rz(2.2255185) q[2];
sx q[2];
rz(-1.4929885) q[2];
sx q[2];
rz(2.3103466) q[2];
rz(-1.6871917) q[3];
sx q[3];
rz(-2.6506861) q[3];
sx q[3];
rz(-0.41674137) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
