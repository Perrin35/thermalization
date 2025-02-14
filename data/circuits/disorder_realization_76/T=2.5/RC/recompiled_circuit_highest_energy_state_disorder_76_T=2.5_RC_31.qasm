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
rz(0.87585706) q[0];
sx q[0];
rz(-0.96867222) q[0];
sx q[0];
rz(2.2157366) q[0];
rz(0.61697382) q[1];
sx q[1];
rz(-0.66520912) q[1];
sx q[1];
rz(-1.8170504) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61875859) q[0];
sx q[0];
rz(-1.6346187) q[0];
sx q[0];
rz(-2.5413949) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.994903) q[2];
sx q[2];
rz(-2.5518199) q[2];
sx q[2];
rz(-1.1663933) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1086667) q[1];
sx q[1];
rz(-0.98331988) q[1];
sx q[1];
rz(2.1290178) q[1];
rz(-pi) q[2];
rz(2.2258513) q[3];
sx q[3];
rz(-0.91980442) q[3];
sx q[3];
rz(-2.8398193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2609743) q[2];
sx q[2];
rz(-1.7558492) q[2];
sx q[2];
rz(2.3495038) q[2];
rz(0.875862) q[3];
sx q[3];
rz(-3.0328817) q[3];
sx q[3];
rz(-2.47827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(2.6302014) q[0];
sx q[0];
rz(-1.8237317) q[0];
sx q[0];
rz(-0.9414916) q[0];
rz(-0.9990274) q[1];
sx q[1];
rz(-2.2143054) q[1];
sx q[1];
rz(-0.082854465) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4946144) q[0];
sx q[0];
rz(-2.4002693) q[0];
sx q[0];
rz(-2.9927918) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.98495667) q[2];
sx q[2];
rz(-0.42030605) q[2];
sx q[2];
rz(-0.46520761) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.99020834) q[1];
sx q[1];
rz(-2.4308204) q[1];
sx q[1];
rz(-2.4113893) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1762455) q[3];
sx q[3];
rz(-2.1248528) q[3];
sx q[3];
rz(-0.29747552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.89722172) q[2];
sx q[2];
rz(-1.3542078) q[2];
sx q[2];
rz(-1.2574035) q[2];
rz(2.2952378) q[3];
sx q[3];
rz(-2.9823163) q[3];
sx q[3];
rz(2.4634821) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9427247) q[0];
sx q[0];
rz(-1.6833479) q[0];
sx q[0];
rz(-0.22920907) q[0];
rz(2.9275059) q[1];
sx q[1];
rz(-2.2702718) q[1];
sx q[1];
rz(-0.3814989) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4320735) q[0];
sx q[0];
rz(-2.1901178) q[0];
sx q[0];
rz(-1.1362057) q[0];
rz(1.7842641) q[2];
sx q[2];
rz(-2.3980283) q[2];
sx q[2];
rz(-1.5876169) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9126622) q[1];
sx q[1];
rz(-0.060256392) q[1];
sx q[1];
rz(2.069239) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5759141) q[3];
sx q[3];
rz(-1.2237142) q[3];
sx q[3];
rz(-2.0818613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4182338) q[2];
sx q[2];
rz(-2.8853719) q[2];
sx q[2];
rz(2.5733433) q[2];
rz(1.7226284) q[3];
sx q[3];
rz(-1.4293554) q[3];
sx q[3];
rz(2.0645781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-3.1133465) q[0];
sx q[0];
rz(-1.1320817) q[0];
sx q[0];
rz(-0.25076732) q[0];
rz(3.057042) q[1];
sx q[1];
rz(-1.0944347) q[1];
sx q[1];
rz(1.9409404) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47287175) q[0];
sx q[0];
rz(-0.78470147) q[0];
sx q[0];
rz(2.0913823) q[0];
x q[1];
rz(1.1293639) q[2];
sx q[2];
rz(-1.7953292) q[2];
sx q[2];
rz(2.6972636) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1657287) q[1];
sx q[1];
rz(-1.670125) q[1];
sx q[1];
rz(-0.71290675) q[1];
x q[2];
rz(-2.8765466) q[3];
sx q[3];
rz(-1.0599905) q[3];
sx q[3];
rz(-3.0731415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.67479977) q[2];
sx q[2];
rz(-1.1226706) q[2];
sx q[2];
rz(1.2910845) q[2];
rz(-2.71991) q[3];
sx q[3];
rz(-2.6899874) q[3];
sx q[3];
rz(1.5765367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6186433) q[0];
sx q[0];
rz(-2.4186501) q[0];
sx q[0];
rz(-1.6408386) q[0];
rz(-0.067642637) q[1];
sx q[1];
rz(-1.5875971) q[1];
sx q[1];
rz(-1.1281475) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4823032) q[0];
sx q[0];
rz(-1.3829675) q[0];
sx q[0];
rz(-1.6085621) q[0];
rz(2.7893547) q[2];
sx q[2];
rz(-1.0718498) q[2];
sx q[2];
rz(-1.3832472) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.50666821) q[1];
sx q[1];
rz(-1.4574058) q[1];
sx q[1];
rz(3.0874814) q[1];
rz(0.87781436) q[3];
sx q[3];
rz(-2.7846309) q[3];
sx q[3];
rz(-1.9138543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.03380123) q[2];
sx q[2];
rz(-1.1148323) q[2];
sx q[2];
rz(-2.4647253) q[2];
rz(0.90773165) q[3];
sx q[3];
rz(-1.2264484) q[3];
sx q[3];
rz(2.7270253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21861976) q[0];
sx q[0];
rz(-0.75160471) q[0];
sx q[0];
rz(2.3722017) q[0];
rz(-1.9006624) q[1];
sx q[1];
rz(-2.2878094) q[1];
sx q[1];
rz(2.9262537) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8148635) q[0];
sx q[0];
rz(-1.0007326) q[0];
sx q[0];
rz(0.22320052) q[0];
rz(1.5211283) q[2];
sx q[2];
rz(-2.81041) q[2];
sx q[2];
rz(2.0803723) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8465855) q[1];
sx q[1];
rz(-2.2071116) q[1];
sx q[1];
rz(-1.8759439) q[1];
rz(-pi) q[2];
rz(1.9086544) q[3];
sx q[3];
rz(-0.90150276) q[3];
sx q[3];
rz(0.026642628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.26663366) q[2];
sx q[2];
rz(-1.5077488) q[2];
sx q[2];
rz(2.5249262) q[2];
rz(2.941361) q[3];
sx q[3];
rz(-0.73035208) q[3];
sx q[3];
rz(2.0088137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1277593) q[0];
sx q[0];
rz(-0.56202373) q[0];
sx q[0];
rz(-0.01509893) q[0];
rz(3.0191782) q[1];
sx q[1];
rz(-1.2950803) q[1];
sx q[1];
rz(-2.1133568) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0434179) q[0];
sx q[0];
rz(-2.2255996) q[0];
sx q[0];
rz(-1.6408987) q[0];
rz(-2.188855) q[2];
sx q[2];
rz(-1.3985123) q[2];
sx q[2];
rz(3.1069482) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.15785698) q[1];
sx q[1];
rz(-0.4326371) q[1];
sx q[1];
rz(3.0745201) q[1];
x q[2];
rz(-0.7042775) q[3];
sx q[3];
rz(-2.7368494) q[3];
sx q[3];
rz(1.2621439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5369109) q[2];
sx q[2];
rz(-1.9313507) q[2];
sx q[2];
rz(0.49986419) q[2];
rz(2.8258421) q[3];
sx q[3];
rz(-2.4707268) q[3];
sx q[3];
rz(2.1226814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4278118) q[0];
sx q[0];
rz(-1.918387) q[0];
sx q[0];
rz(2.1977303) q[0];
rz(-1.4048514) q[1];
sx q[1];
rz(-1.8753139) q[1];
sx q[1];
rz(-0.92165438) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.114678) q[0];
sx q[0];
rz(-1.252838) q[0];
sx q[0];
rz(1.6284031) q[0];
rz(-pi) q[1];
rz(-1.2723075) q[2];
sx q[2];
rz(-1.3278074) q[2];
sx q[2];
rz(1.9858433) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7870175) q[1];
sx q[1];
rz(-2.0817309) q[1];
sx q[1];
rz(2.4035011) q[1];
rz(-pi) q[2];
rz(1.5632024) q[3];
sx q[3];
rz(-2.6287327) q[3];
sx q[3];
rz(-0.22554413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1450119) q[2];
sx q[2];
rz(-1.8700446) q[2];
sx q[2];
rz(-1.5550295) q[2];
rz(-2.0874713) q[3];
sx q[3];
rz(-1.328822) q[3];
sx q[3];
rz(-1.740295) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9911875) q[0];
sx q[0];
rz(-1.0973955) q[0];
sx q[0];
rz(-2.4990668) q[0];
rz(1.7036899) q[1];
sx q[1];
rz(-0.75690401) q[1];
sx q[1];
rz(0.14370758) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8929307) q[0];
sx q[0];
rz(-2.4911059) q[0];
sx q[0];
rz(-2.6269795) q[0];
rz(-pi) q[1];
x q[1];
rz(0.52395487) q[2];
sx q[2];
rz(-0.99159504) q[2];
sx q[2];
rz(-0.29925811) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4776996) q[1];
sx q[1];
rz(-1.8346922) q[1];
sx q[1];
rz(0.66910918) q[1];
x q[2];
rz(-2.4293184) q[3];
sx q[3];
rz(-2.5872719) q[3];
sx q[3];
rz(1.2620776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5370499) q[2];
sx q[2];
rz(-1.9436676) q[2];
sx q[2];
rz(-0.26270467) q[2];
rz(-2.883834) q[3];
sx q[3];
rz(-0.38193211) q[3];
sx q[3];
rz(-0.8530544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2951374) q[0];
sx q[0];
rz(-1.4895804) q[0];
sx q[0];
rz(-1.9211796) q[0];
rz(2.659761) q[1];
sx q[1];
rz(-0.82492963) q[1];
sx q[1];
rz(0.89734546) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8751971) q[0];
sx q[0];
rz(-1.0561854) q[0];
sx q[0];
rz(1.4304016) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4608896) q[2];
sx q[2];
rz(-1.8898003) q[2];
sx q[2];
rz(-1.4665335) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.022696115) q[1];
sx q[1];
rz(-1.5210946) q[1];
sx q[1];
rz(-1.3124052) q[1];
x q[2];
rz(-1.4225716) q[3];
sx q[3];
rz(-0.61398849) q[3];
sx q[3];
rz(2.2999291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.64409488) q[2];
sx q[2];
rz(-0.92659014) q[2];
sx q[2];
rz(0.11478718) q[2];
rz(0.74350205) q[3];
sx q[3];
rz(-1.8758352) q[3];
sx q[3];
rz(-0.44873294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7262309) q[0];
sx q[0];
rz(-2.3804433) q[0];
sx q[0];
rz(-0.95809715) q[0];
rz(1.505898) q[1];
sx q[1];
rz(-2.0151357) q[1];
sx q[1];
rz(-1.9912079) q[1];
rz(-2.220357) q[2];
sx q[2];
rz(-0.9365281) q[2];
sx q[2];
rz(-2.3021883) q[2];
rz(-0.34974364) q[3];
sx q[3];
rz(-0.45827023) q[3];
sx q[3];
rz(-2.6029725) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
