OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7528485) q[0];
sx q[0];
rz(-0.53628439) q[0];
sx q[0];
rz(2.1938238) q[0];
rz(1.8127958) q[1];
sx q[1];
rz(-1.2674018) q[1];
sx q[1];
rz(-1.0277494) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2078903) q[0];
sx q[0];
rz(-0.37651248) q[0];
sx q[0];
rz(3.0794789) q[0];
rz(-pi) q[1];
rz(-1.9704291) q[2];
sx q[2];
rz(-0.52414775) q[2];
sx q[2];
rz(-1.580796) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6537522) q[1];
sx q[1];
rz(-2.2334705) q[1];
sx q[1];
rz(-1.0767879) q[1];
rz(0.62698934) q[3];
sx q[3];
rz(-1.9416182) q[3];
sx q[3];
rz(-2.2771319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1296922) q[2];
sx q[2];
rz(-1.7069858) q[2];
sx q[2];
rz(2.091308) q[2];
rz(-2.0283279) q[3];
sx q[3];
rz(-0.89171019) q[3];
sx q[3];
rz(0.068107001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0691836) q[0];
sx q[0];
rz(-1.2453112) q[0];
sx q[0];
rz(0.29775277) q[0];
rz(-0.61966664) q[1];
sx q[1];
rz(-1.0071808) q[1];
sx q[1];
rz(2.0334977) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71582687) q[0];
sx q[0];
rz(-2.1799488) q[0];
sx q[0];
rz(-2.8858375) q[0];
x q[1];
rz(2.1790702) q[2];
sx q[2];
rz(-1.9307185) q[2];
sx q[2];
rz(2.4134709) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3023492) q[1];
sx q[1];
rz(-1.9560555) q[1];
sx q[1];
rz(1.1883931) q[1];
x q[2];
rz(2.439019) q[3];
sx q[3];
rz(-1.3630023) q[3];
sx q[3];
rz(-2.7489565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6790598) q[2];
sx q[2];
rz(-0.97683895) q[2];
sx q[2];
rz(0.97529808) q[2];
rz(-2.2235928) q[3];
sx q[3];
rz(-1.2851597) q[3];
sx q[3];
rz(-2.8454034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96238962) q[0];
sx q[0];
rz(-0.85150349) q[0];
sx q[0];
rz(0.54291022) q[0];
rz(-0.88223282) q[1];
sx q[1];
rz(-2.0062607) q[1];
sx q[1];
rz(-2.1767445) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1436413) q[0];
sx q[0];
rz(-0.43295857) q[0];
sx q[0];
rz(0.71822449) q[0];
x q[1];
rz(2.8446571) q[2];
sx q[2];
rz(-1.6116217) q[2];
sx q[2];
rz(0.21619913) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1988586) q[1];
sx q[1];
rz(-3.004289) q[1];
sx q[1];
rz(1.0820461) q[1];
rz(-pi) q[2];
rz(-1.1539677) q[3];
sx q[3];
rz(-1.7524476) q[3];
sx q[3];
rz(-2.3895398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.09459153) q[2];
sx q[2];
rz(-2.5307405) q[2];
sx q[2];
rz(-1.1331406) q[2];
rz(2.9099693) q[3];
sx q[3];
rz(-1.8685721) q[3];
sx q[3];
rz(2.384322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27947458) q[0];
sx q[0];
rz(-0.010443895) q[0];
sx q[0];
rz(1.3765155) q[0];
rz(0.51849413) q[1];
sx q[1];
rz(-1.2644178) q[1];
sx q[1];
rz(2.8994697) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0771368) q[0];
sx q[0];
rz(-0.69364871) q[0];
sx q[0];
rz(-0.25607381) q[0];
rz(-pi) q[1];
rz(0.42963117) q[2];
sx q[2];
rz(-1.5589082) q[2];
sx q[2];
rz(-2.5846543) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.33430797) q[1];
sx q[1];
rz(-1.1242928) q[1];
sx q[1];
rz(-1.7590894) q[1];
rz(-pi) q[2];
rz(-1.4297156) q[3];
sx q[3];
rz(-1.1564848) q[3];
sx q[3];
rz(-0.29841081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.018192856) q[2];
sx q[2];
rz(-2.172956) q[2];
sx q[2];
rz(-3.0467395) q[2];
rz(1.3421966) q[3];
sx q[3];
rz(-1.7443402) q[3];
sx q[3];
rz(0.43911394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78385329) q[0];
sx q[0];
rz(-1.8308324) q[0];
sx q[0];
rz(2.7807996) q[0];
rz(-1.7533253) q[1];
sx q[1];
rz(-1.8107982) q[1];
sx q[1];
rz(-2.0070019) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9903588) q[0];
sx q[0];
rz(-0.84999527) q[0];
sx q[0];
rz(-2.4609341) q[0];
x q[1];
rz(0.35933944) q[2];
sx q[2];
rz(-1.1655072) q[2];
sx q[2];
rz(2.5422424) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1281631) q[1];
sx q[1];
rz(-1.1268106) q[1];
sx q[1];
rz(1.7432937) q[1];
rz(0.17825019) q[3];
sx q[3];
rz(-2.6986487) q[3];
sx q[3];
rz(2.6435542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1363042) q[2];
sx q[2];
rz(-1.4207999) q[2];
sx q[2];
rz(-0.57265442) q[2];
rz(-0.92875656) q[3];
sx q[3];
rz(-2.6199665) q[3];
sx q[3];
rz(1.1675534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3271493) q[0];
sx q[0];
rz(-1.0961908) q[0];
sx q[0];
rz(-1.8922528) q[0];
rz(-1.2231187) q[1];
sx q[1];
rz(-1.616281) q[1];
sx q[1];
rz(1.1522326) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1359065) q[0];
sx q[0];
rz(-1.583477) q[0];
sx q[0];
rz(1.6760875) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.639159) q[2];
sx q[2];
rz(-1.2049434) q[2];
sx q[2];
rz(2.1937214) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9705379) q[1];
sx q[1];
rz(-1.4661403) q[1];
sx q[1];
rz(-2.9403951) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96529393) q[3];
sx q[3];
rz(-0.98635841) q[3];
sx q[3];
rz(-2.8402929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.46889177) q[2];
sx q[2];
rz(-1.2142618) q[2];
sx q[2];
rz(-2.1506298) q[2];
rz(0.64783603) q[3];
sx q[3];
rz(-0.9655374) q[3];
sx q[3];
rz(-2.794054) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25794849) q[0];
sx q[0];
rz(-0.22709665) q[0];
sx q[0];
rz(0.062285475) q[0];
rz(2.9557872) q[1];
sx q[1];
rz(-1.4567016) q[1];
sx q[1];
rz(2.7468162) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6900401) q[0];
sx q[0];
rz(-1.8071022) q[0];
sx q[0];
rz(2.7355746) q[0];
rz(2.4418418) q[2];
sx q[2];
rz(-2.664898) q[2];
sx q[2];
rz(-1.2303908) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2784087) q[1];
sx q[1];
rz(-1.3235958) q[1];
sx q[1];
rz(2.7017038) q[1];
rz(0.23468252) q[3];
sx q[3];
rz(-2.9838786) q[3];
sx q[3];
rz(1.371067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4930967) q[2];
sx q[2];
rz(-2.811537) q[2];
sx q[2];
rz(-2.8685692) q[2];
rz(1.3027044) q[3];
sx q[3];
rz(-1.3132934) q[3];
sx q[3];
rz(0.31204143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39772314) q[0];
sx q[0];
rz(-1.1345154) q[0];
sx q[0];
rz(-1.8564818) q[0];
rz(-1.5015191) q[1];
sx q[1];
rz(-1.7506426) q[1];
sx q[1];
rz(-1.8008908) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9047324) q[0];
sx q[0];
rz(-1.8956087) q[0];
sx q[0];
rz(-1.1572641) q[0];
rz(-pi) q[1];
rz(1.5199392) q[2];
sx q[2];
rz(-0.87561456) q[2];
sx q[2];
rz(-2.4620172) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2015842) q[1];
sx q[1];
rz(-1.219698) q[1];
sx q[1];
rz(1.9560948) q[1];
rz(-pi) q[2];
x q[2];
rz(2.826346) q[3];
sx q[3];
rz(-1.1479706) q[3];
sx q[3];
rz(1.3413615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0354707) q[2];
sx q[2];
rz(-0.80646986) q[2];
sx q[2];
rz(-1.0236053) q[2];
rz(-0.18493955) q[3];
sx q[3];
rz(-0.39026323) q[3];
sx q[3];
rz(-2.8997054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(-2.0969365) q[0];
sx q[0];
rz(-2.1450295) q[0];
sx q[0];
rz(-1.5203083) q[0];
rz(2.8114491) q[1];
sx q[1];
rz(-1.9338927) q[1];
sx q[1];
rz(2.3044589) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5846111) q[0];
sx q[0];
rz(-0.71174445) q[0];
sx q[0];
rz(2.8296489) q[0];
rz(-pi) q[1];
rz(-0.57640055) q[2];
sx q[2];
rz(-1.3199558) q[2];
sx q[2];
rz(2.5121411) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3771462) q[1];
sx q[1];
rz(-0.18491491) q[1];
sx q[1];
rz(0.96692337) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0650474) q[3];
sx q[3];
rz(-0.55378434) q[3];
sx q[3];
rz(-1.6365285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5403486) q[2];
sx q[2];
rz(-2.0566172) q[2];
sx q[2];
rz(2.9619651) q[2];
rz(0.99572292) q[3];
sx q[3];
rz(-1.8928173) q[3];
sx q[3];
rz(1.3109591) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3778465) q[0];
sx q[0];
rz(-0.34559956) q[0];
sx q[0];
rz(-1.0572222) q[0];
rz(0.10748848) q[1];
sx q[1];
rz(-1.2534393) q[1];
sx q[1];
rz(0.9799788) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1996795) q[0];
sx q[0];
rz(-1.4451888) q[0];
sx q[0];
rz(3.0903387) q[0];
rz(-pi) q[1];
rz(-2.0699632) q[2];
sx q[2];
rz(-1.1071148) q[2];
sx q[2];
rz(1.0774563) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6085538) q[1];
sx q[1];
rz(-1.3320919) q[1];
sx q[1];
rz(1.1281611) q[1];
rz(-pi) q[2];
rz(0.1006871) q[3];
sx q[3];
rz(-1.8929314) q[3];
sx q[3];
rz(-1.611657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8367299) q[2];
sx q[2];
rz(-1.2048081) q[2];
sx q[2];
rz(-2.1137962) q[2];
rz(1.7547539) q[3];
sx q[3];
rz(-1.3091062) q[3];
sx q[3];
rz(0.28361472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2733611) q[0];
sx q[0];
rz(-2.1049451) q[0];
sx q[0];
rz(2.0275397) q[0];
rz(-0.83203075) q[1];
sx q[1];
rz(-2.6770626) q[1];
sx q[1];
rz(-2.4774036) q[1];
rz(-1.2912512) q[2];
sx q[2];
rz(-0.57689473) q[2];
sx q[2];
rz(2.9070791) q[2];
rz(0.63315331) q[3];
sx q[3];
rz(-0.59637759) q[3];
sx q[3];
rz(-0.37806088) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
