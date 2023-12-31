OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.5503791) q[0];
sx q[0];
rz(3.1449218) q[0];
sx q[0];
rz(9.7700906) q[0];
rz(-1.3357063) q[1];
sx q[1];
rz(-0.3392646) q[1];
sx q[1];
rz(-2.8621434) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0246668) q[0];
sx q[0];
rz(-1.7755839) q[0];
sx q[0];
rz(1.8860399) q[0];
x q[1];
rz(-0.34502132) q[2];
sx q[2];
rz(-2.3292688) q[2];
sx q[2];
rz(1.617384) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.93912032) q[1];
sx q[1];
rz(-1.3807553) q[1];
sx q[1];
rz(0.01357667) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.64198288) q[3];
sx q[3];
rz(-1.5390453) q[3];
sx q[3];
rz(0.25105219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.44148579) q[2];
sx q[2];
rz(-1.6690648) q[2];
sx q[2];
rz(2.1842365) q[2];
rz(0.29933128) q[3];
sx q[3];
rz(-2.744031) q[3];
sx q[3];
rz(2.7295952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9807724) q[0];
sx q[0];
rz(-0.94962025) q[0];
sx q[0];
rz(-2.7278996) q[0];
rz(-1.7970239) q[1];
sx q[1];
rz(-0.78318703) q[1];
sx q[1];
rz(2.5059674) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6247511) q[0];
sx q[0];
rz(-3.0020614) q[0];
sx q[0];
rz(-1.2463039) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0947729) q[2];
sx q[2];
rz(-1.4390107) q[2];
sx q[2];
rz(1.5891967) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7285068) q[1];
sx q[1];
rz(-1.3284849) q[1];
sx q[1];
rz(2.8788024) q[1];
rz(-pi) q[2];
rz(0.73511519) q[3];
sx q[3];
rz(-0.8494091) q[3];
sx q[3];
rz(1.2934367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.7522493) q[2];
sx q[2];
rz(-1.8698591) q[2];
sx q[2];
rz(-2.6129369) q[2];
rz(-1.901249) q[3];
sx q[3];
rz(-2.7820008) q[3];
sx q[3];
rz(-2.8958029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56025958) q[0];
sx q[0];
rz(-1.154705) q[0];
sx q[0];
rz(0.2581968) q[0];
rz(-1.5064346) q[1];
sx q[1];
rz(-0.55748993) q[1];
sx q[1];
rz(2.7071276) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2697849) q[0];
sx q[0];
rz(-1.9315533) q[0];
sx q[0];
rz(-2.0478134) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.065504727) q[2];
sx q[2];
rz(-2.0993877) q[2];
sx q[2];
rz(2.5158109) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8409621) q[1];
sx q[1];
rz(-1.1197487) q[1];
sx q[1];
rz(1.8039963) q[1];
rz(-pi) q[2];
rz(-2.4934019) q[3];
sx q[3];
rz(-0.95578268) q[3];
sx q[3];
rz(-1.0685514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.42052856) q[2];
sx q[2];
rz(-1.3238182) q[2];
sx q[2];
rz(0.055796441) q[2];
rz(-2.2037286) q[3];
sx q[3];
rz(-0.28356975) q[3];
sx q[3];
rz(-0.24308932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.043512251) q[0];
sx q[0];
rz(-0.94399095) q[0];
sx q[0];
rz(-3.0010624) q[0];
rz(2.9699516) q[1];
sx q[1];
rz(-1.834603) q[1];
sx q[1];
rz(-2.8780639) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9517277) q[0];
sx q[0];
rz(-0.57706149) q[0];
sx q[0];
rz(0.18738562) q[0];
rz(-0.61435917) q[2];
sx q[2];
rz(-0.76459568) q[2];
sx q[2];
rz(1.489153) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9159689) q[1];
sx q[1];
rz(-1.9150503) q[1];
sx q[1];
rz(0.72361372) q[1];
rz(-0.88818355) q[3];
sx q[3];
rz(-2.7215951) q[3];
sx q[3];
rz(-2.1995467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.043109743) q[2];
sx q[2];
rz(-0.3442328) q[2];
sx q[2];
rz(-1.8919224) q[2];
rz(1.194687) q[3];
sx q[3];
rz(-2.1134816) q[3];
sx q[3];
rz(2.4627114) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93976218) q[0];
sx q[0];
rz(-1.8562466) q[0];
sx q[0];
rz(0.029065954) q[0];
rz(-1.4020231) q[1];
sx q[1];
rz(-0.35750917) q[1];
sx q[1];
rz(-0.32863858) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11360725) q[0];
sx q[0];
rz(-2.0562045) q[0];
sx q[0];
rz(3.0263682) q[0];
rz(-0.26942307) q[2];
sx q[2];
rz(-1.4002443) q[2];
sx q[2];
rz(1.5829057) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8935834) q[1];
sx q[1];
rz(-1.5072522) q[1];
sx q[1];
rz(0.43339543) q[1];
rz(-pi) q[2];
x q[2];
rz(0.5248431) q[3];
sx q[3];
rz(-2.3224324) q[3];
sx q[3];
rz(0.51845779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0236726) q[2];
sx q[2];
rz(-0.4549883) q[2];
sx q[2];
rz(-0.16061352) q[2];
rz(-1.9953856) q[3];
sx q[3];
rz(-1.8278154) q[3];
sx q[3];
rz(2.5207991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6495431) q[0];
sx q[0];
rz(-2.2821125) q[0];
sx q[0];
rz(0.78654003) q[0];
rz(-0.37711626) q[1];
sx q[1];
rz(-0.84723324) q[1];
sx q[1];
rz(0.4424817) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9839358) q[0];
sx q[0];
rz(-1.929495) q[0];
sx q[0];
rz(-1.1789807) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5673893) q[2];
sx q[2];
rz(-1.310537) q[2];
sx q[2];
rz(-2.0116218) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.46736273) q[1];
sx q[1];
rz(-1.2503887) q[1];
sx q[1];
rz(-0.25735374) q[1];
rz(-pi) q[2];
rz(3.0480012) q[3];
sx q[3];
rz(-1.8851265) q[3];
sx q[3];
rz(2.4623354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4601712) q[2];
sx q[2];
rz(-2.5632016) q[2];
sx q[2];
rz(2.1703413) q[2];
rz(-0.69532895) q[3];
sx q[3];
rz(-0.23614241) q[3];
sx q[3];
rz(-0.42738459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.5044395) q[0];
sx q[0];
rz(-1.210286) q[0];
sx q[0];
rz(-1.0094281) q[0];
rz(2.5908453) q[1];
sx q[1];
rz(-1.8661205) q[1];
sx q[1];
rz(-1.9783463) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79272072) q[0];
sx q[0];
rz(-0.36076818) q[0];
sx q[0];
rz(1.1611206) q[0];
rz(-pi) q[1];
rz(-2.8909056) q[2];
sx q[2];
rz(-0.87826585) q[2];
sx q[2];
rz(-0.78679774) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7904661) q[1];
sx q[1];
rz(-0.27517056) q[1];
sx q[1];
rz(2.5009584) q[1];
x q[2];
rz(2.8323152) q[3];
sx q[3];
rz(-0.89865548) q[3];
sx q[3];
rz(2.6291763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.78272468) q[2];
sx q[2];
rz(-1.1512558) q[2];
sx q[2];
rz(2.3434095) q[2];
rz(2.362137) q[3];
sx q[3];
rz(-2.6051086) q[3];
sx q[3];
rz(1.9688169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1786132) q[0];
sx q[0];
rz(-2.6586752) q[0];
sx q[0];
rz(-0.12284199) q[0];
rz(-0.12610647) q[1];
sx q[1];
rz(-1.6364731) q[1];
sx q[1];
rz(1.925148) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13623304) q[0];
sx q[0];
rz(-1.1598806) q[0];
sx q[0];
rz(-0.67816011) q[0];
rz(1.4037651) q[2];
sx q[2];
rz(-1.1208431) q[2];
sx q[2];
rz(-0.049875967) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.46978894) q[1];
sx q[1];
rz(-0.43699139) q[1];
sx q[1];
rz(0.48204084) q[1];
x q[2];
rz(2.5495278) q[3];
sx q[3];
rz(-1.4148303) q[3];
sx q[3];
rz(-1.1252943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8806261) q[2];
sx q[2];
rz(-0.61769056) q[2];
sx q[2];
rz(-0.043126062) q[2];
rz(-2.96636) q[3];
sx q[3];
rz(-0.8774811) q[3];
sx q[3];
rz(1.5847248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6373428) q[0];
sx q[0];
rz(-0.054333996) q[0];
sx q[0];
rz(2.9328226) q[0];
rz(1.6437221) q[1];
sx q[1];
rz(-1.6521963) q[1];
sx q[1];
rz(-2.0170905) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18692423) q[0];
sx q[0];
rz(-2.9483729) q[0];
sx q[0];
rz(0.78878553) q[0];
rz(1.3051946) q[2];
sx q[2];
rz(-2.2679066) q[2];
sx q[2];
rz(0.12649378) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3325602) q[1];
sx q[1];
rz(-3.0082015) q[1];
sx q[1];
rz(1.5619713) q[1];
rz(-pi) q[2];
x q[2];
rz(0.48891588) q[3];
sx q[3];
rz(-0.93547869) q[3];
sx q[3];
rz(-0.92632252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1407397) q[2];
sx q[2];
rz(-0.72379392) q[2];
sx q[2];
rz(-2.9635584) q[2];
rz(-1.5464276) q[3];
sx q[3];
rz(-1.3083357) q[3];
sx q[3];
rz(2.6509638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.35266018) q[0];
sx q[0];
rz(-1.1505609) q[0];
sx q[0];
rz(2.1066522) q[0];
rz(-2.3433698) q[1];
sx q[1];
rz(-2.1612576) q[1];
sx q[1];
rz(2.9916874) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4231739) q[0];
sx q[0];
rz(-2.383854) q[0];
sx q[0];
rz(-2.5916879) q[0];
x q[1];
rz(0.055585102) q[2];
sx q[2];
rz(-0.28700799) q[2];
sx q[2];
rz(-0.74319786) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.80430154) q[1];
sx q[1];
rz(-1.9807439) q[1];
sx q[1];
rz(1.0211584) q[1];
rz(1.8191765) q[3];
sx q[3];
rz(-2.2349149) q[3];
sx q[3];
rz(0.76457232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4118816) q[2];
sx q[2];
rz(-1.231266) q[2];
sx q[2];
rz(0.40851545) q[2];
rz(-2.216693) q[3];
sx q[3];
rz(-2.3291406) q[3];
sx q[3];
rz(0.48172054) q[3];
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
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9020486) q[0];
sx q[0];
rz(-1.5477381) q[0];
sx q[0];
rz(1.5292194) q[0];
rz(1.7655903) q[1];
sx q[1];
rz(-1.1767495) q[1];
sx q[1];
rz(-1.8935988) q[1];
rz(-0.31870141) q[2];
sx q[2];
rz(-2.3206884) q[2];
sx q[2];
rz(-0.86736292) q[2];
rz(-1.3281214) q[3];
sx q[3];
rz(-0.96488733) q[3];
sx q[3];
rz(0.78028954) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
