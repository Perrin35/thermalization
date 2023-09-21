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
rz(1.8058864) q[1];
sx q[1];
rz(3.4808573) q[1];
sx q[1];
rz(9.1453287) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1301491) q[0];
sx q[0];
rz(-0.37405095) q[0];
sx q[0];
rz(2.1610297) q[0];
rz(-0.7819671) q[2];
sx q[2];
rz(-1.3227533) q[2];
sx q[2];
rz(-2.9458407) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.93912032) q[1];
sx q[1];
rz(-1.7608374) q[1];
sx q[1];
rz(0.01357667) q[1];
x q[2];
rz(-1.6104326) q[3];
sx q[3];
rz(-2.2124024) q[3];
sx q[3];
rz(1.8455781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.44148579) q[2];
sx q[2];
rz(-1.4725279) q[2];
sx q[2];
rz(2.1842365) q[2];
rz(-2.8422614) q[3];
sx q[3];
rz(-2.744031) q[3];
sx q[3];
rz(-0.41199747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9807724) q[0];
sx q[0];
rz(-2.1919724) q[0];
sx q[0];
rz(-0.41369307) q[0];
rz(-1.7970239) q[1];
sx q[1];
rz(-2.3584056) q[1];
sx q[1];
rz(-2.5059674) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6247511) q[0];
sx q[0];
rz(-0.13953129) q[0];
sx q[0];
rz(1.2463039) q[0];
x q[1];
rz(-3.0947729) q[2];
sx q[2];
rz(-1.4390107) q[2];
sx q[2];
rz(1.5891967) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8860652) q[1];
sx q[1];
rz(-0.35554245) q[1];
sx q[1];
rz(-2.38106) q[1];
rz(-pi) q[2];
x q[2];
rz(0.91931822) q[3];
sx q[3];
rz(-0.98005664) q[3];
sx q[3];
rz(-2.7880993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.7522493) q[2];
sx q[2];
rz(-1.2717335) q[2];
sx q[2];
rz(0.52865571) q[2];
rz(1.2403437) q[3];
sx q[3];
rz(-0.35959187) q[3];
sx q[3];
rz(2.8958029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56025958) q[0];
sx q[0];
rz(-1.9868877) q[0];
sx q[0];
rz(2.8833959) q[0];
rz(1.6351581) q[1];
sx q[1];
rz(-0.55748993) q[1];
sx q[1];
rz(2.7071276) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8718078) q[0];
sx q[0];
rz(-1.2100394) q[0];
sx q[0];
rz(-2.0478134) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6824109) q[2];
sx q[2];
rz(-2.6093405) q[2];
sx q[2];
rz(0.75512952) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1669958) q[1];
sx q[1];
rz(-1.3612862) q[1];
sx q[1];
rz(0.4619044) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.86359777) q[3];
sx q[3];
rz(-2.2798385) q[3];
sx q[3];
rz(1.1532702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.42052856) q[2];
sx q[2];
rz(-1.3238182) q[2];
sx q[2];
rz(0.055796441) q[2];
rz(-2.2037286) q[3];
sx q[3];
rz(-2.8580229) q[3];
sx q[3];
rz(0.24308932) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.043512251) q[0];
sx q[0];
rz(-2.1976017) q[0];
sx q[0];
rz(-0.14053024) q[0];
rz(0.17164104) q[1];
sx q[1];
rz(-1.3069897) q[1];
sx q[1];
rz(-2.8780639) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9182501) q[0];
sx q[0];
rz(-1.6726057) q[0];
sx q[0];
rz(0.56901594) q[0];
x q[1];
rz(1.0657004) q[2];
sx q[2];
rz(-0.96955883) q[2];
sx q[2];
rz(0.87841735) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.22562379) q[1];
sx q[1];
rz(-1.2265424) q[1];
sx q[1];
rz(-0.72361372) q[1];
x q[2];
rz(-1.9043546) q[3];
sx q[3];
rz(-1.8309438) q[3];
sx q[3];
rz(3.1317657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.043109743) q[2];
sx q[2];
rz(-0.3442328) q[2];
sx q[2];
rz(-1.2496703) q[2];
rz(-1.194687) q[3];
sx q[3];
rz(-2.1134816) q[3];
sx q[3];
rz(0.67888129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2018305) q[0];
sx q[0];
rz(-1.8562466) q[0];
sx q[0];
rz(-0.029065954) q[0];
rz(-1.7395696) q[1];
sx q[1];
rz(-0.35750917) q[1];
sx q[1];
rz(-2.8129541) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5111361) q[0];
sx q[0];
rz(-1.4689313) q[0];
sx q[0];
rz(1.0826375) q[0];
rz(-pi) q[1];
rz(0.57428898) q[2];
sx q[2];
rz(-0.31775489) q[2];
sx q[2];
rz(0.53900915) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.45916522) q[1];
sx q[1];
rz(-2.703856) q[1];
sx q[1];
rz(-0.15037219) q[1];
x q[2];
rz(2.0629115) q[3];
sx q[3];
rz(-2.2552367) q[3];
sx q[3];
rz(1.9198315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.11792004) q[2];
sx q[2];
rz(-2.6866044) q[2];
sx q[2];
rz(-2.9809791) q[2];
rz(1.9953856) q[3];
sx q[3];
rz(-1.8278154) q[3];
sx q[3];
rz(-2.5207991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6495431) q[0];
sx q[0];
rz(-0.8594802) q[0];
sx q[0];
rz(2.3550526) q[0];
rz(-0.37711626) q[1];
sx q[1];
rz(-2.2943594) q[1];
sx q[1];
rz(2.699111) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0241644) q[0];
sx q[0];
rz(-0.52485835) q[0];
sx q[0];
rz(0.79458046) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2636678) q[2];
sx q[2];
rz(-1.0182292) q[2];
sx q[2];
rz(2.8657258) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9554555) q[1];
sx q[1];
rz(-1.326814) q[1];
sx q[1];
rz(-1.9013491) q[1];
rz(-pi) q[2];
rz(-1.8864185) q[3];
sx q[3];
rz(-1.6597897) q[3];
sx q[3];
rz(-0.86252585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4601712) q[2];
sx q[2];
rz(-2.5632016) q[2];
sx q[2];
rz(-2.1703413) q[2];
rz(0.69532895) q[3];
sx q[3];
rz(-2.9054502) q[3];
sx q[3];
rz(-0.42738459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63715315) q[0];
sx q[0];
rz(-1.9313066) q[0];
sx q[0];
rz(2.1321645) q[0];
rz(-2.5908453) q[1];
sx q[1];
rz(-1.8661205) q[1];
sx q[1];
rz(1.9783463) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9775987) q[0];
sx q[0];
rz(-1.4297276) q[0];
sx q[0];
rz(1.903957) q[0];
rz(-pi) q[1];
rz(1.2802358) q[2];
sx q[2];
rz(-2.4121948) q[2];
sx q[2];
rz(1.1682208) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1313596) q[1];
sx q[1];
rz(-1.3512003) q[1];
sx q[1];
rz(1.7379727) q[1];
rz(-0.30927741) q[3];
sx q[3];
rz(-0.89865548) q[3];
sx q[3];
rz(2.6291763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.78272468) q[2];
sx q[2];
rz(-1.1512558) q[2];
sx q[2];
rz(-0.7981832) q[2];
rz(2.362137) q[3];
sx q[3];
rz(-2.6051086) q[3];
sx q[3];
rz(-1.1727758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1786132) q[0];
sx q[0];
rz(-0.48291746) q[0];
sx q[0];
rz(0.12284199) q[0];
rz(0.12610647) q[1];
sx q[1];
rz(-1.5051196) q[1];
sx q[1];
rz(1.925148) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0053596) q[0];
sx q[0];
rz(-1.1598806) q[0];
sx q[0];
rz(0.67816011) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33151303) q[2];
sx q[2];
rz(-2.6636332) q[2];
sx q[2];
rz(-0.31994672) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6718037) q[1];
sx q[1];
rz(-2.7046013) q[1];
sx q[1];
rz(2.6595518) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5495278) q[3];
sx q[3];
rz(-1.7267623) q[3];
sx q[3];
rz(1.1252943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8806261) q[2];
sx q[2];
rz(-0.61769056) q[2];
sx q[2];
rz(3.0984666) q[2];
rz(-0.17523266) q[3];
sx q[3];
rz(-0.8774811) q[3];
sx q[3];
rz(-1.5847248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-0.50424987) q[0];
sx q[0];
rz(-0.054333996) q[0];
sx q[0];
rz(0.20877008) q[0];
rz(1.6437221) q[1];
sx q[1];
rz(-1.4893963) q[1];
sx q[1];
rz(2.0170905) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97832699) q[0];
sx q[0];
rz(-1.4341337) q[0];
sx q[0];
rz(3.004573) q[0];
rz(-pi) q[1];
rz(1.3051946) q[2];
sx q[2];
rz(-0.87368602) q[2];
sx q[2];
rz(-0.12649378) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3885755) q[1];
sx q[1];
rz(-1.57197) q[1];
sx q[1];
rz(1.4374103) q[1];
x q[2];
rz(2.2665958) q[3];
sx q[3];
rz(-1.1831302) q[3];
sx q[3];
rz(2.8029203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1407397) q[2];
sx q[2];
rz(-0.72379392) q[2];
sx q[2];
rz(0.17803426) q[2];
rz(1.5464276) q[3];
sx q[3];
rz(-1.3083357) q[3];
sx q[3];
rz(0.49062887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35266018) q[0];
sx q[0];
rz(-1.1505609) q[0];
sx q[0];
rz(2.1066522) q[0];
rz(-0.79822284) q[1];
sx q[1];
rz(-0.98033506) q[1];
sx q[1];
rz(-0.14990526) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56652727) q[0];
sx q[0];
rz(-1.9381822) q[0];
sx q[0];
rz(0.67879403) q[0];
rz(-0.055585102) q[2];
sx q[2];
rz(-0.28700799) q[2];
sx q[2];
rz(-2.3983948) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.18969892) q[1];
sx q[1];
rz(-0.67283291) q[1];
sx q[1];
rz(-0.876902) q[1];
rz(0.30431872) q[3];
sx q[3];
rz(-0.70239866) q[3];
sx q[3];
rz(1.1549293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4118816) q[2];
sx q[2];
rz(-1.9103266) q[2];
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
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2395441) q[0];
sx q[0];
rz(-1.5477381) q[0];
sx q[0];
rz(1.5292194) q[0];
rz(1.3760024) q[1];
sx q[1];
rz(-1.9648432) q[1];
sx q[1];
rz(1.2479938) q[1];
rz(1.2462763) q[2];
sx q[2];
rz(-0.80249716) q[2];
sx q[2];
rz(-1.3182166) q[2];
rz(2.8077447) q[3];
sx q[3];
rz(-0.6469938) q[3];
sx q[3];
rz(-2.7713431) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
