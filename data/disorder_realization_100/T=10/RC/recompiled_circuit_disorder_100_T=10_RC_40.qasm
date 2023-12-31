OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.93958062) q[0];
sx q[0];
rz(-0.35020819) q[0];
sx q[0];
rz(2.7749618) q[0];
rz(0.86759138) q[1];
sx q[1];
rz(-2.4974513) q[1];
sx q[1];
rz(1.4555414) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.769387) q[0];
sx q[0];
rz(-1.9908449) q[0];
sx q[0];
rz(0.37680349) q[0];
rz(-pi) q[1];
rz(-0.86906616) q[2];
sx q[2];
rz(-1.7742187) q[2];
sx q[2];
rz(-0.8806526) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9526706) q[1];
sx q[1];
rz(-2.1181207) q[1];
sx q[1];
rz(0.31618936) q[1];
rz(-pi) q[2];
rz(-1.5600169) q[3];
sx q[3];
rz(-1.618715) q[3];
sx q[3];
rz(-0.77424327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.1047487) q[2];
sx q[2];
rz(-1.3354744) q[2];
sx q[2];
rz(2.74995) q[2];
rz(-3.0018905) q[3];
sx q[3];
rz(-2.468686) q[3];
sx q[3];
rz(-0.02123775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.7249107) q[0];
sx q[0];
rz(-1.4044489) q[0];
sx q[0];
rz(-1.2765983) q[0];
rz(2.2712767) q[1];
sx q[1];
rz(-1.5753997) q[1];
sx q[1];
rz(-1.2044027) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6635839) q[0];
sx q[0];
rz(-1.1482129) q[0];
sx q[0];
rz(2.708486) q[0];
rz(0.18560974) q[2];
sx q[2];
rz(-0.87366548) q[2];
sx q[2];
rz(2.6999161) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.345563) q[1];
sx q[1];
rz(-1.9351164) q[1];
sx q[1];
rz(1.6771392) q[1];
x q[2];
rz(2.2435917) q[3];
sx q[3];
rz(-1.3452531) q[3];
sx q[3];
rz(-2.0190092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5045972) q[2];
sx q[2];
rz(-2.6837139) q[2];
sx q[2];
rz(-1.9223928) q[2];
rz(-2.7870264) q[3];
sx q[3];
rz(-2.6597326) q[3];
sx q[3];
rz(1.5725117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5511659) q[0];
sx q[0];
rz(-1.5274436) q[0];
sx q[0];
rz(-2.5235126) q[0];
rz(-0.016050054) q[1];
sx q[1];
rz(-0.87688223) q[1];
sx q[1];
rz(-1.191167) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0389827) q[0];
sx q[0];
rz(-1.257574) q[0];
sx q[0];
rz(0.15977504) q[0];
x q[1];
rz(0.70706681) q[2];
sx q[2];
rz(-2.7756049) q[2];
sx q[2];
rz(2.6235839) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1515854) q[1];
sx q[1];
rz(-1.3602435) q[1];
sx q[1];
rz(1.3596119) q[1];
x q[2];
rz(1.1690087) q[3];
sx q[3];
rz(-2.0766633) q[3];
sx q[3];
rz(-1.5617621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.96737343) q[2];
sx q[2];
rz(-0.9625532) q[2];
sx q[2];
rz(1.013914) q[2];
rz(1.7845456) q[3];
sx q[3];
rz(-0.8299399) q[3];
sx q[3];
rz(-2.1092265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33092609) q[0];
sx q[0];
rz(-1.4241011) q[0];
sx q[0];
rz(2.9372835) q[0];
rz(1.7640242) q[1];
sx q[1];
rz(-1.4148477) q[1];
sx q[1];
rz(-0.92299443) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6072236) q[0];
sx q[0];
rz(-1.9019706) q[0];
sx q[0];
rz(1.2481199) q[0];
x q[1];
rz(0.26206215) q[2];
sx q[2];
rz(-0.045496551) q[2];
sx q[2];
rz(0.87693518) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1624565) q[1];
sx q[1];
rz(-2.7297331) q[1];
sx q[1];
rz(1.3084175) q[1];
rz(0.022425671) q[3];
sx q[3];
rz(-2.2750345) q[3];
sx q[3];
rz(1.1426329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4541645) q[2];
sx q[2];
rz(-1.585008) q[2];
sx q[2];
rz(-2.6320809) q[2];
rz(0.64309684) q[3];
sx q[3];
rz(-1.0984848) q[3];
sx q[3];
rz(-0.96737868) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61280695) q[0];
sx q[0];
rz(-1.2061773) q[0];
sx q[0];
rz(-2.2221785) q[0];
rz(2.1557504) q[1];
sx q[1];
rz(-1.7835833) q[1];
sx q[1];
rz(-1.7808419) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4454942) q[0];
sx q[0];
rz(-1.3514263) q[0];
sx q[0];
rz(1.7646837) q[0];
x q[1];
rz(0.76061337) q[2];
sx q[2];
rz(-1.0906272) q[2];
sx q[2];
rz(-1.7825356) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9779224) q[1];
sx q[1];
rz(-1.1174035) q[1];
sx q[1];
rz(1.8015253) q[1];
rz(2.1702483) q[3];
sx q[3];
rz(-1.506862) q[3];
sx q[3];
rz(-1.5980699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.78648606) q[2];
sx q[2];
rz(-1.3856709) q[2];
sx q[2];
rz(0.11486593) q[2];
rz(2.6857175) q[3];
sx q[3];
rz(-1.3331648) q[3];
sx q[3];
rz(1.8615287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8961261) q[0];
sx q[0];
rz(-1.8969314) q[0];
sx q[0];
rz(-1.2448357) q[0];
rz(-2.4339829) q[1];
sx q[1];
rz(-1.5039624) q[1];
sx q[1];
rz(-2.8663666) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21101418) q[0];
sx q[0];
rz(-2.0654581) q[0];
sx q[0];
rz(2.7633694) q[0];
x q[1];
rz(-2.6456344) q[2];
sx q[2];
rz(-2.1514116) q[2];
sx q[2];
rz(2.1883287) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.59906193) q[1];
sx q[1];
rz(-1.2236569) q[1];
sx q[1];
rz(0.51862006) q[1];
rz(-pi) q[2];
rz(-1.4554548) q[3];
sx q[3];
rz(-0.62386688) q[3];
sx q[3];
rz(1.3184402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.39367166) q[2];
sx q[2];
rz(-1.5449521) q[2];
sx q[2];
rz(0.45864027) q[2];
rz(0.7115055) q[3];
sx q[3];
rz(-2.4578874) q[3];
sx q[3];
rz(2.5512364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8608619) q[0];
sx q[0];
rz(-1.9545398) q[0];
sx q[0];
rz(-1.77805) q[0];
rz(-1.4200312) q[1];
sx q[1];
rz(-2.6224711) q[1];
sx q[1];
rz(-2.4218959) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73368209) q[0];
sx q[0];
rz(-2.1955197) q[0];
sx q[0];
rz(1.2172933) q[0];
x q[1];
rz(-2.4160956) q[2];
sx q[2];
rz(-2.2273387) q[2];
sx q[2];
rz(3.0692435) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.81855782) q[1];
sx q[1];
rz(-2.5948988) q[1];
sx q[1];
rz(-0.085303765) q[1];
rz(-2.8119874) q[3];
sx q[3];
rz(-2.4334987) q[3];
sx q[3];
rz(3.0081188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.482243) q[2];
sx q[2];
rz(-1.5321926) q[2];
sx q[2];
rz(-2.4534295) q[2];
rz(0.041953772) q[3];
sx q[3];
rz(-1.099702) q[3];
sx q[3];
rz(-0.28013128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2280837) q[0];
sx q[0];
rz(-1.1627473) q[0];
sx q[0];
rz(-0.18653175) q[0];
rz(0.60449156) q[1];
sx q[1];
rz(-2.1291321) q[1];
sx q[1];
rz(1.4950745) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.625531) q[0];
sx q[0];
rz(-0.4058668) q[0];
sx q[0];
rz(-2.5748475) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.375884) q[2];
sx q[2];
rz(-0.88390985) q[2];
sx q[2];
rz(-1.9538823) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1284359) q[1];
sx q[1];
rz(-1.6746579) q[1];
sx q[1];
rz(2.8810112) q[1];
rz(-pi) q[2];
x q[2];
rz(0.80963482) q[3];
sx q[3];
rz(-1.0757043) q[3];
sx q[3];
rz(-1.1822869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8994393) q[2];
sx q[2];
rz(-0.95726761) q[2];
sx q[2];
rz(0.014766679) q[2];
rz(1.2773369) q[3];
sx q[3];
rz(-1.3093964) q[3];
sx q[3];
rz(1.7285255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16167851) q[0];
sx q[0];
rz(-2.4665687) q[0];
sx q[0];
rz(-2.263608) q[0];
rz(-0.19628482) q[1];
sx q[1];
rz(-1.1801964) q[1];
sx q[1];
rz(-1.7810129) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1454865) q[0];
sx q[0];
rz(-1.6084451) q[0];
sx q[0];
rz(-2.1426175) q[0];
x q[1];
rz(0.94366818) q[2];
sx q[2];
rz(-1.6086279) q[2];
sx q[2];
rz(0.19247069) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4728002) q[1];
sx q[1];
rz(-1.2392514) q[1];
sx q[1];
rz(-1.214295) q[1];
rz(2.9865773) q[3];
sx q[3];
rz(-1.985637) q[3];
sx q[3];
rz(-0.14839867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.62375162) q[2];
sx q[2];
rz(-1.5072284) q[2];
sx q[2];
rz(-0.8927792) q[2];
rz(0.039285224) q[3];
sx q[3];
rz(-1.6623442) q[3];
sx q[3];
rz(-0.75004309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84353012) q[0];
sx q[0];
rz(-3.1382882) q[0];
sx q[0];
rz(0.046534006) q[0];
rz(-2.2325366) q[1];
sx q[1];
rz(-2.2687056) q[1];
sx q[1];
rz(-2.4216901) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4601446) q[0];
sx q[0];
rz(-1.8143192) q[0];
sx q[0];
rz(0.081865099) q[0];
x q[1];
rz(0.30883046) q[2];
sx q[2];
rz(-1.5656106) q[2];
sx q[2];
rz(-1.3877102) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7457909) q[1];
sx q[1];
rz(-1.5862641) q[1];
sx q[1];
rz(-1.8176816) q[1];
x q[2];
rz(-2.4107237) q[3];
sx q[3];
rz(-1.0167529) q[3];
sx q[3];
rz(2.4027367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7211192) q[2];
sx q[2];
rz(-1.3995918) q[2];
sx q[2];
rz(3.1151248) q[2];
rz(-1.8376393) q[3];
sx q[3];
rz(-0.81726685) q[3];
sx q[3];
rz(1.8855689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3532886) q[0];
sx q[0];
rz(-1.7699387) q[0];
sx q[0];
rz(1.8040245) q[0];
rz(-2.9251255) q[1];
sx q[1];
rz(-1.4420061) q[1];
sx q[1];
rz(-1.9056086) q[1];
rz(0.36454501) q[2];
sx q[2];
rz(-0.74082965) q[2];
sx q[2];
rz(-0.34168591) q[2];
rz(0.032565928) q[3];
sx q[3];
rz(-2.1571772) q[3];
sx q[3];
rz(-1.0909506) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
