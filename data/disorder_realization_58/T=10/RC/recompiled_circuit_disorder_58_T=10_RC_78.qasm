OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9392202) q[0];
sx q[0];
rz(-0.4063172) q[0];
sx q[0];
rz(0.82011861) q[0];
rz(2.7804873) q[1];
sx q[1];
rz(-0.63280025) q[1];
sx q[1];
rz(-0.83067218) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93847371) q[0];
sx q[0];
rz(-1.1955368) q[0];
sx q[0];
rz(1.9012326) q[0];
rz(-2.6248706) q[2];
sx q[2];
rz(-2.3331254) q[2];
sx q[2];
rz(-0.39743039) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.1931397) q[1];
sx q[1];
rz(-1.8857191) q[1];
sx q[1];
rz(-0.83955168) q[1];
rz(-pi) q[2];
rz(-1.9740231) q[3];
sx q[3];
rz(-1.9864169) q[3];
sx q[3];
rz(2.2574539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7044907) q[2];
sx q[2];
rz(-0.4318684) q[2];
sx q[2];
rz(-3.0214156) q[2];
rz(-1.1581356) q[3];
sx q[3];
rz(-1.3985671) q[3];
sx q[3];
rz(-2.2226298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.08081089) q[0];
sx q[0];
rz(-1.7827543) q[0];
sx q[0];
rz(-0.91180116) q[0];
rz(0.78951019) q[1];
sx q[1];
rz(-0.98840886) q[1];
sx q[1];
rz(-0.3266913) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.806843) q[0];
sx q[0];
rz(-2.202889) q[0];
sx q[0];
rz(0.64081162) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8961043) q[2];
sx q[2];
rz(-2.6763958) q[2];
sx q[2];
rz(2.543769) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.82636278) q[1];
sx q[1];
rz(-0.76347199) q[1];
sx q[1];
rz(-2.1748494) q[1];
rz(-pi) q[2];
rz(-2.2585906) q[3];
sx q[3];
rz(-1.6371173) q[3];
sx q[3];
rz(-1.6460713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6313173) q[2];
sx q[2];
rz(-0.77284208) q[2];
sx q[2];
rz(1.2871683) q[2];
rz(-3.0316947) q[3];
sx q[3];
rz(-1.7307614) q[3];
sx q[3];
rz(1.7597642) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.592955) q[0];
sx q[0];
rz(-2.4023963) q[0];
sx q[0];
rz(0.32989311) q[0];
rz(0.27711162) q[1];
sx q[1];
rz(-1.3169293) q[1];
sx q[1];
rz(1.057391) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0111423) q[0];
sx q[0];
rz(-1.5913977) q[0];
sx q[0];
rz(-1.9111454) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.51809394) q[2];
sx q[2];
rz(-0.38803852) q[2];
sx q[2];
rz(-2.5990017) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2705546) q[1];
sx q[1];
rz(-1.241031) q[1];
sx q[1];
rz(2.2689181) q[1];
rz(-pi) q[2];
x q[2];
rz(0.40773817) q[3];
sx q[3];
rz(-0.45555112) q[3];
sx q[3];
rz(1.3726485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.7827591) q[2];
sx q[2];
rz(-1.1153355) q[2];
sx q[2];
rz(-1.7791629) q[2];
rz(-2.5168915) q[3];
sx q[3];
rz(-2.0420859) q[3];
sx q[3];
rz(2.3220298) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7841004) q[0];
sx q[0];
rz(-0.52629137) q[0];
sx q[0];
rz(0.5350565) q[0];
rz(-2.0013981) q[1];
sx q[1];
rz(-1.813872) q[1];
sx q[1];
rz(2.9761956) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9388158) q[0];
sx q[0];
rz(-1.6949777) q[0];
sx q[0];
rz(-1.3846272) q[0];
rz(-pi) q[1];
rz(2.303316) q[2];
sx q[2];
rz(-2.8998313) q[2];
sx q[2];
rz(-1.6844815) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.0488102) q[1];
sx q[1];
rz(-0.65199344) q[1];
sx q[1];
rz(0.55744967) q[1];
x q[2];
rz(1.1449279) q[3];
sx q[3];
rz(-0.59755675) q[3];
sx q[3];
rz(1.715341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.39607221) q[2];
sx q[2];
rz(-1.3362276) q[2];
sx q[2];
rz(-2.9361434) q[2];
rz(2.0139587) q[3];
sx q[3];
rz(-1.1536359) q[3];
sx q[3];
rz(1.2566465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66185343) q[0];
sx q[0];
rz(-2.2275708) q[0];
sx q[0];
rz(-0.35476312) q[0];
rz(-1.154249) q[1];
sx q[1];
rz(-0.92461363) q[1];
sx q[1];
rz(-0.23194557) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1922798) q[0];
sx q[0];
rz(-2.4711907) q[0];
sx q[0];
rz(2.288726) q[0];
rz(-pi) q[1];
rz(-2.4460692) q[2];
sx q[2];
rz(-1.4624274) q[2];
sx q[2];
rz(-2.7930789) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8307454) q[1];
sx q[1];
rz(-1.3535045) q[1];
sx q[1];
rz(2.9160935) q[1];
x q[2];
rz(-2.7730745) q[3];
sx q[3];
rz(-1.5526062) q[3];
sx q[3];
rz(-0.25572488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.140124) q[2];
sx q[2];
rz(-1.3342369) q[2];
sx q[2];
rz(1.9011964) q[2];
rz(-0.59605789) q[3];
sx q[3];
rz(-1.3052992) q[3];
sx q[3];
rz(-1.8381455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0971138) q[0];
sx q[0];
rz(-3.0713186) q[0];
sx q[0];
rz(-0.20275673) q[0];
rz(-2.1525106) q[1];
sx q[1];
rz(-1.443807) q[1];
sx q[1];
rz(-0.99745497) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.049012262) q[0];
sx q[0];
rz(-0.870734) q[0];
sx q[0];
rz(-0.17605619) q[0];
rz(-1.9428271) q[2];
sx q[2];
rz(-2.2673006) q[2];
sx q[2];
rz(-0.55559413) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7194697) q[1];
sx q[1];
rz(-2.649253) q[1];
sx q[1];
rz(1.0455529) q[1];
x q[2];
rz(2.698425) q[3];
sx q[3];
rz(-2.245095) q[3];
sx q[3];
rz(2.5600195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9138907) q[2];
sx q[2];
rz(-1.1662377) q[2];
sx q[2];
rz(2.690199) q[2];
rz(0.40870062) q[3];
sx q[3];
rz(-1.6025851) q[3];
sx q[3];
rz(-1.2020948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.8354427) q[0];
sx q[0];
rz(-0.4040443) q[0];
sx q[0];
rz(2.5174482) q[0];
rz(1.5165326) q[1];
sx q[1];
rz(-0.25696483) q[1];
sx q[1];
rz(-0.61378941) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95425883) q[0];
sx q[0];
rz(-1.3849392) q[0];
sx q[0];
rz(2.3128187) q[0];
x q[1];
rz(1.6348398) q[2];
sx q[2];
rz(-1.8539068) q[2];
sx q[2];
rz(1.7664906) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.92717273) q[1];
sx q[1];
rz(-1.629429) q[1];
sx q[1];
rz(3.0374132) q[1];
x q[2];
rz(2.036014) q[3];
sx q[3];
rz(-1.675616) q[3];
sx q[3];
rz(0.011381586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.43341407) q[2];
sx q[2];
rz(-0.91579473) q[2];
sx q[2];
rz(-1.9667352) q[2];
rz(-2.5332149) q[3];
sx q[3];
rz(-1.404168) q[3];
sx q[3];
rz(-1.4233937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90010086) q[0];
sx q[0];
rz(-1.2986203) q[0];
sx q[0];
rz(-2.6532145) q[0];
rz(1.5178559) q[1];
sx q[1];
rz(-1.3987712) q[1];
sx q[1];
rz(0.98446313) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9280633) q[0];
sx q[0];
rz(-2.5941879) q[0];
sx q[0];
rz(-0.071896032) q[0];
rz(0.77163561) q[2];
sx q[2];
rz(-1.2727591) q[2];
sx q[2];
rz(0.76921295) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9148548) q[1];
sx q[1];
rz(-1.2054772) q[1];
sx q[1];
rz(1.564333) q[1];
x q[2];
rz(-2.5823309) q[3];
sx q[3];
rz(-0.54702938) q[3];
sx q[3];
rz(-0.83838851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.70665923) q[2];
sx q[2];
rz(-1.443053) q[2];
sx q[2];
rz(-0.53517503) q[2];
rz(-2.0914071) q[3];
sx q[3];
rz(-1.340056) q[3];
sx q[3];
rz(2.1323269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
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
rz(-2.500279) q[0];
sx q[0];
rz(-2.2135493) q[0];
sx q[0];
rz(-0.6859268) q[0];
rz(-2.7507239) q[1];
sx q[1];
rz(-2.1978244) q[1];
sx q[1];
rz(2.2156782) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5110916) q[0];
sx q[0];
rz(-1.7112268) q[0];
sx q[0];
rz(-1.7846084) q[0];
rz(-0.7943031) q[2];
sx q[2];
rz(-1.0216733) q[2];
sx q[2];
rz(-1.4505475) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8402108) q[1];
sx q[1];
rz(-1.7173319) q[1];
sx q[1];
rz(2.0135897) q[1];
rz(-pi) q[2];
x q[2];
rz(0.52243201) q[3];
sx q[3];
rz(-0.2345095) q[3];
sx q[3];
rz(1.0684551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.19568504) q[2];
sx q[2];
rz(-0.47161272) q[2];
sx q[2];
rz(-0.53608981) q[2];
rz(0.4195956) q[3];
sx q[3];
rz(-0.59046888) q[3];
sx q[3];
rz(-1.6528116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0855899) q[0];
sx q[0];
rz(-2.7828126) q[0];
sx q[0];
rz(0.39500239) q[0];
rz(1.6292054) q[1];
sx q[1];
rz(-1.869166) q[1];
sx q[1];
rz(1.013247) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2241867) q[0];
sx q[0];
rz(-1.4316443) q[0];
sx q[0];
rz(-1.3361206) q[0];
rz(-1.6147862) q[2];
sx q[2];
rz(-1.9242052) q[2];
sx q[2];
rz(2.8994438) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.71868616) q[1];
sx q[1];
rz(-1.4740853) q[1];
sx q[1];
rz(-0.42214091) q[1];
rz(-pi) q[2];
rz(0.69443955) q[3];
sx q[3];
rz(-0.89530066) q[3];
sx q[3];
rz(-1.6758067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6932678) q[2];
sx q[2];
rz(-1.5193181) q[2];
sx q[2];
rz(2.3821793) q[2];
rz(-1.365472) q[3];
sx q[3];
rz(-2.0335734) q[3];
sx q[3];
rz(-2.571648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7286745) q[0];
sx q[0];
rz(-1.2644132) q[0];
sx q[0];
rz(1.2046474) q[0];
rz(-2.5683174) q[1];
sx q[1];
rz(-1.9075867) q[1];
sx q[1];
rz(1.8214068) q[1];
rz(-1.1258833) q[2];
sx q[2];
rz(-1.5013668) q[2];
sx q[2];
rz(-2.8253386) q[2];
rz(0.063354062) q[3];
sx q[3];
rz(-0.8739211) q[3];
sx q[3];
rz(-2.7779761) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];