OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.206267) q[0];
sx q[0];
rz(-1.6214108) q[0];
sx q[0];
rz(1.4186463) q[0];
rz(-3.0980134) q[1];
sx q[1];
rz(-1.8585304) q[1];
sx q[1];
rz(-0.15951523) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43082008) q[0];
sx q[0];
rz(-0.16834591) q[0];
sx q[0];
rz(1.400186) q[0];
x q[1];
rz(2.5297727) q[2];
sx q[2];
rz(-1.7828336) q[2];
sx q[2];
rz(3.0481047) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.91640305) q[1];
sx q[1];
rz(-1.4328416) q[1];
sx q[1];
rz(-1.0911343) q[1];
rz(-pi) q[2];
rz(-0.84955022) q[3];
sx q[3];
rz(-1.2952943) q[3];
sx q[3];
rz(-1.4412137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.45304766) q[2];
sx q[2];
rz(-1.9828372) q[2];
sx q[2];
rz(-3.0513406) q[2];
rz(0.074617535) q[3];
sx q[3];
rz(-1.4679694) q[3];
sx q[3];
rz(2.727865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8892882) q[0];
sx q[0];
rz(-1.9828718) q[0];
sx q[0];
rz(-1.4537551) q[0];
rz(0.74268913) q[1];
sx q[1];
rz(-1.7748666) q[1];
sx q[1];
rz(-2.3765366) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3111717) q[0];
sx q[0];
rz(-0.97872916) q[0];
sx q[0];
rz(-1.4903846) q[0];
rz(-pi) q[1];
rz(2.3397331) q[2];
sx q[2];
rz(-2.2727786) q[2];
sx q[2];
rz(-0.02846708) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.88095746) q[1];
sx q[1];
rz(-0.85399125) q[1];
sx q[1];
rz(2.6997801) q[1];
rz(-pi) q[2];
rz(0.54799796) q[3];
sx q[3];
rz(-1.573085) q[3];
sx q[3];
rz(1.8016004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.226977) q[2];
sx q[2];
rz(-1.0116297) q[2];
sx q[2];
rz(-1.4545308) q[2];
rz(0.64676532) q[3];
sx q[3];
rz(-2.5819467) q[3];
sx q[3];
rz(-1.8618934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0890559) q[0];
sx q[0];
rz(-1.4707969) q[0];
sx q[0];
rz(-0.045150969) q[0];
rz(-1.2308944) q[1];
sx q[1];
rz(-1.8321593) q[1];
sx q[1];
rz(1.0567573) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3927944) q[0];
sx q[0];
rz(-0.65750376) q[0];
sx q[0];
rz(-2.1211521) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0750404) q[2];
sx q[2];
rz(-0.98598749) q[2];
sx q[2];
rz(2.5677447) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.31683982) q[1];
sx q[1];
rz(-0.72622967) q[1];
sx q[1];
rz(1.9588425) q[1];
rz(-pi) q[2];
rz(1.2892492) q[3];
sx q[3];
rz(-2.4334638) q[3];
sx q[3];
rz(-2.331108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7354273) q[2];
sx q[2];
rz(-2.1801517) q[2];
sx q[2];
rz(0.65262922) q[2];
rz(-1.8485707) q[3];
sx q[3];
rz(-0.76904622) q[3];
sx q[3];
rz(-3.0434216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80766135) q[0];
sx q[0];
rz(-0.46285358) q[0];
sx q[0];
rz(-1.7474784) q[0];
rz(-1.8380503) q[1];
sx q[1];
rz(-1.9775672) q[1];
sx q[1];
rz(2.0533144) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5096579) q[0];
sx q[0];
rz(-2.3096173) q[0];
sx q[0];
rz(-0.37507071) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5456504) q[2];
sx q[2];
rz(-0.94178994) q[2];
sx q[2];
rz(3.0160835) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.478178) q[1];
sx q[1];
rz(-2.3199548) q[1];
sx q[1];
rz(-2.4577702) q[1];
rz(-3.0451891) q[3];
sx q[3];
rz(-2.6402874) q[3];
sx q[3];
rz(0.65602428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.7465705) q[2];
sx q[2];
rz(-1.7638548) q[2];
sx q[2];
rz(-2.6611888) q[2];
rz(0.26614842) q[3];
sx q[3];
rz(-1.9726617) q[3];
sx q[3];
rz(-2.2973785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0323623) q[0];
sx q[0];
rz(-0.6210331) q[0];
sx q[0];
rz(-1.2048298) q[0];
rz(3.0989528) q[1];
sx q[1];
rz(-1.3748906) q[1];
sx q[1];
rz(-2.1991594) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.401817) q[0];
sx q[0];
rz(-1.4352006) q[0];
sx q[0];
rz(1.9996765) q[0];
rz(-pi) q[1];
rz(-2.3187014) q[2];
sx q[2];
rz(-2.3510755) q[2];
sx q[2];
rz(0.86871494) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1548295) q[1];
sx q[1];
rz(-1.0752605) q[1];
sx q[1];
rz(-2.7096728) q[1];
rz(0.14640866) q[3];
sx q[3];
rz(-2.4973828) q[3];
sx q[3];
rz(-0.25800426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.11613906) q[2];
sx q[2];
rz(-2.4567273) q[2];
sx q[2];
rz(-1.2241414) q[2];
rz(-2.4141566) q[3];
sx q[3];
rz(-1.4801315) q[3];
sx q[3];
rz(3.091403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8023476) q[0];
sx q[0];
rz(-0.50528637) q[0];
sx q[0];
rz(-0.25830609) q[0];
rz(-0.91336617) q[1];
sx q[1];
rz(-2.2759627) q[1];
sx q[1];
rz(-2.1736653) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25719698) q[0];
sx q[0];
rz(-1.9564187) q[0];
sx q[0];
rz(2.5325724) q[0];
rz(-pi) q[1];
rz(-1.9759693) q[2];
sx q[2];
rz(-2.3771667) q[2];
sx q[2];
rz(0.38328277) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.55019864) q[1];
sx q[1];
rz(-1.0934783) q[1];
sx q[1];
rz(-1.5309019) q[1];
x q[2];
rz(-2.1139891) q[3];
sx q[3];
rz(-2.1434727) q[3];
sx q[3];
rz(3.0610994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6496868) q[2];
sx q[2];
rz(-0.68444362) q[2];
sx q[2];
rz(-0.41927949) q[2];
rz(2.3573917) q[3];
sx q[3];
rz(-2.1214285) q[3];
sx q[3];
rz(1.7694582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.075031) q[0];
sx q[0];
rz(-3.0176268) q[0];
sx q[0];
rz(3.1228464) q[0];
rz(3.0491414) q[1];
sx q[1];
rz(-1.8094614) q[1];
sx q[1];
rz(-0.094955347) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2128747) q[0];
sx q[0];
rz(-1.2480019) q[0];
sx q[0];
rz(2.5514437) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9906391) q[2];
sx q[2];
rz(-2.125084) q[2];
sx q[2];
rz(1.8009763) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1996205) q[1];
sx q[1];
rz(-0.50040659) q[1];
sx q[1];
rz(-0.54897658) q[1];
rz(-2.9800013) q[3];
sx q[3];
rz(-1.1757903) q[3];
sx q[3];
rz(-0.99564161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2350601) q[2];
sx q[2];
rz(-1.1842714) q[2];
sx q[2];
rz(2.2825784) q[2];
rz(0.61339316) q[3];
sx q[3];
rz(-2.9412013) q[3];
sx q[3];
rz(1.8092417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61619806) q[0];
sx q[0];
rz(-1.2735676) q[0];
sx q[0];
rz(-1.4816351) q[0];
rz(-1.067465) q[1];
sx q[1];
rz(-0.72507247) q[1];
sx q[1];
rz(1.3190837) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7049371) q[0];
sx q[0];
rz(-2.4709513) q[0];
sx q[0];
rz(1.4192307) q[0];
rz(-pi) q[1];
rz(-2.5379792) q[2];
sx q[2];
rz(-2.5381662) q[2];
sx q[2];
rz(0.44925424) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1975648) q[1];
sx q[1];
rz(-2.135471) q[1];
sx q[1];
rz(0.27826635) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4433925) q[3];
sx q[3];
rz(-1.2166598) q[3];
sx q[3];
rz(0.69981258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8760425) q[2];
sx q[2];
rz(-2.7713113) q[2];
sx q[2];
rz(0.041157095) q[2];
rz(2.0187078) q[3];
sx q[3];
rz(-1.6582146) q[3];
sx q[3];
rz(-0.52367228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48162833) q[0];
sx q[0];
rz(-0.97844231) q[0];
sx q[0];
rz(1.6690669) q[0];
rz(0.72498471) q[1];
sx q[1];
rz(-1.1712733) q[1];
sx q[1];
rz(0.57458007) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0313505) q[0];
sx q[0];
rz(-1.6185805) q[0];
sx q[0];
rz(0.33804531) q[0];
rz(2.1442604) q[2];
sx q[2];
rz(-0.49066273) q[2];
sx q[2];
rz(0.034324797) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9159989) q[1];
sx q[1];
rz(-2.288842) q[1];
sx q[1];
rz(1.5313994) q[1];
x q[2];
rz(-0.88317849) q[3];
sx q[3];
rz(-1.6856582) q[3];
sx q[3];
rz(2.3695947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4943115) q[2];
sx q[2];
rz(-2.8936671) q[2];
sx q[2];
rz(-0.9606804) q[2];
rz(2.658355) q[3];
sx q[3];
rz(-1.5944696) q[3];
sx q[3];
rz(1.5620935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1305337) q[0];
sx q[0];
rz(-2.647825) q[0];
sx q[0];
rz(-0.44179398) q[0];
rz(-0.68253016) q[1];
sx q[1];
rz(-1.6923169) q[1];
sx q[1];
rz(-2.0880879) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.080314577) q[0];
sx q[0];
rz(-0.77064454) q[0];
sx q[0];
rz(1.944197) q[0];
rz(-pi) q[1];
rz(3.0526673) q[2];
sx q[2];
rz(-1.94238) q[2];
sx q[2];
rz(-3.0419797) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5805086) q[1];
sx q[1];
rz(-1.2496767) q[1];
sx q[1];
rz(0.91144423) q[1];
rz(2.9656762) q[3];
sx q[3];
rz(-2.7421167) q[3];
sx q[3];
rz(0.3829869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0992451) q[2];
sx q[2];
rz(-1.8727563) q[2];
sx q[2];
rz(-2.919) q[2];
rz(1.7706722) q[3];
sx q[3];
rz(-1.7226847) q[3];
sx q[3];
rz(0.62209779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-1.2400146) q[0];
sx q[0];
rz(-2.1299025) q[0];
sx q[0];
rz(1.0255751) q[0];
rz(1.1657794) q[1];
sx q[1];
rz(-0.60973254) q[1];
sx q[1];
rz(2.9235074) q[1];
rz(2.0596621) q[2];
sx q[2];
rz(-2.0578544) q[2];
sx q[2];
rz(-0.73357972) q[2];
rz(-0.00047666778) q[3];
sx q[3];
rz(-1.2491741) q[3];
sx q[3];
rz(0.76330976) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
