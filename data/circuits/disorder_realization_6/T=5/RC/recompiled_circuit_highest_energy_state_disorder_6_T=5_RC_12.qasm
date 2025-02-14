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
rz(-2.1036086) q[0];
sx q[0];
rz(-1.684364) q[0];
sx q[0];
rz(-1.7888223) q[0];
rz(1.1846722) q[1];
sx q[1];
rz(-0.67263043) q[1];
sx q[1];
rz(1.5157359) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0294581) q[0];
sx q[0];
rz(-1.7687651) q[0];
sx q[0];
rz(-2.9517118) q[0];
rz(-2.6625161) q[2];
sx q[2];
rz(-2.276439) q[2];
sx q[2];
rz(2.3912663) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.63187088) q[1];
sx q[1];
rz(-1.1051105) q[1];
sx q[1];
rz(-1.1759971) q[1];
x q[2];
rz(0.42077371) q[3];
sx q[3];
rz(-1.3075365) q[3];
sx q[3];
rz(-0.0076310633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1379913) q[2];
sx q[2];
rz(-3.0333952) q[2];
sx q[2];
rz(-2.9555964) q[2];
rz(3.1229535) q[3];
sx q[3];
rz(-1.2942856) q[3];
sx q[3];
rz(2.0712461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5828534) q[0];
sx q[0];
rz(-2.5753729) q[0];
sx q[0];
rz(0.3183611) q[0];
rz(-0.63703498) q[1];
sx q[1];
rz(-0.77992264) q[1];
sx q[1];
rz(2.6723518) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56457389) q[0];
sx q[0];
rz(-1.3718318) q[0];
sx q[0];
rz(1.3964723) q[0];
rz(-pi) q[1];
rz(-1.3606684) q[2];
sx q[2];
rz(-2.7581425) q[2];
sx q[2];
rz(2.0166778) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9904203) q[1];
sx q[1];
rz(-1.3440596) q[1];
sx q[1];
rz(2.8137035) q[1];
rz(-pi) q[2];
rz(1.5777228) q[3];
sx q[3];
rz(-2.4516461) q[3];
sx q[3];
rz(1.5010709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3176754) q[2];
sx q[2];
rz(-2.9313512) q[2];
sx q[2];
rz(0.69327411) q[2];
rz(-1.8215826) q[3];
sx q[3];
rz(-1.3258679) q[3];
sx q[3];
rz(1.0953974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1410809) q[0];
sx q[0];
rz(-2.5040099) q[0];
sx q[0];
rz(0.47955036) q[0];
rz(0.50411049) q[1];
sx q[1];
rz(-1.9934374) q[1];
sx q[1];
rz(-3.0832916) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8970015) q[0];
sx q[0];
rz(-2.069733) q[0];
sx q[0];
rz(1.2759491) q[0];
x q[1];
rz(-2.53251) q[2];
sx q[2];
rz(-0.78819617) q[2];
sx q[2];
rz(-0.64465514) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2268277) q[1];
sx q[1];
rz(-1.0221507) q[1];
sx q[1];
rz(0.1424205) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3065058) q[3];
sx q[3];
rz(-1.0796781) q[3];
sx q[3];
rz(2.6437505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7551859) q[2];
sx q[2];
rz(-1.9599954) q[2];
sx q[2];
rz(-0.03726658) q[2];
rz(-1.6218328) q[3];
sx q[3];
rz(-2.683679) q[3];
sx q[3];
rz(0.38145414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94654361) q[0];
sx q[0];
rz(-0.46849546) q[0];
sx q[0];
rz(-0.066548912) q[0];
rz(1.6171803) q[1];
sx q[1];
rz(-1.8981372) q[1];
sx q[1];
rz(2.4066511) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3108705) q[0];
sx q[0];
rz(-1.3390216) q[0];
sx q[0];
rz(-2.5481497) q[0];
rz(-pi) q[1];
x q[1];
rz(0.25937542) q[2];
sx q[2];
rz(-0.7470567) q[2];
sx q[2];
rz(-0.41058985) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0543136) q[1];
sx q[1];
rz(-0.4590408) q[1];
sx q[1];
rz(-2.4695818) q[1];
x q[2];
rz(-2.2913886) q[3];
sx q[3];
rz(-1.461004) q[3];
sx q[3];
rz(1.3444855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.28896004) q[2];
sx q[2];
rz(-1.5776878) q[2];
sx q[2];
rz(0.10870474) q[2];
rz(-2.2147801) q[3];
sx q[3];
rz(-2.1709397) q[3];
sx q[3];
rz(-1.577781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0573334) q[0];
sx q[0];
rz(-2.5045392) q[0];
sx q[0];
rz(1.0908701) q[0];
rz(0.57251656) q[1];
sx q[1];
rz(-1.9520091) q[1];
sx q[1];
rz(-0.82430878) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5474654) q[0];
sx q[0];
rz(-2.6931433) q[0];
sx q[0];
rz(-1.0065824) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9277497) q[2];
sx q[2];
rz(-1.130419) q[2];
sx q[2];
rz(1.7945031) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.2356154) q[1];
sx q[1];
rz(-0.69825497) q[1];
sx q[1];
rz(-2.2119616) q[1];
rz(-pi) q[2];
rz(-0.82309725) q[3];
sx q[3];
rz(-0.56328008) q[3];
sx q[3];
rz(2.7989504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6606286) q[2];
sx q[2];
rz(-2.169256) q[2];
sx q[2];
rz(0.29120293) q[2];
rz(2.5045942) q[3];
sx q[3];
rz(-1.4642508) q[3];
sx q[3];
rz(1.8891107) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3785192) q[0];
sx q[0];
rz(-2.851649) q[0];
sx q[0];
rz(-3.0105403) q[0];
rz(-1.1669) q[1];
sx q[1];
rz(-2.1578372) q[1];
sx q[1];
rz(-3.1189721) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1194099) q[0];
sx q[0];
rz(-1.9227322) q[0];
sx q[0];
rz(-0.58076136) q[0];
x q[1];
rz(2.6220967) q[2];
sx q[2];
rz(-0.24458376) q[2];
sx q[2];
rz(2.1313388) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.80161051) q[1];
sx q[1];
rz(-1.8793072) q[1];
sx q[1];
rz(-2.1478081) q[1];
rz(-pi) q[2];
x q[2];
rz(0.87225391) q[3];
sx q[3];
rz(-2.8917851) q[3];
sx q[3];
rz(-1.8810617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.268198) q[2];
sx q[2];
rz(-0.70316535) q[2];
sx q[2];
rz(-1.084729) q[2];
rz(-1.1449413) q[3];
sx q[3];
rz(-1.2866674) q[3];
sx q[3];
rz(1.0219319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5522083) q[0];
sx q[0];
rz(-1.5188058) q[0];
sx q[0];
rz(2.6680706) q[0];
rz(-1.2208968) q[1];
sx q[1];
rz(-2.7291606) q[1];
sx q[1];
rz(1.268505) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80618033) q[0];
sx q[0];
rz(-2.0783193) q[0];
sx q[0];
rz(2.1730459) q[0];
rz(-pi) q[1];
x q[1];
rz(0.73641554) q[2];
sx q[2];
rz(-2.066095) q[2];
sx q[2];
rz(-2.8335932) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4422677) q[1];
sx q[1];
rz(-2.0751157) q[1];
sx q[1];
rz(1.0902576) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0587993) q[3];
sx q[3];
rz(-1.1361381) q[3];
sx q[3];
rz(2.5744048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1641757) q[2];
sx q[2];
rz(-1.5722534) q[2];
sx q[2];
rz(0.24641985) q[2];
rz(0.94270802) q[3];
sx q[3];
rz(-2.0556512) q[3];
sx q[3];
rz(2.3781618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8134269) q[0];
sx q[0];
rz(-1.8712217) q[0];
sx q[0];
rz(-3.0810007) q[0];
rz(-1.6391485) q[1];
sx q[1];
rz(-1.7785347) q[1];
sx q[1];
rz(-1.5938119) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0133724) q[0];
sx q[0];
rz(-3.1376079) q[0];
sx q[0];
rz(-0.26040034) q[0];
rz(1.3826107) q[2];
sx q[2];
rz(-2.1379092) q[2];
sx q[2];
rz(-0.59394804) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2912078) q[1];
sx q[1];
rz(-2.3602848) q[1];
sx q[1];
rz(-0.92618295) q[1];
x q[2];
rz(2.0261062) q[3];
sx q[3];
rz(-2.2336965) q[3];
sx q[3];
rz(-2.2772706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8760406) q[2];
sx q[2];
rz(-0.92090845) q[2];
sx q[2];
rz(1.457224) q[2];
rz(3.0217116) q[3];
sx q[3];
rz(-0.20457743) q[3];
sx q[3];
rz(1.0286819) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84751713) q[0];
sx q[0];
rz(-0.99440044) q[0];
sx q[0];
rz(-1.1353528) q[0];
rz(3.0382233) q[1];
sx q[1];
rz(-1.4048978) q[1];
sx q[1];
rz(0.9476544) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5948828) q[0];
sx q[0];
rz(-0.93385392) q[0];
sx q[0];
rz(1.4875436) q[0];
x q[1];
rz(2.060076) q[2];
sx q[2];
rz(-2.3913417) q[2];
sx q[2];
rz(0.81847755) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2716213) q[1];
sx q[1];
rz(-2.4050887) q[1];
sx q[1];
rz(-1.9727712) q[1];
x q[2];
rz(-2.3941764) q[3];
sx q[3];
rz(-2.1028165) q[3];
sx q[3];
rz(2.4126787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.4244298) q[2];
sx q[2];
rz(-2.5134176) q[2];
sx q[2];
rz(-2.2372645) q[2];
rz(-1.8782015) q[3];
sx q[3];
rz(-1.9046013) q[3];
sx q[3];
rz(2.625107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.271027) q[0];
sx q[0];
rz(-2.3834383) q[0];
sx q[0];
rz(2.1562321) q[0];
rz(-2.789978) q[1];
sx q[1];
rz(-2.1814587) q[1];
sx q[1];
rz(0.34686372) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3540005) q[0];
sx q[0];
rz(-1.7890325) q[0];
sx q[0];
rz(2.2130593) q[0];
rz(-pi) q[1];
rz(-0.11254452) q[2];
sx q[2];
rz(-1.1273317) q[2];
sx q[2];
rz(0.36641589) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.736944) q[1];
sx q[1];
rz(-2.258806) q[1];
sx q[1];
rz(-1.4555803) q[1];
rz(-pi) q[2];
rz(0.44054359) q[3];
sx q[3];
rz(-1.9048573) q[3];
sx q[3];
rz(-1.3584709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2615307) q[2];
sx q[2];
rz(-2.4226649) q[2];
sx q[2];
rz(3.0066709) q[2];
rz(-3.0123582) q[3];
sx q[3];
rz(-2.7415469) q[3];
sx q[3];
rz(-2.1028886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1226817) q[0];
sx q[0];
rz(-1.2553348) q[0];
sx q[0];
rz(-3.0154764) q[0];
rz(0.46161721) q[1];
sx q[1];
rz(-1.4001662) q[1];
sx q[1];
rz(-2.9495159) q[1];
rz(2.2444637) q[2];
sx q[2];
rz(-1.2012325) q[2];
sx q[2];
rz(2.7129632) q[2];
rz(0.25855385) q[3];
sx q[3];
rz(-1.7333442) q[3];
sx q[3];
rz(2.6556591) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
