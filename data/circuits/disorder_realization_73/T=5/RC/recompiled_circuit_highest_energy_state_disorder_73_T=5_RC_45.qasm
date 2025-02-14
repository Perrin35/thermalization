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
rz(-2.7849164) q[0];
sx q[0];
rz(-1.4465605) q[0];
sx q[0];
rz(-2.2398228) q[0];
rz(-1.8094485) q[1];
sx q[1];
rz(-0.44337115) q[1];
sx q[1];
rz(-2.3763357) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0527406) q[0];
sx q[0];
rz(-0.67702261) q[0];
sx q[0];
rz(0.96816008) q[0];
rz(-pi) q[1];
x q[1];
rz(0.24306007) q[2];
sx q[2];
rz(-1.2420734) q[2];
sx q[2];
rz(-0.49561938) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2697077) q[1];
sx q[1];
rz(-0.34075865) q[1];
sx q[1];
rz(0.12026708) q[1];
rz(-pi) q[2];
x q[2];
rz(0.13007836) q[3];
sx q[3];
rz(-1.8474694) q[3];
sx q[3];
rz(0.94844669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8454933) q[2];
sx q[2];
rz(-1.0881492) q[2];
sx q[2];
rz(-1.4720526) q[2];
rz(1.4269525) q[3];
sx q[3];
rz(-1.497437) q[3];
sx q[3];
rz(-2.6078687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1832837) q[0];
sx q[0];
rz(-1.6338209) q[0];
sx q[0];
rz(-2.2928152) q[0];
rz(-3.1065885) q[1];
sx q[1];
rz(-1.1468381) q[1];
sx q[1];
rz(0.95357198) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6153107) q[0];
sx q[0];
rz(-1.0089951) q[0];
sx q[0];
rz(-1.022382) q[0];
rz(-pi) q[1];
rz(1.9601279) q[2];
sx q[2];
rz(-1.2876533) q[2];
sx q[2];
rz(-2.2623537) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.86480675) q[1];
sx q[1];
rz(-1.3753483) q[1];
sx q[1];
rz(-1.9169743) q[1];
rz(0.75462975) q[3];
sx q[3];
rz(-2.0622232) q[3];
sx q[3];
rz(-1.400465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1635052) q[2];
sx q[2];
rz(-1.3601902) q[2];
sx q[2];
rz(0.84954849) q[2];
rz(-2.9761525) q[3];
sx q[3];
rz(-1.4435507) q[3];
sx q[3];
rz(-0.74976841) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36990944) q[0];
sx q[0];
rz(-0.92215466) q[0];
sx q[0];
rz(0.40192303) q[0];
rz(2.9882714) q[1];
sx q[1];
rz(-2.9647277) q[1];
sx q[1];
rz(-1.1361928) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4118958) q[0];
sx q[0];
rz(-1.9703016) q[0];
sx q[0];
rz(3.1285498) q[0];
rz(-pi) q[1];
rz(1.8449144) q[2];
sx q[2];
rz(-1.1027567) q[2];
sx q[2];
rz(-0.57690358) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3553149) q[1];
sx q[1];
rz(-1.8031562) q[1];
sx q[1];
rz(1.3241598) q[1];
x q[2];
rz(-0.69092423) q[3];
sx q[3];
rz(-0.57583416) q[3];
sx q[3];
rz(-1.761328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.084126964) q[2];
sx q[2];
rz(-0.87856138) q[2];
sx q[2];
rz(-0.055056661) q[2];
rz(1.347524) q[3];
sx q[3];
rz(-1.811458) q[3];
sx q[3];
rz(1.0042892) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36443001) q[0];
sx q[0];
rz(-2.936852) q[0];
sx q[0];
rz(1.5740016) q[0];
rz(0.69152999) q[1];
sx q[1];
rz(-0.66929308) q[1];
sx q[1];
rz(-2.9561668) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5449654) q[0];
sx q[0];
rz(-1.7379679) q[0];
sx q[0];
rz(-2.7388229) q[0];
x q[1];
rz(-2.4302684) q[2];
sx q[2];
rz(-2.5344116) q[2];
sx q[2];
rz(-1.9679697) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.37764964) q[1];
sx q[1];
rz(-2.7353706) q[1];
sx q[1];
rz(-2.2561314) q[1];
rz(-0.54585643) q[3];
sx q[3];
rz(-1.7480769) q[3];
sx q[3];
rz(-1.3455678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4039679) q[2];
sx q[2];
rz(-1.7125968) q[2];
sx q[2];
rz(3.1363623) q[2];
rz(-0.38749203) q[3];
sx q[3];
rz(-2.6369075) q[3];
sx q[3];
rz(-1.3030049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72652793) q[0];
sx q[0];
rz(-0.98353493) q[0];
sx q[0];
rz(0.61597419) q[0];
rz(-1.9527324) q[1];
sx q[1];
rz(-2.523246) q[1];
sx q[1];
rz(0.51330769) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5145288) q[0];
sx q[0];
rz(-1.5634057) q[0];
sx q[0];
rz(-1.8523995) q[0];
rz(-pi) q[1];
rz(1.2890069) q[2];
sx q[2];
rz(-2.2766487) q[2];
sx q[2];
rz(-0.47823804) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.087438093) q[1];
sx q[1];
rz(-2.5028924) q[1];
sx q[1];
rz(0.12104669) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0232459) q[3];
sx q[3];
rz(-2.0353298) q[3];
sx q[3];
rz(-1.9637513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0080228) q[2];
sx q[2];
rz(-1.2761389) q[2];
sx q[2];
rz(-0.7424773) q[2];
rz(1.6978469) q[3];
sx q[3];
rz(-1.7323114) q[3];
sx q[3];
rz(2.0402563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83108574) q[0];
sx q[0];
rz(-2.0948912) q[0];
sx q[0];
rz(0.35980862) q[0];
rz(-0.85044914) q[1];
sx q[1];
rz(-1.1812187) q[1];
sx q[1];
rz(0.46583072) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98679989) q[0];
sx q[0];
rz(-1.354711) q[0];
sx q[0];
rz(-0.84439069) q[0];
rz(0.76325046) q[2];
sx q[2];
rz(-1.6692729) q[2];
sx q[2];
rz(-1.5710448) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9512144) q[1];
sx q[1];
rz(-2.8953027) q[1];
sx q[1];
rz(-0.21065335) q[1];
rz(-pi) q[2];
rz(-0.57329081) q[3];
sx q[3];
rz(-2.0985641) q[3];
sx q[3];
rz(-2.9506369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.74238527) q[2];
sx q[2];
rz(-1.1804322) q[2];
sx q[2];
rz(2.7772389) q[2];
rz(2.897701) q[3];
sx q[3];
rz(-0.049592169) q[3];
sx q[3];
rz(-1.2172788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57109443) q[0];
sx q[0];
rz(-0.97935337) q[0];
sx q[0];
rz(1.1660227) q[0];
rz(-1.0150389) q[1];
sx q[1];
rz(-2.6045585) q[1];
sx q[1];
rz(0.38645116) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58182654) q[0];
sx q[0];
rz(-2.92972) q[0];
sx q[0];
rz(-1.4245524) q[0];
rz(-2.4555444) q[2];
sx q[2];
rz(-0.39271388) q[2];
sx q[2];
rz(-0.60143747) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1931175) q[1];
sx q[1];
rz(-1.5962692) q[1];
sx q[1];
rz(-1.6465882) q[1];
x q[2];
rz(3.0160955) q[3];
sx q[3];
rz(-1.0319796) q[3];
sx q[3];
rz(1.5244689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.40547907) q[2];
sx q[2];
rz(-1.6454641) q[2];
sx q[2];
rz(-2.9138937) q[2];
rz(1.0730216) q[3];
sx q[3];
rz(-0.74879542) q[3];
sx q[3];
rz(-0.29357114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88184083) q[0];
sx q[0];
rz(-2.2417534) q[0];
sx q[0];
rz(-2.7135799) q[0];
rz(1.0182861) q[1];
sx q[1];
rz(-2.7052453) q[1];
sx q[1];
rz(-2.2705618) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7334325) q[0];
sx q[0];
rz(-2.1299612) q[0];
sx q[0];
rz(-0.59691043) q[0];
rz(-1.0592878) q[2];
sx q[2];
rz(-2.7370484) q[2];
sx q[2];
rz(-1.900857) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6884675) q[1];
sx q[1];
rz(-0.18632132) q[1];
sx q[1];
rz(2.2800355) q[1];
rz(-0.68058859) q[3];
sx q[3];
rz(-1.1201123) q[3];
sx q[3];
rz(-2.2623073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.13228664) q[2];
sx q[2];
rz(-0.79557482) q[2];
sx q[2];
rz(1.8629249) q[2];
rz(-2.5631185) q[3];
sx q[3];
rz(-0.63251248) q[3];
sx q[3];
rz(1.9875897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6744252) q[0];
sx q[0];
rz(-2.9245057) q[0];
sx q[0];
rz(2.5417852) q[0];
rz(1.8425875) q[1];
sx q[1];
rz(-1.6103585) q[1];
sx q[1];
rz(-2.4997247) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94980948) q[0];
sx q[0];
rz(-1.2300175) q[0];
sx q[0];
rz(0.50720117) q[0];
rz(-pi) q[1];
x q[1];
rz(0.14417837) q[2];
sx q[2];
rz(-0.81071172) q[2];
sx q[2];
rz(-1.7264896) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.95002) q[1];
sx q[1];
rz(-2.2784462) q[1];
sx q[1];
rz(1.5465082) q[1];
rz(-2.3849726) q[3];
sx q[3];
rz(-1.1706268) q[3];
sx q[3];
rz(-2.4979046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.82015264) q[2];
sx q[2];
rz(-1.7165311) q[2];
sx q[2];
rz(-0.610262) q[2];
rz(-0.52014703) q[3];
sx q[3];
rz(-2.0870233) q[3];
sx q[3];
rz(-0.49351969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.34058061) q[0];
sx q[0];
rz(-1.2435253) q[0];
sx q[0];
rz(2.2035759) q[0];
rz(-2.7998789) q[1];
sx q[1];
rz(-1.7594756) q[1];
sx q[1];
rz(0.15728532) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7559497) q[0];
sx q[0];
rz(-1.3766212) q[0];
sx q[0];
rz(2.5119378) q[0];
rz(-pi) q[1];
rz(1.4862655) q[2];
sx q[2];
rz(-0.8145552) q[2];
sx q[2];
rz(1.5712122) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5648236) q[1];
sx q[1];
rz(-1.2693468) q[1];
sx q[1];
rz(0.5524854) q[1];
rz(-1.1620164) q[3];
sx q[3];
rz(-1.805154) q[3];
sx q[3];
rz(-0.63692585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.21891317) q[2];
sx q[2];
rz(-2.7965386) q[2];
sx q[2];
rz(-1.867713) q[2];
rz(0.0056565469) q[3];
sx q[3];
rz(-1.7001245) q[3];
sx q[3];
rz(2.1605632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2904749) q[0];
sx q[0];
rz(-1.9708451) q[0];
sx q[0];
rz(-3.0921902) q[0];
rz(0.46650096) q[1];
sx q[1];
rz(-1.3460881) q[1];
sx q[1];
rz(0.64429611) q[1];
rz(-0.51085761) q[2];
sx q[2];
rz(-0.96561531) q[2];
sx q[2];
rz(-2.2462788) q[2];
rz(-0.43177035) q[3];
sx q[3];
rz(-1.4960066) q[3];
sx q[3];
rz(-3.0638051) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
