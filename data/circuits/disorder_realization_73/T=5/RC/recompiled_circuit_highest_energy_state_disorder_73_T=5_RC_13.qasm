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
rz(1.6950322) q[0];
sx q[0];
rz(8.5230081) q[0];
rz(1.3321441) q[1];
sx q[1];
rz(-2.6982215) q[1];
sx q[1];
rz(2.3763357) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9741824) q[0];
sx q[0];
rz(-1.2077792) q[0];
sx q[0];
rz(-2.1556751) q[0];
rz(-pi) q[1];
x q[1];
rz(0.95633219) q[2];
sx q[2];
rz(-0.40618375) q[2];
sx q[2];
rz(-1.150591) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.74435788) q[1];
sx q[1];
rz(-1.2325979) q[1];
sx q[1];
rz(-1.5282791) q[1];
rz(-pi) q[2];
rz(1.1422994) q[3];
sx q[3];
rz(-2.836578) q[3];
sx q[3];
rz(1.3950789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.29609933) q[2];
sx q[2];
rz(-1.0881492) q[2];
sx q[2];
rz(1.4720526) q[2];
rz(-1.7146401) q[3];
sx q[3];
rz(-1.497437) q[3];
sx q[3];
rz(-2.6078687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1832837) q[0];
sx q[0];
rz(-1.5077718) q[0];
sx q[0];
rz(-0.84877745) q[0];
rz(-3.1065885) q[1];
sx q[1];
rz(-1.1468381) q[1];
sx q[1];
rz(0.95357198) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7824421) q[0];
sx q[0];
rz(-1.1138565) q[0];
sx q[0];
rz(0.63553973) q[0];
rz(-pi) q[1];
rz(-2.8368901) q[2];
sx q[2];
rz(-1.1977473) q[2];
sx q[2];
rz(-2.3359131) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2767859) q[1];
sx q[1];
rz(-1.7662443) q[1];
sx q[1];
rz(1.2246183) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4783432) q[3];
sx q[3];
rz(-2.2683072) q[3];
sx q[3];
rz(2.8467941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1635052) q[2];
sx q[2];
rz(-1.7814025) q[2];
sx q[2];
rz(-0.84954849) q[2];
rz(0.16544011) q[3];
sx q[3];
rz(-1.4435507) q[3];
sx q[3];
rz(-0.74976841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.36990944) q[0];
sx q[0];
rz(-0.92215466) q[0];
sx q[0];
rz(-0.40192303) q[0];
rz(0.15332128) q[1];
sx q[1];
rz(-0.17686495) q[1];
sx q[1];
rz(-1.1361928) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76321852) q[0];
sx q[0];
rz(-2.741886) q[0];
sx q[0];
rz(1.5399152) q[0];
rz(-0.49164518) q[2];
sx q[2];
rz(-2.6043713) q[2];
sx q[2];
rz(1.1342837) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.043667533) q[1];
sx q[1];
rz(-2.8043724) q[1];
sx q[1];
rz(-0.80101669) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6776776) q[3];
sx q[3];
rz(-1.216421) q[3];
sx q[3];
rz(-2.7257435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.084126964) q[2];
sx q[2];
rz(-0.87856138) q[2];
sx q[2];
rz(-3.086536) q[2];
rz(-1.7940686) q[3];
sx q[3];
rz(-1.3301347) q[3];
sx q[3];
rz(-1.0042892) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7771626) q[0];
sx q[0];
rz(-0.20474064) q[0];
sx q[0];
rz(-1.5675911) q[0];
rz(-0.69152999) q[1];
sx q[1];
rz(-2.4722996) q[1];
sx q[1];
rz(0.18542586) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2381993) q[0];
sx q[0];
rz(-1.1739587) q[0];
sx q[0];
rz(1.7522041) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.71132423) q[2];
sx q[2];
rz(-2.5344116) q[2];
sx q[2];
rz(1.9679697) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1048829) q[1];
sx q[1];
rz(-1.2598902) q[1];
sx q[1];
rz(0.26581146) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.809285) q[3];
sx q[3];
rz(-0.57113591) q[3];
sx q[3];
rz(-2.6337998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7376248) q[2];
sx q[2];
rz(-1.7125968) q[2];
sx q[2];
rz(-0.0052304012) q[2];
rz(-0.38749203) q[3];
sx q[3];
rz(-0.5046851) q[3];
sx q[3];
rz(1.3030049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4150647) q[0];
sx q[0];
rz(-2.1580577) q[0];
sx q[0];
rz(0.61597419) q[0];
rz(1.1888602) q[1];
sx q[1];
rz(-0.6183466) q[1];
sx q[1];
rz(-0.51330769) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6270638) q[0];
sx q[0];
rz(-1.5781869) q[0];
sx q[0];
rz(1.2891931) q[0];
rz(1.8525858) q[2];
sx q[2];
rz(-2.2766487) q[2];
sx q[2];
rz(-2.6633546) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9037902) q[1];
sx q[1];
rz(-2.2040743) q[1];
sx q[1];
rz(-1.6602181) q[1];
rz(-1.1034455) q[3];
sx q[3];
rz(-1.6765521) q[3];
sx q[3];
rz(0.33973628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1335699) q[2];
sx q[2];
rz(-1.8654537) q[2];
sx q[2];
rz(2.3991154) q[2];
rz(1.4437458) q[3];
sx q[3];
rz(-1.7323114) q[3];
sx q[3];
rz(1.1013364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83108574) q[0];
sx q[0];
rz(-2.0948912) q[0];
sx q[0];
rz(2.781784) q[0];
rz(0.85044914) q[1];
sx q[1];
rz(-1.9603739) q[1];
sx q[1];
rz(-2.6757619) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82079262) q[0];
sx q[0];
rz(-2.3893836) q[0];
sx q[0];
rz(1.8899931) q[0];
rz(-2.3783422) q[2];
sx q[2];
rz(-1.6692729) q[2];
sx q[2];
rz(-1.5710448) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9656495) q[1];
sx q[1];
rz(-1.5197943) q[1];
sx q[1];
rz(-2.9005364) q[1];
rz(-2.1773866) q[3];
sx q[3];
rz(-2.0584985) q[3];
sx q[3];
rz(-1.4473947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.74238527) q[2];
sx q[2];
rz(-1.1804322) q[2];
sx q[2];
rz(0.36435374) q[2];
rz(-2.897701) q[3];
sx q[3];
rz(-3.0920005) q[3];
sx q[3];
rz(1.9243139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57109443) q[0];
sx q[0];
rz(-2.1622393) q[0];
sx q[0];
rz(-1.1660227) q[0];
rz(-2.1265538) q[1];
sx q[1];
rz(-2.6045585) q[1];
sx q[1];
rz(-0.38645116) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5597661) q[0];
sx q[0];
rz(-0.21187267) q[0];
sx q[0];
rz(-1.7170402) q[0];
x q[1];
rz(-2.831424) q[2];
sx q[2];
rz(-1.8156689) q[2];
sx q[2];
rz(2.819811) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9484752) q[1];
sx q[1];
rz(-1.5453234) q[1];
sx q[1];
rz(-1.6465882) q[1];
rz(-1.0284958) q[3];
sx q[3];
rz(-1.4631548) q[3];
sx q[3];
rz(0.018317761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7361136) q[2];
sx q[2];
rz(-1.6454641) q[2];
sx q[2];
rz(-2.9138937) q[2];
rz(-2.0685711) q[3];
sx q[3];
rz(-2.3927972) q[3];
sx q[3];
rz(-2.8480215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88184083) q[0];
sx q[0];
rz(-2.2417534) q[0];
sx q[0];
rz(0.42801273) q[0];
rz(1.0182861) q[1];
sx q[1];
rz(-0.4363474) q[1];
sx q[1];
rz(2.2705618) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6414531) q[0];
sx q[0];
rz(-0.79384168) q[0];
sx q[0];
rz(-0.83896933) q[0];
x q[1];
rz(2.0823048) q[2];
sx q[2];
rz(-2.7370484) q[2];
sx q[2];
rz(1.2407357) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4531252) q[1];
sx q[1];
rz(-0.18632132) q[1];
sx q[1];
rz(-0.86155714) q[1];
x q[2];
rz(1.0139129) q[3];
sx q[3];
rz(-0.96864163) q[3];
sx q[3];
rz(0.35246655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.13228664) q[2];
sx q[2];
rz(-0.79557482) q[2];
sx q[2];
rz(-1.2786678) q[2];
rz(-0.57847413) q[3];
sx q[3];
rz(-2.5090802) q[3];
sx q[3];
rz(-1.1540029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6744252) q[0];
sx q[0];
rz(-0.217087) q[0];
sx q[0];
rz(0.59980741) q[0];
rz(-1.8425875) q[1];
sx q[1];
rz(-1.6103585) q[1];
sx q[1];
rz(-0.64186796) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1917832) q[0];
sx q[0];
rz(-1.2300175) q[0];
sx q[0];
rz(-2.6343915) q[0];
x q[1];
rz(-1.4207877) q[2];
sx q[2];
rz(-2.370655) q[2];
sx q[2];
rz(-1.9341759) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.95002) q[1];
sx q[1];
rz(-0.86314647) q[1];
sx q[1];
rz(1.5465082) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0976552) q[3];
sx q[3];
rz(-0.88636413) q[3];
sx q[3];
rz(1.8620644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.82015264) q[2];
sx q[2];
rz(-1.7165311) q[2];
sx q[2];
rz(0.610262) q[2];
rz(-0.52014703) q[3];
sx q[3];
rz(-1.0545694) q[3];
sx q[3];
rz(-2.648073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.801012) q[0];
sx q[0];
rz(-1.2435253) q[0];
sx q[0];
rz(-0.93801671) q[0];
rz(-0.34171379) q[1];
sx q[1];
rz(-1.7594756) q[1];
sx q[1];
rz(-0.15728532) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92625916) q[0];
sx q[0];
rz(-0.65500998) q[0];
sx q[0];
rz(0.32230719) q[0];
x q[1];
rz(-0.089265169) q[2];
sx q[2];
rz(-2.3815739) q[2];
sx q[2];
rz(1.6940728) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.54522911) q[1];
sx q[1];
rz(-2.5198054) q[1];
sx q[1];
rz(-0.53485628) q[1];
x q[2];
rz(1.0299133) q[3];
sx q[3];
rz(-2.6737041) q[3];
sx q[3];
rz(-0.44178007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9226795) q[2];
sx q[2];
rz(-2.7965386) q[2];
sx q[2];
rz(1.2738796) q[2];
rz(3.1359361) q[3];
sx q[3];
rz(-1.4414682) q[3];
sx q[3];
rz(-0.98102942) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
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
rz(2.2412655) q[2];
sx q[2];
rz(-1.1569958) q[2];
sx q[2];
rz(-0.3668084) q[2];
rz(-1.6531108) q[3];
sx q[3];
rz(-2.001279) q[3];
sx q[3];
rz(1.6829987) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
