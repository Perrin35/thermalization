OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72803175) q[0];
sx q[0];
rz(-0.95946884) q[0];
sx q[0];
rz(2.5897107) q[0];
rz(-0.38504398) q[1];
sx q[1];
rz(-1.3621962) q[1];
sx q[1];
rz(0.41866067) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.066025) q[0];
sx q[0];
rz(-0.93187823) q[0];
sx q[0];
rz(-2.781771) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35595591) q[2];
sx q[2];
rz(-0.47544962) q[2];
sx q[2];
rz(-2.8048777) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6578411) q[1];
sx q[1];
rz(-2.9194909) q[1];
sx q[1];
rz(-2.7276843) q[1];
rz(-pi) q[2];
rz(1.5403662) q[3];
sx q[3];
rz(-1.55369) q[3];
sx q[3];
rz(-2.3612983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.60575134) q[2];
sx q[2];
rz(-1.8391106) q[2];
sx q[2];
rz(2.7525986) q[2];
rz(2.244106) q[3];
sx q[3];
rz(-0.56098452) q[3];
sx q[3];
rz(-1.2787904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9101343) q[0];
sx q[0];
rz(-2.1899905) q[0];
sx q[0];
rz(2.8125473) q[0];
rz(-1.344205) q[1];
sx q[1];
rz(-0.75733328) q[1];
sx q[1];
rz(0.12408852) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50590512) q[0];
sx q[0];
rz(-1.3342739) q[0];
sx q[0];
rz(-0.2157477) q[0];
rz(1.7089073) q[2];
sx q[2];
rz(-1.8542854) q[2];
sx q[2];
rz(-2.5140009) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.29247788) q[1];
sx q[1];
rz(-2.5579331) q[1];
sx q[1];
rz(2.2569342) q[1];
rz(-pi) q[2];
rz(-2.9355818) q[3];
sx q[3];
rz(-0.877218) q[3];
sx q[3];
rz(1.3242974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.2825534) q[2];
sx q[2];
rz(-2.6586847) q[2];
sx q[2];
rz(0.60749751) q[2];
rz(2.2533158) q[3];
sx q[3];
rz(-1.5253303) q[3];
sx q[3];
rz(1.7150778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68179503) q[0];
sx q[0];
rz(-1.8495704) q[0];
sx q[0];
rz(-2.9070396) q[0];
rz(1.0598496) q[1];
sx q[1];
rz(-1.1666965) q[1];
sx q[1];
rz(0.92811981) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2560349) q[0];
sx q[0];
rz(-1.9826898) q[0];
sx q[0];
rz(-2.4649863) q[0];
rz(-pi) q[1];
rz(-2.6675148) q[2];
sx q[2];
rz(-1.9454207) q[2];
sx q[2];
rz(1.4371536) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.249284) q[1];
sx q[1];
rz(-1.705348) q[1];
sx q[1];
rz(2.8110912) q[1];
rz(-1.9409402) q[3];
sx q[3];
rz(-1.7142152) q[3];
sx q[3];
rz(-2.7756144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5782535) q[2];
sx q[2];
rz(-1.4207062) q[2];
sx q[2];
rz(1.2467747) q[2];
rz(-2.9747544) q[3];
sx q[3];
rz(-0.92883795) q[3];
sx q[3];
rz(1.5290574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.700915) q[0];
sx q[0];
rz(-3.0755141) q[0];
sx q[0];
rz(-0.75827688) q[0];
rz(-1.5757163) q[1];
sx q[1];
rz(-2.5570452) q[1];
sx q[1];
rz(0.87055269) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93773952) q[0];
sx q[0];
rz(-2.014262) q[0];
sx q[0];
rz(2.4340043) q[0];
rz(1.0494558) q[2];
sx q[2];
rz(-2.0851729) q[2];
sx q[2];
rz(-0.55869049) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.16337559) q[1];
sx q[1];
rz(-1.4184457) q[1];
sx q[1];
rz(1.1635029) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0421889) q[3];
sx q[3];
rz(-2.5113341) q[3];
sx q[3];
rz(-0.46168345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7190651) q[2];
sx q[2];
rz(-2.1336522) q[2];
sx q[2];
rz(2.3818805) q[2];
rz(0.64940137) q[3];
sx q[3];
rz(-1.6067182) q[3];
sx q[3];
rz(-1.5627129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11288697) q[0];
sx q[0];
rz(-0.6539456) q[0];
sx q[0];
rz(1.1829859) q[0];
rz(0.68823367) q[1];
sx q[1];
rz(-1.3166683) q[1];
sx q[1];
rz(0.36929718) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6128176) q[0];
sx q[0];
rz(-1.5433831) q[0];
sx q[0];
rz(2.4779336) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4299116) q[2];
sx q[2];
rz(-1.3698915) q[2];
sx q[2];
rz(2.7679659) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6162045) q[1];
sx q[1];
rz(-1.7553575) q[1];
sx q[1];
rz(-2.2730458) q[1];
rz(1.9019674) q[3];
sx q[3];
rz(-1.2544362) q[3];
sx q[3];
rz(1.4671732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1987622) q[2];
sx q[2];
rz(-1.3254712) q[2];
sx q[2];
rz(-2.14373) q[2];
rz(-0.61740795) q[3];
sx q[3];
rz(-0.95881763) q[3];
sx q[3];
rz(0.82120419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.038641039) q[0];
sx q[0];
rz(-1.8674253) q[0];
sx q[0];
rz(1.8983023) q[0];
rz(-2.359911) q[1];
sx q[1];
rz(-1.2739173) q[1];
sx q[1];
rz(0.26085687) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8600966) q[0];
sx q[0];
rz(-1.3602169) q[0];
sx q[0];
rz(2.6407535) q[0];
rz(-pi) q[1];
rz(2.2736808) q[2];
sx q[2];
rz(-0.54722584) q[2];
sx q[2];
rz(0.99582129) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1697537) q[1];
sx q[1];
rz(-1.1128281) q[1];
sx q[1];
rz(0.16916807) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9294037) q[3];
sx q[3];
rz(-2.1916323) q[3];
sx q[3];
rz(0.60597907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6692052) q[2];
sx q[2];
rz(-2.278625) q[2];
sx q[2];
rz(-0.48259398) q[2];
rz(-0.11416301) q[3];
sx q[3];
rz(-2.0518905) q[3];
sx q[3];
rz(2.193006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77940762) q[0];
sx q[0];
rz(-0.055883378) q[0];
sx q[0];
rz(2.8756397) q[0];
rz(-0.42568046) q[1];
sx q[1];
rz(-1.8961779) q[1];
sx q[1];
rz(-0.46806213) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45586203) q[0];
sx q[0];
rz(-2.1121368) q[0];
sx q[0];
rz(-0.0079519072) q[0];
x q[1];
rz(1.3139901) q[2];
sx q[2];
rz(-0.35126424) q[2];
sx q[2];
rz(-2.6333269) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.148273) q[1];
sx q[1];
rz(-0.56374218) q[1];
sx q[1];
rz(1.0757273) q[1];
rz(-3.1105177) q[3];
sx q[3];
rz(-1.0508176) q[3];
sx q[3];
rz(-1.7620373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.88428664) q[2];
sx q[2];
rz(-1.2904737) q[2];
sx q[2];
rz(-2.5878944) q[2];
rz(-2.2181559) q[3];
sx q[3];
rz(-2.2348576) q[3];
sx q[3];
rz(-3.0807909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(1.6786574) q[0];
sx q[0];
rz(-2.904992) q[0];
sx q[0];
rz(2.4635354) q[0];
rz(-2.1642115) q[1];
sx q[1];
rz(-1.9934318) q[1];
sx q[1];
rz(1.703702) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17409731) q[0];
sx q[0];
rz(-0.50104841) q[0];
sx q[0];
rz(1.1606085) q[0];
rz(-pi) q[1];
rz(-0.0043785574) q[2];
sx q[2];
rz(-2.5776049) q[2];
sx q[2];
rz(-2.0019238) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6594636) q[1];
sx q[1];
rz(-1.8244484) q[1];
sx q[1];
rz(-1.3319904) q[1];
rz(0.74681654) q[3];
sx q[3];
rz(-1.9821315) q[3];
sx q[3];
rz(-0.29838994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4108654) q[2];
sx q[2];
rz(-0.56034708) q[2];
sx q[2];
rz(0.66696683) q[2];
rz(0.11735958) q[3];
sx q[3];
rz(-1.2753692) q[3];
sx q[3];
rz(-1.7608775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9475107) q[0];
sx q[0];
rz(-1.5912594) q[0];
sx q[0];
rz(0.34341735) q[0];
rz(1.4944448) q[1];
sx q[1];
rz(-2.1518555) q[1];
sx q[1];
rz(2.6572878) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0065932) q[0];
sx q[0];
rz(-1.3065225) q[0];
sx q[0];
rz(-2.5359383) q[0];
rz(-pi) q[1];
rz(1.1555919) q[2];
sx q[2];
rz(-1.1744497) q[2];
sx q[2];
rz(1.9981148) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3266294) q[1];
sx q[1];
rz(-0.99363929) q[1];
sx q[1];
rz(-1.3054791) q[1];
x q[2];
rz(-0.72192854) q[3];
sx q[3];
rz(-2.8145104) q[3];
sx q[3];
rz(-0.41072911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1962428) q[2];
sx q[2];
rz(-2.0413155) q[2];
sx q[2];
rz(-0.84570447) q[2];
rz(-1.2196994) q[3];
sx q[3];
rz(-1.2518576) q[3];
sx q[3];
rz(-2.9367101) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6824816) q[0];
sx q[0];
rz(-2.8963608) q[0];
sx q[0];
rz(-2.5526175) q[0];
rz(-0.67539769) q[1];
sx q[1];
rz(-0.94869906) q[1];
sx q[1];
rz(1.677547) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9927514) q[0];
sx q[0];
rz(-2.0939079) q[0];
sx q[0];
rz(2.6639725) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4824994) q[2];
sx q[2];
rz(-0.73511926) q[2];
sx q[2];
rz(-2.5500848) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3879536) q[1];
sx q[1];
rz(-1.6170701) q[1];
sx q[1];
rz(0.58677499) q[1];
rz(-2.0328224) q[3];
sx q[3];
rz(-2.2115876) q[3];
sx q[3];
rz(0.45087157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7318763) q[2];
sx q[2];
rz(-1.7475374) q[2];
sx q[2];
rz(-1.124292) q[2];
rz(2.2935947) q[3];
sx q[3];
rz(-2.336899) q[3];
sx q[3];
rz(3.1112352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.547434) q[0];
sx q[0];
rz(-1.1347329) q[0];
sx q[0];
rz(-1.5190079) q[0];
rz(-0.92319725) q[1];
sx q[1];
rz(-2.1122439) q[1];
sx q[1];
rz(-1.5069638) q[1];
rz(-1.7566924) q[2];
sx q[2];
rz(-1.2658613) q[2];
sx q[2];
rz(-1.801898) q[2];
rz(-0.93529978) q[3];
sx q[3];
rz(-0.37105303) q[3];
sx q[3];
rz(1.1340352) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
