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
rz(2.1821238) q[0];
sx q[0];
rz(9.9766599) q[0];
rz(2.7565487) q[1];
sx q[1];
rz(-1.7793964) q[1];
sx q[1];
rz(-0.41866067) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4872348) q[0];
sx q[0];
rz(-2.4208768) q[0];
sx q[0];
rz(2.0134175) q[0];
rz(-pi) q[1];
rz(1.39327) q[2];
sx q[2];
rz(-1.1273618) q[2];
sx q[2];
rz(2.4088032) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4837515) q[1];
sx q[1];
rz(-0.22210177) q[1];
sx q[1];
rz(-2.7276843) q[1];
rz(1.0585634) q[3];
sx q[3];
rz(-0.034907428) q[3];
sx q[3];
rz(1.839118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5358413) q[2];
sx q[2];
rz(-1.8391106) q[2];
sx q[2];
rz(-2.7525986) q[2];
rz(2.244106) q[3];
sx q[3];
rz(-0.56098452) q[3];
sx q[3];
rz(1.8628023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9101343) q[0];
sx q[0];
rz(-0.95160216) q[0];
sx q[0];
rz(-0.32904539) q[0];
rz(-1.7973876) q[1];
sx q[1];
rz(-2.3842594) q[1];
sx q[1];
rz(0.12408852) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6356875) q[0];
sx q[0];
rz(-1.8073188) q[0];
sx q[0];
rz(-2.9258449) q[0];
x q[1];
rz(-1.7089073) q[2];
sx q[2];
rz(-1.2873073) q[2];
sx q[2];
rz(0.62759179) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8777496) q[1];
sx q[1];
rz(-1.9274492) q[1];
sx q[1];
rz(-2.043173) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.8666919) q[3];
sx q[3];
rz(-1.7287489) q[3];
sx q[3];
rz(-0.37930909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8590392) q[2];
sx q[2];
rz(-0.48290792) q[2];
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
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68179503) q[0];
sx q[0];
rz(-1.2920222) q[0];
sx q[0];
rz(2.9070396) q[0];
rz(1.0598496) q[1];
sx q[1];
rz(-1.1666965) q[1];
sx q[1];
rz(-2.2134728) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88555777) q[0];
sx q[0];
rz(-1.9826898) q[0];
sx q[0];
rz(2.4649863) q[0];
x q[1];
rz(-2.4305827) q[2];
sx q[2];
rz(-2.5463422) q[2];
sx q[2];
rz(2.3884515) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8923087) q[1];
sx q[1];
rz(-1.705348) q[1];
sx q[1];
rz(0.3305015) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9506128) q[3];
sx q[3];
rz(-0.39576021) q[3];
sx q[3];
rz(0.8518962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5633391) q[2];
sx q[2];
rz(-1.7208865) q[2];
sx q[2];
rz(-1.2467747) q[2];
rz(2.9747544) q[3];
sx q[3];
rz(-2.2127547) q[3];
sx q[3];
rz(-1.6125352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4406776) q[0];
sx q[0];
rz(-0.066078521) q[0];
sx q[0];
rz(-0.75827688) q[0];
rz(1.5757163) q[1];
sx q[1];
rz(-0.58454746) q[1];
sx q[1];
rz(0.87055269) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1567932) q[0];
sx q[0];
rz(-2.198303) q[0];
sx q[0];
rz(2.1294562) q[0];
rz(-pi) q[1];
rz(-1.0494558) q[2];
sx q[2];
rz(-2.0851729) q[2];
sx q[2];
rz(-2.5829022) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0724277) q[1];
sx q[1];
rz(-0.43335763) q[1];
sx q[1];
rz(-1.2010203) q[1];
rz(-pi) q[2];
x q[2];
rz(1.643067) q[3];
sx q[3];
rz(-0.94413432) q[3];
sx q[3];
rz(0.58451239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7190651) q[2];
sx q[2];
rz(-1.0079404) q[2];
sx q[2];
rz(-2.3818805) q[2];
rz(2.4921913) q[3];
sx q[3];
rz(-1.6067182) q[3];
sx q[3];
rz(-1.5788797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11288697) q[0];
sx q[0];
rz(-2.4876471) q[0];
sx q[0];
rz(1.1829859) q[0];
rz(0.68823367) q[1];
sx q[1];
rz(-1.3166683) q[1];
sx q[1];
rz(-2.7722955) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6128176) q[0];
sx q[0];
rz(-1.5982096) q[0];
sx q[0];
rz(-2.4779336) q[0];
x q[1];
rz(-0.30226548) q[2];
sx q[2];
rz(-0.73469732) q[2];
sx q[2];
rz(0.96978984) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6162045) q[1];
sx q[1];
rz(-1.3862351) q[1];
sx q[1];
rz(2.2730458) q[1];
rz(-0.7820205) q[3];
sx q[3];
rz(-0.45392515) q[3];
sx q[3];
rz(0.63185121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1987622) q[2];
sx q[2];
rz(-1.8161215) q[2];
sx q[2];
rz(2.14373) q[2];
rz(2.5241847) q[3];
sx q[3];
rz(-2.182775) q[3];
sx q[3];
rz(2.3203885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.038641039) q[0];
sx q[0];
rz(-1.8674253) q[0];
sx q[0];
rz(-1.2432903) q[0];
rz(2.359911) q[1];
sx q[1];
rz(-1.2739173) q[1];
sx q[1];
rz(2.8807358) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2170625) q[0];
sx q[0];
rz(-2.6017761) q[0];
sx q[0];
rz(-0.41882078) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.37520295) q[2];
sx q[2];
rz(-1.979036) q[2];
sx q[2];
rz(-0.21438504) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.971839) q[1];
sx q[1];
rz(-2.0287645) q[1];
sx q[1];
rz(-2.9724246) q[1];
rz(-1.2844159) q[3];
sx q[3];
rz(-0.65154159) q[3];
sx q[3];
rz(-2.1809354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6692052) q[2];
sx q[2];
rz(-0.86296764) q[2];
sx q[2];
rz(0.48259398) q[2];
rz(3.0274296) q[3];
sx q[3];
rz(-1.0897021) q[3];
sx q[3];
rz(-2.193006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77940762) q[0];
sx q[0];
rz(-3.0857093) q[0];
sx q[0];
rz(0.26595297) q[0];
rz(0.42568046) q[1];
sx q[1];
rz(-1.8961779) q[1];
sx q[1];
rz(0.46806213) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6702995) q[0];
sx q[0];
rz(-0.54139304) q[0];
sx q[0];
rz(1.5575717) q[0];
x q[1];
rz(1.2301679) q[2];
sx q[2];
rz(-1.6583031) q[2];
sx q[2];
rz(0.82078314) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7167042) q[1];
sx q[1];
rz(-2.0603097) q[1];
sx q[1];
rz(-0.29177006) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0909833) q[3];
sx q[3];
rz(-1.5438296) q[3];
sx q[3];
rz(-0.20668465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.257306) q[2];
sx q[2];
rz(-1.2904737) q[2];
sx q[2];
rz(-2.5878944) q[2];
rz(0.92343679) q[3];
sx q[3];
rz(-0.90673509) q[3];
sx q[3];
rz(3.0807909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4629352) q[0];
sx q[0];
rz(-2.904992) q[0];
sx q[0];
rz(2.4635354) q[0];
rz(0.97738114) q[1];
sx q[1];
rz(-1.9934318) q[1];
sx q[1];
rz(1.703702) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9674953) q[0];
sx q[0];
rz(-2.6405442) q[0];
sx q[0];
rz(-1.1606085) q[0];
x q[1];
rz(-3.1372141) q[2];
sx q[2];
rz(-2.5776049) q[2];
sx q[2];
rz(2.0019238) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9919093) q[1];
sx q[1];
rz(-1.339773) q[1];
sx q[1];
rz(0.2607338) q[1];
rz(-pi) q[2];
rz(2.1071042) q[3];
sx q[3];
rz(-2.2429129) q[3];
sx q[3];
rz(1.514707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4108654) q[2];
sx q[2];
rz(-2.5812456) q[2];
sx q[2];
rz(2.4746258) q[2];
rz(-0.11735958) q[3];
sx q[3];
rz(-1.2753692) q[3];
sx q[3];
rz(1.7608775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9475107) q[0];
sx q[0];
rz(-1.5912594) q[0];
sx q[0];
rz(-2.7981753) q[0];
rz(1.6471479) q[1];
sx q[1];
rz(-0.98973715) q[1];
sx q[1];
rz(2.6572878) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0750879) q[0];
sx q[0];
rz(-2.4874788) q[0];
sx q[0];
rz(2.6978786) q[0];
x q[1];
rz(-1.9860007) q[2];
sx q[2];
rz(-1.9671429) q[2];
sx q[2];
rz(1.1434778) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2769988) q[1];
sx q[1];
rz(-2.5127257) q[1];
sx q[1];
rz(0.38284812) q[1];
rz(0.24933322) q[3];
sx q[3];
rz(-1.3568546) q[3];
sx q[3];
rz(-2.6765424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1962428) q[2];
sx q[2];
rz(-2.0413155) q[2];
sx q[2];
rz(2.2958882) q[2];
rz(-1.2196994) q[3];
sx q[3];
rz(-1.889735) q[3];
sx q[3];
rz(2.9367101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6824816) q[0];
sx q[0];
rz(-2.8963608) q[0];
sx q[0];
rz(0.58897513) q[0];
rz(-2.466195) q[1];
sx q[1];
rz(-2.1928936) q[1];
sx q[1];
rz(-1.4640456) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1896473) q[0];
sx q[0];
rz(-2.4485817) q[0];
sx q[0];
rz(-0.89784867) q[0];
rz(-pi) q[1];
rz(1.6590933) q[2];
sx q[2];
rz(-0.73511926) q[2];
sx q[2];
rz(0.59150782) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3879536) q[1];
sx q[1];
rz(-1.6170701) q[1];
sx q[1];
rz(-0.58677499) q[1];
rz(-pi) q[2];
rz(-2.4469763) q[3];
sx q[3];
rz(-1.9362078) q[3];
sx q[3];
rz(-1.409274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.4097164) q[2];
sx q[2];
rz(-1.7475374) q[2];
sx q[2];
rz(2.0173006) q[2];
rz(-0.8479979) q[3];
sx q[3];
rz(-2.336899) q[3];
sx q[3];
rz(-0.030357411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
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
rz(2.547434) q[0];
sx q[0];
rz(-1.1347329) q[0];
sx q[0];
rz(-1.5190079) q[0];
rz(-2.2183954) q[1];
sx q[1];
rz(-1.0293488) q[1];
sx q[1];
rz(1.6346288) q[1];
rz(0.53096622) q[2];
sx q[2];
rz(-0.35561564) q[2];
sx q[2];
rz(1.8993062) q[2];
rz(-0.22696678) q[3];
sx q[3];
rz(-1.2746779) q[3];
sx q[3];
rz(-2.6770491) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
