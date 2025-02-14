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
rz(0.25207818) q[0];
sx q[0];
rz(-2.5660388) q[0];
sx q[0];
rz(-0.12230305) q[0];
rz(1.1072493) q[1];
sx q[1];
rz(-0.9613494) q[1];
sx q[1];
rz(-0.038318757) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.092889812) q[0];
sx q[0];
rz(-0.99703353) q[0];
sx q[0];
rz(-0.4452197) q[0];
rz(-pi) q[1];
rz(-1.4517446) q[2];
sx q[2];
rz(-1.6887293) q[2];
sx q[2];
rz(2.8182507) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.029509228) q[1];
sx q[1];
rz(-1.8942299) q[1];
sx q[1];
rz(-0.98677633) q[1];
rz(-0.24738048) q[3];
sx q[3];
rz(-1.5508964) q[3];
sx q[3];
rz(2.6008391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.64500874) q[2];
sx q[2];
rz(-1.4270447) q[2];
sx q[2];
rz(-1.4487779) q[2];
rz(-0.19041666) q[3];
sx q[3];
rz(-1.0551635) q[3];
sx q[3];
rz(2.9922488) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9377624) q[0];
sx q[0];
rz(-1.2685403) q[0];
sx q[0];
rz(0.25800905) q[0];
rz(-1.5915271) q[1];
sx q[1];
rz(-1.9015046) q[1];
sx q[1];
rz(0.16664997) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7014871) q[0];
sx q[0];
rz(-0.99033725) q[0];
sx q[0];
rz(1.6850182) q[0];
x q[1];
rz(-2.7562253) q[2];
sx q[2];
rz(-1.3056985) q[2];
sx q[2];
rz(-2.6814658) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7609561) q[1];
sx q[1];
rz(-1.445862) q[1];
sx q[1];
rz(0.71604587) q[1];
rz(-2.1364501) q[3];
sx q[3];
rz(-2.0453302) q[3];
sx q[3];
rz(2.8090854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0717281) q[2];
sx q[2];
rz(-1.122415) q[2];
sx q[2];
rz(1.9094763) q[2];
rz(2.2169436) q[3];
sx q[3];
rz(-2.2053714) q[3];
sx q[3];
rz(1.1054976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.449618) q[0];
sx q[0];
rz(-1.469935) q[0];
sx q[0];
rz(-0.0019419226) q[0];
rz(-0.034189668) q[1];
sx q[1];
rz(-1.2119774) q[1];
sx q[1];
rz(-1.5431822) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40638518) q[0];
sx q[0];
rz(-1.8282561) q[0];
sx q[0];
rz(-2.6763112) q[0];
x q[1];
rz(2.7091647) q[2];
sx q[2];
rz(-0.78410599) q[2];
sx q[2];
rz(2.6570005) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0405827) q[1];
sx q[1];
rz(-1.9242745) q[1];
sx q[1];
rz(1.6124875) q[1];
rz(-pi) q[2];
rz(1.8887599) q[3];
sx q[3];
rz(-2.468716) q[3];
sx q[3];
rz(1.3313791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.89692846) q[2];
sx q[2];
rz(-2.7074773) q[2];
sx q[2];
rz(-2.1042692) q[2];
rz(1.1050998) q[3];
sx q[3];
rz(-1.5367855) q[3];
sx q[3];
rz(2.0028116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26745519) q[0];
sx q[0];
rz(-2.656811) q[0];
sx q[0];
rz(-1.4111891) q[0];
rz(0.63938582) q[1];
sx q[1];
rz(-1.6849898) q[1];
sx q[1];
rz(-3.0944518) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77098504) q[0];
sx q[0];
rz(-1.8520675) q[0];
sx q[0];
rz(-0.60312834) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7732255) q[2];
sx q[2];
rz(-0.50058156) q[2];
sx q[2];
rz(-2.6427302) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6005391) q[1];
sx q[1];
rz(-0.70912921) q[1];
sx q[1];
rz(2.1313271) q[1];
rz(-2.8465038) q[3];
sx q[3];
rz(-2.2722244) q[3];
sx q[3];
rz(-1.5257344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3186657) q[2];
sx q[2];
rz(-1.1456127) q[2];
sx q[2];
rz(-1.9226496) q[2];
rz(-0.31560358) q[3];
sx q[3];
rz(-0.054840755) q[3];
sx q[3];
rz(0.1489197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8197935) q[0];
sx q[0];
rz(-0.55235523) q[0];
sx q[0];
rz(0.34859443) q[0];
rz(1.2527342) q[1];
sx q[1];
rz(-1.1628954) q[1];
sx q[1];
rz(1.3016275) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5900629) q[0];
sx q[0];
rz(-2.3746535) q[0];
sx q[0];
rz(1.4566896) q[0];
rz(1.1548066) q[2];
sx q[2];
rz(-2.092053) q[2];
sx q[2];
rz(-0.74300569) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0730473) q[1];
sx q[1];
rz(-1.3497735) q[1];
sx q[1];
rz(2.9566492) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0318448) q[3];
sx q[3];
rz(-2.5205118) q[3];
sx q[3];
rz(-0.10824848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9584413) q[2];
sx q[2];
rz(-0.30854598) q[2];
sx q[2];
rz(1.4754254) q[2];
rz(-2.4557377) q[3];
sx q[3];
rz(-0.97969046) q[3];
sx q[3];
rz(-2.6077152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0297861) q[0];
sx q[0];
rz(-2.989558) q[0];
sx q[0];
rz(-1.6492122) q[0];
rz(0.97914186) q[1];
sx q[1];
rz(-1.5359595) q[1];
sx q[1];
rz(1.2243366) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2164152) q[0];
sx q[0];
rz(-1.6013711) q[0];
sx q[0];
rz(-2.9435817) q[0];
rz(2.5923205) q[2];
sx q[2];
rz(-1.3830997) q[2];
sx q[2];
rz(-0.21085462) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8414014) q[1];
sx q[1];
rz(-2.527664) q[1];
sx q[1];
rz(1.5924686) q[1];
rz(-pi) q[2];
rz(2.0322509) q[3];
sx q[3];
rz(-0.68285817) q[3];
sx q[3];
rz(2.7340555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7225723) q[2];
sx q[2];
rz(-1.6442862) q[2];
sx q[2];
rz(0.9160308) q[2];
rz(1.3845059) q[3];
sx q[3];
rz(-1.2576831) q[3];
sx q[3];
rz(-0.70703435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1322121) q[0];
sx q[0];
rz(-0.051932422) q[0];
sx q[0];
rz(0.95426553) q[0];
rz(-0.43680278) q[1];
sx q[1];
rz(-1.5635468) q[1];
sx q[1];
rz(0.11016914) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1230135) q[0];
sx q[0];
rz(-1.2467822) q[0];
sx q[0];
rz(-1.8229681) q[0];
rz(2.4442441) q[2];
sx q[2];
rz(-1.9983872) q[2];
sx q[2];
rz(-0.22874895) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4272144) q[1];
sx q[1];
rz(-2.7958779) q[1];
sx q[1];
rz(1.5835254) q[1];
rz(-pi) q[2];
rz(-1.0694592) q[3];
sx q[3];
rz(-1.3469463) q[3];
sx q[3];
rz(1.4042133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.43180141) q[2];
sx q[2];
rz(-3.1352037) q[2];
sx q[2];
rz(2.1978281) q[2];
rz(-0.56378311) q[3];
sx q[3];
rz(-2.0457334) q[3];
sx q[3];
rz(-2.0311267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7875882) q[0];
sx q[0];
rz(-0.28689757) q[0];
sx q[0];
rz(0.34550825) q[0];
rz(1.0954789) q[1];
sx q[1];
rz(-1.4778719) q[1];
sx q[1];
rz(1.5617721) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6248572) q[0];
sx q[0];
rz(-2.3586732) q[0];
sx q[0];
rz(1.5262414) q[0];
rz(0.061959668) q[2];
sx q[2];
rz(-0.91455787) q[2];
sx q[2];
rz(-1.3861173) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6553538) q[1];
sx q[1];
rz(-1.2062688) q[1];
sx q[1];
rz(0.93185394) q[1];
rz(-pi) q[2];
rz(1.5388266) q[3];
sx q[3];
rz(-1.073146) q[3];
sx q[3];
rz(-2.9835193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.76274189) q[2];
sx q[2];
rz(-0.26143917) q[2];
sx q[2];
rz(1.7108819) q[2];
rz(0.53449574) q[3];
sx q[3];
rz(-1.675324) q[3];
sx q[3];
rz(-2.1889595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5296103) q[0];
sx q[0];
rz(-3.072325) q[0];
sx q[0];
rz(1.8310504) q[0];
rz(2.2546841) q[1];
sx q[1];
rz(-0.84019089) q[1];
sx q[1];
rz(1.2695674) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5909605) q[0];
sx q[0];
rz(-2.7661588) q[0];
sx q[0];
rz(-1.9272789) q[0];
x q[1];
rz(-1.1940895) q[2];
sx q[2];
rz(-0.40222886) q[2];
sx q[2];
rz(-2.5144983) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.15358217) q[1];
sx q[1];
rz(-1.8021823) q[1];
sx q[1];
rz(1.850495) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6107913) q[3];
sx q[3];
rz(-1.3673395) q[3];
sx q[3];
rz(-1.1261765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5294007) q[2];
sx q[2];
rz(-1.3357013) q[2];
sx q[2];
rz(2.9086435) q[2];
rz(-1.285078) q[3];
sx q[3];
rz(-2.464747) q[3];
sx q[3];
rz(-2.7555833) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1583629) q[0];
sx q[0];
rz(-2.5322999) q[0];
sx q[0];
rz(-3.0443211) q[0];
rz(2.3376047) q[1];
sx q[1];
rz(-0.709788) q[1];
sx q[1];
rz(-0.21673094) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9602338) q[0];
sx q[0];
rz(-2.9621819) q[0];
sx q[0];
rz(1.8852194) q[0];
x q[1];
rz(2.9998776) q[2];
sx q[2];
rz(-1.5248305) q[2];
sx q[2];
rz(-1.7739997) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4137553) q[1];
sx q[1];
rz(-2.1247132) q[1];
sx q[1];
rz(0.78110771) q[1];
rz(0.43301591) q[3];
sx q[3];
rz(-2.4320514) q[3];
sx q[3];
rz(0.39329986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.29297605) q[2];
sx q[2];
rz(-0.28102195) q[2];
sx q[2];
rz(-0.49523735) q[2];
rz(1.5826591) q[3];
sx q[3];
rz(-2.0571183) q[3];
sx q[3];
rz(-2.9787279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0611298) q[0];
sx q[0];
rz(-1.6292138) q[0];
sx q[0];
rz(-0.52393352) q[0];
rz(-2.0551266) q[1];
sx q[1];
rz(-0.51849425) q[1];
sx q[1];
rz(0.35324221) q[1];
rz(-1.5389961) q[2];
sx q[2];
rz(-1.563579) q[2];
sx q[2];
rz(-2.210571) q[2];
rz(1.6752406) q[3];
sx q[3];
rz(-0.71487311) q[3];
sx q[3];
rz(-2.5345595) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
