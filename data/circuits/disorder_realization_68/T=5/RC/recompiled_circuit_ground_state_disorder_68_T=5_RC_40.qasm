OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5691583) q[0];
sx q[0];
rz(-0.99826607) q[0];
sx q[0];
rz(-2.7531667) q[0];
rz(-1.2216964) q[1];
sx q[1];
rz(0.95284) q[1];
sx q[1];
rz(9.3508773) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.088053) q[0];
sx q[0];
rz(-2.5754991) q[0];
sx q[0];
rz(-0.38882701) q[0];
rz(-pi) q[1];
rz(2.2227395) q[2];
sx q[2];
rz(-1.5519189) q[2];
sx q[2];
rz(-2.9110661) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2683289) q[1];
sx q[1];
rz(-1.2115098) q[1];
sx q[1];
rz(-2.2993035) q[1];
x q[2];
rz(0.2147712) q[3];
sx q[3];
rz(-1.8146551) q[3];
sx q[3];
rz(-2.7726695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.6834324) q[2];
sx q[2];
rz(-1.4417803) q[2];
sx q[2];
rz(-1.3275576) q[2];
rz(-2.1761927) q[3];
sx q[3];
rz(-0.5894956) q[3];
sx q[3];
rz(-2.086967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4210159) q[0];
sx q[0];
rz(-0.0076616658) q[0];
sx q[0];
rz(-1.3910008) q[0];
rz(1.1945456) q[1];
sx q[1];
rz(-0.60940131) q[1];
sx q[1];
rz(2.4360099) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53853453) q[0];
sx q[0];
rz(-1.0760191) q[0];
sx q[0];
rz(0.81788466) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6383361) q[2];
sx q[2];
rz(-0.76336289) q[2];
sx q[2];
rz(-0.22035519) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0944453) q[1];
sx q[1];
rz(-1.021968) q[1];
sx q[1];
rz(2.7286163) q[1];
rz(-pi) q[2];
rz(-1.0666177) q[3];
sx q[3];
rz(-2.1103577) q[3];
sx q[3];
rz(-2.2062796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0671063) q[2];
sx q[2];
rz(-1.4718082) q[2];
sx q[2];
rz(-0.86522317) q[2];
rz(1.5857961) q[3];
sx q[3];
rz(-0.14388789) q[3];
sx q[3];
rz(-1.5998862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(2.9928352) q[0];
sx q[0];
rz(-2.302763) q[0];
sx q[0];
rz(1.6194153) q[0];
rz(-1.4083699) q[1];
sx q[1];
rz(-2.759628) q[1];
sx q[1];
rz(0.24040374) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5071538) q[0];
sx q[0];
rz(-2.2858862) q[0];
sx q[0];
rz(1.2583744) q[0];
rz(1.1454606) q[2];
sx q[2];
rz(-2.7096203) q[2];
sx q[2];
rz(-1.7329104) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.63096365) q[1];
sx q[1];
rz(-1.6121702) q[1];
sx q[1];
rz(-3.0363185) q[1];
rz(-2.6601276) q[3];
sx q[3];
rz(-1.4507308) q[3];
sx q[3];
rz(2.1019746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0092885) q[2];
sx q[2];
rz(-1.665364) q[2];
sx q[2];
rz(-1.08584) q[2];
rz(0.050431937) q[3];
sx q[3];
rz(-1.4112873) q[3];
sx q[3];
rz(2.1850695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0352935) q[0];
sx q[0];
rz(-1.4426008) q[0];
sx q[0];
rz(-2.371149) q[0];
rz(-2.2162407) q[1];
sx q[1];
rz(-1.4848361) q[1];
sx q[1];
rz(0.83522183) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5073287) q[0];
sx q[0];
rz(-1.5871704) q[0];
sx q[0];
rz(-2.6446675) q[0];
x q[1];
rz(-0.24150924) q[2];
sx q[2];
rz(-1.6577953) q[2];
sx q[2];
rz(2.3149025) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.861117) q[1];
sx q[1];
rz(-1.0442227) q[1];
sx q[1];
rz(3.1269424) q[1];
rz(-pi) q[2];
rz(0.20056574) q[3];
sx q[3];
rz(-1.4576389) q[3];
sx q[3];
rz(3.1386047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.21021065) q[2];
sx q[2];
rz(-1.8385734) q[2];
sx q[2];
rz(-2.4118928) q[2];
rz(-0.99916712) q[3];
sx q[3];
rz(-1.3068643) q[3];
sx q[3];
rz(2.6877747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5900742) q[0];
sx q[0];
rz(-0.22576627) q[0];
sx q[0];
rz(0.99910587) q[0];
rz(-2.0935811) q[1];
sx q[1];
rz(-0.71185714) q[1];
sx q[1];
rz(-0.5447095) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9182129) q[0];
sx q[0];
rz(-1.8210016) q[0];
sx q[0];
rz(1.3859207) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4757521) q[2];
sx q[2];
rz(-2.2975635) q[2];
sx q[2];
rz(-0.11451463) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7919257) q[1];
sx q[1];
rz(-1.3832258) q[1];
sx q[1];
rz(-0.55123185) q[1];
rz(-2.7361511) q[3];
sx q[3];
rz(-1.573994) q[3];
sx q[3];
rz(-2.1759667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5576632) q[2];
sx q[2];
rz(-1.9926535) q[2];
sx q[2];
rz(-1.0168797) q[2];
rz(-2.478638) q[3];
sx q[3];
rz(-0.69279492) q[3];
sx q[3];
rz(-1.3346765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3589631) q[0];
sx q[0];
rz(-2.6435659) q[0];
sx q[0];
rz(1.9899415) q[0];
rz(-0.89371124) q[1];
sx q[1];
rz(-1.3795373) q[1];
sx q[1];
rz(0.11633565) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2710909) q[0];
sx q[0];
rz(-1.4634824) q[0];
sx q[0];
rz(-2.0475325) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6200368) q[2];
sx q[2];
rz(-0.34477012) q[2];
sx q[2];
rz(1.5695968) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.92116881) q[1];
sx q[1];
rz(-1.0225778) q[1];
sx q[1];
rz(2.9698644) q[1];
rz(-1.6508914) q[3];
sx q[3];
rz(-2.3175196) q[3];
sx q[3];
rz(-2.2585845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5748888) q[2];
sx q[2];
rz(-2.7723007) q[2];
sx q[2];
rz(1.0075547) q[2];
rz(1.9762074) q[3];
sx q[3];
rz(-0.27519614) q[3];
sx q[3];
rz(-2.5113441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4646869) q[0];
sx q[0];
rz(-0.46606627) q[0];
sx q[0];
rz(-2.9718072) q[0];
rz(-0.67974293) q[1];
sx q[1];
rz(-2.3092473) q[1];
sx q[1];
rz(-0.92175093) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0089192275) q[0];
sx q[0];
rz(-0.28604315) q[0];
sx q[0];
rz(-0.45222767) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5160732) q[2];
sx q[2];
rz(-1.2642908) q[2];
sx q[2];
rz(0.79105575) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3231735) q[1];
sx q[1];
rz(-1.9193279) q[1];
sx q[1];
rz(1.0437832) q[1];
rz(1.3774817) q[3];
sx q[3];
rz(-1.3018248) q[3];
sx q[3];
rz(1.811816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0035231) q[2];
sx q[2];
rz(-1.6602844) q[2];
sx q[2];
rz(1.8180234) q[2];
rz(2.0924163) q[3];
sx q[3];
rz(-2.0101571) q[3];
sx q[3];
rz(-1.3808892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2915989) q[0];
sx q[0];
rz(-1.1350564) q[0];
sx q[0];
rz(-2.3681613) q[0];
rz(2.6749581) q[1];
sx q[1];
rz(-2.0687053) q[1];
sx q[1];
rz(3.1007865) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.651603) q[0];
sx q[0];
rz(-1.7035653) q[0];
sx q[0];
rz(0.013057166) q[0];
rz(-pi) q[1];
rz(-2.2757024) q[2];
sx q[2];
rz(-0.84176062) q[2];
sx q[2];
rz(1.2129606) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0405214) q[1];
sx q[1];
rz(-0.98993976) q[1];
sx q[1];
rz(3.0809666) q[1];
rz(-pi) q[2];
x q[2];
rz(0.17191903) q[3];
sx q[3];
rz(-2.1287595) q[3];
sx q[3];
rz(-2.362988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1524973) q[2];
sx q[2];
rz(-0.54741198) q[2];
sx q[2];
rz(-3.1067749) q[2];
rz(-3.1133437) q[3];
sx q[3];
rz(-2.0939128) q[3];
sx q[3];
rz(0.69407535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5038274) q[0];
sx q[0];
rz(-2.4389508) q[0];
sx q[0];
rz(1.0850061) q[0];
rz(1.7833692) q[1];
sx q[1];
rz(-0.63301507) q[1];
sx q[1];
rz(-1.0795275) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87303783) q[0];
sx q[0];
rz(-1.3655021) q[0];
sx q[0];
rz(2.1389066) q[0];
rz(-pi) q[1];
rz(2.8582879) q[2];
sx q[2];
rz(-0.46479169) q[2];
sx q[2];
rz(-1.9329485) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.72708817) q[1];
sx q[1];
rz(-1.8782756) q[1];
sx q[1];
rz(-0.25924637) q[1];
rz(-pi) q[2];
rz(-0.50864403) q[3];
sx q[3];
rz(-2.2318342) q[3];
sx q[3];
rz(2.4511372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8980155) q[2];
sx q[2];
rz(-1.5960627) q[2];
sx q[2];
rz(3.0698981) q[2];
rz(2.535635) q[3];
sx q[3];
rz(-2.3809483) q[3];
sx q[3];
rz(3.0683556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83409413) q[0];
sx q[0];
rz(-1.6313169) q[0];
sx q[0];
rz(-2.639005) q[0];
rz(-1.3692299) q[1];
sx q[1];
rz(-1.9160198) q[1];
sx q[1];
rz(-0.1756846) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7384199) q[0];
sx q[0];
rz(-1.4851928) q[0];
sx q[0];
rz(1.9001207) q[0];
rz(-pi) q[1];
rz(-2.530859) q[2];
sx q[2];
rz(-2.1977976) q[2];
sx q[2];
rz(-0.58525733) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8524234) q[1];
sx q[1];
rz(-2.2383177) q[1];
sx q[1];
rz(-0.65594419) q[1];
rz(-pi) q[2];
rz(0.93264742) q[3];
sx q[3];
rz(-1.89092) q[3];
sx q[3];
rz(0.29039106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.48675576) q[2];
sx q[2];
rz(-1.882587) q[2];
sx q[2];
rz(-0.34461018) q[2];
rz(-1.2279855) q[3];
sx q[3];
rz(-2.4495008) q[3];
sx q[3];
rz(-2.1740348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
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
rz(1.6268613) q[0];
sx q[0];
rz(-1.8296965) q[0];
sx q[0];
rz(-2.1924023) q[0];
rz(1.7687891) q[1];
sx q[1];
rz(-0.95284843) q[1];
sx q[1];
rz(1.3292809) q[1];
rz(-0.71277724) q[2];
sx q[2];
rz(-1.7736802) q[2];
sx q[2];
rz(1.6414248) q[2];
rz(2.9459841) q[3];
sx q[3];
rz(-0.8630639) q[3];
sx q[3];
rz(0.56474781) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
