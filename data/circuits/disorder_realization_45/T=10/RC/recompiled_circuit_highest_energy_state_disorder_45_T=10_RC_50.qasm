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
rz(-1.2019914) q[0];
sx q[0];
rz(3.6245873) q[0];
sx q[0];
rz(10.935187) q[0];
rz(-0.13934879) q[1];
sx q[1];
rz(5.723602) q[1];
sx q[1];
rz(8.6948123) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18099526) q[0];
sx q[0];
rz(-0.9708403) q[0];
sx q[0];
rz(0.17790312) q[0];
rz(2.5467186) q[2];
sx q[2];
rz(-1.7584137) q[2];
sx q[2];
rz(-1.3485731) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.12355676) q[1];
sx q[1];
rz(-2.5606321) q[1];
sx q[1];
rz(3.0835129) q[1];
rz(-pi) q[2];
rz(1.8584941) q[3];
sx q[3];
rz(-1.3511063) q[3];
sx q[3];
rz(2.2030597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1175179) q[2];
sx q[2];
rz(-2.8705609) q[2];
sx q[2];
rz(0.81895858) q[2];
rz(-3.135318) q[3];
sx q[3];
rz(-1.2333906) q[3];
sx q[3];
rz(2.1591469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55945021) q[0];
sx q[0];
rz(-1.5511976) q[0];
sx q[0];
rz(-2.5421802) q[0];
rz(0.86743152) q[1];
sx q[1];
rz(-2.1116833) q[1];
sx q[1];
rz(2.3960466) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.065577995) q[0];
sx q[0];
rz(-1.4434555) q[0];
sx q[0];
rz(-2.8498184) q[0];
rz(0.7600766) q[2];
sx q[2];
rz(-1.809263) q[2];
sx q[2];
rz(-2.5395405) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9245431) q[1];
sx q[1];
rz(-0.94247183) q[1];
sx q[1];
rz(1.373686) q[1];
rz(-pi) q[2];
rz(-2.1576513) q[3];
sx q[3];
rz(-2.3592257) q[3];
sx q[3];
rz(-2.7593105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6856689) q[2];
sx q[2];
rz(-0.29752877) q[2];
sx q[2];
rz(1.7342742) q[2];
rz(2.2972441) q[3];
sx q[3];
rz(-0.83365369) q[3];
sx q[3];
rz(2.1048529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5832962) q[0];
sx q[0];
rz(-1.1953657) q[0];
sx q[0];
rz(-1.012828) q[0];
rz(-0.60802513) q[1];
sx q[1];
rz(-1.5826179) q[1];
sx q[1];
rz(1.2695405) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69508775) q[0];
sx q[0];
rz(-1.9534612) q[0];
sx q[0];
rz(-1.0951359) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.059194) q[2];
sx q[2];
rz(-2.2401056) q[2];
sx q[2];
rz(1.5562039) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3568282) q[1];
sx q[1];
rz(-0.85961378) q[1];
sx q[1];
rz(1.8309092) q[1];
rz(-pi) q[2];
rz(0.57165159) q[3];
sx q[3];
rz(-0.56662512) q[3];
sx q[3];
rz(-1.0039312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.515392) q[2];
sx q[2];
rz(-2.0833368) q[2];
sx q[2];
rz(-1.2916279) q[2];
rz(-3.1332704) q[3];
sx q[3];
rz(-2.1885927) q[3];
sx q[3];
rz(-0.21361175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0030647) q[0];
sx q[0];
rz(-0.55532885) q[0];
sx q[0];
rz(1.3767161) q[0];
rz(-1.6142913) q[1];
sx q[1];
rz(-1.9381899) q[1];
sx q[1];
rz(-0.52070224) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10258612) q[0];
sx q[0];
rz(-1.9780567) q[0];
sx q[0];
rz(-0.38577052) q[0];
x q[1];
rz(2.1200373) q[2];
sx q[2];
rz(-2.5168572) q[2];
sx q[2];
rz(0.83323375) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5129052) q[1];
sx q[1];
rz(-0.73343745) q[1];
sx q[1];
rz(0.67104407) q[1];
x q[2];
rz(2.2124452) q[3];
sx q[3];
rz(-2.3722674) q[3];
sx q[3];
rz(-1.820946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.5248096) q[2];
sx q[2];
rz(-0.7889792) q[2];
sx q[2];
rz(1.7899803) q[2];
rz(2.1221519) q[3];
sx q[3];
rz(-2.5717058) q[3];
sx q[3];
rz(-1.1714237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5919507) q[0];
sx q[0];
rz(-1.4627946) q[0];
sx q[0];
rz(0.098966448) q[0];
rz(0.70676604) q[1];
sx q[1];
rz(-2.2711429) q[1];
sx q[1];
rz(-0.4695355) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81585138) q[0];
sx q[0];
rz(-1.5945863) q[0];
sx q[0];
rz(1.7086534) q[0];
rz(-pi) q[1];
rz(-0.37920538) q[2];
sx q[2];
rz(-2.0742886) q[2];
sx q[2];
rz(-0.45631726) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.27891913) q[1];
sx q[1];
rz(-1.1871183) q[1];
sx q[1];
rz(0.63273276) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9088777) q[3];
sx q[3];
rz(-2.5732627) q[3];
sx q[3];
rz(-1.1956904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2609451) q[2];
sx q[2];
rz(-1.3346883) q[2];
sx q[2];
rz(-1.3060695) q[2];
rz(0.4246873) q[3];
sx q[3];
rz(-1.7644019) q[3];
sx q[3];
rz(0.88596058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11923085) q[0];
sx q[0];
rz(-2.5194118) q[0];
sx q[0];
rz(1.3223883) q[0];
rz(2.0139096) q[1];
sx q[1];
rz(-1.6219982) q[1];
sx q[1];
rz(1.7599531) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36395006) q[0];
sx q[0];
rz(-0.76350437) q[0];
sx q[0];
rz(1.6958773) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3076129) q[2];
sx q[2];
rz(-2.4206941) q[2];
sx q[2];
rz(-2.3490459) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0387011) q[1];
sx q[1];
rz(-0.82566092) q[1];
sx q[1];
rz(1.3172512) q[1];
rz(2.9074677) q[3];
sx q[3];
rz(-0.30933274) q[3];
sx q[3];
rz(2.7721289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7606925) q[2];
sx q[2];
rz(-1.6221294) q[2];
sx q[2];
rz(-2.6603928) q[2];
rz(0.5091269) q[3];
sx q[3];
rz(-2.9042518) q[3];
sx q[3];
rz(2.7595162) q[3];
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
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5612438) q[0];
sx q[0];
rz(-2.6600397) q[0];
sx q[0];
rz(-3.0522108) q[0];
rz(-2.0629758) q[1];
sx q[1];
rz(-1.4559454) q[1];
sx q[1];
rz(-2.1481029) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2334796) q[0];
sx q[0];
rz(-2.3589239) q[0];
sx q[0];
rz(-3.0872869) q[0];
rz(0.78623505) q[2];
sx q[2];
rz(-2.1053542) q[2];
sx q[2];
rz(1.3649808) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1878933) q[1];
sx q[1];
rz(-2.7479798) q[1];
sx q[1];
rz(1.1792437) q[1];
rz(-2.4782789) q[3];
sx q[3];
rz(-0.74771008) q[3];
sx q[3];
rz(1.472689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.79954687) q[2];
sx q[2];
rz(-0.61508238) q[2];
sx q[2];
rz(0.40204027) q[2];
rz(-1.0953995) q[3];
sx q[3];
rz(-1.9270555) q[3];
sx q[3];
rz(0.57797617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6680229) q[0];
sx q[0];
rz(-2.3325925) q[0];
sx q[0];
rz(-1.8079669) q[0];
rz(-0.85583055) q[1];
sx q[1];
rz(-1.4060833) q[1];
sx q[1];
rz(-2.2241101) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3429392) q[0];
sx q[0];
rz(-2.9173033) q[0];
sx q[0];
rz(-2.0464315) q[0];
x q[1];
rz(-2.1945004) q[2];
sx q[2];
rz(-1.6674041) q[2];
sx q[2];
rz(-1.7283224) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5405824) q[1];
sx q[1];
rz(-1.2262019) q[1];
sx q[1];
rz(0.2356727) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0794271) q[3];
sx q[3];
rz(-2.1537703) q[3];
sx q[3];
rz(-0.82750083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.63460073) q[2];
sx q[2];
rz(-1.4054106) q[2];
sx q[2];
rz(1.8514006) q[2];
rz(-2.2318132) q[3];
sx q[3];
rz(-1.6981643) q[3];
sx q[3];
rz(-3.1006052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9057587) q[0];
sx q[0];
rz(-1.9941149) q[0];
sx q[0];
rz(-0.27715096) q[0];
rz(1.9174891) q[1];
sx q[1];
rz(-1.6033019) q[1];
sx q[1];
rz(-1.7701497) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0540028) q[0];
sx q[0];
rz(-1.3773019) q[0];
sx q[0];
rz(2.4688379) q[0];
x q[1];
rz(-1.1354228) q[2];
sx q[2];
rz(-2.0757542) q[2];
sx q[2];
rz(-1.4471042) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0321694) q[1];
sx q[1];
rz(-2.210683) q[1];
sx q[1];
rz(0.25856047) q[1];
rz(-pi) q[2];
rz(-1.7684494) q[3];
sx q[3];
rz(-0.73171333) q[3];
sx q[3];
rz(0.55303516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.203043) q[2];
sx q[2];
rz(-2.2800192) q[2];
sx q[2];
rz(2.4533563) q[2];
rz(-2.8271683) q[3];
sx q[3];
rz(-0.36948547) q[3];
sx q[3];
rz(1.2615874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37453434) q[0];
sx q[0];
rz(-2.2735167) q[0];
sx q[0];
rz(0.57149291) q[0];
rz(2.4608965) q[1];
sx q[1];
rz(-1.9711767) q[1];
sx q[1];
rz(-0.64819711) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8777953) q[0];
sx q[0];
rz(-1.5476942) q[0];
sx q[0];
rz(0.013851555) q[0];
rz(-pi) q[1];
rz(0.86096455) q[2];
sx q[2];
rz(-2.2662244) q[2];
sx q[2];
rz(-2.8784213) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1163016) q[1];
sx q[1];
rz(-0.87193438) q[1];
sx q[1];
rz(1.541733) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6246214) q[3];
sx q[3];
rz(-1.3992157) q[3];
sx q[3];
rz(1.3763381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8615243) q[2];
sx q[2];
rz(-1.0736977) q[2];
sx q[2];
rz(1.6746707) q[2];
rz(0.17659771) q[3];
sx q[3];
rz(-0.64591518) q[3];
sx q[3];
rz(1.8471898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27547729) q[0];
sx q[0];
rz(-1.3963516) q[0];
sx q[0];
rz(1.8657952) q[0];
rz(2.7571309) q[1];
sx q[1];
rz(-1.6356331) q[1];
sx q[1];
rz(2.5148139) q[1];
rz(-2.7914417) q[2];
sx q[2];
rz(-1.2731009) q[2];
sx q[2];
rz(0.44558744) q[2];
rz(2.1351142) q[3];
sx q[3];
rz(-1.1075533) q[3];
sx q[3];
rz(1.5816734) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
