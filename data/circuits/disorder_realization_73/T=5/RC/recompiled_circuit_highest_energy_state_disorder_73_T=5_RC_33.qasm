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
rz(0.90176982) q[0];
rz(1.3321441) q[1];
sx q[1];
rz(-2.6982215) q[1];
sx q[1];
rz(2.3763357) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9741824) q[0];
sx q[0];
rz(-1.9338134) q[0];
sx q[0];
rz(-2.1556751) q[0];
rz(2.1852605) q[2];
sx q[2];
rz(-2.7354089) q[2];
sx q[2];
rz(1.9910016) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.81232416) q[1];
sx q[1];
rz(-1.6109038) q[1];
sx q[1];
rz(-0.33848156) q[1];
rz(1.2918858) q[3];
sx q[3];
rz(-1.4456914) q[3];
sx q[3];
rz(-2.5549614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8454933) q[2];
sx q[2];
rz(-1.0881492) q[2];
sx q[2];
rz(1.6695401) q[2];
rz(-1.4269525) q[3];
sx q[3];
rz(-1.6441556) q[3];
sx q[3];
rz(0.53372395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.958309) q[0];
sx q[0];
rz(-1.6338209) q[0];
sx q[0];
rz(-0.84877745) q[0];
rz(0.035004184) q[1];
sx q[1];
rz(-1.1468381) q[1];
sx q[1];
rz(-2.1880207) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67259939) q[0];
sx q[0];
rz(-2.3777588) q[0];
sx q[0];
rz(2.4498778) q[0];
x q[1];
rz(-0.30470253) q[2];
sx q[2];
rz(-1.1977473) q[2];
sx q[2];
rz(2.3359131) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9417022) q[1];
sx q[1];
rz(-0.39559707) q[1];
sx q[1];
rz(-1.0426177) q[1];
x q[2];
rz(-0.75462975) q[3];
sx q[3];
rz(-2.0622232) q[3];
sx q[3];
rz(1.400465) q[3];
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
rz(-0.16544011) q[3];
sx q[3];
rz(-1.6980419) q[3];
sx q[3];
rz(2.3918242) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7716832) q[0];
sx q[0];
rz(-2.219438) q[0];
sx q[0];
rz(0.40192303) q[0];
rz(-0.15332128) q[1];
sx q[1];
rz(-2.9647277) q[1];
sx q[1];
rz(2.0053999) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83602609) q[0];
sx q[0];
rz(-1.5587806) q[0];
sx q[0];
rz(1.1712606) q[0];
rz(-pi) q[1];
rz(-1.8449144) q[2];
sx q[2];
rz(-1.1027567) q[2];
sx q[2];
rz(0.57690358) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.78627777) q[1];
sx q[1];
rz(-1.8031562) q[1];
sx q[1];
rz(1.8174328) q[1];
rz(-pi) q[2];
x q[2];
rz(0.46391507) q[3];
sx q[3];
rz(-1.216421) q[3];
sx q[3];
rz(2.7257435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0574657) q[2];
sx q[2];
rz(-2.2630313) q[2];
sx q[2];
rz(0.055056661) q[2];
rz(1.7940686) q[3];
sx q[3];
rz(-1.811458) q[3];
sx q[3];
rz(-1.0042892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
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
rz(-2.7771626) q[0];
sx q[0];
rz(-0.20474064) q[0];
sx q[0];
rz(-1.5675911) q[0];
rz(-2.4500627) q[1];
sx q[1];
rz(-2.4722996) q[1];
sx q[1];
rz(2.9561668) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7950993) q[0];
sx q[0];
rz(-0.43432626) q[0];
sx q[0];
rz(-2.7350635) q[0];
rz(-0.71132423) q[2];
sx q[2];
rz(-0.60718107) q[2];
sx q[2];
rz(1.173623) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1048829) q[1];
sx q[1];
rz(-1.2598902) q[1];
sx q[1];
rz(-0.26581146) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7774277) q[3];
sx q[3];
rz(-1.0344328) q[3];
sx q[3];
rz(0.11851507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7376248) q[2];
sx q[2];
rz(-1.4289958) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72652793) q[0];
sx q[0];
rz(-2.1580577) q[0];
sx q[0];
rz(2.5256185) q[0];
rz(-1.1888602) q[1];
sx q[1];
rz(-2.523246) q[1];
sx q[1];
rz(-0.51330769) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.05840551) q[0];
sx q[0];
rz(-1.8523916) q[0];
sx q[0];
rz(-0.0076936184) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8262324) q[2];
sx q[2];
rz(-2.3906446) q[2];
sx q[2];
rz(-0.89799228) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0541546) q[1];
sx q[1];
rz(-0.63870027) q[1];
sx q[1];
rz(3.020546) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0381472) q[3];
sx q[3];
rz(-1.6765521) q[3];
sx q[3];
rz(-2.8018564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0080228) q[2];
sx q[2];
rz(-1.2761389) q[2];
sx q[2];
rz(-2.3991154) q[2];
rz(1.4437458) q[3];
sx q[3];
rz(-1.4092813) q[3];
sx q[3];
rz(2.0402563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3105069) q[0];
sx q[0];
rz(-1.0467014) q[0];
sx q[0];
rz(0.35980862) q[0];
rz(2.2911435) q[1];
sx q[1];
rz(-1.1812187) q[1];
sx q[1];
rz(-2.6757619) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7458347) q[0];
sx q[0];
rz(-2.2767249) q[0];
sx q[0];
rz(-2.8559844) q[0];
x q[1];
rz(-0.1419576) q[2];
sx q[2];
rz(-2.3732911) q[2];
sx q[2];
rz(0.10266081) q[2];
rz(-pi) q[3];
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
rz(-pi) q[2];
rz(2.1773866) q[3];
sx q[3];
rz(-1.0830942) q[3];
sx q[3];
rz(-1.4473947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3992074) q[2];
sx q[2];
rz(-1.9611605) q[2];
sx q[2];
rz(0.36435374) q[2];
rz(0.24389167) q[3];
sx q[3];
rz(-0.049592169) q[3];
sx q[3];
rz(1.2172788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57109443) q[0];
sx q[0];
rz(-0.97935337) q[0];
sx q[0];
rz(-1.97557) q[0];
rz(-1.0150389) q[1];
sx q[1];
rz(-2.6045585) q[1];
sx q[1];
rz(-2.7551415) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5597661) q[0];
sx q[0];
rz(-2.92972) q[0];
sx q[0];
rz(1.7170402) q[0];
x q[1];
rz(-0.68604821) q[2];
sx q[2];
rz(-0.39271388) q[2];
sx q[2];
rz(-2.5401552) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4402953) q[1];
sx q[1];
rz(-0.079950182) q[1];
sx q[1];
rz(-1.246212) q[1];
rz(-pi) q[2];
rz(-1.3644045) q[3];
sx q[3];
rz(-2.5897621) q[3];
sx q[3];
rz(1.3760374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7361136) q[2];
sx q[2];
rz(-1.4961286) q[2];
sx q[2];
rz(0.22769895) q[2];
rz(-1.0730216) q[3];
sx q[3];
rz(-0.74879542) q[3];
sx q[3];
rz(0.29357114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
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
rz(-2.7135799) q[0];
rz(-1.0182861) q[1];
sx q[1];
rz(-0.4363474) q[1];
sx q[1];
rz(0.87103081) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7334325) q[0];
sx q[0];
rz(-2.1299612) q[0];
sx q[0];
rz(2.5446822) q[0];
x q[1];
rz(-1.0592878) q[2];
sx q[2];
rz(-2.7370484) q[2];
sx q[2];
rz(-1.900857) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4531252) q[1];
sx q[1];
rz(-0.18632132) q[1];
sx q[1];
rz(-2.2800355) q[1];
rz(-pi) q[2];
rz(2.4860367) q[3];
sx q[3];
rz(-2.3456315) q[3];
sx q[3];
rz(-1.9565005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.13228664) q[2];
sx q[2];
rz(-0.79557482) q[2];
sx q[2];
rz(-1.8629249) q[2];
rz(2.5631185) q[3];
sx q[3];
rz(-0.63251248) q[3];
sx q[3];
rz(-1.9875897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.4671675) q[0];
sx q[0];
rz(-2.9245057) q[0];
sx q[0];
rz(2.5417852) q[0];
rz(1.8425875) q[1];
sx q[1];
rz(-1.5312342) q[1];
sx q[1];
rz(-0.64186796) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80461795) q[0];
sx q[0];
rz(-2.0463159) q[0];
sx q[0];
rz(1.1853976) q[0];
rz(-pi) q[1];
rz(1.4207877) q[2];
sx q[2];
rz(-2.370655) q[2];
sx q[2];
rz(1.9341759) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.36343345) q[1];
sx q[1];
rz(-1.5892519) q[1];
sx q[1];
rz(-0.70779558) q[1];
x q[2];
rz(2.0976552) q[3];
sx q[3];
rz(-2.2552285) q[3];
sx q[3];
rz(-1.2795283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.32144) q[2];
sx q[2];
rz(-1.4250616) q[2];
sx q[2];
rz(0.610262) q[2];
rz(-0.52014703) q[3];
sx q[3];
rz(-2.0870233) q[3];
sx q[3];
rz(-0.49351969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34058061) q[0];
sx q[0];
rz(-1.2435253) q[0];
sx q[0];
rz(0.93801671) q[0];
rz(-2.7998789) q[1];
sx q[1];
rz(-1.7594756) q[1];
sx q[1];
rz(0.15728532) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.324824) q[0];
sx q[0];
rz(-0.95476645) q[0];
sx q[0];
rz(1.3321213) q[0];
rz(-pi) q[1];
rz(-3.0523275) q[2];
sx q[2];
rz(-2.3815739) q[2];
sx q[2];
rz(-1.6940728) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5648236) q[1];
sx q[1];
rz(-1.8722459) q[1];
sx q[1];
rz(0.5524854) q[1];
rz(-pi) q[2];
rz(-0.25453849) q[3];
sx q[3];
rz(-1.9677646) q[3];
sx q[3];
rz(1.0341259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9226795) q[2];
sx q[2];
rz(-2.7965386) q[2];
sx q[2];
rz(-1.2738796) q[2];
rz(-3.1359361) q[3];
sx q[3];
rz(-1.4414682) q[3];
sx q[3];
rz(-2.1605632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8511178) q[0];
sx q[0];
rz(-1.1707476) q[0];
sx q[0];
rz(0.049402417) q[0];
rz(-2.6750917) q[1];
sx q[1];
rz(-1.3460881) q[1];
sx q[1];
rz(0.64429611) q[1];
rz(-0.90032719) q[2];
sx q[2];
rz(-1.1569958) q[2];
sx q[2];
rz(-0.3668084) q[2];
rz(0.1771743) q[3];
sx q[3];
rz(-2.7037947) q[3];
sx q[3];
rz(-1.6537742) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
