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
rz(0.35667625) q[0];
sx q[0];
rz(-1.6950322) q[0];
sx q[0];
rz(2.2398228) q[0];
rz(-1.8094485) q[1];
sx q[1];
rz(-0.44337115) q[1];
sx q[1];
rz(-2.3763357) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0888521) q[0];
sx q[0];
rz(-0.67702261) q[0];
sx q[0];
rz(2.1734326) q[0];
rz(-pi) q[1];
rz(2.1852605) q[2];
sx q[2];
rz(-0.40618375) q[2];
sx q[2];
rz(1.150591) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.87188497) q[1];
sx q[1];
rz(-2.800834) q[1];
sx q[1];
rz(-0.12026708) q[1];
x q[2];
rz(1.8497068) q[3];
sx q[3];
rz(-1.6959012) q[3];
sx q[3];
rz(-2.5549614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.29609933) q[2];
sx q[2];
rz(-1.0881492) q[2];
sx q[2];
rz(-1.4720526) q[2];
rz(-1.4269525) q[3];
sx q[3];
rz(-1.6441556) q[3];
sx q[3];
rz(0.53372395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.958309) q[0];
sx q[0];
rz(-1.6338209) q[0];
sx q[0];
rz(2.2928152) q[0];
rz(0.035004184) q[1];
sx q[1];
rz(-1.9947546) q[1];
sx q[1];
rz(2.1880207) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35915056) q[0];
sx q[0];
rz(-1.1138565) q[0];
sx q[0];
rz(-0.63553973) q[0];
x q[1];
rz(-1.9601279) q[2];
sx q[2];
rz(-1.2876533) q[2];
sx q[2];
rz(2.2623537) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9417022) q[1];
sx q[1];
rz(-2.7459956) q[1];
sx q[1];
rz(-2.0989749) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3869629) q[3];
sx q[3];
rz(-2.0622232) q[3];
sx q[3];
rz(1.400465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1635052) q[2];
sx q[2];
rz(-1.3601902) q[2];
sx q[2];
rz(-0.84954849) q[2];
rz(-0.16544011) q[3];
sx q[3];
rz(-1.4435507) q[3];
sx q[3];
rz(0.74976841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36990944) q[0];
sx q[0];
rz(-2.219438) q[0];
sx q[0];
rz(-2.7396696) q[0];
rz(0.15332128) q[1];
sx q[1];
rz(-2.9647277) q[1];
sx q[1];
rz(-2.0053999) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72969681) q[0];
sx q[0];
rz(-1.9703016) q[0];
sx q[0];
rz(-3.1285498) q[0];
rz(-pi) q[1];
x q[1];
rz(0.49164518) q[2];
sx q[2];
rz(-0.53722135) q[2];
sx q[2];
rz(1.1342837) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.043667533) q[1];
sx q[1];
rz(-0.33722028) q[1];
sx q[1];
rz(-2.340576) q[1];
rz(-pi) q[2];
rz(1.1785169) q[3];
sx q[3];
rz(-1.137737) q[3];
sx q[3];
rz(-2.1585502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.084126964) q[2];
sx q[2];
rz(-0.87856138) q[2];
sx q[2];
rz(0.055056661) q[2];
rz(-1.347524) q[3];
sx q[3];
rz(-1.3301347) q[3];
sx q[3];
rz(1.0042892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7771626) q[0];
sx q[0];
rz(-2.936852) q[0];
sx q[0];
rz(-1.5740016) q[0];
rz(2.4500627) q[1];
sx q[1];
rz(-2.4722996) q[1];
sx q[1];
rz(-2.9561668) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90339336) q[0];
sx q[0];
rz(-1.967634) q[0];
sx q[0];
rz(-1.3893886) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9965956) q[2];
sx q[2];
rz(-1.123872) q[2];
sx q[2];
rz(1.9831234) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1048829) q[1];
sx q[1];
rz(-1.8817025) q[1];
sx q[1];
rz(-2.8757812) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7774277) q[3];
sx q[3];
rz(-2.1071599) q[3];
sx q[3];
rz(0.11851507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7376248) q[2];
sx q[2];
rz(-1.4289958) q[2];
sx q[2];
rz(3.1363623) q[2];
rz(2.7541006) q[3];
sx q[3];
rz(-0.5046851) q[3];
sx q[3];
rz(-1.8385878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4150647) q[0];
sx q[0];
rz(-0.98353493) q[0];
sx q[0];
rz(-0.61597419) q[0];
rz(1.1888602) q[1];
sx q[1];
rz(-0.6183466) q[1];
sx q[1];
rz(2.628285) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5145288) q[0];
sx q[0];
rz(-1.5781869) q[0];
sx q[0];
rz(1.8523995) q[0];
rz(-0.31536021) q[2];
sx q[2];
rz(-0.7509481) q[2];
sx q[2];
rz(2.2436004) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.087438093) q[1];
sx q[1];
rz(-2.5028924) q[1];
sx q[1];
rz(3.020546) q[1];
rz(-pi) q[2];
rz(-1.3393976) q[3];
sx q[3];
rz(-2.6632892) q[3];
sx q[3];
rz(-1.4372642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0080228) q[2];
sx q[2];
rz(-1.2761389) q[2];
sx q[2];
rz(2.3991154) q[2];
rz(-1.6978469) q[3];
sx q[3];
rz(-1.7323114) q[3];
sx q[3];
rz(1.1013364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83108574) q[0];
sx q[0];
rz(-1.0467014) q[0];
sx q[0];
rz(2.781784) q[0];
rz(-0.85044914) q[1];
sx q[1];
rz(-1.1812187) q[1];
sx q[1];
rz(-2.6757619) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98679989) q[0];
sx q[0];
rz(-1.354711) q[0];
sx q[0];
rz(0.84439069) q[0];
rz(-pi) q[1];
rz(2.9996351) q[2];
sx q[2];
rz(-0.76830155) q[2];
sx q[2];
rz(-0.10266081) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.40738525) q[1];
sx q[1];
rz(-1.3300597) q[1];
sx q[1];
rz(1.5182785) q[1];
rz(-pi) q[2];
rz(0.82139246) q[3];
sx q[3];
rz(-2.3830049) q[3];
sx q[3];
rz(0.71739355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3992074) q[2];
sx q[2];
rz(-1.9611605) q[2];
sx q[2];
rz(-0.36435374) q[2];
rz(-2.897701) q[3];
sx q[3];
rz(-3.0920005) q[3];
sx q[3];
rz(-1.2172788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57109443) q[0];
sx q[0];
rz(-0.97935337) q[0];
sx q[0];
rz(-1.1660227) q[0];
rz(-2.1265538) q[1];
sx q[1];
rz(-0.53703419) q[1];
sx q[1];
rz(0.38645116) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58182654) q[0];
sx q[0];
rz(-2.92972) q[0];
sx q[0];
rz(1.4245524) q[0];
x q[1];
rz(-0.68604821) q[2];
sx q[2];
rz(-2.7488788) q[2];
sx q[2];
rz(2.5401552) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7658479) q[1];
sx q[1];
rz(-1.4950291) q[1];
sx q[1];
rz(-3.1160464) q[1];
rz(-pi) q[2];
rz(-0.12549716) q[3];
sx q[3];
rz(-1.0319796) q[3];
sx q[3];
rz(-1.6171238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.40547907) q[2];
sx q[2];
rz(-1.6454641) q[2];
sx q[2];
rz(-0.22769895) q[2];
rz(-2.0685711) q[3];
sx q[3];
rz(-2.3927972) q[3];
sx q[3];
rz(-2.8480215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
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
rz(-1.0182861) q[1];
sx q[1];
rz(-2.7052453) q[1];
sx q[1];
rz(-0.87103081) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6329417) q[0];
sx q[0];
rz(-2.067446) q[0];
sx q[0];
rz(-2.2185241) q[0];
rz(-pi) q[1];
rz(1.0592878) q[2];
sx q[2];
rz(-0.40454421) q[2];
sx q[2];
rz(1.2407357) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.58295137) q[1];
sx q[1];
rz(-1.691733) q[1];
sx q[1];
rz(1.7128829) q[1];
rz(-pi) q[2];
rz(-0.68058859) q[3];
sx q[3];
rz(-1.1201123) q[3];
sx q[3];
rz(0.87928538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.13228664) q[2];
sx q[2];
rz(-2.3460178) q[2];
sx q[2];
rz(1.2786678) q[2];
rz(-2.5631185) q[3];
sx q[3];
rz(-0.63251248) q[3];
sx q[3];
rz(-1.1540029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-1.4671675) q[0];
sx q[0];
rz(-0.217087) q[0];
sx q[0];
rz(-2.5417852) q[0];
rz(1.8425875) q[1];
sx q[1];
rz(-1.5312342) q[1];
sx q[1];
rz(2.4997247) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.079499809) q[0];
sx q[0];
rz(-2.5389517) q[0];
sx q[0];
rz(2.5109765) q[0];
rz(-0.80550216) q[2];
sx q[2];
rz(-1.4664716) q[2];
sx q[2];
rz(0.25539216) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9126666) q[1];
sx q[1];
rz(-0.70799457) q[1];
sx q[1];
rz(0.028381746) q[1];
x q[2];
rz(-0.55223744) q[3];
sx q[3];
rz(-2.304616) q[3];
sx q[3];
rz(2.6058634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.82015264) q[2];
sx q[2];
rz(-1.4250616) q[2];
sx q[2];
rz(-2.5313306) q[2];
rz(2.6214456) q[3];
sx q[3];
rz(-2.0870233) q[3];
sx q[3];
rz(2.648073) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.801012) q[0];
sx q[0];
rz(-1.2435253) q[0];
sx q[0];
rz(-0.93801671) q[0];
rz(2.7998789) q[1];
sx q[1];
rz(-1.382117) q[1];
sx q[1];
rz(-2.9843073) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7559497) q[0];
sx q[0];
rz(-1.7649714) q[0];
sx q[0];
rz(-2.5119378) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6553272) q[2];
sx q[2];
rz(-2.3270375) q[2];
sx q[2];
rz(1.5712122) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5648236) q[1];
sx q[1];
rz(-1.2693468) q[1];
sx q[1];
rz(2.5891073) q[1];
x q[2];
rz(2.1116793) q[3];
sx q[3];
rz(-0.4678886) q[3];
sx q[3];
rz(-0.44178007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.21891317) q[2];
sx q[2];
rz(-2.7965386) q[2];
sx q[2];
rz(-1.867713) q[2];
rz(-0.0056565469) q[3];
sx q[3];
rz(-1.7001245) q[3];
sx q[3];
rz(-2.1605632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2904749) q[0];
sx q[0];
rz(-1.1707476) q[0];
sx q[0];
rz(0.049402417) q[0];
rz(0.46650096) q[1];
sx q[1];
rz(-1.3460881) q[1];
sx q[1];
rz(0.64429611) q[1];
rz(-0.90032719) q[2];
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
