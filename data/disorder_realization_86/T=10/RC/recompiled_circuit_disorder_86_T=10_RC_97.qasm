OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3671626) q[0];
sx q[0];
rz(-2.2280333) q[0];
sx q[0];
rz(1.7295184) q[0];
rz(-2.9867759) q[1];
sx q[1];
rz(5.6875416) q[1];
sx q[1];
rz(7.7653801) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9528113) q[0];
sx q[0];
rz(-2.6872098) q[0];
sx q[0];
rz(-0.36034583) q[0];
x q[1];
rz(0.46967311) q[2];
sx q[2];
rz(-2.854752) q[2];
sx q[2];
rz(-1.4925721) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.60018051) q[1];
sx q[1];
rz(-1.0804847) q[1];
sx q[1];
rz(-2.6580826) q[1];
rz(-pi) q[2];
rz(-1.7300624) q[3];
sx q[3];
rz(-2.5730238) q[3];
sx q[3];
rz(-1.1390613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1564864) q[2];
sx q[2];
rz(-0.50922314) q[2];
sx q[2];
rz(-2.2757754) q[2];
rz(-2.1872897) q[3];
sx q[3];
rz(-1.6031957) q[3];
sx q[3];
rz(-1.8538063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1433379) q[0];
sx q[0];
rz(-1.4366432) q[0];
sx q[0];
rz(3.1153733) q[0];
rz(1.6014618) q[1];
sx q[1];
rz(-1.5427579) q[1];
sx q[1];
rz(-0.96347934) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.026982633) q[0];
sx q[0];
rz(-2.5275505) q[0];
sx q[0];
rz(-1.5747889) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1396779) q[2];
sx q[2];
rz(-1.2895673) q[2];
sx q[2];
rz(2.1525454) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7746437) q[1];
sx q[1];
rz(-1.0732713) q[1];
sx q[1];
rz(-0.36555396) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96879949) q[3];
sx q[3];
rz(-2.7235944) q[3];
sx q[3];
rz(-0.19817643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6271237) q[2];
sx q[2];
rz(-1.1273948) q[2];
sx q[2];
rz(3.0070686) q[2];
rz(-2.3965805) q[3];
sx q[3];
rz(-2.9146505) q[3];
sx q[3];
rz(2.1988595) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9298252) q[0];
sx q[0];
rz(-0.38914248) q[0];
sx q[0];
rz(-2.3441558) q[0];
rz(-1.047661) q[1];
sx q[1];
rz(-0.14973775) q[1];
sx q[1];
rz(2.581596) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6277498) q[0];
sx q[0];
rz(-1.7729513) q[0];
sx q[0];
rz(-1.1802799) q[0];
rz(-pi) q[1];
rz(-0.30324869) q[2];
sx q[2];
rz(-1.5911284) q[2];
sx q[2];
rz(-3.0603527) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.33369505) q[1];
sx q[1];
rz(-1.9296608) q[1];
sx q[1];
rz(-1.3586504) q[1];
rz(-pi) q[2];
x q[2];
rz(0.47645724) q[3];
sx q[3];
rz(-2.0072848) q[3];
sx q[3];
rz(0.95526327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3893163) q[2];
sx q[2];
rz(-1.1976778) q[2];
sx q[2];
rz(0.17253549) q[2];
rz(2.1595188) q[3];
sx q[3];
rz(-1.3970102) q[3];
sx q[3];
rz(2.0836232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.003222) q[0];
sx q[0];
rz(-1.0976185) q[0];
sx q[0];
rz(-0.28451434) q[0];
rz(0.31670397) q[1];
sx q[1];
rz(-2.7088294) q[1];
sx q[1];
rz(-1.2987312) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0239149) q[0];
sx q[0];
rz(-2.0664584) q[0];
sx q[0];
rz(2.8850874) q[0];
rz(-0.0012871731) q[2];
sx q[2];
rz(-2.3372071) q[2];
sx q[2];
rz(3.121701) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.87151566) q[1];
sx q[1];
rz(-1.75711) q[1];
sx q[1];
rz(0.99888505) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.11233791) q[3];
sx q[3];
rz(-1.897052) q[3];
sx q[3];
rz(0.60409594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6056885) q[2];
sx q[2];
rz(-2.9457592) q[2];
sx q[2];
rz(-0.38468012) q[2];
rz(0.7540594) q[3];
sx q[3];
rz(-1.056517) q[3];
sx q[3];
rz(1.6872905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6032747) q[0];
sx q[0];
rz(-0.98709995) q[0];
sx q[0];
rz(-1.7549365) q[0];
rz(-0.23100135) q[1];
sx q[1];
rz(-1.341154) q[1];
sx q[1];
rz(2.8447661) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7291527) q[0];
sx q[0];
rz(-1.531633) q[0];
sx q[0];
rz(2.5705283) q[0];
rz(-pi) q[1];
x q[1];
rz(0.74283959) q[2];
sx q[2];
rz(-1.683871) q[2];
sx q[2];
rz(-2.6807221) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2424803) q[1];
sx q[1];
rz(-2.50933) q[1];
sx q[1];
rz(2.2805023) q[1];
rz(-pi) q[2];
rz(-1.6597219) q[3];
sx q[3];
rz(-0.85907912) q[3];
sx q[3];
rz(-0.54642788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7632873) q[2];
sx q[2];
rz(-1.8313235) q[2];
sx q[2];
rz(2.7491167) q[2];
rz(1.9893507) q[3];
sx q[3];
rz(-0.71458721) q[3];
sx q[3];
rz(2.8241482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1257989) q[0];
sx q[0];
rz(-1.5725461) q[0];
sx q[0];
rz(-2.3902067) q[0];
rz(1.8136576) q[1];
sx q[1];
rz(-1.8782047) q[1];
sx q[1];
rz(-2.5352535) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0814708) q[0];
sx q[0];
rz(-1.524964) q[0];
sx q[0];
rz(1.5541374) q[0];
rz(-pi) q[1];
rz(-2.7251284) q[2];
sx q[2];
rz(-0.93566862) q[2];
sx q[2];
rz(-3.0481899) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6341083) q[1];
sx q[1];
rz(-2.1180696) q[1];
sx q[1];
rz(-2.9299111) q[1];
x q[2];
rz(-1.243152) q[3];
sx q[3];
rz(-2.9122105) q[3];
sx q[3];
rz(0.71654191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5027344) q[2];
sx q[2];
rz(-2.0998462) q[2];
sx q[2];
rz(1.139337) q[2];
rz(-1.4849439) q[3];
sx q[3];
rz(-1.1805725) q[3];
sx q[3];
rz(-3.0373354) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6095603) q[0];
sx q[0];
rz(-2.39344) q[0];
sx q[0];
rz(-2.6334921) q[0];
rz(-1.5787026) q[1];
sx q[1];
rz(-2.0527614) q[1];
sx q[1];
rz(0.79024822) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9433141) q[0];
sx q[0];
rz(-2.4662848) q[0];
sx q[0];
rz(-2.6576463) q[0];
rz(-pi) q[1];
rz(-0.054025606) q[2];
sx q[2];
rz(-0.21394357) q[2];
sx q[2];
rz(2.8097048) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.33752791) q[1];
sx q[1];
rz(-1.3009326) q[1];
sx q[1];
rz(-0.40562628) q[1];
rz(-pi) q[2];
rz(1.948445) q[3];
sx q[3];
rz(-2.1834063) q[3];
sx q[3];
rz(-0.8876422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0885075) q[2];
sx q[2];
rz(-0.44184703) q[2];
sx q[2];
rz(1.7283758) q[2];
rz(-1.6648071) q[3];
sx q[3];
rz(-1.0352742) q[3];
sx q[3];
rz(0.061554519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7168032) q[0];
sx q[0];
rz(-3.1112818) q[0];
sx q[0];
rz(-1.0472263) q[0];
rz(-2.5324902) q[1];
sx q[1];
rz(-1.7276238) q[1];
sx q[1];
rz(-1.75288) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4980709) q[0];
sx q[0];
rz(-2.5812491) q[0];
sx q[0];
rz(1.6269006) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6128204) q[2];
sx q[2];
rz(-1.1985589) q[2];
sx q[2];
rz(-2.1793274) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.74275201) q[1];
sx q[1];
rz(-1.5821777) q[1];
sx q[1];
rz(-0.091817261) q[1];
x q[2];
rz(2.1771168) q[3];
sx q[3];
rz(-1.8165605) q[3];
sx q[3];
rz(2.5442459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1887112) q[2];
sx q[2];
rz(-2.7313576) q[2];
sx q[2];
rz(-0.88225538) q[2];
rz(-1.7404209) q[3];
sx q[3];
rz(-1.1663576) q[3];
sx q[3];
rz(1.9410979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27154487) q[0];
sx q[0];
rz(-0.4168059) q[0];
sx q[0];
rz(-1.4260938) q[0];
rz(0.081461279) q[1];
sx q[1];
rz(-1.1625682) q[1];
sx q[1];
rz(-2.5833599) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0575858) q[0];
sx q[0];
rz(-1.1450197) q[0];
sx q[0];
rz(0.42752479) q[0];
x q[1];
rz(2.0360557) q[2];
sx q[2];
rz(-3.0214587) q[2];
sx q[2];
rz(2.5533954) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7575175) q[1];
sx q[1];
rz(-1.6346524) q[1];
sx q[1];
rz(0.82660316) q[1];
x q[2];
rz(0.48645143) q[3];
sx q[3];
rz(-1.7351741) q[3];
sx q[3];
rz(-2.1742976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2925064) q[2];
sx q[2];
rz(-1.8722653) q[2];
sx q[2];
rz(-0.212184) q[2];
rz(2.9296181) q[3];
sx q[3];
rz(-2.4583355) q[3];
sx q[3];
rz(-1.9395444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6367209) q[0];
sx q[0];
rz(-2.3175406) q[0];
sx q[0];
rz(-1.5378392) q[0];
rz(0.82540712) q[1];
sx q[1];
rz(-0.67276612) q[1];
sx q[1];
rz(-0.5232946) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16738811) q[0];
sx q[0];
rz(-0.61199576) q[0];
sx q[0];
rz(-3.0332855) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.043838219) q[2];
sx q[2];
rz(-0.99928108) q[2];
sx q[2];
rz(-0.29585719) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.20269468) q[1];
sx q[1];
rz(-1.885186) q[1];
sx q[1];
rz(2.7640192) q[1];
rz(-0.2449805) q[3];
sx q[3];
rz(-1.3462726) q[3];
sx q[3];
rz(0.45104879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3045197) q[2];
sx q[2];
rz(-0.7154811) q[2];
sx q[2];
rz(2.8722897) q[2];
rz(-2.6473911) q[3];
sx q[3];
rz(-0.84635693) q[3];
sx q[3];
rz(-1.0860898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3257278) q[0];
sx q[0];
rz(-1.6115191) q[0];
sx q[0];
rz(1.4900526) q[0];
rz(1.6745463) q[1];
sx q[1];
rz(-0.29232262) q[1];
sx q[1];
rz(1.2437337) q[1];
rz(1.3080296) q[2];
sx q[2];
rz(-1.5816214) q[2];
sx q[2];
rz(-2.66253) q[2];
rz(-0.79694637) q[3];
sx q[3];
rz(-1.7399825) q[3];
sx q[3];
rz(0.83132838) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
