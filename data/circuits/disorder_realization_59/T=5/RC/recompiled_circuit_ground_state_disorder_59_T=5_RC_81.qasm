OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6049603) q[0];
sx q[0];
rz(-2.9807615) q[0];
sx q[0];
rz(0.39128006) q[0];
rz(3.0985576) q[1];
sx q[1];
rz(-1.6208836) q[1];
sx q[1];
rz(-2.8032141) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7926377) q[0];
sx q[0];
rz(-0.39039877) q[0];
sx q[0];
rz(2.3304295) q[0];
rz(-pi) q[1];
rz(-1.7994405) q[2];
sx q[2];
rz(-2.9215339) q[2];
sx q[2];
rz(0.30753747) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.30792812) q[1];
sx q[1];
rz(-2.2349155) q[1];
sx q[1];
rz(-0.2425027) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8998221) q[3];
sx q[3];
rz(-0.42754506) q[3];
sx q[3];
rz(-2.0302142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.36838621) q[2];
sx q[2];
rz(-0.97969222) q[2];
sx q[2];
rz(-2.0835908) q[2];
rz(3.0392785) q[3];
sx q[3];
rz(-1.5240069) q[3];
sx q[3];
rz(1.4003096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8232987) q[0];
sx q[0];
rz(-2.3866391) q[0];
sx q[0];
rz(-0.21389432) q[0];
rz(-1.3338044) q[1];
sx q[1];
rz(-0.50024453) q[1];
sx q[1];
rz(-0.28876367) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64015111) q[0];
sx q[0];
rz(-0.94600618) q[0];
sx q[0];
rz(-0.21296091) q[0];
x q[1];
rz(0.66592543) q[2];
sx q[2];
rz(-2.1873432) q[2];
sx q[2];
rz(1.9279955) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.83486356) q[1];
sx q[1];
rz(-0.22278654) q[1];
sx q[1];
rz(1.7053717) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64677644) q[3];
sx q[3];
rz(-2.6138762) q[3];
sx q[3];
rz(2.3257252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5220962) q[2];
sx q[2];
rz(-0.82021004) q[2];
sx q[2];
rz(-1.2843457) q[2];
rz(-1.3075102) q[3];
sx q[3];
rz(-1.3176354) q[3];
sx q[3];
rz(-2.4431084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73806015) q[0];
sx q[0];
rz(-1.0665251) q[0];
sx q[0];
rz(-0.84529483) q[0];
rz(-2.2360133) q[1];
sx q[1];
rz(-0.60494345) q[1];
sx q[1];
rz(-1.1776498) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15525165) q[0];
sx q[0];
rz(-0.097106783) q[0];
sx q[0];
rz(-2.9349021) q[0];
x q[1];
rz(-2.6180116) q[2];
sx q[2];
rz(-2.058002) q[2];
sx q[2];
rz(-2.6139245) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.78598394) q[1];
sx q[1];
rz(-2.5289383) q[1];
sx q[1];
rz(2.775928) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2159141) q[3];
sx q[3];
rz(-1.4429907) q[3];
sx q[3];
rz(2.5296581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5251069) q[2];
sx q[2];
rz(-0.944204) q[2];
sx q[2];
rz(-0.13060972) q[2];
rz(-1.7836001) q[3];
sx q[3];
rz(-1.0087174) q[3];
sx q[3];
rz(1.2880026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65508771) q[0];
sx q[0];
rz(-0.26432744) q[0];
sx q[0];
rz(0.41410145) q[0];
rz(2.2762903) q[1];
sx q[1];
rz(-1.9046141) q[1];
sx q[1];
rz(0.6032595) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8651486) q[0];
sx q[0];
rz(-1.28946) q[0];
sx q[0];
rz(-1.9758609) q[0];
x q[1];
rz(0.19152051) q[2];
sx q[2];
rz(-2.7835803) q[2];
sx q[2];
rz(2.4409118) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8131667) q[1];
sx q[1];
rz(-0.35086497) q[1];
sx q[1];
rz(-3.0163832) q[1];
rz(-3.0217085) q[3];
sx q[3];
rz(-1.2601687) q[3];
sx q[3];
rz(-1.436123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6174751) q[2];
sx q[2];
rz(-1.6092499) q[2];
sx q[2];
rz(0.65940801) q[2];
rz(-0.13255969) q[3];
sx q[3];
rz(-2.5949251) q[3];
sx q[3];
rz(2.4373655) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34683126) q[0];
sx q[0];
rz(-2.1857388) q[0];
sx q[0];
rz(-0.26926789) q[0];
rz(1.5003834) q[1];
sx q[1];
rz(-1.8625448) q[1];
sx q[1];
rz(1.0795116) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0844042) q[0];
sx q[0];
rz(-1.3916172) q[0];
sx q[0];
rz(-2.4898131) q[0];
rz(-pi) q[1];
rz(-2.232021) q[2];
sx q[2];
rz(-0.9098297) q[2];
sx q[2];
rz(1.329042) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0725545) q[1];
sx q[1];
rz(-2.5423706) q[1];
sx q[1];
rz(2.1848382) q[1];
rz(-pi) q[2];
rz(-0.57860878) q[3];
sx q[3];
rz(-2.7670585) q[3];
sx q[3];
rz(-2.7079564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4214804) q[2];
sx q[2];
rz(-2.5613027) q[2];
sx q[2];
rz(0.87023467) q[2];
rz(1.8057711) q[3];
sx q[3];
rz(-1.3355052) q[3];
sx q[3];
rz(-0.67572063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0077591) q[0];
sx q[0];
rz(-3.1133911) q[0];
sx q[0];
rz(-0.61122417) q[0];
rz(-2.4080343) q[1];
sx q[1];
rz(-1.6707784) q[1];
sx q[1];
rz(-0.38331097) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72276536) q[0];
sx q[0];
rz(-1.2929898) q[0];
sx q[0];
rz(-2.4158143) q[0];
rz(2.3141642) q[2];
sx q[2];
rz(-1.7199923) q[2];
sx q[2];
rz(-2.0713446) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.53364175) q[1];
sx q[1];
rz(-1.4429394) q[1];
sx q[1];
rz(-0.92725384) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46054764) q[3];
sx q[3];
rz(-2.3604098) q[3];
sx q[3];
rz(2.6278842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40206259) q[2];
sx q[2];
rz(-0.5548839) q[2];
sx q[2];
rz(-1.0129207) q[2];
rz(-2.5439475) q[3];
sx q[3];
rz(-1.4089855) q[3];
sx q[3];
rz(-2.2156318) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0830072) q[0];
sx q[0];
rz(-0.42917955) q[0];
sx q[0];
rz(-3.10566) q[0];
rz(1.9195456) q[1];
sx q[1];
rz(-0.67600328) q[1];
sx q[1];
rz(-2.8935208) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0029468676) q[0];
sx q[0];
rz(-0.32297501) q[0];
sx q[0];
rz(2.2180411) q[0];
rz(-pi) q[1];
rz(-2.6509097) q[2];
sx q[2];
rz(-1.6958738) q[2];
sx q[2];
rz(-0.016029257) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.386302) q[1];
sx q[1];
rz(-2.7366126) q[1];
sx q[1];
rz(-1.1454357) q[1];
rz(-1.5682427) q[3];
sx q[3];
rz(-1.4658584) q[3];
sx q[3];
rz(1.5401287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0849358) q[2];
sx q[2];
rz(-2.0818905) q[2];
sx q[2];
rz(-0.53209957) q[2];
rz(-0.45446864) q[3];
sx q[3];
rz(-1.5957417) q[3];
sx q[3];
rz(-1.8475629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(0.71659511) q[0];
sx q[0];
rz(-1.0497365) q[0];
sx q[0];
rz(-0.24018921) q[0];
rz(-0.60513085) q[1];
sx q[1];
rz(-1.7938219) q[1];
sx q[1];
rz(-0.64632195) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28531269) q[0];
sx q[0];
rz(-2.3512771) q[0];
sx q[0];
rz(-1.1423443) q[0];
rz(1.4736299) q[2];
sx q[2];
rz(-2.7787913) q[2];
sx q[2];
rz(1.7513315) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.60257116) q[1];
sx q[1];
rz(-0.50316167) q[1];
sx q[1];
rz(-0.43518592) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2924932) q[3];
sx q[3];
rz(-0.67453803) q[3];
sx q[3];
rz(1.9703678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1120844) q[2];
sx q[2];
rz(-1.6479475) q[2];
sx q[2];
rz(0.081681577) q[2];
rz(2.2937842) q[3];
sx q[3];
rz(-0.40926465) q[3];
sx q[3];
rz(0.22296396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34143701) q[0];
sx q[0];
rz(-2.2969022) q[0];
sx q[0];
rz(0.89168125) q[0];
rz(1.910123) q[1];
sx q[1];
rz(-0.48858085) q[1];
sx q[1];
rz(-2.5317392) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.114074) q[0];
sx q[0];
rz(-2.4354012) q[0];
sx q[0];
rz(-0.8154517) q[0];
x q[1];
rz(-3.1270486) q[2];
sx q[2];
rz(-0.78476465) q[2];
sx q[2];
rz(-0.54569983) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.66850502) q[1];
sx q[1];
rz(-0.83946315) q[1];
sx q[1];
rz(1.1727612) q[1];
rz(-pi) q[2];
rz(-3.1248437) q[3];
sx q[3];
rz(-1.2257366) q[3];
sx q[3];
rz(2.6477388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8651198) q[2];
sx q[2];
rz(-2.1742994) q[2];
sx q[2];
rz(-2.8709732) q[2];
rz(-0.5451777) q[3];
sx q[3];
rz(-1.3278278) q[3];
sx q[3];
rz(2.8934208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4973064) q[0];
sx q[0];
rz(-1.8338642) q[0];
sx q[0];
rz(-0.18648952) q[0];
rz(-1.2093774) q[1];
sx q[1];
rz(-1.7861563) q[1];
sx q[1];
rz(0.019651042) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6015897) q[0];
sx q[0];
rz(-0.9603017) q[0];
sx q[0];
rz(-2.6363346) q[0];
rz(-pi) q[1];
rz(2.8803359) q[2];
sx q[2];
rz(-1.3105596) q[2];
sx q[2];
rz(2.2360436) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3557089) q[1];
sx q[1];
rz(-1.9221091) q[1];
sx q[1];
rz(-1.7339049) q[1];
rz(2.4830706) q[3];
sx q[3];
rz(-2.2611475) q[3];
sx q[3];
rz(0.041680574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.43202117) q[2];
sx q[2];
rz(-0.95546466) q[2];
sx q[2];
rz(-1.5891937) q[2];
rz(-3.0324557) q[3];
sx q[3];
rz(-1.8692632) q[3];
sx q[3];
rz(-0.2763589) q[3];
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
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6266262) q[0];
sx q[0];
rz(-1.6635386) q[0];
sx q[0];
rz(-1.2141825) q[0];
rz(-1.4151489) q[1];
sx q[1];
rz(-1.0217923) q[1];
sx q[1];
rz(1.6820977) q[1];
rz(-0.44533163) q[2];
sx q[2];
rz(-1.3117322) q[2];
sx q[2];
rz(-1.7016254) q[2];
rz(0.23701238) q[3];
sx q[3];
rz(-2.7494299) q[3];
sx q[3];
rz(1.3361479) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
