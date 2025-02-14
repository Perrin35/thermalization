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
rz(-1.8981847) q[0];
sx q[0];
rz(-1.564448) q[0];
sx q[0];
rz(2.2301883) q[0];
rz(2.4721594) q[1];
sx q[1];
rz(-2.8652006) q[1];
sx q[1];
rz(-0.82958329) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2612635) q[0];
sx q[0];
rz(-1.1733049) q[0];
sx q[0];
rz(2.7842159) q[0];
rz(-pi) q[1];
rz(-1.6638512) q[2];
sx q[2];
rz(-1.8169799) q[2];
sx q[2];
rz(3.1174146) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1948283) q[1];
sx q[1];
rz(-2.2948537) q[1];
sx q[1];
rz(-1.1652693) q[1];
rz(1.4178278) q[3];
sx q[3];
rz(-2.0024096) q[3];
sx q[3];
rz(-0.6237517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.36898819) q[2];
sx q[2];
rz(-2.172894) q[2];
sx q[2];
rz(-1.9770835) q[2];
rz(0.72689593) q[3];
sx q[3];
rz(-1.8098857) q[3];
sx q[3];
rz(3.0708142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1233391) q[0];
sx q[0];
rz(-0.66272074) q[0];
sx q[0];
rz(0.29139274) q[0];
rz(-0.60028752) q[1];
sx q[1];
rz(-0.95708668) q[1];
sx q[1];
rz(-2.9765863) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24848973) q[0];
sx q[0];
rz(-1.4997673) q[0];
sx q[0];
rz(-1.8336589) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0290452) q[2];
sx q[2];
rz(-1.1243382) q[2];
sx q[2];
rz(-0.017306414) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2875322) q[1];
sx q[1];
rz(-0.89193501) q[1];
sx q[1];
rz(-1.4816585) q[1];
x q[2];
rz(1.6253774) q[3];
sx q[3];
rz(-1.3228058) q[3];
sx q[3];
rz(-0.33948201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0773641) q[2];
sx q[2];
rz(-0.95244971) q[2];
sx q[2];
rz(0.23951086) q[2];
rz(0.32660487) q[3];
sx q[3];
rz(-0.17290792) q[3];
sx q[3];
rz(3.0467564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25642446) q[0];
sx q[0];
rz(-0.435985) q[0];
sx q[0];
rz(-1.5727795) q[0];
rz(2.7478711) q[1];
sx q[1];
rz(-0.55501333) q[1];
sx q[1];
rz(0.47600019) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8825553) q[0];
sx q[0];
rz(-0.60094423) q[0];
sx q[0];
rz(-1.1968234) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4127998) q[2];
sx q[2];
rz(-0.18922986) q[2];
sx q[2];
rz(3.0410724) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4502628) q[1];
sx q[1];
rz(-2.1457377) q[1];
sx q[1];
rz(-2.6460669) q[1];
rz(-2.1952403) q[3];
sx q[3];
rz(-2.7331684) q[3];
sx q[3];
rz(-0.32865903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8126882) q[2];
sx q[2];
rz(-2.4220971) q[2];
sx q[2];
rz(-0.94181124) q[2];
rz(-1.4224667) q[3];
sx q[3];
rz(-0.37426451) q[3];
sx q[3];
rz(0.67371887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0734237) q[0];
sx q[0];
rz(-0.24003679) q[0];
sx q[0];
rz(1.6085251) q[0];
rz(0.34670058) q[1];
sx q[1];
rz(-1.0418714) q[1];
sx q[1];
rz(-2.9516721) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1153206) q[0];
sx q[0];
rz(-2.2824725) q[0];
sx q[0];
rz(1.9188966) q[0];
rz(-pi) q[1];
rz(2.8666977) q[2];
sx q[2];
rz(-1.1091386) q[2];
sx q[2];
rz(1.8732338) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.20224883) q[1];
sx q[1];
rz(-0.96572438) q[1];
sx q[1];
rz(-3.0720965) q[1];
x q[2];
rz(-0.81081049) q[3];
sx q[3];
rz(-0.56104198) q[3];
sx q[3];
rz(1.1277792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3612264) q[2];
sx q[2];
rz(-0.85005886) q[2];
sx q[2];
rz(-0.38193199) q[2];
rz(0.10968883) q[3];
sx q[3];
rz(-2.1083125) q[3];
sx q[3];
rz(-3.0221353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1031951) q[0];
sx q[0];
rz(-1.1829475) q[0];
sx q[0];
rz(-0.84684816) q[0];
rz(3.0021744) q[1];
sx q[1];
rz(-0.81652313) q[1];
sx q[1];
rz(-0.74904186) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7199166) q[0];
sx q[0];
rz(-0.70867175) q[0];
sx q[0];
rz(2.9139374) q[0];
rz(2.6566418) q[2];
sx q[2];
rz(-1.4451988) q[2];
sx q[2];
rz(-1.4664949) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.79758039) q[1];
sx q[1];
rz(-1.7645807) q[1];
sx q[1];
rz(-1.7001726) q[1];
rz(1.3588293) q[3];
sx q[3];
rz(-2.6304768) q[3];
sx q[3];
rz(1.8960475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2060812) q[2];
sx q[2];
rz(-1.8847909) q[2];
sx q[2];
rz(1.067266) q[2];
rz(-2.4308128) q[3];
sx q[3];
rz(-1.658541) q[3];
sx q[3];
rz(1.425364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-2.7988605) q[0];
sx q[0];
rz(-1.4787759) q[0];
sx q[0];
rz(-2.1730098) q[0];
rz(-2.4505278) q[1];
sx q[1];
rz(-1.7993118) q[1];
sx q[1];
rz(-1.8341281) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2742776) q[0];
sx q[0];
rz(-0.68570341) q[0];
sx q[0];
rz(-0.4088485) q[0];
rz(2.8505465) q[2];
sx q[2];
rz(-1.3851162) q[2];
sx q[2];
rz(0.03215511) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.56388748) q[1];
sx q[1];
rz(-1.9576549) q[1];
sx q[1];
rz(-3.0847286) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1627681) q[3];
sx q[3];
rz(-2.8205964) q[3];
sx q[3];
rz(-0.61883607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.4751733) q[2];
sx q[2];
rz(-0.23491493) q[2];
sx q[2];
rz(1.2972181) q[2];
rz(1.746486) q[3];
sx q[3];
rz(-1.8078943) q[3];
sx q[3];
rz(1.5313799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.34201) q[0];
sx q[0];
rz(-2.3260249) q[0];
sx q[0];
rz(2.2014501) q[0];
rz(-0.6483342) q[1];
sx q[1];
rz(-1.9377361) q[1];
sx q[1];
rz(-2.8727093) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9342182) q[0];
sx q[0];
rz(-1.8158127) q[0];
sx q[0];
rz(2.5859358) q[0];
rz(-2.4388695) q[2];
sx q[2];
rz(-1.3303192) q[2];
sx q[2];
rz(0.81288494) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.18792472) q[1];
sx q[1];
rz(-1.3410733) q[1];
sx q[1];
rz(2.1603895) q[1];
rz(-0.65461378) q[3];
sx q[3];
rz(-0.5236434) q[3];
sx q[3];
rz(2.3414827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1342643) q[2];
sx q[2];
rz(-1.8625883) q[2];
sx q[2];
rz(-2.4556665) q[2];
rz(-2.4053597) q[3];
sx q[3];
rz(-2.1015621) q[3];
sx q[3];
rz(2.9161016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69970423) q[0];
sx q[0];
rz(-1.2208953) q[0];
sx q[0];
rz(2.3928483) q[0];
rz(-3.0486095) q[1];
sx q[1];
rz(-2.3817606) q[1];
sx q[1];
rz(2.0704796) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1919928) q[0];
sx q[0];
rz(-1.1566391) q[0];
sx q[0];
rz(-2.5903775) q[0];
rz(0.65437859) q[2];
sx q[2];
rz(-1.9438231) q[2];
sx q[2];
rz(-1.6668863) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.86102137) q[1];
sx q[1];
rz(-1.264466) q[1];
sx q[1];
rz(2.2561399) q[1];
rz(0.052196189) q[3];
sx q[3];
rz(-1.9813271) q[3];
sx q[3];
rz(-0.13971381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6894655) q[2];
sx q[2];
rz(-2.5048544) q[2];
sx q[2];
rz(-2.9316736) q[2];
rz(3.012015) q[3];
sx q[3];
rz(-2.1942873) q[3];
sx q[3];
rz(1.5773704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.121948) q[0];
sx q[0];
rz(-2.8804998) q[0];
sx q[0];
rz(-3.1280532) q[0];
rz(-0.53567046) q[1];
sx q[1];
rz(-2.5745013) q[1];
sx q[1];
rz(0.013669107) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9122152) q[0];
sx q[0];
rz(-1.3446864) q[0];
sx q[0];
rz(-2.3588728) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7510468) q[2];
sx q[2];
rz(-2.0261814) q[2];
sx q[2];
rz(-2.9519338) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7175258) q[1];
sx q[1];
rz(-2.2560511) q[1];
sx q[1];
rz(1.1478893) q[1];
rz(-pi) q[2];
rz(-1.9433653) q[3];
sx q[3];
rz(-0.85597023) q[3];
sx q[3];
rz(-1.7320339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6232264) q[2];
sx q[2];
rz(-1.3806815) q[2];
sx q[2];
rz(1.3631442) q[2];
rz(-0.81950435) q[3];
sx q[3];
rz(-0.47911152) q[3];
sx q[3];
rz(-2.4578186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(2.0226456) q[0];
sx q[0];
rz(-0.85856694) q[0];
sx q[0];
rz(1.6534506) q[0];
rz(0.34291521) q[1];
sx q[1];
rz(-1.5838793) q[1];
sx q[1];
rz(0.75103474) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2011178) q[0];
sx q[0];
rz(-1.507326) q[0];
sx q[0];
rz(1.2918143) q[0];
x q[1];
rz(1.9575167) q[2];
sx q[2];
rz(-2.6474075) q[2];
sx q[2];
rz(1.3665716) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0768834) q[1];
sx q[1];
rz(-0.61184498) q[1];
sx q[1];
rz(0.20670273) q[1];
rz(-1.9225227) q[3];
sx q[3];
rz(-3.0455601) q[3];
sx q[3];
rz(2.072123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.47120961) q[2];
sx q[2];
rz(-0.94380108) q[2];
sx q[2];
rz(2.75441) q[2];
rz(0.47532982) q[3];
sx q[3];
rz(-1.1673704) q[3];
sx q[3];
rz(-0.79132426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8925856) q[0];
sx q[0];
rz(-2.1795166) q[0];
sx q[0];
rz(-2.4695061) q[0];
rz(0.36920209) q[1];
sx q[1];
rz(-0.75405706) q[1];
sx q[1];
rz(-1.3934607) q[1];
rz(-0.21239352) q[2];
sx q[2];
rz(-0.72827374) q[2];
sx q[2];
rz(-0.49401415) q[2];
rz(-1.9561097) q[3];
sx q[3];
rz(-2.8478619) q[3];
sx q[3];
rz(-0.0024402161) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
